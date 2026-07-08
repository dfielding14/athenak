from __future__ import annotations

import os
from pathlib import Path
import shlex
import shutil
import subprocess
import textwrap

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def _function_body(source: str, signature: str) -> str:
    start = source.index(signature)
    brace = source.index("{", start)
    depth = 0
    for index in range(brace, len(source)):
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                return source[brace:index + 1]
    raise AssertionError(f"unterminated function: {signature}")


def test_migration_uses_one_validated_hole_set_and_batched_move() -> None:
    source = (REPO_ROOT / "src/bvals/bvals_part.cpp").read_text(encoding="ascii")
    body = _function_body(
        source, "TaskStatus ParticlesBoundaryValues::RecvAndUnpackPrtcls()"
    )
    assert "BuildParticleCompactionPlan" in body
    assert "particle_compact_survivors" in body
    assert "particle_holes" in body
    assert "nremain_d" not in body
    assert "Kokkos::deep_copy(rdest, rsrc)" not in body


def test_compaction_plan_preserves_exact_survivor_identity(tmp_path: Path) -> None:
    compiler_command = shlex.split(os.environ.get("CXX", "c++"))
    compiler = shutil.which(compiler_command[0])
    if compiler is None:
        pytest.skip("a C++ compiler is required for the compaction planner unit test")
    compiler_command[0] = compiler

    harness = tmp_path / "particle_compaction_test.cpp"
    executable = tmp_path / "particle_compaction_test"
    harness.write_text(
        textwrap.dedent(
            r"""
            #include <algorithm>
            #include <iostream>
            #include <string>
            #include <vector>

            #include "src/bvals/particle_compaction.hpp"

            namespace {

            bool CheckCase(int old_size, const std::vector<int> &send,
                           const std::vector<int> &destroy, int receive_count) {
              particles::ParticleCompactionPlan plan;
              std::string error;
              if (!particles::BuildParticleCompactionPlan(
                      old_size, receive_count, send, destroy, &plan, &error)) {
                std::cerr << "unexpected planning failure: " << error << '\n';
                return false;
              }

              std::vector<int> state(old_size);
              for (int n = 0; n < old_size; ++n) state[n] = n;
              if (plan.final_size > old_size) state.resize(plan.final_size, -9999);
              for (int n = 0; n < receive_count; ++n) {
                const int destination =
                    (n < static_cast<int>(plan.holes.size()))
                        ? plan.holes[n]
                        : old_size + n - static_cast<int>(plan.holes.size());
                state[destination] = old_size + n;
              }
              for (std::size_t n = 0; n < plan.move_sources.size(); ++n) {
                if (plan.move_sources[n] < plan.final_size ||
                    plan.move_destinations[n] >= plan.final_size) {
                  std::cerr << "source and destination ranges overlap\n";
                  return false;
                }
                state[plan.move_destinations[n]] = state[plan.move_sources[n]];
              }
              state.resize(plan.final_size);

              std::vector<bool> removed(old_size, false);
              for (const int index : send) removed[index] = true;
              for (const int index : destroy) removed[index] = true;
              std::vector<int> expected;
              for (int n = 0; n < old_size; ++n) {
                if (!removed[n]) expected.push_back(n);
              }
              for (int n = 0; n < receive_count; ++n) {
                expected.push_back(old_size + n);
              }
              std::sort(state.begin(), state.end());
              std::sort(expected.begin(), expected.end());
              if (state != expected) {
                std::cerr << "survivor identity mismatch\n";
                return false;
              }
              return true;
            }

            }  // namespace

            int main() {
              // Exhaustively assign each original particle to stay/send/destroy.
              for (int old_size = 0; old_size <= 8; ++old_size) {
                int configurations = 1;
                for (int n = 0; n < old_size; ++n) configurations *= 3;
                for (int encoded = 0; encoded < configurations; ++encoded) {
                  int value = encoded;
                  std::vector<int> send;
                  std::vector<int> destroy;
                  for (int particle = 0; particle < old_size; ++particle) {
                    const int fate = value % 3;
                    value /= 3;
                    if (fate == 1) send.push_back(particle);
                    if (fate == 2) destroy.push_back(particle);
                  }
                  for (int receive_count = 0; receive_count <= old_size + 2;
                       ++receive_count) {
                    if (!CheckCase(old_size, send, destroy, receive_count)) return 1;
                    std::reverse(send.begin(), send.end());
                    std::reverse(destroy.begin(), destroy.end());
                    if (!CheckCase(old_size, send, destroy, receive_count)) return 7;
                  }
                }
              }

              particles::ParticleCompactionPlan plan;
              std::string error;
              if (!particles::BuildParticleCompactionPlan(
                      10, 0, {1}, {9}, &plan, &error)) return 2;
              if (plan.move_sources != std::vector<int>{8} ||
                  plan.move_destinations != std::vector<int>{1}) return 3;
              if (particles::BuildParticleCompactionPlan(
                      3, 0, {1}, {1}, &plan, &error)) return 4;
              if (particles::BuildParticleCompactionPlan(
                      3, 0, {3}, {}, &plan, &error)) return 5;
              if (particles::BuildParticleCompactionPlan(
                      -1, 0, {}, {}, &plan, &error)) return 6;
              return 0;
            }
            """
        ),
        encoding="ascii",
    )
    compile_result = subprocess.run(
        [
            *compiler_command,
            "-std=c++17",
            "-Wall",
            "-Wextra",
            "-pedantic",
            f"-I{REPO_ROOT}",
            str(harness),
            "-o",
            str(executable),
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    assert compile_result.returncode == 0, compile_result.stdout

    run_result = subprocess.run(
        [str(executable)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    assert run_result.returncode == 0, run_result.stdout
