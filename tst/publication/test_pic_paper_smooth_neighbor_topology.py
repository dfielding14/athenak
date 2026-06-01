#!/usr/bin/env python3
"""Focused topology checks for paper_smooth periodic receiver routing."""

from __future__ import annotations

from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]


def _neighbor_index(ix: int, iy: int, iz: int, n1: int, n2: int) -> int:
    if abs(ix) + abs(iy) + abs(iz) == 0 or abs(ix * iy * iz) > 1:
        return -1
    if iz == 0:
        if ix * iy == 0:
            return abs(ix) * 2 * (ix + 1) + abs(iy) * 2 * (iy + 5) + n1 + 2 * n2
        return 16 + (ix + 1) + 2 * (iy + 1) + n1
    if ix * iy == 0:
        return (
            24
            + abs(ix) * (ix + 9)
            + abs(iy) * (iy + 17)
            + 2 * (iz + 1)
            + n1
            + 2 * n2
        )
    return 48 + (ix + 1) // 2 + (iy + 1) + 2 * (iz + 1)


def _decode_legal_slots() -> dict[int, tuple[int, int, int]]:
    decoded = {}
    for iz in (-1, 0, 1):
        for iy in (-1, 0, 1):
            for ix in (-1, 0, 1):
                if (ix, iy, iz) == (0, 0, 0):
                    continue
                zeros = (ix == 0) + (iy == 0) + (iz == 0)
                for n1 in range(2 if zeros >= 1 else 1):
                    for n2 in range(2 if zeros == 2 else 1):
                        slot = _neighbor_index(ix, iy, iz, n1, n2)
                        prior = decoded.setdefault(slot, (ix, iy, iz))
                        if prior != (ix, iy, iz):
                            raise AssertionError(
                                f"slot {slot} aliases {prior} and {(ix, iy, iz)}"
                            )
    return decoded


class PicPaperSmoothNeighborTopologyTests(unittest.TestCase):
    def test_all_56_legal_slots_decode_uniquely(self) -> None:
        decoded = _decode_legal_slots()
        self.assertEqual(sorted(decoded), list(range(56)))
        expected_x3_edges = {
            32: (-1, 0, -1),
            33: (-1, 0, -1),
            34: (1, 0, -1),
            35: (1, 0, -1),
            36: (-1, 0, 1),
            37: (-1, 0, 1),
            38: (1, 0, 1),
            39: (1, 0, 1),
            40: (0, -1, -1),
            41: (0, -1, -1),
            42: (0, 1, -1),
            43: (0, 1, -1),
            44: (0, -1, 1),
            45: (0, -1, 1),
            46: (0, 1, 1),
            47: (0, 1, 1),
        }
        self.assertEqual({slot: decoded[slot] for slot in expected_x3_edges},
                         expected_x3_edges)

    def test_cpp_decoder_uses_only_legal_subdivision_arguments(self) -> None:
        source = (
            REPO_ROOT / "src/particles/particles_moments.cpp"
        ).read_text(encoding="ascii")
        self.assertIn(
            "const int n1_max = (ndirections >= 1) ? 1 : 0;",
            source,
        )
        self.assertIn(
            "const int n2_max = (ndirections == 2) ? 1 : 0;",
            source,
        )

    def test_transport_identity_includes_periodic_image_code(self) -> None:
        source = (REPO_ROOT / "src/bvals/bvals_mom.cpp").read_text(encoding="ascii")
        self.assertIn("record.reserved", source)
        self.assertIn("lhs.image_code == rhs.image_code", source)

    def test_record_scaling_matches_standard_deltaf_deposition(self) -> None:
        source = (
            REPO_ROOT / "src/particles/particles_moments.cpp"
        ).read_text(encoding="ascii")
        for snippet in (
            "physical_boundary_scale*deposit_qscale*weight*df_weight*h_qspecies(sp)",
            "vx, vy, vz, physical_boundary_scale*h_pr(IPEBDOT, p)",
            "physical_boundary_scale*df_weight*h_pr(IPDPX, p)",
            "physical_boundary_scale*df_weight*h_pr(IPDPY, p)",
            "physical_boundary_scale*df_weight*h_pr(IPDPZ, p)",
            "physical_boundary_scale*df_weight*h_pr(IPDE, p)",
        ):
            with self.subTest(snippet=snippet):
                self.assertIn(snippet, source)


if __name__ == "__main__":
    unittest.main()
