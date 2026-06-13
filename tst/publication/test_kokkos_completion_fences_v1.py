import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
MAIN = REPO_ROOT / "src/main.cpp"
BVALS = REPO_ROOT / "src/bvals/bvals.cpp"
BOUNDARY_FILES = (
    (
        REPO_ROOT / "src/bvals/bvals_fc.cpp",
        "TaskStatus MeshBoundaryValuesFC::RecvAndUnpackFC",
    ),
    (
        REPO_ROOT / "src/bvals/bvals_cc.cpp",
        "TaskStatus MeshBoundaryValuesCC::RecvAndUnpackCC",
    ),
)


def function_body(path: Path, signature: str) -> str:
    source = path.read_text(encoding="utf-8")
    start = source.index(signature)
    opening_brace = source.index("{", start)
    depth = 0
    for offset, character in enumerate(source[opening_brace:], opening_brace):
        if character == "{":
            depth += 1
        elif character == "}":
            depth -= 1
            if depth == 0:
                return source[opening_brace : offset + 1]
    raise AssertionError(f"unterminated function in {path}")


class KokkosCompletionFenceTests(unittest.TestCase):
    def test_mpi_teardown_precedes_application_and_kokkos_teardown(self) -> None:
        source = MAIN.read_text(encoding="utf-8")
        self.assertIsNotNone(
            re.search(
                r"Kokkos::fence\(\);\s*"
                r"#if MPI_PARALLEL_ENABLED\s*"
                r"MPI_Finalize\(\);\s*"
                r"#endif\s*"
                r"delete pout;\s*"
                r"delete pdriver;\s*"
                r"delete pmesh;\s*"
                r"delete pinput;\s*"
                r"Kokkos::finalize\(\);",
                source,
            )
        )

    def test_mpi_gpu_transport_is_finalized_before_kokkos(self) -> None:
        source = MAIN.read_text(encoding="utf-8")
        helper = function_body(MAIN, "void FinalizeParallelRuntime()")
        self.assertLess(helper.index("MPI_Finalize();"), helper.index("Kokkos::finalize();"))
        self.assertNotRegex(
            source,
            r"Kokkos::finalize\(\);\s*(?:#if MPI_PARALLEL_ENABLED\s*)?"
            r"MPI_Finalize\(\);",
        )

    def test_mpi_communicators_are_not_freed_after_mpi_finalize(self) -> None:
        helper = function_body(BVALS, "void FreeCommunicatorIfMPIActive")
        finalized_check = helper.index("MPI_Finalized(&finalized)")
        communicator_free = helper.index("MPI_Comm_free(comm);")
        self.assertLess(finalized_check, communicator_free)
        self.assertIn("finalized != 0", helper)

        source = BVALS.read_text(encoding="utf-8")
        self.assertEqual(source.count("FreeCommunicatorIfMPIActive(&comm_"), 2)
        self.assertEqual(source.count("FreeCommunicatorIfMPIActive(&mpi_comm_part)"), 1)
        self.assertNotRegex(source, r"MPI_Comm_free\(&(?:comm_vars|comm_flux|mpi_comm_part)\)")

    def test_receive_unpack_finishes_before_task_reports_completion(self) -> None:
        for path, signature in BOUNDARY_FILES:
            with self.subTest(path=path.name):
                body = function_body(path, signature)
                final_unpack = body.rfind('Kokkos::parallel_for("RecvBuff"')
                fence = body.rfind("Kokkos::fence();")
                completion = body.rfind("return TaskStatus::complete;")
                self.assertGreater(final_unpack, 0)
                self.assertGreater(fence, final_unpack)
                self.assertGreater(completion, fence)


if __name__ == "__main__":
    unittest.main()
