#!/usr/bin/env python3
"""Prepare the non-authorizing deterministic Q009 execution matrix.

The checked-in Q009 contract fixes five logical cases and seven scheduler jobs.
This module deterministically prepares aggregate-only decks, launch-contract
syntax and policy review candidates, and explicit launch-prohibited
restart-continuation handoffs in memory.  Its minimal structured syntax
validator and in-process integrity guards detect ordinary accidental drift; they
are not an authority boundary against privileged code in the same CPython
process.  No persistent-publication or cleanup surface is present.  It never
accepts operational bindings, mutates live policy, calls a scheduler, reads
runtime output, produces receipts, performs admission, or grants authority.

Every returned object is a non-authoritative source-local review input,
regardless of its contents.  Any future launch, admission, evidence, or
publication authority must be established by a separately reviewed external
process and receipt chain.
"""

_BUILTIN_FUNCTION_TYPE = [].append.__class__
_FUNCTION_TYPE = (lambda: None).__class__
_TYPE_TYPE = (0).__class__.__class__
_MAPPING_PROXY_TYPE = _TYPE_TYPE(_TYPE_TYPE.__dict__)
_PY_TPFLAGS_HEAPTYPE = 1 << 9


def _is_exact_builtin_function(value, name, module="builtins"):
    return (
        _TYPE_TYPE(value) is _BUILTIN_FUNCTION_TYPE
        and value.__module__ == module
        and value.__name__ == name
    )


def _is_exact_static_builtin_type(value, name, module="builtins"):
    return (
        _TYPE_TYPE(value) is _TYPE_TYPE
        and value.__module__ == module
        and value.__name__ == name
        and value.__flags__ & _PY_TPFLAGS_HEAPTYPE == 0
    )


if not (
    bool is True.__class__
    and bytes is b"".__class__
    and dict is {}.__class__
    and float is (0.0).__class__
    and int is (0).__class__
    and list is [].__class__
    and object is _FUNCTION_TYPE.__base__
    and set is {0}.__class__
    and str is "".__class__
    and tuple is ().__class__
    and type is _TYPE_TYPE
    and _is_exact_static_builtin_type(enumerate, "enumerate")
    and _is_exact_static_builtin_type(range, "range")
    and _is_exact_static_builtin_type(zip, "zip")
    and _is_exact_static_builtin_type(ValueError, "ValueError")
    and _is_exact_static_builtin_type(OSError, "OSError")
    and _is_exact_static_builtin_type(UnicodeDecodeError, "UnicodeDecodeError")
):
    raise 0

for _candidate, _name in (
    (all, "all"),
    (any, "any"),
    (__build_class__, "__build_class__"),
    (callable, "callable"),
    (chr, "chr"),
    (compile, "compile"),
    (exec, "exec"),
    (format, "format"),
    (getattr, "getattr"),
    (isinstance, "isinstance"),
    (len, "len"),
    (open, "open"),
    (ord, "ord"),
    (sorted, "sorted"),
    (sum, "sum"),
):
    if not _is_exact_builtin_function(
        _candidate, _name, "io" if _name == "open" else "builtins"
    ):
        raise 0

try:
    float(b"Q009-invalid-float")
except ValueError:
    pass
else:
    raise 0

_TRUSTED_ALL = all
_TRUSTED_ANY = any
_TRUSTED_BUILD_CLASS = __build_class__
_TRUSTED_BOOL = bool
_TRUSTED_BYTES = bytes
_TRUSTED_CALLABLE = callable
_TRUSTED_CHR = chr
_TRUSTED_COMPILE = compile
_TRUSTED_DICT = dict
_TRUSTED_ENUMERATE = enumerate
_TRUSTED_EXEC = exec
_TRUSTED_FLOAT = float
_TRUSTED_FORMAT = format
_TRUSTED_GETATTR = getattr
_TRUSTED_INT = int
_TRUSTED_ISINSTANCE = isinstance
_TRUSTED_LEN = len
_TRUSTED_LIST = list
_TRUSTED_OS_ERROR = OSError
_TRUSTED_ORD = ord
_TRUSTED_RANGE = range
_TRUSTED_SET = set
_TRUSTED_SORTED = sorted
_TRUSTED_STR = str
_TRUSTED_SUM = sum
_TRUSTED_TYPE = type
_TRUSTED_UNICODE_DECODE_ERROR = UnicodeDecodeError
_TRUSTED_VALUE_ERROR = ValueError
_TRUSTED_ZIP = zip

Callable = list
Mapping = dict
Sequence = list

if type(__file__) is not str:
    raise 0

REPO_ROOT = __file__.rsplit("/", 3)[0]
AUTHORIZED_ORION_ROOT = "/lustre/orion/ast207/proj-shared/dfielding/PIC"
CANONICAL_PROJECT_HOME_ROOT = "/autofs/nccs-svm1_proj/ast207/proj-shared/PIC"
LEGACY_PROJECT_HOME_ROOT = "/ccs/proj/ast207/proj-shared/PIC"
CONTRACT_PATH = "tst/publication/q009_dynamic_amr_load_balance_registered_pilot_contract_v1.json"
CONTRACT_SHA256 = "2567bd9490e1e643c37cc63973307eead271a781b6ebbe8abc60edfefbb71450"
CONTROL_PLANE_LAUNCH_SYNTAX_REFERENCE_PATH = (
    "tst/publication/frontier_control_plane/control_plane_common.py"
)
CONTRACT_VALIDATOR_PATH = (
    "tst/publication/q009_dynamic_amr_load_balance_registered_pilot_contract_v1.py"
)
MINIMAL_LAUNCH_SYNTAX_VALIDATOR_ID = (
    "_q009_preparation_minimal_structured_launch_syntax"
)
CAMPAIGN = "q009_dynamic_amr_load_balance_registered_pilot_v1"
SUCCESSOR_ID = "q009_deterministic_execution_preparation_only_successor_v2"
SCHEMA_VERSION = 2
RECORD_TYPE = "q009_deterministic_execution_preparation_only_manifest"
LAUNCH_RECORD_TYPE = "q009_dynamic_amr_launch_review_candidate"
POLICY_RECORD_TYPE = "q009_dynamic_amr_policy_slice_review_candidate"
HANDOFF_RECORD_TYPE = "q009_dynamic_amr_restart_continuation_launch_prohibited_handoff"
QUALIFICATION_EFFECT = (
    "source_local_execution_preparation_only_no_launch_no_policy_no_q009_"
    "no_dynamic_amr_science_no_receipt_no_admission_no_q011_no_publication_authority"
)
AUTHORITY_SEMANTICS = (
    "none_source_local_review_input_only_never_launch_admission_evidence_or_"
    "publication_authority_regardless_of_content"
)
EVIDENCE_CLASS = "q009_dynamic_amr_load_balance_registered_engineering_pilot"
PHYSICAL_MODE = "q009_dynamic_amr_load_balance_aggregate_only"
RUNTIME_PROFILE = "frontier_minimum_supported"
SUBMISSION_ID_TEMPLATE = "{submission_id}"
MAXIMUM_ATTEMPTS_PER_CASE = 1
MAXIMUM_RETRIES_PER_CASE = 0
MAXIMUM_LIVE_Q009_SUBMISSIONS = 1

AUTHORIZATION_BOUNDARY = _MAPPING_PROXY_TYPE(
    {
        "launch_authorized": False,
        "scheduler_submission_authorized": False,
        "policy_mutation_authorized": False,
        "frontier_execution_authorized": False,
        "receipt_production_authorized": False,
        "admission_authorized": False,
        "operational_acceptance_authorized": False,
        "runtime_telemetry_claim_authorized": False,
        "q009_qualified": False,
        "dynamic_amr_science_qualified": False,
        "exact_amr_conservation_claim_authorized": False,
        "q011_production_authorized": False,
        "scientific_claim_authorized": False,
        "publication_authorized": False,
    }
)
THREAT_MODEL = _MAPPING_PROXY_TYPE(
    {
        "contract": "trusted_cpython_process_source_local_preparation_only",
        "privileged_same_process_python_trusted": True,
        "privileged_same_process_python_resisted": False,
        "closure_or_code_object_mutation_resisted": False,
        "mapping_proxy_backing_dict_recovery_resisted": False,
        "os_or_process_isolation_present": False,
        "external_sha256_receipt_verifier_present": False,
        "in_process_integrity_guards_role": "ordinary_accidental_drift_detection_only",
        "all_public_outputs_authoritative": False,
        "all_public_outputs_can_grant_authority": False,
    }
)
ABSENT_CAPABILITIES = {
    "operational_binding_input_implemented": False,
    "exact_launch_candidate_validator_implemented": False,
    "runtime_telemetry_implemented": False,
    "pairwise_lb_on_off_acceptance_implemented": False,
    "trusted_unique_rst_restart_launching_implemented": False,
    "hardened_external_completion_root_inventory_receipts_implemented": False,
    "receipt_producer_validator_implemented": False,
    "admission_bridge_implemented": False,
    "persistent_review_bundle_publication_implemented": False,
}
UNRESOLVED_OPERATIONAL_BINDINGS = (
    "final_clean_candidate_manifest_sha256",
    "final_executable_sha256",
    "installed_environment_profile_sha256",
    "installed_job_script_sha256",
    "runtime_telemetry_analysis_script_sha256",
    "pairwise_acceptance_analysis_script_sha256",
    "separately_promoted_registered_science_authorization_id",
)

COMMON_RUNTIME_BLOCKERS = (
    "final_clean_candidate_source_archive_executable_and_environment_bound",
    "paired_installed_control_plane_generation_bound_and_independently_verified",
    "separate_reviewed_registered_science_policy_slice_promoted_for_exact_stage",
    "project_wide_node_hour_ledger_other_campaign_reserves_and_storage_bound",
    "fresh_empty_user_queue_proved_immediately_before_each_serial_submission",
    "one_logical_attempt_per_case_and_zero_retries_enforced",
    "contract_complete_event_per_rank_per_level_runtime_telemetry_implemented_and_bound",
    "pairwise_lb_on_off_operational_acceptance_implemented_and_bound",
    "hardened_external_completion_root_and_inventory_receipts_installed_and_bound",
    "forbidden_physical_field_and_particle_inspection_not_performed",
)
RESTART_CONTINUATION_BLOCKER = (
    "trusted_unique_rst_restart_snapshot_launching_not_implemented_or_bound"
)
MIGRATION_DESIGN_BLOCKER = (
    "migration_cases_require_redesign_to_produce_contract_required_repeated_"
    "refine_derefine_and_inter_rank_migration_events"
)
PERSISTENT_REVIEW_BUNDLE_BLOCKER = (
    "persistent_review_bundle_publication_requires_kernel_enforced_immutability_"
    "not_provided_by_the_same_principal_mutable_source_tree"
)


PreparationError = ValueError


def _make_runtime_guard():
    trusted_all = _TRUSTED_ALL
    trusted_any = _TRUSTED_ANY
    trusted_bool = _TRUSTED_BOOL
    trusted_bytes = _TRUSTED_BYTES
    trusted_callable = _TRUSTED_CALLABLE
    trusted_chr = _TRUSTED_CHR
    trusted_dict = _TRUSTED_DICT
    trusted_enumerate = _TRUSTED_ENUMERATE
    trusted_float = _TRUSTED_FLOAT
    trusted_format = _TRUSTED_FORMAT
    trusted_getattr = _TRUSTED_GETATTR
    trusted_int = _TRUSTED_INT
    trusted_isinstance = _TRUSTED_ISINSTANCE
    trusted_len = _TRUSTED_LEN
    trusted_list = _TRUSTED_LIST
    trusted_object = object
    trusted_ord = _TRUSTED_ORD
    trusted_os_error = _TRUSTED_OS_ERROR
    trusted_range = _TRUSTED_RANGE
    trusted_set = _TRUSTED_SET
    trusted_sorted = _TRUSTED_SORTED
    trusted_str = _TRUSTED_STR
    trusted_sum = _TRUSTED_SUM
    trusted_tuple = tuple
    trusted_type = _TRUSTED_TYPE
    trusted_unicode_decode_error = _TRUSTED_UNICODE_DECODE_ERROR
    trusted_value_error = _TRUSTED_VALUE_ERROR
    trusted_zip = _TRUSTED_ZIP

    def require_trusted_runtime() -> None:
        if not (
            all is trusted_all
            and any is trusted_any
            and bool is trusted_bool
            and bytes is trusted_bytes
            and callable is trusted_callable
            and chr is trusted_chr
            and dict is trusted_dict
            and enumerate is trusted_enumerate
            and float is trusted_float
            and format is trusted_format
            and getattr is trusted_getattr
            and int is trusted_int
            and isinstance is trusted_isinstance
            and len is trusted_len
            and list is trusted_list
            and object is trusted_object
            and ord is trusted_ord
            and OSError is trusted_os_error
            and range is trusted_range
            and set is trusted_set
            and sorted is trusted_sorted
            and str is trusted_str
            and sum is trusted_sum
            and tuple is trusted_tuple
            and type is trusted_type
            and UnicodeDecodeError is trusted_unicode_decode_error
            and ValueError is trusted_value_error
            and zip is trusted_zip
        ):
            raise trusted_value_error(
                "Q009 preparation rejected poisoned runtime primitives"
            )

    return require_trusted_runtime


_require_trusted_runtime = _make_runtime_guard()
del _make_runtime_guard


def _require(condition: bool, message: str) -> None:
    _require_trusted_runtime()
    if type(condition) is not bool or type(message) is not str:
        raise PreparationError("Q009 preparation condition is not exact")
    if not condition:
        raise PreparationError(message)


def _strict_equal(left: object, right: object) -> bool:
    _require_trusted_runtime()
    if not _is_exact_json_graph(left) or not _is_exact_json_graph(right):
        return False
    pending = [(left, right)]
    while pending:
        left_member, right_member = pending.pop()
        kind = type(left_member)
        if kind is not type(right_member):
            return False
        if kind is dict:
            if len(left_member) != len(right_member):
                return False
            for key, value in left_member.items():
                if key not in right_member:
                    return False
                pending.append((value, right_member[key]))
            continue
        if kind is list:
            if len(left_member) != len(right_member):
                return False
            for index in range(len(left_member)):
                pending.append((left_member[index], right_member[index]))
            continue
        if left_member != right_member:
            return False
    return True


def _is_exact_json_graph(value: object) -> bool:
    _require_trusted_runtime()
    pending = [value]
    containers = []
    while pending:
        member = pending.pop()
        kind = type(member)
        if member is None or kind is bool or kind is int or kind is str:
            continue
        if kind is float:
            if not (
                member == member
                and -1.7976931348623157e308 <= member <= 1.7976931348623157e308
            ):
                return False
            continue
        if kind is list:
            for prior in containers:
                if member is prior:
                    return False
            containers.append(member)
            pending.extend(member)
            continue
        if kind is dict:
            for prior in containers:
                if member is prior:
                    return False
            containers.append(member)
            for key, value in member.items():
                if type(key) is not str:
                    return False
                pending.append(value)
            continue
        return False
    return True


def _is_exact_file_graph(value: object) -> bool:
    _require_trusted_runtime()
    if type(value) is not dict:
        return False
    for path, payload in value.items():
        if type(path) is not str or type(payload) is not bytes:
            return False
    return True


_JSON_HEX = "0123456789abcdef"
_JSON_ESCAPES = {
    '"': '\\"',
    "\\": "\\\\",
    "\b": "\\b",
    "\f": "\\f",
    "\n": "\\n",
    "\r": "\\r",
    "\t": "\\t",
}


def _json_hex4(value: int) -> str:
    return "".join(_JSON_HEX[(value >> shift) & 15] for shift in (12, 8, 4, 0))


def _json_string(value: str) -> str:
    encoded = ['"']
    for character in value:
        if character in _JSON_ESCAPES:
            encoded.append(_JSON_ESCAPES[character])
            continue
        codepoint = ord(character)
        if 0x20 <= codepoint <= 0x7E:
            encoded.append(character)
        elif codepoint <= 0xFFFF:
            encoded.append("\\u" + _json_hex4(codepoint))
        else:
            adjusted = codepoint - 0x10000
            encoded.append("\\u" + _json_hex4(0xD800 + (adjusted >> 10)))
            encoded.append("\\u" + _json_hex4(0xDC00 + (adjusted & 0x3FF)))
    encoded.append('"')
    return "".join(encoded)


def _canonical_json(value: object) -> str:
    kind = type(value)
    if value is None:
        return "null"
    if value is True:
        return "true"
    if value is False:
        return "false"
    if kind is int:
        return str(value)
    if kind is float:
        _require(
            value == value
            and -1.7976931348623157e308 <= value <= 1.7976931348623157e308,
            "nonfinite JSON number",
        )
        rendered = format(value, ".17g")
        if "." not in rendered and "e" not in rendered:
            rendered += ".0"
        return rendered
    if kind is str:
        return _json_string(value)
    if kind is list:
        return "[" + ",".join(_canonical_json(item) for item in value) + "]"
    if kind is dict:
        _require(
            all(type(key) is str for key in value),
            "Q009 preparation JSON object key is not a string",
        )
        return "{" + ",".join(
            _json_string(key) + ":" + _canonical_json(value[key])
            for key in sorted(value)
        ) + "}"
    raise PreparationError("Q009 preparation contains noncanonical JSON")


def _json_bytes(value: object) -> bytes:
    _require_trusted_runtime()
    _require(_is_exact_json_graph(value), "Q009 preparation is not an exact JSON graph")
    return (_canonical_json(value) + "\n").encode("utf-8")


_SHA256_INITIAL = (
    0x6A09E667,
    0xBB67AE85,
    0x3C6EF372,
    0xA54FF53A,
    0x510E527F,
    0x9B05688C,
    0x1F83D9AB,
    0x5BE0CD19,
)
_SHA256_ROUND = (
    0x428A2F98, 0x71374491, 0xB5C0FBCF, 0xE9B5DBA5,
    0x3956C25B, 0x59F111F1, 0x923F82A4, 0xAB1C5ED5,
    0xD807AA98, 0x12835B01, 0x243185BE, 0x550C7DC3,
    0x72BE5D74, 0x80DEB1FE, 0x9BDC06A7, 0xC19BF174,
    0xE49B69C1, 0xEFBE4786, 0x0FC19DC6, 0x240CA1CC,
    0x2DE92C6F, 0x4A7484AA, 0x5CB0A9DC, 0x76F988DA,
    0x983E5152, 0xA831C66D, 0xB00327C8, 0xBF597FC7,
    0xC6E00BF3, 0xD5A79147, 0x06CA6351, 0x14292967,
    0x27B70A85, 0x2E1B2138, 0x4D2C6DFC, 0x53380D13,
    0x650A7354, 0x766A0ABB, 0x81C2C92E, 0x92722C85,
    0xA2BFE8A1, 0xA81A664B, 0xC24B8B70, 0xC76C51A3,
    0xD192E819, 0xD6990624, 0xF40E3585, 0x106AA070,
    0x19A4C116, 0x1E376C08, 0x2748774C, 0x34B0BCB5,
    0x391C0CB3, 0x4ED8AA4A, 0x5B9CCA4F, 0x682E6FF3,
    0x748F82EE, 0x78A5636F, 0x84C87814, 0x8CC70208,
    0x90BEFFFA, 0xA4506CEB, 0xBEF9A3F7, 0xC67178F2,
)


def _rotate_right(value: int, count: int) -> int:
    return ((value >> count) | (value << (32 - count))) & 0xFFFFFFFF


def _sha256_bytes(payload: bytes) -> str:
    """Return SHA-256 without trusting an importable hashing module."""
    _require_trusted_runtime()
    _require(type(payload) is bytes, "SHA-256 input must be exact bytes")
    message = list(payload)
    bit_length = len(message) * 8
    message.append(0x80)
    while len(message) % 64 != 56:
        message.append(0)
    message.extend(bit_length.to_bytes(8, "big"))
    state = list(_SHA256_INITIAL)
    for offset in range(0, len(message), 64):
        words = [
            int.from_bytes(message[index : index + 4], "big")
            for index in range(offset, offset + 64, 4)
        ]
        for index in range(16, 64):
            sigma0 = (
                _rotate_right(words[index - 15], 7)
                ^ _rotate_right(words[index - 15], 18)
                ^ (words[index - 15] >> 3)
            )
            sigma1 = (
                _rotate_right(words[index - 2], 17)
                ^ _rotate_right(words[index - 2], 19)
                ^ (words[index - 2] >> 10)
            )
            words.append(
                (words[index - 16] + sigma0 + words[index - 7] + sigma1)
                & 0xFFFFFFFF
            )
        a, b, c, d, e, f, g, h = state
        for index in range(64):
            sum1 = _rotate_right(e, 6) ^ _rotate_right(e, 11) ^ _rotate_right(e, 25)
            choose = (e & f) ^ ((~e) & g)
            temporary1 = (h + sum1 + choose + _SHA256_ROUND[index] + words[index]) & 0xFFFFFFFF
            sum0 = _rotate_right(a, 2) ^ _rotate_right(a, 13) ^ _rotate_right(a, 22)
            majority = (a & b) ^ (a & c) ^ (b & c)
            temporary2 = (sum0 + majority) & 0xFFFFFFFF
            h, g, f, e, d, c, b, a = (
                g,
                f,
                e,
                (d + temporary1) & 0xFFFFFFFF,
                c,
                b,
                a,
                (temporary1 + temporary2) & 0xFFFFFFFF,
            )
        state = [
            (value + addition) & 0xFFFFFFFF
            for value, addition in zip(state, (a, b, c, d, e, f, g, h))
        ]
    return "".join(f"{value:08x}" for value in state)


def _binding(path: str, payload: bytes) -> dict[str, object]:
    return {"path": path, "sha256": _sha256_bytes(payload), "byte_count": len(payload)}


def _stable_repository_bytes(relative: str, *, label: str) -> bytes:
    """Capture descriptor-bound repository bytes using no importable module."""
    _require_trusted_runtime()
    _require(
        type(relative) is str and type(label) is str,
        "repository capture requires exact string inputs",
    )
    parts = relative.split("/")
    _require(
        relative
        and not relative.startswith("/")
        and "/".join(parts) == relative
        and all(part not in {"", ".", ".."} for part in parts),
        f"{label}: unsafe repository-relative path",
    )
    reader = open
    _require(
        _is_exact_builtin_function(reader, "open", "io"),
        f"{label}: rejected poisoned read primitive",
    )
    try:
        with reader(REPO_ROOT + "/" + relative, "rb") as source:
            payload = source.read()
    except OSError as error:
        raise PreparationError(f"{label}: unavailable below repository root") from error
    _require(payload.__class__ is bytes, f"{label}: did not produce exact bytes")
    return payload


def _captured_source_binding(
    relative: str, payload: bytes, *, loading_role: str
) -> dict[str, object]:
    return {
        "path": relative,
        "sha256": _sha256_bytes(payload),
        "byte_count": len(payload),
        "loading_role": loading_role,
    }


def _decode_json(payload: bytes, *, label: str) -> object:
    """Decode exact JSON without trusting an importable parser."""
    _require_trusted_runtime()
    _require(
        type(payload) is bytes and type(label) is str,
        "JSON decoding requires exact byte and string inputs",
    )
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise PreparationError(label + ": JSON is not UTF-8") from error
    index = 0

    def skip_space() -> None:
        nonlocal index
        while index < len(text) and text[index] in " \t\r\n":
            index += 1

    def parse_string() -> str:
        nonlocal index
        _require(index < len(text) and text[index] == '"', label + ": expected JSON string")
        index += 1
        result: list[str] = []
        escapes = {
            '"': '"',
            "\\": "\\",
            "/": "/",
            "b": "\b",
            "f": "\f",
            "n": "\n",
            "r": "\r",
            "t": "\t",
        }
        while index < len(text):
            character = text[index]
            index += 1
            if character == '"':
                return "".join(result)
            _require(ord(character) >= 0x20, label + ": control character in JSON string")
            if character != "\\":
                result.append(character)
                continue
            _require(index < len(text), label + ": truncated JSON escape")
            escape = text[index]
            index += 1
            if escape in escapes:
                result.append(escapes[escape])
                continue
            _require(
                escape == "u" and index + 4 <= len(text),
                label + ": malformed JSON escape",
            )
            digits = text[index : index + 4]
            _require(
                all(character in _JSON_HEX or character in "ABCDEF" for character in digits),
                label + ": malformed JSON Unicode escape",
            )
            index += 4
            result.append(chr(int(digits, 16)))
        raise PreparationError(label + ": unterminated JSON string")

    def parse_number() -> object:
        nonlocal index
        start = index
        while index < len(text) and text[index] in "-+0123456789.eE":
            index += 1
        token = text[start:index]
        _require(token and token not in {"-", "+"}, label + ": malformed JSON number")
        try:
            if "." in token or "e" in token or "E" in token:
                value = float(token)
                _require(
                    value == value and -1.7976931348623157e308 <= value <= 1.7976931348623157e308,
                    label + ": nonfinite JSON number",
                )
                return value
            return int(token)
        except ValueError as error:
            raise PreparationError(label + ": malformed JSON number") from error

    def parse_value() -> object:
        nonlocal index
        skip_space()
        _require(index < len(text), label + ": truncated JSON")
        character = text[index]
        if character == '"':
            return parse_string()
        if character == "{":
            index += 1
            result: dict[str, object] = {}
            skip_space()
            if index < len(text) and text[index] == "}":
                index += 1
                return result
            while True:
                skip_space()
                key = parse_string()
                _require(key not in result, label + ": duplicate JSON object key")
                skip_space()
                _require(index < len(text) and text[index] == ":", label + ": expected JSON colon")
                index += 1
                result[key] = parse_value()
                skip_space()
                _require(index < len(text), label + ": truncated JSON object")
                separator = text[index]
                index += 1
                if separator == "}":
                    return result
                _require(separator == ",", label + ": expected JSON object separator")
        if character == "[":
            index += 1
            result_list: list[object] = []
            skip_space()
            if index < len(text) and text[index] == "]":
                index += 1
                return result_list
            while True:
                result_list.append(parse_value())
                skip_space()
                _require(index < len(text), label + ": truncated JSON array")
                separator = text[index]
                index += 1
                if separator == "]":
                    return result_list
                _require(separator == ",", label + ": expected JSON array separator")
        for literal, value in (("true", True), ("false", False), ("null", None)):
            if text.startswith(literal, index):
                index += len(literal)
                return value
        return parse_number()

    value = parse_value()
    skip_space()
    _require(index == len(text), label + ": trailing JSON bytes")
    return value


def _load_frozen_contract() -> dict[str, object]:
    payload = _stable_repository_bytes(CONTRACT_PATH, label="Q009 exact contract JSON")
    _require(_sha256_bytes(payload) == CONTRACT_SHA256, "Q009 exact contract JSON digest drifted")
    value = _decode_json(payload, label="Q009 exact contract JSON")
    _require(
        value.__class__ is dict
        and value.get("contract_id") == "q009_dynamic_amr_load_balance_registered_pilot_contract_v1"
        and value.get("status") == "preregistered_non_authorizing_execution_blocked"
        and value.get("authority_boundary").__class__ is dict
        and all(member is False for member in value["authority_boundary"].values()),
        "Q009 exact contract JSON boundary drifted",
    )
    return value


_MINIMAL_LAUNCH_SYNTAX_VALIDATOR_SOURCE = r'''
TRUSTED_LAUNCH_EXECUTOR = "trusted_trampoline_athena_argv_v1"
_LOWER_ALNUM = "abcdefghijklmnopqrstuvwxyz0123456789"
_LOWER_ID_REST = _LOWER_ALNUM + "_-"
_ASCII_ALNUM_UNDERSCORE = (
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_"
)
_HEX = "0123456789abcdef"
_JSON_ESCAPES = (
    ('"', '\\"'),
    ("\\", "\\\\"),
    ("\b", "\\b"),
    ("\f", "\\f"),
    ("\n", "\\n"),
    ("\r", "\\r"),
    ("\t", "\\t"),
)
_SHA256_INITIAL = (
    0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
    0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19,
)
_SHA256_ROUND = (
    0x428A2F98, 0x71374491, 0xB5C0FBCF, 0xE9B5DBA5,
    0x3956C25B, 0x59F111F1, 0x923F82A4, 0xAB1C5ED5,
    0xD807AA98, 0x12835B01, 0x243185BE, 0x550C7DC3,
    0x72BE5D74, 0x80DEB1FE, 0x9BDC06A7, 0xC19BF174,
    0xE49B69C1, 0xEFBE4786, 0x0FC19DC6, 0x240CA1CC,
    0x2DE92C6F, 0x4A7484AA, 0x5CB0A9DC, 0x76F988DA,
    0x983E5152, 0xA831C66D, 0xB00327C8, 0xBF597FC7,
    0xC6E00BF3, 0xD5A79147, 0x06CA6351, 0x14292967,
    0x27B70A85, 0x2E1B2138, 0x4D2C6DFC, 0x53380D13,
    0x650A7354, 0x766A0ABB, 0x81C2C92E, 0x92722C85,
    0xA2BFE8A1, 0xA81A664B, 0xC24B8B70, 0xC76C51A3,
    0xD192E819, 0xD6990624, 0xF40E3585, 0x106AA070,
    0x19A4C116, 0x1E376C08, 0x2748774C, 0x34B0BCB5,
    0x391C0CB3, 0x4ED8AA4A, 0x5B9CCA4F, 0x682E6FF3,
    0x748F82EE, 0x78A5636F, 0x84C87814, 0x8CC70208,
    0x90BEFFFA, 0xA4506CEB, 0xBEF9A3F7, 0xC67178F2,
)


def _rotate_right(value, count):
    return ((value >> count) | (value << (32 - count))) & 0xFFFFFFFF


def _sha256_hex(payload):
    message = list(payload)
    bit_length = len(message) * 8
    message.append(0x80)
    while len(message) % 64 != 56:
        message.append(0)
    message.extend(bit_length.to_bytes(8, "big"))
    state = list(_SHA256_INITIAL)
    for offset in range(0, len(message), 64):
        words = [
            int.from_bytes(message[index:index + 4], "big")
            for index in range(offset, offset + 64, 4)
        ]
        for index in range(16, 64):
            sigma0 = (
                _rotate_right(words[index - 15], 7)
                ^ _rotate_right(words[index - 15], 18)
                ^ (words[index - 15] >> 3)
            )
            sigma1 = (
                _rotate_right(words[index - 2], 17)
                ^ _rotate_right(words[index - 2], 19)
                ^ (words[index - 2] >> 10)
            )
            words.append(
                (words[index - 16] + sigma0 + words[index - 7] + sigma1)
                & 0xFFFFFFFF
            )
        a, b, c, d, e, f, g, h = state
        for index in range(64):
            sum1 = _rotate_right(e, 6) ^ _rotate_right(e, 11) ^ _rotate_right(e, 25)
            choose = (e & f) ^ ((~e) & g)
            temporary1 = (
                h + sum1 + choose + _SHA256_ROUND[index] + words[index]
            ) & 0xFFFFFFFF
            sum0 = _rotate_right(a, 2) ^ _rotate_right(a, 13) ^ _rotate_right(a, 22)
            majority = (a & b) ^ (a & c) ^ (b & c)
            temporary2 = (sum0 + majority) & 0xFFFFFFFF
            h, g, f, e, d, c, b, a = (
                g,
                f,
                e,
                (d + temporary1) & 0xFFFFFFFF,
                c,
                b,
                a,
                (temporary1 + temporary2) & 0xFFFFFFFF,
            )
        state = [
            (value + addition) & 0xFFFFFFFF
            for value, addition in zip(state, (a, b, c, d, e, f, g, h))
        ]
    return "".join(f"{value:08x}" for value in state)


def _lower_identifier(value):
    return (
        type(value) is str
        and 1 <= len(value) <= 64
        and value[0] in _LOWER_ALNUM
        and all(character in _LOWER_ID_REST for character in value[1:])
    )


def _analysis_role(value):
    prefix = "analysis-script-"
    suffix = value[len(prefix):] if type(value) is str and value.startswith(prefix) else ""
    return len(suffix) == 3 and all(character in "0123456789" for character in suffix)


def _relative_artifact_path(value, field):
    if type(value) is not str:
        raise ValueError(field + " must be an exact JSON string")
    text = value
    parts = text.split("/")
    if (
        not text
        or text.startswith("/")
        or any(part in {"", ".", ".."} for part in parts)
    ):
        raise ValueError(field + " must be a non-empty relative artifact path")
    return text


def _override_key(value):
    parts = value.split("/")
    return (
        len(parts) >= 2
        and all(
            part
            and all(character in _ASCII_ALNUM_UNDERSCORE for character in part)
            for part in parts
        )
    )


def _hex4(value):
    return "".join(_HEX[(value >> shift) & 15] for shift in [12, 8, 4, 0])


def _json_string(value):
    encoded = ['"']
    for character in value:
        escaped = None
        for source, replacement in _JSON_ESCAPES:
            if character == source:
                escaped = replacement
                break
        if escaped is not None:
            encoded.append(escaped)
            continue
        codepoint = ord(character)
        if 0x20 <= codepoint <= 0x7e:
            encoded.append(character)
        elif codepoint <= 0xffff:
            encoded.append("\\u" + _hex4(codepoint))
        else:
            adjusted = codepoint - 0x10000
            encoded.append("\\u" + _hex4(0xd800 + (adjusted >> 10)))
            encoded.append("\\u" + _hex4(0xdc00 + (adjusted & 0x3ff)))
    encoded.append('"')
    return "".join(encoded)


def _canonical_json(value):
    if value is None:
        return "null"
    if value is True:
        return "true"
    if value is False:
        return "false"
    if type(value) is int:
        return str(value)
    if type(value) is str:
        return _json_string(value)
    if type(value) is list:
        return "[" + ",".join(_canonical_json(item) for item in value) + "]"
    if type(value) is dict and all(type(key) is str for key in value):
        return "{" + ",".join(
            _json_string(key) + ":" + _canonical_json(value[key])
            for key in sorted(value)
        ) + "}"
    raise ValueError("Launch contract contains a noncanonical JSON value")


def validate_launch_contract(value):
    if type(value) is not dict or set(value) != {
        "schema_version",
        "executor",
        "pre_actions",
        "actions",
        "post_actions",
    }:
        raise ValueError(
            "Launch contract must contain only schema_version, executor, "
            "pre_actions, actions and post_actions"
        )
    if type(value.get("schema_version")) is not int or value.get("schema_version") != 1:
        raise ValueError("Unsupported launch-contract schema")
    if value.get("executor") != TRUSTED_LAUNCH_EXECUTOR:
        raise ValueError("Launch contract does not select the trusted Athena executor")
    actions = value.get("actions")
    if type(actions) is not list or not 1 <= len(actions) <= 16:
        raise ValueError("Launch contract requires between one and sixteen Athena actions")
    identifiers = set()
    for phase in ["pre_actions", "post_actions"]:
        bounded_actions = value.get(phase)
        if type(bounded_actions) is not list or len(bounded_actions) > 16:
            raise ValueError("Launch contract " + phase + " must contain at most sixteen actions")
        for action in bounded_actions:
            if type(action) is not dict:
                raise ValueError("Launch-contract " + phase + " action is malformed")
            identifier = action.get("action_id")
            if not _lower_identifier(identifier):
                raise ValueError("Launch-contract " + phase + " action ID is malformed")
            if identifier in identifiers:
                raise ValueError("Duplicate launch action ID: " + identifier)
            identifiers.add(identifier)
            kind = action.get("kind")
            if kind == "snapshot_sha256":
                if set(action) != {
                    "action_id",
                    "kind",
                    "snapshot_role",
                    "output_artifact",
                }:
                    raise ValueError("Launch-contract " + phase + " snapshot action is malformed")
                role = action.get("snapshot_role")
                if role not in {
                    "job-script",
                    "executable",
                    "input-deck",
                    "environment-profile",
                    "timeout-margin",
                    "queue-snapshot",
                } and not _analysis_role(role):
                    raise ValueError("Launch-contract " + phase + " snapshot role is malformed")
                _relative_artifact_path(action.get("output_artifact"), "output_artifact")
            elif kind == "artifact_sha256":
                if set(action) != {
                    "action_id",
                    "kind",
                    "artifact",
                    "output_artifact",
                }:
                    raise ValueError("Launch-contract " + phase + " artifact action is malformed")
                artifact = _relative_artifact_path(action.get("artifact"), "artifact")
                output = _relative_artifact_path(
                    action.get("output_artifact"), "output_artifact"
                )
                if artifact == output:
                    raise ValueError("Artifact checksum output must differ from its input")
            elif kind == "artifact_nonempty":
                if set(action) != {"action_id", "kind", "artifact"}:
                    raise ValueError("Launch-contract " + phase + " assertion action is malformed")
                _relative_artifact_path(action.get("artifact"), "artifact")
            else:
                raise ValueError(
                    "Launch-contract " + phase + " accepts only bounded built-in actions"
                )
    for action in actions:
        if type(action) is not dict or set(action) != {
            "action_id",
            "kind",
            "resources",
            "arguments",
            "stdout_artifact",
            "stderr_artifact",
        }:
            raise ValueError("Launch action has unexpected or missing fields")
        identifier = action.get("action_id")
        if not _lower_identifier(identifier):
            raise ValueError("Launch action ID is malformed")
        if identifier in identifiers:
            raise ValueError("Duplicate launch action ID: " + identifier)
        identifiers.add(identifier)
        if action.get("kind") != "athena":
            raise ValueError("Trusted trampoline accepts only Athena actions")
        resources = action.get("resources")
        if type(resources) is not dict or set(resources) != {
            "nodes",
            "tasks",
            "cpus_per_task",
            "gpus_per_task",
            "gpu_bind",
        }:
            raise ValueError("Launch-action resources are malformed")
        for field in ["nodes", "tasks", "cpus_per_task", "gpus_per_task"]:
            number = resources.get(field)
            if type(number) is not int or number <= 0:
                raise ValueError("Launch-action resources." + field + " must be a positive integer")
        if resources.get("gpu_bind") != "closest":
            raise ValueError("Launch-action resources.gpu_bind must be closest")
        arguments = action.get("arguments")
        if type(arguments) is not list:
            raise ValueError("Launch-action arguments must be an array")
        argument_index = 0
        while argument_index < len(arguments):
            argument = arguments[argument_index]
            if type(argument) is not dict or len(argument) != 1:
                raise ValueError("Launch-action argument must be one structured token")
            if "literal" in argument:
                literal = argument["literal"]
                if (
                    type(literal) is not str
                    or not literal
                    or "\0" in literal
                    or len(literal) > 4096
                ):
                    raise ValueError("Launch-action literal is malformed")
                if literal in {"-i", "-d"}:
                    argument_index += 1
                    if argument_index >= len(arguments):
                        raise ValueError("Launch-action " + literal + " requires one value")
                    structured = arguments[argument_index]
                    if literal == "-i" and structured != {"snapshot_role": "input-deck"}:
                        raise ValueError("Launch-action -i requires the frozen input deck")
                    if literal == "-d":
                        if (
                            type(structured) is not dict
                            or set(structured) != {"artifact_directory"}
                        ):
                            raise ValueError(
                                "Launch-action -d requires one artifact directory"
                            )
                        _relative_artifact_path(
                            structured["artifact_directory"], "artifact_directory"
                        )
                elif literal == "-r" or literal.startswith("-r="):
                    raise ValueError(
                        "Launch-action restart input is not authorized before trusted "
                        "restart snapshots are implemented"
                    )
                elif literal in {"-n", "-c"}:
                    pass
                elif literal.startswith("-"):
                    raise ValueError("Launch-action CLI flag is not authorized: " + literal)
                else:
                    key, separator, override = literal.partition("=")
                    if (
                        not separator
                        or not _override_key(key)
                        or not override
                        or "/" in override
                        or "\\" in override
                        or ".." in override
                    ):
                        raise ValueError(
                            "Launch-action Athena override is not authorized: " + literal
                        )
            else:
                raise ValueError(
                    "Launch-action structured path token must immediately follow -i or -d"
                )
            argument_index += 1
        stdout = _relative_artifact_path(action.get("stdout_artifact"), "stdout_artifact")
        stderr = _relative_artifact_path(action.get("stderr_artifact"), "stderr_artifact")
        if stdout == stderr:
            raise ValueError("Launch-action stdout and stderr artifacts must differ")
    return value


def launch_contract_sha256(value):
    contract = validate_launch_contract(value)
    payload = _canonical_json(contract).encode("utf-8")
    return _sha256_hex(payload)
'''


def _minimal_launch_syntax_surfaces() -> tuple[
    str,
    Callable[[object], dict[str, object]],
    Callable[[object], str],
    Callable[[], None],
]:
    runtime_guard = _require_trusted_runtime
    runtime_guard_code = runtime_guard.__code__
    runtime_globals = runtime_guard.__globals__
    function_type = _FUNCTION_TYPE
    require_condition = _require
    require_condition_code = require_condition.__code__
    error_type = _TRUSTED_VALUE_ERROR
    trusted_all = _TRUSTED_ALL
    trusted_dict = _TRUSTED_DICT
    trusted_len = _TRUSTED_LEN
    trusted_sorted = _TRUSTED_SORTED
    trusted_str = _TRUSTED_STR
    trusted_type = _TRUSTED_TYPE

    def require_exact_runtime() -> None:
        if not (
            runtime_guard.__code__ is runtime_guard_code
            and runtime_guard.__globals__ is runtime_globals
            and runtime_globals.get("_require_trusted_runtime") is runtime_guard
            and runtime_globals.get("_require") is require_condition
            and require_condition.__code__ is require_condition_code
            and require_condition.__globals__ is runtime_globals
            and runtime_globals.get("PreparationError") is error_type
        ):
            raise error_type("Q009 preparation runtime guard drifted")
        runtime_guard()

    require_exact_runtime_code = require_exact_runtime.__code__
    safe_builtins = _MAPPING_PROXY_TYPE(
        {
            "ValueError": _TRUSTED_VALUE_ERROR,
            "all": _TRUSTED_ALL,
            "any": _TRUSTED_ANY,
            "bool": _TRUSTED_BOOL,
            "dict": _TRUSTED_DICT,
            "int": _TRUSTED_INT,
            "len": _TRUSTED_LEN,
            "list": _TRUSTED_LIST,
            "ord": _TRUSTED_ORD,
            "range": _TRUSTED_RANGE,
            "set": _TRUSTED_SET,
            "sorted": _TRUSTED_SORTED,
            "str": _TRUSTED_STR,
            "type": _TRUSTED_TYPE,
            "zip": _TRUSTED_ZIP,
        }
    )
    namespace: dict[str, object] = {
        "__builtins__": safe_builtins,
        "__name__": MINIMAL_LAUNCH_SYNTAX_VALIDATOR_ID,
    }
    _TRUSTED_EXEC(
        _TRUSTED_COMPILE(
            _MINIMAL_LAUNCH_SYNTAX_VALIDATOR_SOURCE,
            f"<{MINIMAL_LAUNCH_SYNTAX_VALIDATOR_ID}>",
            "exec",
            dont_inherit=True,
        ),
        namespace,
    )
    validator = namespace["validate_launch_contract"]
    digest = namespace["launch_contract_sha256"]
    _require(
        _TRUSTED_CALLABLE(validator)
        and _TRUSTED_CALLABLE(digest)
        and validator.__module__ == MINIMAL_LAUNCH_SYNTAX_VALIDATOR_ID
        and digest.__module__ == MINIMAL_LAUNCH_SYNTAX_VALIDATOR_ID,
        "minimal launch-contract syntax implementation is malformed",
    )
    function_records = tuple(
        sorted(
            (
                name,
                value,
                value.__code__,
            )
            for name, value in namespace.items()
            if _TRUSTED_TYPE(value) is _FUNCTION_TYPE
        )
    )
    constant_records = tuple(
        sorted(
            (name, value)
            for name, value in namespace.items()
            if name != "__builtins__" and _TRUSTED_TYPE(value) is not _FUNCTION_TYPE
        )
    )
    namespace_keys = tuple(sorted(namespace))

    def require_exact_validator_state() -> None:
        if require_exact_runtime.__code__ is not require_exact_runtime_code:
            raise error_type("minimal launch-contract runtime guard drifted")
        require_exact_runtime()
        require_condition(
            trusted_type(namespace) is trusted_dict
            and trusted_len(namespace) == trusted_len(namespace_keys)
            and trusted_all(trusted_type(key) is trusted_str for key in namespace)
            and tuple(trusted_sorted(namespace)) == namespace_keys
            and namespace.get("__builtins__") is safe_builtins,
            "minimal launch-contract validator namespace drifted",
        )
        for name, expected, code in function_records:
            current = namespace.get(name)
            require_condition(
                current is expected
                and trusted_type(current) is function_type
                and current.__code__ is code
                and current.__globals__ is namespace
                and current.__builtins__ is safe_builtins,
                "minimal launch-contract validator code drifted",
            )
        for name, expected in constant_records:
            require_condition(
                namespace.get(name) is expected,
                "minimal launch-contract validator constant drifted",
            )

    require_exact_validator_state()
    return (
        _TRUSTED_STR(namespace["TRUSTED_LAUNCH_EXECUTOR"]),
        validator,
        digest,
        require_exact_validator_state,
    )


(
    MINIMAL_TRUSTED_LAUNCH_EXECUTOR,
    MINIMAL_VALIDATE_LAUNCH_CONTRACT,
    MINIMAL_LAUNCH_CONTRACT_SHA256,
    REQUIRE_EXACT_MINIMAL_VALIDATOR_STATE,
) = _minimal_launch_syntax_surfaces()
del _TRUSTED_COMPILE, _TRUSTED_EXEC, _minimal_launch_syntax_surfaces


def _load_captured_control_plane_launch_syntax() -> tuple[
    str,
    Callable[[object], dict[str, object]],
    Callable[[object], str],
    dict[str, object],
]:
    """Bind full source as reference and return the isolated minimal projection."""
    payload = _stable_repository_bytes(
        CONTROL_PLANE_LAUNCH_SYNTAX_REFERENCE_PATH,
        label="captured control-plane launch-syntax reference",
    )
    _require(
        b'def validate_launch_contract(value: object)' in payload
        and b'def launch_contract_sha256(value: object)' in payload
        and MINIMAL_TRUSTED_LAUNCH_EXECUTOR.encode("ascii") in payload,
        "captured control-plane launch-syntax reference is malformed",
    )
    bindings = {
        "control_plane_common_reference": _captured_source_binding(
            CONTROL_PLANE_LAUNCH_SYNTAX_REFERENCE_PATH,
            payload,
            loading_role="captured_verified_reference_not_executed",
        )
    }
    return (
        MINIMAL_TRUSTED_LAUNCH_EXECUTOR,
        MINIMAL_VALIDATE_LAUNCH_CONTRACT,
        MINIMAL_LAUNCH_CONTRACT_SHA256,
        bindings,
    )


CAPTURED_CONTRACT_VALIDATOR_BINDING = _captured_source_binding(
    CONTRACT_VALIDATOR_PATH,
    _stable_repository_bytes(
        CONTRACT_VALIDATOR_PATH,
        label="captured Q009 contract validator reference",
    ),
    loading_role="captured_verified_contract_validator_reference_not_executed",
)
(
    TRUSTED_LAUNCH_EXECUTOR,
    _INTERNAL_VALIDATE_LAUNCH_CONTRACT,
    _INTERNAL_LAUNCH_CONTRACT_SHA256,
    CAPTURED_CONTROL_PLANE_REFERENCE_BINDINGS,
) = (
    _load_captured_control_plane_launch_syntax()
)


def _make_public_launch_syntax_surfaces(
    validator,
    digest,
    require_exact_validator_state,
    exact_json_graph,
):
    exact_json_code = exact_json_graph.__code__
    validator_state_code = require_exact_validator_state.__code__
    error_type = _TRUSTED_VALUE_ERROR

    def validate_launch_contract(value):
        if require_exact_validator_state.__code__ is not validator_state_code:
            raise error_type("minimal launch-contract validator guard drifted")
        require_exact_validator_state()
        if exact_json_graph.__code__ is not exact_json_code or not exact_json_graph(value):
            raise error_type("Launch contract must be an exact JSON graph")
        return validator(value)

    def launch_contract_sha256(value):
        if require_exact_validator_state.__code__ is not validator_state_code:
            raise error_type("minimal launch-contract validator guard drifted")
        require_exact_validator_state()
        if exact_json_graph.__code__ is not exact_json_code or not exact_json_graph(value):
            raise error_type("Launch contract must be an exact JSON graph")
        return digest(value)

    return validate_launch_contract, launch_contract_sha256, require_exact_validator_state


(
    validate_launch_contract,
    launch_contract_sha256,
    _PUBLIC_REQUIRE_EXACT_VALIDATOR_STATE,
) = _make_public_launch_syntax_surfaces(
    _INTERNAL_VALIDATE_LAUNCH_CONTRACT,
    _INTERNAL_LAUNCH_CONTRACT_SHA256,
    REQUIRE_EXACT_MINIMAL_VALIDATOR_STATE,
    _is_exact_json_graph,
)
del (
    _INTERNAL_LAUNCH_CONTRACT_SHA256,
    _INTERNAL_VALIDATE_LAUNCH_CONTRACT,
    MINIMAL_LAUNCH_CONTRACT_SHA256,
    MINIMAL_VALIDATE_LAUNCH_CONTRACT,
    REQUIRE_EXACT_MINIMAL_VALIDATOR_STATE,
    _load_captured_control_plane_launch_syntax,
    _make_public_launch_syntax_surfaces,
)


def _parse_deck(text: str) -> tuple[list[str], dict[str, dict[str, str]]]:
    _require(type(text) is str, "deck parser requires an exact string")
    order: list[str] = []
    blocks: dict[str, dict[str, str]] = {}
    current: str | None = None
    for line_number, raw in enumerate(text.splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            current = line[1:-1].strip()
            _require(current and current not in blocks, f"template:{line_number}: duplicate block")
            order.append(current)
            blocks[current] = {}
            continue
        name, separator, value = line.partition("=")
        name = name.strip()
        value = value.strip()
        _require(
            current is not None
            and separator == "="
            and bool(name)
            and bool(value)
            and all(not character.isspace() and character != "=" for character in name),
            f"template:{line_number}: malformed line",
        )
        _require(name not in blocks[current] and bool(value), f"template:{line_number}: duplicate parameter")
        blocks[current][name] = value
    return order, blocks


def _render_deck(order: Sequence[str], blocks: Mapping[str, Mapping[str, str]], *, case_id: str, stage_id: str) -> str:
    lines = [
        "# Q009 REGISTERED ENGINEERING PILOT CANDIDATE: aggregate-only, non-authorizing.",
        f"# logical_case={case_id} scheduler_stage={stage_id}",
        "# Physical field and particle outputs are intentionally disabled.",
        "",
    ]
    for block in order:
        if block not in blocks:
            continue
        lines.append(f"<{block}>")
        for name, value in blocks[block].items():
            lines.append(f"{name:<38} = {value}")
        lines.append("")
    return "\n".join(lines)


def _template_payload(case: Mapping[str, object]) -> tuple[bytes, str]:
    frozen = _load_frozen_contract()
    binding = frozen["source_templates"][str(case["template_id"])]
    payload = _stable_repository_bytes(
        str(binding["path"]), label=f"Q009 source template {case['template_id']}"
    )
    _require(_sha256_bytes(payload) == binding["sha256"], "Q009 source template digest drifted")
    try:
        return payload, payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise PreparationError("Q009 source template is not UTF-8") from error


def _set(blocks: dict[str, dict[str, str]], block: str, name: str, value: str) -> None:
    _require(block in blocks and name in blocks[block], f"deck override references missing {block}/{name}")
    blocks[block][name] = value


def _remove_outputs(order: list[str], blocks: dict[str, dict[str, str]]) -> None:
    for block in list(order):
        if block.startswith("output"):
            order.remove(block)
            del blocks[block]


def _stage_ids(case: Mapping[str, object]) -> tuple[str, ...]:
    return (
        ("checkpoint", "continuation")
        if case["restart_profile_id"] == "checkpoint_before45_resume_after45"
        else ("main",)
    )


def _materialized_deck(case: Mapping[str, object], stage: str) -> tuple[bytes, dict[str, dict[str, str]]]:
    _, text = _template_payload(case)
    order, blocks = _parse_deck(text)
    _remove_outputs(order, blocks)
    case_id = str(case["case_id"])
    stage_id = f"{case_id}-{stage}"
    _set(blocks, "job", "basename", stage_id.replace("-", "_"))
    _set(
        blocks,
        "particles",
        "pic_load_balance_cost_per_particle",
        format(float(case["load_balance"]["pic_load_balance_cost_per_particle"]), ".17g"),
    )
    _set(blocks, "particles", "pic_random_seed", str(case["random_seed"]))
    if case["template_id"] == "q011_section54_production_preparation":
        _set(blocks, "problem", "ps_inject_seed", str(case["random_seed"]))
        _set(blocks, "problem", "ps_seed_noise_seed", str(case["random_seed"]))
    geometry = str(case["geometry_profile_id"])
    if geometry == "reduced_transverse_shock":
        _set(blocks, "mesh", "nx2", "20")
        _set(blocks, "mesh", "x2max", "240.0")
        _set(blocks, "time", "tlim", "44.0" if stage == "checkpoint" else "50.0")
    elif geometry == "full_section54_geometry":
        exact = _load_frozen_contract()["profile_definitions"]["geometry_profiles"][
            geometry
        ]
        for key, value in exact["exact_geometry"].items():
            block, name = key.split("/", 1)
            _set(blocks, block, name, str(value))
        _set(blocks, "time", "nlim", "512")
    elif geometry != "reduced_lifetime":
        raise PreparationError(f"unsupported Q009 geometry profile: {geometry}")
    if stage == "checkpoint":
        blocks["output1"] = {
            "file_type": "rst",
            "id": "q009_checkpoint",
            "dt": "44.0",
            "single_file_per_rank": "false",
        }
        order.append("output1")
    _validate_deck(case, stage, blocks)
    return _render_deck(order, blocks, case_id=case_id, stage_id=stage_id).encode("utf-8"), blocks


def _validate_deck(
    case: Mapping[str, object], stage: str, blocks: Mapping[str, Mapping[str, str]]
) -> None:
    outputs = {name: values for name, values in blocks.items() if name.startswith("output")}
    if stage == "checkpoint":
        _require(
            outputs == {
                "output1": {
                    "file_type": "rst",
                    "id": "q009_checkpoint",
                    "dt": "44.0",
                    "single_file_per_rank": "false",
                }
            },
            "checkpoint deck output contract drifted",
        )
    else:
        _require(not outputs, "aggregate-only Q009 deck retained a physical output")
    _require(
        blocks["particles"]["pic_random_seed"] == str(case["random_seed"])
        and float(blocks["particles"]["pic_load_balance_cost_per_particle"])
        == float(case["load_balance"]["pic_load_balance_cost_per_particle"]),
        "Q009 deck seed or load-balance weight drifted",
    )
    geometry = case["geometry_profile_id"]
    if geometry == "reduced_transverse_shock":
        _require(
            blocks["mesh"]["nx2"] == "20"
            and blocks["mesh"]["x2max"] == "240.0"
            and blocks["time"]["tlim"] == ("44.0" if stage == "checkpoint" else "50.0"),
            "reduced-transverse Q009 geometry drifted",
        )
    if geometry == "full_section54_geometry":
        _require(
            blocks["mesh"]["nx1"] == "4000"
            and blocks["mesh"]["nx2"] == "260"
            and blocks["mesh"]["x1max"] == "48000.0"
            and blocks["mesh"]["x2max"] == "3120.0"
            and blocks["meshblock"]["nx1"] == "20"
            and blocks["meshblock"]["nx2"] == "20"
            and blocks["time"]["nlim"] == "512",
            "held-out full geometry drifted",
        )


def _launch_contract(case: Mapping[str, object], stage: str) -> dict[str, object]:
    resources = case["resources"]
    action_id = f"q009-{case['sequence']}-{stage}"
    value = {
        "schema_version": 1,
        "executor": TRUSTED_LAUNCH_EXECUTOR,
        "pre_actions": [
            {
                "action_id": "sha-input-deck",
                "kind": "snapshot_sha256",
                "snapshot_role": "input-deck",
                "output_artifact": "bindings/input_deck.sha256",
            },
            {
                "action_id": "sha-executable",
                "kind": "snapshot_sha256",
                "snapshot_role": "executable",
                "output_artifact": "bindings/executable.sha256",
            },
        ],
        "actions": [
            {
                "action_id": action_id,
                "kind": "athena",
                "resources": {
                    "nodes": resources["nodes"],
                    "tasks": resources["mpi_ranks"],
                    "cpus_per_task": resources["cpus_per_rank"],
                    "gpus_per_task": resources["gpus_per_rank"],
                    "gpu_bind": "closest",
                },
                "arguments": [
                    {"literal": "-i"},
                    {"snapshot_role": "input-deck"},
                    {"literal": "-d"},
                    {"artifact_directory": "raw"},
                ],
                "stdout_artifact": "athena_stdout.txt",
                "stderr_artifact": "athena_stderr.txt",
            }
        ],
        "post_actions": [
            {
                "action_id": "require-stdout",
                "kind": "artifact_nonempty",
                "artifact": "athena_stdout.txt",
            },
            {
                "action_id": "sha-stdout",
                "kind": "artifact_sha256",
                "artifact": "athena_stdout.txt",
                "output_artifact": "bindings/athena_stdout.sha256",
            },
        ],
    }
    return validate_launch_contract(value)


def _artifact_root_template() -> str:
    return f"{AUTHORIZED_ORION_ROOT}/runs/{CAMPAIGN}/{SUBMISSION_ID_TEMPLATE}"


def _review_candidate_id(case_id: str, stage: str) -> str:
    return f"q009-review-{case_id}-{stage}-v2"


def _prior_stage_ids(
    cases: Sequence[Mapping[str, object]], case_index: int, stage_index: int
) -> list[str]:
    result: list[str] = []
    for prior_case in cases[:case_index]:
        result.extend(f"{prior_case['case_id']}:{stage}" for stage in _stage_ids(prior_case))
    result.extend(
        f"{cases[case_index]['case_id']}:{stage}"
        for stage in _stage_ids(cases[case_index])[:stage_index]
    )
    return result


def _launch_candidate(
    *,
    cases: Sequence[Mapping[str, object]],
    case_index: int,
    stage_index: int,
    deck: Mapping[str, object],
) -> dict[str, object]:
    case = cases[case_index]
    stage = _stage_ids(case)[stage_index]
    launchable = stage != "continuation"
    launch = _launch_contract(case, stage) if launchable else None
    case_id = str(case["case_id"])
    stage_id = f"{case_id}:{stage}"
    blockers = list(COMMON_RUNTIME_BLOCKERS)
    nonpromotable_blockers: list[str] = []
    if case["template_id"] == "q009_coupled_boundary_lifetime":
        blockers.append(MIGRATION_DESIGN_BLOCKER)
        nonpromotable_blockers.append(MIGRATION_DESIGN_BLOCKER)
    if not launchable:
        blockers.append(RESTART_CONTINUATION_BLOCKER)
        nonpromotable_blockers.append(RESTART_CONTINUATION_BLOCKER)
    return {
        "record_type": LAUNCH_RECORD_TYPE,
        "schema_version": SCHEMA_VERSION,
        "successor_id": SUCCESSOR_ID,
        "qualification_effect": QUALIFICATION_EFFECT,
        "authority_semantics": AUTHORITY_SEMANTICS,
        "preparation_role": "source_local_non_operational_review_only",
        "identity": {
            "campaign": CAMPAIGN,
            "logical_case_id": case_id,
            "case_sequence": case["sequence"],
            "scheduler_stage": stage,
            "stage_id": stage_id,
            "stage_sequence_within_case": stage_index + 1,
            "scheduler_jobs_per_logical_attempt": case["resources"][
                "scheduler_jobs_per_attempt"
            ],
            "maximum_logical_attempts": MAXIMUM_ATTEMPTS_PER_CASE,
            "maximum_retries": MAXIMUM_RETRIES_PER_CASE,
            "review_candidate_id": _review_candidate_id(case_id, stage),
        },
        "future_operational_topology_review_only": {
            "artifact_dir_template": _artifact_root_template(),
            "raw_root_template": f"{_artifact_root_template()}/raw",
            "analysis_root_template": f"{_artifact_root_template()}/analysis",
            "artifact_inventory_template": f"{_artifact_root_template()}/artifact_inventory.json",
            "required_shape": "runs/<campaign>/<submission-id>",
            "campaign_is_single_path_segment": True,
            "case_id_is_bound_as_test_id_not_as_an_extra_run_path_component": True,
        },
        "logical_case_contract": dict(case),
        "deck": dict(deck),
        "operational_bindings_accepted": False,
        "unresolved_operational_bindings": list(UNRESOLVED_OPERATIONAL_BINDINGS),
        "prior_stage_closures_required": _prior_stage_ids(cases, case_index, stage_index),
        "launch_contract_candidate": launch,
        "launch_contract_sha256": launch_contract_sha256(launch) if launch is not None else None,
        "launch_contract_syntax_validated_by_minimal_projection": launchable,
        "restart_continuation_handoff_required": not launchable,
        "nonpromotable_design_blockers": nonpromotable_blockers,
        "resource_ceiling": {
            "nodes": case["resources"]["nodes"],
            "mpi_ranks": case["resources"]["mpi_ranks"],
            "gpus_per_node": case["resources"]["gpus_per_node"],
            "gpus_per_rank": case["resources"]["gpus_per_rank"],
            "cpus_per_rank": case["resources"]["cpus_per_rank"],
            "walltime_seconds": case["resources"]["walltime_seconds_per_scheduler_job"],
            "logical_case_maximum_node_hours": case["ceilings"]["maximum_node_hours"],
            "logical_case_maximum_artifact_storage_gib": case["ceilings"][
                "maximum_artifact_storage_gib"
            ],
        },
        "required_blockers": blockers,
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }


def _restart_handoff(launch: Mapping[str, object]) -> dict[str, object] | None:
    if not launch["restart_continuation_handoff_required"]:
        return None
    identity = launch["identity"]
    case = launch["logical_case_contract"]
    return {
        "record_type": HANDOFF_RECORD_TYPE,
        "schema_version": SCHEMA_VERSION,
        "successor_id": SUCCESSOR_ID,
        "status": "launch_prohibited_pending_trusted_restart_snapshot_successor",
        "qualification_effect": QUALIFICATION_EFFECT,
        "authority_semantics": AUTHORITY_SEMANTICS,
        "campaign": CAMPAIGN,
        "logical_case_id": identity["logical_case_id"],
        "scheduler_stage": identity["scheduler_stage"],
        "review_candidate_id": identity["review_candidate_id"],
        "checkpoint_source_stage_id": f"{identity['logical_case_id']}:checkpoint",
        "required_checkpoint_time_window": {
            "minimum_inclusive": 40.0,
            "maximum_exclusive": 45.0,
        },
        "required_resumed_terminal_time_window": {
            "minimum_exclusive": 45.0,
            "maximum_inclusive": 50.0,
        },
        "required_checkpoint_bindings": [
            "sealed_checkpoint_artifact_inventory_sha256",
            "exact_checkpoint_artifact_sha256",
            "checkpoint_observed_time",
            "checkpoint_observed_cycle",
            "pre_restart_state_fingerprint",
            "pre_restart_integer_counter_digest",
            "checkpoint_scheduler_job_id",
        ],
        "future_argv_shape": [
            "-r",
            "<trusted_restart_snapshot_role>",
            "-d",
            "raw",
            "time/tlim=50.0",
        ],
        "minimal_projection_syntax_validation_result": (
            "minimal_projection_rejects_restart_input_by_design"
        ),
        "logical_case_maximum_node_hours": case["ceilings"]["maximum_node_hours"],
        "required_blocker": RESTART_CONTINUATION_BLOCKER,
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }


def _policy_candidate(launch: Mapping[str, object]) -> dict[str, object]:
    identity = launch["identity"]
    case = launch["logical_case_contract"]
    policy_slice = {
        "review_candidate_id": identity["review_candidate_id"],
        "status": (
            "blocked_not_promotable_design_or_launch_syntax_missing"
            if launch["nonpromotable_design_blockers"]
            else "review_only_missing_all_operational_bindings"
        ),
        "structurally_live_policy_compatible": False,
        "campaign": CAMPAIGN,
        "test_id": str(identity["logical_case_id"]).replace("-", "_")
        + "_"
        + str(identity["scheduler_stage"]),
        "evidence_class": EVIDENCE_CLASS,
        "physical_mode": PHYSICAL_MODE,
        "runtime_profile": RUNTIME_PROFILE,
        "selected_qos": case["resources"]["qos"],
        "registered_short_nonproduction": False,
        "maximum_nodes": case["resources"]["nodes"],
        "maximum_walltime_seconds": case["resources"]["walltime_seconds_per_scheduler_job"],
        "maximum_attempts": 1,
        "input_deck_sha256": launch["deck"]["sha256"],
        "launch_contract_sha256": launch["launch_contract_sha256"],
        "unresolved_operational_bindings": list(UNRESOLVED_OPERATIONAL_BINDINGS),
    }
    return {
        "record_type": POLICY_RECORD_TYPE,
        "schema_version": SCHEMA_VERSION,
        "successor_id": SUCCESSOR_ID,
        "qualification_effect": QUALIFICATION_EFFECT,
        "authority_semantics": AUTHORITY_SEMANTICS,
        "logical_case_id": identity["logical_case_id"],
        "scheduler_stage": identity["scheduler_stage"],
        "policy_slice_candidate": policy_slice,
        "logical_attempt_guard": {
            "maximum_logical_attempts_per_case": 1,
            "maximum_retries_per_case": 0,
            "scheduler_jobs_in_this_logical_case": case["resources"][
                "scheduler_jobs_per_attempt"
            ],
            "serial_prior_stage_closures_required": launch["prior_stage_closures_required"],
        },
        "live_policy_mutation_authorized": False,
        "required_blockers": launch["required_blockers"],
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }


def _budget(cases: Sequence[Mapping[str, object]], stage_count: int) -> dict[str, object]:
    node_hours = sum(float(case["ceilings"]["maximum_node_hours"]) for case in cases)
    storage = sum(float(case["ceilings"]["maximum_artifact_storage_gib"]) for case in cases)
    _require(
        node_hours == 100.0 and storage == 200.0,
        "Q009 contract budget or storage sum drifted",
    )
    return {
        "logical_case_count": len(cases),
        "scheduler_stage_count": stage_count,
        "maximum_live_q009_submissions": MAXIMUM_LIVE_Q009_SUBMISSIONS,
        "maximum_logical_attempts_per_case": MAXIMUM_ATTEMPTS_PER_CASE,
        "maximum_retries_per_case": MAXIMUM_RETRIES_PER_CASE,
        "contract_hard_cap_node_hours": node_hours,
        "contract_hard_cap_artifact_storage_gib": storage,
        "project_frontier_cap_node_hours": 10000.0,
        "fresh_live_ledger_other_campaign_reserve_and_storage_preflight_required": True,
        "empty_user_queue_required_before_every_serial_stage": True,
        "ledger_mutation_authorized": False,
    }


def _prepare_materialization_once() -> tuple[dict[str, object], dict[str, bytes]]:
    """Build all deterministic Q009 decks and review candidates in memory."""
    _require_trusted_runtime()
    frozen = _load_frozen_contract()
    cases = [dict(case) for case in frozen["pilot_matrix"]]
    files: dict[str, bytes] = {}
    records: list[dict[str, object]] = []
    policy_slices: list[dict[str, object]] = []
    stage_count = 0
    for case_index, case in enumerate(cases):
        case_id = str(case["case_id"])
        for stage_index, stage in enumerate(_stage_ids(case)):
            stage_count += 1
            deck_payload, _ = _materialized_deck(case, stage)
            deck_path = f"decks/{case_id}/{stage}.athinput"
            files[deck_path] = deck_payload
            deck = _binding(deck_path, deck_payload)
            launch = _launch_candidate(
                cases=cases,
                case_index=case_index,
                stage_index=stage_index,
                deck=deck,
            )
            policy = _policy_candidate(launch)
            handoff = _restart_handoff(launch)
            launch_path = f"launch_candidates/{case_id}/{stage}.json"
            policy_path = f"policy_candidates/{case_id}/{stage}.json"
            launch_payload = _json_bytes(launch)
            policy_payload = _json_bytes(policy)
            files[launch_path] = launch_payload
            files[policy_path] = policy_payload
            record: dict[str, object] = {
                "logical_case_id": case_id,
                "case_sequence": case["sequence"],
                "scheduler_stage": stage,
                "global_stage_sequence": stage_count,
                "deck": deck,
                "launch_candidate": _binding(launch_path, launch_payload),
                "policy_candidate": _binding(policy_path, policy_payload),
                "launch_contract_syntax_candidate_present": launch[
                    "launch_contract_syntax_validated_by_minimal_projection"
                ],
            }
            if handoff is not None:
                handoff_path = f"restart_handoffs/{case_id}.json"
                handoff_payload = _json_bytes(handoff)
                files[handoff_path] = handoff_payload
                record["restart_continuation_handoff"] = _binding(
                    handoff_path, handoff_payload
                )
            records.append(record)
            policy_slices.append(policy["policy_slice_candidate"])
    budget = _budget(cases, stage_count)
    budget_payload = _json_bytes(budget)
    files["batch_budget_accounting_input.json"] = budget_payload
    policy_fragment = {
        "record_type": "q009_dynamic_amr_registered_policy_slice_review_fragment",
        "schema_version": SCHEMA_VERSION,
        "successor_id": SUCCESSOR_ID,
        "status": "review_fragment_only_not_live_policy",
        "qualification_effect": QUALIFICATION_EFFECT,
        "authority_semantics": AUTHORITY_SEMANTICS,
        "logical_case_count": len(cases),
        "scheduler_stage_count": stage_count,
        "non_operational_policy_slice_review_candidates": policy_slices,
        "restart_continuation_slice_count_blocked_not_promotable": 2,
        "migration_slice_count_blocked_not_promotable": 2,
        "live_policy_mutation_authorized": False,
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }
    policy_fragment_payload = _json_bytes(policy_fragment)
    files["registered_policy_slice_review_fragment.json"] = policy_fragment_payload
    manifest = {
        "record_type": RECORD_TYPE,
        "schema_version": SCHEMA_VERSION,
        "successor_id": SUCCESSOR_ID,
        "status": "source_local_preparation_complete_execution_blocked",
        "qualification_effect": QUALIFICATION_EFFECT,
        "authority_semantics": AUTHORITY_SEMANTICS,
        "threat_model": dict(THREAT_MODEL),
        "operational_binding_input_supported": False,
        "exact_launch_candidate_validation_supported": False,
        "unresolved_operational_bindings": list(UNRESOLVED_OPERATIONAL_BINDINGS),
        "source_bindings": {
            "q009_contract": _binding(
                CONTRACT_PATH,
                _stable_repository_bytes(
                    CONTRACT_PATH,
                    label="Q009 contract source binding",
                ),
            ),
            "q009_contract_validator_reference_not_executed": dict(
                CAPTURED_CONTRACT_VALIDATOR_BINDING
            ),
            "captured_control_plane_common_reference_not_executed": dict(
                CAPTURED_CONTROL_PLANE_REFERENCE_BINDINGS[
                    "control_plane_common_reference"
                ]
            ),
            "execution_preparation_source": _binding(
                "tst/publication/q009_deterministic_execution_preparation_only_successor_v2.py",
                _stable_repository_bytes(
                    "tst/publication/q009_deterministic_execution_preparation_only_successor_v2.py",
                    label="Q009 execution preparation source",
                ),
            ),
        },
        "logical_case_count": len(cases),
        "scheduler_stage_count": stage_count,
        "stage_records": records,
        "batch_budget_accounting_input": _binding(
            "batch_budget_accounting_input.json", budget_payload
        ),
        "registered_policy_slice_review_fragment": _binding(
            "registered_policy_slice_review_fragment.json", policy_fragment_payload
        ),
        "minimal_structured_syntax_compatibility": {
            "run_topology": "runs/<campaign>/<submission-id>",
            "minimal_validator_identity": MINIMAL_LAUNCH_SYNTAX_VALIDATOR_ID,
            "captured_control_plane_reference_role": (
                "captured_verified_reference_not_executed"
            ),
            "exposed_minimal_surfaces": [
                "validate_launch_contract",
                "launch_contract_sha256",
            ],
            "operational_control_plane_module_exposed": False,
            "aggregate_only_main_and_checkpoint_launch_contract_syntax_validates": (
                True
            ),
            "transitive_sys_modules_imports_used_by_minimal_validator": False,
            "importable_standard_modules_used_by_preparation": False,
            "exact_contract_json_parsed_by_internal_digest_bound_decoder": True,
            "captured_contract_validator_executed": False,
            "minimal_validator_uses_captured_builtin_primitives": True,
            "minimal_validator_builtin_mapping_rejects_ordinary_assignment": True,
            "minimal_validator_builtin_mapping_recoverable_and_mutable_by_privileged_same_process_python": (
                True
            ),
            "minimal_validator_namespace_and_code_identity_checks_present": True,
            "runtime_builtin_primitive_drift_fails_closed": True,
            "privileged_same_process_can_synchronize_guard_roots": True,
            "spoofed_heap_type_static_builtin_wrappers_fail_before_execution": True,
            "reachable_internal_guard_code_identity_checks_present": True,
            "retained_mutable_trusted_primitive_globals_present": False,
            "post_freeze_construction_or_repository_read_helpers_exposed": False,
            "post_freeze_mutable_primitive_verification_roots_exposed": False,
            "retained_module_local_callable_count": 7,
            "retained_runtime_guard_uses_identity_pinned_non_write_non_exec_primitives": (
                True
            ),
            "in_process_integrity_guards_are_authority_boundary": False,
            "in_process_integrity_guards_role": (
                "ordinary_accidental_drift_detection_only"
            ),
            "non_exact_json_rejected_before_coercion_or_protocol_calls": True,
            "cyclic_or_aliased_container_graphs_rejected_as_non_json": True,
            "hostile_loader_file_rejected_before_protocol_calls": True,
            "preverification_imports_or_class_definitions_present": False,
            "import_time_preparation_graph_cached": True,
            "import_time_preparation_graph_authoritative": False,
            "ordinary_module_global_drift_changes_cached_preparation": False,
            "privileged_same_process_closure_or_mapping_proxy_mutation_outside_contract": (
                True
            ),
            "validator_process_isolation_present": False,
            "authoritative_output_or_receipt_produced": False,
            "all_public_outputs_non_authoritative_review_inputs": True,
            "authority_requires_separate_external_verifier_and_admission": True,
            "operational_bindings_accepted": False,
            "exact_launch_candidate_validator_present": False,
            "trusted_restart_input_supported": False,
            "trusted_unique_rst_restart_launching_supported": False,
            "contract_complete_event_per_rank_per_level_telemetry_present": False,
            "pairwise_lb_on_off_acceptance_present": False,
            "hardened_external_completion_root_inventory_receipts_present": False,
            "receipt_producer_validator_present": False,
            "admission_bridge_present": False,
            "persistent_review_bundle_publication_supported": False,
            "persistent_review_bundle_publication_api_present": False,
            "filesystem_write_or_destructive_cleanup_path_present": False,
            "persistent_review_bundle_blocker": PERSISTENT_REVIEW_BUNDLE_BLOCKER,
            "reduced_migration_profile_single_job_design_closes_required_events": False,
            "restart_continuation_handoffs_fail_closed": True,
            "restart_continuation_blocker": RESTART_CONTINUATION_BLOCKER,
            "migration_design_blocker": MIGRATION_DESIGN_BLOCKER,
        },
        "absent_capabilities": dict(ABSENT_CAPABILITIES),
        "inspection_boundary": {
            "physical_field_outputs_materialized": False,
            "particle_outputs_materialized": False,
            "restart_payload_inspection_authorized": False,
            "runtime_output_inspection_authorized": False,
            "runtime_telemetry_implemented": False,
        },
        "execution_boundary": {
            "empty_user_queue_assumed": False,
            "live_policy_complete": False,
            "scheduler_calls_authorized": False,
            "launches_jobs": False,
            **AUTHORIZATION_BOUNDARY,
        },
    }
    return manifest, files


def _make_frozen_materialization_surfaces(
    manifest_bytes,
    file_items,
    decoder,
    exact_json_graph,
    exact_file_graph,
    strict_equal,
    require_exact_validator_state,
):
    decoder_code = decoder.__code__
    exact_json_code = exact_json_graph.__code__
    exact_file_code = exact_file_graph.__code__
    strict_equal_code = strict_equal.__code__
    strict_equal_globals = strict_equal.__globals__
    decoder_globals = decoder.__globals__
    json_hex = decoder_globals.get("_JSON_HEX")
    validator_state_code = require_exact_validator_state.__code__
    error_type = _TRUSTED_VALUE_ERROR

    def require_exact_surface_state() -> None:
        if require_exact_validator_state.__code__ is not validator_state_code:
            raise error_type("Q009 frozen preparation validator guard drifted")
        require_exact_validator_state()
        if not (
            decoder.__code__ is decoder_code
            and decoder.__globals__ is decoder_globals
            and decoder_globals.get("_JSON_HEX") is json_hex
            and exact_json_graph.__code__ is exact_json_code
            and exact_file_graph.__code__ is exact_file_code
            and strict_equal.__code__ is strict_equal_code
            and strict_equal.__globals__ is strict_equal_globals
            and strict_equal_globals.get("_is_exact_json_graph") is exact_json_graph
        ):
            raise error_type("Q009 frozen preparation surface drifted")

    surface_state_code = require_exact_surface_state.__code__

    def build_materialization():
        """Return a fresh copy of the cached non-authoritative Q009 review graph."""
        if require_exact_surface_state.__code__ is not surface_state_code:
            raise error_type("Q009 frozen preparation guard drifted")
        require_exact_surface_state()
        manifest = decoder(manifest_bytes, label="frozen Q009 preparation manifest")
        if not exact_json_graph(manifest):
            raise error_type("frozen Q009 preparation manifest is not exact JSON")
        return manifest, {path: payload for path, payload in file_items}

    def validate_materialization(manifest, files):
        """Compare against the cached non-authoritative Q009 review graph."""
        if require_exact_surface_state.__code__ is not surface_state_code:
            raise error_type("Q009 frozen preparation guard drifted")
        require_exact_surface_state()
        if not exact_json_graph(manifest):
            raise error_type("Q009 preparation manifest is not an exact JSON graph")
        if not exact_file_graph(files):
            raise error_type("Q009 preparation files are not an exact byte graph")
        expected_manifest = decoder(
            manifest_bytes, label="frozen Q009 preparation manifest"
        )
        if not strict_equal(manifest, expected_manifest):
            raise error_type("Q009 preparation manifest drifted")
        if len(files) != len(file_items):
            raise error_type("Q009 preparation inventory drifted")
        for path, payload in file_items:
            if path not in files or files[path] != payload:
                raise error_type("Q009 preparation member drifted: " + path)
        return expected_manifest, {path: payload for path, payload in file_items}

    return build_materialization, validate_materialization


_FROZEN_MANIFEST, _FROZEN_FILES = _prepare_materialization_once()
_FROZEN_MANIFEST_BYTES = _json_bytes(_FROZEN_MANIFEST)
_FROZEN_FILE_ITEMS = tuple(sorted(_FROZEN_FILES.items()))
build_materialization, validate_materialization = _make_frozen_materialization_surfaces(
    _FROZEN_MANIFEST_BYTES,
    _FROZEN_FILE_ITEMS,
    _decode_json,
    _is_exact_json_graph,
    _is_exact_file_graph,
    _strict_equal,
    _PUBLIC_REQUIRE_EXACT_VALIDATOR_STATE,
)
del (
    _FROZEN_FILES,
    _FROZEN_MANIFEST,
    _FROZEN_FILE_ITEMS,
    _FROZEN_MANIFEST_BYTES,
    _PUBLIC_REQUIRE_EXACT_VALIDATOR_STATE,
    _make_frozen_materialization_surfaces,
    _prepare_materialization_once,
)
del (
    Callable,
    Mapping,
    Sequence,
    _BUILTIN_FUNCTION_TYPE,
    _FUNCTION_TYPE,
    _MAPPING_PROXY_TYPE,
    _PY_TPFLAGS_HEAPTYPE,
    _TYPE_TYPE,
    _TRUSTED_ALL,
    _TRUSTED_ANY,
    _TRUSTED_BOOL,
    _TRUSTED_BUILD_CLASS,
    _TRUSTED_BYTES,
    _TRUSTED_CALLABLE,
    _TRUSTED_CHR,
    _TRUSTED_DICT,
    _TRUSTED_ENUMERATE,
    _TRUSTED_FLOAT,
    _TRUSTED_FORMAT,
    _TRUSTED_GETATTR,
    _TRUSTED_INT,
    _TRUSTED_ISINSTANCE,
    _TRUSTED_LEN,
    _TRUSTED_LIST,
    _TRUSTED_ORD,
    _TRUSTED_OS_ERROR,
    _TRUSTED_RANGE,
    _TRUSTED_SET,
    _TRUSTED_SORTED,
    _TRUSTED_STR,
    _TRUSTED_SUM,
    _TRUSTED_TYPE,
    _TRUSTED_UNICODE_DECODE_ERROR,
    _TRUSTED_VALUE_ERROR,
    _TRUSTED_ZIP,
)
del (
    _JSON_ESCAPES,
    _MINIMAL_LAUNCH_SYNTAX_VALIDATOR_SOURCE,
    _SHA256_INITIAL,
    _SHA256_ROUND,
    _artifact_root_template,
    _binding,
    _budget,
    _canonical_json,
    _captured_source_binding,
    _decode_json,
    _is_exact_builtin_function,
    _is_exact_file_graph,
    _is_exact_static_builtin_type,
    _json_bytes,
    _json_hex4,
    _json_string,
    _launch_candidate,
    _launch_contract,
    _load_frozen_contract,
    _materialized_deck,
    _parse_deck,
    _policy_candidate,
    _prior_stage_ids,
    _remove_outputs,
    _render_deck,
    _restart_handoff,
    _review_candidate_id,
    _rotate_right,
    _set,
    _sha256_bytes,
    _stable_repository_bytes,
    _stage_ids,
    _strict_equal,
    _template_payload,
    _validate_deck,
)
