#!/bin/bash
#SBATCH --account=AST207
#SBATCH --partition=batch
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --job-name=pic-q011-build-freeze
#SBATCH --output=/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/%x.%j.log

set -euo pipefail

export PATH=/usr/bin:/bin
unset BASH_ENV ENV GIT_ALTERNATE_OBJECT_DIRECTORIES GIT_CONFIG_COUNT
unset GIT_CONFIG_GLOBAL GIT_CONFIG_KEY_0 GIT_CONFIG_NOSYSTEM GIT_CONFIG_SYSTEM
unset GIT_CONFIG_VALUE_0 GIT_DIR GIT_INDEX_FILE GIT_OBJECT_DIRECTORY GIT_WORK_TREE
unset PYTHONHOME PYTHONINSPECT PYTHONPATH PYTHONSTARTUP PYTHONUSERBASE
unset TMPDIR

PIC_ROOT=/lustre/orion/ast207/proj-shared/dfielding/PIC
SRC_DIR=/ccs/home/dfielding/athenak-pic
USAGE="usage: sbatch $0 FULL_GIT_COMMIT INSTALLED_CONTROL_PLANE_DIGEST EXPECTED_ACTIVE_POLICY_SHA256 EXPECTED_ACTIVE_PROMOTION_SHA256 EXPECTED_AUTHORIZED_FREEZE_MANIFEST EXPECTED_AUTHORIZED_FREEZE_MANIFEST_SHA256 EXPECTED_AUTHORIZED_FREEZE_BUILD_CONTROLLER"
[[ "$#" -eq 7 ]] || {
  printf '%s\n' "$USAGE" >&2
  exit 1
}

EXPECTED_GIT_COMMIT="$1"
CONTROL_PLANE_VERSION="$2"
EXPECTED_ACTIVE_POLICY_SHA256="$3"
EXPECTED_ACTIVE_PROMOTION_SHA256="$4"
EXPECTED_AUTHORIZED_FREEZE_MANIFEST="$5"
EXPECTED_AUTHORIZED_FREEZE_MANIFEST_SHA256="$6"
EXPECTED_AUTHORIZED_FREEZE_BUILD_CONTROLLER="$7"
CONTROL_PLANE_DIR="${PIC_ROOT}/control_plane/${CONTROL_PLANE_VERSION}"
ENV_FILE="${CONTROL_PLANE_DIR}/frontier_pic_environment.sh"
PYTHON=/opt/cray/pe/python/3.11.7/bin/python3
CONTROL_PLANE=("$PYTHON" -I "${CONTROL_PLANE_DIR}/run_control_plane.py")
CONFIG=hip-mpi-release-paper-pic

[[ "$EXPECTED_GIT_COMMIT" =~ ^[0-9a-f]{40}$ ]] || {
  printf 'Expected one full lowercase Git commit, got: %s\n' "$EXPECTED_GIT_COMMIT" >&2
  exit 1
}
[[ "$CONTROL_PLANE_VERSION" =~ ^[0-9a-f]{64}$ ]] || {
  printf 'Expected one installed control-plane digest, got: %s\n' \
    "$CONTROL_PLANE_VERSION" >&2
  exit 1
}
[[ "$EXPECTED_ACTIVE_POLICY_SHA256" =~ ^[0-9a-f]{64}$ ]] || {
  printf 'Malformed expected active-policy digest: %s\n' \
    "$EXPECTED_ACTIVE_POLICY_SHA256" >&2
  exit 1
}
[[ "$EXPECTED_ACTIVE_PROMOTION_SHA256" =~ ^[0-9a-f]{64}$ ]] || {
  printf 'Malformed expected active-promotion digest: %s\n' \
    "$EXPECTED_ACTIVE_PROMOTION_SHA256" >&2
  exit 1
}
EXPECTED_AUTHORIZED_FREEZE_MANIFEST_RE="^${PIC_ROOT}/clean_candidates/[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}/clean_candidate_manifest\\.json$"
[[ "$EXPECTED_AUTHORIZED_FREEZE_MANIFEST" =~ $EXPECTED_AUTHORIZED_FREEZE_MANIFEST_RE ]] || {
  printf 'Malformed expected authorized-freeze manifest path: %s\n' \
    "$EXPECTED_AUTHORIZED_FREEZE_MANIFEST" >&2
  exit 1
}
[[ "$EXPECTED_AUTHORIZED_FREEZE_MANIFEST_SHA256" =~ ^[0-9a-f]{64}$ ]] || {
  printf 'Malformed expected authorized-freeze manifest digest: %s\n' \
    "$EXPECTED_AUTHORIZED_FREEZE_MANIFEST_SHA256" >&2
  exit 1
}
[[ "$EXPECTED_AUTHORIZED_FREEZE_BUILD_CONTROLLER" =~ ^[0-9a-f]{64}$ ]] || {
  printf 'Malformed expected authorized-freeze build controller: %s\n' \
    "$EXPECTED_AUTHORIZED_FREEZE_BUILD_CONTROLLER" >&2
  exit 1
}

VERIFY_ACTIVE_GENERATION=(
  "${CONTROL_PLANE[@]}" promote_active_policy.py
  --verify-active-launch-prohibited-generation
  --expected-control-plane-version "$CONTROL_PLANE_VERSION"
  --expected-active-policy-sha256 "$EXPECTED_ACTIVE_POLICY_SHA256"
  --expected-active-promotion-sha256 "$EXPECTED_ACTIVE_PROMOTION_SHA256"
  --expected-authorized-freeze-manifest "$EXPECTED_AUTHORIZED_FREEZE_MANIFEST"
  --expected-authorized-freeze-manifest-sha256 \
    "$EXPECTED_AUTHORIZED_FREEZE_MANIFEST_SHA256"
  --expected-authorized-freeze-build-controller \
    "$EXPECTED_AUTHORIZED_FREEZE_BUILD_CONTROLLER"
)

cd "$SRC_DIR"
SOURCE_STATUS=$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null \
    status --porcelain --untracked-files=all
)
test -z "$SOURCE_STATUS"
test "$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null rev-parse HEAD
)" = "$EXPECTED_GIT_COMMIT"
test "$(
  /usr/bin/env -i HOME=/ LANG=C LC_ALL=C PATH=/usr/bin:/bin \
    /usr/bin/git -c core.fsmonitor=false -c core.hooksPath=/dev/null \
    rev-parse origin/PIC
)" = "$EXPECTED_GIT_COMMIT"

"${CONTROL_PLANE[@]}" validate_and_reserve_frontier_job.py verify-control-plane >/dev/null
export HOME=/
export USER="$(/usr/bin/id -un)"
source /opt/cray/pe/lmod/lmod/init/profile
unset BASH_ENV ENV
source "$ENV_FILE" || exit $?

"${VERIFY_ACTIVE_GENERATION[@]}" >/dev/null
printf 'source_commit=%s\n' "$EXPECTED_GIT_COMMIT"
printf 'control_plane_version=%s\n' "$CONTROL_PLANE_VERSION"
printf 'expected_active_policy_sha256=%s\n' "$EXPECTED_ACTIVE_POLICY_SHA256"
printf 'expected_active_promotion_sha256=%s\n' "$EXPECTED_ACTIVE_PROMOTION_SHA256"
printf 'expected_authorized_freeze_manifest=%s\n' \
  "$EXPECTED_AUTHORIZED_FREEZE_MANIFEST"
printf 'expected_authorized_freeze_manifest_sha256=%s\n' \
  "$EXPECTED_AUTHORIZED_FREEZE_MANIFEST_SHA256"
printf 'expected_authorized_freeze_build_controller=%s\n' \
  "$EXPECTED_AUTHORIZED_FREEZE_BUILD_CONTROLLER"

"${CONTROL_PLANE[@]}" write_orion_build_profile.py \
  --source-root "$SRC_DIR" \
  --expected-git-commit "$EXPECTED_GIT_COMMIT" \
  --profile-id "$CONFIG"

COMMIT12="${EXPECTED_GIT_COMMIT:0:12}"
BIN_DIR="${PIC_ROOT}/bin/${COMMIT12}/${CONFIG}"

FREEZE_BINDING=$("${CONTROL_PLANE[@]}" create_clean_candidate_freeze.py \
  --source-root "$SRC_DIR" \
  --executable "${BIN_DIR}/athena" \
  --build-profile "${BIN_DIR}/build_profile.json" \
  --build-profile-id "$CONFIG" \
  --prepared-artifact-inventory \
    tst/publication/frontier_control_plane/prepared_pic_artifact_inventory.json)

read -r CLEAN_CANDIDATE_MANIFEST CLEAN_CANDIDATE_MANIFEST_SHA256 EXTRA \
  <<< "$FREEZE_BINDING"
[[ -n "$CLEAN_CANDIDATE_MANIFEST" && -z "${EXTRA:-}" ]] || {
  printf 'Malformed clean-candidate freeze binding: %s\n' "$FREEZE_BINDING" >&2
  exit 1
}
[[ "$CLEAN_CANDIDATE_MANIFEST_SHA256" =~ ^[0-9a-f]{64}$ ]] || {
  printf 'Malformed clean-candidate manifest digest: %s\n' \
    "$CLEAN_CANDIDATE_MANIFEST_SHA256" >&2
  exit 1
}

printf 'clean_candidate_manifest=%s\n' "$CLEAN_CANDIDATE_MANIFEST"
printf 'clean_candidate_manifest_sha256=%s\n' "$CLEAN_CANDIDATE_MANIFEST_SHA256"

"${CONTROL_PLANE[@]}" revalidate_clean_candidate.py \
  --manifest "$CLEAN_CANDIDATE_MANIFEST" \
  --expected-manifest-sha256 "$CLEAN_CANDIDATE_MANIFEST_SHA256" \
  --expected-git-commit "$EXPECTED_GIT_COMMIT" \
  --expected-receipt-control-plane-version "$CONTROL_PLANE_VERSION"

"${VERIFY_ACTIVE_GENERATION[@]}" >/dev/null
