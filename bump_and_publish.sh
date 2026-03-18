#!/usr/bin/env bash
set -euo pipefail

die() {
    echo "error: $*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage:
  ./bump_and_publish.sh <version> [--publish] [--dry-run]
  ./bump_and_publish.sh [--publish] [--dry-run] <version>

Options:
  --publish  Publish sshash-lib first, then sshash, after bumping and committing
  --dry-run  Show what would be done without modifying files, creating commits, tags, or publishing
  -h, --help Show this help message
EOF
}

print_cmd() {
    printf '+'
    printf ' %q' "$@"
    printf '\n'
}

run() {
    print_cmd "$@"
    if [[ "$DRY_RUN" == true ]]; then
        return 0
    fi
    "$@"
}

VERSION=""
PUBLISH=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --publish)
            PUBLISH=true
            ;;
        --dry-run)
            DRY_RUN=true
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            die "unknown option: $1"
            ;;
        *)
            if [[ -n "$VERSION" ]]; then
                die "version specified more than once"
            fi
            VERSION="$1"
            ;;
    esac
    shift
done

[[ -n "$VERSION" ]] || {
    usage
    exit 1
}

if ! [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+([+-][0-9A-Za-z.-]+)*$ ]]; then
    die "version must look like X.Y.Z, optionally with prerelease/build suffixes"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ROOT_CARGO="Cargo.toml"
LOCKFILE="Cargo.lock"
TAG="v${VERSION}"
LIB_CRATE="sshash-lib"
CLI_CRATE="sshash"

[[ -f "$ROOT_CARGO" ]] || die "not found: $ROOT_CARGO"
[[ -f "$LOCKFILE" ]] || die "not found: $LOCKFILE"

CURRENT_VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$ROOT_CARGO" | head -1)"
CURRENT_LIB_DEP_VERSION="$(sed -n 's/^sshash-lib = { version = "\(.*\)", path = "crates\/sshash-lib" }/\1/p' "$ROOT_CARGO")"

[[ -n "$CURRENT_VERSION" ]] || die "could not determine current workspace version from $ROOT_CARGO"
[[ -n "$CURRENT_LIB_DEP_VERSION" ]] || die "could not determine sshash-lib workspace dependency version from $ROOT_CARGO"

if [[ "$CURRENT_VERSION" == "$VERSION" ]]; then
    die "workspace version is already $VERSION"
fi

if git rev-parse "$TAG" >/dev/null 2>&1; then
    die "tag $TAG already exists"
fi

if [[ -n "$(git status --porcelain)" ]]; then
    die "working tree is not clean; commit or stash existing changes first"
fi

echo "Current workspace version : $CURRENT_VERSION"
echo "Current sshash-lib dep    : $CURRENT_LIB_DEP_VERSION"
echo "New workspace version     : $VERSION"
echo "Tag                       : $TAG"
if [[ "$PUBLISH" == true ]]; then
    echo "Publish crates            : yes"
else
    echo "Publish crates            : no"
fi
if [[ "$DRY_RUN" == true ]]; then
    echo "Dry-run                   : yes"
else
    echo "Dry-run                   : no"
fi
echo

echo "Updating $ROOT_CARGO"
echo "  [workspace.package] version: $CURRENT_VERSION -> $VERSION"
echo "  [workspace.dependencies].sshash-lib version: $CURRENT_LIB_DEP_VERSION -> $VERSION"

if [[ "$DRY_RUN" == false ]]; then
    sed -i.bak "1,/^version = /s/^version = \".*\"/version = \"${VERSION}\"/" "$ROOT_CARGO"
    rm -f "${ROOT_CARGO}.bak"

    sed -i.bak "s/^sshash-lib = { version = \".*\", path = \"crates\\/sshash-lib\" }/sshash-lib = { version = \"${VERSION}\", path = \"crates\\/sshash-lib\" }/" "$ROOT_CARGO"
    rm -f "${ROOT_CARGO}.bak"
fi

UPDATED_VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$ROOT_CARGO" | head -1)"
UPDATED_LIB_DEP_VERSION="$(sed -n 's/^sshash-lib = { version = "\(.*\)", path = "crates\/sshash-lib" }/\1/p' "$ROOT_CARGO")"

if [[ "$DRY_RUN" == false ]]; then
    [[ "$UPDATED_VERSION" == "$VERSION" ]] || die "workspace version update failed"
    [[ "$UPDATED_LIB_DEP_VERSION" == "$VERSION" ]] || die "sshash-lib workspace dependency update failed"
else
    echo "Dry-run: would rewrite $ROOT_CARGO and refresh $LOCKFILE"
fi

run cargo check -p "$LIB_CRATE" -p "$CLI_CRATE" -q
run git add "$ROOT_CARGO" "$LOCKFILE"
run git commit -m "chore(release): bump Rust crates to v${VERSION}"

if [[ "$PUBLISH" == true ]]; then
    run cargo publish -p "$LIB_CRATE"

    if [[ "$DRY_RUN" == true ]]; then
        echo "Dry-run: would wait for ${LIB_CRATE} v${VERSION} to appear on crates.io before publishing ${CLI_CRATE}"
    else
        echo "Waiting for ${LIB_CRATE} v${VERSION} to appear on crates.io..."
        MAX_ATTEMPTS=30
        ATTEMPT=0
        while (( ATTEMPT < MAX_ATTEMPTS )); do
            if cargo search "$LIB_CRATE" 2>/dev/null | grep -q "^${LIB_CRATE} = \"${VERSION}\""; then
                echo "${LIB_CRATE} v${VERSION} is available"
                break
            fi
            ATTEMPT=$((ATTEMPT + 1))
            echo "  attempt ${ATTEMPT}/${MAX_ATTEMPTS}: not indexed yet, waiting 10s..."
            sleep 10
        done

        if (( ATTEMPT >= MAX_ATTEMPTS )); then
            die "timed out waiting for ${LIB_CRATE} v${VERSION} on crates.io"
        fi
    fi

    run cargo publish -p "$CLI_CRATE"
fi

run git tag -a "$TAG" -m "Release ${VERSION}"
run git push origin HEAD
run git push origin "$TAG"

if [[ "$DRY_RUN" == true ]]; then
    echo
    echo "Dry-run complete"
else
    echo
    echo "Release bump complete for v${VERSION}"
fi
