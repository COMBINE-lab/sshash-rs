#!/usr/bin/env bash
# release-py.sh — bump sshash-py version, commit, tag, and optionally push.
#
# Usage:
#   ./release-py.sh <version>          # bump, commit, tag (no push)
#   ./release-py.sh <version> --push   # bump, commit, tag, then push
#
# Example:
#   ./release-py.sh 0.2.0 --push

set -euo pipefail

# ── helpers ──────────────────────────────────────────────────────────────────

die() { echo "error: $*" >&2; exit 1; }

usage() {
    echo "usage: $0 <version> [--push]"
    echo "  version  semver string, e.g. 1.2.3"
    exit 1
}

# ── args ─────────────────────────────────────────────────────────────────────

[[ $# -lt 1 ]] && usage

VERSION="$1"
PUSH=false
[[ "${2-}" == "--push" ]] && PUSH=true

# Basic semver sanity check
if ! [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    die "version must be in X.Y.Z format, got '$VERSION'"
fi

TAG="sshash-py-v${VERSION}"

CARGO_TOML="crates/sshash-py/Cargo.toml"
PYPROJECT="crates/sshash-py/pyproject.toml"

# ── pre-flight checks ─────────────────────────────────────────────────────────

cd "$(dirname "$0")"   # run from repo root regardless of cwd

[[ -f "$CARGO_TOML" ]]  || die "not found: $CARGO_TOML"
[[ -f "$PYPROJECT" ]]   || die "not found: $PYPROJECT"

if git rev-parse "$TAG" &>/dev/null; then
    die "tag $TAG already exists"
fi

# Warn about uncommitted changes (but allow staged changes from this script)
if ! git diff --quiet; then
    die "working tree has unstaged changes — commit or stash them first"
fi

# ── current versions (for display) ───────────────────────────────────────────

OLD_CARGO=$(grep '^version' "$CARGO_TOML" | head -1 | sed 's/.*"\(.*\)"/\1/')
OLD_PY=$(grep '^version' "$PYPROJECT"    | head -1 | sed 's/.*"\(.*\)"/\1/')

echo "sshash-py current versions:"
echo "  Cargo.toml  : $OLD_CARGO"
echo "  pyproject.toml: $OLD_PY"
echo ""
echo "Bumping to $VERSION  (tag: $TAG)"
echo ""

# ── update version strings ───────────────────────────────────────────────────

# Cargo.toml: replace the first `version = "..."` line
sed -i.bak "1,/^version = /s/^version = \"[^\"]*\"/version = \"${VERSION}\"/" "$CARGO_TOML"
rm "$CARGO_TOML.bak"

# pyproject.toml: same approach
sed -i.bak "1,/^version = /s/^version = \"[^\"]*\"/version = \"${VERSION}\"/" "$PYPROJECT"
rm "$PYPROJECT.bak"

# Refresh Cargo.lock so it reflects the new version
cargo check -p sshash-py -q

# ── verify the edits look right ───────────────────────────────────────────────

NEW_CARGO=$(grep '^version' "$CARGO_TOML" | head -1 | sed 's/.*"\(.*\)"/\1/')
NEW_PY=$(grep '^version' "$PYPROJECT"    | head -1 | sed 's/.*"\(.*\)"/\1/')

[[ "$NEW_CARGO" == "$VERSION" ]] || die "Cargo.toml update failed (got $NEW_CARGO)"
[[ "$NEW_PY"    == "$VERSION" ]] || die "pyproject.toml update failed (got $NEW_PY)"

echo "Updated:"
echo "  Cargo.toml    : $OLD_CARGO → $NEW_CARGO"
echo "  pyproject.toml: $OLD_PY → $NEW_PY"
echo ""

# ── commit and tag ────────────────────────────────────────────────────────────

git add "$CARGO_TOML" "$PYPROJECT" Cargo.lock
git commit -m "chore(sshash-py): release v${VERSION}"
git tag -a "$TAG" -m "sshash-py v${VERSION}"

echo "Created commit and tag $TAG"

# ── optional push ─────────────────────────────────────────────────────────────

if $PUSH; then
    echo ""
    git push origin HEAD
    git push origin "$TAG"
    echo "Pushed commit and tag to origin."
else
    echo ""
    echo "To push when ready:"
    echo "  git push origin HEAD && git push origin $TAG"
fi
