#!/usr/bin/env bash
#
# merge-high-order-micelle-to-main.sh
#
# Merges the high-order-micelle branch into main and pushes main to origin.
#
# Usage:
#   ./merge-high-order-micelle-to-main.sh

set -euo pipefail

FEATURE_BRANCH="high-order-micelle"
TARGET_BRANCH="main"

git rev-parse --is-inside-work-tree > /dev/null 2>&1 || {
  echo "Error: not inside a git repository." >&2
  exit 1
}

# Make sure we have the latest refs from origin
git fetch origin

# Update the feature branch first so main gets the latest committed work
git checkout "$FEATURE_BRANCH"
git pull origin "$FEATURE_BRANCH" --ff-only || true

# Switch to main and update it
git checkout "$TARGET_BRANCH"
git pull origin "$TARGET_BRANCH" --ff-only

# Merge the feature branch in
echo "Merging '$FEATURE_BRANCH' into '$TARGET_BRANCH'..."
git merge --no-ff "$FEATURE_BRANCH" -m "Merge $FEATURE_BRANCH into $TARGET_BRANCH"

# Push the result
git push origin "$TARGET_BRANCH"

echo "Done. '$FEATURE_BRANCH' merged into '$TARGET_BRANCH' and pushed."
