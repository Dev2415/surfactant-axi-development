#!/usr/bin/env bash
#
# push-high-order-micelle.sh
#
# Commits any pending changes and pushes the high-order-micelle branch to origin.
#
# Usage:
#   ./push-high-order-micelle.sh "commit message here"
#
# If no commit message is given, a default timestamped message is used.

set -euo pipefail

BRANCH="high-order-micelle"
COMMIT_MSG="${1:-Update on $(date '+%Y-%m-%d %H:%M:%S')}"

# Make sure we're inside a git repo
git rev-parse --is-inside-work-tree > /dev/null 2>&1 || {
  echo "Error: not inside a git repository." >&2
  exit 1
}

CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)

if [ "$CURRENT_BRANCH" != "$BRANCH" ]; then
  echo "Switching from '$CURRENT_BRANCH' to '$BRANCH'..."
  git checkout "$BRANCH"
fi

# Stage everything (new, modified, deleted)
git add -A

# Only commit if there's something staged
if ! git diff --cached --quiet; then
  git commit -m "$COMMIT_MSG"
  echo "Committed: $COMMIT_MSG"
else
  echo "No changes to commit."
fi

# Push, setting upstream if this is the first push of the branch
if git rev-parse --abbrev-ref --symbolic-full-name "@{u}" > /dev/null 2>&1; then
  git push
else
  git push -u origin "$BRANCH"
fi

echo "Done. '$BRANCH' pushed to origin."
