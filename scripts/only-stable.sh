#! /bin/bash

set -e

BRANCH=`git rev-parse --abbrev-ref HEAD`

if [[ "$BRANCH" != "stable" ]]; then
    echo "Please switch to branch 'stable' for executing '$@'." >&2
    exit 1
fi
