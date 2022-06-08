#! /bin/bash

set -e

BRANCH=`git rev-parse --abbrev-ref HEAD`

if [[ "$BRANCH" != "main" ]]; then
    echo "Please switch to branch 'main' for executing '$@'." >&2
    exit 1
fi
