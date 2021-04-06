#! /bin/bash

set -e
XLSX="$HOME/OneDrive - Stanford/COVDB-DatabaseFiles/SARS2MutationComments.xlsx"
CSV="/tmp/sars2-mutation-comments.csv"
cd $(dirname $0)/..

pipenv run xlsx2csv "${XLSX}" "${CSV}" -n "ForExporting (don't edit)" -f %Y-%m-%d
# this probably shouldn't belong here, it's just a trade-off since the mutation comments
# will evolve so often and we don't want to release a big Sierra version everyday
pipenv run python scripts/mutcmt2json.py /tmp/sars2-mutation-comments.csv downloads/mutation-comments/$(date +'%Y%m%d').json $(date +'%Y%m%d')
rm -rf downloads/mutation-comments/latest.json 2>/dev/null || true
ln -vs $(date +'%Y%m%d').json downloads/mutation-comments/latest.json
