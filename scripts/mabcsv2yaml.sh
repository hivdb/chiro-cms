#! /bin/bash

set -e
XLSXDIR="$HOME/OneDrive - Stanford/COVDB-DatabaseFiles"
CSVDIR="$HOME/MAbSummary-CSVFILES"
cd $(dirname $0)/..

if [[ "$1" == "--update-refid-lookup" ]]; then
    echo "Updating RefID lookup table (${CSVDIR}/articles.json)..."
    pipenv run xlsx2csv "${XLSXDIR}/References.xlsx" "${CSVDIR}/articles.csv" -n tblReferences -f %Y-%m-%d
    pipenv run python3 scripts/build_refid_doi_lookup.py "${CSVDIR}/articles.csv" "${CSVDIR}/articles.json"
fi

pipenv run xlsx2csv "${XLSXDIR}/MAbSummary.xlsx" "${CSVDIR}/MAbs.csv" -n MAbs -f %Y-%m-%d

rm resources/mab-table-data.yml 2>/dev/null || true

pipenv run python3 scripts/mabcsv2yaml.py \
    ${CSVDIR}/MAbs.csv \
    --refid-lookup ${CSVDIR}/articles.json \
    resources/mabs-table-data.yml

touch resources/
