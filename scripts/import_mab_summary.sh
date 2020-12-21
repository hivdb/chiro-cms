#! /bin/bash

set -e
XLSXDIR="$HOME/OneDrive - Stanford/COVDB-DatabaseFiles"
CSVDIR="$HOME/MAbSummary-CSVFILES"
cd $(dirname $0)/..

if [[ "$1" != "--no-update-refid-lookup" ]]; then
    echo "Updating RefID lookup table (resources/refid_lookup.json)..."
    pipenv run xlsx2csv "${XLSXDIR}/References.xlsx" "${CSVDIR}/articles.csv" -n tblReferences -f %Y-%m-%d
    pipenv run python3 scripts/build_refid_doi_lookup.py "${CSVDIR}/articles.csv" "resources/refid_lookup.json"
fi

pipenv run xlsx2csv "${XLSXDIR}/MAbSummary.xlsx" "${CSVDIR}/MAbs.csv" -n MAbs -f %Y-%m-%d
pipenv run xlsx2csv "${XLSXDIR}/MAbSummary.xlsx" "${CSVDIR}/MAbStructures.csv" -n "Structure meta"

rm resources/mab-table-data.yml 2>/dev/null || true

pipenv run python3 scripts/import_mab_summary.py \
    ${CSVDIR}/MAbs.csv \
    --refid-lookup resources/refid_lookup.json \
    resources/mabs-table.yml

pipenv run python3 scripts/import_epitopes.py \
     ${CSVDIR}/MAbStructures.csv \
     resources/mutannots/spike

touch resources/
