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

pipenv run xlsx2csv "${XLSXDIR}/MAbSummary.xlsx" "${CSVDIR}/publishedMAbs.csv" -n publishedMAbs -f %Y-%m-%d
# pipenv run xlsx2csv "${XLSXDIR}/MAbSummary.xlsx" "${CSVDIR}/sdAbs.csv" -n sdAbs -f %Y-%m-%d
pipenv run xlsx2csv "${XLSXDIR}/MAbSummary.xlsx" "${CSVDIR}/animalModelMAbs.csv" -n animalModelMAbs -f %Y-%m-%d
pipenv run xlsx2csv "${XLSXDIR}/MAbSummary.xlsx" "${CSVDIR}/clinicalTrialMAbs.csv" -n clinicalTrialMAbs -f %Y-%m-%d
pipenv run xlsx2csv "${XLSXDIR}/MAbSummary.xlsx" "${CSVDIR}/sequences.csv" -n Sequences -f %Y-%m-%d

rm resources/mab-table-data.yml 2>/dev/null || true

pipenv run python3 scripts/mabcsv2yaml.py \
    ${CSVDIR}/publishedMAbs.csv \
    --alignment ${CSVDIR}/sequences.csv \
    --refid-lookup ${CSVDIR}/articles.json \
    resources/published-mabs.yml

pipenv run python3 scripts/mabcsv2yaml.py \
    ${CSVDIR}/animalModelMAbs.csv \
    --refid-lookup ${CSVDIR}/articles.json \
    resources/animal-model-mabs.yml

pipenv run python3 scripts/mabcsv2yaml.py \
    ${CSVDIR}/clinicalTrialMAbs.csv \
    --refid-lookup ${CSVDIR}/articles.json \
    resources/clinical-trial-mabs.yml

touch resources/
