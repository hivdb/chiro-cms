#! /bin/bash

set -e
XLSXDIR="$HOME/OneDrive - Stanford/COVDB-DatabaseFiles"
CSVDIR="$HOME/MAbSummary-CSVFILES"

mkdir -p "${CSVDIR}"

cd $(dirname $0)/..

pipenv run xlsx2csv "${XLSXDIR}/MAbSummary.xlsx" "${CSVDIR}/MAbs.csv" -n "__pending_removal_MAbs"
pipenv run xlsx2csv "${XLSXDIR}/MAbSummary.xlsx" "${CSVDIR}/Structures.csv" -n "Structure meta"
pipenv run xlsx2csv "${XLSXDIR}/MAbSummary.xlsx" "${CSVDIR}/Sequences.csv" -n "Sequences"


pipenv run python3 scripts/check-structure-sequence.py ${CSVDIR}
