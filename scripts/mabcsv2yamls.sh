#! /bin/bash

set -e

cd $(dirname $0)/..

rm -rf resources/mab-table-data

pipenv run python3 scripts/mabcsv2yamls.py \
    ~/Dropbox/Coronavirus/DatabaseFiles/MabSummarySep10-CSVFILES/neutralizing-monoclonal-antibodies-described-in-published-papers.csv \
    ~/Dropbox/Coronavirus/DatabaseFiles/mab-seq-compare/Mab-aligned-simple.csv \
    resources/mab-table-data
