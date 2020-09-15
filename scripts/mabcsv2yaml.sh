#! /bin/bash

set -e

cd $(dirname $0)/..

rm resources/mab-table-data.yml 2>/dev/null || true

pipenv run python3 scripts/mabcsv2yaml.py \
    ~/Dropbox/Coronavirus/DatabaseFiles/MabSummary-CSVFILES/publishedMAbs.csv \
    ~/Dropbox/Coronavirus/DatabaseFiles/mab-seq-compare/Mab-aligned-simple.csv \
    resources/published-mabs-data.yml

pipenv run python3 scripts/mabcsv2yaml.py \
    ~/Dropbox/Coronavirus/DatabaseFiles/MabSummary-CSVFILES/animalModelMAbs.csv \
    __none__ \
    resources/animal-model-mabs-data.yml

pipenv run python3 scripts/mabcsv2yaml.py \
    ~/Dropbox/Coronavirus/DatabaseFiles/MabSummary-CSVFILES/clinicalTrialMAbs.csv \
    __none__ \
    resources/clinical-trial-mabs-data.yml
