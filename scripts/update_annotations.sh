#! /bin/bash

set -e

cd $(dirname $0)/..

pipenv run python ./scripts/update_annotations.py scripts/update_annotations_config.yml
pipenv run python ./scripts/clean_annotations.py resources/mutannot-spike.yml
pipenv run python ./scripts/clean_annotations.py resources/mutannot-rdrp.yml
