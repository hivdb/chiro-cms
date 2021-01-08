#! /bin/bash

set -e

cd $(dirname $0)/..

pipenv run python scripts/import_gluepcnt.py RdRP S
