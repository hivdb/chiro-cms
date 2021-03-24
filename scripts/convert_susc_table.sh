#!/bin/sh

cp ../covid-drdb-reports/susceptibility-data_files/table/*.json ./resources/susc-table/
pipenv run python ./scripts/convert_susc_table.py ./resources/susc-table/table_summary.json ./resources/susc-table/table_summary.json
pipenv run python ./scripts/convert_susc_table.py ./resources/susc-table/table_mab.json ./resources/susc-table/table_mab.json
pipenv run python ./scripts/convert_susc_table.py ./resources/susc-table/table_plasma.json ./resources/susc-table/table_plasma.json
pipenv run python ./scripts/convert_susc_table.py ./resources/susc-table/table_plasma_variant.json ./resources/susc-table/table_plasma_variant.json
pipenv run python ./scripts/convert_susc_table.py ./resources/susc-table/table_plasma_muts.json ./resources/susc-table/table_plasma_muts.json
pipenv run python ./scripts/convert_susc_table.py ./resources/susc-table/table_mab_variant.json ./resources/susc-table/table_mab_variant.json
pipenv run python ./scripts/convert_susc_table.py ./resources/susc-table/table_mab_muts.json ./resources/susc-table/table_mab_muts.json


pipenv run python ./scripts/json2yaml.py ./resources/susc-table/table_summary.json resources/susc-table/table_summary.yml
pipenv run python ./scripts/json2yaml.py ./resources/susc-table/table_mab.json resources/susc-table/table_mab.yml
pipenv run python ./scripts/json2yaml.py ./resources/susc-table/table_plasma.json resources/susc-table/table_plasma.yml
pipenv run python ./scripts/json2yaml.py ./resources/susc-table/table_plasma_variant.json resources/susc-table/table_plasma_variant.yml
pipenv run python ./scripts/json2yaml.py ./resources/susc-table/table_plasma_muts.json resources/susc-table/table_plasma_muts.yml
pipenv run python ./scripts/json2yaml.py ./resources/susc-table/table_mab_variant.json resources/susc-table/table_mab_variant.yml
pipenv run python ./scripts/json2yaml.py ./resources/susc-table/table_mab_muts.json resources/susc-table/table_mab_muts.yml
