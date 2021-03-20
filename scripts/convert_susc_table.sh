#!/bin/sh

cp ../covid-drdb-reports/susceptibility-data_files/table/*.json ./resources/susc-table/
pipenv run python ./scripts/convert_susc_table.py ./resources/susc-table/table3.json ./resources/susc-table/table3.json
pipenv run python ./scripts/convert_susc_table.py ./resources/susc-table/table4.json ./resources/susc-table/table4.json
pipenv run python ./scripts/convert_susc_table.py ./resources/susc-table/tableS3.json ./resources/susc-table/tableS3.json
pipenv run python ./scripts/convert_susc_table.py ./resources/susc-table/tableS4.json ./resources/susc-table/tableS4.json
pipenv run python ./scripts/convert_susc_table.py ./resources/susc-table/tableS5.json ./resources/susc-table/tableS5.json
pipenv run python ./scripts/convert_susc_table.py ./resources/susc-table/tableS6.json ./resources/susc-table/tableS6.json
pipenv run python ./scripts/convert_susc_table.py ./resources/susc-table/tableS7.json ./resources/susc-table/tableS7.json


pipenv run python ./scripts/json2yaml.py ./resources/susc-table/table3.json resources/susc-table/table3.yml
pipenv run python ./scripts/json2yaml.py ./resources/susc-table/table4.json resources/susc-table/table4.yml
pipenv run python ./scripts/json2yaml.py ./resources/susc-table/tableS3.json resources/susc-table/tableS3.yml
pipenv run python ./scripts/json2yaml.py ./resources/susc-table/tableS4.json resources/susc-table/tableS4.yml
pipenv run python ./scripts/json2yaml.py ./resources/susc-table/tableS5.json resources/susc-table/tableS5.yml
pipenv run python ./scripts/json2yaml.py ./resources/susc-table/tableS6.json resources/susc-table/tableS6.yml
pipenv run python ./scripts/json2yaml.py ./resources/susc-table/tableS7.json resources/susc-table/tableS7.yml
