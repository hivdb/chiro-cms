requirements.txt: Pipfile.lock
	@pipenv lock --requirements > requirements.txt

build-docker: requirements.txt
	@docker pull python:3
	@docker build . --no-cache -t hivdb/chiro-cms-builder:latest

push-docker: build-docker
	@docker push hivdb/chiro-cms-builder:latest

pull-docker:
	@docker pull hivdb/chiro-cms-builder:latest

_docker-build:
	@rm -rf build/
	@docker run \
		--mount type=bind,source=$(PWD),target=/app \
		--workdir /app --rm -it \
		hivdb/chiro-cms-builder:latest \
		python build.py

_fast-build:
	@rm -rf build/
	@pipenv run python build.py

local/timestamp/%:
	@mkdir -p $$(dirname $@)
	@rm $$(dirname $@)/* 2>/dev/null || true
	touch $@

resources/sierra-sars2/outbreak.info/lineages.json: % : local/timestamp/%.$(shell date +"%Y%m%d")
	@mkdir -p local/
	curl -sSL -H 'authorization: Bearer 0ed52bbfb6c79d1fd8e9c6f267f9b6311c885a4c4c6f037d6ab7b3a40d586ad0' -o "local/lineages.json" "https://api.outbreak.info/genomics/lineage?name=*"
	python3.9 -m "json.tool" --indent=2 local/lineages.json > resources/sierra-sars2/outbreak.info/lineages.json
	@test "$$(jq .success resources/sierra-sars2/outbreak.info/lineages.json)" = "true" || (echo "success != true" && false)

build: $(shell find . -type f -not -path "./.git*" -a -not -path "*.swp" -a -not -path "*.swo" -a -not -path "*/.DS_Store" -a -not -path "*/.gradle/*" -a -not -path "*/build/*" -a -not -path "*/build_gz/*" -a -not -path "*.log" -a -not -path "*/local/*" | sed 's#\([| ]\)#\\\1#g') build.py build_plugins/*.py resources/sierra-sars2/outbreak.info/lineages.json
	@test -e $(shell which pipenv) && make _fast-build || make _docker-build

downloads/resistance-mutations/latest.json: downloads/resistance-mutations/3cl.tsv downloads/resistance-mutations/rdrp.tsv downloads/resistance-mutations/latest.tsv scripts/sars2_drms_to_json.py
	@scripts/sars2_drms_to_json.py

build_gz: build
	@scripts/build_gz.sh

deploy-dev: build_gz
	@docker run \
		--mount type=bind,source=$(HOME)/.aws,target=/root/.aws,readonly \
		--mount type=bind,source=$(PWD),target=/app,readonly \
		--rm -it hivdb/chiro-cms-builder:latest \
		aws s3 sync /app/build_gz s3://cms.hivdb.org/chiro-dev --delete --content-encoding gzip

deploy-dev2: build_gz
	@docker run \
		--mount type=bind,source=$(HOME)/.aws,target=/root/.aws,readonly \
		--mount type=bind,source=$(PWD),target=/app,readonly \
		--rm -it hivdb/chiro-cms-builder:latest \
		aws s3 sync /app/build_gz s3://cms.hivdb.org/chiro-dev2 --delete --content-encoding gzip

deploy-prod: build_gz
	@scripts/only-main.sh make deploy-prod
	@docker run \
		--mount type=bind,source=$(HOME)/.aws,target=/root/.aws,readonly \
		--mount type=bind,source=$(PWD),target=/app,readonly \
		--rm -it hivdb/chiro-cms-builder:latest \
		aws s3 sync /app/build_gz s3://cms.hivdb.org/chiro-prod \
		--delete --content-encoding gzip --cache-control max-age=600

deploy-all: deploy-dev deploy-dev2 deploy-prod

.PHONY: deploy-* _*
