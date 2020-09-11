requirements.txt: Pipfile.lock
	@pipenv freeze --requirements > requirements.txt

build-docker:
	@docker pull python:3
	@docker build . --no-cache -t hivdb/chiro-cms-builder:latest

push-docker: build-docker
	@docker push hivdb/chiro-cms-builder:latest

pull-docker:
	@docker pull hivdb/chiro-cms-builder:latest

build: build.py pages images resources downloads
	@rm -rf build/
	@docker run \
		--mount type=bind,source=$(PWD),target=/app \
		--workdir /app --rm -it \
		hivdb/chiro-cms-builder:latest \
		python build.py

deploy-dev: build
	@docker run \
		--mount type=bind,source=$(HOME)/.aws,target=/root/.aws,readonly \
		--mount type=bind,source=$(PWD),target=/app,readonly \
		--rm -it hivdb/chiro-cms-builder:latest \
		aws s3 sync /app/build s3://cms.hivdb.org/chiro-dev --delete

deploy-prod: build
	@docker run \
		--mount type=bind,source=$(HOME)/.aws,target=/root/.aws,readonly \
		--mount type=bind,source=$(PWD),target=/app,readonly \
		--rm -it hivdb/chiro-cms-builder:latest \
		aws s3 sync /app/build s3://cms.hivdb.org/chiro-prod \
		--delete --cache-control max-age=600

deploy-all: deploy-dev deploy-prod

.PHONY: deploy-*
