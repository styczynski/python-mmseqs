ALL_SOURCES := $(shell find ./src -iname *.h -o -iname *.cpp -o -iname *.hpp -o -iname *.c)

lint:
	poetry run black .
	poetry run isort .
	poetry run flakehell lint

install:
	poetry install

update:
	poetry update

test:
	mkdir -p test_reports && poetry run pytest -n 4 --ignore lib --junitxml=test_reports/latest.xml --cov=biosnake --cov-report html && mv htmlcov test_reports/coverage_latest_html

publish:
	poetry build && poetry run s3pypi --bucket pypi.covidgenomics.com --private --region eu-west-1 --dist-path dist

format: lint
	docker run -it -v "$(CURDIR)":/workdir -w /workdir unibeautify/clang-format -style=Google -i $(ALL_SOURCES)

docs:
	rm -rfd docs && cd .doc && rm -rfd _build && poetry run make html && mv _build/html ../docs && cp CNAME ../docs && touch ../docs/.nojekyll

alpine-build:
	mkdir -p alpine_dist > /dev/null 2> /dev/null
	docker build -t covid-genomics/biosnake-alpine-builder:1.0 .
	docker run -it -v "$(CURDIR)/alpine_dist":/usr/src/biosnake/dist covid-genomics/biosnake-alpine-builder:1.0

.PHONY: docs
