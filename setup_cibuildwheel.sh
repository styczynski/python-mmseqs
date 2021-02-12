#!/bin/bash

echo "Install Poetry..."
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python

echo "Prepare cibuildwheel global dependencies..."
poetry export -f requirements.txt --output requirements.txt
pip install -r requirements.txt
