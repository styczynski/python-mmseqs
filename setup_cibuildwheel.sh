#!/bin/bash

echo "Install Poetry..."
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python
export PATH="$PATH:$HOME/.poetry/bin"
source $HOME/.poetry/env
echo "Poetry installed."

echo "Export Poetry requirements..."
poetry export -f requirements.txt --output requirements.txt
echo "Install everything..."
pip install -r requirements.txt
