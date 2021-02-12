#!/bin/bash

echo "Prepare cibuildwheel global dependencies..."
poetry export -f requirements.txt --output requirements.txt
pip install -r requirements.txt
