#!/bin/bash
set -o errexit

python -m pytest test "$@"
