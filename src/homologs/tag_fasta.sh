#!/usr/bin/env bash

set -euo pipefail
TAG=$1

sed "s/^>/>$TAG/"