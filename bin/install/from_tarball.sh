#!/usr/bin/env bash

mkdir -p src/
pushd src/

# Fastcodeml 
wget \
    --continue \
    --output-document FastCodeML-1.1.0.tar.gz \
    ftp://ftp.vital-it.ch/tools/FastCodeML/FastCodeML-1.1.0.tar.gz
tar xvf FastCodeML-1.1.0.tar.gz
cp FastCodeML-1.1.0/fast ../bin/

popd