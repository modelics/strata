#!/usr/bin/env bash

# Author: Shashwat Sharma
# Created on Nov. 07, 2021

# Prerequisites:
# latex, latexmk, sphinx

sphinx-build -b html ./source ./build/html
sphinx-build -b latex ./source ./build/latex
make -C ./build/latex
cp ./build/latex/strata.pdf ./
ln -sf ./build/html/index.html ./strata.html


