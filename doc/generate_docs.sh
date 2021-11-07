#!/usr/bin/env bash

# Author: Shashwat Sharma
# Created on Nov. 07, 2021

# Prerequisites:
# latex, latexmk, sphinx

make clean
make html
make latexpdf
cp ./build/latex/Strata.pdf ./
ln -sf ./build/html/contents.html ./Strata.html
