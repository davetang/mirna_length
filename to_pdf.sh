#!/bin/bash

file=`basename $1 .tex`

latex $file
bibtex $file
latex $file
pdflatex $file

rm -rf $file.blg $file.bbl $file.dvi $file.aux $file.log $file.out
