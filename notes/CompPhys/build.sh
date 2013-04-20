#!/bin/bash
rm -v *.aux *.dvi *.log *.pdf *.ps *.spl *.bbl *.out

F=main

latex $F.tex
bibtex $F.aux
latex $F.tex
latex $F.tex
dvips -o $F.ps $F.dvi
ps2pdf -sOutputFile=CompPhys_BW.pdf \
       -sColorConversionStrategy=Gray -dProcessColorModel=/DeviceGray $F.ps
       
#ps2pdf -sOutputFile=${F}.pdf $F.ps

