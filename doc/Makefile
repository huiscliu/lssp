
doc = lssp.pdf

SRC_TEX   = $(wildcard *.tex)
TEXFILES = $(SRC_TEX)

LATEX = latex
DVIPDF = dvipdf

default: $(doc)

lssp.pdf: lssp.dvi
	$(DVIPDF) lssp.dvi

lssp.dvi: $(TEXFILES)
	$(LATEX) lssp.tex
	$(LATEX) lssp.tex

.PHONY: clean

clean:
	@rm -fv *.aux lssp.dvi *.log *.out *.toc lssp.idx lssp.ind lssp.ilg
