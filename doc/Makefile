.SECONDARY:

LATEX = pdflatex

doc: mock_zeldovich.pdf linear_zeldovich.pdf

%.pdf: %.tex
	$(LATEX) $<
	@grep -e 'may have changed' $(basename $<).log && rm $@; $(MAKE) $@ || echo "done."


.PHONY: clean check

clean:
	rm -f *.bbl *.blg *.bak *.dvi *.log *.aux

check:
	aspell -c --lang=en_GB linear_zeldovich.tex
