.PHONY: all
all: fbm.pdf

.PHONY: clean
clean:
	latexrun -O .tex --clean-all

.PHONY: distclean
distclean: clean
	rm -f fbm.tex

.PHONY: force
fbm.pdf: fbm.tex force
fbm.tex: fbm.nw

.SUFFIXES: .nw .tex .pdf
.nw.tex:
	noweave -delay -index $< > $@
.tex.pdf:
	latexrun -O .tex -o $@ $<
