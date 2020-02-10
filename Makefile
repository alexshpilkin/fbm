.PHONY: all clean
.SUFFIXES: .nw .tex

all: fbm.tex

clean:
	rm -f fbm.tex

.nw.tex:
	noweave -delay -index $< > $@
