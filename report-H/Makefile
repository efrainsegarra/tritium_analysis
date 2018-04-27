.PHONY: all clean

all: short_report.pdf

short_report.pdf: short_report.tex
	pdflatex $^

clean:
	rm -rf *.aux *.log short_report.pdf
