html: ntt32 prt
	echo "HTML build complete."

pdf: ntt32-pdf prt-pdf
	echo "PDF build complete."

ntt32: ntt32.md ntt.yaml style.css
	pandoc ntt32.md --metadata-file=ntt.yaml --mathjax --standalone --number-sections --css=style.css --syntax-highlighting=tango -f markdown -t html -o ntt32.html

prt: prt.md poly-math.yaml style.css
	pandoc prt.md --metadata-file=poly-math.yaml --mathjax --standalone --number-sections --css=style.css --syntax-highlighting=tango -f markdown -t html -o prt.html

ntt32-pdf: ntt32.md ntt.yaml
	pandoc ntt32.md --log=pandoc.log --pdf-engine=xelatex --metadata-file=ntt.yaml --syntax-highlighting=tango -o ntt32.pdf

prt-pdf: prt.md poly-math.yaml
	pandoc prt.md --log=pandoc.log --pdf-engine=xelatex --metadata-file=poly-math.yaml --syntax-highlighting=tango -o prt.pdf

clean:
	rm -f ntt32.html
	rm -f ntt32.pdf
	rm -f *.aux
	rm -f *.log
	rm -f *.out
