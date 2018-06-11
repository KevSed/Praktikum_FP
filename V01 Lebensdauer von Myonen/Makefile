all: plots latex

latex: FORCE | build
	@lualatex --halt-on-error --interaction=nonstopmode --output-directory=build --draftmode main.tex
	@biber build/main.bcf
	@lualatex --halt-on-error --interaction=nonstopmode --output-directory=build main.tex

FORCE:

build:
	@mkdir -p build/

clean:
	@rm -rf build/

plots:
	@python analysis/verzoegerung.py
	@python analysis/kalibrierung.py
	@python analysis/spektrum.py
