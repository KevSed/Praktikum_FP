all: build/main.pdf

build/main.pdf: FORCE | build
	@lualatex --halt-on-error --interaction=nonstopmode --output-directory=build --draftmode main.tex
	@biber build/main.bcf
	@lualatex --halt-on-error --interaction=nonstopmode --output-directory=build main.tex

FORCE:

build:
	mkdir -p build/
 
clean:
	rm -rf build