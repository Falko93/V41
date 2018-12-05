all: build/main.pdf

fast: build/fastmain.pdf


TeXOptions = -lualatex \
			 -interaction=nonstopmode \
			 -halt-on-error \
			 -output-directory=build
                                                                                
build/main.pdf: FORCE | build
	latexmk $(TeXOptions) main.tex

build/fastmain.pdf: FORCE | build
	lualatex -output-directory=build main.tex 
	
FORCE:

build:
	mkdir -p build/

clean:
	rm -rf build
