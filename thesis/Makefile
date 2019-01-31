ifeq (,$(shell sh -c 'cygpath --version 2> /dev/null'))
  # Unix
  program = okular
else
  # Windows mit MSys2/Cygwin
  program = sumatraPDF
endif

all: build/thesis.pdf


TeXOptions = -pdf \
			 -interaction=nonstopmode \
			 -halt-on-error \
			 -output-directory=build \
			 -shell-escape

build/thesis.pdf: FORCE | build
	latexmk $(TeXOptions) thesis.tex

FORCE:

build:
	mkdir -p build/

clean:
	rm -rf build

once: build
	pdflatex --output-directory=build --shell-escape thesis.tex

see:
	$(program) build/thesis.pdf
