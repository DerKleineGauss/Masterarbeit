ifeq (,$(shell sh -c 'cygpath --version 2> /dev/null'))
  # Unix
  program = okular
else
  # Windows mit MSys2/Cygwin
  program = sumatraPDF
endif

all: build/thesis.pdf python/scaling.py


TeXOptions = -pdf \
			 -interaction=nonstopmode \
			 -halt-on-error \
			 -output-directory=build \
			 -shell-escape

python/scaling.py: FORCE
	python python/scaling.py

build/thesis.pdf: FORCE | build
	latexmk $(TeXOptions) thesis.tex

FORCE:

build:
	mkdir -p build/

clean:
	rm -rf build

once: build python/scaling.py
	pdflatex --output-directory=build --shell-escape thesis.tex

see:
	$(program) build/thesis.pdf
