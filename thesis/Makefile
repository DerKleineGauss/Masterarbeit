ifeq (,$(shell sh -c 'cygpath --version 2> /dev/null'))
  # Unix
	pwd := $$(pwd)
	translate = $1
  program = okular
else
  # Windows mit MSys2/Cygwin
	pwd := $$(cygpath -m "$$(pwd)")
	translate = $(shell echo '$1' | sed 's/:/;/g')
  program = sumatraPDF
endif

all: build/thesis.pdf python/scaling.py

thoughts: build/thoughts.pdf


TeXOptions = -pdf \
			 -interaction=nonstopmode \
			 -halt-on-error \
			 -output-directory=build \
			 -shell-escape

python: python/PythonSkript.py python/scaling.py python/meine-matplotlibrc

python/scaling.py: FORCE
	python python/scaling.py

python/PythonSkript.py: FORCE
	TEXINPUTS="$(call translate,$(pwd):)" python python/PythonSkript.py

build/thesis.pdf: FORCE | build
	latexmk $(TeXOptions) thesis.tex

build/thoughts.pdf: FORCE | build
	latexmk $(TeXOptions) thoughts.tex

FORCE:

build:
	mkdir -p build/

clean:
	rm -rf build

once: build python/scaling.py
	pdflatex --output-directory=build --shell-escape thesis.tex

see:
	$(program) build/thesis.pdf

seet:
	$(program) build/thoughts.pdf
