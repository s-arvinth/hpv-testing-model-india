# Common tasks. Run `make help` for the list.
PY ?= python
WORKERS ?=
W = $(if $(WORKERS),--workers $(WORKERS))

.PHONY: help install inputs examples test test-full paper paper-simulate paper-tables paper-figures clean

help:              ## list targets
	@grep -E '^[a-z-]+:.*##' $(MAKEFILE_LIST) | sed 's/:.*##/\t/'

install:           ## install the package with maps and test extras
	$(PY) -m pip install -e ".[maps,dev]"

inputs:            ## rebuild data/processed/units from data/raw
	$(PY) -m hpvcc prepare-inputs

examples:          ## run every script in examples/ (outputs in results/examples)
	@for f in examples/0*.py; do echo "== $$f"; $(PY) $$f || exit 1; done

test:              ## fast checks (~1 min)
	$(PY) -m pytest -m "not slow"

test-full:         ## + regenerated ensembles and LaTeX tables vs the reference ones
	$(PY) -m pytest

paper:             ## reproduce the paper: simulations, tables, figures (~15 min on 40 cores)
	$(PY) -m hpvcc $(W) paper all

paper-simulate:    ## the paper's simulations only
	$(PY) -m hpvcc $(W) paper simulate

paper-tables:      ## the paper's tables
	$(PY) -m hpvcc paper tables

paper-figures:     ## the paper's figures
	$(PY) -m hpvcc paper figure fig3 fig4 fig5 fig6 figS1 figS2 figS3 figS4 \
	    figS11 figS12 figS13 figS14 figS15 figS16

clean:             ## delete results/
	rm -rf results
