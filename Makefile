SHELL := bash
.ONESHELL:
.SHELLFLAGS := -euio pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

### Variables ###
# Tools
PYTEST           = pytest
BASH             = bash
CONDA            = conda
PYTHON           = python3.8
SNAKEMAKE        = snakemake
CONDA_ACTIVATE   = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

# Paths
TEST_COMMON      = scripts/common_script_vass.py
TEST_PIPELINE    = scripts/prepare_pipeline.py
SNAKE_FILE       = Snakefile
ENV_YAML         = envs/workflow.yaml
GWASCAT_PATH     = tests/annotations/test.gwascat.tsv
DBNSFP_PATH      = tests/annotations/test.dbNSFP.vcf.gz
GENESETS_PATH    = tests/annotations/test.gmt
VCF_PATH         = tests/vcfs/empty

# Arguments
ENV_NAME         = vcf-annotate-snpeff-snpsift
SNAKE_THREADS    = 1
SNPEFF_ARGS      = '\-noGenome'

# Recipes
default: all-unit-tests

# Environment building through conda
conda-install:
	${CONDA_ACTIVATE} base && \
	${CONDA} env create --file ${ENV_YAML} --force && \
	${CONDA} activate ${ENV_NAME}

### UNIT TESTS ###
# Running all unit-tests (one for each python scripts)
all-unit-tests:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTEST} -v ${TEST_COMMON} ${TEST_PIPELINE}
.PHONY: all-unit-tests

# Running all unit test (on common_script_vass.py only)
common-tests:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTEST} -vv ${TEST_COMMON}
.PHONY: common-tests

### Continuous Integration Tests ###
# Running snakemake on test datasets
test-conda-report.html:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTHON} ${TEST_PIPELINE} --force --recursive --debug ${VCF_PATH} ${GWASCAT_PATH} ${GENESETS_PATH} ${DBNSFP_PATH} --snpeff-extra ${SNPEFF_ARGS} --workdir ${PWD}/tests &&
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --forceall --printshellcmds --reason --directory ${PWD}/tests && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --directory ${PWD}/tests --report test-conda-report.html

# Running snakemake on test datasets with singularity flag raised on
test-singularity-report.html:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTHON} ${TEST_PIPELINE} --force --recursive --debug ${VCF_PATH} ${GWASCAT_PATH} ${GENESETS_PATH} ${DBNSFP_PATH} --snpeff-extra ${SNPEFF_ARGS} --workdir ${PWD}/tests &&
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --forceall --printshellcmds --reason --directory ${PWD}/tests --use-singularity&& \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --directory ${PWD}/tests --report singularity-tests.html

# Cleaning Snakemake outputs
clean:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --force --configfile ${PWD}/tests/config.yaml --use-singularity --directory ${PWD}/tests --delete-all-output
.PHONY: clean


# Display pipeline graph
workflow.png:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --force --configfile ${PWD}/tests/config.yaml --directory ${PWD}/tests --rulegraph | dot -T png > workflow.png

example.png:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --force --configfile ${PWD}/tests/config.yaml --directory ${PWD}/tests --dag | dot -T png > example.png
