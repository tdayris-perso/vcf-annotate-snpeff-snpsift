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
PYTHON           = python3.7
SNAKEMAKE        = snakemake
CONDA_ACTIVATE   = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate

# Paths
TEST_CONFIG      = scripts/prepare_config.py
TEST_DESIGN      = scripts/prepare_design.py
SNAKE_FILE       = Snakefile
ENV_YAML         = envs/workflow.yaml
GWASCAT_PATH     = annotations/test.gwascat.tsv
DBNSFP_PATH      = annotations/test.dbNSFP.tsv.gz
GENESETS_PATH    = annotations/test.gmt
VCF_PATH         = tests/vcfs/

# Arguments
ENV_NAME         = vcf-annotate-snpeff-snpsift
SNAKE_THREADS    = 1
SNPEFF_ARGS      = '\-noGenome'

# Recipes
default: all-unit-tests

# Environment building through conda
conda-tests:
	${CONDA_ACTIVATE} base && \
	${CONDA} env create --file ${ENV_YAML} --force && \
	${CONDA} activate ${ENV_NAME}

### UNIT TESTS ###
# Running all unit-tests (one for each python scripts)
all-unit-tests:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTEST} -v ${TEST_CONFIG} ${TEST_DESIGN}
.PHONY: all-unit-tests

# Running all unit test (on prepare_config.py only)
config-tests:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTEST} -v ${TEST_CONFIG} && \
	${PYTHON} ${TEST_CONFIG} ${GWASCAT_PATH} ${GENESETS_PATH} ${DBNSFP_PATH} --workdir ${PWD}/tests --debug
.PHONY: config-tests

# Running all unit test (on prepare_design.py only)
design-tests:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTEST} -v ${TEST_DESIGN} && \
	${PYTHON} ${TEST_DESIGN} ${VCF_PATH} -o ${PWD}/tests/design.tsv --debug
.PHONY: design-tests

### Continuous Integration Tests ###
# Running snakemake on test datasets
ci-tests:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTHON} ${TEST_DESIGN} ${VCF_PATH} -o ${PWD}/tests/design.tsv --debug && \
	${PYTHON} ${TEST_CONFIG} ${TRANSCRIPT_PATH} ${GWASCAT_PATH} ${GENESETS_PATH} ${DBNSFP_PATH} --snpeff-extra ${SNPEFF_ARGS} --workdir ${PWD}/tests --debug  && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --forceall --printshellcmds --reason --directory ${PWD}/tests && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --directory ${PWD}/tests --report
.PHONY: ci-tests

# Running snakemake on test datasets with singularity flag raised on
singularity-tests:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTHON} ${TEST_DESIGN} ${VCF_PATH} -o ${PWD}/tests/design.tsv --debug && \
	${PYTHON} ${TEST_CONFIG} ${TRANSCRIPT_PATH} ${GWASCAT_PATH} ${GENESETS_PATH} ${DBNSFP_PATH} --snpeff-extra ${SNPEFF_ARGS} --workdir ${PWD}/tests --debug && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --forceall --printshellcmds --reason --directory ${PWD}/tests --use-singularity && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --directory ${PWD}/tests --report
.PHONY: singularity-tests

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
