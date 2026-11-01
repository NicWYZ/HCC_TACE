# Reproduce the full analysis from data/raw/liver.csv.
# Requires the package to be installed (pip install -e .) and R with MatchIt/cobalt.

NB_RUN = jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=3600

.PHONY: all cohort matching analysis exploratory clean

all: cohort matching analysis

cohort:
	$(NB_RUN) notebooks/01_cohort_construction.ipynb

matching:
	Rscript R/psm_matching.R

analysis:
	$(NB_RUN) notebooks/02_covariate_balance.ipynb
	$(NB_RUN) notebooks/03_tace_vs_systemic.ipynb
	$(NB_RUN) notebooks/04_three_group_comparison.ipynb
	$(NB_RUN) notebooks/05_stage_subgroups.ipynb
	$(NB_RUN) notebooks/06_interval_optimization.ipynb

exploratory:
	$(NB_RUN) notebooks/exploratory/mice_imputation_all.ipynb
	$(NB_RUN) notebooks/exploratory/mice_imputation_predictive.ipynb

clean:
	rm -rf results/figures/* results/tables/*
