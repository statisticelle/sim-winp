# sim-winp
R code to reproduce simulation study for:
Generalized pairwise comparisons for cluster randomized trials using win fractions

Steps for running the simulation study are as follows:
1. scenarios.R - Input parameter values of interest and generate all scenarios for data generation/investigation. Save to scenarios.csv. 
2. main.R - Load simulation scripts and other dependencies (stored in dependencies.R), import scenarios.csv, and simulate all scenarios using simScenarios() function.

Scripts within src/
1. help.R - Helper functions (e.g., logit, correlation matrix construction, identification of parameter values for given win probability).
3. genResponse.R - Function to generate bivariate responses given simulation scenario parameters.
4. simScenarios.R - Function which runs all simulation scenarios, and returns results. 

Files within out/
1. scenarios.csv - Dataframe of scenarios output by scenarios.R
2. simResults.csv - Dataframe of results output by simScenarios() within main.R
