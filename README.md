# Immunosuppression-Project

# Survival Analysis
To generate the results for survival analysis:

1. Download perform_survival_analysis.R, liver_dat_cleaned.csv, and immuno_dat_cleaned.csv into the same folder. 
2. Change the working directory provided in setwd() to your current working directory. 
3. Source the R file to generate the patient survival curves and graft survival curves for the 3 age groups (18 - 49, 50 - 64, and 65+) as well as a pvalue.txt file containing the p-values from the log-rank test in LaTeX format. 
4. The extra discard.txt file can be discarded.

# Propensity Score Matching
To generate the results for propensity score matching: 
1. Download propensity_score_matching.R, liver_dat_large.csv, and immuno_dat_small.csv into the same folder.
2. Change the working directory provided in setwd() to your current working directory.
3. Source the R file to generate the stacked histograms (stored in stacked_histograms.pdf) and the kernel density plots (stored in kernel_density_plots.pdf), two pdf files which provide data visualization of the propensity scores from the logistic regression classifier.
