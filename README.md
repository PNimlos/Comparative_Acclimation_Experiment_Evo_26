# Comparative_Acclimation_Experiment_Evo_26
R code for running analyses from Nimlos et al. 2026 (doi).
We acclimated six taxa of Peromyscus mice spanning a gradient of native elevations (0-4,300m) to laboratory-generated hypobaria to test whether flexible responses to hypoxia have evolved predictably. See paper for detailed methods.

Gene_Expr_and_Flex_modeling.R contains scripts necessary to model the expression and flexibility of expression across native elevation, including expression normalization, GLS and PGLS modeling, and visualization through ggplot2. Input data for these scripts can be found in both acclim_data_all.csv and cpm_counts.csv.

Physio_trait_modeling.R uses the same framework to test for a relationship between various physiological traits (and their flexibility) and native elevation. Input data is only acclim_data_all.csv. 
