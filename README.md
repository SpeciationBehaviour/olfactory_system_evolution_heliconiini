# R Analysis scripts for *Evolution of the Olfactory system during the radiation of Heliconiini butterflies*
Analysis scripts for "Evolution of the Olfactory system during the radiation of Heliconiini butterflies"

These documents contains scripts for Toh et al., 2025 “Evolution of the Olfactory system during the radiation of Heliconiini butterflies”. All starting data used in these analyses can be found in the supplementary materials under “Heliconiini_olfactory_evolution_supplementary_tables.xlsx”

This repository contains 6 annotated scripts.

1) Scripts for phylogenetically controlled comparative analysis of neuropil volumes and olfactory receptor numbers of the Heliconiini

a. MCMCglmm phylogenetic.R - tests for broad phylogenetic relationships as well as pairwise comparisons of neuropil volumes in the Heliconiini. This script also contains analyses which tests for the strength of the phylogenetic signal (lambda) for the different neuropils.

b. MCMCglmm ecology and environment.R- tests for significant environmental and ecological variables affecting neuropil volume

c. MCMCglmm OR analysis.R- tests whether neuropil volumes of each species are affected by the number of olfactory receptors , by using BUSCO scores or N50 scores as a control. The script also contains analyses which tests for significant environmental and ecological variables affecting the number of olfactory receptors for each species.

2) Likelihood ratio testing of phylogenetic signal lambda

a. likelihood ratio testing of phylogenetic (lambda) values.R- tests for the significance of the phylogenetic signal (lambda) by using a likelihood ratio test method where a model with the observed lambda is compared against one that has no phylogenetic signal.

3) Scripts to test for differences in neuropil volumes of wild vs insectary individuals of Heliconius butterflies

a. lmer.R- This analysis runs a linear mixed model to test for the significant variables affecting variation in neuropil volume

b. SMATR.R- This analysis tests for differences in slope between two groups (wild vs insectary) and subsequently the presence of elevational or major-axis shifts if slopes are not significantly different
