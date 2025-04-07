##This analysis tests for the significance of the phylogenetic signal (lambda) by using a likelihood ratio test method where a model with the observed lambda is compared against one that has no phylogenetic signal.

#load libraries
library(nlme)
library(ape)
library(phytools)
library(car)
#read tree
tree = read.nexus("Heliconiini_full sp name.trees")
#view sequence of species in tree
tree$tip.label
#make ultrametric tree
tree<-force.ultrametric(tree, method=c("nnls"),)

#create data and data.frame, used average species value 
ALtrait_data<-c("Philaethria dido"=7.119276567, "Drydaula phaetusa"=7.148002722, "Podotricha telesephe"=6.825322527, "Dryas iulia"=7.111033408, "Dione juno"=6.982604004,
                "Agraulis vanillae"=7.068104737, "E.aliphera"=6.883686154, "E.tales"=7.093350697, "E.lybia"=7.266286463, "E.isabella"=7.078000132, "E.vibilia"=7.051488232,
                "H.telesephe"=7.159154815, "H.clysonymous"=7.249237711, "H.hecalesia"=7.235096707, "H.himera"=7.126808816, "H.erato cyrbia"=7.059489175, "H.erato demophoon"=7.0737583,
                "H.charithonia"=6.986123759, "H.demeter"=7.113703652, "H.eratosignis"=7.101468315, "H.ricini"=6.990585608, "H.sara"=7.110131952, "H.eleuchia"=7.135341913,
                "H.sapho"=7.062662931, "H.hewitsoni"=7.238464026, "H.aoede"=7.029598007, "H.wallacei"=7.092290226, "H.burneyi"=7.148820858, "H.ismenius"=7.074963278,
                "H.numata"=7.068583103, "H.ethila"=7.110758125, "H.atthis"=7.121642384, "H.hecale"=7.084968396, "H.m.amaryllis"=7.056737822, "H.m.melpomene"=7.141849524,
                "H.melpomene rosina"=7.087413634, "H.timareta"=7.115678685, "H.cydno chioneus"=7.120734001, "H.pachinus"=7.123662127, "H.cydno galanthus"=7.133448075,
                "H.doris"=7.060702249)

ALHtrait_data<-c("Philaethria dido"=6.222278556, "Drydaula phaetusa"=6.341083326, "Podotricha telesephe"=5.973327249, "Dryas iulia"=6.109449433, "Dione juno"=5.942759914,
                 "Agraulis vanillae"=6.145363321, "E.aliphera"=5.918132161, "E.tales"=6.103000561, "E.lybia"=6.437419495, "E.isabella"=5.950027322, "E.vibilia"=6.2101426,
                 "H.telesephe"=6.194579013, "H.clysonymous"=6.253210839, "H.hecalesia"=6.172675056, "H.himera"=6.538781955, "H.erato cyrbia"=6.495044474, "H.erato demophoon"=6.277432669,
                 "H.charithonia"=6.070267557, "H.demeter"=6.112344308, "H.eratosignis"=6.033856018, "H.ricini"=5.95022368, "H.sara"=6.23019231, "H.eleuchia"=6.474459307,
                 "H.sapho"=6.320564005, "H.hewitsoni"=6.308292323, "H.aoede"=6.302742385, "H.wallacei"=6.163304891, "H.burneyi"=6.207338981, "H.ismenius"=6.204557087,
                 "H.numata"=6.031209353, "H.ethila"=6.132451588, "H.atthis"=6.286391681, "H.hecale"=6.247348695, "H.m.amaryllis"=6.026447284, "H.m.melpomene"=6.160147161,
                 "H.melpomene rosina"=6.271890579, "H.timareta"=6.245365466, "H.cydno chioneus"=6.269561479, "H.pachinus"=6.260129635, "H.cydno galanthus"=6.186305014,
                 "H.doris"=5.996802499)
GLtrait_data<-c("Philaethria dido"=7.06040732, "Drydaula phaetusa"=7.065317594, "Podotricha telesephe"=6.909591556, "Dryas iulia"=7.065451558, "Dione juno"=6.941056155,
                "Agraulis vanillae"=7.030764747, "E.aliphera"=6.833927074, "E.tales"=7.046507643, "E.lybia"=7.196574524, "E.isabella"=7.044386944, "E.vibilia"=7.01162721,
                "H.telesephe"=7.124243092, "H.clysonymous"=7.20303651, "H.hecalesia"=7.205040644, "H.himera"=6.997089762, "H.erato cyrbia"=6.921251436, "H.erato demophoon"=6.998122412,
                "H.charithonia"=6.936705339, "H.demeter"=7.068096938, "H.eratosignis"=7.062612588, "H.ricini"=6.949089713, "H.sara"=7.064303957, "H.eleuchia"=7.028364246,
                "H.sapho"=6.975897645, "H.hewitsoni"=7.184206635, "H.aoede"=6.966633918, "H.wallacei"=7.037874719, "H.burneyi"=7.096046567, "H.ismenius"=7.012094757,
                "H.numata"=7.026786472, "H.ethila"=7.062521974, "H.atthis"=7.053032443, "H.hecale"=7.0167626, "H.m.amaryllis"=7.014219151, "H.m.melpomene"=7.080819566,
                "H.melpomene rosina"=7.015335777, "H.timareta"=7.05279568, "H.cydno chioneus"=7.054791416, "H.pachinus"=7.059713082, "H.cydno galanthus"=7.081399906,
                "H.doris"=7.016914918)
rCBRtrait_data<-c("Philaethria dido"=8.081347311, "Drydaula phaetusa"=8.11380331, "Podotricha telesephe"=7.932016278, "Dryas iulia"=8.119077252, "Dione juno"=7.914671288,
                  "Agraulis vanillae"=7.992079705, "E.aliphera"=7.717954248, "E.tales"=7.904020848, "E.lybia"=7.971481584, "E.isabella"=7.952059071, "E.vibilia"=7.927167604,
                  "H.telesephe"=8.262657918, "H.clysonymous"=8.348795177, "H.hecalesia"=8.240285543, "H.himera"=8.069644847, "H.erato cyrbia"=8.068035936, "H.erato demophoon"=8.03380339,
                  "H.charithonia"=8.02350345, "H.demeter"=8.11668086, "H.eratosignis"=8.085650145, "H.ricini"=7.96534579, "H.sara"=7.946939863, "H.eleuchia"=8.034349994,
                  "H.sapho"=8.029733523, "H.hewitsoni"=8.181971474, "H.aoede"=7.998166221, "H.wallacei"=8.089269234, "H.burneyi"=8.170905623, "H.ismenius"=8.139512877,
                  "H.numata"=8.061643999, "H.ethila"=8.155425137, "H.atthis"=8.17755147, "H.hecale"=8.202257883, "H.m.amaryllis"=8.05478235, "H.m.melpomene"=8.17747411,
                  "H.melpomene rosina"=8.091831569, "H.timareta"=8.12730993, "H.cydno chioneus"=8.120980954, "H.pachinus"=8.152750038, "H.cydno galanthus"=8.154499068,
                  "H.doris"=8.045951153)

#create models
ALrCBR_df<-data.frame(
  Species = names(ALtrait_data),
  Antennal_lobe = ALtrait_data,
  rCBR = rCBRtrait_data
)

ALHrCBR_df<-data.frame(
  Species = names(ALHtrait_data),
  Antennal_lobe_hub = ALHtrait_data,
  rCBR = rCBRtrait_data
)

GLrCBR_df<-data.frame(
  Species = names(GLtrait_data),
  Glomeruli = GLtrait_data,
  rCBR = rCBRtrait_data
)


bothtrait_df<-data.frame(
  Species = names(GLtrait_data),
  Glomeruli = GLtrait_data,
  Hub = ALHtrait_data
)

#lambda obtained from previous MCMC analysed (see script MCMCglmm phylogenetic)
AL_lambda_mean=0.7101038
ALH_lambda_mean=0.671386
GL_lambda_mean=0.698
both_lambda_mean=0.6053521
#test AL 
# Fit the model with the observed lambda
ALmodel_lambda <- gls(Antennal_lobe ~ rCBR, correlation = corPagel(value = AL_lambda_mean, phy = tree), data = ALrCBR_df, method = "ML")

# Fit the model assuming lambda = 0 (no phylogenetic signal), convergence issues with lambda at 0, so set a small nonzero starting value
ALmodel_no_lambda <- gls(Antennal_lobe ~ rCBR, correlation = corPagel(value = 0.01, phy = tree), data = ALrCBR_df, method = "ML")


# Compare the log-likelihoods
AL_LRT_stat <- -2 * (logLik(ALmodel_no_lambda) - logLik(ALmodel_lambda))

# Compute p-value using chi-square distribution with 1 degree of freedom
AL_p_value <- pchisq(AL_LRT_stat, df = 1, lower.tail = FALSE)

# Print results
cat("Likelihood Ratio Test statistic:", AL_LRT_stat, "\n")
cat("p-value:", AL_p_value, "\n")



#test ALH 
# Fit the model with the observed lambda
ALHmodel_lambda <- gls(Antennal_lobe_hub ~ rCBR, correlation = corPagel(value = ALH_lambda_mean, phy = tree), data = ALHrCBR_df, method = "REML")
# Fit the model assuming lambda = 0 (no phylogenetic signal)
ALHmodel_no_lambda <- gls(Antennal_lobe_hub ~ rCBR, correlation = corPagel(value = 0, phy = tree), data = ALHrCBR_df, method = "ML")
# Compare the log-likelihoods
ALH_LRT_stat <- -2 * (logLik(ALHmodel_no_lambda) - logLik(ALHmodel_lambda))

# Compute p-value using chi-square distribution with 1 degree of freedom
ALH_p_value <- pchisq(ALH_LRT_stat, df = 1, lower.tail = FALSE)

# Print results
cat("Likelihood Ratio Test statistic:", ALH_LRT_stat, "\n")
cat("p-value:", ALH_p_value, "\n")


#test GL 
# Fit the model with the observed lambda
GLmodel_lambda <- gls(Glomeruli ~ rCBR, correlation = corPagel(value = GL_lambda_mean, phy = tree), data = GLrCBR_df, method = "ML")
# Fit the model assuming lambda = 0 (no phylogenetic signal)
GLmodel_no_lambda <- gls(Glomeruli ~ rCBR, correlation = corPagel(value = 0, phy = tree), data = GLrCBR_df, method = "ML")
# Compare the log-likelihoods
GL_LRT_stat <- -2 * (logLik(GLmodel_no_lambda) - logLik(GLmodel_lambda))

# Compute p-value using chi-square distribution with 1 degree of freedom
GL_p_value <- pchisq(GL_LRT_stat, df = 1, lower.tail = FALSE)

# Print results
cat("Likelihood Ratio Test statistic:", GL_LRT_stat, "\n")
cat("p-value:", GL_p_value, "\n")


#test GL~ALH 
# Fit the model with the observed lambda
bothmodel_lambda <- gls(Glomeruli ~ Hub, correlation = corPagel(value = both_lambda_mean, phy = tree), data = bothtrait_df, method = "REML")
# Fit the model assuming lambda = 0 (no phylogenetic signal)
bothmodel_no_lambda <- gls(Glomeruli ~ Hub, correlation = corPagel(value = 0, phy = tree), data = bothtrait_df, method = "ML")
# Compare the log-likelihoods
both_LRT_stat <- -2 * (logLik(bothmodel_no_lambda) - logLik(bothmodel_lambda))

# Compute p-value using chi-square distribution with 1 degree of freedom
both_p_value <- pchisq(both_LRT_stat, df = 1, lower.tail = FALSE)

# Print results
cat("Likelihood Ratio Test statistic:", both_LRT_stat, "\n")
cat("p-value:", both_p_value, "\n")
