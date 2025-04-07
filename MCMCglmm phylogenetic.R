##This analysis tests for broad phylogenetic relationships as well as pairwise comparisons of neuropil volumes in the Heliconiini while taking into account the phylogenetic relationship amongst species.
##This analysis also tests for the strength of the phylogenetic signal (lambda) for the different neuropils.  

#run MCMCglmm models 
#load libraries
library(MCMCglmm)
library(phytools)
#read in metadata 
data<-read.csv("Table S8-metadata.csv",header=T)
#read in filtered data for pairwise comparison
filtered_data<-read.csv("Table S9-metadata-filtered.csv",header=T)
#structuring data for analysis
#note: sp.app= species name abbrevation, sim.clade=simplified clade name 
data<-data[,-c(21:24)]
data$sim.clade<-as.factor(data$sim.clade)
data$Heliconius<-as.factor(data$Heliconius)
data$sex<-as.factor(data$sex)
data$sp.abb<-as.factor(data$sp.abb)
data$genus<-as.factor(data$genus)
filtered_data$sim.clade<-as.factor(filtered_data$sim.clade)
filtered_data$Heliconius<-as.factor(filtered_data$Heliconius)
filtered_data$sex<-as.factor(filtered_data$sex)
filtered_data$sp.abb<-as.factor(filtered_data$sp.abb)
filtered_data$genus<-as.factor(filtered_data$genus)
#omit NAs
data<-na.omit(data)
filtered_data<-na.omit(filtered_data)
#add logs to AL ALH, GL and rCBR
data$AL.log <-log10(data$AL)
data$ALH.log<-log10(data$ALH)
data$GL.log<-log10(data$GL)
data$rCBR.log<-log10(data$rCBR) 
filtered_data$AL.log <-log10(filtered_data$AL)
filtered_data$ALH.log<-log10(filtered_data$ALH)
filtered_data$GL.log<-log10(filtered_data$GL)
filtered_data$rCBR.log<-log10(filtered_data$rCBR) 
#check data structure 
str(data)
str(filtered_data)

#filter data for Heliconius genus testing 
Helionly<-subset(data, Heliconius=="y")
#drop groups that are not in used 
Helionly<-droplevels(Helionly)

## set priors, essential for Bayesian analysis
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
#read in phylogenetic tree of Helliconiini
tree = read.nexus("Heliconiini.trees")
#make ultrametric tree
tree<-force.ultrametric(tree, method=c("nnls"),)
#make inverse of tree/covariance matrix
inv.phylo<-inverseA(tree,nodes="TIPS",scale=TRUE)
#AL  Heliconius vs non-Heliconius
#define model, with interaction with sex
sm1<-AL.log~rCBR.log+Heliconius*sex
## create models
smodel1<-MCMCglmm(sm1, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)

smodel2<-MCMCglmm(sm1, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)
##model checks
## check convergence
gelman.diag(mcmc.list(smodel1$Sol,smodel2$Sol))
gelman.diag(mcmc.list(smodel1$VCV,smodel2$VCV))
#visual representation of the output from gelman.diag convergence check
plot(smodel1$VCV)
plot(smodel1$Sol)
plot(smodel2$VCV)
plot(smodel2$Sol)
##check autocorrelation
autocorr(smodel1$Sol)
autocorr(smodel1$VCV)
autocorr(smodel2$Sol)
autocorr(smodel2$VCV)
## summary
summary(smodel1)

#model simplification (-interaction)
#no effect of Heliconius*sex
#define model 
sm1.1<-AL.log~rCBR.log+Heliconius+sex
## create models
smodel1.1<-MCMCglmm(sm1.1, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)

smodel1.2<-MCMCglmm(sm1.1, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)
##model checks
## check convergence
gelman.diag(mcmc.list(smodel1.1$Sol,smodel1.2$Sol))
gelman.diag(mcmc.list(smodel1.1$VCV,smodel1.2$VCV))

#visual representation of the output from gelman.diag convergence check
plot(smodel1.1$VCV)
plot(smodel1.1$Sol)

plot(smodel1.2$VCV)
plot(smodel1.2$Sol)
##check autocorrelation
autocorr(smodel1.1$Sol)
autocorr(smodel1.1$VCV)
autocorr(smodel1.2$Sol)
autocorr(smodel1.2$VCV)
## summary
summary(smodel1.1)

#model simplification (-Heliconius), test for sex effects
#define model
sm2<-AL.log~rCBR.log+sex
##create models
smodel3<-MCMCglmm(sm2, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)

smodel4<-MCMCglmm(sm2, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)
##model checks
## check convergence
gelman.diag(mcmc.list(smodel3$Sol,smodel4$Sol))
gelman.diag(mcmc.list(smodel3$VCV,smodel4$VCV))
#visual representation of the output from gelman.diag convergence check
plot(smodel3$VCV)
plot(smodel3$Sol)
plot(smodel4$VCV)
plot(smodel4$Sol)

##check autocorrelation
autocorr(smodel3$Sol)
autocorr(smodel3$VCV)
autocorr(smodel4$Sol)
autocorr(smodel4$VCV)
## summary
summary(smodel3)
#no sex differences! 

#investigate major clade differences 
#define model
sm3<-AL.log~rCBR.log+genus*sex
##create models
smodel5<-MCMCglmm(sm3, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)

smodel6<-MCMCglmm(sm3, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)
##model checks
## check convergence
gelman.diag(mcmc.list(smodel5$Sol,smodel6$Sol))
gelman.diag(mcmc.list(smodel5$VCV,smodel6$VCV))
#visual representation of the output from gelman.diag convergence check
plot(smodel5$VCV)
plot(smodel5$Sol)
plot(smodel6$VCV)
plot(smodel6$Sol)
##check autocorrelation
autocorr(smodel5$Sol)
autocorr(smodel5$VCV)
autocorr(smodel6$Sol)
autocorr(smodel6$VCV)
#summary
summary(smodel5)

#model simplification (-interaction)
#define model
sm4<-AL.log~rCBR.log+genus+sex
##create models
smodel7<-MCMCglmm(sm4, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)

smodel8<-MCMCglmm(sm4, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)
##model checks
# check convergence
gelman.diag(mcmc.list(smodel7$Sol,smodel8$Sol))
gelman.diag(mcmc.list(smodel7$VCV,smodel8$VCV))
#visual representation of the output from gelman.diag convergence check
plot(smodel7$VCV)
plot(smodel7$Sol)
plot(smodel8$VCV)
plot(smodel8$Sol)
##check autocorrelation
autocorr(smodel7$Sol)
autocorr(smodel7$VCV)
autocorr(smodel8$Sol)
autocorr(smodel8$VCV)
#summary
summary(smodel7)
#no significant effect of genus
#model simplification (-genus), same as sex model above 

##AL clade differences within Heliconius 
#define model
gmodel1<-AL.log~rCBR.log+sim.clade*sex
##create models 
gmodel1.1<-MCMCglmm(gmodel1, random=~sp.abb, family=c("gaussian"), 
                    ginverse=list(sp.abb=inv.phylo$Ainv),
                    prior=prior, data=Helionly,
                    nitt=500000,burnin=10000,thin=500)

gmodel1.2<-MCMCglmm(gmodel1, random=~sp.abb, family=c("gaussian"), 
                    ginverse=list(sp.abb=inv.phylo$Ainv),
                    prior=prior, data=Helionly,
                    nitt=500000,burnin=10000,thin=500)
##model checks
## check convergence
gelman.diag(mcmc.list(gmodel1.1$Sol,gmodel1.2$Sol))
gelman.diag(mcmc.list(gmodel1.1$VCV,gmodel1.2$VCV))
#visual representation of the output from gelman.diag convergence check
plot(gmodel.1$VCV)
plot(gmodel.1$Sol)
plot(gmodel.2$VCV)
plot(gmodel.2$Sol)
##check autocorrelation
autocorr(gmodel.1$Sol)
autocorr(gmodel.1$VCV)
autocorr(gmodel.2$Sol)
autocorr(gmodel.2$VCV)
#summary
summary(gmodel1.1)
#interaction not significant
#model simplification (-interaction)
#define model
gmodel2<-AL.log~rCBR.log+sim.clade+sex
#create models
gmodel2.1<-MCMCglmm(gmodel2, random=~sp.abb, family=c("gaussian"), 
                    ginverse=list(sp.abb=inv.phylo$Ainv),
                    prior=prior, data=Helionly,
                    nitt=500000,burnin=10000,thin=500)

gmodel2.2<-MCMCglmm(gmodel2, random=~sp.abb, family=c("gaussian"), 
                    ginverse=list(sp.abb=inv.phylo$Ainv),
                    prior=prior, data=Helionly,
                    nitt=500000,burnin=10000,thin=500)
## check convergence
gelman.diag(mcmc.list(gmodel2.1$Sol,gmodel2.2$Sol))
gelman.diag(mcmc.list(gmodel2.1$VCV,gmodel2.2$VCV))  
#visual representation of the output from gelman.diag convergence check
plot(gmodel2.1$VCV)
plot(gmodel2.1$Sol)
plot(gmodel2.2$VCV)
plot(gmodel2.2$Sol)
##check autocorrelation
autocorr(gmodel2.1$Sol)
autocorr(gmodel2.1$VCV)
autocorr(gmodel2.2$Sol)
autocorr(gmodel2.2$VCV)
#summary
summary(gmodel2.1)
#sim.clade not significant 
#model simplification (-sim.clade), same as sex model above 
##no difference detected for all 3 phylogenetic groupings, run NULL model to get coefficients for plotting 
#define model
nullmodel<-AL.log~rCBR.log
#create models
nullmodel1<-MCMCglmm(nullmodel, random=~sp.abb, family=c("gaussian"), 
                     ginverse=list(sp.abb=inv.phylo$Ainv),
                     prior=prior, data=data,
                     nitt=500000,burnin=10000,thin=500)

nullmodel2<-MCMCglmm(nullmodel, random=~sp.abb, family=c("gaussian"), 
                     ginverse=list(sp.abb=inv.phylo$Ainv),
                     prior=prior, data=data,
                     nitt=500000,burnin=10000,thin=500)
#model checks 
# check convergence
gelman.diag(mcmc.list(nullmodel1$Sol,nullmodel2$Sol))
gelman.diag(mcmc.list(nullmodel1$VCV,nullmodel2$VCV))
#visual representation of the output from gelman.diag convergence check
plot(nullmodel1$VCV)
plot(nullmodel1$Sol)
plot(nullmodel2$VCV)
plot(nullmodel2$Sol)
#check autocorrelation
autocorr(nullmodel1$Sol)
autocorr(nullmodel1$VCV)
autocorr(nullmodel2$Sol)
autocorr(nullmodel2$VCV)

# summary
summary(nullmodel1)


#extract phylogenetic signal lambda (mean)
lambdaAL<-nullmodel1$VCV[,'sp.abb']/(nullmodel1$VCV[,'sp.abb']+nullmodel1$VCV[,'units'])
mean(lambdaAL)

#compare significance of the phylogenetic signal 
#run model without phylogeny 
nophylo1<-MCMCglmm(nullmodel, random=~units, family=c("gaussian"), 
                   prior=prior, data=data,
                   nitt=500000,burnin=10000,thin=500)

# Compute log-likelihoods
logL_full <- nullmodel1$DIC  # DIC for full model
logL_null <- nophylo1$DIC  # DIC for no phylo model
# Compute likelihood ratio test statistic
LRT_stat <- logL_null - logL_full
p_value <- pchisq(LRT_stat, df = 1, lower.tail = FALSE)
# Print results
cat("Likelihood Ratio Test statistic:", LRT_stat, "\n")
cat("p-value:", p_value, "\n")

#next test for Credible Interval (CI) of lambda
AL_lambda_CI <- quantile(lambdaAL, probs = c(0.025, 0.975))
# Print results
cat("Lambda mean:", mean(lambdaAL), "\n")
cat("95% credible interval:", AL_lambda_CI, "\n")

##species level pairwise comparisons (AL)
#use previous defined null model and let random effect be species, use filtered data 
#create models
nullmodel3<-MCMCglmm(nullmodel, random=~sp.abb, family=c("gaussian"), 
                     ginverse=list(sp.abb=inv.phylo$Ainv),
                     prior=prior, data=filtered_data,
                     nitt=500000,burnin=10000,thin=500,pr=T)

nullmodel4<-MCMCglmm(nullmodel, random=~sp.abb, family=c("gaussian"), 
                     ginverse=list(sp.abb=inv.phylo$Ainv),
                     prior=prior, data=filtered_data,
                     nitt=500000,burnin=10000,thin=500,pr=T)

#model checks 
# check convergence
gelman.diag(mcmc.list(nullmodel3$Sol,nullmodel4$Sol))
gelman.diag(mcmc.list(nullmodel3$VCV,nullmodel4$VCV))

#visual representation of the output from gelman.diag convergence check
plot(nullmodel3$VCV)
plot(nullmodel3$Sol)
plot(nullmodel4$VCV)
plot(nullmodel4$Sol)
#check autocorrelation
autocorr(nullmodel3$Sol)
autocorr(nullmodel3$VCV)
autocorr(nullmodel4$Sol)
autocorr(nullmodel4$VCV)

# summary
summary(nullmodel3)
# Extract posterior distribution for the random effects
random_effects <- nullmodel3$Sol
# Extract the random effects for species
random_effects_species <- random_effects[, grep("sp.abb", colnames(random_effects))]
# Convert to mcmc object
random_effects_species_mcmc <- as.mcmc(random_effects_species)
# Define a function for pairwise comparisons with a threshold
pairwise_comparisons <- function(effects_matrix, threshold = 0.975) {
  n <- ncol(effects_matrix)
  comparisons <- matrix(NA, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      diff_effect <- effects_matrix[, i] - effects_matrix[, j]
      prob_ij <- mean(diff_effect > 0)
      comparisons[i, j] <- prob_ij
      comparisons[j, i] <- 1 - prob_ij
    }
  }
  colnames(comparisons) <- rownames(comparisons) <- paste("Species", 1:n)
  comparisons_adjusted <- ifelse(comparisons >= threshold | comparisons <= (1 - threshold), comparisons, NA)
  return(list(comparisons = comparisons, adjusted = comparisons_adjusted))
}
# pairwise comparison and apply the function with a higher threshold to account for false discoveries
comparison_results <- pairwise_comparisons(random_effects_species_mcmc, threshold = 0.99)
#extract comparisons 
comparison_matrix <- comparison_results$comparisons
#save as csv.
write.csv(comparison_matrix, "AL pairwise_comparisons all values.csv")

##ALH  Heliconius vs non-Heliconius
#define model
ALH1<-ALH.log~rCBR.log+Heliconius*sex
## create model
ALHmodel1<-MCMCglmm(ALH1, random=~sp.abb, family=c("gaussian"), 
                    ginverse=list(sp.abb=inv.phylo$Ainv),
                    prior=prior, data=data,
                    nitt=500000,burnin=10000,thin=500)

ALHmodel2<-MCMCglmm(ALH1, random=~sp.abb, family=c("gaussian"), 
                    ginverse=list(sp.abb=inv.phylo$Ainv),
                    prior=prior, data=data,
                    nitt=500000,burnin=10000,thin=500)
#model checks
## check convergence
gelman.diag(mcmc.list(ALHmodel1$Sol,ALHmodel2$Sol))
gelman.diag(mcmc.list(ALHmodel1$VCV,ALHmodel2$VCV))
#visual representation of the output from gelman.diag convergence check
plot(ALHmodel1$VCV)
plot(ALHmodel1$Sol)
plot(ALHmodel2$VCV)
plot(ALHmodel2$Sol)
##check autocorrelation
autocorr(ALHmodel1$Sol)
autocorr(ALHmodel1$VCV)
autocorr(ALHmodel2$Sol)
autocorr(ALHmodel2$VCV)
# summary
summary(ALHmodel1)

#model simplification (-interaction)
#define model
ALH1.1<-ALH.log~rCBR.log+Heliconius+sex
#create models
ALHmodel1.1<-MCMCglmm(ALH1.1, random=~sp.abb, family=c("gaussian"), 
                    ginverse=list(sp.abb=inv.phylo$Ainv),
                    prior=prior, data=data,
                    nitt=500000,burnin=10000,thin=500)

ALHmodel1.2<-MCMCglmm(ALH1.1, random=~sp.abb, family=c("gaussian"), 
                    ginverse=list(sp.abb=inv.phylo$Ainv),
                    prior=prior, data=data,
                    nitt=500000,burnin=10000,thin=500)
#model checks
## check convergence
gelman.diag(mcmc.list(ALHmodel1.1$Sol,ALHmodel1.2$Sol))
gelman.diag(mcmc.list(ALHmodel1.1$VCV,ALHmodel1.2$VCV))
#visual representation of the output from gelman.diag convergence check
plot(ALHmodel1.1$VCV)
plot(ALHmodel1.1$Sol)
plot(ALHmodel1.2$VCV)
plot(ALHmodel1.2$Sol)
##check autocorrelation
autocorr(ALHmodel1.1$Sol)
autocorr(ALHmodel1.1$VCV)
autocorr(ALHmodel1.2$Sol)
autocorr(ALHmodel1.2$VCV)
# summary
summary(ALHmodel1.1)
#Heliconius term not significant
#model simplification (-Heliconius)
#ALH just testing for sex differences
#define model
ALH2<-ALH.log~rCBR.log+sex
#create models
ALHmodel3<-MCMCglmm(ALH2, random=~sp.abb, family=c("gaussian"), 
                    ginverse=list(sp.abb=inv.phylo$Ainv),
                    prior=prior, data=data,
                    nitt=500000,burnin=10000,thin=500)

ALHmodel4<-MCMCglmm(ALH2, random=~sp.abb, family=c("gaussian"), 
                    ginverse=list(sp.abb=inv.phylo$Ainv),
                    prior=prior, data=data,
                    nitt=500000,burnin=10000,thin=500)
#model checks
# check convergence
gelman.diag(mcmc.list(ALHmodel3$Sol,ALHmodel4$Sol))
gelman.diag(mcmc.list(ALHmodel3$VCV,ALHmodel4$VCV))
#visual representation of the output from gelman.diag convergence check
plot(ALHmodel3$VCV)
plot(ALHmodel3$Sol)
plot(ALHmodel4$VCV)
plot(ALHmodel4$Sol)
#check autocorrelation
autocorr(ALHmodel3$Sol)
autocorr(ALHmodel3$VCV)
autocorr(ALHmodel4$Sol)
autocorr(ALHmodel4$VCV)
# summary
summary(ALHmodel4)
#no sex differences 

##ALH Heliconiinii genera differences, investigate major clade differences 
#define model
ALH3<-ALH.log~rCBR.log+genus*sex
#create models
ALHmodel5<-MCMCglmm(ALH3, random=~sp.abb, family=c("gaussian"), 
                    ginverse=list(sp.abb=inv.phylo$Ainv),
                    prior=prior, data=data,
                    nitt=500000,burnin=10000,thin=500)

ALHmodel6<-MCMCglmm(ALH3, random=~sp.abb, family=c("gaussian"), 
                    ginverse=list(sp.abb=inv.phylo$Ainv),
                    prior=prior, data=data,
                    nitt=500000,burnin=10000,thin=500)
#model checks
## check convergence
gelman.diag(mcmc.list(ALHmodel5$Sol,ALHmodel6$Sol))
gelman.diag(mcmc.list(ALHmodel5$VCV,ALHmodel6$VCV))
#visual representation of the output from gelman.diag convergence check
plot(ALHmodel5$VCV)
plot(ALHmodel5$Sol)
plot(ALHmodel6$VCV)
plot(ALHmodel6$Sol)
##check autocorrelation
autocorr(ALHmodel5$Sol)
autocorr(ALHmodel5$VCV)
autocorr(ALHmodel6$Sol)
autocorr(ALHmodel6$VCV)
## summary
summary(ALHmodel6)

#model simplification (-interaction)
#define model
ALH3.1<-ALH.log~rCBR.log+genus+sex
#create models
ALHmodel5.1<-MCMCglmm(ALH3.1, random=~sp.abb, family=c("gaussian"), 
                    ginverse=list(sp.abb=inv.phylo$Ainv),
                    prior=prior, data=data,
                    nitt=500000,burnin=10000,thin=500)

ALHmodel6.1<-MCMCglmm(ALH3.1, random=~sp.abb, family=c("gaussian"), 
                    ginverse=list(sp.abb=inv.phylo$Ainv),
                    prior=prior, data=data,
                    nitt=500000,burnin=10000,thin=500)

#model checks
## check convergence
gelman.diag(mcmc.list(ALHmodel5.1$Sol,ALHmodel6.1$Sol))
gelman.diag(mcmc.list(ALHmodel5.1$VCV,ALHmodel6.1$VCV))
#visual representation of the output from gelman.diag convergence check
plot(ALHmodel5.1$VCV)
plot(ALHmodel5.1$Sol)
plot(ALHmodel6.1$VCV)
plot(ALHmodel6.1$Sol)
##check autocorrelation
autocorr(ALHmodel5.1$Sol)
autocorr(ALHmodel5.1$VCV)
autocorr(ALHmodel6.1$Sol)
autocorr(ALHmodel6.1$VCV)
## summary
summary(ALHmodel6.1)

#genus term not significant
#model simplification (-genus), same as sex model above

##clade differences within Heliconius
#define model
ALH4<-ALH.log~rCBR.log+sim.clade*sex
#create models 
ALHmodel7<-MCMCglmm(ALH4, random=~sp.abb, family=c("gaussian"), 
                    ginverse=list(sp.abb=inv.phylo$Ainv),
                    prior=prior, data=Helionly,
                    nitt=500000,burnin=10000,thin=500)

ALHmodel8<-MCMCglmm(ALH4, random=~sp.abb, family=c("gaussian"), 
                    ginverse=list(sp.abb=inv.phylo$Ainv),
                    prior=prior, data=Helionly,
                    nitt=500000,burnin=10000,thin=500)
#model checks
## check convergence
gelman.diag(mcmc.list(ALHmodel7$Sol,ALHmodel8$Sol))
gelman.diag(mcmc.list(ALHmodel7$VCV,ALHmodel8$VCV))

#visual representation of the output from gelman.diag convergence check
plot(ALHmodel7$VCV)
plot(ALHmodel7$Sol)
plot(ALHmodel8$VCV)
plot(ALHmodel8$Sol)
#check autocorrelation
autocorr(ALHmodel7$Sol)
autocorr(ALHmodel7$VCV)
autocorr(ALHmodel8$Sol)
autocorr(ALHmodel8$VCV)
# summary
summary(ALHmodel7)

#model simplification (-interaction)
#define model 
ALH5<-ALH.log~rCBR.log+sim.clade+sex
#create models
ALHmodel9<-MCMCglmm(ALH5, random=~sp.abb, family=c("gaussian"), 
                    ginverse=list(sp.abb=inv.phylo$Ainv),
                    prior=prior, data=Helionly,
                    nitt=500000,burnin=10000,thin=500)

ALHmodel10<-MCMCglmm(ALH5, random=~sp.abb, family=c("gaussian"), 
                     ginverse=list(sp.abb=inv.phylo$Ainv),
                     prior=prior, data=Helionly,
                     nitt=500000,burnin=10000,thin=500)
#model checks
## check convergence
gelman.diag(mcmc.list(ALHmodel9$Sol,ALHmodel10$Sol))
gelman.diag(mcmc.list(ALHmodel9$VCV,ALHmodel10$VCV))
#visual representation of the output from gelman.diag convergence check
plot(ALHmodel9$VCV)
plot(ALHmodel9$Sol)
plot(ALHmodel10$VCV)
plot(ALHmodel10$Sol)
#check autocorrelation
autocorr(ALHmodel9$Sol)
autocorr(ALHmodel9$VCV)
autocorr(ALHmodel10$Sol)
autocorr(ALHmodel10$VCV)
# summary
summary(ALHmodel9)
#model simplification (-sex)
#define model 
ALH6<-ALH.log~rCBR.log+sim.clade
#create models 
ALHmodel11<-MCMCglmm(ALH6, random=~sp.abb, family=c("gaussian"), 
                     ginverse=list(sp.abb=inv.phylo$Ainv),
                     prior=prior, data=Helionly,
                     nitt=500000,burnin=10000,thin=500)

ALHmodel12<-MCMCglmm(ALH6, random=~sp.abb, family=c("gaussian"), 
                     ginverse=list(sp.abb=inv.phylo$Ainv),
                     prior=prior, data=Helionly,
                     nitt=500000,burnin=10000,thin=500)
#model checks
## check convergence
gelman.diag(mcmc.list(ALHmodel11$Sol,ALHmodel12$Sol))
gelman.diag(mcmc.list(ALHmodel11$VCV,ALHmodel12$VCV))
#visual representation of the output from gelman.diag convergence check
plot(ALHmodel11$VCV)
plot(ALHmodel11$Sol)
plot(ALHmodel12$VCV)
plot(ALHmodel12$Sol)
#check autocorrelation
autocorr(ALHmodel11$Sol)
autocorr(ALHmodel11$VCV)
autocorr(ALHmodel12$Sol)
autocorr(ALHmodel12$VCV)

# summary
summary(ALHmodel11)

#sim.clade is not significant 
#model simplification (-sim.clade)
##ALH null model
#define model 
ALHnull<-ALH.log~rCBR.log
#create models 
ALHnull1<-MCMCglmm(ALHnull, random=~sp.abb, family=c("gaussian"), 
                   ginverse=list(sp.abb=inv.phylo$Ainv),
                   prior=prior, data=data,
                   nitt=500000,burnin=10000,thin=500)

ALHnull2<-MCMCglmm(ALHnull, random=~sp.abb, family=c("gaussian"), 
                   ginverse=list(sp.abb=inv.phylo$Ainv),
                   prior=prior, data=data,
                   nitt=500000,burnin=10000,thin=500)


#model checks
## check convergence
gelman.diag(mcmc.list(ALHnull1$Sol,ALHnull2$Sol))
gelman.diag(mcmc.list(ALHnull1$VCV,ALHnull2$VCV))
#visual representation of the output from gelman.diag convergence check
plot(ALHnull1$VCV)
plot(ALHnull1$Sol)
plot(ALHnull2$VCV)
plot(ALHnull2$Sol)
#check autocorrelation
autocorr(ALHnull1$Sol)
autocorr(ALHnull1$VCV)
autocorr(ALHnull2$Sol)
autocorr(ALHnull2$VCV)
# summary
summary(ALHnull1)
#extract phylogenetic signal lambda (mean)
lambdaALH<-ALHnull1$VCV[,'sp.abb']/(ALHnull1$VCV[,'sp.abb']+ALHnull1$VCV[,'units'])
mean(lambdaALH)

#compare significance of the phylogenetic signal 
#run model without phylogeny 
ALHnophylo<-MCMCglmm(ALHnull, random=~units, family=c("gaussian"), 
                   prior=prior, data=data,
                   nitt=500000,burnin=10000,thin=500)

# Compute log-likelihoods
ALHlogL_full <- ALHnull1$DIC  # DIC for full model
ALHlogL_null <- ALHnophylo$DIC  # DIC for no phylo model
# Compute likelihood ratio test statistic
ALH_LRT_stat <- ALHlogL_null - ALHlogL_full
ALH_p_value <- pchisq(ALH_LRT_stat, df = 1, lower.tail = FALSE)
# Print results
cat("Likelihood Ratio Test statistic:", ALH_LRT_stat, "\n")
cat("p-value:", ALH_p_value, "\n")

#next test
ALH_lambda_CI <- quantile(lambdaALH, probs = c(0.025, 0.975))
# Print results
cat("Lambda mean:", mean(lambdaALH), "\n")
cat("95% credible interval:", ALH_lambda_CI, "\n")

#ALH species level comparisons 
#used previously defined null model
#create model 
ALHnull3<-MCMCglmm(ALHnull, random=~sp.abb, family=c("gaussian"), 
                   ginverse=list(sp.abb=inv.phylo$Ainv),
                   prior=prior, data=filtered_data,
                   nitt=500000,burnin=10000,thin=500,pr=T)

ALHnull4<-MCMCglmm(ALHnull, random=~sp.abb, family=c("gaussian"), 
                   ginverse=list(sp.abb=inv.phylo$Ainv),
                   prior=prior, data=filtered_data,
                   nitt=500000,burnin=10000,thin=500,pr=T)
#model checks
# check convergence
gelman.diag(mcmc.list(ALHnull3$Sol,ALHnull4$Sol))
gelman.diag(mcmc.list(ALHnull3$VCV,ALHnull4$VCV))
#visual representation of the output from gelman.diag convergence check
plot(ALHnull3$VCV)
plot(ALHnull3$Sol)
plot(ALHnull4$VCV)
plot(ALHnull4$Sol)
#check autocorrelation
autocorr(ALHnull3$Sol)
autocorr(ALHnull3$VCV)
autocorr(ALHnull4$Sol)
autocorr(ALHnull4$VCV)
# summary
summary(ALHnull3)

##pairwise comparison of ALH volumes 
# Extract posterior distribution for the random effects
ALH_random_effects <- ALHnull3$Sol
# Extract the random effects for species
ALH_random_effects_species <- ALH_random_effects[, grep("sp.abb", colnames(ALH_random_effects))]
# Convert to mcmc object
ALH_random_effects_species_mcmc <- as.mcmc(ALH_random_effects_species)
# pairwise comparison and apply the function with a higher threshold to account for false discoveries
ALH_comparison_results <- pairwise_comparisons(ALH_random_effects_species_mcmc, threshold = 0.99)
#extract comparisons 
ALH_comparison_matrix <- ALH_comparison_results$comparisons
#save as csv.
write.csv(ALH_comparison_matrix, "ALH pairwise_comparisons all values.csv")

##GL  Heliconius vs non-Heliconius
#define model
GL1<-GL.log~rCBR.log+Heliconius*sex
#create models 
GLmodel1<-MCMCglmm(GL1, random=~sp.abb, family=c("gaussian"), 
                   ginverse=list(sp.abb=inv.phylo$Ainv),
                   prior=prior, data=data,
                   nitt=500000,burnin=10000,thin=500)

GLmodel2<-MCMCglmm(GL1, random=~sp.abb, family=c("gaussian"), 
                   ginverse=list(sp.abb=inv.phylo$Ainv),
                   prior=prior, data=data,
                   nitt=500000,burnin=10000,thin=500)
#model checks
# check convergence
gelman.diag(mcmc.list(GLmodel1$Sol,GLmodel2$Sol))
gelman.diag(mcmc.list(GLmodel1$VCV,GLmodel2$VCV))
#visual representation of the output from gelman.diag convergence check
plot(GLmodel1$VCV)
plot(GLmodel1$Sol)
plot(GLmodel2$VCV)
plot(GLmodel2$Sol)
#check autocorrelation
autocorr(GLmodel1$Sol)
autocorr(GLmodel1$VCV)
autocorr(GLmodel2$Sol)
autocorr(GLmodel2$VCV)
# summary
summary(GLmodel1)

#model simplification (-interaction)
#define model
GL1.1<-GL.log~rCBR.log+Heliconius+sex
#create models 
GLmodel1.1<-MCMCglmm(GL1.1, random=~sp.abb, family=c("gaussian"), 
                   ginverse=list(sp.abb=inv.phylo$Ainv),
                   prior=prior, data=data,
                   nitt=500000,burnin=10000,thin=500)

GLmodel1.2<-MCMCglmm(GL1.1, random=~sp.abb, family=c("gaussian"), 
                   ginverse=list(sp.abb=inv.phylo$Ainv),
                   prior=prior, data=data,
                   nitt=500000,burnin=10000,thin=500)
#model checks
## check convergence
gelman.diag(mcmc.list(GLmodel1.1$Sol,GLmodel1.2$Sol))
gelman.diag(mcmc.list(GLmodel1.1$VCV,GLmodel1.2$VCV))
#visual representation of the output from gelman.diag convergence check
plot(GLmodel1.1$VCV)
plot(GLmodel1.1$Sol)
plot(GLmodel1.2$VCV)
plot(GLmodel1.2$Sol)
#check autocorrelation
autocorr(GLmodel1.1$Sol)
autocorr(GLmodel1.1$VCV)
autocorr(GLmodel1.2$Sol)
autocorr(GLmodel1.2$VCV)
#summary
summary(GLmodel1.1)

#model simplification (-sex)
#define model
GL1.2<-GL.log~rCBR.log+Heliconius
#create models 
GLmodel1.3<-MCMCglmm(GL1.2, random=~sp.abb, family=c("gaussian"), 
                     ginverse=list(sp.abb=inv.phylo$Ainv),
                     prior=prior, data=data,
                     nitt=500000,burnin=10000,thin=500)

GLmodel1.4<-MCMCglmm(GL1.2, random=~sp.abb, family=c("gaussian"), 
                     ginverse=list(sp.abb=inv.phylo$Ainv),
                     prior=prior, data=data,
                     nitt=500000,burnin=10000,thin=500)
#model checks
## check convergence
gelman.diag(mcmc.list(GLmodel1.3$Sol,GLmodel1.4$Sol))
gelman.diag(mcmc.list(GLmodel1.3$VCV,GLmodel1.4$VCV))
#visual representation of the output from gelman.diag convergence check
plot(GLmodel1.3$VCV)
plot(GLmodel1.3$Sol)
plot(GLmodel1.4$VCV)
plot(GLmodel1.4$Sol)
#check autocorrelation
autocorr(GLmodel1.3$Sol)
autocorr(GLmodel1.3$VCV)
autocorr(GLmodel1.4$Sol)
autocorr(GLmodel1.4$VCV)
#summary
summary(GLmodel1.4)
#Heliconius term not significant 
#model simplification (-Heliconius)
#null model GL
#define model
GLnull<-GL.log~rCBR.log
#create models
GLnull1<-MCMCglmm(GLnull, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)

GLnull2<-MCMCglmm(GLnull, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)
#model checks
## check convergence
gelman.diag(mcmc.list(GLnull1$Sol,GLnull2$Sol))
gelman.diag(mcmc.list(GLnull1$VCV,GLnull2$VCV))
#visual representation of the output from gelman.diag convergence check
plot(GLnull1$VCV)
plot(GLnull1$Sol)
plot(GLnull2$VCV)
plot(GLnull2$Sol)
#check autocorrelation
autocorr(GLnull1$Sol)
autocorr(GLnull1$VCV)
autocorr(GLnull2$Sol)
autocorr(GLnull2$VCV)
#summary
summary(GLnull2)

##GL genera within Heliconiini
#define model 
GL2<-GL.log~rCBR.log+genus*sex
#create models 
GLmodel3<-MCMCglmm(GL2, random=~sp.abb, family=c("gaussian"), 
                   ginverse=list(sp.abb=inv.phylo$Ainv),
                   prior=prior, data=data,
                   nitt=500000,burnin=10000,thin=500)

GLmodel4<-MCMCglmm(GL2, random=~sp.abb, family=c("gaussian"), 
                   ginverse=list(sp.abb=inv.phylo$Ainv),
                   prior=prior, data=data,
                   nitt=500000,burnin=10000,thin=500)
#model checks
# check convergence
gelman.diag(mcmc.list(GLmodel3$Sol,GLmodel4$Sol))
gelman.diag(mcmc.list(GLmodel3$VCV,GLmodel4$VCV))
#visual representation of the output from gelman.diag convergence check
plot(GLmodel3$VCV)
plot(GLmodel3$Sol)
plot(GLmodel4$VCV)
plot(GLmodel4$Sol)
#check autocorrelation
autocorr(GLmodel3$Sol)
autocorr(GLmodel3$VCV)
autocorr(GLmodel4$Sol)
autocorr(GLmodel4$VCV)
#summary 
summary(GLmodel3)

#model simplification (-interaction)
#define model 
GL2.1<-GL.log~rCBR.log+genus+sex
#create models 
GLmodel3.1<-MCMCglmm(GL2.1, random=~sp.abb, family=c("gaussian"), 
                   ginverse=list(sp.abb=inv.phylo$Ainv),
                   prior=prior, data=data,
                   nitt=500000,burnin=10000,thin=500)

GLmodel4.1<-MCMCglmm(GL2.1, random=~sp.abb, family=c("gaussian"), 
                   ginverse=list(sp.abb=inv.phylo$Ainv),
                   prior=prior, data=data,
                   nitt=500000,burnin=10000,thin=500)
#model checks
# check convergence
gelman.diag(mcmc.list(GLmodel3.1$Sol,GLmodel4.1$Sol))
gelman.diag(mcmc.list(GLmodel3.1$VCV,GLmodel4.1$VCV))
#visual representation of the output from gelman.diag convergence check
plot(GLmodel3.1$VCV)
plot(GLmodel3.1$Sol)
plot(GLmodel4.1$VCV)
plot(GLmodel4.1$Sol)
#check autocorrelation
autocorr(GLmodel3.1$Sol)
autocorr(GLmodel3.1$VCV)
autocorr(GLmodel4.1$Sol)
autocorr(GLmodel4.1$VCV)
#summary 
summary(GLmodel3.1)

#model simplification (-genus)
#define model 
GL4<-GL.log~rCBR.log+sex
#create models 
GLmodel7<-MCMCglmm(GL4, random=~sp.abb, family=c("gaussian"), 
                   ginverse=list(sp.abb=inv.phylo$Ainv),
                   prior=prior, data=data,
                   nitt=500000,burnin=10000,thin=500)

GLmodel8<-MCMCglmm(GL4, random=~sp.abb, family=c("gaussian"), 
                   ginverse=list(sp.abb=inv.phylo$Ainv),
                   prior=prior, data=data,
                   nitt=500000,burnin=10000,thin=500)
#model checks
# check convergence
gelman.diag(mcmc.list(GLmodel7$Sol,GLmodel8$Sol))
gelman.diag(mcmc.list(GLmodel7$VCV,GLmodel8$VCV))
#visual representation of the output from gelman.diag convergence check
plot(GLmodel7$VCV)
plot(GLmodel7$Sol)
plot(GLmodel8$VCV)
plot(GLmodel8$Sol)
#check autocorrelation
autocorr(GLmodel7$Sol)
autocorr(GLmodel7$VCV)
autocorr(GLmodel8$Sol)
autocorr(GLmodel8$VCV)
#summary 
summary(GLmodel7)

#no sex difference 

##GL clade within Heliconius 
#define model 
GL3<-GL.log~rCBR.log+sim.clade*sex
#create models 
GLmodel5<-MCMCglmm(GL3, random=~sp.abb, family=c("gaussian"), 
                   ginverse=list(sp.abb=inv.phylo$Ainv),
                   prior=prior, data=Helionly,
                   nitt=500000,burnin=10000,thin=500)

GLmodel6<-MCMCglmm(GL3, random=~sp.abb, family=c("gaussian"), 
                   ginverse=list(sp.abb=inv.phylo$Ainv),
                   prior=prior, data=Helionly,
                   nitt=500000,burnin=10000,thin=500)
#model checks
# check convergence
gelman.diag(mcmc.list(GLmodel5$Sol,GLmodel6$Sol))
gelman.diag(mcmc.list(GLmodel5$VCV,GLmodel6$VCV))
#visual representation of the output from gelman.diag convergence check
plot(GLmodel5$VCV)
plot(GLmodel5$Sol)
plot(GLmodel6$VCV)
plot(GLmodel6$Sol)
#check autocorrelation
autocorr(GLmodel5$Sol)
autocorr(GLmodel5$VCV)
autocorr(GLmodel6$Sol)
autocorr(GLmodel6$VCV)
#summary 
summary(GLmodel5)

#model simplification (-interaction)
#define model 
GL3.1<-GL.log~rCBR.log+sim.clade+sex
#create models 
GLmodel5.1<-MCMCglmm(GL3.1, random=~sp.abb, family=c("gaussian"), 
                     ginverse=list(sp.abb=inv.phylo$Ainv),
                     prior=prior, data=Helionly,
                     nitt=500000,burnin=10000,thin=500)

GLmodel6.1<-MCMCglmm(GL3.1, random=~sp.abb, family=c("gaussian"), 
                     ginverse=list(sp.abb=inv.phylo$Ainv),
                     prior=prior, data=Helionly,
                     nitt=500000,burnin=10000,thin=500)
#model checks
# check convergence
gelman.diag(mcmc.list(GLmodel5.1$Sol,GLmodel6.1$Sol))
gelman.diag(mcmc.list(GLmodel5.1$VCV,GLmodel6.1$VCV))
#visual representation of the output from gelman.diag convergence check
plot(GLmodel5.1$VCV)
plot(GLmodel5.1$Sol)
plot(GLmodel6.1$VCV)
plot(GLmodel6.1$Sol)
#check autocorrelation
autocorr(GLmodel5.1$Sol)
autocorr(GLmodel5.1$VCV)
autocorr(GLmodel6.1$Sol)
autocorr(GLmodel6.1$VCV)
#summary 
summary(GLmodel5.1)

#model simplification (-sim.clade), same as sex model above 


#extract phylogenetic signal lambda (mean)
lambdaGL<-GLnull2$VCV[,'sp.abb']/(GLnull2$VCV[,'sp.abb']+GLnull2$VCV[,'units'])
mean(lambdaGL)

#compare significance of the phylogenetic signal 
#run model without phylogeny 
GLnophylo<-MCMCglmm(GLnull, random=~units, family=c("gaussian"), 
                     prior=prior, data=data,
                     nitt=500000,burnin=10000,thin=500)

# Compute log-likelihoods
GLlogL_full <- GLnull2$DIC  # DIC for full model
GLlogL_null <- GLnophylo$DIC  # DIC for no phylo model
# Compute likelihood ratio test statistic
GL_LRT_stat <- GLlogL_null - GLlogL_full
GL_p_value <- pchisq(GL_LRT_stat, df = 1, lower.tail = FALSE)
# Print results
cat("Likelihood Ratio Test statistic:", GL_LRT_stat, "\n")
cat("p-value:", GL_p_value, "\n")

#next test
GL_lambda_CI <- quantile(lambdaGL, probs = c(0.025, 0.975))
# Print results
cat("Lambda mean:", mean(lambdaGL), "\n")
cat("95% credible interval:", GL_lambda_CI, "\n")

##ALH species level comparisons 
#used previously defined null model
#create model 
GLnull3<-MCMCglmm(GLnull, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=filtered_data,
                  nitt=500000,burnin=10000,thin=500,pr=T)
GLnull4<-MCMCglmm(GLnull, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=filtered_data,
                  nitt=500000,burnin=10000,thin=500,pr=T)
#model checks
# check convergence
gelman.diag(mcmc.list(GLnull3$Sol,GLnull4$Sol))
gelman.diag(mcmc.list(GLnull3$VCV,GLnull4$VCV))
#visual representation of the output from gelman.diag convergence check
plot(GLnull3$VCV)
plot(GLnull3$Sol)
plot(GLnull4$VCV)
plot(GLnull4$Sol)
#check autocorrelation
autocorr(GLnull3$Sol)
autocorr(GLnull3$VCV)
autocorr(GLnull4$Sol)
autocorr(GLnull4$VCV)
#summary 
summary(GLnull3)
# Extract posterior samples for the random effects
GL_random_effects <- GLnull3$Sol
# Extract the random effects for species
GL_random_effects_species <- GL_random_effects[, grep("sp.abb", colnames(GL_random_effects))]
# Convert to mcmc object
GL_random_effects_species_mcmc <- as.mcmc(GL_random_effects_species)
# pairwise comparison and apply the function with a higher threshold to account for false discoveries
GL_comparison_results <- pairwise_comparisons(GL_random_effects_species_mcmc, threshold = 0.99)
#extract comparisons
GL_comparison_matrix <- GL_comparison_results$comparisons
#save as csv.
write.csv(GL_comparison_matrix, "GL pairwise_comparisons all values.csv")

##GL vs ALH Heliconiini vs non-Heliconiini
#define model
comp1<-GL.log~ALH.log+Heliconius*sex
#create models 
comp1.1<-MCMCglmm(comp1, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)

comp1.2<-MCMCglmm(comp1, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)
#model checks
# check convergence
gelman.diag(mcmc.list(comp1.1$Sol,comp1.2$Sol))
gelman.diag(mcmc.list(comp1.1$VCV,comp1.2$VCV))
#visual representation of the output from gelman.diag convergence check
plot(comp1.1$VCV)
plot(comp1.1$Sol)
plot(comp1.2$VCV)
plot(comp1.2$Sol)
#check autocorrelation
autocorr(comp1.1$Sol)
autocorr(comp1.1$VCV)
autocorr(comp1.2$Sol)
autocorr(comp1.2$VCV)
#summary
summary(comp1.2)

#model simplification (-interaction)
#define model 
comp1.1<-GL.log~ALH.log+Heliconius+sex
#create models 
comp1.3<-MCMCglmm(comp1.1, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)

comp1.4<-MCMCglmm(comp1.1, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)
#model checks
# check convergence
gelman.diag(mcmc.list(comp1.3$Sol,comp1.4$Sol))
gelman.diag(mcmc.list(comp1.3$VCV,comp1.4$VCV))
#visual representation of the output from gelman.diag convergence check
plot(comp1.3$VCV)
plot(comp1.3$Sol)
plot(comp1.4$VCV)
plot(comp1.4$Sol)
#check autocorrelation
autocorr(comp1.3$Sol)
autocorr(comp1.3$VCV)
autocorr(comp1.4$Sol)
autocorr(comp1.4$VCV)
#summary
summary(comp1.3)

#model simplification (-Heliconius)
#define model 
comp4<-GL.log~ALH.log+sex
#create models
comp4.1<-MCMCglmm(comp4, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)

comp4.2<-MCMCglmm(comp4, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)
#model checks
# check convergence
gelman.diag(mcmc.list(comp4.1$Sol,comp4.2$Sol))
gelman.diag(mcmc.list(comp4.1$VCV,comp4.2$VCV))
#visual representation of the output from gelman.diag convergence check
plot(comp4.1$VCV)
plot(comp4.1$Sol)
plot(comp4.2$VCV)
plot(comp4.2$Sol)
#check autocorrelation
autocorr(comp4.1$Sol)
autocorr(comp4.1$VCV)
autocorr(comp4.2$Sol)
autocorr(comp4.2$VCV)

#summary
summary(comp4.2)
#sex is significant

##GL vs ALH genera within Heliconiini
#define model 
comp2<-GL.log~ALH.log+genus*sex
#create models
comp2.1<-MCMCglmm(comp2, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)

comp2.2<-MCMCglmm(comp2, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)
#model checks
# check convergence
gelman.diag(mcmc.list(comp2.1$Sol,comp2.2$Sol))
gelman.diag(mcmc.list(comp2.1$VCV,comp2.2$VCV))
#visual representation of the output from gelman.diag convergence check
plot(comp2.1$VCV)
plot(comp2.1$Sol)
plot(comp2.2$VCV)
plot(comp2.2$Sol)

#check autocorrelation
autocorr(comp2.1$Sol)
autocorr(comp2.1$VCV)
autocorr(comp2.2$Sol)
autocorr(comp2.2$VCV)

#summary
summary(comp2.1)

#model simplification (-interaction)
#define model 
comp2.1<-GL.log~ALH.log+genus+sex
#create models
comp2.3<-MCMCglmm(comp2.1, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)

comp2.4<-MCMCglmm(comp2.1, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)
#model checks
# check convergence
gelman.diag(mcmc.list(comp2.3$Sol,comp2.4$Sol))
gelman.diag(mcmc.list(comp2.3$VCV,comp2.4$VCV))

#visual representation of the output from gelman.diag convergence check
plot(comp2.3$VCV)
plot(comp2.3$Sol)
plot(comp2.4$VCV)
plot(comp2.4$Sol)

#check autocorrelation
autocorr(comp2.3$Sol)
autocorr(comp2.3$VCV)
autocorr(comp2.4$Sol)
autocorr(comp2.4$VCV)

#summary
summary(comp2.3)

#model simplification (-genus), same as sex model above 

##GL vs ALH clades within Heliconius
#define model
comp3<-GL.log~ALH.log+sim.clade*sex
#create models
comp3.1<-MCMCglmm(comp3, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=Helionly,
                  nitt=500000,burnin=10000,thin=500)

comp3.2<-MCMCglmm(comp3, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=Helionly,
                  nitt=500000,burnin=10000,thin=500)
#model checks
# check convergence
gelman.diag(mcmc.list(comp3.1$Sol,comp3.2$Sol))
gelman.diag(mcmc.list(comp3.1$VCV,comp3.2$VCV))
#visual representation of the output from gelman.diag convergence check
plot(comp3.1$VCV)
plot(comp3.1$Sol)
plot(comp3.2$VCV)
plot(comp3.2$Sol)
#check autocorrelation
autocorr(comp3.1$Sol)
autocorr(comp3.1$VCV)
autocorr(comp3.2$Sol)
autocorr(comp3.2$VCV)
#summary 
summary(comp3.1)

#model simplification (-interaction)
#define model
comp5<-GL.log~ALH.log+sim.clade+sex
#create models 
comp5.1<-MCMCglmm(comp5, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=Helionly,
                  nitt=500000,burnin=10000,thin=500)

comp5.2<-MCMCglmm(comp5, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=Helionly,
                  nitt=500000,burnin=10000,thin=500)
#model checks
# check convergence
gelman.diag(mcmc.list(comp5.1$Sol,comp5.2$Sol))
gelman.diag(mcmc.list(comp5.1$VCV,comp5.2$VCV))
#visual representation of the output from gelman.diag convergence check
plot(comp5.1$VCV)
plot(comp5.1$Sol)
plot(comp5.2$VCV)
plot(comp5.2$Sol)
#check autocorrelation
autocorr(comp5.1$Sol)
autocorr(comp5.1$VCV)
autocorr(comp5.2$Sol)
autocorr(comp5.2$VCV)
#summary 
summary(comp5.1)
#model simplification (-sim.clade), same as sex model above

##null model
#define model
comp6<-GL.log~ALH.log
#create models
comp6.1<-MCMCglmm(comp6, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)

comp6.2<-MCMCglmm(comp6, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=data,
                  nitt=500000,burnin=10000,thin=500)
#model checks
# check convergence
gelman.diag(mcmc.list(comp6.1$Sol,comp6.2$Sol))
gelman.diag(mcmc.list(comp6.1$VCV,comp6.2$VCV))
#visual representation of the output from gelman.diag convergence check
plot(comp6.1$VCV)
plot(comp6.1$Sol)
plot(comp6.2$VCV)
plot(comp6.2$Sol)
#check autocorrelation
autocorr(comp6.1$Sol)
autocorr(comp6.1$VCV)
autocorr(comp6.2$Sol)
autocorr(comp6.2$VCV)
#summary
summary(comp6.2)

#extract phylogenetic signal lambda (mean)
lambdacomp<-comp6.2$VCV[,'sp.abb']/(comp6.2$VCV[,'sp.abb']+comp6.2$VCV[,'units'])
mean(lambdacomp)

#run model without phylogeny 
bothnophylo<-MCMCglmm(comp6, random=~units, family=c("gaussian"), 
                    prior=prior, data=data,
                    nitt=500000,burnin=10000,thin=500)

# Compute log-likelihoods
bothlogL_full <- comp6.2$DIC  # DIC for full model
bothlogL_null <- bothnophylo$DIC  # DIC for no phylo model
# Compute likelihood ratio test statistic
both_LRT_stat <- bothlogL_null - bothlogL_full
both_p_value <- pchisq(both_LRT_stat, df = 1, lower.tail = FALSE)
# Print results
cat("Likelihood Ratio Test statistic:", both_LRT_stat, "\n")
cat("p-value:", both_p_value, "\n")

#next test
both_lambda_CI <- quantile(lambdacomp, probs = c(0.025, 0.975))
# Print results
cat("Lambda mean:", mean(lambdacomp), "\n")
cat("95% credible interval:", both_lambda_CI, "\n")

#GL~ALH species level comparisons 
#used previously defined null model
#create model 
comp6.3<-MCMCglmm(comp6, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=filtered_data,
                  nitt=500000,burnin=10000,thin=500,pr=T)
comp6.4<-MCMCglmm(comp6, random=~sp.abb, family=c("gaussian"), 
                  ginverse=list(sp.abb=inv.phylo$Ainv),
                  prior=prior, data=filtered_data,
                  nitt=500000,burnin=10000,thin=500,pr=T)
#model checks
# check convergence
gelman.diag(mcmc.list(comp6.3$Sol,comp6.4$Sol))
gelman.diag(mcmc.list(comp6.3$VCV,comp6.4$VCV))
#visual representation of the output from gelman.diag convergence check
plot(comp6.3$VCV)
plot(comp6.3$Sol)
plot(comp6.4$VCV)
plot(comp6.4$Sol)
#check autocorrelation
autocorr(comp6.3$Sol)
autocorr(comp6.3$VCV)
autocorr(comp6.4$Sol)
autocorr(comp6.4$VCV)
#summary
summary(comp6.3)

# Extract posterior samples for the random effects
comp_random_effects <- comp6.3$Sol
# Extract the random effects for species
comp_random_effects_species <- comp_random_effects[, grep("sp.abb", colnames(comp_random_effects))]
# Convert to mcmc object
comp_random_effects_species_mcmc <- as.mcmc(comp_random_effects_species)
# pairwise comparison and apply the function with a higher threshold to account for false discoveries
comp_comparison_results <- pairwise_comparisons(comp_random_effects_species_mcmc, threshold = 0.99)
#extract comparisons 
comp_comparison_matrix <- comp_comparison_results$comparisons
#save as csv.
write.csv(comp_comparison_matrix, "GL~ALH pairwise_comparisons all values.csv")