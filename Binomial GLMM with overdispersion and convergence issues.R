######################################################################
######                Fitting a GLMM, trouble shooting          ######
######                Dealing with overdispersion, optimisers   ######
######################################################################

##################################################################################################
### Section 1: Information about the script -----
##################################################################################################


############ A. Script information      -----------------------

# Script Content:         This script works with data from Hischier et al (2023) to illusrate an appraoch to fitting a generalized linear mixed model and trouble-shooting issues with the model. The data come from an experiment of seed germination of 10 plant species with 3 different treatments at two elevations (1400 m.a.s.l. and 2000 m.a.s.l.) on Calanda Mountain (Switzerland)

# Data required:         HischierEtal_DataForAnalysis_230629.txt

# Authors of this script: Jake Alexander 
# Main contact (and email): jake.alexander@usys.ethz.ch
# Date created:           2023.09.15
# Last changes:           2023.09.15


############ B. Abbreviations of factor levels used in the dataset (if relevant) --------------------------------


#### Lowland species:
# Pla.med         Plantago media
# Sca.col         Scabiosa columbaria
# Med.lup         Medicago lupulina
# Bra.pin         Brachypodium pinnatum
# Bro.ere         Bromus erectus

#### Highland species:
# Pla.atr         Plantago atrata
# Sca.luc         Scabiosa lucida
# Lot.alp         Lotus alpinus
# Poa.alp         Poa alpinus
# Ses.cae         Sesleria caerulea

# Low elevation site (at 1400 m.a.s.l.)         Nesselboden
# High elevation site (at 2000 ma.s.l.)         Calanda


##################################################################################################
### Section 2: Package installation, functions and data import -----
##################################################################################################


############ Working directory ############

setwd("/Users/alexander/polybox2/Shared/Paper Hischier et al./Data and code")


############ Packages ############

# Packages ----------------------------------------------------------------

library("car")   					#for Anova on glmer objects (establishment models)
library("MuMIn")					#for AICc
library("lme4")						#for glmer
library("dfoptim")        # to compare different optimisers

############ Functions ############

#function to check a glm(m) for overdispersion
# see https://github.com/lme4/lme4/issues/220
overdisp_fun <- function(model) {
  ## number of variance parameters in
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# -------------------------------------------------------------------------


############ Read-in data ############

seed.germination.max_data <- read.table("HischierEtal_DataForAnalysis_230629.txt", header=T)
head(seed.germination.max_data)


##################################################################################################
### Section 3: Data analysis - establishment (main analysis) -----
##################################################################################################

#this binomial model has the germination rate as the response variable, expressed as "successes" (number of germinants observed) and "failures" (non-germinating seeds). The explanatory variables are the experimental site (2 levels, high or low), the origin of the focal species (2 levels, high or low) and the experimental treatment (3 levels, which aren't important here). Observations were made on 10 species, with 3 replicates per species within an experimental plot. We therefore fit random effects for species and plot to account for the nestedness of observations within these factors

germination.focalorigin_glmm <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * focalorigin 
                                      + (1|species) + (1|plot), family = "binomial", data = seed.germination.max_data)
summary(germination.focalorigin_glmm) 
overdisp_fun(germination.focalorigin_glmm) #the dispersion ratio is 2.12 > 1.2, so model is overdispersed. To deal with this, we refit the model with an individual-level random effect

germination.focalorigin_glmm <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * focalorigin 
                                      + (1|species) + (1|plot)+ (1|toothpick), family = "binomial", data = seed.germination.max_data)
# this yields a warning message that the model has failed to converge
# We can use "allFit" to compare different optimisers to see if this helps. For more information: (?convergence ; https://joshua-nugent.github.io/allFit/)


diff_optims <- allFit(germination.focalorigin_glmm, maxfun=1e5) #function fits the same model with multiple optimisers to see which ones work
diff_optims_OK <- diff_optims[sapply(diff_optims, is, "merMod")]
lapply(diff_optims_OK, function(x) x@optinfo$conv$lme4$messages) #this shows us the warning messages for each method. We see that bobyqa and two others converged 

#refit our model using the bobyqa optimiser
germination.focalorigin_glmm2 <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * focalorigin 
                                       + (1|species) + (1|plot) + (1|toothpick), family = "binomial", data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa")) 
isSingular(germination.focalorigin_glmm2, tol = 1e-4)  # to test model for singularity - none
overdisp_fun(germination.focalorigin_glmm2) #dispersion is now <1.2
Anova(germination.focalorigin_glmm2) #anova table (from "car" library) with P-values for each variable


## this led me to use the "bobyqa" optimiser for all models

germination.focalorigin_glmm <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * focalorigin 
                                      + (1|species) + (1|plot), family = "binomial", data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa"))
summary(germination.focalorigin_glmm) 
overdisp_fun(germination.focalorigin_glmm)
isSingular(germination.focalorigin_glmm, tol = 1e-4)                             # to test model for singularity
#model is overdispersed, add observation level random effect

germination.focalorigin_glmm2 <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * focalorigin 
                                       + (1|species) + (1|plot) + (1|toothpick), family = "binomial", data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa"))
summary(germination.focalorigin_glmm2)  # converges
isSingular(germination.focalorigin_glmm2, tol = 1e-4)                             # to test model for singularity
overdisp_fun(germination.focalorigin_glmm2)
Anova(germination.focalorigin_glmm2)

# Analysis of Deviance Table (Type II Wald chisquare tests)

# Response: cbind(seedlings, meanseeds - seedlings)
# Chisq Df Pr(>Chisq)    
# treatment                  45.5562  2  1.281e-10 ***
# site                       15.6150  1  7.764e-05 ***
# focalorigin                 2.3832  1    0.12265    
# treatment:site             26.4347  2  1.819e-06 ***
# treatment:focalorigin      28.8728  2  5.375e-07 ***
# site:focalorigin            1.0275  1    0.31074    
# treatment:site:focalorigin  7.7974  2    0.02027 *  

AICc(germination.focalorigin_glmm2)	#5699.066 AICc (AIC corrected for small sample size) of the final model


