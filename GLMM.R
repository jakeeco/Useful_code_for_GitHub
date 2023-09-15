


setwd("/Users/alexander/polybox2/Shared/Paper Hischier et al./Data and code")

seed.germination.max_data <- read.table("HischierEtal_DataForAnalysis_230629.txt", header=T)
head(seed.germination.max_data)


# First 


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

AICc(germination.focalorigin_glmm2)	#5699.066
aovtab <- as.data.frame(Anova(germination.focalorigin_glmm2))
aovtab$Chisq <- round(aovtab$Chisq, 2)
aovtab[,3] <- round(aovtab[,3], 3)
aovtab[,3][aovtab[,3] < 0.001] <- "<0.001"

aovtab.origin <- aovtab



## ASIDE: fit of this model using default optimiser returned a convergence warning. I use "allFit" to compare different optimisers: (?convergence ; https://joshua-nugent.github.io/allFit/)
## this led me to use the "bobyqa" optimiser for all models
germination.focalorigin_glmm2 <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * focalorigin 
                                       + (1|species) + (1|plot) + (1|toothpick), family = "binomial", data = seed.germination.max_data)
library(dfoptim)
diff_optims <- allFit(germination.focalorigin_glmm2, maxfun=1e5)
diff_optims_OK <- diff_optims[sapply(diff_optims, is, "merMod")]
lapply(diff_optims_OK, function(x) x@optinfo$conv$lme4$messages)
germination.focalorigin_glmm2 <- glmer(cbind(seedlings,meanseeds-seedlings) ~ treatment * site * focalorigin 
                                       + (1|species) + (1|plot) + (1|toothpick), family = "binomial", data = seed.germination.max_data, control = glmerControl(optimizer = "bobyqa"))
Anova(germination.focalorigin_glmm2)

