#this code re-runs the mixed-effect model to estimate ERR for all sampels in Tom's original data, 
#
#code mostly by Tom Crellen (tomcrellen@gmail.com /  thomas.crellen@bdi.ox.ac.uk) and Martin Walker (mwalker@rvc.ac.uk)
#some edits by James Cotton (jc17@sanger.ac.uk)
#
#Script for analysis of S. mansoni learance data
#exact results will differ from published due to Markov error in the MCMC, but probably a tiny discrepancy
#
library(MCMCglmm)
library(ggplot2)
library(scales)
library(plyr)

ug <- read.table('uganda.smansoni_minimised.model.txt',header=T)

source("margeff.R")
source("ERR-dist-func.R")
source("intercept_mixed.R")

#Create treatment intensity categories
ug$treat_intensity <- with(ug, ifelse(school== "Kocoge", "Low",
                                    ifelse(school == "Bukoba", "Medium",
                                           ifelse(school == "Bukagabo", "Medium","High")))) 
ug$treat_intensity <- as.factor(ug$treat_intensity)

#Create weight categories
ug$weight_cat <- with(ug, ifelse(weight >26, "High",
                               ifelse(weight <22, "Low", "Medium")))
ug$weight_cat <- as.factor(ug$weight_cat)

#Create binary variable for before or after treatment
ug$treatment <- as.character(ug$time)
ug$treatment[ug$time=="BT"] <- 0
ug$treatment[ug$time=="AT"] <- 1
ug$treatment <- as.numeric(ug$treatment)

#Ensure count is numeric
ug$count<- as.numeric(as.character(ug$count))

#Remove NAs (here from weight)
ug <- ug[!is.na(ug$weight),]

#Create dataframe for Sm (note: these DO NOT include "uninfected individuals")
sm <- subset(ug,  treated == 1 & count != "NA" & days_total_post > 0 & days_total_pre > 0 & sm_pos == 1)

#Drop unused factors (schools) and reset factor order
sm$school <- factor(sm$school)
sm$school <- factor(sm$school, 
             levels = c("Bwondha", "Bugoto", "Musubi", 
             "Kocoge", "Bukoba", "Bukagabo"))
sm$treat_intensity <- factor(sm$treat_intensity, 
                      levels = c("Low", "Medium", "High"))

#Add co-infection with STH covariate for Sm dataframe
sm$sth <- with(sm, ifelse((hk_pos + as_pos + tr_pos) > 0, 1, 0))

# fixed effects for Smansoni
full_form_sm <- formula(count ~ treatment*treat_intensity+treatment*weight_cat + treatment*sex + treatment*sth)
# random effects
ranef <- formula( ~ us(1+treatment):UID + us(1+treatment):school)
# observation level random effects which permits overdispersion among egg counts sharing a common set of covariates
overdisp = formula(~idh(1):units)

# prior distributions 
prior <- list(R = list(V = 1, nu = 0.002),
              G = list(G1 = list(V = diag(c(1,1)), nu=0.002), 
                       G2 = list(V = diag(c(1,1)), nu=0.002)))

#######################
#MCMC Model
#######################

# run Markov chain for "nitt" iterations with a "burnin"
nitt <- 100000; burnin <- 2500
# fit the model to Smansoni
fullmix_sm <- MCMCglmm(fixed = full_form_sm, random = ranef, prior = prior, 
                    rcov = overdisp, family = "poisson", data = sm, 
                    thin = 10, verbose = TRUE, nitt = nitt, 
                    burnin = burnin, pr=TRUE, pl=TRUE)

out_UID_sm <- margeff(modlist = list(fullmix_sm), data = sm, nm = "treatment", marg = "UID") 


egg_bl_sm <- margeff_intercept(modlist = list(fullmix_sm), data = sm, nm = "treatment", marg = "UID") 
egg_bl_sm$post <- egg_bl_sm$mn-(egg_bl_sm$mn * out_UID_sm$mn)


out_sm_ERRdist_IND <- ERRdist(modlist = list(fullmix_sm), data = sm, nm = "treatment", marg = "UID")



#order individuals by ascending posterior mean
dfplot <- out_UID_sm[order(out_UID_sm$mn),]
# setting x-axis position to be 1:num of individuals
dfplot$x <- seq(1, nrow(dfplot))

extra <- unique(sm[,c(1,2,4,3)])
colnames(extra)[2] <- "id"

dfplot_extra <- join(dfplot,extra,by="id")

Sample_lists <- read.delim("Sample_classification.txt")
PID_to_CID  <- read.delim("Sample_IDs.txt")
sample.data <- data.frame(patient_ID=Sample_lists$patient_ID,classification=Sample_lists$classification)
sample.data.CID <- join(sample.data,PID_to_CID)



good_ids <- dfplot_extra$id[dfplot_extra$CID %in% unique(sample.data.CID$CID[sample.data$classification == "Good_clearers"])]
poor_ids <- dfplot_extra$id[dfplot_extra$CID %in% unique(sample.data.CID$CID[sample.data$classification == "Pre-treatment" | sample.data$classification == "Post-treatment" ])]


good_ERRdist <- out_sm_ERRdist_IND[out_sm_ERRdist_IND$id %in% good_ids,]
poor_ERRdist <- out_sm_ERRdist_IND[out_sm_ERRdist_IND$id %in% poor_ids,]

#generate means and medians per sample category for each MCMC iteration.
#MEAN
mean(colMeans(good_ERRdist[,-ncol(good_ERRdist)]))
HPDinterval(as.mcmc(colMeans(good_ERRdist[,-ncol(good_ERRdist)])))       
mean(colMeans(poor_ERRdist[,-ncol(poor_ERRdist)]))
HPDinterval(as.mcmc(colMeans(poor_ERRdist[,-ncol(poor_ERRdist)])))          
#MEDIAN
mean(apply(good_ERRdist[,-ncol(good_ERRdist)],2,median))
HPDinterval(as.mcmc(apply(good_ERRdist[,-ncol(good_ERRdist)],2,median)))        
mean(apply(poor_ERRdist[,-ncol(poor_ERRdist)],2,median))
HPDinterval(as.mcmc(apply(poor_ERRdist[,-ncol(poor_ERRdist)],2,median)))         

## this data frame contains mean, median and 95% CIs for each individual 
all_samples_ERR_summary <- dfplot_extra[dfplot_extra$CID %in% unique(sample.data.CID$CID),]
