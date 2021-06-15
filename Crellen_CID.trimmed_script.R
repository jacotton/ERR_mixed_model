library(MCMCglmm)
source("/Users/tc13/Desktop/margeff.R")

#Read in data
d <- read.table("uganda_input_model.txt", header = TRUE)

#Create District Column
d$district <- ifelse(d$school == "Kocoge", "Tororo", "Mayuge")

#Create "Treatment Intensity" Column
d$treat_intensity <- with(d, ifelse(school== "Kocoge", "Low",
                                    ifelse(school == "Bukoba", "Medium",
                                           ifelse(school == "Bukagabo", "Medium","High")))) 
d$treat_intensity <- as.factor(d$treat_intensity)

#Create weight categories
d$weight_cat <- with(d, ifelse(weight >26, "High",
                               ifelse(weight <22, "Low", "Medium")))
d$weight_cat <- as.factor(d$weight_cat)

#Create binary variable for before or after treatment
d$treatment <- as.character(d$time)
d$treatment[d$time=="BT"] <- 0
d$treatment[d$time=="AT"] <- 1
d$treatment <- as.numeric(d$treatment)

# fixed effects
full_form <- formula(count ~ treatment*treat_intensity+treatment*weight_cat + treatment*sex)

# random effects
ranef <- formula( ~ us(1+treatment):numID +  us(1+treatment):school)

# observation level random effects which permits overdispersion among egg counts sharing a common set of covariates
overdisp = formula(~idh(1):units)

# prior distributions 
prior <- list(R = list(V = 1, nu = 0.002),
               G = list(G1 = list(V = diag(c(1,1)), nu=0.002), 
                        G2 = list(V = diag(c(1,1)), nu=0.002)))

# run Markov chain for "nitt" iterations with a "burnin"
nitt <- 100000; burnin <- 2500

# fit the model
fullmix <- MCMCglmm(fixed = full_form, random = ranef, prior = prior, 
                    rcov = overdisp, family = "poisson", data = d, 
                    thin = 10, verbose = FALSE, nitt = 
                    nitt, burnin = burnin, pr=TRUE, pl=TRUE)

# marginalize over treatment intensity group
out_treat <- margeff(modlist = list(fullmix), data = d, nm = "treatment", marg = "treat_intensity") 

# marginalize over school
out_sc <- margeff(modlist = list(fullmix), data = d, nm = "treatment", marg = "school") 

# marginalize over weight category
out_wt <- margeff(modlist = list(fullmix), data = d, nm = "treatment", marg = "weight_cat") 

# "marginalize" (it's not really marginalizing for individuals as 
## ERR is defined as an individual level effect) over individuals
out_id <- margeff(modlist = list(fullmix), data = d, nm = "treatment", marg = "numID") 

#Exctract intercept values - ie baseline egg counts from model
pre_id <- margeff_intercept(modlist = list(fullmix), data = d, marg = "numID") 

