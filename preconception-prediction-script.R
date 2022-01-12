#------------------------------------------------------------------------------#
#   Sample Script: Model Selection, Model Fitting, and Model Assessment        #
#------------------------------------------------------------------------------#

library(mice)
source("calcHUM-IPW.R")
source("selection-and-prediction-functions.R")
source("data-cleaning-and-bootstrapping-functions.R")

### Note:
  # The data file used in this script (containing de-identified pregnancy outcome 
  # data from the EAGeR trial) is available upon reasonable request. This script
  # produces: (1) the predicted preconception risk profiles for all individuals
  # in the EAGeR dataset; (2) the optimistic HUM; and (3) the optimism-corrected
  # HUM.
  # To obtain the bootstrapped CIs for (2) and (3), re-run this script on 50
  # bootstrapped datasets to obtain the bootstrap distribution for the optimistic
  # HUM and optimism-corrected HUM.
###


#-------------------------------------------------------------------------------
# Reading in the Dataset
#-------------------------------------------------------------------------------

# Reading in the original dataset 
#   Replace with bootstrapped dataset for bootstrapped CI
ppt <- read.csv(file="eager.csv")


#-------------------------------------------------------------------------------
# Two-Stage Imputation of Covariates & Pregnancy Outcomes
#-------------------------------------------------------------------------------

# Storing study id numbers with constraints on the outcome imputation
#   Imputed outcome must be either pre-term or full-term birth
id.lb <- ppt$study.id[which(ppt$final.status=="live birth")] 
#   Women with a known pregnancy, but unknown final outcome status
id.no.w <- ppt$study.id[which(ppt$final.status=="withdrawal" & ppt$clin.conf.preg==1)]
#   Women with an unknown pregnancy status
id.no.ppt <- ppt$study.id[which(ppt$clin.conf.preg==0 & ppt$cycle.clin.preg < 6)]

# Creating additional variables to assist with two-stage outcome imputation
ppt$i.1 <- rep(0, nrow(ppt)); ppt$i.2 <- rep(0, nrow(ppt)); ppt$i.3 <- rep(0, nrow(ppt))
ppt$i.4 <- rep(0, nrow(ppt)); ppt$i.5 <- rep(0, nrow(ppt)); ppt$i.6 <- rep(0, nrow(ppt))
ppt$r.1 <- rep(NA, nrow(ppt)); ppt$r.2 <- rep(NA, nrow(ppt)); ppt$r.3 <- rep(NA, nrow(ppt))
ppt$r.4 <- rep(NA, nrow(ppt)); ppt$r.5 <- rep(NA, nrow(ppt)); ppt$r.6 <- rep(NA, nrow(ppt))
ppt$cycle.clin.preg.impute <- ppt$cycle.clin.preg
ppt$cycle.clin.preg.impute[ppt$study.id %in% id.no.ppt] <- NA
ppt$cycle.clin.preg.impute[ppt$cycle.clin.preg.impute==6 & ppt$clin.conf.preg==0] <- 7
ppt$r.1 <- ppt$cycle.clin.preg.impute
ppt$i.1[which(ppt$cycle.clin.preg.impute==1 & ppt$clin.conf.preg==1)] <- 1
ppt$r.2 <- ppt$cycle.clin.preg.impute
ppt$r.2[which(ppt$cycle.clin.preg.impute==1)] <- NA
ppt$i.2[which((ppt$cycle.clin.preg.impute==2 & ppt$clin.conf.preg==1) | (ppt$cycle.clin.preg==1 & ppt$clin.conf.preg==0))] <- 1
ppt$r.3 <- ppt$cycle.clin.preg.impute
ppt$r.3[which(ppt$cycle.clin.preg.impute %in% c(1, 2))] <- NA
ppt$i.3[which((ppt$cycle.clin.preg.impute==3 & ppt$clin.conf.preg==1) | (ppt$cycle.clin.preg==2 & ppt$clin.conf.preg==0))] <- 1
ppt$r.4 <- ppt$cycle.clin.preg.impute
ppt$r.4[which(ppt$cycle.clin.preg.impute %in% c(1, 2, 3))] <- NA
ppt$i.4[which((ppt$cycle.clin.preg.impute==4 & ppt$clin.conf.preg==1) | (ppt$cycle.clin.preg==3 & ppt$clin.conf.preg==0))] <- 1
ppt$r.5 <- ppt$cycle.clin.preg.impute
ppt$r.5[which(ppt$cycle.clin.preg.impute %in% c(1, 2, 3, 4))] <- NA
ppt$i.5[which((ppt$cycle.clin.preg.impute==5 & ppt$clin.conf.preg==1) | (ppt$cycle.clin.preg==4 & ppt$clin.conf.preg==0))] <- 1
ppt$r.6 <- ppt$cycle.clin.preg.impute
ppt$r.6[which(ppt$cycle.clin.preg.impute %in% c(1, 2, 3, 4, 5))] <- NA
ppt$i.6[which((ppt$cycle.clin.preg.impute==6 & ppt$clin.conf.preg==1) | (ppt$cycle.clin.preg==5 & ppt$clin.conf.preg==0))] <- 1

# Designating factor variables with appropriate levels
ppt$income <- as.factor(ppt$income)
ppt$smoking.ind <- as.factor(ppt$smoking.ind)
ppt$hypertension <- as.factor(ppt$hypertension)
ppt$clin.conf.preg <- as.factor(ppt$clin.conf.preg)
ppt$white <- ifelse(ppt$race=="", NA, ifelse(ppt$race=="Non-Hispanic white", 1, 0))
ppt$white <- as.factor(ppt$white)
ppt$r.1 <- as.factor(ppt$r.1)
ppt$r.2 <- as.factor(ppt$r.2)
ppt$r.3 <- as.factor(ppt$r.3)
ppt$r.4 <- as.factor(ppt$r.4)
ppt$r.5 <- as.factor(ppt$r.5)
ppt$r.6 <- as.factor(ppt$r.6)
ppt$i.1 <- as.factor(ppt$i.1)
ppt$i.2 <- as.factor(ppt$i.2)
ppt$i.3 <- as.factor(ppt$i.3)
ppt$i.4 <- as.factor(ppt$i.4)
ppt$i.5 <- as.factor(ppt$i.5)
ppt$i.6 <- as.factor(ppt$i.6)
ppt$preterm.birth <- as.factor(ppt$preterm.birth)
ppt$final.status[ppt$cycle.clin.preg==6 & ppt$clin.conf.preg==0] <- "efuwp"
ppt$final.status[ppt$study.id %in% union(union(id.lb, id.no.w), id.no.ppt)] <- NA
ppt$final.status <- droplevels(as.factor(ppt$final.status))
ppt$valid <- ifelse(is.na(ppt$final.status) | 
                      (ppt$final.status=="withdrawal") |
                      (ppt$final.status=="efuwp" & ppt$cycle.clin.preg != 6), 0, 1)

# Creating additional factor variable to assist with imputations with constraints
ppt$final.no.withdrawal <- ppt$final.status
ppt$final.no.withdrawal[ppt$final.status=="withdrawal" | ppt$final.status=="efuwp"] <- NA
ppt$final.no.withdrawal <- droplevels(as.factor(ppt$final.no.withdrawal))

# Specifying imputation method 
meth <- c("", "", "", "", "norm", "logreg", "", "", "polyreg", "sample", "sample", "sample", "", "", "sample", "sample", 
          "polr", "sample", "sample", "sample", "sample", "sample", "logreg", "logreg", "logreg", "", "", 
          "", "", "", "", "", "", "polyreg", "polyreg", "polyreg", "polyreg", "polyreg", "polyreg", 
          "~I(as.numeric(as.character(i.1))*as.numeric(as.character(r.1)) + as.numeric(as.character(i.2))*as.numeric(as.character(r.2)) + as.numeric(as.character(i.3))*as.numeric(as.character(r.3)) + as.numeric(as.character(i.4))*as.numeric(as.character(r.4)) + as.numeric(as.character(i.5))*as.numeric(as.character(r.5)) + as.numeric(as.character(i.6))*as.numeric(as.character(r.6)))", "", "polyreg")

# Creating the initial predictor matrix
init <- mice(ppt[, c(2, 3, 4, 7, 9, 19, 30, 31, 33, 42, 49, 53, 54, 57, 58, 59, 60, 65, 68, 69, 71, 76, 94, 109, 
                     110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127)],
             max = 0, seed = 1123)
vs <- init$visitSequence[c(1, 2, 3, 4, 7, 8, 10, 11, 12, 14, 15, 16, 18, 19, 20, 21, 22, 6, 5, 17, 13, 24, 25, 
                           26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 23, 9, 42, 41)]
pred.mat <- init$predictorMatrix

# Updating the specification of the predictors in the imputation model
pred.mat[c('final.no.withdrawal'), ] <- pred.mat[c('final.status'), ]
pred.mat[c('preterm.birth'), ] <- pred.mat[c('final.status'),]
pred.mat[c('cycle.clin.preg.impute'), ] <- pred.mat[c('cycle.clin.preg'), ]
pred.mat[, c('cycle.clin.preg.impute')] <- pred.mat[, c('cycle.clin.preg')]
pred.mat[c('study.id', 'aspirin', 'age', 'n.previous.losses', 'n.live.birth.NIH', 'cycles.trying.prior.rand',
           'income', 'parous', 'cycle.clin.preg', 'clin.conf.preg', 'i.1', 'i.2', 'i.3', 'i.4', 'i.5', 'i.6', 'valid'),] <- 0
pred.mat[c(1, 2, 3, 4, 6, 7, 8, 9, 13, 14, 17, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42),
         c('SHBG', 'Albumin', 'insulin2', 'waist2hip','baseline.hsCRP', 'T3', 'ApoB', 'C.Peptide', 'Uric.Acid', 'Thyroglobulin')] <- 0
pred.mat[, c('study.id', 'preterm.birth', 'cycle.clin.preg', 'clin.conf.preg', 'i.1', 'i.2', 'i.3', 'i.4', 
             'i.5', 'i.6', 'r.1', 'r.2', 'r.3', 'r.4', 'r.5', 'r.6', 'valid', 'final.no.withdrawal')] <- 0
pred.mat[c('preterm.birth', 'r.1', 'r.2', 'r.3', 'r.4', 'r.5', 'r.6'), 
         c('aspirin', 'age', 'n.previous.losses', 'BMI', 'n.live.birth.NIH', 'cycles.trying.prior.rand',
           'income', 'parous', 'edu', 'smoking.ind', 'hypertension', 'white')] <- 1
pred.mat[c('preterm.birth'), c('cycle.clin.preg.impute')] <- 1
pred.mat[c('r.1', 'r.2', 'r.3', 'r.4', 'r.5', 'r.6'), c('final.status')] <- 1
pred.mat[c('r.1', 'r.2', 'r.3', 'r.4', 'r.5', 'r.6'), c('cycle.clin.preg.impute')] <- 0
pred.mat[, c('parous')] <- 0

# Running the multiple imputation
imp <- mice(ppt[, c(2, 3, 4, 7, 9, 19, 30, 31, 33, 42, 49, 53, 54, 57, 58, 59, 60, 65, 68, 69, 71, 76, 94, 109, 
                    110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127)],
            m = 10, pred = pred.mat, visitSequence = vs, method = meth, seed = 1123)

# Returning the imputed datasets
ppt.1 <- clean_imputed_data(complete(imp, 1))
ppt.2 <- clean_imputed_data(complete(imp, 2))
ppt.3 <- clean_imputed_data(complete(imp, 3))
ppt.4 <- clean_imputed_data(complete(imp, 4))
ppt.5 <- clean_imputed_data(complete(imp, 5))
ppt.6 <- clean_imputed_data(complete(imp, 6))
ppt.7 <- clean_imputed_dataa(complete(imp, 7))
ppt.8 <- clean_imputed_data(complete(imp, 8))
ppt.9 <- clean_imputed_data(complete(imp, 9))
ppt.10 <- clean_imputed_data(complete(imp, 10))


#-------------------------------------------------------------------------------
# Model Selection 
#-------------------------------------------------------------------------------

# Selecting the first-stage, second-stage, and IPW models
stacked <- rbind(ppt.1, ppt.2, ppt.3, ppt.4, ppt.5, ppt.6, ppt.7, ppt.8, ppt.9, ppt.10)
selected.models <- model_selection(stacked, first.stage.y = c("cycle.clin.preg", "clin.conf.preg"),
                                   first.stage.var = list("max"=c("aspirin", "age", "BMI", "cycles.trying.prior.rand", "income",
                                                            "parous", "college", "smoking.ind", "hypertension", "n.previous.losses", "white"),
                                                          "min"=c("n.previous.losses", "white")),
                                   second.stage.selection = "clin.conf.preg.impute", 
                                   second.stage.y = "final.status",
                                   second.stage.var = list("max"=c("aspirin", "age", "BMI", "cycles.trying.prior.rand", "income",
                                                                   "parous", "college", "smoking.ind", "hypertension", "n.previous.losses", "white"),
                                                           "min"=c("n.previous.losses", "white")),
                                   ipw.y = "valid", 
                                   ipw.var = list("max"=c("aspirin", "age", "BMI", "cycles.trying.prior.rand", "income",
                                                          "parous", "college", "smoking.ind", "hypertension", "n.previous.losses", "white"),
                                                  "min"=c("1")),
                                   M=10)


#-------------------------------------------------------------------------------
# Model Fitting
#-------------------------------------------------------------------------------

# Fitting the models on the imputed data: first stage (time to event)
est.coef <- fit_stage_1(selected.models[[1]], ppt.1, "cycle.clin.preg", "clin.conf.preg")
est.coef <- rbind(est.coef, fit_stage_1(selected.models.old[[1]], ppt.2, "cycle.clin.preg", "clin.conf.preg"))
est.coef <- rbind(est.coef, fit_stage_1(selected.models.old[[1]], ppt.3, "cycle.clin.preg", "clin.conf.preg"))
est.coef <- rbind(est.coef, fit_stage_1(selected.models.old[[1]], ppt.4, "cycle.clin.preg", "clin.conf.preg"))
est.coef <- rbind(est.coef, fit_stage_1(selected.models.old[[1]], ppt.5, "cycle.clin.preg", "clin.conf.preg"))
est.coef <- rbind(est.coef, fit_stage_1(selected.models.old[[1]], ppt.6, "cycle.clin.preg", "clin.conf.preg"))
est.coef <- rbind(est.coef, fit_stage_1(selected.models.old[[1]], ppt.7, "cycle.clin.preg", "clin.conf.preg"))
est.coef <- rbind(est.coef, fit_stage_1(selected.models.old[[1]], ppt.8, "cycle.clin.preg", "clin.conf.preg"))
est.coef <- rbind(est.coef, fit_stage_1(selected.models.old[[1]], ppt.9, "cycle.clin.preg", "clin.conf.preg"))
est.coef <- rbind(est.coef, fit_stage_1(selected.models.old[[1]], ppt.10, "cycle.clin.preg", "clin.conf.preg"))
ppt.coef <- apply(est.coef[seq(1, nrow(est.coef), by=2), ], 2, mean)
var.w <- apply(est.coef[seq(2, nrow(est.coef), by=2), ], 2, mean)
var.b <- split(est.coef[seq(1, nrow(est.coef), by=2), ], 1:nrow(est.coef[seq(1, nrow(est.coef), by=2), ]))
var.b <- lapply(var.b, FUN=function(x)as.matrix(x - ppt.coef) %*% t(as.matrix(x-ppt.coef)))
var.b <- Reduce("+", var.b)/9
var.ppt.coef <- var.w + (1 + 1/10)*diag(var.b)

# Fitting the models on the imputed data: second stage (pregnancy outcome)
est.coef <- fit_stage_2(selected.models[[2]], ppt.1, "clin.conf.preg.impute")
est.coef <- rbind(est.coef, fit_stage_2(selected.models[[2]], ppt.2, "clin.conf.preg.impute"))
est.coef <- rbind(est.coef, fit_stage_2(selected.models[[2]], ppt.3, "clin.conf.preg.impute"))
est.coef <- rbind(est.coef, fit_stage_2(selected.models[[2]], ppt.4, "clin.conf.preg.impute"))
est.coef <- rbind(est.coef, fit_stage_2(selected.models[[2]], ppt.5, "clin.conf.preg.impute"))
est.coef <- rbind(est.coef, fit_stage_2(selected.models[[2]], ppt.6, "clin.conf.preg.impute"))
est.coef <- rbind(est.coef, fit_stage_2(selected.models[[2]], ppt.7, "clin.conf.preg.impute"))
est.coef <- rbind(est.coef, fit_stage_2(selected.models[[2]], ppt.8, "clin.conf.preg.impute"))
est.coef <- rbind(est.coef, fit_stage_2(selected.models[[2]], ppt.9, "clin.conf.preg.impute"))
est.coef <- rbind(est.coef, fit_stage_2(selected.models[[2]], ppt.10, "clin.conf.preg.impute"))
out.coef <- matrix(NA, nrow=2, ncol=ncol(est.coef))
out.coef[1,] <- apply(est.coef[seq(1, nrow(est.coef), by=4), ], 2, mean)
out.coef[2,] <- apply(est.coef[seq(2, nrow(est.coef), by=4), ], 2, mean)
colnames(out.coef) <- colnames(est.coef)
var.w <- c(apply(est.coef[seq(3, nrow(est.coef), by=4), ], 2, mean), apply(est.coef[seq(4, nrow(est.coef), by=4), ], 2, mean))
var.b <- split(as.data.frame(est.coef[c(seq(1, nrow(est.coef), by=4), seq(2, nrow(est.coef), by=4)), ]), rep(seq(1, nrow(est.coef[seq(1, nrow(est.coef), by=4),])), by=2))
var.b <- lapply(var.b, FUN=function(x)unlist(c(x[1,], x[2,])))
mean.vec <- unlist(c(out.coef[1,], out.coef[2,]))
var.b <- lapply(var.b, FUN=function(x)as.matrix(x - mean.vec) %*% t(as.matrix(x-mean.vec)))
var.b <- Reduce("+", var.b)/9
var.out.coef <- var.w + (1 + 1/10)*diag(var.b)
var.out.coef <- matrix(var.out.coef, nrow=2, byrow=T)
colnames(var.out.coef) <- colnames(est.coef)
sd.out.coef <- sqrt(var.out.coef)

# Fitting the models on the imputed data: verification IP weights
est.coef <- fit_stage_3(selected.models[[3]], ppt.1)
est.coef <- rbind(est.coef, fit_stage_3(selected.models[[3]], ppt.2))
est.coef <- rbind(est.coef, fit_stage_3(selected.models[[3]], ppt.3))
est.coef <- rbind(est.coef, fit_stage_3(selected.models[[3]], ppt.4))
est.coef <- rbind(est.coef, fit_stage_3(selected.models[[3]], ppt.5))
est.coef <- rbind(est.coef, fit_stage_3(selected.models[[3]], ppt.6))
est.coef <- rbind(est.coef, fit_stage_3(selected.models[[3]], ppt.7))
est.coef <- rbind(est.coef, fit_stage_3(selected.models[[3]], ppt.8))
est.coef <- rbind(est.coef, fit_stage_3(selected.models[[3]], ppt.9))
est.coef <- rbind(est.coef, fit_stage_3(selected.models[[3]], ppt.10))
ipw.coef <- apply(est.coef[seq(1, nrow(est.coef), by=2), ], 2, mean)
var.w <- apply(est.coef[seq(2, nrow(est.coef), by=2), ], 2, mean)
var.b <- split(est.coef[seq(1, nrow(est.coef), by=2), ], 1:nrow(est.coef[seq(1, nrow(est.coef), by=2), ]))
var.b <- lapply(var.b, FUN=function(x)as.matrix(x - ipw.coef) %*% t(as.matrix(x - ipw.coef)))
var.b <- Reduce("+", var.b)/9
var.ipw.coef <- var.w + (1 + 1/10)*diag(var.b)


#-------------------------------------------------------------------------------
# Obtaining Preconception Risk Profiles
#-------------------------------------------------------------------------------

# Cleaning the pregnancy dataset
ppt <- clean_ppt_data(ppt)

# Predicted probabilities: clinically confirmed pregnancy
ppt.prob <- pred_stage_1(ppt, selected.models[[1]], ppt.coef)

# Predicted probabilities: pregnancy outcomes
out.prob <- pred_stage_2(ppt, selected.models[[2]], out.coef)

# Predicted probabilities: validation weights
ipw <- pred_r(ppt, selected.models[[3]], ipw.coef)

# Calculating pregnancy outcome risk profiles
profiles <- t(apply(cbind(ppt.prob[, ncol(ppt.prob)], out.prob[, c(2, 3, 1)]), 1, risk_profile))
profiles <- as.data.frame(profiles)
names(profiles) <- c("prob.no.ccp", "prob.clin.loss", "prob.preterm", "prob.fullterm")

# Re-ordering risk profiles
profiles <- profiles[, c(2, 1, 4, 3)]


#-------------------------------------------------------------------------------
# Optimistic, Verification-Bias-Adjusted HUM
#-------------------------------------------------------------------------------

# True outcome status
y <- ppt$final.status
y[ppt$final.status=="withdrawal"] <- NA
y[ppt$final.status=="live birth"] <- NA
y[ppt$final.status=="efuwp" & ppt$cycle.clin.preg != 6] <- NA
y <- droplevels(y)
y <- as.numeric(y)

# Removing observations with NAs
rows <- c(which(!complete.cases(profiles)), which(is.na(y)), which(is.na(ipw)))
y <- y[-rows]
profiles <- profiles[-rows,]
risk.mat <- as.matrix(cbind(profiles, y))
weights <- ipw[-rows]

# Calculating the optimistic IPW-adjusted HUM
opt.hum <- calcHUM(risk.mat, weights); print(opt.hum)


#-------------------------------------------------------------------------------
# Optimism Correction
#-------------------------------------------------------------------------------

# Estimating optimism in the model performance assessment
optimism <- data.frame()
for (i in 1:50){
  boot <- ppt[sample.int(nrow(ppt), replace = TRUE), ]
  temp <- hum_boot(data, boot, selected.models)
  optimism <- rbind(optimism, temp)
}

# Optimism-adjusted HUM
adj.hum <- opt.hum - mean(optimism$boot.hum - optimism$orig.hum)
