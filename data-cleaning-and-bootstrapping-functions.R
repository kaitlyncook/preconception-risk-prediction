#------------------------------------------------------------------------------#
#     Functions for Data Cleaning and Conducting the Bootstrapped Analysis     #
#------------------------------------------------------------------------------#

library(mice)
source("calcHUM-IPW.R")
source("selection-and-prediction-functions.R")


#-------------------------------------------------------------------------------
# Data Cleaning Functions
#-------------------------------------------------------------------------------

### Cleaning the analysis dataset
  # data                dataframe         analysis dataset     
clean_ppt_data <- function(data){
  # Creating a college indicator
  data$college <- ifelse(data$edu %in% c("College degree", "Postgraduate degree"), 1, 0)
  
  # Cleaning indicator variables
  data$white <- as.numeric(as.character(data$white))
  data$hypertension <- as.numeric(as.character(data$hypertension))
  data$smoking.ind <- as.numeric(as.character(data$smoking.ind)) 
  data$income <- as.factor(data$income)
  
  # Cleaning final outcome status
  data$final.status <- relevel(data$final.status, ref="full-term birth")
  
  return(data)
}

### Cleaning the imputed datasets
  # data                dataframe         imputed dataset 
clean_imputed_data <- function(data){
  # Updating both the imputed time and imputed event variable indicators
  data$clin.conf.preg.impute <- ifelse(data$cycle.clin.preg.impute == 7, 0, 1)
  data$cycle.clin.preg.impute[data$cycle.clin.preg.impute==7] <- 6
  
  # Cleaning indicator variables
  data$white <- as.numeric(as.character(data$white))
  data$hypertension <- as.numeric(as.character(data$hypertension))
  data$smoking.ind <- as.numeric(as.character(data$smoking.ind))
  data$income <- as.factor(data$income)
  
  # Determining final outcome status
  data$final.status <- data$final.no.withdrawal
  data$final.status <- as.character(data$final.status)
  data$final.status[data$clin.conf.preg.impute==0] <- "efuwp"
  data$final.status[data$study.id %in% id.lb] <- ifelse(data$preterm.birth[data$study.id %in% id.lb]==1, "preterm birth", "full-term birth")
  data$final.status <- as.factor(data$final.status)
  data$final.status <- relevel(data$final.status, ref="full-term birth")
  
  # Creating a college indicator
  data$college <- ifelse(data$edu %in% c("College degree", "Postgraduate degree"), 1, 0)
  
  return(data)
}


#-------------------------------------------------------------------------------
# Function for Bootstrapped Optimism Analysis
#-------------------------------------------------------------------------------

### Obtaining a single estimate of the optimism in the HUM
  # data                dataframe         original analysis dataset 
  # boot                dataframe         bootstrapped analysis dataset
  # selected.models     list              selected model formulations for the first-stage, second-stage, and IPW models
hum_boot <- function(data, boot, selected.models){
  
  # Conducting the two-stage imputation procedure on the bootstrapped data
  meth <- c("", "", "", "", "norm", "logreg", "", "", "polyreg", "sample", "sample", "sample", "", "", "sample", "sample", 
            "polr", "sample", "sample", "sample", "sample", "sample", "logreg", "logreg", "logreg", "", "", 
            "", "", "", "", "", "", "polyreg", "polyreg", "polyreg", "polyreg", "polyreg", "polyreg", 
            "~I(as.numeric(as.character(i.1))*as.numeric(as.character(r.1)) + as.numeric(as.character(i.2))*as.numeric(as.character(r.2)) + as.numeric(as.character(i.3))*as.numeric(as.character(r.3)) + as.numeric(as.character(i.4))*as.numeric(as.character(r.4)) + as.numeric(as.character(i.5))*as.numeric(as.character(r.5)) + as.numeric(as.character(i.6))*as.numeric(as.character(r.6)))", "", "polyreg")
  init <- mice(boot[, c(2, 3, 4, 7, 9, 19, 30, 31, 33, 42, 49, 53, 54, 57, 58, 59, 60, 65, 68, 69, 71, 76, 94, 109, 
                        110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127)],
               max=0)
  vs <- init$visitSequence[c(1, 2, 3, 4, 7, 8, 10, 11, 12, 14, 15, 16, 18, 19, 20, 21, 22, 6, 5, 17, 13, 24, 25, 
                             26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 23, 9, 42, 41)]
  pred.mat <- init$predictorMatrix
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
  imp <- mice(boot[, c(2, 3, 4, 7, 9, 19, 30, 31, 33, 42, 49, 53, 54, 57, 58, 59, 60, 65, 68, 69, 71, 76, 94, 109, 
                       110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127)],
              m = 10, pred = pred.mat, visitSequence = vs, method = meth)
  
  # Cleaning the imputed bootstrapped datasets
  ppt.1 <- clean_imputed_data(complete(imp, 1))
  ppt.2 <- clean_imputed_data(complete(imp, 2))
  ppt.3 <- clean_imputed_data(complete(imp, 3))
  ppt.4 <- clean_imputed_data(complete(imp, 4))
  ppt.5 <- clean_imputed_data(complete(imp, 5))
  ppt.6 <- clean_imputed_data(complete(imp, 6))
  ppt.7 <- clean_imputed_data(complete(imp, 7))
  ppt.8 <- clean_imputed_data(complete(imp, 8))
  ppt.9 <- clean_imputed_data(complete(imp, 9))
  ppt.10 <- clean_imputed_data(complete(imp, 10))
  
  # Fitting the models on the imputed bootstrapped data: first stage (time to event)
  est.coef <- fit_stage_1(selected.models[[1]], ppt.1, "cycle.clin.preg", "clin.conf.preg")
  est.coef <- rbind(est.coef, fit_stage_1(selected.models[[1]], ppt.2, "cycle.clin.preg", "clin.conf.preg"))
  est.coef <- rbind(est.coef, fit_stage_1(selected.models[[1]], ppt.3, "cycle.clin.preg", "clin.conf.preg"))
  est.coef <- rbind(est.coef, fit_stage_1(selected.models[[1]], ppt.4, "cycle.clin.preg", "clin.conf.preg"))
  est.coef <- rbind(est.coef, fit_stage_1(selected.models[[1]], ppt.5, "cycle.clin.preg", "clin.conf.preg"))
  est.coef <- rbind(est.coef, fit_stage_1(selected.models[[1]], ppt.6, "cycle.clin.preg", "clin.conf.preg"))
  est.coef <- rbind(est.coef, fit_stage_1(selected.models[[1]], ppt.7, "cycle.clin.preg", "clin.conf.preg"))
  est.coef <- rbind(est.coef, fit_stage_1(selected.models[[1]], ppt.8, "cycle.clin.preg", "clin.conf.preg"))
  est.coef <- rbind(est.coef, fit_stage_1(selected.models[[1]], ppt.9, "cycle.clin.preg", "clin.conf.preg"))
  est.coef <- rbind(est.coef, fit_stage_1(selected.models[[1]], ppt.10, "cycle.clin.preg", "clin.conf.preg"))
  ppt.coef <- apply(est.coef[seq(1, nrow(est.coef), by=2), ], 2, mean)
  
  # Fitting the models on the imputed bootstrapped data: second stage (pregnancy outcome)
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
  
  # Fitting the models on the imputed bootstrapped data: verification IP weights
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
  
  # In-sample model evaluation: bootstrapped analysis dataset
  boot <- clean_ppt_data(boot)
  ppt.prob <- pred_stage_1(boot, selected.models[[1]], ppt.coef)
  out.prob <- pred_stage_2(boot, selected.models[[2]], out.coef)
  ipw <- pred_r(boot, selected.models[[3]], ipw.coef)
  profiles <- t(apply(cbind(ppt.prob[, ncol(ppt.prob)], out.prob[, c(2, 3, 1)]), 1, risk_profile))
  profiles <- as.data.frame(profiles)
  names(profiles) <- c("prob.no.ccp", "prob.clin.loss", "prob.preterm", "prob.fullterm")
  profiles <- profiles[, c(2, 1, 4, 3)]
  y <- boot$final.status
  y[boot$final.status=="withdrawal"] <- NA
  y[boot$final.status=="live birth"] <- NA
  y[boot$final.status=="efuwp" & boot$cycle.clin.preg != 6] <- NA
  y <- droplevels(y)
  y <- as.numeric(y)
  rows <- c(which(!complete.cases(profiles)), which(is.na(y)), which(is.na(ipw)))
  y <- y[-rows]
  profiles <- profiles[-rows,]
  risk.mat <- as.matrix(cbind(profiles, y))
  weights <- ipw[-rows]
  boot.hum <- calcHUM(risk.mat, weights)
  
  # Out-of-sample model evaluation: original analysis dataset
  data <- clean_ppt_data(data)
  
  ppt.prob <- pred_stage_1(data, selected.models[[1]], ppt.coef)
  out.prob <- pred_stage_2(data, selected.models[[2]], out.coef)
  ipw <- pred_r(data, selected.models[[3]], ipw.coef)
  profiles <- t(apply(cbind(ppt.prob[, ncol(ppt.prob)], out.prob[, c(2, 3, 1)]), 1, risk_profile))
  profiles <- as.data.frame(profiles)
  names(profiles) <- c("prob.no.ccp", "prob.clin.loss", "prob.preterm", "prob.fullterm")
  profiles <- profiles[, c(2, 1, 4, 3)]
  y <- data$final.status
  y[data$final.status=="withdrawal"] <- NA
  y[data$final.status=="live birth"] <- NA
  y[data$final.status=="efuwp" & data$cycle.clin.preg != 6] <- NA
  y <- droplevels(y)
  y <- as.numeric(y)
  rows <- c(which(!complete.cases(profiles)), which(is.na(y)), which(is.na(ipw)))
  y <- y[-rows]
  profiles <- profiles[-rows,]
  risk.mat <- as.matrix(cbind(profiles, y))
  weights <- ipw[-rows]
  orig.hum <- calcHUM(risk.mat, weights)
  
  return(data.frame(boot.hum=boot.hum, orig.hum=orig.hum))
}

