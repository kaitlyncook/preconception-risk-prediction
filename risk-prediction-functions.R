#------------------------------------------------------------------------------#
#           Functions for Model Selection, Risk Prediction and                 #
#                   Validation for Pregnancy Prediction                        #
#------------------------------------------------------------------------------#

library(discSurv)
library(nnet)
library(glmnet)
library(MASS)
library(mice)
source("calcHUM-IPW.R")

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
# Model Selection Function
#-------------------------------------------------------------------------------

### Determining the selected prediction and weighting models
  # stacked                 dataframe         single data frame containing all M imputed datasets
  # first.stage.y           character         vector containing the column names for the first-stage outcome: (time, delta)
  # first.stage.var         list              list containing (1) a vector of column names for the maximal set of variables in the first-stage model
  #                                                           (2) a vector of column names for the minimal set of variables in the first-stage model
  # second.stage.selection  character         column name of the selection event for the second-stage outcome/model
  # second.stage.y          character         column name for the second-stage outcome
  # second.stage.var        list              list containing (1) a vector of column names for the maximal set of variables in the second-stage model
  #                                                           (2) a vector of column names for the minimal set of variables in the second-stage model
  # ipw.y                   vector            vector containing either (i) the column names of R_1 and R_2 or (ii) the column name of R
  # ipw.var                 list              list containing either (i) two further lists, each with two elements corresponding to the column names for the maximal and minimal sets of variables for the IPW models for R_1 and R_2
  #                                                                  (ii) two vectors of column names for the maximal and minimal sets of variables for the single IPW model of R
  # M                       numeric           number of imputed datasets
model_selection <- function(stacked, first.stage.y, first.stage.var, second.stage.selection, second.stage.y, second.stage.var, ipw.y, ipw.var, M){ 
  
  # Model selection for the first-stage outcome
  dat.long <- dataLong(stacked, timeColumn=first.stage.y[1], censColumn=first.stage.y[2])
  max.formula.1 <- as.formula(paste0("y ~ timeInt + ", paste0(first.stage.var[[1]], collapse=" + ")))
  min.formula.1 <- as.formula(paste0("y ~ timeInt + ", paste0(first.stage.var[[2]], collapse=" + ")))
  max.model.1 <- glm(max.formula.1, family=binomial, data=dat.long, weights=rep(1/M, nrow(dat.long)))
  min.model.1 <- glm(min.formula.1, family=binomial, data=dat.long, weights=rep(1/M, nrow(dat.long)))
  selected.1 <- stepAIC(max.model.1, scope=list("upper"=max.model.1, "lower"=min.model.1), direction="backward", trace=0)
  selected.1.formula <- selected.1$formula
  
  # Model selection for the second-stage outcome
  data.2 <- stacked[stacked[, second.stage.selection]==1, ]
  data.2[, second.stage.y] <- droplevels(data.2[, second.stage.y])
  max.formula.2 <- as.formula(paste0(second.stage.y, " ~ ", paste0(second.stage.var[[1]], collapse=" + ")))
  min.formula.2 <- as.formula(paste0(second.stage.y, " ~ ", paste0(second.stage.var[[2]], collapse = " + ")))
  max.model.2 <- multinom(max.formula.2, data=data.2, trace=FALSE, weights=rep(1/M, nrow(data.2)))
  min.model.2 <- multinom(min.formula.2, data=data.2, trace=FALSE, weights=rep(1/M, nrow(data.2)))
  selected.2 <- stepAIC(max.model.2, scope=list("upper"=max.model.2, "lower"=min.model.2), direction="backward", trace=0)
  selected.2.formula <- as.formula(strsplit(as.character(selected.2$call), "=")[[2]])
  
  # Model selection for the inverse probability weights
  if (length(ipw.y)==1){
    max.formula.3 <- as.formula(paste0(ipw.y, " ~ ", paste0(ipw.var[[1]], collapse=" + ")))
    min.formula.3 <- as.formula(paste0(ipw.y, " ~ ", paste0(ipw.var[[2]], collapse=" + ")))
    max.model.3 <- glm(max.formula.3, family=binomial, data=stacked, weights=rep(1/M, nrow(stacked)))
    min.model.3 <- glm(min.formula.3, family=binomial, data=stacked, weights=rep(1/M, nrow(stacked)))
    selected.3 <- stepAIC(max.model.3, scope=list("upper"=max.model.3, "lower"=min.model.3), direction="backward", trace=0)
    selected.3.formula <- selected.3$formula
  } else if (length(ipw.y)==2){
    selected.3.formula <- list()
    for (i in 1:2){
      max.formula.3 <- as.formula(paste0(ipw.y[i], " ~ ", paste0(ipw.var[[i]][[1]], collapse=" + ")))
      min.formula.3 <- as.formula(paste0(ipw.y[i], " ~ ", paste0(ipw.var[[i]][[2]], collapse=" + ")))
      if (i == 1){
        data <- stacked
      } else if (i == 2){
        data <- stacked[ipw.y[1]==1, ]
      }
      max.model.3 <- glm(max.formula.3, family=binomial, data=data, weights=rep(1/M, nrow(stacked)))
      min.model.3 <- glm(min.formula.3, family=binomial, data=data, weights=rep(1/M, nrow(stacked)))
      selected.3 <- stepAIC(max.model.3, scope=list("upper"=max.model.3, "lower"=min.model.3), direction="backward", trace=0)
      selected.3.formula[[i]] <- selected.3$formula
    }
  }
  
  # Returning the model formulas
  return(list("firstStage"=selected.1.formula,
              "secondStage"=selected.2.formula,
              "IPW"=selected.3.formula))
}


#-------------------------------------------------------------------------------
# Model Fitting Functions
#-------------------------------------------------------------------------------

### Fitting the selected discrete time-to-event model
  # formula             formula           selected first-stage model
  # data                dataframe         stacked imputed datasets
  # time                character         column name of the time variable
  # delta               character         column name of the event indicator
fit_stage_1 <- function(formula, data, time, delta){
  dat.long <- dataLong(data, timeColumn=time, censColumn=delta)
  
  # Fitting a discrete-time logistic regression model
  dat.glm <- glm(formula, family=binomial, data=dat.long)
  
  return(rbind(coef(dat.glm), diag(vcov(dat.glm))))
}

### Fitting the selected pregnancy outcome model
  # formula             formula           selected first-stage model
  # data                dataframe         stacked imputed datasets
  # selection.event     character         column name of the selection event for the second-stage model
fit_stage_2 <- function(formula, data, selection.event){
  out.data <- data[data[, selection.event]==1, ]
  out.data$final.status <- droplevels(out.data$final.status)
  
  # Fitting the multinomial logistic regression model
  dat.mn <- multinom(formula, data=out.data, trace=FALSE, Hess=TRUE)
  
  return(rbind(coef(dat.mn), tryCatch(matrix(diag(solve(dat.mn$Hessian)), nrow=2, byrow=T), error=function(e) matrix(rep(NA, length(coef(dat.mn))), nrow=2))))
}

### Fitting the selected IPW model
  # formula             list              either (i) a list containing formula objects for the selected IPW models for R_1 and R_2 or (ii) a formula object for the selected IPW model for R
  # data                dataframe         stacked imputed datasets       
fit_stage_3 <- function(formula, data){
  
  if (!is.list(formula)){
    # Fitting the IPW validation model
    dat.ipw <- glm(formula, family=binomial, data=data)
    out <- rbind(coef(dat.ipw), diag(vcov(dat.ipw)))
  } else {
    out <- list()
    # Fitting the IPW validation model
    dat.ipw.1 <- glm(formula[[1]], family=binomial, data=data)
    out[[1]] <- rbind(coef(dat.ipw.1), diag(vcov(dat.ipw.1)))
    # Fitting the IPW validation model
    dat.ipw.2 <- glm(formula[[2]], family=binomial, data=data[data[, formula[[1]][[2]]]==1,])
    out[[2]] <- rbind(coef(dat.ipw.2), diag(vcov(dat.ipw.2)))
  }
  
  return(out)
}

#-------------------------------------------------------------
# Risk Prediction/IPW Estimation Functions
#-------------------------------------------------------------

### Obtaining fitted first-stage probabilities
  # data                dataframe         dataset for prediction
  # formula             formula           selected first-stage prediction model
  # coef                numeric           vector of fitted regression coefficients for the first-stage prediction model
pred_stage_1 <- function(data, formula, coef){
  
  # Expanding the dataset to include discrete-time indicators
  tau <- length(grep("timeInt", names(coef))) + 1
  n <- nrow(data)
  data <- data[rep(seq_len(n), each=tau), ]
  time.dat <- data.frame(X = rep(1:tau, n))
  for (i in 2:tau){
    temp <- ifelse(time.dat$X == i, 1, 0)
    time.dat <- cbind(time.dat, temp)
  }
  colnames(time.dat)[2:tau] <- sapply(2:tau, FUN=function(x)paste0("timeInt", x))
  time.dat$X <- NULL
  data <- cbind(data, time.dat)
  
  # Creating the X matrix
  update.fm <- paste0("~ ", paste0(colnames(time.dat), collapse=" + "), " + . - timeInt")
  X <- model.matrix(update(formula[c(1, 3)], as.formula(update.fm)), model.frame(~ ., data, na.action=na.pass))
  
  # Estimated probabilities
  linear <- as.matrix(X) %*% coef
  X <- as.data.frame(X)
  X$ppt.cond.prob <- exp(linear)/(1 + exp(linear))
  
  # Predicted probabilities: joint
  X$ppt.prob <- rep(NA, nrow(X))
  X$ppt.prob[seq(1, nrow(X), by=tau)] <- X$ppt.cond.prob[seq(1, nrow(X), by=tau)]
  for (i in 2:tau){
    temp <- 1
    for (j in 1:(i-1)){
      temp <- temp*(1-X$ppt.cond.prob[seq(j, nrow(X), by=tau)])
    }
    X$ppt.prob[seq(i, nrow(X), by=tau)] <- temp*X$ppt.cond.prob[seq(i, nrow(X), by=tau)]
  }
  
  # Output matrix
  out <- matrix(NA, nrow=n, ncol=tau*2)
  for (i in 1:tau){
    out[, i] <- X$ppt.prob[seq(i, nrow(X), by=tau)]
  }
  out[, tau+1] <- out[,1]
  for (i in 2:tau){
    out[, tau + i] <- apply(out[, 1:i], 1, sum)
  }
  name.1 <- sapply(1:tau, FUN=function(x)paste0("pred.first.ppt.at.", x))
  name.2 <- sapply(1:tau, FUN=function(x)paste0("pred.first.ppt.by.", x))
  colnames(out) <- c(name.1, name.2)
  
  return(out)
}

### Obtaining fitted second-stage probabilities
  # data                dataframe         dataset for prediction
  # formula             formula           selected second-stage prediction model
  # coef                numeric           vector of fitted regression coefficients for the second-stage prediction model
pred_stage_2 <- function(data, formula, coef){
  
  # Creating the X matrix
  X <- model.matrix(formula[c(1, 3)], model.frame(~ ., data, na.action=na.pass))

  # Linear predictor: clinical loss vs. full-term birth
  linear.1 <- as.matrix(X) %*% coef[1,]
  
  # Linear predictor: preterm birth vs. full-term birth
  linear.2 <- as.matrix(X) %*% coef[2,]
  
  # Output matrix
  out <- matrix(NA, nrow=nrow(data), ncol=3)
  out[, 1] <- 1/(1 + exp(linear.1) + exp(linear.2))
  out[, 2] <- exp(linear.1)/(1 + exp(linear.1) + exp(linear.2))
  out[, 3] <- exp(linear.2)/(1 + exp(linear.1) + exp(linear.2))
  colnames(out) <- c("cond.prob.full", "cond.prob.clin", "cond.prob.pre")
  
  return(out)
}

### Obtaining fitted probability of complete data
  # data                dataframe         dataset for prediction
  # formula             list              either (i) a list containing the two selected IPW models for R_1 and R_2 or (ii) a formula object containing the IPW model for R
  # coef                list              either (i) a list containing the corresponding fitted regression coefficients or (ii) a vector of fitted regression coefficients for the IPW model
pred_r <- function(data, formula, coef){
  
  if (!is.list(formula)){
    # Creating the X matrix
    X <- model.matrix(formula[c(1, 3)], model.frame(~ ., data, na.action=na.pass))
    
    # Estimated probabilities
    linear <- as.matrix(X) %*% coef
    out <- exp(linear)/(1 + exp(linear))
  } else{
    out <- list()
    for (i in 1:length(formula)){
      # Creating the X matrix
      X <- model.matrix(formula[[i]][c(1, 3)], model.frame(~ ., data, na.action=na.pass))
      
      # Estimated probabilities
      linear <- as.matrix(X) %*% coef
      out[[i]] <- exp(linear)/(1 + exp(linear))
    }
    out <- out[[2]]*out[[1]]
  }
  
  return(out)
}

### Constructing final preconception risk profile
  # x                   numeric           vector containing all first- and second-stage predictions 
risk_profile <- function(x){
  prob.np <- 1-x[1]
  prob.clin <- x[2]*x[1]
  prob.pre <- x[3]*x[1]
  prob.full <- x[4]*x[1]
  
  out <- c(prob.np, prob.clin, prob.pre, prob.full)
  
  return(out)
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


