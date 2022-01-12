#------------------------------------------------------------------------------#
#           Functions for Model Selection, Risk Prediction and                 #
#                   Validation for Pregnancy Prediction                        #
#------------------------------------------------------------------------------#

library(discSurv)
library(nnet)
library(glmnet)
library(MASS)


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
