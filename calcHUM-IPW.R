#-------------------------------------------------------------#
#     Implementing calcHUM in C with IPW-Based Adjustment     #
#                  for Verification Bias                      #
#-------------------------------------------------------------#

library(gtools)

### Note:
  # The C implementation for calculating the unadjusted HUM was written by
  # Kyu Ha Lee (Harvard T.H. Chan School of Public Health) and 
  # Catherine Lee (Kaiser Permanente Division of Research)
###


#-------------------------------------------------------------
# Loading the C functions/interface
#-------------------------------------------------------------

dyn.load("loop.so")


#-------------------------------------------------------------
# Wrapper Function
#-------------------------------------------------------------

calcHUM <- function(riskMat, probVec=NULL){
  
  seq <- permutations(4, 4, 1:4)
  n <- nrow(riskMat)
  
  if (is.null(probVec)){
    probVec <- rep(1, n)
  }
  
  # Creating a list that holds index of subjects in categories 1, 2, 3, 4
  indCat <- list()
  for (c in 1:4){
    indCat[[c]] <- which(riskMat[,5] == c)
  }
  Ns <- unlist(lapply(indCat, length))
  
  # Computing the HUM using C implementation
  CR <- 0
  denom <- 0
  val <- .C("loop",
            Data = as.double(as.matrix(riskMat)),
            Prob = as.double(probVec),
            Seq = as.double(as.matrix(seq)),
            Cat1 = as.double(indCat[[1]]),
            Cat2 = as.double(indCat[[2]]),
            Cat3 = as.double(indCat[[3]]),
            Cat4 = as.double(indCat[[4]]),
            n = as.integer(nrow(riskMat)),
            p = as.integer(ncol(riskMat)),
            I = as.integer(Ns[1]),
            J = as.integer(Ns[2]),
            K = as.integer(Ns[3]),
            L = as.integer(Ns[4]),  
            CRval = as.double(CR),
            Denom = as.double(denom))  	  	
  
  CR <- val$CRval; denom <- val$Denom
  print(CR)
  print(denom)
  value <- CR/denom
  return(value)
}
