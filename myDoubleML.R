#######################################################################################################################################################
#  This is modified from "Double/Debiased Machine Learning of Treatment and Causal Parameters",  AER May 2017     
# Data source: Yannis Bilias, "Sequential Testing of Duration Data: The Case of Pennsylvania 'Reemployment Bonus' Experiment", 
# Journal of Applied Econometrics, Vol. 15, No. 6, 2000, pp. 575-594

# Description of the data set taken from Bilias (2000):

# The 23 variables (columns) of the datafile utilized in the article may be described as follows:

# abdt:       chronological time of enrollment of each claimant in the Pennsylvania reemployment bonus experiment.
# tg:         indicates the treatment group (bonus amount - qualification period) of each claimant. 
# inuidur1:   a measure of length (in weeks) of the first spell ofunemployment
# inuidur2:   a second measure for the length (in weeks) of 
# female:     dummy variable; it indicates if the claimant's sex is female (=1) or male (=0).
# black:      dummy variable; it  indicates a person of black race (=1).
# hispanic:   dummy variable; it  indicates a person of hispanic race (=1).
# othrace:    dummy variable; it  indicates a non-white, non-black, not-hispanic person (=1).
# dep:        the number of dependents of each claimant;
# q1-q6:      six dummy variables indicating the quarter of experiment  during which each claimant enrolled.
# recall:     takes the value of 1 if the claimant answered ``yes'' when was asked if he/she had any expectation to be recalled
# agelt35:    takes the value of 1 if the claimant's age is less  than 35 and 0 otherwise.
# agegt54:    takes the value of 1 if the claimant's age is more than 54 and 0 otherwise.
# durable:    it takes the value of 1 if the occupation  of the claimant was in the sector of durable manufacturing and 0 otherwise.
# nondurable: it takes the value of 1 if the occupation of the claimant was in the sector of nondurable manufacturing and 0 otherwise.
# lusd:       it takes the value of 1 if the claimant filed  in Coatesville, Reading, or Lancaster and 0 otherwise.
#             These three sites were considered to be located in areas characterized by low unemployment rate and short duration of unemployment.
# husd:       it takes the value of 1 if the claimant filed in Lewistown, Pittston, or Scranton and 0 otherwise.
#             These three sites were considered to be located in areas characterized by high unemployment rate and short duration of unemployment.
# muld:       it takes the value of 1 if the claimant filed in Philadelphia-North, Philadelphia-Uptown, McKeesport, Erie, or Butler and 0 otherwise.
#             These three sites were considered to be located in areas characterized by moderate unemployment rate and long duration of unemployment."
#######################################################################################################################################################

###################### Loading packages ###########################
library(doParallel)
library(doSNOW)
library(foreign)
library(gbm)
library(glmnet)
library(hdm)
library(MASS)
library(matrixStats)
library(mnormt)
library(nnet)
library(quadprog)
library(quantreg)
library(randomForest)
library(rpart)
library(sandwich)

################ Loading functions and Data ########################


#rm(list = ls())  # Clear everything out so we're starting clean
source("myML_Funs.R")  
source("myMomentFuns.R")  
options(warn=-1)
#set.seed(1210);
cl <- makeCluster(12, outfile="")
registerDoSNOW(cl)
Penn<- as.data.frame(read.table("penn_jae.dat", header=T ));

########################### Sample Construction ######################

index       <- (Penn$tg==0) | (Penn$tg==4)
data        <- Penn[index,]
data$tg[(data$tg==4)] <- 1
data$dep  <- as.factor(data$dep)
data$inuidur1 <- log(data$inuidur1)

################################ Inputs ##############################

# Outcome Variable
y      <- "inuidur1";

# Treatment Indicator
d      <- "tg";    

# Controls
x      <- "female+black+othrace+ dep+q2+q3+q4+q5+q6+agelt35+agegt54+durable+lusd+husd"         # use this for tree-based methods like forests and boosted trees
xl     <- "(female+black+othrace+dep+q2+q3+q4+q5+q6+agelt35+agegt54+durable+lusd+husd)^2";     # use this for rlasso etc.

# Method names: Boosting, Nnet, RLasso, PostRLasso, Forest, Trees, Ridge, Lasso, Elnet, Ensemble

Boosting     <- list(bag.fraction = .5, train.fraction = 1.0, interaction.depth=2, n.trees=1000, shrinkage=.01, n.cores=1, cv.folds=5, verbose = FALSE, clas_dist= 'adaboost', reg_dist='gaussian')
Forest       <- list(clas_nodesize=1, reg_nodesize=5, ntree=1000, na.action=na.omit, replace=TRUE)
RLasso       <- list(penalty = list(homoscedastic = FALSE, X.dependent.lambda =FALSE, lambda.start = NULL, c = 1.1), intercept = TRUE)
Nnet         <- list(size=2,  maxit=1000, decay=0.02, MaxNWts=10000,  trace=FALSE)
Trees        <- list(reg_method="anova", clas_method="class")

arguments    <- list(Boosting=Boosting, Forest=Forest, RLasso=RLasso, Nnet=Nnet, Trees=Trees)

ensemble     <- list(methods=c("RLasso", "Boosting", "Forest", "Nnet"))              # specify the methods for the ensemble estimation
methods      <- c("RLasso","Trees", "Forest", "Boosting", "Nnet", "Ensemble")        # method names to be estimated
ite          <- 2                                                                    # number of iteration

################################ Estimation ##################################################

r <- foreach(k = 1:ite, .combine='rbind', .inorder=FALSE, .packages=c('MASS','randomForest','neuralnet','gbm', 'sandwich', 'hdm', 'nnet', 'rpart','glmnet')) %dopar% { 
    
    # table 1-2, Panel B, 2 fold
    res <- DoubleML(data, y, d, x, xl, methods=methods, nfold=2, est="plinear", arguments=arguments, ensemble=ensemble, silent=FALSE) 
    
    # table 1-2, Panel B, 5 fold
    #res <- DoubleML(data, y, d, x, xl, methods=methods, nfold=5, est="plinear", arguments=arguments, ensemble=ensemble, silent=FALSE) 
    
    # table 1-2, Panel A, 2 fold
    #res <- DoubleML(data, y, d, x, xl, methods=methods, nfold=2, est="interactive", arguments=arguments, ensemble=ensemble, silent=FALSE, trim=c(0.01,0.99)) 
    
    # table 1-2, Panel A, 5 fold
    #res <- DoubleML(data, y, d, x, xl, methods=methods, nfold=5, est="interactive", arguments=arguments, ensemble=ensemble, silent=FALSE, trim=c(0.01,0.99))  
    
    data.frame(t(res[1,]), t(res[2,]))
}

################################ Compute Output Table ########################################

result           <- matrix(0,4, length(methods)+1)
colnames(result) <- cbind(t(methods), "best")
rownames(result) <- cbind("Mean ATE", "se", "Median ATE", "se")

result[1,]       <- colMeans(r[,1:(length(methods)+1)])
result[2,]       <- (sqrt(colSums(r[,(length(methods)+2):ncol(r)]^2+(r[,1:(length(methods)+1)] - colMeans(r[,1:(length(methods)+1)]))^2)/ite))
result[3,]       <- colQuantiles(r[,1:(length(methods)+1)], probs=0.5)
result[4,]       <- colQuantiles(sqrt(r[,(length(methods)+2):ncol(r)]^2+(r[,1:(length(methods)+1)] - colQuantiles(r[,1:(length(methods)+1)], probs=0.5))^2), probs=0.5)

print(result)






