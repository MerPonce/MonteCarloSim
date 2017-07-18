remove(list=ls())
library(aod)
library(Boruta)
library(caret)
library(causalTree)
library(class)
library(randomForest)
library(speff2trial)
library(varSelRF)

#set.seed(12345678)
obs <- 1000
treatEffect <- 0
covariateEffect<-1

normalize <- function(x) {
    num <- x - min(x)
    denom <- max(x) - min(x)
    return (num/denom)
}

mySimFunc<-function(obs, treatEffect, covariateEffect){
        treatment <- 1*(runif(obs) < 0.5)
        x1 <- rnorm(obs)
        x3 <- rnorm(obs)
        x8 <- rnorm(obs)
        x2 <- 0.2*x1 + 0.98*rnorm(obs)
        x4 <- 1*(runif(obs) < 0.3)
        x5 <- 0.1*x1 + 0.2*x3 +0.97*rnorm(obs)
        x6 <- 1*(runif(obs) < 0.5)
        x7 <- 0.1*x3 + 0.99*rnorm(obs)
        ## x1-x4 are important covariates and x5-x8 are unimportant covariates
        b1 <- 1*covariateEffect
        b2 <- 1*covariateEffect
        b3 <- 1*covariateEffect
        b4 <- 1*covariateEffect
        b5 <- 0.01*covariateEffect
        b6 <- 0.01*covariateEffect
        b7 <- 0.01*covariateEffect
        b8 <- 0.01*covariateEffect
        xb <- b1*x1 + b2*x2 + b3*x3 + b4*x4 +b5*x5 +b6*x6 + b7*x7 + b8*x8 + treatEffect*treatment
        
        ## DGP
        prob <- exp(xb)/(1+exp(xb))
        y <- rbinom(obs,1,prob) 
        myData <- data.frame(y,x1,x2,x3,x4,x5,x6,x7,x8,treatment)
        
        ## OLS model
        controlMod <- glm(y~x1+x2+x3+x4+x5+x6+x7+x8, data=myData[treatment==0,], family="binomial")
        treatMod <- glm(y~x1+x2+x3+x4+x5+x6+x7+x8, data=myData[treatment==1,], family="binomial")
        
        covariate <- cbind(x1,x2,x3,x4,x5,x6,x7,x8)
        colnames(covariate) <- c("x1","x2","x3","x4","x5","x6","x7","x8")
        expit <- function(linpred){
            exp(linpred)/(1+exp(linpred))
        }
        controlPred <- as.numeric(expit(cbind(rep(1,obs),covariate) %*% controlMod$coef) >0.5)
        treatPred <- as.numeric(expit(cbind(rep(1,obs),covariate) %*% treatMod$coef) > 0.5)
        
        ## Show the prediction result
        #table(controlPred>0.5, y)
        #table(treatPred>0.5, y)
        accControlMod <- sum((controlPred>0.5)==y)/length(y)
        accTreatMod <- sum((treatPred>0.5)==y)/length(y)
        
        semiPara <- speff(y ~ 1, endpoint="dichotomous", data=myData, trt.id="treatment",
                          endCtrlPre=controlPred, endTreatPre=treatPred)
        
        ## Random forest model
        controlModRF <- randomForest(as.factor(y)~.-treatment, data=myData[treatment==0,], importance=TRUE, ntree=50)
        treatModRF <- randomForest(as.factor(y)~.-treatment, data=myData[treatment==1,], importance=TRUE, ntree=50)
        
        
        #varImpPlot(controlModRF)
        #varImpPlot(treatModRF)
        
        ## Select important covariates from random forest model and then create the prediction
        #borutaControl <- Boruta(as.factor(y)~.-treatment, data = myData[treatment==0,], doTrace = 0)
        #borutaTreat <- Boruta(as.factor(y)~.-treatment, data = myData[treatment==1,], doTrace = 0)
        
        #varCtrlRF<-getSelectedAttributes(borutaControl)
        #varTreatRF<-getSelectedAttributes(borutaTreat)
        #ctrlFormula<-paste("y~", varCtrlRF[1], sep="")
        #i<-2
        #while(i <= length(varCtrlRF)) {
        #    ctrlFormula<-paste(ctrlFormula, varCtrlRF[i], sep="+")
        #    i=i+1
        #}
        #treatFormula<-paste("y~",varTreatRF[1], sep="")
        #i<-2
        #while (i <= length(varTreatRF)) {
        #    treatFormula<-paste(treatFormula, varTreatRF[i], sep="+")
        #    i = i+1
        #}
        
        #tempData<-myData
        #tempData[,1]<-data.frame(apply(tempData[1], 2, as.factor))
        
        #if (length(varCtrlRF)==0) {
        #    ctrlFormula<-paste("y~x1+x2+x3+x4+x5+x6+x7+x8")
        #    varCtrlRF<-c("x1", "x2","x3","x4","x5","x6","x7","x8")
        #}
        #if (length(varTreatRF)==0) {
        #    treatFormula<-paste("y~x1+x2+x3+x4+x5+x6+x7+x8")
        #    varTreatRF<-c("x1", "x2","x3","x4","x5","x6","x7","x8")
        #}
        #print(ctrlFormula)
        #print(treatFormula)
        #controlModRF <- glm(as.formula(ctrlFormula), data=tempData[treatment==0,],family="binomial")
        #treatModRF <- glm(as.formula(treatFormula), data=tempData[treatment==1,], family="binomial")
        
        #controlPredRF <- as.numeric(expit(cbind(rep(1,obs),covariate[,varCtrlRF]) %*% controlModRF$coef) >0.5)
        #treatPredRF <- as.numeric(expit(cbind(rep(1,obs),covariate[,varTreatRF]) %*% treatModRF$coef) >0.5)
        #semiParaRF <- speff(y ~ 1, endpoint="dichotomous", data=myData, trt.id="treatment", 
        #                    endCtrlPre=controlPredRF, endTreatPre=treatPredRF)
        
        controlPredRF <- predict(controlModRF, myData)
        treatPredRF <- predict(treatModRF, myData)
        
        #table(controlPredRF, y)
        accControlModRF <-sum(as.numeric(as.character(controlPredRF))==y)/length(y)
        #table(treatPredRF, y)
        accTreatModRF<-sum(as.numeric(as.character(treatPredRF))==y)/length(y)
        
        
        semiParaRF <- speff(y ~ 1, endpoint="dichotomous", data=myData, trt.id="treatment", 
                          endCtrlPre=as.numeric(levels(controlPredRF)[controlPredRF]), 
                          endTreatPre=as.numeric(levels(treatPredRF)[treatPredRF]))
        # KNN
        knnDF<-as.data.frame(apply(as.matrix(covariate),2,normalize))
        
        controlPredKnn<-knn(train=knnDF[treatment==0,], test=knnDF, 
                           as.factor(y[treatment==0]), k = 3, prob=TRUE)
        treatPredKnn<-knn(train=knnDF[treatment==1,], test=knnDF, 
                         as.factor(y[treatment==1]), k = 3, prob=TRUE)
        semiParaKnn <- speff(y ~ 1, endpoint="dichotomous", data=myData, trt.id="treatment", 
                            endCtrlPre=as.numeric(levels(controlPredKnn)[controlPredKnn]), 
                            endTreatPre=as.numeric(levels(treatPredKnn)[treatPredKnn]))
        accControlModKnn <- sum(as.numeric(as.character(controlPredKnn))==y)/length(y)
        accTreatModKnn <- sum(as.numeric(as.character(treatPredKnn))==y)/length(y)
                     
        #summary(semiPara)
        #summary(semiParaRF)
        te<-semiPara$coef[,3]
        teRF<-semiParaRF$coef[,3]
        teKnn<-semiParaKnn$coef[,3]
        se <-sqrt(semiPara$varbeta)
        seRF <-sqrt(semiParaRF$varbeta)
        seKnn <-sqrt(semiParaKnn$varbeta)
        ## return the p value of speff
        pval <- 2*pnorm(-abs(semiPara$coef[,3])/sqrt(semiPara$varbeta))
        pvalRF <- 2*pnorm(-abs(semiParaRF$coef[,3])/sqrt(semiParaRF$varbeta))
        pvalKnn <- 2*pnorm(-abs(semiParaKnn$coef[,3])/sqrt(semiParaKnn$varbeta))
        ## return the p value of recovering true treatment effect
        #pval <- 2*pnorm(-abs(semiPara$coef[,3]-c(treatEffect,treatEffect))/sqrt(semiPara$varbeta))
        #pvalRF <- 2*pnorm(-abs(semiParaRF$coef[,3]-c(treatEffect,treatEffect))/sqrt(semiParaRF$varbeta))
        result <-list("naiveTE"=te[1],"olsTE" =te[2], "rfTE"=teRF[2], 
                      "naiveSE"=se[1], "olsSE"=se[2], "rfSE"=seRF[2], 
                      "pvalNaive"=pval[1], "pvalOLS"=pval[2], "pvalRF"=pvalRF[2],
                      "accControl"=accControlMod, "accTreat"=accTreatMod, 
                      "accControlRF"=accControlModRF, "accTreatRF"=accTreatModRF,
                      "knnTE"=teKnn[2], "knnSE"=seKnn[2], "pvalKnn"=pvalKnn[2], 
                      "accControlKnn"= accControlModKnn, "accTreatKnn"=accTreatModKnn)
        return(result)
}

simNum <- 1000
mySim <- function(simNum, obs, treatEffect, covariateEffect, alpha = 0.05){
    result = replicate(simNum, mySimFunc(obs, treatEffect, covariateEffect))
    # Compute average slope
    size<-length(result)/simNum
    
    temp<-result
    temp[7,]<-temp[7,] < alpha
    temp[8,]<-temp[8,] < alpha
    temp[9,]<-temp[9,] < alpha
    temp[16,]<-temp[16,] < alpha
    
    listMeans<-function(data, size){
        means<-c()
        for (i in 1:size) {
            means<-c(means, mean(unlist(data[i,]), na.rm=TRUE))
        }
        return(means)
    }
    means<-listMeans(temp, size)
    
    listSd<-function(data, size) {
        SD<-c()
        for (i in 1:size) {
            SD<-c(SD, sd(unlist(data[i,]), na.rm=TRUE)/sqrt(sum(ifelse(is.na(result[i,]),0,1))))
        }
        return(SD)
    }
    SDs <- listSd(temp, size)
    
    simResult <- list("naiveTE"=means[1], "olsTE"=means[2], "rfTE"=means[3], 
                      "naiveSE"=means[4], "olsSE"=means[5], "rfSE"=means[6], 
                      "naivePower"=means[7], "olsPower"=means[8], "rfPower"=means[9],
                      "accControl"=means[10], "accTreat"=means[11], "accControlRF"=means[12], "accTreatRF"=means[13],
                      "knnTE"=means[14], "knnSE"=means[15], "knnPower"=means[16], "accControlKnn"=means[17],"accTreatKnn"=means[18],
                      "se_of_naiveTE"=SDs[1], "se_of_olsTE"=SDs[2], "se_of_rfTE"=SDs[3], 
                      "se_of_naiveSE"=SDs[4], "se_of_olsSE"=SDs[5], "se_of_rfSE"=SDs[6],
                      "se_of_naivePower"=SDs[7], "se_of_olsPower"=SDs[8], "se_of_rfPower"=SDs[9], 
                      "se_of_accControl"=SDs[10], "se_of_accTreat"=SDs[11], "se_of_accControlRF"=SDs[12], "se_of_accTreatRF"=SDs[13])
    
    return(simResult)
}

library(parallel)

alpha<-0.05
grid <- expand.grid(covariateEffect = seq(0,2,0.5), treatEffect = seq(0, 1, 0.2))
start.time <- Sys.time()
simulatedResults <- mcmapply(mySim, covariateEffect = grid$covariateEffect, treatEffect = grid$treatEffect,
                             # set fixed parameters
                             MoreArgs = list(simNum,obs, alpha),
                             # distribute work over CPU cores
                             mc.cores = 1)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
# Output results (transposed for clarity)
simulatedResults <- cbind(grid, data.frame(t(simulatedResults)))
#write.csv(as.matrix(simulatedResults[1:28]), paste0("simulatedResults", "_simNUM", simNum, "_obs", obs,".csv"), row.names = FALSE)
#test<-matrix(simulatedResults[29], nrow=simNum, ncol=13)
#write.csv(matrix(simulatedResults[29], nrow=obs, ncol=13), paste0("simulatedResultsDetails", "_simNUM", simNum, "_obs", obs,".csv"), row.names = FALSE)
#simulatedResults<-read.csv("simulatedResults_simNum10000_obs1000.csv")

library(ggplot2)
df<-as.data.frame(simulatedResults)
powerVersusTreatment <- ggplot() + geom_line(data=df, aes(x=treatEffect, y=unlist(naivePower), color ="Naive")) +
    #geom_errorbar(data=df,aes(x=treatEffect, ymin=unlist(naivePower)-unlist(se_of_naivePower), ymax=unlist(naivePower)+unlist(se_of_naivePower)), width=.1, color ="red") +
    geom_line(data=df, aes(x=treatEffect, y=unlist(olsPower), color ="OLS")) +
    #geom_errorbar(data=df,aes(x=treatEffect, ymin=unlist(olsPower)-unlist(se_of_olsPower), ymax=unlist(olsPower)+unlist(se_of_olsPower)), width=.1, color ="green") +
    geom_line(data=df, aes(x=treatEffect, y=unlist(rfPower), color ="RF")) + 
    #geom_errorbar(data=df,aes(x=treatEffect, ymin=unlist(rfPower)-unlist(se_of_rfPower), ymax=unlist(rfPower)+unlist(se_of_rfPower)), width=.1, color ="blue") +
    geom_line(data=df, aes(x=treatEffect, y=unlist(knnPower), color ="KNN")) + 
    facet_grid(covariateEffect~.) + xlab('Treatment Effect') + ylab('Power')
seVersusTreatment <- ggplot() + geom_line(data=df, aes(x=treatEffect, y=unlist(naiveSE), color ="Naive")) +
    #geom_errorbar(data=df,aes(x=treatEffect, ymin=unlist(naiveSE)-unlist(se_of_naiveSE), ymax=unlist(naiveSE)+unlist(se_of_naiveSE)), width=.1, color ="red") +
    geom_line(data=df, aes(x=treatEffect, y=unlist(olsSE), color ="OLS")) +
    #geom_errorbar(data=df,aes(x=treatEffect, ymin=unlist(olsSE)-unlist(se_of_olsSE), ymax=unlist(olsSE)+unlist(se_of_olsSE)), width=.1, color ="green") +
    geom_line(data=df, aes(x=treatEffect, y=unlist(rfSE), color ="RF")) + 
    #geom_errorbar(data=df,aes(x=treatEffect, ymin=unlist(rfSE)-unlist(se_of_rfSE), ymax=unlist(rfSE)+unlist(se_of_rfSE)), width=.1, color ="blue") +
    geom_line(data=df, aes(x=treatEffect, y=unlist(knnSE), color ="KNN")) + 
    facet_grid(covariateEffect~.) + xlab('Treatment Effect') + ylab('Standard error')
teVersusTreatment <- ggplot() + geom_line(data=df, aes(x=treatEffect, y=unlist(naiveTE), color ="Naive")) +
    #geom_errorbar(data=df,aes(x=treatEffect, ymin=unlist(naiveTE)-unlist(se_of_naiveTE), ymax=unlist(naiveTE)+unlist(se_of_naiveTE)), width=.1, color ="red") +
    geom_line(data=df, aes(x=treatEffect, y=unlist(olsTE), color ="OLS")) +
    #geom_errorbar(data=df,aes(x=treatEffect, ymin=unlist(olsTE)-unlist(se_of_olsTE), ymax=unlist(olsTE)+unlist(se_of_olsTE)), width=.1, color ="green") +
    geom_line(data=df, aes(x=treatEffect, y=unlist(rfTE), color ="RF")) + 
    #geom_errorbar(data=df,aes(x=treatEffect, ymin=unlist(rfTE)-unlist(se_of_rfTE), ymax=unlist(rfTE)+unlist(se_of_rfTE)), width=.1, color ="blue") +
    geom_line(data=df, aes(x=treatEffect, y=unlist(knnTE), color ="KNN")) + 
    facet_grid(covariateEffect~.) + xlab('Treatment Effect') + ylab('Estimated Treatment Effect')
accVersusTreatment <- ggplot() + geom_line(data=df, aes(x=treatEffect, y=unlist(accControl), color ="OLS")) +
    geom_line(data=df, aes(x=treatEffect, y=unlist(accControlRF), color ="RF")) +
    geom_line(data=df, aes(x=treatEffect, y=unlist(accControlKnn), color ="KNN")) + 
    facet_grid(covariateEffect~.) + xlab('Treatment Effect') + ylab('Accuracy')

powerVersusTreatment
seVersusTreatment
teVersusTreatment
accVersusTreatment

#ggsave(filename=paste0("powerVersusTreatment_simNum", simNum, "_obs", obs, ".pdf"), plot=powerVersusTreatment)
#ggsave(filename=paste0("teVersusTreatment_simNum", simNum, "_obs", obs, ".pdf"), plot=teVersusTreatment)
#ggsave(filename=paste0("seVersusTreatment_simNum", simNum, "_obs", obs, ".pdf"), plot=seVersusTreatment)
#ggsave(filename=paste0("accVersusTreatment_simNum", simNum, "_obs", obs, ".pdf"), plot=accVersusTreatment)
 
