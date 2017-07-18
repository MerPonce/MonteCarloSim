DoubleML <- function(data, y, d, xx, xL, methods, nfold, est, arguments, ensemble, silent=FALSE, trim) {
    K         <- nfold
    TE        <- matrix(0,1,length(methods)+1)
    STE       <- matrix(0,1,length(methods)+1)
    result    <- matrix(0,2,length(methods)+1)
    MSE1      <- matrix(0,length(methods), K)
    MSE2      <- matrix(0,length(methods), K)
    MSE3      <- matrix(0,length(methods), K)
    cond.comp <- matrix(list(),length(methods),K)
    
    binary    <- as.numeric(checkBinary(data[,d]))
    split     <-runif(nrow(data))
    cvgroup   <- as.numeric(cut(split, quantile(split, probs = seq(0,1,1/K)),include.lowest = TRUE))
    
    for (k in 1:length(methods)) {
        if (silent==FALSE) {
            cat(methods[k], '\n')
        }
        
        if (any(c("RLasso", "PostRLasso", "Ridge", "Lasso", "Elnet") ==methods[k])){
            x=xL
        } else {
            x=xx
        }
        
        for(j in 1:K) {
            if (silent==FALSE) {
                cat(' fold',j,'\n')
            }
            
            ii <-cvgroup == j
            nii <- cvgroup != j
            
            datause <-as.data.frame(data[nii,])
            dataout <-as.data.frame(data[ii,])
            
            if(est=="interactive" && (length(methods)>0)) {
                if (methods[k]=="Ensemble") {
                    cond.comp[[k,j]] <- ensembleF(datause=datause, dataout=dataout,y,d,x,methods[k],plinear=0,xL,binary,arguments,ensemble)
                } else {
                    cond.comp[[k,j]] <- cond_comp(datause=datause, dataout=dataout,y,d,x,methods[k], plinear=0,xL,binary,arguments)
                }
                MSE1[k,j] <- cond.comp[[k,j]]$err.yz0
                MSE2[k,j] <- cond.comp[[k,j]]$err.zy1
                MSE3[k,j] <- cond.comp[[k,j]]$err.z
                
                drop <- which(cond.comp[[k,j]]$mz_x>trim[1] & cond.comp[[k,j]]$mz_x<trim[2])
                mz_x <- cond.comp[[k,j]]$mz_x[drop]
                my_z1x <- cond.comp[[k,j]]$my_z1x[drop]
                my_z0x <- cond.comp[[k,j]]$my_z0x[drop]
                
                yout <- dataout[drop,y]
                dout <- dataout[drop,d]
                
                TE[1,k] <- ATE(yout,dout,my_z1x,my_z0x,mz_x)/k + TE[1,k]
                STE[1,k] <- (1/(K^2))*((SE.ATE(yout, dout, my_z1x, my_z0x, mz_x))^2) + STE[1,k];
            }
            
            if(est=="plinear" && (length(methods)>0)) {
                if(methods[k]=="Ensemble") { 
                    cond.comp[[k,j]] <- ensembleF(datause=datause, dataout=dataout, y, d, x, methods[k], plinear=1, xL, binary, arguments, ensemble)
                } else {
                    cond.comp[[k,j]] <- cond_comp(datause=datause, dataout=dataout, y, d, x, methods[k], plinear=1, xL, binary, arguments)
                }   
                MSE1[k,j] <- cond.comp[[k,j]]$err.y
                MSE2[k,j] <- cond.comp[[k,j]]$err.z
                lm.fit.ry <- lm(as.matrix(cond.comp[[k,j]]$ry) ~ as.matrix(cond.comp[[k,j]]$rz)-1)
                ate <- lm.fit.ry$coef
                HCV.coefs <- vcovHC(lm.fit.ry, type = 'HC')
            }
        }
    }
}

DoubleML <- function(data, y, d, xx, xL, methods, nfold, est, arguments, ensemble, silent=FALSE, trim){
                
                
                
                
                
                STE[1,k]               <- (1/(K^2))*(diag(HCV.coefs)) +  STE[1,k] 
                TE[1,k]                <- ate/K + TE[1,k] ;
                
    
    
    if(est=="interactive"){
        
        min1 <- which.min(rowMeans(MSE1))
        min2 <- which.min(rowMeans(MSE2))
        min3 <- which.min(rowMeans(MSE3))
        
        if(silent==FALSE){
            cat('  best methods for E[Y|X, D=0]:',methods[min1],'\n')
            cat('  best methods for E[Y|X, D=1]:',methods[min2],'\n')
            cat('  best methods for E[D|X]:',methods[min3],'\n')
        }
        
    }
    
    if(est=="plinear"){
        
        min1 <- which.min(rowMeans(MSE1))
        min2 <- which.min(rowMeans(MSE2))
        
        if(silent==FALSE){   
            cat('  best methods for E[Y|X]:',methods[min1],'\n')
            cat('  best methods for E[D|X]:',methods[min2],'\n')
        }    
    }
    
    for(j in 1:K){  
        
        ii = cvgroup == j
        nii = cvgroup != j
        
        datause = as.data.frame(data[nii,])
        dataout = as.data.frame(data[ii,])  
        
        if(est=="interactive"){
            
            drop                   <- which(cond.comp[[min3,j]]$mz_x>trim[1] & cond.comp[[min3,j]]$mz_x<trim[2])      
            mz_x                   <- cond.comp[[min1,j]]$mz_x[drop]
            my_z1x                 <- cond.comp[[min2,j]]$my_z1x[drop]
            my_z0x                 <- cond.comp[[min3,j]]$my_z0x[drop]
            yout                   <- dataout[drop,y]
            dout                   <- dataout[drop,d]
            
            TE[1,(k+1)]            <- ATE(yout, dout, my_z1x, my_z0x, mz_x)/K + TE[1,(k+1)];
            STE[1,(k+1)]           <- (1/(K^2))*((SE.ATE(yout, dout, my_z1x, my_z0x, mz_x))^2) + STE[1,(k+1)];
            
        }
        
        if(est=="plinear"){
            
            lm.fit.ry              <- lm(as.matrix(cond.comp[[min1,j]]$ry) ~ as.matrix(cond.comp[[min2,j]]$rz)-1);
            ate                    <- lm.fit.ry$coef;
            HCV.coefs              <- vcovHC(lm.fit.ry, type = 'HC');
            STE[1,(k+1)]           <- (1/(K^2))*(diag(HCV.coefs)) +  STE[1,(k+1)] 
            TE[1,(k+1)]            <- ate/K + TE[1,(k+1)] ;
            
        }
    }
    
    colnames(result)   <- c(methods, "best") 
    rownames(MSE1)     <- c(methods) 
    rownames(MSE2)     <- c(methods) 
    rownames(MSE3)     <- c(methods) 
    rownames(result)   <- c("ATE", "se")
    result[1,]         <- colMeans(TE)
    result[2,]         <- sqrt(STE)
    
    if(est=="plinear"){   
        table <- rbind(result, cbind(rbind(rowMeans(MSE1), rowMeans(MSE2)), c(NA,NA)))    
        rownames(table)[3:4]   <- c("MSE[Y|X]", "MSE[D|X]") 
    }
    
    if(est=="interactive"){    
        table <- rbind(result, cbind(rbind(rowMeans(MSE1), rowMeans(MSE2), rowMeans(MSE3)), c(NA,NA, NA)))    
        rownames(table)[3:5]   <- c("MSE[Y|X, D=0]", "MSE[Y|X, D=1]", "MSE[D|X]")
    }  
    
    colnames(table)[length(methods)+1] = "best"
    return(table)
}  
