#########################################################################################################
#  Program:     Functions for estimating moments for using Machine Learning Methods.                    #
#  References: "Double/Debiased Machine Learning of Treatment and Causal Parameters",  AER P&P 2017     #
#              "Double Machine Learning for Treatment and Causal Parameters",  Arxiv 2016               #
#  by V.Chernozhukov, D. Chetverikov, M. Demirer, E. Duflo, C. Hansen, W. Newey                         #           
#########################################################################################################

source("ML_Functions.R") 

DoubleML <- function(data, y, d, xx, xL, methods, nfold, est, arguments, ensemble, silent=FALSE, trim){
  
  K         <- nfold
  
  TE        <- matrix(0,1,(length(methods)+1))
  STE       <- matrix(0,1,(length(methods)+1))
  result    <- matrix(0,2,(length(methods)+1))
  MSE1      <- matrix(0,length(methods),K)
  MSE2      <- matrix(0,length(methods),K)
  MSE3      <- matrix(0,length(methods),K)
  cond.comp <- matrix(list(),length(methods),K)
  
  binary    <- as.numeric(checkBinary(data[,d]))
  
  split     <- runif(nrow(data))
  cvgroup   <- as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/K)),include.lowest = TRUE))  
  
  for(k in 1:length(methods)){   
    
    if(silent==FALSE){
      cat(methods[k],'\n')
    }
    
    if (any(c("RLasso", "PostRLasso", "Ridge", "Lasso", "Elnet")==methods[k])){
      x=xL
    } else {
      x=xx
    }
    
    for(j in 1:K){   
      
      if(silent==FALSE){
        cat('  fold',j,'\n')
      }
      
      ii  <- cvgroup == j
      nii <- cvgroup != j
      
      datause <- as.data.frame(data[nii,])
      dataout <- as.data.frame(data[ii,])  
      
      if(est=="interactive" && (length(methods)>0)){
        
        if(methods[k]=="Ensemble") { cond.comp[[k,j]] <- ensembleF(datause=datause, dataout=dataout, y, d, x, methods[k], plinear=0, xL, binary, arguments, ensemble)}
        else{                       cond.comp[[k,j]] <- cond_comp(datause=datause, dataout=dataout, y, d, x, methods[k], plinear=0, xL, binary, arguments);  }
        
        MSE1[k,j]               <- cond.comp[[k,j]]$err.yz0
        MSE2[k,j]               <- cond.comp[[k,j]]$err.yz1
        MSE3[k,j]               <- cond.comp[[k,j]]$err.z
        
        drop                   <- which(cond.comp[[k,j]]$mz_x>trim[1] & cond.comp[[k,j]]$mz_x<trim[2])      
        mz_x                   <- cond.comp[[k,j]]$mz_x[drop]
        my_z1x                 <- cond.comp[[k,j]]$my_z1x[drop]
        my_z0x                 <- cond.comp[[k,j]]$my_z0x[drop]
        yout                   <- dataout[drop,y]
        dout                   <- dataout[drop,d]
        
        TE[1,k]                <- ATE(yout, dout, my_z1x, my_z0x, mz_x)/K + TE[1,k];
        STE[1,k]               <- (1/(K^2))*((SE.ATE(yout, dout, my_z1x, my_z0x, mz_x))^2) + STE[1,k];
      }
      
      if(est=="plinear" && (length(methods)>0)){
        
        if(methods[k]=="Ensemble") { cond.comp[[k,j]] <- ensembleF(datause=datause, dataout=dataout, y, d, x, methods[k], plinear=1, xL, binary, arguments, ensemble)}
        else{                       cond.comp[[k,j]] <- cond_comp(datause=datause, dataout=dataout, y, d, x, methods[k], plinear=1, xL, binary, arguments)}
        
        MSE1[k,j]              <- cond.comp[[k,j]]$err.y
        MSE2[k,j]              <- cond.comp[[k,j]]$err.z
        
        lm.fit.ry              <- lm(as.matrix(cond.comp[[k,j]]$ry) ~ as.matrix(cond.comp[[k,j]]$rz)-1);
        ate                    <- lm.fit.ry$coef;
        HCV.coefs              <- vcovHC(lm.fit.ry, type = 'HC');
        STE[1,k]               <- (1/(K^2))*(diag(HCV.coefs)) +  STE[1,k] 
        TE[1,k]                <- ate/K + TE[1,k] ;
        
      }
    }  
  }
  
  
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

cond_comp <- function(datause, dataout, y, d, x, method, plinear,xL, binary, arguments){
  
  form_y   <- y
  form_d   <- d
  form_x   <- x
  form_xL  <- xL
  ind_u    <- which(datause[,d]==1)
  ind_o    <- which(dataout[,d]==1)
  err.yz1  <- NULL
  err.yz0  <- NULL
  my_z1x   <- NULL
  my_z0x   <- NULL
  fit.yz1  <- NULL
  fit.yz0  <- NULL
  
  ########################## Boosted  Trees ###################################################;
  
  if(method=="Boosting")
  {
    
    option <- arguments[[method]]
    arg    <- option
    arg[which(names(arg) %in% c("clas_dist","reg_dist"))] <-  NULL
    
    if(plinear==0){
      
      fit.yz1        <- boost(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, distribution=option[['reg_dist']], option=arg)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      my_z1x         <- predict(fit.yz1$model, n.trees=fit.yz1$best, dataout, type="response") 
      
      fit.yz0        <- boost(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_y, distribution=option[['reg_dist']], option=arg)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      my_z0x         <- predict(fit.yz0$model, n.trees=fit.yz0$best, dataout, type="response") 
    }
    
    if(binary==1){
      fit.z          <- boost(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_d, distribution=option[['clas_dist']], option=arg)
      mis.z          <- error(fit.z$yhatout, dataout[,d])$mis
    }
    
    if(binary==0){
      fit.z          <- boost(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_d, distribution=option[['reg_dist']], option=arg)
      mis.z          <- NA
    }
    
    fit.y            <- boost(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_y, distribution=option[['reg_dist']], option=arg)
    
  }  
  
  
  ########################## Neural Network(Nnet Package) ###################################################;   
  
  
  if(method=="Nnet"){
    
    option <- arguments[[method]]
    arg    <- option
    
    if(plinear==0){
      
      fit.yz1        <- nnetF(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, arg=arg)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      dataouts       <- dataout
      dataouts[,!fit.yz1$f] <- as.data.frame(scale(dataouts[,!fit.yz1$f], center = fit.yz1$min, scale = fit.yz1$max - fit.yz1$min))
      my_z1x         <- predict(fit.yz1$model, dataouts)*(fit.yz1$max[fit.yz1$k]-fit.yz1$min[fit.yz1$k])+fit.yz1$min[fit.yz1$k] 
      
      fit.yz0        <- nnetF(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_y, arg=arg)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      dataouts       <- dataout
      dataouts[,!fit.yz0$f] <- as.data.frame(scale(dataouts[,!fit.yz0$f], center = fit.yz0$min, scale = fit.yz0$max - fit.yz0$min))
      my_z0x         <- predict(fit.yz0$model, dataouts)*(fit.yz0$max[fit.yz0$k]-fit.yz0$min[fit.yz0$k])+fit.yz0$min[fit.yz0$k] 
    }
    
    if(binary==1){
      fit.z          <- nnetF(datause=datause, dataout=dataout, form_x=form_x, form_y=form_d, clas=TRUE, arg=arg)
      mis.z          <- error(fit.z$yhatout, dataout[,d])$mis
    }
    
    if(binary==0){
      fit.z          <- nnetF(datause=datause, dataout=dataout, form_x=form_x, form_y=form_d, clas=FALSE, arg=arg)
      mis.z          <- NA
    }
    
    fit.y          <- nnetF(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_y, arg=arg)
    
  } 
  
  ########################## Lasso and Post Lasso(Hdm Package) ###################################################;    
  
  if(method=="RLasso" || method=="PostRLasso"){
    
    post = FALSE
    if(method=="PostRLasso"){ post=TRUE }
    
    option    <- arguments[[method]]
    arg       <- option
    
    if(plinear==0){
      
      fit.yz1        <- rlassoF(datause=datause[ind_u,], dataout=dataout[ind_o,],  form_x, form_y, post, arg=arg)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      my_z1x         <- predict(fit.yz1$model, newdata=formC(form_y, form_x, dataout)$x , type="response") 
      
      fit.yz0        <- rlassoF(datause=datause[-ind_u,], dataout=dataout[-ind_o,], form_x, form_y, post, arg=arg)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      my_z0x         <- predict(fit.yz0$model, newdata=formC(form_y, form_x, dataout)$x, type="response")   
      
    }
    
    if(binary==1){
      fit.z          <- rlassoF(datause=datause, dataout=dataout,  form_x, form_d, post, logit=TRUE, arg=arg)
      mis.z          <- error(fit.z$yhatout, dataout[,d])$mis
    }
    
    if(binary==0){
      fit.z          <- rlassoF(datause=datause, dataout=dataout,  form_x, form_d, post, logit=FALSE, arg=arg)
      mis.z          <- NA
    }   
    
    fit.y          <- rlassoF(datause=datause, dataout=dataout,  form_x, form_y, post, arg=arg)
  }    
  
  
  ########################## Lasso and Post Lasso(Glmnet) Package) ###################################################;    
  
  if(method=="Ridge" || method=="Lasso" || method=="Elnet"){
    
    if(method=="Ridge"){ alp=0 }
    if(method=="Lasso"){ alp=1 }
    if(method=="Elnet"){ alp=0.5 }
    
    option    <- arguments[[method]]
    arg       <- option
    arg[which(names(arg) %in% c("s"))] <-  NULL
    s         <- option[['s']]
    
    if(plinear==0){
      
      fit.yz1        <- lassoF(datause=datause[ind_u,], dataout=dataout[ind_o,],  form_x, form_y, alp=alp, arg=arg, s=s)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      fit.p          <- lm(as.formula(paste(form_y, "~", form_x)),  x = TRUE, y = TRUE, data=dataout);
      my_z1x         <- predict(fit.yz1$model, newx=fit.p$x[,-1] ) 
      
      fit.yz0        <- lassoF(datause=datause[-ind_u,], dataout=dataout[-ind_o,], form_x, form_y, alp=alp, arg=arg, s=s)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      fit.p          <- lm(as.formula(paste(form_y, "~", form_x)),  x = TRUE, y = TRUE, data=dataout);
      my_z0x         <- predict(fit.yz0$model,  newx=fit.p$x[,-1])   
      
    }
    
    if(binary==1){
      fit.z          <- lassoF(datause=datause, dataout=dataout,  form_x, form_d, logit=TRUE, alp=alp, arg=arg, s=s)
      mis.z          <- error(fit.z$yhatout, dataout[,d])$mis
    }
    
    if(binary==0){
      fit.z          <- lassoF(datause=datause, dataout=dataout,  form_x, form_d, logit=FALSE, alp=alp, arg=arg, s=s)
      mis.z          <- NA
    }   
    
    fit.y            <- lassoF(datause=datause, dataout=dataout,  form_x, form_y, alp=alp, arg=arg, s=s)
    
  }    
  
  ############# Random Forest ###################################################;
  
  if(method=="Forest" | method=="TForest"){
    
    tune = FALSE
    if(method=="TForest"){tune=TRUE}
    
    option    <- arguments[[method]]
    
    arg       <- option
    arg[which(names(arg) %in% c("clas_nodesize","reg_nodesize"))] <-  NULL
    
    if(plinear==0){
      
      fit.yz1        <- RF(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, nodesize=option[["reg_nodesize"]], arg=arg, tune=tune)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      my_z1x         <- predict(fit.yz1$model, dataout, type="response") 
      
      fit.yz0        <- RF(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_y, nodesize=option[["reg_nodesize"]], arg=arg, tune=tune)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      my_z0x         <- predict(fit.yz0$model, dataout, type="response")
      
    }
    
    if(binary==1){
      fit.z          <- RF(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_d, nodesize=option[["clas_nodesize"]], arg=arg, reg=TRUE, tune=tune)
      mis.z          <- error(as.numeric(fit.z$yhatout), dataout[,y])$mis
    }
    
    if(binary==0){
      fit.z          <- RF(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_d,nodesize=option[["reg_nodesize"]], arg=arg, reg=TRUE, tune=tune)
      mis.z          <- NA
    }   
    
    fit.y          <- RF(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_y, nodesize=option[["reg_nodesize"]],  arg=arg, tune=tune)
    
  }
  
  ########################## Regression Trees ###################################################;     
  
  if(method=="Trees"){
    
    option    <- arguments[[method]]
    arg       <- option
    arg[which(names(arg) %in% c("reg_method","clas_method"))] <-  NULL
    
    if(plinear==0){
      
      fit.yz1        <- tree(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, method=option[["reg_method"]], arg=arg)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      my_z1x         <- predict(fit.yz1$model, dataout) 
      
      fit.yz0        <- tree(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_y, method=option[["reg_method"]], arg=arg)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      my_z0x         <- predict(fit.yz0$model,dataout)   
      
    }
    
    if(binary==1){
      
      fit.z          <- tree(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_d, method=option[["clas_method"]], arg=arg)
      mis.z          <- error(as.numeric(fit.z$yhatout), dataout[,y])$mis
    }
    
    if(binary==0){
      
      fit.z          <- tree(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_d, method=option[["reg_method"]], arg=arg)
      mis.z          <- NA
    }        
    
    fit.y          <- tree(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_y, method=option[["cont_method"]], arg=arg)
    
  }
  
  err.z          <- error(fit.z$yhatout, dataout[,d])$err
  mz_x           <- fit.z$yhatout       
  rz             <- fit.z$resout 
  err.z          <- error(fit.z$yhatout, dataout[,d])$err 
  
  ry             <- fit.y$resout
  my_x           <- fit.y$yhatout
  err.y          <- error(fit.y$yhatout, dataout[,y])$err
  
  return(list(my_z1x=my_z1x, mz_x= mz_x, my_z0x=my_z0x, my_x = my_x, err.z = err.z,  err.yz0= err.yz0,  err.yz1=err.yz1, mis.z=mis.z, ry=ry , rz=rz, err.y=err.y,  fit.y=fit.y, fit.yz1out= fit.yz1$yhatout,  fit.yz0out= fit.yz0$yhatout));
  
}  

ensembleF <- function(datause, dataout, y, d, xx, method, plinear, xL, binary, arguments, ensemble){
  
  K         <- 5
  k         <- length(ensemble[['methods']])
  fits      <- vector("list", k)
  method    <- ensemble[['methods']]
  ind_u     <- which(datause[,d]==1)
  ind_o     <- which(dataout[,d]==1)
  
  split     <- runif(nrow(datause))
  cvgroup   <- as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/K)),include.lowest = TRUE))  
  
  if(k<4)  {  lst <- lapply(numeric(k), function(x) as.vector(seq(0,1,0.01))) }
  if(k==4) {  lst <- lapply(numeric(k), function(x) as.vector(seq(0,1,0.02))) }
  if(k==5) {  lst <- lapply(numeric(k), function(x) as.vector(seq(0,1,0.04))) }
  if(k==6) {  lst <- lapply(numeric(k), function(x) as.vector(seq(0,1,0.1)))  }
  
  gr      <- as.matrix(expand.grid(lst))
  weight  <- gr[rowSums(gr)==1,]
  

  if(plinear==1){
    
    errorM  <- array(0,dim=c(nrow(weight),5,2))
    pred2   <- array(0,dim=c(nrow(dataout),2,k))  
    
    for(j in 1:K){
      
      ii   <- cvgroup == j
      nii  <- cvgroup != j
      
      datause1 <- as.data.frame(datause[nii,])
      datause2 <- as.data.frame(datause[ii,])  
      pred1    <- array(0,dim=c(nrow(datause2),2,k))  
      
      for(i in 1:k){
        
        if (any(c("RLasso", "PostRLasso", "Ridge", "Lasso", "Elnet")==method[i])){
          x=xL
        } else {
          x=xx
        }
        
        fits[[i]]   <- cond_comp(datause=datause1, dataout=datause2, y, d, x, method[i], plinear, xL, binary, arguments);
        
        pred1[,1,i] <- fits[[i]][['my_x']]
        pred1[,2,i] <- fits[[i]][['mz_x']]
        
      }
      
      for(p in 1:nrow(weight)){
        
        errorM[p,j,1] <- error(pred1[,1,] %*% (weight[p,]), datause2[,y])$err 
        errorM[p,j,2] <- error(pred1[,2,] %*% (weight[p,]), datause2[,d])$err 
      }
    }
    
    min1 <- which.min(as.matrix(rowSums(errorM[,,1])))
    min2 <- which.min(as.matrix(rowSums(errorM[,,2])))
    
    for(i in 1:k){
      
      if (any(c("RLasso", "PostRLasso", "Ridge", "Lasso", "Elnet")==method[i])){
        x=xL
      } else {
        x=xx
      }
      
      fits[[i]] <-cond_comp(datause=datause, dataout=dataout, y, d, x, method[i], plinear, xL, binary, arguments);
      
      pred2[,1,i] <- fits[[i]][['my_x']]
      pred2[,2,i] <- fits[[i]][['mz_x']]
      
    }
    
    fit.y <- pred2[,1,] %*% (weight[min1,])
    fit.z <- pred2[,2,] %*% (weight[min2,])
    
    ry <- dataout[,y] - fit.y
    rz <- dataout[,d] - fit.z
    
    err.y <- error(fit.y, dataout[,y])$err 
    err.z <- error(fit.z, dataout[,d])$err 
    
    return(list(err.z=err.z, err.y=err.y, ry=ry, rz=rz));
    
  }
  
  
  if(plinear==0){
    
    errorM  <- array(0,dim=c(nrow(weight),5,3))
    pred2   <- array(0,dim=c(nrow(dataout),3,k))  
    pred3   <- matrix(0,nrow(dataout[ind_o,]),k)
    pred4   <- matrix(0,nrow(dataout[-ind_o,]),k)
    
    for(j in 1:K){
      
      ii = cvgroup == j
      nii = cvgroup != j
      
      datause1 = as.data.frame(datause[nii,])
      datause2 = as.data.frame(datause[ii,])  
      
      ind2_u   <- which(datause1[,d]==1)
      ind2_o   <- which(datause2[,d]==1)
      
      pred1  <- array(0,dim=c(nrow(datause2),3,k))  

      for(i in 1:k){
        
        if (any(c("RLasso", "PostRLasso", "Ridge", "Lasso", "Elnet")==method[i])){
          x=xL
        } else {
          x=xx
        }
        
        fits[[i]] <-cond_comp(datause=datause1, dataout=datause2, y, d, x, method[i], plinear, xL, binary, arguments);
        
        pred1[,1,i] <- fits[[i]][['my_z1x']]
        pred1[,2,i] <- fits[[i]][['my_z0x']]
        pred1[,3,i] <- fits[[i]][['mz_x']]
      }
      
      for(p in 1:nrow(weight)){
        
        errorM[p,j,1] <- error(pred1[,1,] %*% (weight[p,]), datause2[ind2_o,y])$err 
        errorM[p,j,2] <- error(pred1[,2,] %*% (weight[p,]), datause2[-ind2_o,d])$err 
        errorM[p,j,3] <- error(pred1[,3,] %*% (weight[p,]), datause2[,d])$err 
      }
    }
    
    min1 <- which.min(as.matrix(rowSums(errorM[,,1])))
    min2 <- which.min(as.matrix(rowSums(errorM[,,2])))
    min3 <- which.min(as.matrix(rowSums(errorM[,,3])))
    
    for(i in 1:k){
      
      if (any(c("RLasso", "PostRLasso", "Ridge", "Lasso", "Elnet")==method[i])){
        x=xL
      } else {
        x=xx
      }
      
      fits[[i]]   <- cond_comp(datause=datause, dataout=dataout, y, d, x, method[i], plinear, xL, binary, arguments);
      
      pred2[,1,i] <- fits[[i]][['my_z1x']]
      pred2[,2,i] <- fits[[i]][['my_z0x']]
      pred2[,3,i] <- fits[[i]][['mz_x']]
      pred3[,i]   <- fits[[i]][['fit.yz1out']]
      pred4[,i]   <- fits[[i]][['fit.yz0out']]
      
    }
    
    fit.y1x    <- pred2[,1,] %*% (weight[min1,])
    fit.y0x    <- pred2[,2,] %*% (weight[min2,])
    fit.zx     <- pred2[,3,] %*% (weight[min3,])
    fit.yz1out <- pred3 %*% (weight[min1,])
    fit.yz0out <- pred4 %*% (weight[min2,])
    
    err.y1x <- error(fit.yz1out, dataout[ind_o,y])$err 
    err.y0x <- error(fit.yz0out, dataout[-ind_o,y])$err 
    err.z   <- error(fit.zx, dataout[,d])$err 
    
    return(list(err.z=err.z, err.yz1=err.y1x, err.yz0=err.y0x,  my_z1x=fit.y1x, my_z0x=fit.y0x,mz_x=fit.zx ));
    
  }
}



















