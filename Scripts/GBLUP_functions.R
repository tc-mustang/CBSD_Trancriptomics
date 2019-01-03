################################### Functions  ##################################################################

# Wrapper for multikernel emmreml 3.0 (V5: Used in Cross-validation function) ------------------------------
emmremlMultiKernel.CrossVal.mod<-function(data,train,test,trait,Klist){
  require(EMMREML)
  data1 <- data[which(data$CLONE %in% union(train,test)),c("CLONE",trait,paste(trait,"ebv",sep="."),paste(trait,"wt",sep="."))]
  data1 <- data1[!is.na(data1[,trait]),]
  rownames(data1) <- data1$CLONE
  train = train[train %in% data1$CLONE]
  test = test[test %in% data1$CLONE]
  for(kernels in 1:length(Klist)){ Klist[[kernels]]<-Klist[[kernels]][union(train,test),union(train,test)] }
  datatrain <- data1[train,]
  datatest <- data1[test,]
  datatrain$CLONE <- factor(as.character(datatrain$CLONE),levels=rownames(Klist[[1]]))
  Z = as.matrix(model.matrix(~datatrain$CLONE-1))
  y = datatrain[,trait]
  
  ################ ACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
  X = model.matrix(~1,data=datatrain)
  weight<-sqrt(data1[train,paste(trait,"wt",sep=".")])
  Z = Z/weight
  X = X/weight
  y = y/weight
  Zlist = list()
  for(kernels in 1:length(Klist)){ Zlist[[paste("Z",kernels,sep="")]] <- Z }
  if(length(Klist) == 1){ 
    funout <- emmreml(y=y, X=X, Z=Zlist[[1]], K=Klist[[1]]) } else {
      funout <- emmremlMultiKernel(y=y, X=X, Z=Zlist, K=Klist, varbetahat=F,varuhat=F, PEVuhat=F, test=F) }
  uhatmat<-matrix(funout$uhat, ncol=length(Klist))
  rownames(uhatmat)<-rownames(Klist[[1]])
  totvals<-rowSums(uhatmat,na.rm=T)
  names(totvals)<-rownames(Klist[[1]])
  funout[["Train"]]<-train
  funout[["Test"]]<-test
  funout[["Trait"]]<-trait    
  funout[["uhatmat"]]<-uhatmat 
  funout[["TotVals"]]<-totvals 
  funout[["Dataset"]]<-data1
  return(funout)
}



# Fold cross-validation function for emmreml multikernel 3.0 (V3: saves loglik + BLUPs for each rep) --------------------------------------------------------------------
FoldCrossValidation.V3.emmreml <- function(dataset, trait, genoID, Klist, nFolds, nRepeats){
  data<-dataset[which(dataset[ ,genoID] %in% rownames(Klist[[1]])),]
  data<-data[!is.na(data[,trait]),]
  rownames(data)<-data$CLONE
  for(kernels in 1:length(Klist)){
    K<-Klist[[kernels]]
    K<-K[rownames(K) %in% data[ ,genoID], rownames(K) %in% data[ ,genoID]]
    Klist[[kernels]]<-K }
  nInd <- dim(data)[1] 
  accuracies<-data.frame()
  varcomps<-data.frame()
  blups<-data.frame()
  for (rep in 1:nRepeats){ 
    print(paste("Rep ",rep," of ",nRepeats,sep=""))
    folds <- sample(rep(1:nFolds, length.out=nInd))
    BLUPSthisRep<-data.frame()
    varcomps.thisrep<-data.frame()
    for (fold in 1:nFolds){
      print(paste("Fold ",fold," of ",nFolds,sep=""))
      indInFold <- which(folds == fold)
      indNotInFold <- which(folds != fold)
      ValidSet<-data[indInFold,genoID]
      TrainSet<-data[indNotInFold,genoID]
      out<-emmremlMultiKernel.CrossVal.mod(data,train=TrainSet,test=ValidSet,trait,Klist=Klist)
      out[["AIC"]]<-2*length(Klist)-2*out$loglik
      if(length(Klist)>1){
        colnames(out$uhatmat)<-paste("K",1:length(Klist),sep="")
        BLUPSthisFold<-data.frame(CLONE=rownames(out$uhatmat[ValidSet,]))
        BLUPSthisFold<-cbind(BLUPSthisFold,out$uhatmat[ValidSet,])
        BLUPSthisRep<-rbind(BLUPSthisRep,BLUPSthisFold)
        names(out$weights)<-colnames(out$uhatmat)        
        varcomps.thisfold<-data.frame(Trait=trait,Rep=rep,Fold=fold,Vu=out$Vu, Ve=out$Ve, loglik=out$loglik, AIC=out$AIC)
        varcomps.thisfold<-cbind(varcomps.thisfold,t(out$weights)) } 
      else { 
        BLUPSthisFold<-data.frame(CLONE=names(out$uhatmat[ValidSet,]))
        BLUPSthisFold<-cbind(BLUPSthisFold,out$uhatmat[ValidSet,])
        colnames(BLUPSthisFold)<-c(genoID,"K1")
        BLUPSthisRep<-rbind(BLUPSthisRep,BLUPSthisFold)
        varcomps.thisfold<-data.frame(Trait=trait,Rep=rep,Fold=fold,Vu=out$Vu, Ve=out$Ve, loglik=out$loglik, AIC=out$AIC) 
      } 
      varcomps.thisrep<-rbind(varcomps.thisrep,varcomps.thisfold)   
    }
    BLUPSthisRep[,"Rep"]<-rep
    # Calc accuracy after predicting each fold to complete the rep
    if(length(Klist)>1){
      BLUPSthisRep[,"TotGen"]<-rowSums(BLUPSthisRep[,c(paste("K",1:length(Klist),sep=""))],na.rm=T)
      BLUPSthisRep<-merge(BLUPSthisRep,data[,c("CLONE",paste(trait,".ebv",sep=""))],by="CLONE")
      blup.names<-c(paste("K",1:length(Klist),sep=""),"TotGen")
      accuracies.thisrep<-vector()
      for(tacos in 1:length(blup.names)){
        acc = cor(BLUPSthisRep[,blup.names[tacos]],BLUPSthisRep[,paste(trait,".ebv",sep="")], use="complete.obs")
        accuracies.thisrep<-c(accuracies.thisrep,acc) }
      names(accuracies.thisrep)<-paste("Acc.",blup.names,sep="")
      AccuracyThisRep<-data.frame(Trait=trait,Rep=rep)
      AccuracyThisRep<-cbind(AccuracyThisRep,t(accuracies.thisrep))
      accuracies<-rbind(accuracies,AccuracyThisRep) 
      varcomps<-rbind(varcomps,varcomps.thisrep)
    } else {  
      BLUPSthisRep<-merge(BLUPSthisRep,data[,c("CLONE",paste(trait,".ebv",sep=""))],by="CLONE")
      blup.names<-c(paste("K",1:length(Klist),sep=""))
      accuracies.thisrep<-vector()
      for(tacos in 1:length(blup.names)){
        acc = cor(BLUPSthisRep[,blup.names[tacos]],BLUPSthisRep[,paste(trait,".ebv",sep="")], use="complete.obs")
        accuracies.thisrep<-c(accuracies.thisrep,acc) }
      names(accuracies.thisrep)<-paste("Acc.",blup.names,sep="")
      AccuracyThisRep<-data.frame(Trait=trait,Rep=rep)
      AccuracyThisRep<-cbind(AccuracyThisRep,t(accuracies.thisrep))
      accuracies<-rbind(accuracies,AccuracyThisRep) 
      varcomps<-rbind(varcomps,varcomps.thisrep)
    }
    blups<-rbind(blups,BLUPSthisRep)
  }
  return(list(accuracies=accuracies,varcomps=varcomps, blups=blups))
}

