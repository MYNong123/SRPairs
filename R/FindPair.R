#' FindStabRvsPair
#'
#' @param ControlFrame Expression matrix of control group.
#' @param TreatFrame Expression matrix of treat group.
#' @param HighPercent Percentage of samples with higher expression of Gene i than Gene j in total samples.
#' @param LowPercent Percentage of samples with lower expression of Gene i than Gene j in total samples.
#'
#' @return The stable reversal gene pairs were returned by inputting the control group and the experimental group.
#' @export
#'
#' @examples FindStabRvsPair(ControlFrame,TreatFrame,HighPercent=0.9,ncore=1,showbar=T)
FindStabRvsPair<-function(ControlFrame,TreatFrame,HighPercent=0.9,ncore=1,showbar=T){
  #Entity operation function.
  LowPercent=1-HighPercent
  library(doSNOW)
  print("WeChat official account: shengxinjy123")
  if (ncore==1) {
    #Read the gene expression file of Control group.
    sampleNum=ncol(ControlFrame)
    ControlData=data.frame()
    print("Calculation of control.")
    
    for(i in 1:(nrow(ControlFrame)-1)){
      for(j in (i+1):nrow(ControlFrame)){
        pair=ifelse(ControlFrame[i,]<ControlFrame[j,],0,1)
        pairRatio=sum(pair)/sampleNum
        if(pairRatio>HighPercent|pairRatio==HighPercent){
          StablePair=paste0(rownames(ControlFrame)[i],"&",rownames(ControlFrame)[j])      #Stable Gi>Gj gene pair.
          ControlData=as.matrix(rbind(ControlData,c(StablePair,1)))
        }
        if(pairRatio<LowPercent){
          StablePair=paste0(rownames(ControlFrame)[i],"&",rownames(ControlFrame)[j])      #Stable Gi<Gj gene pair.
          ControlData=as.matrix(rbind(ControlData,c(StablePair,0)))
        }
      }
      print(paste0(Sys.time()," ","Control:There are ",(nrow(ControlFrame)-(i+1))," data left unprocessed."))
    }
    
    colnames(ControlData)=c("StablePair1","Index1")
    
    #Read the gene expression file of Control group.
    sampleNum2=ncol(TreatFrame)
    TreatData=data.frame()
    
    print("Calculation of treat.")
    
    for(i in 1:(nrow(TreatFrame)-1)){
      for(j in (i+1):nrow(TreatFrame)){
        pair=ifelse(TreatFrame[i,]<TreatFrame[j,],0,1)
        pairRatio=sum(pair)/sampleNum2
        if(pairRatio>HighPercent|pairRatio==HighPercent){
          StablePair=paste0(rownames(TreatFrame)[i],"&",rownames(TreatFrame)[j])      #Stable Gi>Gj gene pair.
          TreatData=as.matrix(rbind(TreatData,c(StablePair,1)))
        }
        if(pairRatio<LowPercent){
          StablePair=paste0(rownames(TreatFrame)[i],"&",rownames(TreatFrame)[j])      #Stable Gi<Gj gene pair.
          TreatData=as.matrix(rbind(TreatData,c(StablePair,0)))
        }
      }
      print(paste0(Sys.time()," ","Treat:There are ",(nrow(TreatFrame)-(i+1))," data left unprocessed."))
    }
    
    colnames(TreatData)=c("StablePair2","Index2")
    
    #Search for stable reversal gene pairs
    rownames(ControlData)=ControlData[,1]
    rownames(TreatData)=TreatData[,1]
    genemix=intersect(ControlData[,1],TreatData[,1])
    Treatsect=TreatData[genemix,]
    Controlsect=ControlData[genemix,]
    InputData=cbind(Controlsect,Treatsect)
    InputData=InputData[,c(-1,-3)]
    InputData<-as.data.frame(InputData)
    InputData$Index3<-ifelse(InputData$Index1==InputData$Index2,0,1)
    InputPair=InputData[which(InputData$Index3==1),]
    colnames(InputPair)=c("Control","Treat","Pair")
    print("WeChat official account: shengxinjy123")
    
    return(InputPair)
  } else if (showbar==T) {
    
    #Read the gene expression file of Control group.
    sampleNum=ncol(ControlFrame)
    ControlData=data.frame()
    
    print("Calculation of control.")
    
    fun1<-function(x){
      for(i in x){
        for(j in (i+1):nrow(ControlFrame)){
          pair=ifelse(ControlFrame[i,]<ControlFrame[j,],0,1)
          pairRatio=sum(pair)/sampleNum
          if(pairRatio>HighPercent|pairRatio==HighPercent){
            StablePair=paste0(rownames(ControlFrame)[i],"&",rownames(ControlFrame)[j])      #Stable Gi>Gj gene pair.
            ControlData=as.matrix(rbind(ControlData,c(StablePair,1)))
          }
          if(pairRatio<LowPercent){
            StablePair=paste0(rownames(ControlFrame)[i],"&",rownames(ControlFrame)[j])      #Stable Gi<Gj gene pair.
            ControlData=as.matrix(rbind(ControlData,c(StablePair,0)))
          }
        }
        return(ControlData)
        print(paste0(Sys.time()," ","Control:There are ",(nrow(ControlFrame)-(i+1))," data left unprocessed."))
      }
    }
    
    cl <- makeCluster(getOption("cl.cores", ncore))

    clusterExport(cl,"sampleNum",envir = environment())
    clusterExport(cl,"ControlData",envir = environment())
    clusterExport(cl,"ControlFrame",envir = environment())
    clusterExport(cl,"HighPercent",envir = environment())
    clusterExport(cl,"LowPercent",envir = environment())
    
    registerDoSNOW(cl)
    
    pb <- txtProgressBar(max=(nrow(ControlFrame)-1), style=3, char = "*",)
    on.exit(close(pb))
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    
    library(doParallel)

    res<-foreach(i=1:(nrow(ControlFrame)-1), .options.snow=opts) %dopar% {
      Sys.sleep(0.1)
      fun1(i)
    }
    
    on.exit(close(pb))
    
    stopCluster(cl)
    
    for (i in 1:(nrow(ControlFrame)-1)) {
      if (ncol(res[[i]])==2){
        colnames(res[[i]])<-c("StablePair1","Index1")
        ControlData=rbind(ControlData,res[[i]])
      }
    }
    
    print(" ")
    print("Calculation of treat.")
    
    #Read the gene expression file of Control group.
    sampleNum2=ncol(TreatFrame)
    TreatData=data.frame()
    
    fun2<-function(x){
      for(i in x){
        for(j in (i+1):nrow(TreatFrame)){
          pair=ifelse(TreatFrame[i,]<TreatFrame[j,],0,1)
          pairRatio=sum(pair)/sampleNum2
          if(pairRatio>HighPercent|pairRatio==HighPercent){
            StablePair=paste0(rownames(TreatFrame)[i],"&",rownames(TreatFrame)[j])      #Stable Gi>Gj gene pair.
            TreatData=as.matrix(rbind(TreatData,c(StablePair,1)))
          }
          if(pairRatio<LowPercent){
            StablePair=paste0(rownames(TreatFrame)[i],"&",rownames(TreatFrame)[j])      #Stable Gi<Gj gene pair.
            TreatData=as.matrix(rbind(TreatData,c(StablePair,0)))
          }
        }
        return(TreatData)
        print(paste0(Sys.time()," ","Treat:There are ",(nrow(TreatFrame)-(i+1))," data left unprocessed."))
      }
    }
    
    cl <- makeCluster(getOption("cl.cores", ncore))

    clusterExport(cl,"sampleNum2",envir = environment())
    clusterExport(cl,"TreatData",envir = environment())
    clusterExport(cl,"TreatFrame",envir = environment())
    clusterExport(cl,"HighPercent",envir = environment())
    clusterExport(cl,"LowPercent",envir = environment())
    
    registerDoSNOW(cl)
    
    pb <- txtProgressBar(max=(nrow(TreatFrame)-1), style=3, char = "*",)
    on.exit(close(pb))
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    
    res<-foreach(i=1:(nrow(TreatFrame)-1), .options.snow=opts) %dopar% {
      Sys.sleep(0.1)
      fun2(i)
    }
    
    on.exit(close(pb))
    
    stopCluster(cl)
    
    for (i in 1:(nrow(TreatFrame)-1)) {
      if (ncol(res[[i]])==2){
        colnames(res[[i]])<-c("StablePair2","Index2")
        TreatData=rbind(TreatData,res[[i]])
      }
    }
    
    #Search for stable reversal gene pairs
    rownames(ControlData)=ControlData[,1]
    rownames(TreatData)=TreatData[,1]
    genemix=intersect(ControlData[,1],TreatData[,1])
    Treatsect=TreatData[genemix,]
    Controlsect=ControlData[genemix,]
    InputData=cbind(Controlsect,Treatsect)
    InputData=InputData[,c(-1,-3)]
    InputData<-as.data.frame(InputData)
    InputData$Index3<-ifelse(InputData$Index1==InputData$Index2,0,1)
    InputPair=InputData[which(InputData$Index3==1),]
    
    colnames(InputPair)=c("Control","Treat","Pair")
    
    return(InputPair)
    print("WeChat official account: shengxinjy123")
  } else {
    library(parallel)
    #Read the gene expression file of Control group.
    sampleNum=ncol(ControlFrame)
    ControlData=data.frame()
    
    print("Calculation of control.")
    
    fun1<-function(x){
      for(i in x){
        for(j in (i+1):nrow(ControlFrame)){
          pair=ifelse(ControlFrame[i,]<ControlFrame[j,],0,1)
          pairRatio=sum(pair)/sampleNum
          if(pairRatio>HighPercent|pairRatio==HighPercent){
            StablePair=paste0(rownames(ControlFrame)[i],"&",rownames(ControlFrame)[j])      #Stable Gi>Gj gene pair.
            ControlData=as.matrix(rbind(ControlData,c(StablePair,1)))
          }
          if(pairRatio<LowPercent){
            StablePair=paste0(rownames(ControlFrame)[i],"&",rownames(ControlFrame)[j])      #Stable Gi<Gj gene pair.
            ControlData=as.matrix(rbind(ControlData,c(StablePair,0)))
          }
        }
        return(ControlData)
        
      }
    }
    
    cl <- makeCluster(getOption("cl.cores", ncore))
    
    clusterExport(cl,"sampleNum",envir = environment())
    clusterExport(cl,"ControlData",envir = environment())
    clusterExport(cl,"ControlFrame",envir = environment())
    clusterExport(cl,"HighPercent",envir = environment())
    clusterExport(cl,"LowPercent",envir = environment())
    
    res <- parLapply(cl,1:(nrow(ControlFrame)-1), fun1)

    stopCluster(cl)
    
    for (i in 1:(nrow(ControlFrame)-1)) {
      if (ncol(res[[i]])==2){
        colnames(res[[i]])<-c("StablePair1","Index1")
        ControlData=rbind(ControlData,res[[i]])
      }
    }
    
    #Read the gene expression file of Control group.
    sampleNum2=ncol(TreatFrame)
    TreatData=data.frame()
    
    print("Calculation of treat.")
    
    fun2<-function(x){
      for(i in x){
        for(j in (i+1):nrow(TreatFrame)){
          pair=ifelse(TreatFrame[i,]<TreatFrame[j,],0,1)
          pairRatio=sum(pair)/sampleNum2
          if(pairRatio>HighPercent|pairRatio==HighPercent){
            StablePair=paste0(rownames(TreatFrame)[i],"&",rownames(TreatFrame)[j])      #Stable Gi>Gj gene pair.
            TreatData=as.matrix(rbind(TreatData,c(StablePair,1)))
          }
          if(pairRatio<LowPercent){
            StablePair=paste0(rownames(TreatFrame)[i],"&",rownames(TreatFrame)[j])      #Stable Gi<Gj gene pair.
            TreatData=as.matrix(rbind(TreatData,c(StablePair,0)))
          }
        }
        return(TreatData)
        print(paste0(Sys.time()," ","Treat:There are ",(nrow(TreatFrame)-(i+1))," data left unprocessed."))
      }
    }
    
    cl <- makeCluster(getOption("cl.cores", ncore))
    
    clusterExport(cl,"sampleNum2",envir = environment())
    clusterExport(cl,"TreatData",envir = environment())
    clusterExport(cl,"TreatFrame",envir = environment())
    clusterExport(cl,"HighPercent",envir = environment())
    clusterExport(cl,"LowPercent",envir = environment())
    
    res <- parLapply(cl,1:(nrow(TreatFrame)-1), fun2)
    
    stopCluster(cl)
    
    for (i in 1:(nrow(TreatFrame)-1)) {
      if (ncol(res[[i]])==2){
        colnames(res[[i]])<-c("StablePair2","Index2")
        TreatData=rbind(TreatData,res[[i]])
      }
    }
    
    #Search for stable reversal gene pairs
    rownames(ControlData)=ControlData[,1]
    rownames(TreatData)=TreatData[,1]
    genemix=intersect(ControlData[,1],TreatData[,1])
    Treatsect=TreatData[genemix,]
    Controlsect=ControlData[genemix,]
    InputData=cbind(Controlsect,Treatsect)
    InputData=InputData[,c(-1,-3)]
    InputData<-as.data.frame(InputData)
    InputData$Index3<-ifelse(InputData$Index1==InputData$Index2,0,1)
    InputPair=InputData[which(InputData$Index3==1),]
    colnames(InputPair)=c("Control","Treat","Pair")
    
    return(InputPair)
    print("WeChat official account: shengxinjy123")
  }
}
