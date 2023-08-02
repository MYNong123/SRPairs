#' testRunTime
#'
#' @param ControlFrame Expression matrix of control group.
#' @param TreatFrame Expression matrix of treat group.
#' @param HighPercent Percentage of samples with higher expression of Gene i than Gene j in total samples.
#' @param LowPercent Percentage of samples with lower expression of Gene i than Gene j in total samples.
#'
#' @return The approximate time required for this computer to return the data.
#' @export
#'
#' @examples testRunTime(ControlFrame,TreatFrame,HighPercent=0.9,ncore=1)
testRunTime<-function(ControlFrame,TreatFrame,HighPercent=0.9,ncore=1){
  #Test the run time required for your computer.
  LowPercent=1-HighPercent
  #Read the gene expression file of Control group.
  sampleNum=ncol(ControlFrame)
  ControlData=data.frame()

  i=1
  Time1=system.time(
    for(j in (i+1):nrow(ControlFrame)){
      pair=ifelse(ControlFrame[i,]>ControlFrame[j,],1,0)
      pairRatio=sum(pair)/sampleNum
      if(pairRatio>HighPercent){
        StablePair=paste0(rownames(ControlFrame)[i],"&",rownames(ControlFrame)[j])      #Stable Gi>Gj gene pair.
        ControlData=as.matrix(rbind(ControlData,c(StablePair,1)))
      }
      if(pairRatio<LowPercent){
        StablePair=paste0(rownames(ControlFrame)[i],"&",rownames(ControlFrame)[j])      #Stable Gi<Gj gene pair.
        ControlData=as.matrix(rbind(ControlData,c(StablePair,0)))
      }
    })
  print(paste0("Control:The calculation time takes ",round(Time1[3]*nrow(ControlFrame))/3600," hours."))

  #Read the gene expression file of Control group.
  sampleNum2=ncol(TreatFrame)
  TreatData=data.frame()

  i=1
  Time2=system.time(
    for(j in (i+1):nrow(TreatFrame)){
      pair=ifelse(TreatFrame[i,]>TreatFrame[j,],1,0)
      pairRatio=sum(pair)/sampleNum2
      if(pairRatio>HighPercent){
        StablePair=paste0(rownames(TreatFrame)[i],"&",rownames(TreatFrame)[j])      #Stable Gi>Gj gene pair.
        TreatData=as.matrix(rbind(TreatData,c(StablePair,1)))
      }
      if(pairRatio<LowPercent){
        StablePair=paste0(rownames(TreatFrame)[i],"&",rownames(TreatFrame)[j])      #Stable Gi<Gj gene pair.
        TreatData=as.matrix(rbind(TreatData,c(StablePair,0)))
      }
    })
  print(paste0("Treat:The calculation time takes ",round(Time2[3]*nrow(TreatFrame))/3600*ncore," hours."))
  print("WeChat official account: shengxinjy123")
}
