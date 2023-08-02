#' CalStabRvsPair
#'
#' @param FindPairFrame The Pair Frame.
#' @param ControlFrame Expression matrix of control group.
#' @param TreatFrame Expression matrix of treat group.
#'
#' @return Sort out the sample-gene pair matrix by reversing the gene pair matrix.
#' @export
#'
#' @examples CalStabRvsPair(FindPairFrame,ControlFrame,TreatFrame)
CalStabRvsPair<-function(FindPairFrame,ControlFrame,TreatFrame){
  library(stringr)
  print("'stringr' package must be loaded! Otherwise, this function cannot be run!")
  InputPairs = FindPairFrame
  FacPair<-rownames(InputPairs)
  GenePair <- t(data.frame(str_split(FacPair,'\\&')))
  colnames(GenePair) <- c('Gene1','Gene2')
  rownames(GenePair) <- rownames(InputPairs)

  PairGeneIndex<-data.frame()

  for (i in 1:nrow(GenePair)) {
    a=as.vector(GenePair[i,])
    b=rbind(cbind(ControlFrame[a[1],],TreatFrame[a[1],]),cbind(ControlFrame[a[2],],TreatFrame[a[2],]))
    b[3,]<-ifelse(b[1,]>b[2,],1,0)
    b<-b[-1:-2,]
    row.names(b)=paste0(a[1],"|",a[2])
    PairGeneIndex=rbind(PairGeneIndex,b)
  }
  print("WeChat official account: shengxinjy123")

  return(PairGeneIndex)
}
