#proc operant data:
setwd('~/Box Sync/skinner/projects_analyses/Project SAPAP3/ca')
library(R.matlab)
binsize = .1 # start at .1 seconds
if(!file.exists("operantData.rdata")){
  dx<-readMat("operantData.mat")
  save(dx,file = "operantData.rdata")
} else { load (file = "operantData.rdata")}

#import this recursive function to deal with matlab stuff
rc_mat<-function(matdt){
  mx_dt<-lapply(drop(matdt), drop)
  if( any(sapply(mx_dt,is.list)) ) {
    mx_dt<-c(mx_dt[!sapply(mx_dt,is.list)],lapply(mx_dt[sapply(mx_dt,is.list)],rc_mat))
  } else {return(mx_dt)}
}

all_sumdt<-lapply(1:dim(dx$operantData)[3],function(mx){
  mx_dt<-rc_mat(dx$operantData[,,mx])
  #message(mx)
  sum_ls<-lapply(1:length(mx_dt$discreteBehavior$allTrialStart$timeStamp),function(trialn){
    #message(trialn)
    startn <- mx_dt$discreteBehavior$allTrialStart$timeStamp[trialn]
    endn <- ifelse(trialn<length(mx_dt$discreteBehavior$allTrialStart$timeStamp), mx_dt$discreteBehavior$allTrialStart$timeStamp[trialn+1],Inf)
    zt<-lapply(mx_dt$discreteBehavior[!grepl(paste("allTrialStart","name",sep = "|"),names(mx_dt$discreteBehavior))],function(tx){
      tps<-tx$timeStamp[tx$timeStamp >= startn & tx$timeStamp < endn]
    }) 
    
    dfx<-data.frame(Trial=trialn,StartTime=startn)
    for (nx in names(zt)){
      dfx[[nx]]<-ifelse(length(zt[[nx]])>1,list(zt[[nx]]),zt[[nx]])
    }  
    dfx$name<-mx_dt$name
    dfx$treatment<-mx_dt$treatment
    dfx$mouseID<-mx_dt$mouseID
    
    if(length(zt$allMagazineEntry)<1){zt$allMagazineEntry<- -999;noentrytag<-T}
    edf<-do.call(rbind,lapply(1:length(zt$allMagazineEntry),function(entryc){
      #print(entryc)
      enter <-zt$allMagazineEntry[entryc]
      nextenter <- ifelse(entryc<length(zt$allMagazineEntry), zt$allMagazineEntry[entryc+1],Inf)
      re<-lapply(zt[-grep("allMagazineEntry",names(zt))],function(xe){
        ge<-xe[xe>=enter & xe < nextenter]
        if(length(ge)<1){ge<-NA}
        return(ge)
      })
      re_df<-do.call(rbind,lapply(c("allCorrectResponses","allIncorrectResponses"),function(rx){
        if( any(!is.na(re[[rx]])) ){
        dfax<-data.frame(responsetimepoint=re[[rx]],respondtype=rx)
        } else {
        dfax<-NULL  
        }
        }))
      if(is.null(re_df)){re_df<-data.frame(responsetimepoint=NA,respondtype="NoResponses")}
      re_df$respondtype<-gsub("all","",re_df$respondtype)
      re_df$RewardDeliberyTime<-re$allRewardDelilvery
      re_df$RewardRetrievalTime<-re$rewardRetrieval
      re_df$MagazineEntryTime<-enter
      re_df$MagazineEntryCount<-entryc
      re_df$Trial<-trialn
      re_df$TrialStartTime<-startn
      re_df$name<-mx_dt$name
      re_df$treatment<-mx_dt$treatment
      re_df$mouseID<-mx_dt$mouseID
      return(re_df)
    })
    )
    return(list(timpointdf=dfx,entrydf=edf))
  })
  
  return(list(
    entrydf = do.call(rbind,lapply(sum_ls,function(ab){ab$entrydf})),
    timpointdf = do.call(rbind,lapply(sum_ls,function(ab){ab$timpointdf}))
  ))
})

MiceOperantsData<-list(
  Entry_df = do.call(rbind,lapply(all_sumdt,function(ab){ab$entrydf})),
  Timpoint_df = do.call(rbind,lapply(all_sumdt,function(ab){ab$timpointdf}))
)

save(MiceOperantData,file = "MiceOperantsData.rdata")
