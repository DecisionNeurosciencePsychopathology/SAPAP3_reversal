#proc operant data:
setwd('~/Box/skinner/projects_analyses/Project SAPAP3/ca')
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

#Do the bins;
allmice_bin<-lapply(all_sumdt,function(mx1){
  mx1$timpointdf$endoftrial<-mx1$timpointdf$rewardRetrieval
  mx1$timpointdf$endoftrial[is.na(mx1$timpointdf$endoftrial)]<-mx1$timpointdf$allRewardDelilvery[is.na(mx1$timpointdf$endoftrial)]
  
  iti_timepoint<-round(as.numeric(unlist(sapply(2:length(mx1$timpointdf$endoftrial),function(lx){seq(from=lag(mx1$timpointdf$endoftrial,1)[lx],to=mx1$timpointdf$StartTime[lx],by = binsize) 
  }))),1)
  
  dfx<-rbind(
  data.frame(type="Correct",value=unlist(mx1$timpointdf$allCorrectResponses),stringsAsFactors = F),
  data.frame(type="Incorrect",value=unlist(mx1$timpointdf$allIncorrectResponses),stringsAsFactors = F),
  data.frame(type="MagazineEntry",value=unlist(mx1$timpointdf$allMagazineEntry),stringsAsFactors = F),
  data.frame(type="RewardDelivery",value=unlist(mx1$timpointdf$allRewardDelilvery),stringsAsFactors = F),
  data.frame(type="RewardRetrieval",value=unlist(mx1$timpointdf$rewardRetrieval),stringsAsFactors = F)
  )
  dfx<-dfx[!is.na(dfx$value),]
  
  rdx<-data.frame(timepoint=seq(from=min(dfx$value),to=max(dfx$value),by=binsize),action=0,type=NA,magzineentry=0,rewardretrieval=0)
  rtx<-do.call(rbind,lapply(c("Correct","Incorrect"),function(tx){
    rdx$action[sapply(dfx$value[dfx$type %in% c(tx)],function(x){which.min(abs(x-rdx$timepoint))})]<-1
    rdx$type[sapply(dfx$value[dfx$type %in% c(tx)],function(x){which.min(abs(x-rdx$timepoint))})]<-tx
    rdx$magzineentry[sapply(dfx$value[dfx$type %in% c("MagazineEntry")],function(x){which.min(abs(x-rdx$timepoint))})]<-1
    #rdx$rewardretrieval[sapply(dfx$value[dfx$type %in% c("RewardRetrieval")],function(x){which.min(abs(x-rdx$timepoint))})]<-1
    rdx$type[is.na(rdx$type)]<-tx
    return(rdx)
  }))
  rtx<-rtx[order(rtx$timepoint),]
  
  rtx$action[rtx$timepoint %in% iti_timepoint]<-NA
  
  rtx$magzineentry[rtx$type=="Incorrect"]<-NA
  #rtx$rewardretrieval[rtx$type=="Incorrect"]<-NA #Do we need to do that? probably, if it's incorrect then obv no reward?
  rtx$trial<-NA
  rtx$rewardretrieval<-0 
  for(ti in mx1$timpointdf$Trial){
    rtx$trial[rtx$timepoint>=mx1$timpointdf$StartTime[ti]]<-ti
    rtx$rewardretrieval[rtx$timepoint==mx1$timpointdf$rewardRetrieval[ti]]<-1
  }
  
  rtx<-cbind(unique(mx1$timpointdf[c("name","mouseID","treatment")]),rtx)
  return(rtx)
})

MiceOperantsData<-list(
  Entry_df = do.call(rbind,lapply(all_sumdt,function(ab){ab$entrydf})),
  Timpoint_df = do.call(rbind,lapply(all_sumdt,function(ab){ab$timpointdf})),
  mice_bin_df = do.call(rbind,allmice_bin)
)

save(MiceOperantsData,file = "MiceOperantsData.rdata")
