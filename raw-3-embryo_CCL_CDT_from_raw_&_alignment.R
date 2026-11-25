


############################################################################
#					extract raw data
############################################################################

focal_projectFolderPath<-r"(E:\canalizationENV)"

##

setwd(paste0(focal_projectFolderPath,"/data/basic"))

load(paste0(focal_projectFolderPath,"/data/basic/q1_cell350_alignment_xyz.RData"))
load(paste0(focal_projectFolderPath,"/data/alignment/embryo_duration_prcomp.RData"))

###

#ccl_mt

ccl_mt<-read.csv(paste0(focal_projectFolderPath,"/raw_data/Raw_data_RNAi/temporal/CCL_RNAi.csv"),header=F)
ccl_mt_colnames<-as.character(ccl_mt[1,-1])
ccl_mt_rownames<-as.character(ccl_mt[-1,1])
ccl_mt<-read.csv(paste0(focal_projectFolderPath,"/raw_data/Raw_data_RNAi/temporal/CCL_RNAi.csv"))
ccl_mt<-as.matrix(ccl_mt[,-1])
colnames(ccl_mt)<-ccl_mt_colnames
rownames(ccl_mt)<-ccl_mt_rownames

save(ccl_mt,file="ccl_mt.RData")

#ccl_wt

ccl_wt<-read.csv(paste0(focal_projectFolderPath,"/raw_data/Raw_data_Control/temporal/CCL_control.csv"),header=F)
ccl_wt_colnames<-as.character(ccl_wt[1,-1])
ccl_wt_rownames<-as.character(ccl_wt[-1,1])
ccl_wt<-read.csv(paste0(focal_projectFolderPath,"/raw_data/Raw_data_Control/temporal/CCL_control.csv"))
ccl_wt<-as.matrix(ccl_wt[,-1])
colnames(ccl_wt)<-ccl_wt_colnames
rownames(ccl_wt)<-ccl_wt_rownames

save(ccl_wt,file="ccl_wt.RData")

#ccl_env

ccl_env<-read.csv(paste0(focal_projectFolderPath,"/raw_data/Raw_data_Env/temporal/CCL_Env.csv"),header=F)
ccl_env_colnames<-as.character(ccl_env[1,-1])
ccl_env_rownames<-as.character(ccl_env[-1,1])
ccl_env<-read.csv(paste0(focal_projectFolderPath,"/raw_data/Raw_data_Env/temporal/CCL_Env.csv"))
ccl_env<-as.matrix(ccl_env[,-1])
colnames(ccl_env)<-ccl_env_colnames
rownames(ccl_env)<-ccl_env_rownames

save(ccl_env,file="ccl_env.RData")

#par(mfrow=c(1,2))
#plot(apply(ccl_wt,1,mean,na.rm=T),apply(ccl_mt,1,mean,na.rm=T))
#abline(0,1)
#plot(apply(ccl_wt,1,mean,na.rm=T),apply(ccl_env,1,mean,na.rm=T))
#abline(0,1)

dataTemp<-cbind(ccl_wt,ccl_mt,ccl_env)
ccl_raw_nonNA<-t(dataTemp[,apply(dataTemp,2,function(x){ifelse(length(which(is.na(x)))>0,FALSE,TRUE)})])
ccl_raw_nonNA<-ccl_raw_nonNA[na.omit(match(rownames(q1_cell350_alignment_xyz),rownames(ccl_raw_nonNA))),]

colnames(ccl_raw_nonNA)<-paste0(colnames(ccl_raw_nonNA),"_ccl")

save(ccl_raw_nonNA,file=paste0(focal_projectFolderPath,"/data/basic/ccl_raw_nonNA.RData"))


### alignment 

# scaling factor by wt58

embryo_duration_prcomp<-sapply(1:dim(ccl_raw_nonNA)[1],function(i){
	temp<-prcomp(t(ccl_raw_nonNA[c(58,i),]),center=FALSE)$rotation[,1]
	return(temp[2]/temp[1])
})


save(embryo_duration_prcomp,file=paste0(focal_projectFolderPath,"/data/alignment/embryo_duration_prcomp.RData"))

# scaling factor by wtMean

wtMean<-apply(ccl_raw_nonNA[1:105,],2,mean)

embryo_duration_prcomp_wtMean<-sapply(1:dim(ccl_raw_nonNA)[1],function(i){
	focalDT<-t(ccl_raw_nonNA[c(58,i),])
	focalDT[,1]<-wtMean
	temp<-prcomp(focalDT,center=FALSE)$rotation[,1]
	return(temp[2]/temp[1])
})

save(embryo_duration_prcomp_wtMean,file=paste0(focal_projectFolderPath,"/data/alignment/embryo_duration_prcomp_wtMean.RData"))

# scaling

ccl_raw_nonNA_scaling<-t(sapply(1:dim(ccl_raw_nonNA)[1],function(i){as.numeric(ccl_raw_nonNA[i,]/embryo_duration_prcomp[i])}))
rownames(ccl_raw_nonNA_scaling)<-rownames(ccl_raw_nonNA)
colnames(ccl_raw_nonNA_scaling)<-colnames(ccl_raw_nonNA)

ccl_raw_nonNA_scaling_noise<-scale(ccl_raw_nonNA_scaling,scale=apply(ccl_raw_nonNA_scaling,2,mean))

save(ccl_raw_nonNA_scaling,file=paste0(focal_projectFolderPath,"/data/basic/ccl_raw_nonNA_scaling.RData"))
save(ccl_raw_nonNA_scaling_noise,file=paste0(focal_projectFolderPath,"/data/basic/ccl_raw_nonNA_scaling_noise.RData"))

ccl_raw_nonNA_scaling_mean<-apply(ccl_raw_nonNA_scaling,2,mean)
save(ccl_raw_nonNA_scaling_mean,file=paste0(focal_projectFolderPath,"/data/basic/ccl_raw_nonNA_scaling_mean.RData"))
ccl_raw_nonNA_scaling_scale<-apply(ccl_raw_nonNA_scaling,2,sd)
save(ccl_raw_nonNA_scaling_scale,file=paste0(focal_projectFolderPath,"/data/basic/ccl_raw_nonNA_scaling_scale.RData"))


############################################################################################3

############################################################################################3

#cdt_wt

cdt_wt<-read.csv(paste0(focal_projectFolderPath,"/raw_data/Raw_data_Control/temporal/cell_division_time_wt.csv"),header=F)
cdt_wt_colnames<-as.character(cdt_wt[1,-1])
cdt_wt_rownames<-as.character(cdt_wt[-1,1])
cdt_wt<-read.csv(paste0(focal_projectFolderPath,"/raw_data/Raw_data_Control/temporal/cell_division_time_wt.csv"))
cdt_wt<-as.matrix(cdt_wt[,-1])
colnames(cdt_wt)<-cdt_wt_colnames
rownames(cdt_wt)<-cdt_wt_rownames


#write.csv(cdt_wt,paste0(focal_projectFolderPath,"/raw_data/Raw_data_Control/temporal/cdt_control.csv"))

save(cdt_wt,file="cdt_wt.RData")



#cdt_mt (artificial modification two RNAi embryos)

cdt_mt<-read.csv(paste0(focal_projectFolderPath,"/raw_data/Raw_data_RNAi/temporal/cell_division_time_mt_revise.csv"),header=F)
cdt_mt_colnames<-as.character(cdt_mt[1,-1])
cdt_mt_rownames<-as.character(cdt_mt[-1,1])
cdt_mt<-read.csv(paste0(focal_projectFolderPath,"/raw_data/Raw_data_RNAi/temporal/cell_division_time_mt_revise.csv"))
cdt_mt<-as.matrix(cdt_mt[,-1])
colnames(cdt_mt)<-cdt_mt_colnames
rownames(cdt_mt)<-cdt_mt_rownames


cdt_mt<-cdt_mt[match(rownames(cdt_wt),rownames(cdt_mt)),]
cdt_mt<-cdt_mt[,match(colnames(ccl_mt),colnames(cdt_mt))]
#write.csv(cdt_mt,paste0(focal_projectFolderPath,"/raw_data/Raw_data_RNAi/temporal/cdt_RNAi.csv"))

save(cdt_mt,file="cdt_mt.RData")



#cdt_env

cdt_env<-read.table(paste0(focal_projectFolderPath,"/raw_data/Raw_data_Env/temporal/Cell_division_timing_matrix_environment.txt"),header=F,sep="\t")
cdt_env_colnames<-as.character(cdt_env[1,-1])
cdt_env_rownames<-as.character(cdt_env[-1,1])
cdt_env<-read.table(paste0(focal_projectFolderPath,"/raw_data/Raw_data_Env/temporal/Cell_division_timing_matrix_environment.txt"),header=T,sep="\t")
cdt_env<-as.matrix(cdt_env[,-1])
colnames(cdt_env)<-cdt_env_colnames
rownames(cdt_env)<-cdt_env_rownames

cdt_env<-cdt_env[match(rownames(cdt_wt),rownames(cdt_env)),]
cdt_env<-cdt_env[,match(colnames(ccl_env),colnames(cdt_env))]
#write.csv(cdt_env,paste0(focal_projectFolderPath,"/raw_data/Raw_data_Env/temporal/cdt_Env.csv"))

save(cdt_env,file="cdt_env.RData")



# combine 

dataTemp<-cbind(cdt_wt,cdt_mt,cdt_env)
cdt_raw_nonNA<-t(dataTemp[,apply(dataTemp,2,function(x){ifelse(length(which(is.na(x)))>0,FALSE,TRUE)})])
cdt_raw_nonNA<-cdt_raw_nonNA[na.omit(match(rownames(q1_cell350_alignment_xyz),rownames(cdt_raw_nonNA))),]

colnames(cdt_raw_nonNA)<-paste0(colnames(cdt_raw_nonNA),"_cdt")

save(cdt_raw_nonNA,file=paste0(focal_projectFolderPath,"/data/basic/cdt_raw_nonNA.RData"))



# adoping the same scaling factor in CCL scaling 

# scaling

cdt_raw_nonNA_scaling<-t(sapply(1:dim(cdt_raw_nonNA)[1],function(i){as.numeric(cdt_raw_nonNA[i,]/embryo_duration_prcomp[i])}))
rownames(cdt_raw_nonNA_scaling)<-rownames(cdt_raw_nonNA)
colnames(cdt_raw_nonNA_scaling)<-colnames(cdt_raw_nonNA)


save(cdt_raw_nonNA_scaling,file=paste0(focal_projectFolderPath,"/data/basic/cdt_raw_nonNA_scaling.RData"))









