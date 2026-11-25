




### load Rscripts and RData 

loadResources<-function(resourceVector,projectFolderPath=focal_projectFolderPath){
	
	library(tools)
	
	projectFolderPath<-gsub("\\\\", "/", projectFolderPath)
	projectFolder<-basename(projectFolderPath)
	
	resourceVector<-sapply(resourceVector,function(rs){
		rs_folders<-unlist(strsplit(rs,"/"))
		return(gsub("//","/",paste0(projectFolderPath,"/",paste0(rs_folders[(which(rs_folders==projectFolder)+1):length(rs_folders)],collapse="/"))))
	})
	
	sapply(1:length(resourceVector),function(i){
		if(file_ext(resourceVector[i])=="R"){
			source(resourceVector[i],chdir=T)
			return(resourceVector[i])
		}
		if(file_ext(resourceVector[i])=="RData"){
			load(resourceVector[i],verbose=T,envir=globalenv())
			return(resourceVector[i])
		}

	})


}



#####################################################################################

#####################################################################################


library(ggplot2)
library(glmnet)
library(yardstick)


focal_projectFolderPath<-r"(E:\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/figs"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/embryo_lethal.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/cvglm_wormbase.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/cvglm_wormbase_withoutCCL.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/wormbase_pheno_mat_select.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/residual_feature_comb.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/onto_phe.RData"
	),projectFolderPath=focal_projectFolderPath)



###

embryo_lethal<-as.numeric(embryo_lethal)

###############################################################################

### with CCL 

# ROC curve and AUC by logistic regression under machine learning framework 

auc_summary<-data.frame(phenoID=colnames(wormbase_pheno_mat_select),
	AUC_cvm=sapply(cvglm_wormbase,function(x){x$cvm[which(x$lambda==x$lambda.min)]}),
	AUC_cvsd=sapply(cvglm_wormbase,function(x){x$cvsd[which(x$lambda==x$lambda.min)]}),
	phenoAnnot=onto_phe[match(colnames(wormbase_pheno_mat_select),onto_phe[,1]),2]
)

save(auc_summary,file=paste0(focal_projectFolderPath,"/data/prediction/auc_summary.RData"))

roc_list<-lapply(1:length(cvglm_wormbase),function(i){
	focalPheno<-wormbase_pheno_mat_select[which(!is.na(wormbase_pheno_mat_select[,i])),i]
	return(roc.glmnet(cvglm_wormbase[[i]]$fit.preval, newy = focalPheno,family="binomial")[[which(cvglm_wormbase[[i]]$lambda==cvglm_wormbase[[i]]$lambda.min)]])
})

names(roc_list)<-names(cvglm_wormbase)


dataTemp<-data.frame(do.call(rbind,roc_list),
	phenoAnnot=factor(unlist(lapply(1:20,function(i){rep(auc_summary$phenoAnnot[i],dim(roc_list[[i]])[1])})),level=auc_summary$phenoAnnot[order(auc_summary$AUC_cvm,decreasing=T)]))

ggplot(dataTemp,aes(x=FPR,y=TPR))+
geom_line()+
geom_abline(slope=1,intercept=0,lty=2)+
facet_wrap(vars(phenoAnnot),nrow=5)+
scale_x_continuous(breaks=c(0,0.5,1))+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(title="20 phenotypes curated from Wormbase")+
theme(strip.background=element_blank())


ggsave("Prediction_wormbase_ROCs.pdf",width=unit(6.4,"in"),height=unit(8,"in"))

paste("AUC=",paste(round(auc_summary$AUC_cvm[order(auc_summary$AUC_cvm,decreasing=T)],2),round(auc_summary$AUC_cvsd[order(auc_summary$AUC_cvm,decreasing=T)],2),sep="±"),sep="")

paste("n=",sapply(1:20,function(i){length(which(!is.na(wormbase_pheno_mat_select[,order(auc_summary$AUC_cvm,decreasing=T)[i]])))}),sep="")


### slow growth and dumpy

dataTemp2<-dataTemp[dataTemp$phenoAnnot=="slow growth",]

ggplot(dataTemp2,aes(x=FPR,y=TPR))+
geom_line()+
geom_abline(slope=1,intercept=0,lty=2)+
scale_x_continuous(breaks=c(0,0.5,1))+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(title="slow growth")

ggsave("Prediction_wormbase_ROC_slow_growth.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


dataTemp3<-dataTemp[dataTemp$phenoAnnot=="dumpy",]

ggplot(dataTemp3,aes(x=FPR,y=TPR))+
geom_line()+
geom_abline(slope=1,intercept=0,lty=2)+
scale_x_continuous(breaks=c(0,0.5,1))+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(title="dumpy")

ggsave("Prediction_wormbase_ROC_dumpy.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



### without CCL 

# ROC curve and AUC by logistic regression under machine learning framework 

auc_summary_withoutCCL<-data.frame(phenoID=colnames(wormbase_pheno_mat_select),
	AUC_cvm=sapply(cvglm_wormbase_withoutCCL,function(x){x$cvm[which(x$lambda==x$lambda.min)]}),
	AUC_cvsd=sapply(cvglm_wormbase_withoutCCL,function(x){x$cvsd[which(x$lambda==x$lambda.min)]}),
	phenoAnnot=onto_phe[match(colnames(wormbase_pheno_mat_select),onto_phe[,1]),2]
)

save(auc_summary_withoutCCL,file=paste0(focal_projectFolderPath,"/data/prediction/auc_summary_withoutCCL.RData"))

roc_list_withoutCCL<-lapply(1:length(cvglm_wormbase_withoutCCL),function(i){
	focalPheno<-wormbase_pheno_mat_select[which(!is.na(wormbase_pheno_mat_select[,i])),i]
	return(roc.glmnet(cvglm_wormbase_withoutCCL[[i]]$fit.preval, newy = focalPheno,family="binomial")[[which(cvglm_wormbase_withoutCCL[[i]]$lambda==cvglm_wormbase_withoutCCL[[i]]$lambda.min)]])
})

names(roc_list_withoutCCL)<-names(cvglm_wormbase_withoutCCL)


dataTemp<-data.frame(do.call(rbind,roc_list_withoutCCL),
	phenoAnnot=factor(unlist(lapply(1:20,function(i){rep(auc_summary_withoutCCL$phenoAnnot[i],dim(roc_list_withoutCCL[[i]])[1])})),level=auc_summary_withoutCCL$phenoAnnot[order(auc_summary_withoutCCL$AUC_cvm,decreasing=T)]))

ggplot(dataTemp,aes(x=FPR,y=TPR))+
geom_line()+
geom_abline(slope=1,intercept=0,lty=2)+
facet_wrap(vars(phenoAnnot),nrow=5)+
scale_x_continuous(breaks=c(0,0.5,1))+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(title="20 phenotypes curated from Wormbase")+
theme(strip.background=element_blank())


ggsave("Prediction_wormbase_ROCs_withoutCCL.pdf",width=unit(6.4,"in"),height=unit(8,"in"))

paste("AUC=",paste(round(auc_summary_withoutCCL$AUC_cvm[order(auc_summary_withoutCCL$AUC_cvm,decreasing=T)],2),round(auc_summary_withoutCCL$AUC_cvsd[order(auc_summary_withoutCCL$AUC_cvm,decreasing=T)],2),sep="±"),sep="")

paste("n=",sapply(1:20,function(i){length(which(!is.na(wormbase_pheno_mat_select[,order(auc_summary_withoutCCL$AUC_cvm,decreasing=T)[i]])))}),sep="")








