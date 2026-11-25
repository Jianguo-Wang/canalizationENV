


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


#######################################################################

#######################################################################


library(mnormt)
library(ggplot2)



focal_projectFolderPath<-r"(F:\wangj\Desktop\temp\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/data/canalGene"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_duration_prcomp.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling_noise.RData"
	),projectFolderPath=focal_projectFolderPath)



	
####################################################################################

###CCL 

feature_names<-unlist(strsplit(colnames(ccl_raw_nonNA_scaling_noise),"_ccl"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]

#

X_mo<-ccl_raw_nonNA_scaling_noise[1:1928,match(mo_sis_pair$mo,feature_names)]
X_da1<-ccl_raw_nonNA_scaling_noise[1:1928,match(mo_sis_pair$da1,feature_names)]
X_da2<-ccl_raw_nonNA_scaling_noise[1:1928,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo


### mutants with deviation from population mean of WT


embryo_geneName<-do.call(rbind,strsplit(rownames(delta_X),"_"))[,1]

#

ccl_raw_nonNA_scaling_mean_WT<-apply(ccl_raw_nonNA_scaling[1:105,],2,mean)

dist_with_WT_mean<-sqrt(sapply(1:1928,function(i){sum((ccl_raw_nonNA_scaling_mean_WT-ccl_raw_nonNA_scaling[i,])^2)}))


#1-pnorm(max(dist_with_WT_mean[1:105]),mean(dist_with_WT_mean[1:105]),sd(dist_with_WT_mean[1:105]))=4.032075e-07<0.05/(1928-105)
#FDR=length(which(dist_with_WT_mean[1:105]>cutoff_1))/105<1%

cutoff_1<-qnorm(1-0.05/(1928-105),mean(dist_with_WT_mean[1:105]),sd(dist_with_WT_mean[1:105]))


#

sigEmbryoInd_dist_with_WT_mean<-which(dist_with_WT_mean>cutoff_1)

sigEmbryoInd_dist_with_WT_mean_summary=data.frame(
	embryoInd=sigEmbryoInd_dist_with_WT_mean,
	embryoID=rownames(delta_X)[sigEmbryoInd_dist_with_WT_mean],
	geneName=embryo_geneName[sigEmbryoInd_dist_with_WT_mean],
	totalEmbryoNum=as.numeric(table(embryo_geneName)[match(embryo_geneName[sigEmbryoInd_dist_with_WT_mean],names(table(embryo_geneName)))]),
	sigEmbryoNum=as.numeric(table(embryo_geneName[sigEmbryoInd_dist_with_WT_mean])[match(embryo_geneName[sigEmbryoInd_dist_with_WT_mean],names(table(embryo_geneName[sigEmbryoInd_dist_with_WT_mean])))]),
	dist_with_WT_mean=dist_with_WT_mean[sigEmbryoInd_dist_with_WT_mean]
)

sigEmbryoInd_dist_with_WT_mean_summary$ratio<-round(sigEmbryoInd_dist_with_WT_mean_summary$sigEmbryoNum/sigEmbryoInd_dist_with_WT_mean_summary$totalEmbryoNum,2)

sigEmbryoInd_dist_with_WT_mean_summary_replicate<-sigEmbryoInd_dist_with_WT_mean_summary[sigEmbryoInd_dist_with_WT_mean_summary$ratio>0.5,]
sigEmbryoInd_dist_with_WT_mean_summary_replicate<-sigEmbryoInd_dist_with_WT_mean_summary_replicate[order(sigEmbryoInd_dist_with_WT_mean_summary_replicate$ratio,sigEmbryoInd_dist_with_WT_mean_summary_replicate$totalEmbryoNum,decreasing=T),]

sigGeneName_dist_with_WT_mean<-unique(sigEmbryoInd_dist_with_WT_mean_summary_replicate$geneName)


##

dataTemp<-data.frame(Dist=dist_with_WT_mean,type=c(rep("WT",105),rep("MT",1928-105)))

ggplot(dataTemp,aes(x=Dist,col=type,fill=type))+geom_density(alpha=0.5,size=0.3)+
geom_vline(xintercept=cutoff_1,lty=2,lwd=0.5)+
labs(x="Distance from WT average",y="Density",title="")

ggsave("temporal_mutant_distance_from_average.pdf",width=unit(2.4,"in"),height=unit(2,"in"))



### mother model -1

delta_X_pred_mo<-sapply(1:dim(delta_X)[2],function(i){
	return(cbind(rep(1,dim(delta_X)[1]),X_mo[,i])%*%
		summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[,1])
})


Res2_mo_1<-sapply(1:1928,function(i){sum((delta_X[i,]-delta_X_pred_mo[i,])^2)/sum(delta_X[i,]^2)})
sigEmbryoInd_dist_with_WT_mean_summary_replicate$Res2_mo_1<-round(Res2_mo_1[sigEmbryoInd_dist_with_WT_mean_summary_replicate$embryoInd],2)

#FDR=length(which(Res2_mo_1[1:105]>cutoff_2))/105<1%

cutoff_2<-qnorm(1-0.05/length(sigEmbryoInd_dist_with_WT_mean_summary_replicate$embryoInd),mean(Res2_mo_1[1:105]),sd(Res2_mo_1[1:105]))


# mutant with mother negative feedback defect

sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect<-sigEmbryoInd_dist_with_WT_mean_summary_replicate[sigEmbryoInd_dist_with_WT_mean_summary_replicate$Res2_mo_1>cutoff_2,]


##

dataTemp<-data.frame(Res2=Res2_mo_1[c(1:105,sigEmbryoInd_dist_with_WT_mean_summary_replicate$embryoInd)],type=c(rep("WT",105),rep("MT",length(sigEmbryoInd_dist_with_WT_mean_summary_replicate$embryoInd))))

ggplot(dataTemp,aes(x=Res2,col=type,fill=type))+geom_density(alpha=0.5,size=0.3)+
geom_vline(xintercept=cutoff_2,lty=2,lwd=0.5)+
labs(x="SSE/SST",y="Density",title=expression(paste(italic(Delta),italic(L),"~",italic(L)[M])))

ggsave("temporal_mutant_SSE_SST_mother_defect.pdf",width=unit(2.4,"in"),height=unit(2,"in"))


sigGeneName_mother_defect<-unique(sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect$geneName)

write.csv(sigGeneName_mother_defect,"sigGeneName_temporal.csv")

### temporal_mutant_res

temporal_mutant_res<-data.frame(
	EmbryoID=sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect$embryoID,
	Gene_name=sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect$geneName,
	Distance_from_average=sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect$dist_with_WT_mean,
	Distance_from_average_pvalue=1-pnorm(sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect$dist_with_WT_mean,mean(dist_with_WT_mean[1:105]),sd(dist_with_WT_mean[1:105])),
	SSE_SST=sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect$Res2_mo_1,
	SSE_SST_pvalue=1-pnorm(sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect$Res2_mo_1,mean(Res2_mo_1[1:105]),sd(Res2_mo_1[1:105]))
)


write.csv(temporal_mutant_res,"Tab.S1-temporal_mutant_res.csv",row.names=FALSE)


############################################################################################
































