




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

library(ggplot2)
library(reshape)
library(gplots)

focal_projectFolderPath<-r"(F:\wangj\Desktop\temp\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/figs"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"E:/papers/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling_noise.RData"
	),projectFolderPath=focal_projectFolderPath)


# xyz

wt_mean<-apply(q1_cell350_alignment_xyz[1:105,],2,mean)
mt_mean<-apply(q1_cell350_alignment_xyz[106:1928,],2,mean)
env_mean<-apply(q1_cell350_alignment_xyz[1929:2039,],2,mean)

wt_sd<-apply(q1_cell350_alignment_xyz[1:105,],2,sd)
mt_sd<-apply(q1_cell350_alignment_xyz[106:1928,],2,sd)
env_sd<-apply(q1_cell350_alignment_xyz[1929:2039,],2,sd)

### CCL 

wt_mean_CCL<-apply(ccl_raw_nonNA_scaling[1:105,],2,mean)
mt_mean_CCL<-apply(ccl_raw_nonNA_scaling[106:1928,],2,mean)
env_mean_CCL<-apply(ccl_raw_nonNA_scaling[1929:2039,],2,mean)

wt_sd_CCL<-apply(ccl_raw_nonNA_scaling[1:105,],2,sd)
mt_sd_CCL<-apply(ccl_raw_nonNA_scaling[106:1928,],2,sd)
env_sd_CCL<-apply(ccl_raw_nonNA_scaling[1929:2039,],2,sd)


### xyz: mean : wt-env

dataTemp<-data.frame(Control=wt_mean,Perturbation=env_mean)
ggplot(dataTemp,aes(x=Control,y=Perturbation))+
geom_point(size=0.5,stroke=0.2,col=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp),breaks=c(-100,0,100))+
scale_y_continuous(limits=range(dataTemp),breaks=c(-100,0,100))+
labs(x="Control embryo",y="Environment embryo",title="Mean of xyz")

ggsave("common_expectation_xyz_mean_wt_Vs_env.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### xyz: mean : wt-RNAi

dataTemp<-data.frame(Control=wt_mean,Perturbation=mt_mean)
ggplot(dataTemp,aes(x=Control,y=Perturbation))+
geom_point(size=0.5,stroke=0.2,col=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp),breaks=c(-100,0,100))+
scale_y_continuous(limits=range(dataTemp),breaks=c(-100,0,100))+
labs(x="Control embryo",y="RNAi embryo",title="Mean of xyz")

ggsave("common_expectation_xyz_mean_wt_Vs_RNAi.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

### xyz: sd : wt-env

dataTemp<-data.frame(Control=wt_sd,Perturbation=env_sd)
ggplot(dataTemp,aes(x=Control,y=Perturbation))+
geom_point(size=0.5,stroke=0.2,col=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp),breaks=c(5,10,15))+
scale_y_continuous(limits=range(dataTemp),breaks=c(5,10,15))+
labs(x="Control embryo",y="Environment embryo",title="SD of xyz")

ggsave("common_expectation_xyz_sd_wt_Vs_env.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

### xyz: sd : wt-RNAi

dataTemp<-data.frame(Control=wt_sd,Perturbation=mt_sd)
ggplot(dataTemp,aes(x=Control,y=Perturbation))+
geom_point(size=0.5,stroke=0.2,col=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp),breaks=c(5,10,15))+
scale_y_continuous(limits=range(dataTemp),breaks=c(5,10,15))+
labs(x="Control embryo",y="RNAi embryo",title="SD of xyz")

ggsave("common_expectation_xyz_sd_wt_Vs_RNAi.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### CCL: mean : wt-env

dataTemp<-data.frame(Control=wt_mean_CCL,Perturbation=env_mean_CCL)
ggplot(dataTemp,aes(x=Control,y=Perturbation))+
geom_point(size=0.5,stroke=0.2,col=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp),breaks=c(20,50,80))+
scale_y_continuous(limits=range(dataTemp),breaks=c(20,50,80))+
labs(x="Control embryo",y="Environment embryo",title="Mean of CCL")

ggsave("common_expectation_CCL_mean_wt_Vs_env.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

### CCL: mean : wt-RNAi

dataTemp<-data.frame(Control=wt_mean_CCL,Perturbation=mt_mean_CCL)
ggplot(dataTemp,aes(x=Control,y=Perturbation))+
geom_point(size=0.5,stroke=0.2,col=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp),breaks=c(20,50,80))+
scale_y_continuous(limits=range(dataTemp),breaks=c(20,50,80))+
labs(x="Control embryo",y="RNAi embryo",title="Mean of CCL")

ggsave("common_expectation_CCL_mean_wt_Vs_RNAi.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

### CCL: sd : wt-env

dataTemp<-data.frame(Control=wt_sd_CCL,Perturbation=env_sd_CCL)
ggplot(dataTemp,aes(x=Control,y=Perturbation))+
geom_point(size=0.5,stroke=0.2,col=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp),breaks=c(1,3,5))+
scale_y_continuous(limits=range(dataTemp),breaks=c(1,3,5))+
labs(x="Control embryo",y="Environment embryo",title="SD of CCL")

ggsave("common_expectation_CCL_sd_wt_Vs_env.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

### CCL: sd : wt-RNAi

dataTemp<-data.frame(Control=wt_sd_CCL,Perturbation=mt_sd_CCL)
ggplot(dataTemp,aes(x=Control,y=Perturbation))+
geom_point(size=0.5,stroke=0.2,col=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp),breaks=c(1,3,5))+
scale_y_continuous(limits=range(dataTemp),breaks=c(1,3,5))+
labs(x="Control embryo",y="RNAi embryo",title="SD of CCL")

ggsave("common_expectation_CCL_sd_wt_Vs_RNAi.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))





par(mfrow=c(2,4))
plot(wt_mean,mt_mean)
abline(0,1)
plot(wt_mean,env_mean)
abline(0,1)
plot(wt_sd,mt_sd)
abline(0,1)
plot(wt_sd,env_sd)
abline(0,1)

plot(wt_mean_CCL,mt_mean_CCL)
abline(0,1)
plot(wt_mean_CCL,env_mean_CCL)
abline(0,1)
plot(wt_sd_CCL,mt_sd_CCL)
abline(0,1)
plot(wt_sd_CCL,env_sd_CCL)
abline(0,1)


PC12<-princomp(q1_cell350_alignment_xyz[,1:714])$scores[,1:2]
col1<-c(rep("grey",105),rep("#F8766d",1823),rep("#00BFC4",111))

plot(PC12[,1],PC12[,2],col=col1)

dataTemp<-data.frame(PC1=PC12[,1],PC2=PC12[,2],type=c(rep("wt",105),rep("mt",1823),rep("env",111)))

ggplot(dataTemp,aes(x=PC1,y=PC2,col=type))+
geom_point(size=0.5,alpha=0.5)+
geom_density2d()


dataTemp<-data.frame(PC1=PC12[,1],PC2=PC12[,2],type=c(rep("wt",105),rep("mt",1823),rep("env",111)))
ggplot(dataTemp,aes(x=PC1,y=PC2))+
geom_point(size=0.5,alpha=0.5)+
geom_density2d()+
facet_wrap(vars(type))

ggplot(dataTemp,aes(col=type,x=PC1))+geom_density(alpha=0.75)
ggplot(dataTemp,aes(col=type,x=PC2))+geom_density(alpha=0.75)






