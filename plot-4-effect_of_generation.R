

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

focal_projectFolderPath<-r"(E:\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/figs"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"E:/papers/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_size.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mo_da1_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mo_da1_pvalue.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mo_da1_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mo_da1_pvalue.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mo_da1_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mo_da1_pvalue.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_duration_prcomp.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling_noise.RData",
	"E:/papers/canalizationENV/data/temporal/ccl_model_mother_R.RData",
	"E:/papers/canalizationENV/data/temporal/ccl_model_mother_coef.RData"
	),projectFolderPath=focal_projectFolderPath)


#################################################################################################

### generation 

mo_sis_pair_generation<-sapply(match(do.call(rbind,strsplit(names(xyz_x_model_mother_coef),"_"))[,2],lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[1,]

dataTemp<-data.frame(generation=rep(mo_sis_pair_generation,3),k=c(xyz_x_model_mother_coef,xyz_y_model_mother_coef,xyz_z_model_mother_coef),coo=c(rep("x",706),rep("y",706),rep("z",706)))

dataTemp2<-data.frame(generation=rep(5:10,3),
k=c(sapply(5:10,function(i){median(xyz_x_model_mother_coef[mo_sis_pair_generation==i])}),
	sapply(5:10,function(i){median(xyz_y_model_mother_coef[mo_sis_pair_generation==i])}),
	sapply(5:10,function(i){median(xyz_z_model_mother_coef[mo_sis_pair_generation==i])})),
coo=c(rep("x",6),rep("y",6),rep("z",6)))

##
ggplot(dataTemp,aes(x=generation,y=k))+
geom_boxplot(data=dataTemp,aes(group=generation),width=0.5,fill=NA,col="DarkGray",outlier.shape=NA,size=0.3)+
geom_jitter(height=0,width=0.1,stroke=0.1,size=0.2,col="#F8766D",alpha=0.7)+
geom_line(data=dataTemp2,aes(x=generation,y=k),lwd=0.3)+
#geom_smooth(formula=y~x,method="loess",lwd=0.2)+
geom_hline(yintercept=c(-1,0),lty=2,lwd=0.2,col="grey")+
scale_y_continuous(breaks=c(-1,0))+
scale_x_continuous(breaks=c(5:10))+
facet_wrap(nrow=3,vars(coo),strip.position="left", labeller = as_labeller(c(x="kx", y="ky", z="kz")))+
labs(x="Generation (daughter cell)",y=NULL,title="")+
theme(strip.background = element_blank(), strip.placement = "outside")

ggsave("spatial_generations_kx_ky_kz.pdf",width=unit(1.7217,"in"),height=unit(3.4,"in"))


##

#x 

ggplot(dataTemp[dataTemp$coo=="x",],aes(x=generation,y=k))+
geom_boxplot(data=dataTemp[dataTemp$coo=="x",],aes(group=generation),width=0.5,fill=NA,col="DarkGray",outlier.shape=NA,size=0.3)+
geom_jitter(height=0,width=0.15,stroke=0.2,size=0.2,col="#F8766D",alpha=0.7)+
geom_line(data=dataTemp2[dataTemp2$coo=="x",],aes(x=generation,y=k),lwd=0.3)+
#geom_smooth(formula=y~x,method="loess",lwd=0.2)+
geom_hline(yintercept=c(-1,0),lty=2,lwd=0.2,col="grey")+
scale_y_continuous(breaks=c(-1,0),limits=range(dataTemp$k))+
scale_x_continuous(breaks=c(5:10))+
labs(x="Generation (daughter cell)",y=expression(italic(k)[x]),title="")

ggsave("spatial_generations_kx.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


#y 

ggplot(dataTemp[dataTemp$coo=="y",],aes(x=generation,y=k))+
geom_boxplot(data=dataTemp[dataTemp$coo=="y",],aes(group=generation),width=0.5,fill=NA,col="DarkGray",outlier.shape=NA,size=0.3)+
geom_jitter(height=0,width=0.1,stroke=0.1,size=0.2,col="#F8766D",alpha=0.7)+
geom_line(data=dataTemp2[dataTemp2$coo=="y",],aes(x=generation,y=k),lwd=0.3)+
#geom_smooth(formula=y~x,method="loess",lwd=0.2)+
geom_hline(yintercept=c(-1,0),lty=2,lwd=0.2,col="grey")+
scale_y_continuous(breaks=c(-1,0),limits=range(dataTemp$k))+
scale_x_continuous(breaks=c(5:10))+
labs(x="Generation (daughter cell)",y=expression(italic(k)[y]),title="")

ggsave("spatial_generations_ky.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


#z  

ggplot(dataTemp[dataTemp$coo=="z",],aes(x=generation,y=k))+
geom_boxplot(data=dataTemp[dataTemp$coo=="z",],aes(group=generation),width=0.5,fill=NA,col="DarkGray",outlier.shape=NA,size=0.3)+
geom_jitter(height=0,width=0.1,stroke=0.1,size=0.2,col="#F8766D",alpha=0.7)+
geom_line(data=dataTemp2[dataTemp2$coo=="z",],aes(x=generation,y=k),lwd=0.3)+
#geom_smooth(formula=y~x,method="loess",lwd=0.2)+
geom_hline(yintercept=c(-1,0),lty=2,lwd=0.2,col="grey")+
scale_y_continuous(breaks=c(-1,0),limits=range(dataTemp$k))+
scale_x_continuous(breaks=c(5:10))+
labs(x="Generation (daughter cell)",y=expression(italic(k)[z]),title="")

ggsave("spatial_generations_kz.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

#### CCL 


mo_sis_pair_generation<-sapply(match(do.call(rbind,strsplit(names(ccl_model_mother_coef),"_"))[,2],lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[1,]

dataTemp<-data.frame(generation=mo_sis_pair_generation,k=ccl_model_mother_coef)

dataTemp2<-data.frame(generation=c(5:9),
k=sapply(5:9,function(i){median(ccl_model_mother_coef[mo_sis_pair_generation==i])}))

ggplot(dataTemp,aes(x=generation,y=k))+
geom_boxplot(data=dataTemp,aes(group=generation),width=0.5,fill=NA,col="DarkGray",outlier.shape=NA,size=0.3)+
geom_jitter(height=0,width=0.1,stroke=0.1,size=0.2,col="#F8766D",alpha=0.7)+
geom_line(data=dataTemp2,aes(x=generation,y=k),lwd=0.3)+
geom_hline(yintercept=c(-1,0),lty=2,lwd=0.2,col="grey")+
scale_y_continuous(breaks=c(-1,0),limits=range(c(0,dataTemp$k)))+
scale_x_continuous(breaks=c(5:9))+
labs(x="Generation (daughter cell)",y=expression(italic(k)[CCL]),title="")

ggsave("temporal_generations_kCCL.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))









