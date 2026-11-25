


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

#

identy_score2<-function(x,y){
	return(sum((x+y)^2)/2/sum(x^2+y^2))
}


#######################################################################

#######################################################################

library(ggplot2)
library(reshape)
library(gplots)
library(scatterplot3d)

library(extrafont)

#font_import()
loadfonts()


focal_projectFolderPath<-r"(E:\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/figs"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"E:/papers/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_duration_prcomp.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_duration_prcomp_wtMean.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/temporal/ccl_model_mother_coef.RData"
	),projectFolderPath=focal_projectFolderPath)


#


feature_names<-unlist(strsplit(colnames(ccl_raw_nonNA_scaling),"_ccl"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]


##

X_mo<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$mo,feature_names)]
X_da1<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$da1,feature_names)]
X_da2<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo



###########################################################################################


### alignment process of ctr_emb17 

which(rownames(ccl_raw_nonNA_scaling)=="HMP-2_emb6")

dataTemp<-data.frame(ref=rep(c(ccl_raw_nonNA[58,]),2),emp_example=c(c(ccl_raw_nonNA[832,]),c(ccl_raw_nonNA_scaling[832,])),type=c(rep("raw",353),rep("scaled",353)))

#

ggplot(dataTemp[dataTemp$type=="raw",],aes(x=ref,y=emp_example))+
geom_point(size=0.1,stroke=0.1,col=colors()[182])+
geom_abline(slope=1,intercept=0,lty=2,lwd=0.3)+
geom_abline(slope=embryo_duration_prcomp[832],intercept=0,col="red",lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,1:2]),breaks=c(20,60,100))+
scale_y_continuous(limits=range(dataTemp[,1:2]),breaks=c(20,60,100))+
labs(x="CCL of ctr_emb58",y="CCL of HMP-2_emb6",title="before scaling")

ggsave("temporal_HMP-2_emb6_scaling_factor_before.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

#

ggplot(dataTemp[dataTemp$type=="scaled",],aes(x=ref,y=emp_example))+
geom_point(size=0.1,stroke=0.1,col=colors()[182])+
geom_abline(slope=1,intercept=0,lty=2,lwd=0.3)+
geom_abline(slope=embryo_duration_prcomp[832],intercept=0,col="red",lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,1:2]),breaks=c(20,60,100))+
scale_y_continuous(limits=range(dataTemp[,1:2]),breaks=c(20,60,100))+
labs(x="CCL of ctr_emb58",y="CCL of HMP-2_emb6",title="after scaling")

ggsave("temporal_HMP-2_emb6_scaling_factor_after.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

	
### quality evaluation

# distribution of IS  
 

IS_CCL<-sapply(1:2039,function(i){
	return(identy_score2(c(ccl_raw_nonNA_scaling[i,]),c(ccl_raw_nonNA_scaling[58,])))
})



dataTemp<-data.frame(IS=IS_CCL)
	
ggplot(dataTemp,aes(IS))+
geom_density()+
scale_x_continuous(breaks=c(0.985,1))+
scale_y_continuous(breaks=c(0,3000))+
labs(x="IS",y="Density",title="alignment quality")


ggsave("temporal_alignment_quality_IS_distribution.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



### robustness of scaling factor 

# compare between scaling factor obtained based on wt58 and wtMean

dataTemp<-data.frame(SF=embryo_duration_prcomp,SF_wtMean=embryo_duration_prcomp_wtMean)

ggplot(dataTemp,aes(x=SF,y=SF_wtMean))+
geom_point(size=0.7)+
geom_smooth(formula=y~0+x,method="lm",lwd=0.7)+
labs(x="SF estimated by ctr_emb58",y="SF estimated by average of control embryos",title="robustness of SF")


ggsave("temporal_alignment_quality_robustness_of_scaling_factor.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


##################################################################3


## all with good quality

poorInd<-which(IS_CCL<0.9)












