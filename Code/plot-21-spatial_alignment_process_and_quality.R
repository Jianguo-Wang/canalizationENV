





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
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_size.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_size_wtMean.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/interpolation/coo_list_q1_cell350.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/coo_list_q1_cell350_translation.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/coo_list_q1_cell350_translation_rotation.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/coo_list_q1_cell350_translation_rotation_scaling.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/spatial/xyz_x_model_mother_coef.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/spatial/xyz_y_model_mother_coef.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/spatial/xyz_z_model_mother_coef.RData"
	),projectFolderPath=focal_projectFolderPath)


#


feature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise)[1:714],"_x"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]


###########################################################################################


### alignment process of ctr_emb17 

color1<-rep(colors()[192],714)

# raw 

pdf("spatial_ctr_emb17_raw.pdf", width=3.2, height=3.4)

scatterplot3d(coo_list_q1_cell350[[17]],xlab="x",ylab="y",zlab="z",cex.symbols=0.1,color=color1,pch=1,xlim=c(0,320),ylim=c(0,320),zlim=c(0,320),main="raw coordinates")

dev.off()

# translation 

pdf("spatial_ctr_emb17_raw_translation.pdf", width=3.2, height=3.4)

scatterplot3d(coo_list_q1_cell350_translation[[17]],xlab="x",ylab="y",zlab="z",cex.symbols=0.1,color=color1,pch=1,xlim=c(-160,160),ylim=c(-160,160),zlim=c(-160,160),main="translation")

dev.off()

# rotation

pdf("spatial_ctr_emb17_raw_translation_rotation.pdf", width=3.2, height=3.4)

scatterplot3d(coo_list_q1_cell350_translation_rotation[[17]],xlab="x",ylab="y",zlab="z",cex.symbols=0.1,color=color1,pch=1,xlim=c(-160,160),ylim=c(-160,160),zlim=c(-160,160),main="rotation")

dev.off()

#scaling 

pdf("spatial_ctr_emb17_raw_translation_rotation_scaling.pdf", width=3.2, height=3.4)

scatterplot3d(coo_list_q1_cell350_translation_rotation_scaling[[17]],xlab="x",ylab="y",zlab="z",cex.symbols=0.1,color=color1,pch=1,xlim=c(-160,160),ylim=c(-160,160),zlim=c(-160,160),main="scaling")

dev.off()


### scaling example 

# before

dataTemp<-data.frame(x_ref=coo_list_q1_cell350_translation_rotation[[58]][,1],
	x_example=coo_list_q1_cell350_translation_rotation[[17]][,1])

ggplot(dataTemp,aes(x=x_ref,y=x_example))+
geom_point(size=0.1,stroke=0.1,col=colors()[182])+
geom_abline(slope=1,intercept=0,lty=2,lwd=0.3)+
geom_abline(slope=embryo_size[17,1],intercept=0,col="red",lwd=0.3)+
scale_x_continuous(limits=range(dataTemp),breaks=c(-100,0,100))+
scale_y_continuous(limits=range(dataTemp),breaks=c(-100,0,100))+
labs(x="ctr_emb58",y="ctr_emb17",title="before scaling")


ggsave("spatial_ctr_emb17_scaling_factor_before.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

# after

dataTemp<-data.frame(x_ref=coo_list_q1_cell350_translation_rotation_scaling[[58]][,1],
	x_example=coo_list_q1_cell350_translation_rotation_scaling[[17]][,1])

ggplot(dataTemp,aes(x=x_ref,y=x_example))+
geom_point(size=0.1,stroke=0.1,col=colors()[182])+
geom_abline(slope=1,intercept=0,lty=2,lwd=0.3)+
geom_abline(slope=embryo_size[17,1],intercept=0,col="red",lwd=0.3)+
scale_x_continuous(limits=range(dataTemp),breaks=c(-100,0,100))+
scale_y_continuous(limits=range(dataTemp),breaks=c(-100,0,100))+
labs(x="ctr_emb58",y="ctr_emb17",title="after scaling")


ggsave("spatial_ctr_emb17_scaling_factor_after.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


	
### quality evaluation

# definition

dataTemp<-data.frame(x_ref=coo_list_q1_cell350_translation_rotation_scaling[[58]][,1],
	x_example=coo_list_q1_cell350_translation_rotation_scaling[[17]][,1])

ggplot(dataTemp,aes(x=x_ref,y=x_example))+
geom_point(size=0.1,stroke=0.1,col=colors()[182])+
geom_abline(slope=1,intercept=0,lty=2,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp),breaks=c(-100,0,100))+
scale_y_continuous(limits=range(dataTemp),breaks=c(-100,0,100))+
labs(x="ctr_emb58",y="ctr_emb17",title="definition of IS")


ggsave("spatial_ctr_emb17_IS_definition.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

# distribution of IS  
 

IS_x<-sapply(1:2039,function(i){
	return(identy_score2(coo_list_q1_cell350_translation_rotation_scaling[[i]][,1],coo_list_q1_cell350_translation_rotation_scaling[[58]][,1]))
})

IS_y<-sapply(1:2039,function(i){
	return(identy_score2(coo_list_q1_cell350_translation_rotation_scaling[[i]][,2],coo_list_q1_cell350_translation_rotation_scaling[[58]][,2]))
})

IS_z<-sapply(1:2039,function(i){
	return(identy_score2(coo_list_q1_cell350_translation_rotation_scaling[[i]][,3],coo_list_q1_cell350_translation_rotation_scaling[[58]][,3]))
})



dataTemp<-data.frame(IS=c(IS_x,IS_y,IS_z),
	coo=c(rep("x",length(IS_x)),rep("y",length(IS_y)),rep("z",length(IS_z))))
	
ggplot(dataTemp,aes(IS,col=coo))+
geom_density(size=0.2)+
scale_x_continuous(breaks=c(0.6,0.8,1))+
scale_y_continuous(breaks=c(0,150,300))+
theme(legend.position=c(0.2,0.8))+
labs(x="IS",y="Density",title="alignment quality")


ggsave("spatial_alignment_quality_IS_distribution.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



### discrepance between coordinates 


R_x_x<-sapply(1:2039,function(i){cor(coo_list_q1_cell350_translation_rotation_scaling[[58]][,1],coo_list_q1_cell350_translation_rotation_scaling[[i]][,1])})
R_x_y<-sapply(1:2039,function(i){cor(coo_list_q1_cell350_translation_rotation_scaling[[58]][,1],coo_list_q1_cell350_translation_rotation_scaling[[i]][,2])})
R_x_z<-sapply(1:2039,function(i){cor(coo_list_q1_cell350_translation_rotation_scaling[[58]][,1],coo_list_q1_cell350_translation_rotation_scaling[[i]][,3])})

R_y_x<-sapply(1:2039,function(i){cor(coo_list_q1_cell350_translation_rotation_scaling[[58]][,2],coo_list_q1_cell350_translation_rotation_scaling[[i]][,1])})
R_y_y<-sapply(1:2039,function(i){cor(coo_list_q1_cell350_translation_rotation_scaling[[58]][,2],coo_list_q1_cell350_translation_rotation_scaling[[i]][,2])})
R_y_z<-sapply(1:2039,function(i){cor(coo_list_q1_cell350_translation_rotation_scaling[[58]][,2],coo_list_q1_cell350_translation_rotation_scaling[[i]][,3])})

R_z_x<-sapply(1:2039,function(i){cor(coo_list_q1_cell350_translation_rotation_scaling[[58]][,3],coo_list_q1_cell350_translation_rotation_scaling[[i]][,1])})
R_z_y<-sapply(1:2039,function(i){cor(coo_list_q1_cell350_translation_rotation_scaling[[58]][,3],coo_list_q1_cell350_translation_rotation_scaling[[i]][,2])})
R_z_z<-sapply(1:2039,function(i){cor(coo_list_q1_cell350_translation_rotation_scaling[[58]][,3],coo_list_q1_cell350_translation_rotation_scaling[[i]][,3])})


#

dataTemp<-data.frame(type=as.character(rep(1,2039*9)),
	R=c(R_x_x,R_x_y,R_x_z,
	R_y_x,R_y_y,R_y_z,
	R_z_x,R_z_y,R_z_z),
	coo_ref=c(rep("reference: x",2039*3),rep("reference: y",2039*3),rep("reference: z",2039*3)),
	coo_focal=rep(c(rep("x",2039),rep("y",2039),rep("z",2039)),3))


ggplot(dataTemp,aes(x=coo_focal,y=R))+
geom_boxplot(outlier.size=0.1,outlier.stroke=0.1,width=0.5,size=0.3)+
geom_hline(yintercept=0,lty=2,lwd=0.3)+
facet_wrap(vars(coo_ref),strip.position="right",nrow=3)+
labs(x="coordinates of embryos",y=expression(paste("Pearson's ",italic(R))),title="discrepancy between coordinates after alignment")+
theme(panel.border=element_rect(fill=NA))


ggsave("spatial_alignment_quality_discrepancy_between_coordinates.pdf",width=unit(3.2,"in"),height=unit(3.4,"in"))


### robustness of scaling factor 
# compare between scaling factor obtained based on wt58 and wtMean

dataTemp<-data.frame(SF=c(embryo_size),SF_wtMean=c(embryo_size_wtMean),coo=c(rep("x",2039),rep("y",2039),rep("z",2039)))

ggplot(dataTemp,aes(x=SF,y=SF_wtMean))+
geom_point(size=0.7)+
geom_smooth(formula=y~0+x,method="lm",lwd=0.7)+
facet_wrap(vars(coo),nrow=2,scales="free")+
labs(x="SF estimated by ctr_emb58",y="SF estimated by average of control embryos",title="robustness of SF")


ggsave("spatial_alignment_quality_robustness_of_scaling_factor.pdf",width=unit(3.2,"in"),height=unit(3.4,"in"))


###############################################################333

###

q1_cell350_alignment_xyz_noise_x<-q1_cell350_alignment_xyz_noise[,1:714]
q1_cell350_alignment_xyz_noise_y<-q1_cell350_alignment_xyz_noise[,(1+714):(714*2)]
q1_cell350_alignment_xyz_noise_z<-q1_cell350_alignment_xyz_noise[,(1+714*2):(714*3)]

#

feature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise[,1:714]),"_x"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]

##

poorInd<-which(IS_x<0.9|IS_y<0.9|IS_z<0.9)

save(poorInd,file=paste0(focal_projectFolderPath,"/data/alignment/poorInd.RData"))

## kx reestimated by removing poor-quality embryos 

X_mo<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo

xyz_x_model_mother_coef_quality<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[-poorInd,i]~X_mo[-poorInd,i]))$coefficients[2,1]})


## ky reestimated by removing poor-quality embryos 


X_mo<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo

xyz_y_model_mother_coef_quality<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[-poorInd,i]~X_mo[-poorInd,i]))$coefficients[2,1]})


## kz reestimated by removing poor-quality embryos 


X_mo<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo

xyz_z_model_mother_coef_quality<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[-poorInd,i]~X_mo[-poorInd,i]))$coefficients[2,1]})


###

dataTemp<-data.frame(k=c(xyz_x_model_mother_coef,xyz_y_model_mother_coef,xyz_z_model_mother_coef),
	k_q=c(xyz_x_model_mother_coef_quality,xyz_y_model_mother_coef_quality,xyz_z_model_mother_coef_quality),
	coo=c(rep("x",706),rep("y",706),rep("z",706)))

ggplot(dataTemp,aes(x=k,y=k_q,col=coo))+
geom_point(size=0.2,stroke=0.2)+
geom_abline(slope=1,inteercept=0,lwd=0.3)+
labs(x="all embryos",y="embryos with IS>=0.9",title=expression(italic(k)," comparison"))+
scale_x_continuous(limits=range(dataTemp[,1:2]),breaks=c(-1,0))+
scale_y_continuous(limits=range(dataTemp[,1:2]),breaks=c(-1,0))+
theme(legend.position=c(0.2,0.8))

ggsave("spatial_remove_bad_quality_embryos.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))






