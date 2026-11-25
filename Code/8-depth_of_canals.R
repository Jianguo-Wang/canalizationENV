

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


####################################################################################

library(ggplot2)
library(glmnet)
library(bayesmeta)


focal_projectFolderPath<-r"(E:\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/data/depth"))


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
	"F:/wangj/Desktop/temp/canalizationENV/data/spatial/xyz_x_model_mother_coef.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/spatial/xyz_y_model_mother_coef.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/spatial/xyz_z_model_mother_coef.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/cvglm_residual_training.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/cvglm_residual.RData"
	),projectFolderPath=focal_projectFolderPath)


embryo_lethal<-as.numeric(embryo_lethal)

	
	
############################################################################################



# depth


xyz_depth_abs<-t(sapply(1:dim(q1_cell350_alignment_xyz_noise)[2],function(i){

	print(i)
	
	noise_abs<-abs(q1_cell350_alignment_xyz_noise[,i])
	
	fit<-glm(embryo_lethal~noise_abs,family=binomial(link='logit'))

	d<-ifelse(coef(fit)[2]>0&(-coef(fit)[1]/coef(fit)[2])<=max(abs(noise_abs)),-coef(fit)[1]/coef(fit)[2],max(abs(noise_abs)))
	
	return(c(d,summary(fit)$coefficients[2,c(1:2,4)],summary(fit)$coefficients[1,c(1:2,4)]))
	
}))

colnames(xyz_depth_abs)[1]<-"depth"
colnames(xyz_depth_abs)[2]<-"beta1"
colnames(xyz_depth_abs)[5]<-"beta0"

rownames(xyz_depth_abs)<-colnames(q1_cell350_alignment_xyz_noise)

save(xyz_depth_abs,file=paste0(focal_projectFolderPath,"/data/depth/xyz_depth_abs.RData"))


depth_pos<-t(sapply(1:dim(q1_cell350_alignment_xyz_noise)[2],function(i){

	print(i)
	
	posInd<-which(q1_cell350_alignment_xyz_noise[,i]>0)
	noise_pos<-q1_cell350_alignment_xyz_noise[posInd,i]
	
	fit<-glm(embryo_lethal[posInd]~noise_pos,family=binomial(link='logit'))

	d<-ifelse(coef(fit)[2]>0&(-coef(fit)[1]/coef(fit)[2])<=max(abs(noise_pos)),-coef(fit)[1]/coef(fit)[2],max(abs(noise_pos)))
	
	return(c(d,summary(fit)$coefficients[2,c(1:2,4)],summary(fit)$coefficients[1,c(1:2,4)]))
	
}))


depth_neg<-t(sapply(1:dim(q1_cell350_alignment_xyz_noise)[2],function(i){

	print(i)
	
	negInd<-which(q1_cell350_alignment_xyz_noise[,i]<0)
	
	noise_neg<--q1_cell350_alignment_xyz_noise[negInd,i]
	
	fit<-glm(embryo_lethal[negInd]~noise_neg,family=binomial(link='logit'))

	d<-ifelse(coef(fit)[2]>0&(-coef(fit)[1]/coef(fit)[2])<=max(abs(noise_neg)),-coef(fit)[1]/coef(fit)[2],max(abs(noise_neg)))
	
	return(c(d,summary(fit)$coefficients[2,c(1:2,4)],summary(fit)$coefficients[1,c(1:2,4)]))
	
}))


#

pos_neg_test<-sapply(1:dim(depth_pos)[1],function(i){
	
	print(i)

	location_pos<--depth_pos[i,5]/depth_pos[i,2]
	location_neg<--depth_neg[i,5]/depth_neg[i,2]
	
	scale_pos<-1/depth_pos[i,2]
	scale_neg<-1/depth_neg[i,2]
	
	return(min(c(plogis(location_neg,location=location_pos,scale=scale_pos),1-plogis(location_neg,location=location_pos,scale=scale_pos))))
	
})

depth_summary<-data.frame(xyz_depth_abs=xyz_depth_abs[,1],depth_pos=depth_pos[,1],depth_neg=depth_neg[,1],p_abs=xyz_depth_abs[,4],p_pos=depth_pos[,4],p_neg=depth_neg[,4],p_diff=pos_neg_test)

save(depth_summary,file=paste0(focal_projectFolderPath,"/data/depth/depth_summary.RData"))


### ccl depth 


ccl_depth_abs<-t(sapply(1:dim(ccl_raw_nonNA_scaling_noise)[2],function(i){

	print(i)
	
	noise_abs<-abs(ccl_raw_nonNA_scaling_noise[,i])
	
	fit<-glm(embryo_lethal~noise_abs,family=binomial(link='logit'))

	d<-ifelse(coef(fit)[2]>0&(-coef(fit)[1]/coef(fit)[2])<=max(abs(noise_abs)),-coef(fit)[1]/coef(fit)[2],max(abs(noise_abs)))
	
	return(c(d,summary(fit)$coefficients[2,c(1:2,4)],summary(fit)$coefficients[1,c(1:2,4)]))
	
}))

colnames(ccl_depth_abs)[1]<-"depth"
colnames(ccl_depth_abs)[2]<-"beta1"
colnames(ccl_depth_abs)[5]<-"beta0"

rownames(ccl_depth_abs)<-colnames(ccl_raw_nonNA_scaling_noise)

save(ccl_depth_abs,file=paste0(focal_projectFolderPath,"/data/depth/ccl_depth_abs.RData"))









