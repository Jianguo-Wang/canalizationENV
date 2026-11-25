



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
	"F:/wangj/Desktop/temp/canalizationENV/functions/RowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_size.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/alignment/rotation_mat_list.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling_noise.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_duration_prcomp.RData",
	"E:/papers/canalizationENV/data/temporal/ccl_model_mother_coef.RData"
	),projectFolderPath=focal_projectFolderPath)


#

feature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise)[1:714],"_x"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]


###############################################################################################

#### 
embryo_measure_var<-t(sapply(1:length(rotation_mat_list),function(i){

	return(t(c(1,1,4.5)^2/12)%*%(rotation_mat_list[[i]]^2))
	
}))


### x

X_mo<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo


#

k_dilution_x<-sapply(1:dim(X_mo)[2],function(i){
	cov_xy<-cov(X_mo[,i],X_da1[,i])
	var_x<-var(X_mo[,i])
	var_dilution<-sum(embryo_measure_var[,1]/c(embryo_size[,1]^2))/1926
	return(cov_xy/(var_x-var_dilution))

})-1


### y

X_mo<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$mo,feature_names)+714]
X_da1<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da1,feature_names)+714]
X_da2<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da2,feature_names)+714]

delta_X<-X_da1-X_mo


#

k_dilution_y<-sapply(1:dim(X_mo)[2],function(i){
	cov_xy<-cov(X_mo[,i],X_da1[,i])
	var_x<-var(X_mo[,i])
	var_dilution<-sum(embryo_measure_var[,2]/c(embryo_size[,2]^2))/1926
	return(cov_xy/(var_x-var_dilution))

})-1


### z

X_mo<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$mo,feature_names)+714*2]
X_da1<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da1,feature_names)+714*2]
X_da2<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da2,feature_names)+714*2]

delta_X<-X_da1-X_mo


#

k_dilution_z<-sapply(1:dim(X_mo)[2],function(i){
	cov_xy<-cov(X_mo[,i],X_da1[,i])
	var_x<-var(X_mo[,i])
	var_dilution<-sum(embryo_measure_var[,3]/c(embryo_size[,3]^2))/1926
	return(cov_xy/(var_x-var_dilution))

})-1

save(k_dilution_x,file=paste0(focal_projectFolderPath,"/data/confounding/k_dilution_x.RData"))
save(k_dilution_y,file=paste0(focal_projectFolderPath,"/data/confounding/k_dilution_y.RData"))
save(k_dilution_z,file=paste0(focal_projectFolderPath,"/data/confounding/k_dilution_z.RData"))


### CCL 


feature_names<-unlist(strsplit(colnames(ccl_raw_nonNA_scaling_noise),"_ccl"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]


X_mo<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$mo,feature_names)]
X_da1<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$da1,feature_names)]
X_da2<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo


X_mo_pre<-ccl_raw_nonNA_scaling[,match(mo_sis_pair$mo,feature_names)]
X_da1_pre<-ccl_raw_nonNA_scaling[,match(mo_sis_pair$da1,feature_names)]
X_da2_pre<-ccl_raw_nonNA_scaling[,match(mo_sis_pair$da2,feature_names)]

delta_X_pre<-X_da1_pre-X_mo_pre



k_regression=ccl_model_mother_coef


k_dilution_CCL=sapply(1:dim(X_mo)[2],function(i){
	cov_xy<-cov(X_mo[,i],X_da1[,i])
	var_x<-var(X_mo[,i])
	var_y<-var(X_da1[,i])
	
	#var_dilution<-sum(0.5^2/12/c(embryo_duration_prcomp^2))/1926/(mean(X_mo_pre[,i])^2)
	var_dilution<-sum(2.5^2/12/c(embryo_duration_prcomp^2))/1926/(mean(X_mo_pre[,i])^2)
	
	beta1<-cov_xy/(var_x-var_dilution)
	
	return(beta1)
	
})-1

save(k_dilution_CCL,file=paste0(focal_projectFolderPath,"/data/confounding/k_dilution_CCL.RData"))



############################################################################################

par(mfrow=c(2,2))


plot(k_dilution_x,xyz_x_model_mother_coef,ylab="k_x")
abline(0,1)
cor(k_dilution_x,xyz_x_model_mother_coef)

plot(k_dilution_y,xyz_y_model_mother_coef,ylab="k_y")
abline(0,1)
cor(k_dilution_y,xyz_y_model_mother_coef)

plot(k_dilution_z,xyz_z_model_mother_coef,ylab="k_z")
abline(0,1)
cor(k_dilution_z,xyz_z_model_mother_coef)


### x 


dataTemp<-data.frame(dilution=k_dilution_x,original=xyz_x_model_mother_coef)
ggplot(dataTemp,aes(x=dilution,y=original))+
geom_point(size=0.2,stroke=0.1)+
geom_abline(slope=1,intercept=0,lwd=0.3)+
scale_x_continuous(breaks=c(-1,-0.5,0))+
scale_y_continuous(breaks=c(-1,-0.5,0))+
labs(x=expression(italic(k)[dilution]),y=expression(italic(k)[original]),title=expression(italic(k)[x]))

ggsave("effect_dilution_effect_x.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### y 


dataTemp<-data.frame(dilution=k_dilution_y,original=xyz_y_model_mother_coef)
ggplot(dataTemp,aes(x=dilution,y=original))+
geom_point(size=0.2,stroke=0.1)+
geom_abline(slope=1,intercept=0,lwd=0.3)+
scale_x_continuous(breaks=c(-1,-0.5,0))+
scale_y_continuous(breaks=c(-1,-0.5,0))+
labs(x=expression(italic(k)[dilution]),y=expression(italic(k)[original]),title=expression(italic(k)[y]))

ggsave("effect_dilution_effect_y.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### z


dataTemp<-data.frame(dilution=k_dilution_z,original=xyz_z_model_mother_coef)
ggplot(dataTemp,aes(x=dilution,y=original))+
geom_point(size=0.2,stroke=0.1)+
geom_abline(slope=1,intercept=0,lwd=0.3)+
scale_x_continuous(breaks=c(-1,-0.5,0),limits=c(min(dataTemp),max(dataTemp)))+
scale_y_continuous(breaks=c(-1,-0.5,0),limits=c(min(dataTemp),max(dataTemp)))+
labs(x=expression(italic(k)[dilution]),y=expression(italic(k)[original]),title=expression(italic(k)[z]))

ggsave("effect_dilution_effect_z.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### CCL


dataTemp<-data.frame(dilution=k_dilution_CCL,original=ccl_model_mother_coef)
ggplot(dataTemp,aes(x=dilution,y=original))+
geom_point(size=0.2,stroke=0.1)+
geom_abline(slope=1,intercept=0,lwd=0.3)+
scale_x_continuous(breaks=c(-1,-0.5,0),limits=c(min(dataTemp),max(dataTemp)))+
scale_y_continuous(breaks=c(-1,-0.5,0),limits=c(min(dataTemp),max(dataTemp)))+
labs(x=expression(italic(k)[dilution]),y=expression(italic(k)[original]),title=expression(italic(k)[CCL]))

ggsave("effect_dilution_effect_CCL.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


############################################################################################

### xyz 


dataTemp<-data.frame(dilution=c(k_dilution_x,k_dilution_y,k_dilution_z),
	original=c(xyz_x_model_mother_coef,xyz_y_model_mother_coef,xyz_z_model_mother_coef),
	coo=c(rep("x",length(k_dilution_x)),rep("y",length(k_dilution_y)),rep("z",length(k_dilution_z))))
	
ggplot(dataTemp,aes(x=dilution,y=original,col=coo))+
geom_point(size=0.2,stroke=0.1,show.legend=FALSE)+
geom_abline(slope=1,intercept=0,lwd=0.3,lty=2)+
scale_x_continuous(breaks=c(-1,-0.5,0),limits=range(dataTemp[,1:2]))+
scale_y_continuous(breaks=c(-1,-0.5,0),limits=range(dataTemp[,1:2]))+
labs(x=expression(italic(k)[dilution]),y=expression(italic(k)[original]),title="")

ggsave("effect_dilution_effect_xyz.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



