





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
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_size.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_R.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling_noise.RData",
	"E:/papers/canalizationENV/data/temporal/ccl_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_R_shuffle.RData",	
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_coef_shuffle.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mo_da1_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mo_da1_pvalue.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mo_da1_R_shuffle.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mo_da1_pvalue_shuffle.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_R_shuffle.RData",	
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_coef_shuffle.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mo_da1_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mo_da1_pvalue.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mo_da1_R_shuffle.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mo_da1_pvalue_shuffle.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_R_shuffle.RData",	
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_coef_shuffle.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mo_da1_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mo_da1_pvalue.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mo_da1_R_shuffle.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mo_da1_pvalue_shuffle.RData"	
	),projectFolderPath=focal_projectFolderPath)


# xyz : b=r*Sy/Sx<Sy/Sx


feature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise)[1:714],"_x"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]

#x

X_mo_x<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$mo,feature_names)]
X_da1_x<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da1,feature_names)]

delta_X_x<-X_da1_x-X_mo_x

#y

X_mo_y<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$mo,feature_names)+714]
X_da1_y<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da1,feature_names)+714]

delta_X_y<-X_da1_y-X_mo_y

#z

X_mo_z<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$mo,feature_names)+714*2]
X_da1_z<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da1,feature_names)+714*2]

delta_X_z<-X_da1_z-X_mo_z

# combine 

X_mo_xyz<-cbind(X_mo_x,X_mo_y,X_mo_z)
X_da1_xyz<-cbind(X_da1_x,X_da1_y,X_da1_z)

plot(apply(X_mo_xyz,2,sd),apply(X_da1_xyz,2,sd))
abline(0,1)
t.test(apply(X_mo_xyz,2,sd),apply(X_da1_xyz,2,sd),paired=T)

#

plot(apply(X_da1_x,2,sd)/apply(X_mo_x,2,sd),xyz_x_model_mother_coef+1)
abline(0,1)
plot(apply(X_da1_y,2,sd)/apply(X_mo_y,2,sd),xyz_y_model_mother_coef+1)
abline(0,1)
plot(apply(X_da1_z,2,sd)/apply(X_mo_z,2,sd),xyz_z_model_mother_coef+1)
abline(0,1)



### CCL 

feature_names_CCL<-unlist(strsplit(colnames(ccl_raw_nonNA_scaling_noise),"_ccl"))

temp1<-unlist(lapply(feature_names_CCL,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names_CCL,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names_CCL,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]

#x

X_mo_CCL<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$mo,feature_names_CCL)]
X_da1_CCL<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$da1,feature_names_CCL)]

delta_X_CCL<-X_da1_CCL-X_mo_CCL


plot(apply(X_mo_CCL,2,sd),apply(X_da1_CCL,2,sd))
abline(0,1)
t.test(apply(X_mo_CCL,2,sd),apply(X_da1_CCL,2,sd),paired=T)

plot(apply(X_da1_CCL,2,sd)/apply(X_mo_CCL,2,sd),ccl_model_mother_coef+1)
abline(0,1)

####


dataTemp<-data.frame(k=c(xyz_x_model_mother_coef,xyz_y_model_mother_coef,xyz_z_model_mother_coef,ccl_model_mother_coef),
	uplimit=c(apply(X_da1_x,2,sd)/apply(X_mo_x,2,sd),apply(X_da1_y,2,sd)/apply(X_mo_y,2,sd),apply(X_da1_z,2,sd)/apply(X_mo_z,2,sd),apply(X_da1_CCL,2,sd)/apply(X_mo_CCL,2,sd)),
	feature=factor(c(rep("x",706),rep("y",706),rep("z",706),rep("CCL",345)),level=c("x","y","z","CCL")))
	
ggplot(dataTemp,aes(x=uplimit,y=k+1))+
geom_point(size=0.5,col=colors()[182],stroke=0.2)+
geom_abline(slope=1,intercept=0,lty=2,lwd=0.3)+
geom_vline(xintercept=1,lty=2,lwd=0.3)+
geom_hline(yintercept=1,lty=2,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,1:2]))+
scale_y_continuous(limits=range(dataTemp[,1:2]))+
facet_wrap(vars(feature))+
labs(x=expression(paste("Upper limit (",σ[D]/σ[M],")")),y=expression(paste(italic(k),"+1")))

ggsave("effect_mo_da_sd_cv.pdf",width=unit(4,"in"),height=unit(4,"in"))


	
ggplot(dataTemp[dataTemp$feature!="CCL",],aes(x=uplimit,y=k+1))+
geom_point(size=0.5,col=colors()[182],stroke=0.2)+
geom_abline(slope=1,intercept=0,lty=2,lwd=0.3)+
geom_vline(xintercept=1,lty=2,lwd=0.3)+
geom_hline(yintercept=1,lty=2,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,1:2]))+
scale_y_continuous(limits=range(dataTemp[,1:2]))+
facet_wrap(vars(feature))+
labs(x=expression(paste("Upper limit (",σ[D]/σ[M],")")),y=expression(paste(italic(k),"+1")))

ggsave("effect_mo_da_sd_kxyz.pdf",width=unit(6,"in"),height=unit(2.3,"in"))


	
ggplot(dataTemp[dataTemp$feature=="CCL",],aes(x=uplimit,y=k+1))+
geom_point(size=0.2,col=colors()[182],stroke=0.2)+
geom_abline(slope=1,intercept=0,lty=2,lwd=0.3)+
geom_vline(xintercept=1,lty=2,lwd=0.3)+
geom_hline(yintercept=1,lty=2,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,1:2]))+
scale_y_continuous(limits=range(dataTemp[,1:2]))+
labs(x=expression(paste("Upper limit (",σ[D]/σ[M],")")),y=expression(paste(italic(k),"+1")),title="CCL")

ggsave("effect_mo_da_sd_kCCL.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))






