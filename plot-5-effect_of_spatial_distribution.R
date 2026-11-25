

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
	"E:/papers/canalizationENV/data/alignment/coo_list_q1_cell350_translation_rotation_scaling.RData",
	"E:/papers/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"E:/papers/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_coef.RData"
	),projectFolderPath=focal_projectFolderPath)



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


#################################################################################################




coo_x<-apply(q1_cell350_alignment_xyz[1:105,1:714],2,mean)
coo_y<-apply(q1_cell350_alignment_xyz[1:105,c(1:714)+714],2,mean)
coo_z<-apply(q1_cell350_alignment_xyz[1:105,c(1:714)+714*2],2,mean)

WT_POS_xyz<-data.frame(x=(coo_x[match(mo_sis_pair$mo,feature_names)]+coo_x[match(mo_sis_pair$da1,feature_names)])/2,
	y=(coo_y[match(mo_sis_pair$mo,feature_names)]+coo_y[match(mo_sis_pair$da1,feature_names)])/2,
	z=(coo_z[match(mo_sis_pair$mo,feature_names)]+coo_z[match(mo_sis_pair$da1,feature_names)])/2,
	kx=xyz_x_model_mother_coef,
	ky=xyz_y_model_mother_coef,
	kz=xyz_z_model_mother_coef)

save(WT_POS_xyz,file=paste0(focal_projectFolderPath,"/data/basic/WT_POS_xyz.RData"))

####################################################################################################

### spatial patterns

#kx~x

ggplot(WT_POS_xyz,aes(x=x,y=kx))+geom_point(size=0.1,stroke=0.1,col="red",alpha=0.5)+
geom_smooth(lwd=0.2)+
scale_x_continuous(breaks=c(-100,0,100))+
scale_y_continuous(limits=range(WT_POS_xyz[,4]),breaks=c(-1,-0.5,0))+
labs(x="x",y=expression(italic(k)[x]),title="")

ggsave("spatial_pattern_x_kx.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


#ky~z

ggplot(WT_POS_xyz,aes(x=z,y=ky))+geom_point(size=0.1,stroke=0.1,col="red",alpha=0.5)+
geom_smooth(lwd=0.2)+
scale_x_continuous(breaks=c(-25,0,25))+
scale_y_continuous(limits=range(WT_POS_xyz[,5]),breaks=c(-1,-0.5,0))+
labs(x="z",y=expression(italic(k)[y]),title="")

ggsave("spatial_pattern_z_ky.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



