

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
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",	
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_coef.RData"
	),projectFolderPath=focal_projectFolderPath)


###

TissueFC_group<-read.csv(paste0(focal_projectFolderPath,"/data/annotation/TissueFC_group.csv"))


##############################################################################################


feature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise)[1:714],"_x"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]

### 

dataTemp<-data.frame(mo_sis_pair,kx=xyz_x_model_mother_coef,ky=xyz_y_model_mother_coef,kz=xyz_z_model_mother_coef)

dataTemp$mo_tissue<-sapply(mo_sis_pair[,1],function(x){ifelse(length(which(TissueFC_group[,1]==x))==0,"plu",TissueFC_group[which(TissueFC_group[,1]==x),2])})
dataTemp$da1_tissue<-sapply(mo_sis_pair[,2],function(x){ifelse(length(which(TissueFC_group[,1]==x))==0,"plu",TissueFC_group[which(TissueFC_group[,1]==x),2])})
dataTemp$da2_tissue<-sapply(mo_sis_pair[,3],function(x){ifelse(length(which(TissueFC_group[,1]==x))==0,"plu",TissueFC_group[which(TissueFC_group[,1]==x),2])})

dataTemp$generation<-sapply(dataTemp[,2],function(x){rowColRec(which(lineage_info==x),dim(lineage_info)[1],dim(lineage_info)[2])[1]})

dataTemp[which(dataTemp$mo_tissue=="plu"|(dataTemp$mo_tissue==dataTemp$da2_tissue&dataTemp$mo_tissue==dataTemp$da2_tissue)),]
dataTemp[-which(dataTemp$mo_tissue=="plu"|(dataTemp$mo_tissue==dataTemp$da2_tissue&dataTemp$mo_tissue==dataTemp$da2_tissue)),]


#consider sister

dataTemp2<-data.frame(dataTemp[((dataTemp$da1_tissue=="plu"&dataTemp$da2_tissue!="plu")|(dataTemp$da1_tissue!="plu"&dataTemp$da2_tissue=="plu"))&dataTemp$generation==10,])

mos<-names(table(dataTemp2$mo))

wilcox.test(dataTemp2[dataTemp2$da1_tissue=="plu",4][match(mos,dataTemp2[dataTemp2$da1_tissue=="plu",1])],
dataTemp2[dataTemp2$da1_tissue!="plu",4][match(mos,dataTemp2[dataTemp2$da1_tissue!="plu",1])],paired=T)

wilcox.test(dataTemp2[dataTemp2$da1_tissue=="plu",5][match(mos,dataTemp2[dataTemp2$da1_tissue=="plu",1])],
dataTemp2[dataTemp2$da1_tissue!="plu",5][match(mos,dataTemp2[dataTemp2$da1_tissue!="plu",1])],paired=T)

wilcox.test(dataTemp2[dataTemp2$da1_tissue=="plu",6][match(mos,dataTemp2[dataTemp2$da1_tissue=="plu",1])],
dataTemp2[dataTemp2$da1_tissue!="plu",6][match(mos,dataTemp2[dataTemp2$da1_tissue!="plu",1])],paired=T)


#consider non-sister 

dataTemp2<-dataTemp[dataTemp$generation==10&dataTemp$da1_tissue!="plu",]
dataTemp3<-dataTemp[dataTemp$generation==10&dataTemp$da1_tissue=="plu",]
dataTemp3<-dataTemp3[-unique(na.omit(match(dataTemp2[,1],dataTemp3[,1]))),]

wilcox.test(dataTemp2[,4],dataTemp3[,4])
wilcox.test(dataTemp2[,5],dataTemp3[,5])
wilcox.test(dataTemp2[,6],dataTemp3[,6])













	