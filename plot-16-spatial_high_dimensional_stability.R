


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
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_full_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_full_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_full_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_full_R2.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_full_R2.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_full_R2.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_block_descendant_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_block_descendant_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_block_descendant_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_block_descendant_R2.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_block_descendant_R2.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_block_descendant_R2.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_M1M11_full_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_M1M11_full_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_M1M11_full_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_M1M11g9_full_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_M1M11g9_full_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_M1M11g9_full_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_M1M11_full_R2.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_M1M11_full_R2.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_M1M11_full_R2.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_M1M11g9_full_R2.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_M1M11g9_full_R2.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_M1M11g9_full_R2.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_nested_without_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_nested_without_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_nested_without_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_nested_without_sister_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_nested_without_sister_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_nested_without_sister_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_sister_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_sister_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_sister_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_sister_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_sister_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_sister_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_sister_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_sister_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_sister_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_sister_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_sister_R.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_sister_R.RData"
	),projectFolderPath=focal_projectFolderPath)


###


feature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise)[1:714],"_x"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]


### coefficients of relatives 


# coef : x


full_coef_mother_x<-sapply(1:dim(mo_sis_pair)[1],function(i){
	return(setdiff(xyz_x_model_full_coef[which(rownames(xyz_x_model_full_coef)==paste0(mo_sis_pair$mo[i],"_x")),i],0))
})


full_coef_sister_x<-sapply(1:dim(mo_sis_pair)[1],function(i){
	temp<-setdiff(xyz_x_model_full_coef[which(rownames(xyz_x_model_full_coef)==paste0(mo_sis_pair$da2[i],"_x")),i],0)
	return(ifelse(is.null(temp),NA,temp))
})



full_coef_aunt_x<-sapply(1:dim(mo_sis_pair)[1],function(i){
	focal_aunt<-mo_sis_pair$da2[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	focal_aunt_ind<-which(rownames(xyz_x_model_full_coef)==paste0(focal_aunt,"_x"))
	if(length(focal_aunt_ind)==0){
		return(NA)
	}else{
		temp<-setdiff(xyz_x_model_full_coef[focal_aunt_ind,i],0)
		return(ifelse(is.null(temp),NA,temp))	
	}
})


full_coef_grandmother_x<-sapply(1:dim(mo_sis_pair)[1],function(i){
	focal_grandmother<-mo_sis_pair$mo[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	focal_grandmother_ind<-which(rownames(xyz_x_model_full_coef)==paste0(focal_grandmother,"_x"))
	if(length(focal_grandmother_ind)==0){
		return(NA)
	}else{
		temp<-setdiff(xyz_x_model_full_coef[focal_grandmother_ind,i],0)
		return(ifelse(is.null(temp),NA,temp))	
	}
})



full_coef_cousin_x<-unlist(lapply(1:dim(mo_sis_pair)[1],function(i){
	focal_aunt<-mo_sis_pair$da2[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	if(is.na(focal_aunt)){
		return(NA)
	}else{
		tempInd<-which(mo_sis_pair$mo==focal_aunt)
		if(length(tempInd)==0){
			return(NA)
		}else{
			focal_cousin<-unique(c(mo_sis_pair$da1[tempInd],mo_sis_pair$da2[tempInd]))
			focal_cousin_ind<-match(paste0(focal_cousin,"_x"),rownames(xyz_x_model_full_coef))
			if(length(focal_cousin_ind)==0){
				return(NA)
			}else{
				temp<-setdiff(xyz_x_model_full_coef[focal_cousin_ind,i],0)
				if(is.null(temp)){
					return(NA)
				}else{
					return(temp)
				}
			}			
		}
	}
}))


full_coef_other_x<-unlist(lapply(1:dim(mo_sis_pair)[1],function(i){

	ind_mo<-which(rownames(xyz_x_model_full_coef)==paste0(mo_sis_pair$mo[i],"_x"))
	ind_sis<-which(rownames(xyz_x_model_full_coef)==paste0(mo_sis_pair$da2[i],"_x"))
	
	focal_aunt<-mo_sis_pair$da2[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	ind_aunt<-which(rownames(xyz_x_model_full_coef)==paste0(focal_aunt,"_x"))
	
	focal_grandmother<-mo_sis_pair$mo[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	ind_grandmother<-which(rownames(xyz_x_model_full_coef)==paste0(focal_grandmother,"_x"))	
	
	ind_cousin<-NA
	
	if(!is.na(focal_aunt)){
		tempInd<-which(mo_sis_pair$mo==focal_aunt)
		if(length(tempInd)!=0){
			focal_cousin<-unique(c(mo_sis_pair$da1[tempInd],mo_sis_pair$da2[tempInd]))
			ind_cousin<-unlist(lapply(paste0(focal_cousin,"_x"),function(x){which(rownames(xyz_x_model_full_coef)==x)}))
		}
	}
	
	ind_non_relative<-setdiff(2:707,c(ind_mo,ind_sis,ind_aunt,ind_grandmother,ind_cousin))
	temp<-setdiff(xyz_x_model_full_coef[ind_non_relative,i],0)
	
	return(temp)

}))


# coef : y 


full_coef_mother_y<-sapply(1:dim(mo_sis_pair)[1],function(i){
	return(setdiff(xyz_y_model_full_coef[which(rownames(xyz_y_model_full_coef)==paste0(mo_sis_pair$mo[i],"_y")),i],0))
})


full_coef_sister_y<-sapply(1:dim(mo_sis_pair)[1],function(i){
	temp<-setdiff(xyz_y_model_full_coef[which(rownames(xyz_y_model_full_coef)==paste0(mo_sis_pair$da2[i],"_y")),i],0)
	return(ifelse(is.null(temp),NA,temp))
})




full_coef_aunt_y<-sapply(1:dim(mo_sis_pair)[1],function(i){
	focal_aunt<-mo_sis_pair$da2[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	focal_aunt_ind<-which(rownames(xyz_y_model_full_coef)==paste0(focal_aunt,"_y"))
	if(length(focal_aunt_ind)==0){
		return(NA)
	}else{
		temp<-setdiff(xyz_y_model_full_coef[focal_aunt_ind,i],0)
		return(ifelse(is.null(temp),NA,temp))	
	}
})


full_coef_grandmother_y<-sapply(1:dim(mo_sis_pair)[1],function(i){
	focal_grandmother<-mo_sis_pair$mo[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	focal_grandmother_ind<-which(rownames(xyz_y_model_full_coef)==paste0(focal_grandmother,"_y"))
	if(length(focal_grandmother_ind)==0){
		return(NA)
	}else{
		temp<-setdiff(xyz_y_model_full_coef[focal_grandmother_ind,i],0)
		return(ifelse(is.null(temp),NA,temp))	
	}
})




full_coef_cousin_y<-unlist(lapply(1:dim(mo_sis_pair)[1],function(i){
	focal_aunt<-mo_sis_pair$da2[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	if(is.na(focal_aunt)){
		return(NA)
	}else{
		tempInd<-which(mo_sis_pair$mo==focal_aunt)
		if(length(tempInd)==0){
			return(NA)
		}else{
			focal_cousin<-unique(c(mo_sis_pair$da1[tempInd],mo_sis_pair$da2[tempInd]))
			focal_cousin_ind<-match(paste0(focal_cousin,"_y"),rownames(xyz_y_model_full_coef))
			if(length(focal_cousin_ind)==0){
				return(NA)
			}else{
				temp<-setdiff(xyz_y_model_full_coef[focal_cousin_ind,i],0)
				if(is.null(temp)){
					return(NA)
				}else{
					return(temp)
				}
			}			
		}
	}
}))


full_coef_other_y<-unlist(lapply(1:dim(mo_sis_pair)[1],function(i){

	ind_mo<-which(rownames(xyz_y_model_full_coef)==paste0(mo_sis_pair$mo[i],"_y"))
	ind_sis<-which(rownames(xyz_y_model_full_coef)==paste0(mo_sis_pair$da2[i],"_y"))
	
	focal_aunt<-mo_sis_pair$da2[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	ind_aunt<-which(rownames(xyz_y_model_full_coef)==paste0(focal_aunt,"_y"))
	
	focal_grandmother<-mo_sis_pair$mo[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	ind_grandmother<-which(rownames(xyz_y_model_full_coef)==paste0(focal_grandmother,"_y"))	
	
	ind_cousin<-NA
	
	if(!is.na(focal_aunt)){
		tempInd<-which(mo_sis_pair$mo==focal_aunt)
		if(length(tempInd)!=0){
			focal_cousin<-unique(c(mo_sis_pair$da1[tempInd],mo_sis_pair$da2[tempInd]))
			ind_cousin<-unlist(lapply(paste0(focal_cousin,"_y"),function(x){which(rownames(xyz_y_model_full_coef)==x)}))
		}
	}
	
	ind_non_relative<-setdiff(2:707,c(ind_mo,ind_sis,ind_aunt,ind_grandmother,ind_cousin))
	temp<-setdiff(xyz_y_model_full_coef[ind_non_relative,i],0)
	
	return(temp)

}))



# coef : z 


full_coef_mother_z<-sapply(1:dim(mo_sis_pair)[1],function(i){
	return(setdiff(xyz_z_model_full_coef[which(rownames(xyz_z_model_full_coef)==paste0(mo_sis_pair$mo[i],"_z")),i],0))
})


full_coef_sister_z<-sapply(1:dim(mo_sis_pair)[1],function(i){
	temp<-setdiff(xyz_z_model_full_coef[which(rownames(xyz_z_model_full_coef)==paste0(mo_sis_pair$da2[i],"_z")),i],0)
	return(ifelse(is.null(temp),NA,temp))
})


full_coef_aunt_z<-sapply(1:dim(mo_sis_pair)[1],function(i){
	focal_aunt<-mo_sis_pair$da2[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	focal_aunt_ind<-which(rownames(xyz_z_model_full_coef)==paste0(focal_aunt,"_z"))
	if(length(focal_aunt_ind)==0){
		return(NA)
	}else{
		temp<-setdiff(xyz_z_model_full_coef[focal_aunt_ind,i],0)
		return(ifelse(is.null(temp),NA,temp))	
	}
})


full_coef_grandmother_z<-sapply(1:dim(mo_sis_pair)[1],function(i){
	focal_grandmother<-mo_sis_pair$mo[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	focal_grandmother_ind<-which(rownames(xyz_z_model_full_coef)==paste0(focal_grandmother,"_z"))
	if(length(focal_grandmother_ind)==0){
		return(NA)
	}else{
		temp<-setdiff(xyz_z_model_full_coef[focal_grandmother_ind,i],0)
		return(ifelse(is.null(temp),NA,temp))	
	}
})


full_coef_cousin_z<-unlist(lapply(1:dim(mo_sis_pair)[1],function(i){
	focal_aunt<-mo_sis_pair$da2[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	if(is.na(focal_aunt)){
		return(NA)
	}else{
		tempInd<-which(mo_sis_pair$mo==focal_aunt)
		if(length(tempInd)==0){
			return(NA)
		}else{
			focal_cousin<-unique(c(mo_sis_pair$da1[tempInd],mo_sis_pair$da2[tempInd]))
			focal_cousin_ind<-match(paste0(focal_cousin,"_z"),rownames(xyz_z_model_full_coef))
			if(length(focal_cousin_ind)==0){
				return(NA)
			}else{
				temp<-setdiff(xyz_z_model_full_coef[focal_cousin_ind,i],0)
				if(is.null(temp)){
					return(NA)
				}else{
					return(temp)
				}
			}			
		}
	}
}))


full_coef_other_z<-unlist(lapply(1:dim(mo_sis_pair)[1],function(i){

	ind_mo<-which(rownames(xyz_z_model_full_coef)==paste0(mo_sis_pair$mo[i],"_z"))
	ind_sis<-which(rownames(xyz_z_model_full_coef)==paste0(mo_sis_pair$da2[i],"_z"))
	
	focal_aunt<-mo_sis_pair$da2[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	ind_aunt<-which(rownames(xyz_z_model_full_coef)==paste0(focal_aunt,"_z"))
	
	focal_grandmother<-mo_sis_pair$mo[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	ind_grandmother<-which(rownames(xyz_z_model_full_coef)==paste0(focal_grandmother,"_z"))	
	
	ind_cousin<-NA
	
	if(!is.na(focal_aunt)){
		tempInd<-which(mo_sis_pair$mo==focal_aunt)
		if(length(tempInd)!=0){
			focal_cousin<-unique(c(mo_sis_pair$da1[tempInd],mo_sis_pair$da2[tempInd]))
			ind_cousin<-unlist(lapply(paste0(focal_cousin,"_z"),function(x){which(rownames(xyz_z_model_full_coef)==x)}))
		}
	}
	
	ind_non_relative<-setdiff(2:707,c(ind_mo,ind_sis,ind_aunt,ind_grandmother,ind_cousin))
	temp<-setdiff(xyz_z_model_full_coef[ind_non_relative,i],0)
	
	return(temp)

}))


###

coef_df<-data.frame(mo_sis_pair,mo_x=diag(xyz_x_model_full_coef[-1,]),mo_y=diag(xyz_y_model_full_coef[-1,]),mo_z=diag(xyz_z_model_full_coef[-1,]),
	sis_x=sapply(1:dim(mo_sis_pair)[1],function(i){xyz_x_model_full_coef[match(paste0(mo_sis_pair[i,3],"_x"),rownames(xyz_x_model_full_coef)),i]}),
	sis_y=sapply(1:dim(mo_sis_pair)[1],function(i){xyz_y_model_full_coef[match(paste0(mo_sis_pair[i,3],"_y"),rownames(xyz_y_model_full_coef)),i]}),
	sis_z=sapply(1:dim(mo_sis_pair)[1],function(i){xyz_z_model_full_coef[match(paste0(mo_sis_pair[i,3],"_z"),rownames(xyz_z_model_full_coef)),i]})
	)


### R2 

#### x


#

X_mo<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo


R2_mother_x<-xyz_x_model_mother_R^2
R2_sister_x<-xyz_x_model_sister_R^2
R2_mother_sister_x<-xyz_x_model_mother_sister_R^2

pred_full_x<-cbind(rep(1,2039),X_mo)%*%xyz_x_model_full_coef
R2_full_x<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],pred_full_x[,i])})^2

pred_without_mother_x<-cbind(rep(1,2039),X_mo)%*%xyz_x_model_nested_without_mother_coef
R2_without_mother_x<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],pred_without_mother_x[,i])^2})

pred_without_sister_x<-cbind(rep(1,2039),X_mo)%*%xyz_x_model_nested_without_sister_coef
R2_without_sister_x<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],pred_without_sister_x[,i])^2})

with_sister_ind_x<-which(!is.na(R2_mother_sister_x))


#### y

#

X_mo<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$mo,feature_names)+714]
X_da1<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da1,feature_names)+714]
X_da2<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da2,feature_names)+714]

delta_X<-X_da1-X_mo


R2_mother_y<-xyz_y_model_mother_R^2
R2_sister_y<-xyz_y_model_sister_R^2
R2_mother_sister_y<-xyz_y_model_mother_sister_R^2

pred_full_y<-cbind(rep(1,2039),X_mo)%*%xyz_y_model_full_coef
R2_full_y<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],pred_full_y[,i])})^2

pred_without_mother_y<-cbind(rep(1,2039),X_mo)%*%xyz_y_model_nested_without_mother_coef
R2_without_mother_y<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],pred_without_mother_y[,i])^2})

pred_without_sister_y<-cbind(rep(1,2039),X_mo)%*%xyz_y_model_nested_without_sister_coef
R2_without_sister_y<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],pred_without_sister_y[,i])^2})

with_sister_ind_y<-which(!is.na(R2_mother_sister_y))

#### z

#

X_mo<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$mo,feature_names)+714*2]
X_da1<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da1,feature_names)+714*2]
X_da2<-q1_cell350_alignment_xyz_noise[,match(mo_sis_pair$da2,feature_names)+714*2]

delta_X<-X_da1-X_mo


R2_mother_z<-xyz_z_model_mother_R^2
R2_sister_z<-xyz_z_model_sister_R^2
R2_mother_sister_z<-xyz_z_model_mother_sister_R^2

pred_full_z<-cbind(rep(1,2039),X_mo)%*%xyz_z_model_full_coef
R2_full_z<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],pred_full_z[,i])})^2

pred_without_mother_z<-cbind(rep(1,2039),X_mo)%*%xyz_z_model_nested_without_mother_coef
R2_without_mother_z<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],pred_without_mother_z[,i])^2})

pred_without_sister_z<-cbind(rep(1,2039),X_mo)%*%xyz_z_model_nested_without_sister_coef
R2_without_sister_z<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],pred_without_sister_z[,i])^2})

with_sister_ind_z<-which(!is.na(R2_mother_sister_z))



#######################################################################################################

par(mfrow=c(2,3))
plot(xyz_x_model_mother_coef,diag(xyz_x_model_full_coef[-1,]),pch=16)
abline(0,1)
plot(xyz_y_model_mother_coef,diag(xyz_y_model_full_coef[-1,]),pch=16)
abline(0,1)
plot(xyz_z_model_mother_coef,diag(xyz_z_model_full_coef[-1,]),pch=16)
abline(0,1)
cor.test(xyz_x_model_mother_coef,diag(xyz_x_model_full_coef[-1,]))
cor.test(xyz_y_model_mother_coef,diag(xyz_y_model_full_coef[-1,]))
cor.test(xyz_z_model_mother_coef,diag(xyz_z_model_full_coef[-1,]))

plot(coef_df$mo_x,coef_df$sis_x,pch=16)
plot(coef_df$mo_y,coef_df$sis_y,pch=16)
plot(coef_df$mo_z,coef_df$sis_z,pch=16)
cor.test(coef_df$mo_x,coef_df$sis_x)
cor.test(coef_df$mo_y,coef_df$sis_y)
cor.test(coef_df$mo_z,coef_df$sis_z)

dataTemp<-data.frame(coef_mo=c(coef_df$mo_x,coef_df$mo_y,coef_df$mo_z),
	coef_sis=c(coef_df$sis_x,coef_df$sis_y,coef_df$sis_z),
	coo=c(rep("x",706),rep("y",706),rep("z",706)))

#x 
	
ggplot(na.omit(dataTemp[dataTemp$coo=="x",]),aes(x=coef_mo,y=coef_sis))+
geom_point(size=0.5,stroke=0.5,show.legend=FALSE,col=colors()[182])+
geom_smooth(method="lm",show.legend=FALSE)+
labs(x="Coefficients of mother",y="Coefficients of sister",title="x")

ggsave("spatial_full_model_coef_mo_sis_comparison_x.pdf",width=unit(2,"in"),height=unit(2.1,"in"))

#y 
	
ggplot(na.omit(dataTemp[dataTemp$coo=="y",]),aes(x=coef_mo,y=coef_sis))+
geom_point(size=0.5,stroke=0.5,show.legend=FALSE,col=colors()[182])+
geom_smooth(method="lm",show.legend=FALSE)+
labs(x="Coefficients of mother",y="Coefficients of sister",title="y")

ggsave("spatial_full_model_coef_mo_sis_comparison_y.pdf",width=unit(2,"in"),height=unit(2.1,"in"))

#z 
	
ggplot(na.omit(dataTemp[dataTemp$coo=="z",]),aes(x=coef_mo,y=coef_sis))+
geom_point(size=0.5,stroke=0.5,show.legend=FALSE,col=colors()[182])+
geom_smooth(method="lm",show.legend=FALSE)+
labs(x="Coefficients of mother",y="Coefficients of sister",title="z")

ggsave("spatial_full_model_coef_mo_sis_comparison_z.pdf",width=unit(2,"in"),height=unit(2.1,"in"))


### coef of x
da_cell<-unlist(strsplit(rownames(t(xyz_x_model_full_coef[-1,])),"_x"))
orderInd<-na.omit(match(unlist(lapply(1:10,function(i){names(table(lineage_info[i,]))})),da_cell))

temp<-t(xyz_x_model_full_coef[-1,])
temp<-temp[orderInd,orderInd]

colnames(temp)<-unlist(strsplit(colnames(temp),"_x"))
rownames(temp)<-unlist(strsplit(rownames(temp),"_x"))

#nonSisInd<-which(!sapply(colnames(temp),function(x){ifelse(rownames(temp)==cell_mo_da1_da2[match(x,cell_mo_da1_da2[,2]),3],TRUE,FALSE)}))

rownames(temp)<-c(1:706)
colnames(temp)<-c(1:706)


temp_melt<-melt(temp,c("daughter","mother"))

temp_melt<-temp_melt[which(temp_melt$value!=0),]
temp_melt[which(temp_melt[,3]>0),3]<-temp_melt[which(temp_melt[,3]>0),3]*2

ggplot(temp_melt,aes(x=mother,y=daughter,col=value))+geom_point(size=0.2,stroke=0.05)+
scale_colour_gradient2(low="blue",high="red",mid="white",midpoint=0)+
scale_x_continuous(breaks=c(1,706))+
scale_y_continuous(breaks=c(1,706))+
labs(x=expression(italic(X)[M]),y=expression(paste(italic(Delta),italic(X))),title=expression(italic(K)[x]))


ggsave("spatial_full_model_coef_Kx.pdf",width=unit(2.6,"in"),height=unit(2.1,"in"))


### coef of y

temp<-t(xyz_y_model_full_coef[-1,])
temp<-temp[orderInd,orderInd]

colnames(temp)<-unlist(strsplit(colnames(temp),"_y"))
rownames(temp)<-unlist(strsplit(rownames(temp),"_y"))
rownames(temp)<-c(1:706)
colnames(temp)<-c(1:706)


temp_melt<-melt(temp,c("daughter","mother"))

temp_melt<-temp_melt[which(temp_melt$value!=0),]
temp_melt[which(temp_melt[,3]>0),3]<-temp_melt[which(temp_melt[,3]>0),3]*2

ggplot(temp_melt,aes(x=mother,y=daughter,col=value))+geom_point(size=0.2,stroke=0.05)+
scale_colour_gradient2(low="blue",high="red",mid="white",midpoint=0)+
scale_x_continuous(breaks=c(1,706))+
scale_y_continuous(breaks=c(1,706))+
labs(x=expression(italic(Y)[M]),y=expression(paste(italic(Delta),italic(Y))),title=expression(italic(K)[y]))


ggsave("spatial_full_model_coef_Ky.pdf",width=unit(2.6,"in"),height=unit(2.1,"in"))


### coef of z

temp<-t(xyz_z_model_full_coef[-1,])
temp<-temp[orderInd,orderInd]

colnames(temp)<-unlist(strsplit(colnames(temp),"_z"))
rownames(temp)<-unlist(strsplit(rownames(temp),"_z"))
rownames(temp)<-c(1:706)
colnames(temp)<-c(1:706)


temp_melt<-melt(temp,c("daughter","mother"))

temp_melt<-temp_melt[which(temp_melt$value!=0),]
temp_melt[which(temp_melt[,3]>0),3]<-temp_melt[which(temp_melt[,3]>0),3]*2

ggplot(temp_melt,aes(x=mother,y=daughter,col=value))+geom_point(size=0.2,stroke=0.05)+
scale_colour_gradient2(low="blue",high="red",mid="white",midpoint=0)+
scale_x_continuous(breaks=c(1,706))+
scale_y_continuous(breaks=c(1,706))+
labs(x=expression(italic(Z)[M]),y=expression(paste(italic(Delta),italic(Z))),title=expression(italic(K)[z]))


ggsave("spatial_full_model_coef_Kz.pdf",width=unit(2.6,"in"),height=unit(2.1,"in"))



### coefficents of relatives 


dataTemp<-data.frame(Coef=c(full_coef_mother_x,full_coef_sister_x,full_coef_aunt_x,full_coef_cousin_x,full_coef_grandmother_x,full_coef_other_x,
							full_coef_mother_y,full_coef_sister_y,full_coef_aunt_y,full_coef_cousin_y,full_coef_grandmother_y,full_coef_other_y,
							full_coef_mother_z,full_coef_sister_z,full_coef_aunt_z,full_coef_cousin_z,full_coef_grandmother_z,full_coef_other_z),
					 Coo=factor(c(rep("x",length(full_coef_mother_x)+length(full_coef_sister_x)+length(full_coef_aunt_x)+length(full_coef_cousin_x)+length(full_coef_other_x)+length(full_coef_grandmother_x)),
								  rep("y",length(full_coef_mother_y)+length(full_coef_sister_y)+length(full_coef_aunt_y)+length(full_coef_cousin_y)+length(full_coef_other_y)+length(full_coef_grandmother_y)),
								  rep("z",length(full_coef_mother_z)+length(full_coef_sister_z)+length(full_coef_aunt_z)+length(full_coef_cousin_z)+length(full_coef_other_z)+length(full_coef_grandmother_z))),level=c("x","y","z")),
					 Relative=factor(c(rep("Mother",length(full_coef_mother_x)),rep("Sister",length(full_coef_sister_x)),rep("Aunt",length(full_coef_aunt_x)),rep("Cousin",length(full_coef_cousin_x)),rep("Other",length(full_coef_other_x)),rep("Grandmother",length(full_coef_grandmother_x)),
									   rep("Mother",length(full_coef_mother_y)),rep("Sister",length(full_coef_sister_y)),rep("Aunt",length(full_coef_aunt_y)),rep("Cousin",length(full_coef_cousin_y)),rep("Other",length(full_coef_other_y)),rep("Grandmother",length(full_coef_grandmother_y)),
									   rep("Mother",length(full_coef_mother_z)),rep("Sister",length(full_coef_sister_z)),rep("Aunt",length(full_coef_aunt_z)),rep("Cousin",length(full_coef_cousin_z)),rep("Other",length(full_coef_other_z)),rep("Grandmother",length(full_coef_grandmother_z))),
									   level=c("Other","Grandmother","Aunt","Cousin","Sister","Mother")))


#

ggplot(dataTemp,aes(x=Coef,y=Relative,fill=Coo))+
geom_boxplot(show.legend=FALSE,col="black",width=0.6,size=0.25,outlier.size=0.02,outlier.stroke=0.1,alpha=0.5)+
geom_vline(xintercept=0,lty=1,col="red",lwd=0.2)+
scale_x_continuous(breaks=c(-1,0,1))+
labs(x="Coefficients",y="",title="")+
theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.5))


ggsave("spatial_ODE_relative_coefficients.pdf",width=unit(2.25,"in"),height=unit(2.1,"in"))

# legend 

ggplot(dataTemp,aes(x=Coef,y=Relative,fill=Coo))+
geom_boxplot(show.legend=TRUE,col="black",width=0.6,size=0.25,outlier.size=0.02,outlier.stroke=0.1,alpha=0.5)+
geom_vline(xintercept=0,lty=1,col="red",lwd=0.2)+
scale_x_continuous(breaks=c(-1,0,1))+
labs(x="Coefficients",y="",title="")+
theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.5))


ggsave("spatial_ODE_relative_coefficients_legend.pdf",width=unit(2.25,"in"),height=unit(2.1,"in"))



### R2 comparison : full ; mother ; sister ; mother and sister

dataTemp<-data.frame(R2=c(R2_full_x[with_sister_ind_x],R2_mother_x[with_sister_ind_x],R2_sister_x[with_sister_ind_x],R2_mother_sister_x[with_sister_ind_x],
						  R2_full_y[with_sister_ind_x],R2_mother_y[with_sister_ind_x],R2_sister_y[with_sister_ind_x],R2_mother_sister_y[with_sister_ind_x],
						  R2_full_z[with_sister_ind_x],R2_mother_z[with_sister_ind_x],R2_sister_z[with_sister_ind_x],R2_mother_sister_z[with_sister_ind_x]),
					 Coo=factor(c(rep("x",length(with_sister_ind_x)*4),rep("y",length(with_sister_ind_x)*4),rep("z",length(with_sister_ind_x)*4)),level=c("x","y","z")),
					 Model=factor(c(rep("Full",length(with_sister_ind_x)),rep("Mother",length(with_sister_ind_x)),rep("Sister",length(with_sister_ind_x)),rep("Mother and sister",length(with_sister_ind_x)),
								    rep("Full",length(with_sister_ind_x)),rep("Mother",length(with_sister_ind_x)),rep("Sister",length(with_sister_ind_x)),rep("Mother and sister",length(with_sister_ind_x)),
									rep("Full",length(with_sister_ind_x)),rep("Mother",length(with_sister_ind_x)),rep("Sister",length(with_sister_ind_x)),rep("Mother and sister",length(with_sister_ind_x))),
					 level=c("Mother and sister","Sister","Mother","Full")))

#

ggplot(dataTemp,aes(x=R2,y=Model,fill=Coo))+
geom_boxplot(show.legend=FALSE,width=0.7,col="black",outlier.shape=NA,size=0.2,alpha=0.5)+
scale_x_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.5))+
labs(x=expression(paste("Model performance (",italic(R)^2,")")),y="",title="")


ggsave("spatial_ODE_R2_mother_sister.pdf",width=unit(2.2,"in"),height=unit(2.1,"in"))




### R2 comparison : full ; no mother ; no sister

dataTemp<-data.frame(R2=c(R2_full_x[with_sister_ind_x],R2_without_mother_x[with_sister_ind_x],R2_without_sister_x[with_sister_ind_x],
						  R2_full_y[with_sister_ind_x],R2_without_mother_y[with_sister_ind_x],R2_without_sister_y[with_sister_ind_x],
						  R2_full_z[with_sister_ind_x],R2_without_mother_z[with_sister_ind_x],R2_without_sister_z[with_sister_ind_x]),
					 Coo=factor(c(rep("x",length(with_sister_ind_x)*3),rep("y",length(with_sister_ind_x)*3),rep("z",length(with_sister_ind_x)*3)),level=c("x","y","z")),
					 Model=factor(c(rep("Full",length(with_sister_ind_x)),rep("No mother",length(with_sister_ind_x)),rep("No sister",length(with_sister_ind_x)),
								    rep("Full",length(with_sister_ind_x)),rep("No mother",length(with_sister_ind_x)),rep("No sister",length(with_sister_ind_x)),
									rep("Full",length(with_sister_ind_x)),rep("No mother",length(with_sister_ind_x)),rep("No sister",length(with_sister_ind_x))),
					 level=c("No sister","No mother","Full")))


#

ggplot(dataTemp,aes(x=R2,y=Model,fill=Coo))+
geom_boxplot(show.legend=FALSE,width=0.7,col="black",outlier.shape=NA,size=0.2,alpha=0.5)+
scale_x_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.5))+
labs(x=expression(paste("Model performance (",italic(R)^2,")")),y="",title="")


ggsave("spatial_ODE_R2_without_mother_sister.pdf",width=unit(1.8,"in"),height=unit(1.6,"in"))




### eigenvalues of full model : xyz


Re_full_x<-Re(eigen(t(xyz_x_model_full_coef[-1,]))$values)
Re_full_y<-Re(eigen(t(xyz_y_model_full_coef[-1,]))$values)
Re_full_z<-Re(eigen(t(xyz_z_model_full_coef[-1,]))$values)


Im_full_x<-Im(eigen(t(xyz_x_model_full_coef[-1,]))$values)
Im_full_y<-Im(eigen(t(xyz_y_model_full_coef[-1,]))$values)
Im_full_z<-Im(eigen(t(xyz_z_model_full_coef[-1,]))$values)



dataTemp<-data.frame(Re=c(Re_full_x,Re_full_y,Re_full_z),
					 Im=c(Im_full_x,Im_full_y,Im_full_z),
					 Coo=factor(c(rep("x",length(Re_full_x)),
								  rep("y",length(Re_full_y)),
								  rep("z",length(Re_full_z))),level=c("x","y","z")))


ggplot(dataTemp,aes(x=Re,y=Im,col=Coo))+
geom_point(width=0.2,size=0.2,stroke=0.1,show.legend=FALSE)+
geom_vline(xintercept=0,col="black",lwd=0.2)+
geom_hline(yintercept=0,col="black",lwd=0.2)+
geom_vline(xintercept=max(dataTemp[,1]),lty=2,col="black",lwd=0.2)+
scale_x_continuous(limits=c(-1.3,1.3),breaks=c(-1,0,1))+
scale_y_continuous(breaks=c(-0.3,0,0.3))+
labs(x="Real",y="Imaginary",title=expression(paste("Eigenvalues of ",italic(K))))

ggsave("spatial_ODE_xyz_eigenvalue_full.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


#legend


ggplot(dataTemp,aes(x=Re,y=Im,col=Coo))+
geom_point(width=0.2,size=0.2,stroke=0.1)+
geom_vline(xintercept=0,col="black",lwd=0.2)+
geom_hline(yintercept=0,col="black",lwd=0.2)+
geom_vline(xintercept=max(dataTemp[,1]),lty=2,col="black",lwd=0.2)+
scale_x_continuous(limits=c(-1.3,1.3),breaks=c(-1,0,1))+
scale_y_continuous(breaks=c(-0.3,0,0.3))+
labs(x="Real",y="Imaginary",title=expression(paste("Eigenvalues of ",italic(K))))


ggsave("spatial_ODE_xyz_eigenvalue_full_legend.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))



### eigenvalue of no mother model : xyz

Re_no_mother_x<-Re(eigen(t(xyz_x_model_nested_without_mother_coef[-1,]))$values)
Re_no_mother_y<-Re(eigen(t(xyz_y_model_nested_without_mother_coef[-1,]))$values)
Re_no_mother_z<-Re(eigen(t(xyz_z_model_nested_without_mother_coef[-1,]))$values)

Im_no_mother_x<-Im(eigen(t(xyz_x_model_nested_without_mother_coef[-1,]))$values)
Im_no_mother_y<-Im(eigen(t(xyz_y_model_nested_without_mother_coef[-1,]))$values)
Im_no_mother_z<-Im(eigen(t(xyz_z_model_nested_without_mother_coef[-1,]))$values)

dataTemp<-data.frame(Re=c(Re_no_mother_x,Re_no_mother_y,Re_no_mother_z),
					 Im=c(Im_no_mother_x,Im_no_mother_y,Im_no_mother_z),
					 Coo=factor(c(rep("x",length(Re_no_mother_x)),
								  rep("y",length(Re_no_mother_y)),
								  rep("z",length(Re_no_mother_z))),level=c("x","y","z")))

ggplot(dataTemp,aes(x=Re,y=Im,col=Coo))+
geom_point(width=0.2,size=0.2,stroke=0.1,show.legend=FALSE)+
geom_vline(xintercept=0,col="black",lwd=0.2)+
geom_hline(yintercept=0,col="black",lwd=0.2)+
geom_vline(xintercept=max(dataTemp[,1]),lty=2,col="black",lwd=0.2)+
scale_x_continuous(limits=c(-1.3,1.3),breaks=c(-1,0,1))+
scale_y_continuous(breaks=c(-0.3,0,0.3))+
labs(x="Real",y="Imaginary",title="No mother")


ggsave("spatial_ODE_xyz_eigenvalue_no_mother.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


### eigenvalue of no sister model : xyz

Re_no_sister_x<-Re(eigen(t(xyz_x_model_nested_without_sister_coef[-1,]))$values)
Re_no_sister_y<-Re(eigen(t(xyz_y_model_nested_without_sister_coef[-1,]))$values)
Re_no_sister_z<-Re(eigen(t(xyz_z_model_nested_without_sister_coef[-1,]))$values)

Im_no_sister_x<-Im(eigen(t(xyz_x_model_nested_without_sister_coef[-1,]))$values)
Im_no_sister_y<-Im(eigen(t(xyz_y_model_nested_without_sister_coef[-1,]))$values)
Im_no_sister_z<-Im(eigen(t(xyz_z_model_nested_without_sister_coef[-1,]))$values)

dataTemp<-data.frame(Re=c(Re_no_sister_x,Re_no_sister_y,Re_no_sister_z),
					 Im=c(Im_no_sister_x,Im_no_sister_y,Im_no_sister_z),
					 Coo=factor(c(rep("x",length(Re_no_sister_x)),
								  rep("y",length(Re_no_sister_y)),
								  rep("z",length(Re_no_sister_z))),level=c("x","y","z")))

ggplot(dataTemp,aes(x=Re,y=Im,col=Coo))+
geom_point(width=0.2,size=0.2,stroke=0.1,show.legend=FALSE)+
geom_vline(xintercept=0,col="black",lwd=0.2)+
geom_hline(yintercept=0,col="black",lwd=0.2)+
geom_vline(xintercept=max(dataTemp[,1]),lty=2,col="black",lwd=0.2)+
scale_x_continuous(limits=c(-1.3,1.3),breaks=c(-1,0,1))+
scale_y_continuous(breaks=c(-0.3,0,0.3))+
labs(x="Real",y="Imaginary",title="No sister")


ggsave("spatial_ODE_xyz_eigenvalue_no_sister.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


### cell contact 

cell_contact<-read.csv(paste0(focal_projectFolderPath,"/data/annotation/effective_connect_list.csv"))

cell_contact_pair<-do.call(rbind,strsplit(cell_contact[,1],"->"))

sort(table(c(apply(cell_contact_pair,1,function(x){paste0(x[1],"_",x[2])}),apply(cell_contact_pair,1,function(x){paste0(x[2],"_",x[1])}))))


cell_contact_pair2<-do.call(c,lapply(unique(c(cell_contact_pair)),function(x){
paste0(x,"_",c(cell_contact_pair[which(cell_contact_pair[,2]==x),1],cell_contact_pair[which(cell_contact_pair[,1]==x),2]))}))



#x

x_full_coef_df<-do.call(rbind,lapply(1:dim(xyz_x_model_full_coef)[2],function(i){
	da_name<-colnames(xyz_x_model_full_coef)[i]
	mo_name<-rownames(xyz_x_model_full_coef)[i+1]
	focalInd<-setdiff(2:dim(xyz_x_model_full_coef)[1],setdiff(union(which(duplicated(rownames(xyz_x_model_full_coef))),
		union(which(rownames(xyz_x_model_full_coef)==da_name),which(rownames(xyz_x_model_full_coef)==mo_name))),i+1))
		
	focalCell<-unlist(strsplit(rownames(xyz_x_model_full_coef)[focalInd],"_x"))
	dataTemp<-data.frame(cell_pair=paste0(unlist(strsplit(colnames(xyz_x_model_full_coef)[i],"_x")),"_",focalCell),
		coef=xyz_x_model_full_coef[focalInd,i])
	
	dataTemp$contact<-rep(0,dim(dataTemp)[1])
	dataTemp$contact[na.omit(match(cell_contact_pair2,dataTemp$cell_pair))]=1
	dataTemp$contact<-factor(dataTemp$contact,level=c(1,0))

	dataTemp$mo_sis<-rep("other",dim(dataTemp)[1])
	dataTemp$mo_sis[na.omit(match(apply(cell_mo_da1_da2,1,function(x){paste0(x[2],"_",x[1])}),dataTemp$cell_pair))]<-"mo"
	dataTemp$mo_sis[na.omit(match(apply(cell_mo_da1_da2,1,function(x){paste0(x[2],"_",x[3])}),dataTemp$cell_pair))]<-"sis"
	
	
	rownames(dataTemp)<-NULL
	return(dataTemp)

}))


#y

y_full_coef_df<-do.call(rbind,lapply(1:dim(xyz_y_model_full_coef)[2],function(i){
	da_name<-colnames(xyz_y_model_full_coef)[i]
	mo_name<-rownames(xyz_y_model_full_coef)[i+1]
	focalInd<-setdiff(2:dim(xyz_y_model_full_coef)[1],setdiff(union(which(duplicated(rownames(xyz_y_model_full_coef))),
		union(which(rownames(xyz_y_model_full_coef)==da_name),which(rownames(xyz_y_model_full_coef)==mo_name))),i+1))
		
	focalCell<-unlist(strsplit(rownames(xyz_y_model_full_coef)[focalInd],"_y"))
	dataTemp<-data.frame(cell_pair=paste0(unlist(strsplit(colnames(xyz_y_model_full_coef)[i],"_y")),"_",focalCell),
		coef=xyz_y_model_full_coef[focalInd,i])
	
	dataTemp$contact<-rep(0,dim(dataTemp)[1])
	dataTemp$contact[na.omit(match(cell_contact_pair2,dataTemp$cell_pair))]=1
	dataTemp$contact<-factor(dataTemp$contact,level=c(1,0))

	dataTemp$mo_sis<-rep("other",dim(dataTemp)[1])
	dataTemp$mo_sis[na.omit(match(apply(cell_mo_da1_da2,1,function(x){paste0(x[2],"_",x[1])}),dataTemp$cell_pair))]<-"mo"
	dataTemp$mo_sis[na.omit(match(apply(cell_mo_da1_da2,1,function(x){paste0(x[2],"_",x[3])}),dataTemp$cell_pair))]<-"sis"
	
	
	rownames(dataTemp)<-NULL
	return(dataTemp)

}))


#z

z_full_coef_df<-do.call(rbind,lapply(1:dim(xyz_z_model_full_coef)[2],function(i){
	da_name<-colnames(xyz_z_model_full_coef)[i]
	mo_name<-rownames(xyz_z_model_full_coef)[i+1]
	focalInd<-setdiff(2:dim(xyz_z_model_full_coef)[1],setdiff(union(which(duplicated(rownames(xyz_z_model_full_coef))),
		union(which(rownames(xyz_z_model_full_coef)==da_name),which(rownames(xyz_z_model_full_coef)==mo_name))),i+1))
		
	focalCell<-unlist(strsplit(rownames(xyz_z_model_full_coef)[focalInd],"_z"))
	dataTemp<-data.frame(cell_pair=paste0(unlist(strsplit(colnames(xyz_z_model_full_coef)[i],"_z")),"_",focalCell),
		coef=xyz_z_model_full_coef[focalInd,i])
	
	dataTemp$contact<-rep(0,dim(dataTemp)[1])
	dataTemp$contact[na.omit(match(cell_contact_pair2,dataTemp$cell_pair))]=1
	dataTemp$contact<-factor(dataTemp$contact,level=c(1,0))

	dataTemp$mo_sis<-rep("other",dim(dataTemp)[1])
	dataTemp$mo_sis[na.omit(match(apply(cell_mo_da1_da2,1,function(x){paste0(x[2],"_",x[1])}),dataTemp$cell_pair))]<-"mo"
	dataTemp$mo_sis[na.omit(match(apply(cell_mo_da1_da2,1,function(x){paste0(x[2],"_",x[3])}),dataTemp$cell_pair))]<-"sis"
	
	
	rownames(dataTemp)<-NULL
	return(dataTemp)

}))



# proportion of non-zero coefficients 

temp=x_full_coef_df[x_full_coef_df$mo_sis=="other",]
c1=length(which(temp$coef!=0&temp$contact==0))/length(which(temp$contact==0))
c2=length(which(temp$coef!=0&temp$contact==1))/length(which(temp$contact==1))

c1_se=sqrt(c1*(1-c1)/length(which(temp$contact==0)))
c2_se=sqrt(c2*(1-c2)/length(which(temp$contact==1)))

temp=y_full_coef_df[y_full_coef_df$mo_sis=="other",]
c3=length(which(temp$coef!=0&temp$contact==0))/length(which(temp$contact==0))
c4=length(which(temp$coef!=0&temp$contact==1))/length(which(temp$contact==1))

c3_se=sqrt(c3*(1-c3)/length(which(temp$contact==0)))
c4_se=sqrt(c4*(1-c4)/length(which(temp$contact==1)))

temp=z_full_coef_df[z_full_coef_df$mo_sis=="other",]
c5=length(which(temp$coef!=0&temp$contact==0))/length(which(temp$contact==0))
c6=length(which(temp$coef!=0&temp$contact==1))/length(which(temp$contact==1))

c5_se=sqrt(c5*(1-c5)/length(which(temp$contact==0)))
c6_se=sqrt(c6*(1-c6)/length(which(temp$contact==1)))

#dataTemp<-data.frame(count=c(c2,c1,c4,c3,c6,c5),
#	contact=c(c(0,1),c(0.2,1.2),c(0.4,1.4)),
#	coo=c("x","x","y","y","z","z"))

#ggplot(dataTemp,aes(x=contact,y=count,col=coo,group=coo))+geom_point()+geom_line()+
#scale_x_continuous(breaks=c(0.2,1.2))+
#scale_y_continuous(limits=c(0.2,0.8),breaks=c(0.2,0.5,0.8))+
#labs(x="Cell type",y="Proportion of\nnon-zero coefficients",title="")+
#theme(legend.position=c(0.8,0.8))

#ggsave("spatial_full_contact_non_zero_count.pdf",width=unit(1.2,"in"),height=unit(1.7,"in"))

dataTemp<-data.frame(proportion=c(c2,c1,c4,c3,c6,c5),
	proportion_se=c(c2_se,c1_se,c4_se,c3_se,c6_se,c5_se)*2,
	type=rep(c("contact","no contact"),3),
	coo=c("x","x","y","y","z","z"))
	
ggplot(dataTemp,aes(x=type,y=proportion,col=coo,fill=coo,ymin=proportion-proportion_se,ymax=proportion+proportion_se))+
geom_bar(position=position_dodge2(padding=0.3),stat="identity",width=0.65,show.legend=FALSE)+
geom_errorbar(position=position_dodge2(padding=0.3),stat="identity",width=0.65,lwd=0.2,show.legend=FALSE)+
scale_y_continuous(breaks=c(0,0.3,0.6))+
labs(x="Cell type",y="Proportion of\nnon-zero coefficients",title="")

ggsave("spatial_full_contact_non_zero_proportion.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))

# legend

ggplot(dataTemp,aes(x=type,y=proportion,col=coo,fill=coo))+
geom_bar(position=position_dodge2(padding=0.3),stat="identity",width=0.65,show.legend=TRUE)+
scale_y_continuous(breaks=c(0,0.3,0.6))+
labs(x="Cell type",y="Proportion of\nnon-zero coefficients",title="")

ggsave("spatial_full_contact_non_zero_proportion_legend.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


# effect size of coefficients 

dataTemp<-rbind(x_full_coef_df[x_full_coef_df$mo_sis=="other"&x_full_coef_df$coef!=0,],
	y_full_coef_df[y_full_coef_df$mo_sis=="other"&y_full_coef_df$coef!=0,],
	z_full_coef_df[z_full_coef_df$mo_sis=="other"&z_full_coef_df$coef!=0,])

dataTemp$coo=c(rep("italic(K)[x]",length(which(x_full_coef_df$mo_sis=="other"&x_full_coef_df$coef!=0))),
	rep("italic(K)[y]",length(which(y_full_coef_df$mo_sis=="other"&y_full_coef_df$coef!=0))),
	rep("italic(K)[z]",length(which(z_full_coef_df$mo_sis=="other"&z_full_coef_df$coef!=0))))

dataTemp$contact<-factor(ifelse(dataTemp$contact==1,"contact","no contact"),level=c("contact","no contact"))
	
ggplot(dataTemp,aes(x=abs(coef),lty=contact,col=coo))+geom_density(show.legend=FALSE)+
scale_x_log10(limits=c(1e-5,1),breaks=c(1e-4,1e-2,1))+
scale_y_continuous(breaks=c(0,0.4,0.8))+
theme(legend.position=c(0.1,0.8),strip.background=element_blank())+
facet_wrap(vars(coo),labeller=label_parsed)+
labs(x=expression(paste("Absolute value of non-zero coefficients in ",italic(K))),y="Density")
	
ggsave("spatial_full_contact_non_zero_distribution_xyz.pdf",width=unit(3.4,"in"),height=unit(1.6,"in"))

	
	
max(abs(x_full_coef_df[x_full_coef_df$mo_sis=="other"&x_full_coef_df$coef!=0,]$coef))

ggplot(x_full_coef_df[x_full_coef_df$mo_sis=="other"&x_full_coef_df$coef!=0,],aes(lty=contact,x=abs(coef)))+geom_density()+
scale_x_log10(limits=c(1e-5,1),breaks=c(1e-4,1e-2,1))+
scale_y_continuous(breaks=c(0,0.4,0.8))+
theme(legend.position=c(1e-1,0.8))+
labs(x="Absolute coefficient",y="Density",title=expression(italic(K)[x]))

ggsave("spatial_full_contact_non_zero_distribution_x.pdf",width=unit(1.2,"in"),height=unit(1.7,"in"))


max(abs(y_full_coef_df[y_full_coef_df$mo_sis=="other"&y_full_coef_df$coef!=0,]$coef))

ggplot(y_full_coef_df[y_full_coef_df$mo_sis=="other"&y_full_coef_df$coef!=0,],aes(lty=contact,x=abs(coef)))+geom_density()+
scale_x_log10(limits=c(1e-5,1),breaks=c(1e-4,1e-2,1))+
scale_y_continuous(breaks=c(0,0.4,0.8))+
theme(legend.position=c(1e-1,0.8))+
labs(x="Absolute coefficient",y="Density",title=expression(italic(K)[y]))

ggsave("spatial_full_contact_non_zero_distribution_y.pdf",width=unit(1.2,"in"),height=unit(1.7,"in"))


max(abs(z_full_coef_df[z_full_coef_df$mo_sis=="other"&z_full_coef_df$coef!=0,]$coef))

ggplot(z_full_coef_df[z_full_coef_df$mo_sis=="other"&z_full_coef_df$coef!=0,],aes(lty=contact,x=abs(coef)))+geom_density()+
scale_x_log10(limits=c(1e-5,1),breaks=c(1e-4,1e-2,1))+
scale_y_continuous(breaks=c(0,0.4,0.8))+
theme(legend.position=c(1e-1,0.8))+
labs(x="Absolute coefficient",y="Density",title=expression(italic(K)[z]))

ggsave("spatial_full_contact_non_zero_distribution_z.pdf",width=unit(1.2,"in"),height=unit(1.7,"in"))


##########################################################################################

### block descendant 


## R2 compare with full model 

dataTemp<-data.frame(xyz_x_model_full_R2,xyz_y_model_full_R2,xyz_z_model_full_R2,
	xyz_x_model_block_descendant_R2,xyz_y_model_block_descendant_R2,xyz_z_model_block_descendant_R2)

#x 

ggplot(dataTemp,aes(x=xyz_x_model_full_R2,y=xyz_x_model_block_descendant_R2))+
geom_point(size=0.3,col=colors()[182])+
geom_abline(slope=1,intercept=0,lwd=0.3)+
labs(x="full model",y="descendant blocked",title=expression(italic(R)[x]^2))

ggsave("spatial_descendant_block_R2_compare_x.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

#y

ggplot(dataTemp,aes(x=xyz_y_model_full_R2,y=xyz_y_model_block_descendant_R2))+
geom_point(size=0.3,col=colors()[182])+
geom_abline(slope=1,intercept=0,lwd=0.3)+
labs(x="full model",y="descendant blocked",title=expression(italic(R)[y]^2))

ggsave("spatial_descendant_block_R2_compare_y.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

#z

ggplot(dataTemp,aes(x=xyz_z_model_full_R2,y=xyz_z_model_block_descendant_R2))+
geom_point(size=0.3,col=colors()[182])+
geom_abline(slope=1,intercept=0,lwd=0.3)+
labs(x="full model",y="descendant blocked",title=expression(italic(R)[z]^2))

ggsave("spatial_descendant_block_R2_compare_z.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


## eigenvalue of descendant locked 

#

Re_block_descendant_x<-Re(eigen(t(xyz_x_model_block_descendant_coef[-1,]))$values)
Re_block_descendant_y<-Re(eigen(t(xyz_y_model_block_descendant_coef[-1,]))$values)
Re_block_descendant_z<-Re(eigen(t(xyz_z_model_block_descendant_coef[-1,]))$values)


Im_block_descendant_x<-Im(eigen(t(xyz_x_model_block_descendant_coef[-1,]))$values)
Im_block_descendant_y<-Im(eigen(t(xyz_y_model_block_descendant_coef[-1,]))$values)
Im_block_descendant_z<-Im(eigen(t(xyz_z_model_block_descendant_coef[-1,]))$values)



dataTemp<-data.frame(Re=c(Re_block_descendant_x,Re_block_descendant_y,Re_block_descendant_z),
					 Im=c(Im_block_descendant_x,Im_block_descendant_y,Im_block_descendant_z),
					 Coo=factor(c(rep("x",length(Re_block_descendant_x)),
								  rep("y",length(Re_block_descendant_y)),
								  rep("z",length(Re_block_descendant_z))),level=c("x","y","z")))


ggplot(dataTemp,aes(x=Re,y=Im,col=Coo))+
geom_point(width=0.2,size=0.2,stroke=0.1,show.legend=FALSE)+
geom_vline(xintercept=0,col="black",lwd=0.2)+
geom_hline(yintercept=0,col="black",lwd=0.2)+
geom_vline(xintercept=max(dataTemp[,1]),lty=2,col="black",lwd=0.2)+
scale_x_continuous(limits=c(-1.3,1.3),breaks=c(-1,0,1))+
scale_y_continuous(breaks=c(-0.3,0,0.3))+
labs(x="Real",y="Imaginary",title=expression(paste(italic(K)," of descendant blocked")))

ggsave("spatial_ODE_xyz_eigenvalue_block_descendant.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### in one generation 

# g10

Re_M1M11_full_x<-Re(eigen(t(xyz_x_M1M11_full_coef[-1,]))$values)
Re_M1M11_full_y<-Re(eigen(t(xyz_y_M1M11_full_coef[-1,]))$values)
Re_M1M11_full_z<-Re(eigen(t(xyz_z_M1M11_full_coef[-1,]))$values)


Im_M1M11_full_x<-Im(eigen(t(xyz_x_M1M11_full_coef[-1,]))$values)
Im_M1M11_full_y<-Im(eigen(t(xyz_y_M1M11_full_coef[-1,]))$values)
Im_M1M11_full_z<-Im(eigen(t(xyz_z_M1M11_full_coef[-1,]))$values)



dataTemp<-data.frame(Re=c(Re_M1M11_full_x,Re_M1M11_full_y,Re_M1M11_full_z),
					 Im=c(Im_M1M11_full_x,Im_M1M11_full_y,Im_M1M11_full_z),
					 Coo=factor(c(rep("x",length(Re_M1M11_full_x)),
								  rep("y",length(Re_M1M11_full_y)),
								  rep("z",length(Re_M1M11_full_z))),level=c("x","y","z")))


ggplot(dataTemp,aes(x=Re,y=Im,col=Coo))+
geom_point(width=0.2,size=0.2,stroke=0.1,show.legend=FALSE)+
geom_vline(xintercept=0,col="black",lwd=0.2)+
geom_hline(yintercept=0,col="black",lwd=0.2)+
geom_vline(xintercept=max(dataTemp[,1]),lty=2,col="black",lwd=0.2)+
scale_x_continuous(limits=c(-1.3,1.3),breaks=c(-1,0,1))+
scale_y_continuous(breaks=c(-0.3,0,0.3))+
labs(x="Real",y="Imaginary",title=expression(paste(italic(K)," in one generation")))

ggsave("spatial_ODE_xyz_eigenvalue_one_generation_g10.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



# g9

Re_M1M11g9_full_x<-Re(eigen(t(xyz_x_M1M11g9_full_coef[-1,]))$values)
Re_M1M11g9_full_y<-Re(eigen(t(xyz_y_M1M11g9_full_coef[-1,]))$values)
Re_M1M11g9_full_z<-Re(eigen(t(xyz_z_M1M11g9_full_coef[-1,]))$values)


Im_M1M11g9_full_x<-Im(eigen(t(xyz_x_M1M11g9_full_coef[-1,]))$values)
Im_M1M11g9_full_y<-Im(eigen(t(xyz_y_M1M11g9_full_coef[-1,]))$values)
Im_M1M11g9_full_z<-Im(eigen(t(xyz_z_M1M11g9_full_coef[-1,]))$values)



dataTemp<-data.frame(Re=c(Re_M1M11g9_full_x,Re_M1M11g9_full_y,Re_M1M11g9_full_z),
					 Im=c(Im_M1M11g9_full_x,Im_M1M11g9_full_y,Im_M1M11g9_full_z),
					 Coo=factor(c(rep("x",length(Re_M1M11g9_full_x)),
								  rep("y",length(Re_M1M11g9_full_y)),
								  rep("z",length(Re_M1M11g9_full_z))),level=c("x","y","z")))


ggplot(dataTemp,aes(x=Re,y=Im,col=Coo))+
geom_point(width=0.2,size=0.2,stroke=0.1,show.legend=FALSE)+
geom_vline(xintercept=0,col="black",lwd=0.2)+
geom_hline(yintercept=0,col="black",lwd=0.2)+
geom_vline(xintercept=max(dataTemp[,1]),lty=2,col="black",lwd=0.2)+
scale_x_continuous(limits=c(-1.3,1.3),breaks=c(-1,0,1))+
scale_y_continuous(breaks=c(-0.3,0,0.3))+
labs(x="Real",y="Imaginary",title=expression(paste(italic(K)," in one generation")))

ggsave("spatial_ODE_xyz_eigenvalue_one_generation_g9.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))




























