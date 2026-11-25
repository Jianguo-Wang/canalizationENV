


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


### translation 

translation<-function(focal_coo_list,focal_translation_mat){
	focal_coo_list_translation<-lapply(1:length(focal_coo_list),function(i){
	
		temp<-cbind(focal_coo_list[[i]][,1]-focal_translation_mat[i,1],
			focal_coo_list[[i]][,2]-focal_translation_mat[i,2],
			focal_coo_list[[i]][,3]-focal_translation_mat[i,3])
			
		colnames(temp)<-c("X","Y","Z")
		rownames(temp)<-rownames(focal_coo_list[[i]])
		
		return(temp)
	})
	
	names(focal_coo_list_translation)<-names(focal_coo_list)
	return(focal_coo_list_translation)
	
}


### rotation 

rotation<-function(focal_coo_list_translation,focal_rotation_mat){
	focal_coo_list_translation_rotation<-lapply(1:length(focal_coo_list_translation),function(i){
		temp<-focal_coo_list_translation[[i]]%*%focal_rotation_mat[[i]]
		colnames(temp)<-c("X","Y","Z")
		rownames(temp)<-rownames(focal_coo_list_translation[[i]])
		return(temp)
	})
	
	names(focal_coo_list_translation_rotation)<-names(focal_coo_list_translation)
	
	return(focal_coo_list_translation_rotation)
	
}

### scaling  

scaling<-function(focal_coo_list_translation_rotation,focal_embryo_size){
	focal_coo_list_translation_rotation_scaling<-lapply(1:length(focal_coo_list_translation_rotation),function(i){
		temp<-cbind(focal_coo_list_translation_rotation[[i]][,1]/focal_embryo_size[i,1],
			focal_coo_list_translation_rotation[[i]][,2]/focal_embryo_size[i,2],
			focal_coo_list_translation_rotation[[i]][,3]/focal_embryo_size[i,3])
		
		colnames(temp)<-c("X","Y","Z")
		rownames(temp)<-rownames(focal_coo_list_translation_rotation[[i]])
		return(temp)
	})
	
	names(focal_coo_list_translation_rotation_scaling)<-names(focal_coo_list_translation_rotation)
	
	return(focal_coo_list_translation_rotation_scaling)
	
}

### coo to xyz  

coo2xyz<-function(focal_coo){
	focal_coo2xyz<-t(sapply(1:length(focal_coo),function(i){c(focal_coo[[i]][,1],focal_coo[[i]][,2],focal_coo[[i]][,3])}))
	colnames(focal_coo2xyz)<-c(paste0(rownames(focal_coo[[1]]),"_x"),paste0(rownames(focal_coo[[1]]),"_y"),paste0(rownames(focal_coo[[1]]),"_z"))
	rownames(focal_coo2xyz)<-names(focal_coo)
	
	return(focal_coo2xyz)
	
}


####################################################################################################
#							for initial stage
####################################################################################################


focal_projectFolderPath<-r"(E:\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/data/alignment"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/data/interpolation/coo_list_q1_cell350.RData",
				"F:/wangj/Desktop/temp/canalizationENV/data/interpolation/coo_list_q2_cell350.RData",
				"F:/wangj/Desktop/temp/canalizationENV/data/interpolation/coo_list_q3_cell350.RData",
				"F:/wangj/Desktop/temp/canalizationENV/data/interpolation/coo_list_q4_cell350.RData",	
				"F:/wangj/Desktop/temp/canalizationENV/data/interpolation/coo_list_q5_cell350.RData",
				"F:/wangj/Desktop/temp/canalizationENV/data/interpolation/coo_list_q6_cell350.RData",
				"F:/wangj/Desktop/temp/canalizationENV/data/interpolation/coo_list_q7_cell350.RData",
				"F:/wangj/Desktop/temp/canalizationENV/data/interpolation/coo_list_q8_cell350.RData",
				"F:/wangj/Desktop/temp/canalizationENV/data/interpolation/coo_list_q9_cell350.RData",
				"F:/wangj/Desktop/temp/canalizationENV/data/interpolation/coo_list_q10_cell350.RData",
				"F:/wangj/Desktop/temp/canalizationENV/data/interpolation/coo_list_q11_cell350.RData"					
	),projectFolderPath=focal_projectFolderPath)


### translation matrix 

translation_mat<-t(sapply(coo_list_q1_cell350,function(x){
	return(apply(x,2,mean))
}))

save(translation_mat,file="translation_mat.RData")


# translation 

coo_list_q1_cell350_translation<-translation(coo_list_q1_cell350,translation_mat)

save(coo_list_q1_cell350_translation,file="coo_list_q1_cell350_translation.RData")


### rotation of template 

x2_y2_z2_sum_wt<-sapply(coo_list_q1_cell350_translation[1:105],function(x){(max(x[,1])-min(x[,1]))^2+(max(x[,2])-min(x[,2]))^2+(max(x[,3])-min(x[,3]))^2})

max_emb_ind<-which(x2_y2_z2_sum_wt==max(x2_y2_z2_sum_wt))

template=coo_list_q1_cell350_translation[[max_emb_ind]]

template_pca<-princomp(template)
template_scores<-template_pca$scores
template_scores[,2]<--template_scores[,2]
template_scores[,3]<--template_scores[,3]


# rotation matrix 

rotation_mat_list<-lapply(1:length(coo_list_q1_cell350_translation),function(i){
	print(i)
	
	focal_embryo<-coo_list_q1_cell350_translation[[i]]
	focal_embryo_pca<-princomp(focal_embryo)
	focal_embryo_scores<-focal_embryo_pca$scores
	focal_embryo_loadings<-focal_embryo_pca$loadings
	
	
	cor_temp<-cor(focal_embryo_scores,template_scores)
	
	ind1<-which.max(abs(cor_temp[,1]))
	ind2<-which.max(abs(cor_temp[,2]))
	ind3<-which.max(abs(cor_temp[,3]))
	

	if(ind1!=ind2&ind1!=ind3&ind2!=ind3){
		focal_embryo_loadings_refine<-cbind(focal_embryo_loadings[,ind1]*sign(cor(focal_embryo_scores[,ind1],template_scores[,1])),
			focal_embryo_loadings[,ind2]*sign(cor(focal_embryo_scores[,ind2],template_scores[,2])),
			focal_embryo_loadings[,ind3]*sign(cor(focal_embryo_scores[,ind3],template_scores[,3])))
			
	}else{
		focal_embryo_loadings_refine<-cbind(focal_embryo_loadings[,1]*sign(cor(focal_embryo_scores[,1],template_scores[,1])),
			focal_embryo_loadings[,2]*sign(cor(focal_embryo_scores[,2],template_scores[,2])),
			focal_embryo_loadings[,3]*sign(cor(focal_embryo_scores[,3],template_scores[,3])))
			
	}
	
	rownames(focal_embryo_loadings_refine)<-c("X","Y","Z")
	colnames(focal_embryo_loadings_refine)<-c("Comp.1","Comp.2","Comp.3")
	
	return(focal_embryo_loadings_refine)
})

names(rotation_mat_list)<-names(coo_list_q1_cell350_translation)

save(rotation_mat_list,file="rotation_mat_list.RData")


# rotation

coo_list_q1_cell350_translation_rotation<-rotation(coo_list_q1_cell350_translation,rotation_mat_list)

save(coo_list_q1_cell350_translation_rotation,file="coo_list_q1_cell350_translation_rotation.RData")


### scaling factor


embryo_size<-t(sapply(coo_list_q1_cell350_translation_rotation,function(x){
	x_pca<-princomp(cbind(coo_list_q1_cell350_translation_rotation[[58]][,1],x[,1]))
	x_slope<-x_pca$loadings[2,1]/x_pca$loadings[1,1]
	
	y_pca<-princomp(cbind(coo_list_q1_cell350_translation_rotation[[58]][,2],x[,2]))
	y_slope<-y_pca$loadings[2,1]/y_pca$loadings[1,1]	
	
	z_pca<-princomp(cbind(coo_list_q1_cell350_translation_rotation[[58]][,3],x[,3]))
	z_slope<-z_pca$loadings[2,1]/z_pca$loadings[1,1]	
	
	return(c(x_slope,y_slope,z_slope))
}))

colnames(embryo_size)<-c("x_size","y_size","z_size")

save(embryo_size,file="embryo_size.RData")


# scaling factor by wtMean 

coo_wtMean<-cbind(X=apply(do.call(cbind,lapply(coo_list_q1_cell350_translation_rotation[1:105],function(x){x[,1]})),1,mean),
	Y=apply(do.call(cbind,lapply(coo_list_q1_cell350_translation_rotation[1:105],function(x){x[,2]})),1,mean),
	Z=apply(do.call(cbind,lapply(coo_list_q1_cell350_translation_rotation[1:105],function(x){x[,3]})),1,mean))


embryo_size_wtMean<-t(sapply(coo_list_q1_cell350_translation_rotation,function(x){
	x_pca<-princomp(cbind(coo_wtMean[,1],x[,1]))
	x_slope<-x_pca$loadings[2,1]/x_pca$loadings[1,1]
	
	y_pca<-princomp(cbind(coo_wtMean[,2],x[,2]))
	y_slope<-y_pca$loadings[2,1]/y_pca$loadings[1,1]	
	
	z_pca<-princomp(cbind(coo_wtMean[,3],x[,3]))
	z_slope<-z_pca$loadings[2,1]/z_pca$loadings[1,1]	
	
	return(c(x_slope,y_slope,z_slope))
}))

colnames(embryo_size_wtMean)<-c("x_size","y_size","z_size")

save(embryo_size_wtMean,file="embryo_size_wtMean.RData")



# scaling 

coo_list_q1_cell350_translation_rotation_scaling<-scaling(coo_list_q1_cell350_translation_rotation,embryo_size)

save(coo_list_q1_cell350_translation_rotation_scaling,file="coo_list_q1_cell350_translation_rotation_scaling.RData")



### obtain noise 

q1_cell350_translation_rotation_xyz<-coo2xyz(coo_list_q1_cell350_translation_rotation)
q1_cell350_alignment_xyz<-coo2xyz(coo_list_q1_cell350_translation_rotation_scaling)
q1_cell350_alignment_xyz_noise<-scale(q1_cell350_alignment_xyz,center=T,scale=F)


save(q1_cell350_translation_rotation_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q1_cell350_translation_rotation_xyz.RData"))
save(q1_cell350_alignment_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q1_cell350_alignment_xyz.RData"))
save(q1_cell350_alignment_xyz_noise,file=paste0(focal_projectFolderPath,"/data/basic/q1_cell350_alignment_xyz_noise.RData"))

q1_cell350_alignment_xyz_mean=apply(q1_cell350_alignment_xyz,2,mean,na.rm=T)
save(q1_cell350_alignment_xyz_mean,file=paste0(focal_projectFolderPath,"/data/basic/q1_cell350_alignment_xyz_mean.RData"))


####################################################################################################
#							for final stage
####################################################################################################


# translation 

coo_list_q11_cell350_translation<-translation(coo_list_q11_cell350,translation_mat)

save(coo_list_q11_cell350_translation,file="coo_list_q11_cell350_translation.RData")

	
# rotation

coo_list_q11_cell350_translation_rotation<-rotation(coo_list_q11_cell350_translation,rotation_mat_list)

save(coo_list_q11_cell350_translation_rotation,file="coo_list_q11_cell350_translation_rotation.RData")


# scaling 

coo_list_q11_cell350_translation_rotation_scaling<-scaling(coo_list_q11_cell350_translation_rotation,embryo_size)

save(coo_list_q11_cell350_translation_rotation_scaling,file="coo_list_q11_cell350_translation_rotation_scaling.RData")

# obtain noise 

q11_cell350_translation_rotation_xyz<-coo2xyz(coo_list_q11_cell350_translation_rotation)

q11_cell350_alignment_xyz<-coo2xyz(coo_list_q11_cell350_translation_rotation_scaling)
q11_cell350_alignment_xyz_noise<-scale(q11_cell350_alignment_xyz,center=T,scale=F)

save(q11_cell350_translation_rotation_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q11_cell350_translation_rotation_xyz.RData"))
save(q11_cell350_alignment_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q11_cell350_alignment_xyz.RData"))
save(q11_cell350_alignment_xyz_noise,file=paste0(focal_projectFolderPath,"/data/basic/q11_cell350_alignment_xyz_noise.RData"))


q11_cell350_alignment_xyz_mean=apply(q11_cell350_alignment_xyz,2,mean,na.rm=T)
save(q11_cell350_alignment_xyz_mean,file=paste0(focal_projectFolderPath,"/data/basic/q11_cell350_alignment_xyz_mean.RData"))



####################################################################################################
#							for 2 stage
####################################################################################################


# translation 

coo_list_q2_cell350_translation<-translation(coo_list_q2_cell350,translation_mat)

save(coo_list_q2_cell350_translation,file="coo_list_q2_cell350_translation.RData")

	
# rotation

coo_list_q2_cell350_translation_rotation<-rotation(coo_list_q2_cell350_translation,rotation_mat_list)

save(coo_list_q2_cell350_translation_rotation,file="coo_list_q2_cell350_translation_rotation.RData")


# scaling 

coo_list_q2_cell350_translation_rotation_scaling<-scaling(coo_list_q2_cell350_translation_rotation,embryo_size)

save(coo_list_q2_cell350_translation_rotation_scaling,file="coo_list_q2_cell350_translation_rotation_scaling.RData")

# obtain noise 

q2_cell350_translation_rotation_xyz<-coo2xyz(coo_list_q2_cell350_translation_rotation)

q2_cell350_alignment_xyz<-coo2xyz(coo_list_q2_cell350_translation_rotation_scaling)
q2_cell350_alignment_xyz_noise<-scale(q2_cell350_alignment_xyz,center=T,scale=F)

save(q2_cell350_translation_rotation_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q2_cell350_translation_rotation_xyz.RData"))
save(q2_cell350_alignment_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q2_cell350_alignment_xyz.RData"))
save(q2_cell350_alignment_xyz_noise,file=paste0(focal_projectFolderPath,"/data/basic/q2_cell350_alignment_xyz_noise.RData"))


####################################################################################################
#							for 3 stage
####################################################################################################


# translation 

coo_list_q3_cell350_translation<-translation(coo_list_q3_cell350,translation_mat)

save(coo_list_q3_cell350_translation,file="coo_list_q3_cell350_translation.RData")

	
# rotation

coo_list_q3_cell350_translation_rotation<-rotation(coo_list_q3_cell350_translation,rotation_mat_list)

save(coo_list_q3_cell350_translation_rotation,file="coo_list_q3_cell350_translation_rotation.RData")


# scaling 

coo_list_q3_cell350_translation_rotation_scaling<-scaling(coo_list_q3_cell350_translation_rotation,embryo_size)

save(coo_list_q3_cell350_translation_rotation_scaling,file="coo_list_q3_cell350_translation_rotation_scaling.RData")

# obtain noise 

q3_cell350_translation_rotation_xyz<-coo2xyz(coo_list_q3_cell350_translation_rotation)

q3_cell350_alignment_xyz<-coo2xyz(coo_list_q3_cell350_translation_rotation_scaling)
q3_cell350_alignment_xyz_noise<-scale(q3_cell350_alignment_xyz,center=T,scale=F)

save(q3_cell350_translation_rotation_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q3_cell350_translation_rotation_xyz.RData"))
save(q3_cell350_alignment_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q3_cell350_alignment_xyz.RData"))
save(q3_cell350_alignment_xyz_noise,file=paste0(focal_projectFolderPath,"/data/basic/q3_cell350_alignment_xyz_noise.RData"))


####################################################################################################
#							for 4 stage
####################################################################################################


# translation 

coo_list_q4_cell350_translation<-translation(coo_list_q4_cell350,translation_mat)

save(coo_list_q4_cell350_translation,file="coo_list_q4_cell350_translation.RData")

	
# rotation

coo_list_q4_cell350_translation_rotation<-rotation(coo_list_q4_cell350_translation,rotation_mat_list)

save(coo_list_q4_cell350_translation_rotation,file="coo_list_q4_cell350_translation_rotation.RData")


# scaling 

coo_list_q4_cell350_translation_rotation_scaling<-scaling(coo_list_q4_cell350_translation_rotation,embryo_size)

save(coo_list_q4_cell350_translation_rotation_scaling,file="coo_list_q4_cell350_translation_rotation_scaling.RData")

# obtain noise 

q4_cell350_translation_rotation_xyz<-coo2xyz(coo_list_q4_cell350_translation_rotation)

q4_cell350_alignment_xyz<-coo2xyz(coo_list_q4_cell350_translation_rotation_scaling)
q4_cell350_alignment_xyz_noise<-scale(q4_cell350_alignment_xyz,center=T,scale=F)

save(q4_cell350_translation_rotation_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q4_cell350_translation_rotation_xyz.RData"))
save(q4_cell350_alignment_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q4_cell350_alignment_xyz.RData"))
save(q4_cell350_alignment_xyz_noise,file=paste0(focal_projectFolderPath,"/data/basic/q4_cell350_alignment_xyz_noise.RData"))


####################################################################################################
#							for 5 stage
####################################################################################################


# translation 

coo_list_q5_cell350_translation<-translation(coo_list_q5_cell350,translation_mat)

save(coo_list_q5_cell350_translation,file="coo_list_q5_cell350_translation.RData")

	
# rotation

coo_list_q5_cell350_translation_rotation<-rotation(coo_list_q5_cell350_translation,rotation_mat_list)

save(coo_list_q5_cell350_translation_rotation,file="coo_list_q5_cell350_translation_rotation.RData")


# scaling 

coo_list_q5_cell350_translation_rotation_scaling<-scaling(coo_list_q5_cell350_translation_rotation,embryo_size)

save(coo_list_q5_cell350_translation_rotation_scaling,file="coo_list_q5_cell350_translation_rotation_scaling.RData")

# obtain noise 

q5_cell350_translation_rotation_xyz<-coo2xyz(coo_list_q5_cell350_translation_rotation)

q5_cell350_alignment_xyz<-coo2xyz(coo_list_q5_cell350_translation_rotation_scaling)
q5_cell350_alignment_xyz_noise<-scale(q5_cell350_alignment_xyz,center=T,scale=F)

save(q5_cell350_translation_rotation_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q5_cell350_translation_rotation_xyz.RData"))
save(q5_cell350_alignment_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q5_cell350_alignment_xyz.RData"))
save(q5_cell350_alignment_xyz_noise,file=paste0(focal_projectFolderPath,"/data/basic/q5_cell350_alignment_xyz_noise.RData"))


####################################################################################################
#							for 6 stage
####################################################################################################


# translation 

coo_list_q6_cell350_translation<-translation(coo_list_q6_cell350,translation_mat)

save(coo_list_q6_cell350_translation,file="coo_list_q6_cell350_translation.RData")

	
# rotation

coo_list_q6_cell350_translation_rotation<-rotation(coo_list_q6_cell350_translation,rotation_mat_list)

save(coo_list_q6_cell350_translation_rotation,file="coo_list_q6_cell350_translation_rotation.RData")


# scaling 

coo_list_q6_cell350_translation_rotation_scaling<-scaling(coo_list_q6_cell350_translation_rotation,embryo_size)

save(coo_list_q6_cell350_translation_rotation_scaling,file="coo_list_q6_cell350_translation_rotation_scaling.RData")

# obtain noise 

q6_cell350_translation_rotation_xyz<-coo2xyz(coo_list_q6_cell350_translation_rotation)

q6_cell350_alignment_xyz<-coo2xyz(coo_list_q6_cell350_translation_rotation_scaling)
q6_cell350_alignment_xyz_noise<-scale(q6_cell350_alignment_xyz,center=T,scale=F)

save(q6_cell350_translation_rotation_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q6_cell350_translation_rotation_xyz.RData"))
save(q6_cell350_alignment_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q6_cell350_alignment_xyz.RData"))
save(q6_cell350_alignment_xyz_noise,file=paste0(focal_projectFolderPath,"/data/basic/q6_cell350_alignment_xyz_noise.RData"))


####################################################################################################
#							for 7 stage
####################################################################################################


# translation 

coo_list_q7_cell350_translation<-translation(coo_list_q7_cell350,translation_mat)

save(coo_list_q7_cell350_translation,file="coo_list_q7_cell350_translation.RData")

	
# rotation

coo_list_q7_cell350_translation_rotation<-rotation(coo_list_q7_cell350_translation,rotation_mat_list)

save(coo_list_q7_cell350_translation_rotation,file="coo_list_q7_cell350_translation_rotation.RData")


# scaling 

coo_list_q7_cell350_translation_rotation_scaling<-scaling(coo_list_q7_cell350_translation_rotation,embryo_size)

save(coo_list_q7_cell350_translation_rotation_scaling,file="coo_list_q7_cell350_translation_rotation_scaling.RData")

# obtain noise 

q7_cell350_translation_rotation_xyz<-coo2xyz(coo_list_q7_cell350_translation_rotation)

q7_cell350_alignment_xyz<-coo2xyz(coo_list_q7_cell350_translation_rotation_scaling)
q7_cell350_alignment_xyz_noise<-scale(q7_cell350_alignment_xyz,center=T,scale=F)

save(q7_cell350_translation_rotation_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q7_cell350_translation_rotation_xyz.RData"))
save(q7_cell350_alignment_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q7_cell350_alignment_xyz.RData"))
save(q7_cell350_alignment_xyz_noise,file=paste0(focal_projectFolderPath,"/data/basic/q7_cell350_alignment_xyz_noise.RData"))


####################################################################################################
#							for 8 stage
####################################################################################################


# translation 

coo_list_q8_cell350_translation<-translation(coo_list_q8_cell350,translation_mat)

save(coo_list_q8_cell350_translation,file="coo_list_q8_cell350_translation.RData")

	
# rotation

coo_list_q8_cell350_translation_rotation<-rotation(coo_list_q8_cell350_translation,rotation_mat_list)

save(coo_list_q8_cell350_translation_rotation,file="coo_list_q8_cell350_translation_rotation.RData")


# scaling 

coo_list_q8_cell350_translation_rotation_scaling<-scaling(coo_list_q8_cell350_translation_rotation,embryo_size)

save(coo_list_q8_cell350_translation_rotation_scaling,file="coo_list_q8_cell350_translation_rotation_scaling.RData")

# obtain noise 

q8_cell350_translation_rotation_xyz<-coo2xyz(coo_list_q8_cell350_translation_rotation)

q8_cell350_alignment_xyz<-coo2xyz(coo_list_q8_cell350_translation_rotation_scaling)
q8_cell350_alignment_xyz_noise<-scale(q8_cell350_alignment_xyz,center=T,scale=F)

save(q8_cell350_translation_rotation_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q8_cell350_translation_rotation_xyz.RData"))
save(q8_cell350_alignment_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q8_cell350_alignment_xyz.RData"))
save(q8_cell350_alignment_xyz_noise,file=paste0(focal_projectFolderPath,"/data/basic/q8_cell350_alignment_xyz_noise.RData"))


####################################################################################################
#							for 9 stage
####################################################################################################


# translation 

coo_list_q9_cell350_translation<-translation(coo_list_q9_cell350,translation_mat)

save(coo_list_q9_cell350_translation,file="coo_list_q9_cell350_translation.RData")

	
# rotation

coo_list_q9_cell350_translation_rotation<-rotation(coo_list_q9_cell350_translation,rotation_mat_list)

save(coo_list_q9_cell350_translation_rotation,file="coo_list_q9_cell350_translation_rotation.RData")


# scaling 

coo_list_q9_cell350_translation_rotation_scaling<-scaling(coo_list_q9_cell350_translation_rotation,embryo_size)

save(coo_list_q9_cell350_translation_rotation_scaling,file="coo_list_q9_cell350_translation_rotation_scaling.RData")

# obtain noise 

q9_cell350_translation_rotation_xyz<-coo2xyz(coo_list_q9_cell350_translation_rotation)

q9_cell350_alignment_xyz<-coo2xyz(coo_list_q9_cell350_translation_rotation_scaling)
q9_cell350_alignment_xyz_noise<-scale(q9_cell350_alignment_xyz,center=T,scale=F)

save(q9_cell350_translation_rotation_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q9_cell350_translation_rotation_xyz.RData"))
save(q9_cell350_alignment_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q9_cell350_alignment_xyz.RData"))
save(q9_cell350_alignment_xyz_noise,file=paste0(focal_projectFolderPath,"/data/basic/q9_cell350_alignment_xyz_noise.RData"))


####################################################################################################
#							for 10 stage
####################################################################################################


# translation 

coo_list_q10_cell350_translation<-translation(coo_list_q10_cell350,translation_mat)

save(coo_list_q10_cell350_translation,file="coo_list_q10_cell350_translation.RData")

	
# rotation

coo_list_q10_cell350_translation_rotation<-rotation(coo_list_q10_cell350_translation,rotation_mat_list)

save(coo_list_q10_cell350_translation_rotation,file="coo_list_q10_cell350_translation_rotation.RData")


# scaling 

coo_list_q10_cell350_translation_rotation_scaling<-scaling(coo_list_q10_cell350_translation_rotation,embryo_size)

save(coo_list_q10_cell350_translation_rotation_scaling,file="coo_list_q10_cell350_translation_rotation_scaling.RData")

# obtain noise 

q10_cell350_translation_rotation_xyz<-coo2xyz(coo_list_q10_cell350_translation_rotation)

q10_cell350_alignment_xyz<-coo2xyz(coo_list_q10_cell350_translation_rotation_scaling)
q10_cell350_alignment_xyz_noise<-scale(q10_cell350_alignment_xyz,center=T,scale=F)

save(q10_cell350_translation_rotation_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q10_cell350_translation_rotation_xyz.RData"))
save(q10_cell350_alignment_xyz,file=paste0(focal_projectFolderPath,"/data/basic/q10_cell350_alignment_xyz.RData"))
save(q10_cell350_alignment_xyz_noise,file=paste0(focal_projectFolderPath,"/data/basic/q10_cell350_alignment_xyz_noise.RData"))




#############################################################################################################################
									END
#############################################################################################################################






















