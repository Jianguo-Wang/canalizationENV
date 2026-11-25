

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


##########################################################################
                           reverse transform into random raw data 
##########################################################################

### functions 

#

transform_scaling_reverse<-function(focal_data_x,focal_data_y,focal_data_z){
	temp_x<-focal_data_x
	temp_y<-focal_data_y
	temp_z<-focal_data_z

	for(i in 1:dim(focal_data_x)[1]){
		temp_x[i,]<-focal_data_x[i,]*embryo_size[i,1]
	}
	
	for(i in 1:dim(focal_data_y)[1]){
		temp_y[i,]<-focal_data_y[i,]*embryo_size[i,2]
	}
	
	for(i in 1:dim(focal_data_z)[1]){
		temp_z[i,]<-focal_data_z[i,]*embryo_size[i,3]
	}
	
	return(list(temp_x,temp_y,temp_z))	
}

#

transform_mean_reverse<-function(focal_data_x,focal_data_y,focal_data_z){
	temp_x<-focal_data_x
	temp_y<-focal_data_y
	temp_z<-focal_data_z

	for(i in 1:dim(focal_data_x)[2]){
		temp_x[,i]<-focal_data_x[,i]+mean(q1_cell350_alignment_xyz_x[,i])
	}
	
	for(i in 1:dim(focal_data_y)[2]){
		temp_y[,i]<-focal_data_y[,i]+mean(q1_cell350_alignment_xyz_y[,i])
	}
	
	for(i in 1:dim(focal_data_z)[2]){
		temp_z[,i]<-focal_data_z[,i]+mean(q1_cell350_alignment_xyz_z[,i])
	}
	
	return(list(temp_x,temp_y,temp_z))
		
}

#

transform_rotate_reverse<-function(focal_data_x,focal_data_y,focal_data_z){
	temp_x<-focal_data_x
	temp_y<-focal_data_y
	temp_z<-focal_data_z

	for(i in 1:dim(focal_data_x)[1]){
		temp<-cbind(focal_data_x[i,],focal_data_y[i,],focal_data_z[i,])%*%solve(rotation_mat_list[[i]])
		temp_x[i,]<-temp[,1]
		temp_y[i,]<-temp[,2]
		temp_z[i,]<-temp[,3]	
	}
	
	return(list(temp_x,temp_y,temp_z))
}

#

transform_translate_reverse<-function(focal_data_x,focal_data_y,focal_data_z){
	temp_x<-focal_data_x
	temp_y<-focal_data_y
	temp_z<-focal_data_z
	
	
	for(i in 1:dim(focal_data_x)[1]){
		temp_x[i,]<-focal_data_x[i,]+translation_mat[i,1]
	}
	
	for(i in 1:dim(focal_data_y)[1]){
		temp_y[i,]<-focal_data_y[i,]+translation_mat[i,2]
	}
	
	for(i in 1:dim(focal_data_z)[1]){
		temp_z[i,]<-focal_data_z[i,]+translation_mat[i,3]
	}
	
	return(list(temp_x,temp_y,temp_z))
	
}


######################################################################
                      translation -> rotation -> scaling -> mean
######################################################################



#

transform_translate<-function(focal_data_x,focal_data_y,focal_data_z){
	temp_x<-focal_data_x
	temp_y<-focal_data_y
	temp_z<-focal_data_z
	
	
	for(i in 1:dim(focal_data_x)[1]){
		temp_x[i,]<-focal_data_x[i,]-mean(focal_data_x[i,])
	}
	
	for(i in 1:dim(focal_data_y)[1]){
		temp_y[i,]<-focal_data_y[i,]-mean(focal_data_y[i,])
	}
	
	for(i in 1:dim(focal_data_z)[1]){
		temp_z[i,]<-focal_data_z[i,]-mean(focal_data_z[i,])
	}
	
	return(list(temp_x,temp_y,temp_z))
	
}

#

transform_rotate<-function(focal_data_x,focal_data_y,focal_data_z,tempInd){
	temp_x<-focal_data_x
	temp_y<-focal_data_y
	temp_z<-focal_data_z
	
	template=cbind(focal_data_x[tempInd,],focal_data_y[tempInd,],focal_data_z[tempInd,])
	template_pca<-princomp(template)
	template_rot<-template_pca$scores
	
	for(i in 1:dim(focal_data_x)[1]){
		focal_data<-cbind(focal_data_x[i,],focal_data_y[i,],focal_data_z[i,])
		focal_pca<-princomp(focal_data)
		
		print(i)
		temp<-focal_pca$scores
		cor_temp<-cor(temp,template_rot)		
		
		ind1<-which.max(abs(cor_temp[,1]))
		ind2<-which.max(abs(cor_temp[,2]))
		ind3<-which.max(abs(cor_temp[,3]))
			
		if(ind1!=ind2&ind1!=ind3&ind2!=ind3){
			temp_x[i,]<-temp[,ind1]*sign(cor(temp[,ind1],template_rot[,1]))
			temp_y[i,]<-temp[,ind2]*sign(cor(temp[,ind2],template_rot[,2]))
			temp_z[i,]<-temp[,ind3]*sign(cor(temp[,ind3],template_rot[,3]))			
		}else{
			temp_x[i,]<-temp[,1]*sign(cor(temp[,1],template_rot[,1]))
			temp_y[i,]<-temp[,2]*sign(cor(temp[,2],template_rot[,2]))
			temp_z[i,]<-temp[,3]*sign(cor(temp[,3],template_rot[,3]))
		}

	}
	
	return(list(temp_x,temp_y,temp_z))
}

#

transform_scaling<-function(focal_data_x,focal_data_y,focal_data_z,tempInd){
	temp_x<-focal_data_x
	temp_y<-focal_data_y
	temp_z<-focal_data_z

	for(i in 1:dim(focal_data_x)[1]){
		x_pca<-princomp(cbind(focal_data_x[tempInd,],focal_data_x[i,]))
		x_slope<-x_pca$loadings[2,1]/x_pca$loadings[1,1]
		temp_x[i,]<-focal_data_x[i,]/x_slope
	}
	
	for(i in 1:dim(focal_data_y)[1]){
		y_pca<-princomp(cbind(focal_data_y[tempInd,],focal_data_y[i,]))
		y_slope<-y_pca$loadings[2,1]/y_pca$loadings[1,1]
		temp_y[i,]<-focal_data_y[i,]/y_slope
	}
	
	for(i in 1:dim(focal_data_z)[1]){
		z_pca<-princomp(cbind(focal_data_z[tempInd,],focal_data_z[i,]))
		z_slope<-z_pca$loadings[2,1]/z_pca$loadings[1,1]
		temp_z[i,]<-focal_data_z[i,]/z_slope
	}
	
	return(list(temp_x,temp_y,temp_z))	
}

#

transform_mean<-function(focal_data_x,focal_data_y,focal_data_z){
	temp_x<-scale(focal_data_x,scale=FALSE)
	temp_y<-scale(focal_data_y,scale=FALSE)
	temp_z<-scale(focal_data_z,scale=FALSE)

	return(list(temp_x,temp_y,temp_z))
}

#######################################################################

#######################################################################


library(ggplot2)


focal_projectFolderPath<-r"(F:\wangj\Desktop\temp\canalizationENV)"

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
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/translation_mat.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/rotation_mat_list.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_coef.RData"
),projectFolderPath=focal_projectFolderPath)




#

q1_cell350_alignment_xyz_noise_x<-q1_cell350_alignment_xyz_noise[,1:714]
q1_cell350_alignment_xyz_noise_y<-q1_cell350_alignment_xyz_noise[,(1+714):(714*2)]
q1_cell350_alignment_xyz_noise_z<-q1_cell350_alignment_xyz_noise[,(1+714*2):(714*3)]

q1_cell350_alignment_xyz_x<-q1_cell350_alignment_xyz[,1:714]
q1_cell350_alignment_xyz_y<-q1_cell350_alignment_xyz[,(1+714):(714*2)]
q1_cell350_alignment_xyz_z<-q1_cell350_alignment_xyz[,(1+714*2):(714*3)]

#

feature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise[,1:714]),"_x"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]


#
cell_name<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise[,1:714]),"_x"))
cell_gen<-sapply(match(cell_name,lineage_info),function(x){rowColRec(x,dim(lineage_info)[1],dim(lineage_info)[2])[1]})

################################################################################
                                x
################################################################################

X_mo<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da1,feature_names)]


###

focal_x<-q1_cell350_alignment_xyz_noise_x

for(i in 5:10){
	focal_da_s<-cell_name[cell_gen==i]
	for(j in 1:length(focal_da_s)){
		
		focal_da<-focal_da_s[j]
		focal_mo<-cell_mo_da1_da2[which(cell_mo_da1_da2[,2]==focal_da),1]
		
		focal_mo_da<-paste0(focal_mo,"_",focal_da)
		focal_mo_da_k<-xyz_x_model_mother_coef[focal_mo_da]
		
		focal_mo_ind<-match(focal_mo,cell_name)
		focal_da_ind<-match(focal_da,cell_name)
		
		set.seed(j)
		focal_x[,focal_da_ind]<-focal_x[,focal_mo_ind]*(focal_mo_da_k+1)+sample(residuals(lm(q1_cell350_alignment_xyz_noise_x[,focal_da_ind]~q1_cell350_alignment_xyz_noise_x[,focal_mo_ind])))
		
	}
	
}


plot(apply(q1_cell350_alignment_xyz_noise_x,2,sd),apply(focal_x,2,sd))
abline(0,1)

q1_cell350_alignment_xyz_noise_x_rand<-focal_x

#save(q1_cell350_alignment_xyz_noise_x_rand,file=paste0(focal_projectFolderPath,"/data/confounding/q1_cell350_alignment_xyz_noise_x_rand.RData"))

################################################################################
                                y
################################################################################

X_mo<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da1,feature_names)]


###

focal_y<-q1_cell350_alignment_xyz_noise_y

for(i in 5:10){
	focal_da_s<-cell_name[cell_gen==i]
	for(j in 1:length(focal_da_s)){
		
		focal_da<-focal_da_s[j]
		focal_mo<-cell_mo_da1_da2[which(cell_mo_da1_da2[,2]==focal_da),1]
		
		focal_mo_da<-paste0(focal_mo,"_",focal_da)
		focal_mo_da_k<-xyz_y_model_mother_coef[focal_mo_da]
		
		focal_mo_ind<-match(focal_mo,cell_name)
		focal_da_ind<-match(focal_da,cell_name)
		
		set.seed(j)
		focal_y[,focal_da_ind]<-focal_y[,focal_mo_ind]*(focal_mo_da_k+1)+sample(residuals(lm(q1_cell350_alignment_xyz_noise_y[,focal_da_ind]~q1_cell350_alignment_xyz_noise_y[,focal_mo_ind])))
		
	}
	
}


plot(apply(q1_cell350_alignment_xyz_noise_y,2,sd),apply(focal_y,2,sd))
abline(0,1)

q1_cell350_alignment_xyz_noise_y_rand<-focal_y

#save(q1_cell350_alignment_xyz_noise_y_rand,file=paste0(focal_projectFolderPath,"/data/confounding/q1_cell350_alignment_xyz_noise_y_rand.RData"))


################################################################################
                                z
################################################################################

X_mo<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da1,feature_names)]


###

focal_z<-q1_cell350_alignment_xyz_noise_z

for(i in 5:10){
	focal_da_s<-cell_name[cell_gen==i]
	for(j in 1:length(focal_da_s)){
		
		focal_da<-focal_da_s[j]
		focal_mo<-cell_mo_da1_da2[which(cell_mo_da1_da2[,2]==focal_da),1]
		
		focal_mo_da<-paste0(focal_mo,"_",focal_da)
		focal_mo_da_k<-xyz_z_model_mother_coef[focal_mo_da]
		
		focal_mo_ind<-match(focal_mo,cell_name)
		focal_da_ind<-match(focal_da,cell_name)
		
		set.seed(j)
		focal_z[,focal_da_ind]<-focal_z[,focal_mo_ind]*(focal_mo_da_k+1)+sample(residuals(lm(q1_cell350_alignment_xyz_noise_z[,focal_da_ind]~q1_cell350_alignment_xyz_noise_z[,focal_mo_ind])))
		
	}
	
}


plot(apply(q1_cell350_alignment_xyz_noise_z,2,sd),apply(focal_z,2,sd))
abline(0,1)

q1_cell350_alignment_xyz_noise_z_rand<-focal_z

#save(q1_cell350_alignment_xyz_noise_z_rand,file=paste0(focal_projectFolderPath,"/data/confounding/q1_cell350_alignment_xyz_noise_z_rand.RData"))

###################################################################################################################
							reverse process -> alignment
####################################################################################################################

###reverse

temp<-transform_mean_reverse(q1_cell350_alignment_xyz_noise_x_rand,q1_cell350_alignment_xyz_noise_y_rand,q1_cell350_alignment_xyz_noise_z_rand)
temp2<-transform_scaling_reverse(temp[[1]],temp[[2]],temp[[3]])
temp3<-transform_rotate_reverse(temp2[[1]],temp2[[2]],temp2[[3]])
temp4<-transform_translate_reverse(temp3[[1]],temp3[[2]],temp3[[3]])

rand_raw_x<-temp4[[1]]
rand_raw_y<-temp4[[2]]
rand_raw_z<-temp4[[3]]

#save(rand_raw_x,file=paste0(focal_projectFolderPath,"/data/confounding/rand_raw_x.RData"))
#save(rand_raw_y,file=paste0(focal_projectFolderPath,"/data/confounding/rand_raw_y.RData"))
#save(rand_raw_z,file=paste0(focal_projectFolderPath,"/data/confounding/rand_raw_z.RData"))

cor(cbind(temp3[[1]][10,],temp3[[2]][10,],temp3[[3]][10,]),cbind(temp2[[1]][10,],temp2[[2]][10,],temp2[[3]][10,]))
cor(cbind(temp4[[1]][10,],temp4[[2]][10,],temp4[[3]][10,]),cbind(temp3[[1]][10,],temp3[[2]][10,],temp3[[3]][10,]))

cor(cbind(rand_raw_x[1,],rand_raw_y[1,],rand_raw_z[1,]),cbind(rand_raw_x[2,],rand_raw_y[2,],rand_raw_z[2,]))


###alignment

temp<-transform_translate(rand_raw_x,rand_raw_y,rand_raw_z)
temp2<-transform_rotate(temp[[1]],temp[[2]],temp[[3]],tempInd=1)
temp3<-transform_scaling(temp2[[1]],temp2[[2]],temp2[[3]],tempInd=1)
temp4<-transform_mean(temp3[[1]],temp3[[2]],temp3[[3]])

rand_raw_x_process<-temp4[[1]]
rand_raw_y_process<-temp4[[2]]
rand_raw_z_process<-temp4[[3]]

###estimate k

#x 

X_mo<-rand_raw_x_process[,match(mo_sis_pair$mo,feature_names)]
X_da1<-rand_raw_x_process[,match(mo_sis_pair$da1,feature_names)]
delta_X<-X_da1-X_mo 
rand_k_x<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})

#y

X_mo<-rand_raw_y_process[,match(mo_sis_pair$mo,feature_names)]
X_da1<-rand_raw_y_process[,match(mo_sis_pair$da1,feature_names)]
delta_X<-X_da1-X_mo 
rand_k_y<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})

#z

X_mo<-rand_raw_z_process[,match(mo_sis_pair$mo,feature_names)]
X_da1<-rand_raw_z_process[,match(mo_sis_pair$da1,feature_names)]
delta_X<-X_da1-X_mo 
rand_k_z<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})

#save(rand_k_x,file=paste0(focal_projectFolderPath,"/data/confounding/rand_k_x.RData"))
#save(rand_k_y,file=paste0(focal_projectFolderPath,"/data/confounding/rand_k_y.RData"))
#save(rand_k_z,file=paste0(focal_projectFolderPath,"/data/confounding/rand_k_z.RData"))


###

par(mfrow=c(2,2))
plot(xyz_x_model_mother_coef,rand_k_x,main="R2=0.99")
cor(xyz_x_model_mother_coef,rand_k_x)
abline(0,1)
plot(xyz_y_model_mother_coef,rand_k_y,main="R2=0.99")
cor(xyz_y_model_mother_coef,rand_k_y)
abline(0,1)
plot(xyz_z_model_mother_coef,rand_k_z,main="R2=0.99")
cor(xyz_z_model_mother_coef,rand_k_z)
abline(0,1)



dataTemp<-data.frame(k_original=c(xyz_x_model_mother_coef,xyz_y_model_mother_coef,xyz_z_model_mother_coef),
	k_random=c(rand_k_x,rand_k_y,rand_k_z),
	coo=c(rep("kx",length(rand_k_x)),rep("ky",length(rand_k_y)),rep("kz",length(rand_k_z))))

ggplot(dataTemp,aes(x=k_original,y=k_random,col=coo,fill=coo))+
geom_point(alpha=0.5,size=0.1,stroke=0.1)+
geom_abline(slope=1,intercept=0,lty=2,lwd=0.3)+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])),breaks=c(-1,0))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])),breaks=c(-1,0))+
labs(x=expression(italic(k)[theoretical]),y=expression(italic(k)[estimated]),title="Simulated embryos")+
theme(legend.position=c(0.2,0.8))

ggsave("effect_of_alignment_in_simulated_embryos.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



###################################


