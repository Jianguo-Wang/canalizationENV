

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


focal_projectFolderPath<-r"(F:\wangj\Desktop\temp\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/figs"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/RowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData",		
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling_noise.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_duration_prcomp.RData",
	"E:/papers/canalizationENV/data/temporal/ccl_model_mother_coef.RData"
	),projectFolderPath=focal_projectFolderPath)


#########################################################################################

### 

feature_names<-unlist(strsplit(colnames(ccl_raw_nonNA_scaling_noise),"_ccl"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]

cell_8_mo<-intersect(na.omit(unique(c(lineage_info[8,]))),mo_sis_pair$mo)
cell_7_mo<-intersect(na.omit(unique(c(lineage_info[7,]))),mo_sis_pair$mo)
cell_6_mo<-intersect(na.omit(unique(c(lineage_info[6,]))),mo_sis_pair$mo)
cell_5_mo<-intersect(na.omit(unique(c(lineage_info[5,]))),mo_sis_pair$mo)
cell_4_mo<-intersect(na.omit(unique(c(lineage_info[4,]))),mo_sis_pair$mo)


##################################################################################################
                                           _CCL
##################################################################################################

num_of_embryo<-dim(ccl_raw_nonNA_scaling_noise)[1]


## 8

accuracy_8<-c()

cell_8_mo_ind<-unlist(lapply(cell_8_mo,function(x){which(mo_sis_pair$mo==x)}))


for(i in 1:50){

	focal_x<-ccl_raw_nonNA_scaling_noise

	print(i)
	
	set.seed(1)
	total_error<-t(sapply(sample(1:(num_of_embryo*length(cell_8_mo)),round(num_of_embryo*length(cell_8_mo)*0.1*i/50)),rowColRec,num_of_embryo,length(cell_8_mo)))
	error_mo<-lapply(1:num_of_embryo,function(ii){
		focalInd<-which(total_error[,1]==ii)
		if(length(focalInd)==0){
			return(NA)
		}else{
			return(cell_8_mo[total_error[focalInd,2]])
		}
	})
	
	for(j in 1:num_of_embryo){
		if(!is.na(error_mo[[j]][1])){
			focal_error_da1<-mo_sis_pair$da1[unlist(lapply(error_mo[[j]],function(x){which(mo_sis_pair$mo==x)}))]
			focal_error_da2<-mo_sis_pair$da2[unlist(lapply(error_mo[[j]],function(x){which(mo_sis_pair$mo==x)}))]
		
			focal_error_da1_da_ind<-unlist(lapply(focal_error_da1,function(x){
				temp<-which(mo_sis_pair$mo==x)
				if(length(temp)==0){
					temp<-c(NA,NA)
				}
				return(temp)
			}))
			focal_error_da2_da_ind<-unlist(lapply(focal_error_da2,function(x){
				temp<-which(mo_sis_pair$mo==x)
				if(length(temp)==0){
					temp<-c(NA,NA)
				}
				return(temp)
			}))
		

			focal_error_da1_ind<-c(match(focal_error_da1,feature_names),focal_error_da1_da_ind)
			focal_error_da2_ind<-c(match(focal_error_da2,feature_names),focal_error_da2_da_ind)
			
			focal_error_da1_da2_nonNAInd<-which(!is.na(focal_error_da1_ind)&!is.na(focal_error_da2_ind))
			focal_error_da1_da2_NAInd<-which(is.na(focal_error_da1_ind)|is.na(focal_error_da2_ind))
		
			temp<-focal_x[j,focal_error_da1_ind[focal_error_da1_da2_nonNAInd]]
			focal_x[j,focal_error_da1_ind[focal_error_da1_da2_nonNAInd]]<-focal_x[j,focal_error_da2_ind[focal_error_da1_da2_nonNAInd]]
			focal_x[j,focal_error_da2_ind[focal_error_da1_da2_nonNAInd]]<-temp
		
			focal_error_NA_ind<-na.omit(c(focal_error_da1_ind[focal_error_da1_da2_NAInd],focal_error_da2_ind[focal_error_da1_da2_NAInd]))
		
			focal_x[j,focal_error_NA_ind]<-NA	
		
		}	
	}
	
	X_mo<-focal_x[,match(mo_sis_pair$mo,feature_names)]
	X_da1<-focal_x[,match(mo_sis_pair$da1,feature_names)]
	delta_X<-X_da1-X_mo
	
	ccl_model_mother_coef_error<-ccl_model_mother_coef
	ccl_model_mother_coef_error[c(cell_8_mo_ind)]<-sapply(c(cell_8_mo_ind),function(iii){summary(lm(delta_X[,iii]~X_mo[,iii]))$coefficients[2,1]})
	
	accuracy_8[i]<-cor(ccl_model_mother_coef_error,ccl_model_mother_coef)
	
	plot(accuracy_8)
}




## 7

accuracy_7<-c()

cell_7_mo_ind<-unlist(lapply(cell_7_mo,function(x){which(mo_sis_pair$mo==x)}))


for(i in 1:50){

	focal_x<-ccl_raw_nonNA_scaling_noise

	print(i)
	
	set.seed(1)
	total_error<-t(sapply(sample(1:(num_of_embryo*length(cell_7_mo)),round(num_of_embryo*length(cell_7_mo)*0.1*i/50)),rowColRec,num_of_embryo,length(cell_7_mo)))
	error_mo<-lapply(1:num_of_embryo,function(ii){
		focalInd<-which(total_error[,1]==ii)
		if(length(focalInd)==0){
			return(NA)
		}else{
			return(cell_7_mo[total_error[focalInd,2]])
		}
	})
	
	for(j in 1:num_of_embryo){
		if(!is.na(error_mo[[j]][1])){
		
			focal_error_da1<-mo_sis_pair$da1[unlist(lapply(error_mo[[j]],function(x){which(mo_sis_pair$mo==x)}))]
			focal_error_da2<-mo_sis_pair$da2[unlist(lapply(error_mo[[j]],function(x){which(mo_sis_pair$mo==x)}))]
		
			focal_error_da1_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da1,function(x){which(mo_sis_pair$mo==x)}))]
			focal_error_da2_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da2,function(x){which(mo_sis_pair$mo==x)}))]		
		
			focal_error_da1_da_da_ind<-unlist(lapply(focal_error_da1_da,function(x){
				temp<-which(mo_sis_pair$mo==x)
				if(length(temp)==0){
					temp<-c(NA,NA)
				}
				return(temp)
			}))
			focal_error_da2_da_da_ind<-unlist(lapply(focal_error_da2_da,function(x){
				temp<-which(mo_sis_pair$mo==x)
				if(length(temp)==0){
					temp<-c(NA,NA)
				}
				return(temp)
			}))
		
		

			focal_error_da1_ind<-c(match(c(focal_error_da1,focal_error_da1_da),feature_names),focal_error_da1_da_da_ind)
			focal_error_da2_ind<-c(match(c(focal_error_da2,focal_error_da2_da),feature_names),focal_error_da2_da_da_ind)
		
			focal_error_da1_da2_nonNAInd<-which(!is.na(focal_error_da1_ind)&!is.na(focal_error_da2_ind))
			focal_error_da1_da2_NAInd<-which(is.na(focal_error_da1_ind)|is.na(focal_error_da2_ind))
		
			temp<-focal_x[j,focal_error_da1_ind[focal_error_da1_da2_nonNAInd]]
			focal_x[j,focal_error_da1_ind[focal_error_da1_da2_nonNAInd]]<-focal_x[j,focal_error_da2_ind[focal_error_da1_da2_nonNAInd]]
			focal_x[j,focal_error_da2_ind[focal_error_da1_da2_nonNAInd]]<-temp
		
			focal_error_NA_ind<-na.omit(c(focal_error_da1_ind[focal_error_da1_da2_NAInd],focal_error_da2_ind[focal_error_da1_da2_NAInd]))
		
			focal_x[j,focal_error_NA_ind]<-NA
		}
	}
	
	X_mo<-focal_x[,match(mo_sis_pair$mo,feature_names)]
	X_da1<-focal_x[,match(mo_sis_pair$da1,feature_names)]
	delta_X<-X_da1-X_mo
	
	ccl_model_mother_coef_error<-ccl_model_mother_coef
	ccl_model_mother_coef_error[c(cell_7_mo_ind,cell_8_mo_ind)]<-sapply(c(cell_7_mo_ind,cell_8_mo_ind),function(iii){summary(lm(delta_X[,iii]~X_mo[,iii]))$coefficients[2,1]})
	
	accuracy_7[i]<-cor(ccl_model_mother_coef_error,ccl_model_mother_coef)
	
	plot(accuracy_7)
}




## 6

accuracy_6<-c()

cell_6_mo_ind<-unlist(lapply(cell_6_mo,function(x){which(mo_sis_pair$mo==x)}))


for(i in 1:50){

	focal_x<-ccl_raw_nonNA_scaling_noise

	print(i)
	
	set.seed(1)
	total_error<-t(sapply(sample(1:(num_of_embryo*length(cell_6_mo)),round(num_of_embryo*length(cell_6_mo)*0.1*i/50)),rowColRec,num_of_embryo,length(cell_6_mo)))
	error_mo<-lapply(1:num_of_embryo,function(ii){
		focalInd<-which(total_error[,1]==ii)
		if(length(focalInd)==0){
			return(NA)
		}else{
			return(cell_6_mo[total_error[focalInd,2]])
		}
	})
	
	for(j in 1:num_of_embryo){
		if(!is.na(error_mo[[j]][1])){
		
			focal_error_da1<-mo_sis_pair$da1[unlist(lapply(error_mo[[j]],function(x){which(mo_sis_pair$mo==x)}))]
			focal_error_da2<-mo_sis_pair$da2[unlist(lapply(error_mo[[j]],function(x){which(mo_sis_pair$mo==x)}))]
		
			focal_error_da1_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da1,function(x){which(mo_sis_pair$mo==x)}))]
			focal_error_da2_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da2,function(x){which(mo_sis_pair$mo==x)}))]	
		
			focal_error_da1_da_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da1_da,function(x){which(mo_sis_pair$mo==x)}))]
			focal_error_da2_da_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da2_da,function(x){which(mo_sis_pair$mo==x)}))]	
		
			
			focal_error_da1_da_da_da_ind<-unlist(lapply(focal_error_da1_da_da,function(x){
				temp<-which(mo_sis_pair$mo==x)
				if(length(temp)==0){
					temp<-c(NA,NA)
				}
				return(temp)
			}))
			focal_error_da2_da_da_da_ind<-unlist(lapply(focal_error_da2_da_da,function(x){
				temp<-which(mo_sis_pair$mo==x)
				if(length(temp)==0){
					temp<-c(NA,NA)
				}
				return(temp)
			}))
		
		

			focal_error_da1_ind<-c(match(c(focal_error_da1,focal_error_da1_da,focal_error_da1_da_da),feature_names),focal_error_da1_da_da_da_ind)
			focal_error_da2_ind<-c(match(c(focal_error_da2,focal_error_da2_da,focal_error_da2_da_da),feature_names),focal_error_da2_da_da_da_ind)
		
			focal_error_da1_da2_nonNAInd<-which(!is.na(focal_error_da1_ind)&!is.na(focal_error_da2_ind))
			focal_error_da1_da2_NAInd<-which(is.na(focal_error_da1_ind)|is.na(focal_error_da2_ind))
		
			temp<-focal_x[j,focal_error_da1_ind[focal_error_da1_da2_nonNAInd]]
			focal_x[j,focal_error_da1_ind[focal_error_da1_da2_nonNAInd]]<-focal_x[j,focal_error_da2_ind[focal_error_da1_da2_nonNAInd]]
			focal_x[j,focal_error_da2_ind[focal_error_da1_da2_nonNAInd]]<-temp
		
			focal_error_NA_ind<-na.omit(c(focal_error_da1_ind[focal_error_da1_da2_NAInd],focal_error_da2_ind[focal_error_da1_da2_NAInd]))
		
			focal_x[j,focal_error_NA_ind]<-NA
		}
	}
	
	X_mo<-focal_x[,match(mo_sis_pair$mo,feature_names)]
	X_da1<-focal_x[,match(mo_sis_pair$da1,feature_names)]
	delta_X<-X_da1-X_mo
	
	ccl_model_mother_coef_error<-ccl_model_mother_coef
	ccl_model_mother_coef_error[c(cell_6_mo_ind,cell_7_mo_ind,cell_8_mo_ind)]<-sapply(c(cell_6_mo_ind,cell_7_mo_ind,cell_8_mo_ind),function(iii){summary(lm(delta_X[,iii]~X_mo[,iii]))$coefficients[2,1]})
	
	accuracy_6[i]<-cor(ccl_model_mother_coef_error,ccl_model_mother_coef)
	
	plot(accuracy_6)
}



#

accuracy_5<-c()

cell_5_mo_ind<-unlist(lapply(cell_5_mo,function(x){which(mo_sis_pair$mo==x)}))


for(i in 1:50){

	focal_x<-ccl_raw_nonNA_scaling_noise

	print(i)
	
	set.seed(1)
	total_error<-t(sapply(sample(1:(num_of_embryo*length(cell_5_mo)),round(num_of_embryo*length(cell_5_mo)*0.1*i/50)),rowColRec,num_of_embryo,length(cell_5_mo)))
	error_mo<-lapply(1:num_of_embryo,function(ii){
		focalInd<-which(total_error[,1]==ii)
		if(length(focalInd)==0){
			return(NA)
		}else{
			return(cell_5_mo[total_error[focalInd,2]])
		}
	})
	
	for(j in 1:num_of_embryo){
		if(!is.na(error_mo[[j]][1])){
		
			focal_error_da1<-mo_sis_pair$da1[unlist(lapply(error_mo[[j]],function(x){which(mo_sis_pair$mo==x)}))]
			focal_error_da2<-mo_sis_pair$da2[unlist(lapply(error_mo[[j]],function(x){which(mo_sis_pair$mo==x)}))]
		
			focal_error_da1_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da1,function(x){which(mo_sis_pair$mo==x)}))]
			focal_error_da2_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da2,function(x){which(mo_sis_pair$mo==x)}))]	
		
			focal_error_da1_da_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da1_da,function(x){which(mo_sis_pair$mo==x)}))]
			focal_error_da2_da_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da2_da,function(x){which(mo_sis_pair$mo==x)}))]	
		
			focal_error_da1_da_da_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da1_da_da,function(x){which(mo_sis_pair$mo==x)}))]
			focal_error_da2_da_da_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da2_da_da,function(x){which(mo_sis_pair$mo==x)}))]	
				
			focal_error_da1_da_da_da_da_ind<-unlist(lapply(focal_error_da1_da_da_da,function(x){
				temp<-which(mo_sis_pair$mo==x)
				if(length(temp)==0){
					temp<-c(NA,NA)
				}
				return(temp)
			}))
			focal_error_da2_da_da_da_da_ind<-unlist(lapply(focal_error_da2_da_da_da,function(x){
				temp<-which(mo_sis_pair$mo==x)
				if(length(temp)==0){
					temp<-c(NA,NA)
				}
				return(temp)
			}))
		
		

			focal_error_da1_ind<-c(match(c(focal_error_da1,focal_error_da1_da,focal_error_da1_da_da,focal_error_da1_da_da_da),feature_names),focal_error_da1_da_da_da_da_ind)
			focal_error_da2_ind<-c(match(c(focal_error_da2,focal_error_da2_da,focal_error_da2_da_da,focal_error_da1_da_da_da),feature_names),focal_error_da2_da_da_da_da_ind)
		
			focal_error_da1_da2_nonNAInd<-which(!is.na(focal_error_da1_ind)&!is.na(focal_error_da2_ind))
			focal_error_da1_da2_NAInd<-which(is.na(focal_error_da1_ind)|is.na(focal_error_da2_ind))
		
			temp<-focal_x[j,focal_error_da1_ind[focal_error_da1_da2_nonNAInd]]
			focal_x[j,focal_error_da1_ind[focal_error_da1_da2_nonNAInd]]<-focal_x[j,focal_error_da2_ind[focal_error_da1_da2_nonNAInd]]
			focal_x[j,focal_error_da2_ind[focal_error_da1_da2_nonNAInd]]<-temp
		
			focal_error_NA_ind<-na.omit(c(focal_error_da1_ind[focal_error_da1_da2_NAInd],focal_error_da2_ind[focal_error_da1_da2_NAInd]))
		
			focal_x[j,focal_error_NA_ind]<-NA
		}
	}
	
	X_mo<-focal_x[,match(mo_sis_pair$mo,feature_names)]
	X_da1<-focal_x[,match(mo_sis_pair$da1,feature_names)]
	delta_X<-X_da1-X_mo
	
	ccl_model_mother_coef_error<-ccl_model_mother_coef
	ccl_model_mother_coef_error[c(cell_5_mo_ind,cell_6_mo_ind,cell_7_mo_ind,cell_8_mo_ind)]<-sapply(c(cell_5_mo_ind,cell_6_mo_ind,cell_7_mo_ind,cell_8_mo_ind),function(iii){summary(lm(delta_X[,iii]~X_mo[,iii]))$coefficients[2,1]})
	
	accuracy_5[i]<-cor(ccl_model_mother_coef_error,ccl_model_mother_coef)
	
	plot(accuracy_5)
}



#

accuracy_4<-c()

cell_4_mo_ind<-unlist(lapply(cell_4_mo,function(x){which(mo_sis_pair$mo==x)}))


for(i in 1:50){

	focal_x<-ccl_raw_nonNA_scaling_noise

	print(i)
	
	set.seed(1)
	total_error<-t(sapply(sample(1:(num_of_embryo*length(cell_4_mo)),round(num_of_embryo*length(cell_4_mo)*0.1*i/50)),rowColRec,num_of_embryo,length(cell_4_mo)))
	error_mo<-lapply(1:num_of_embryo,function(ii){
		focalInd<-which(total_error[,1]==ii)
		if(length(focalInd)==0){
			return(NA)
		}else{
			return(cell_4_mo[total_error[focalInd,2]])
		}
	})
	
	for(j in 1:num_of_embryo){
		if(!is.na(error_mo[[j]][1])){
		
			focal_error_da1<-mo_sis_pair$da1[unlist(lapply(error_mo[[j]],function(x){which(mo_sis_pair$mo==x)}))]
			focal_error_da2<-mo_sis_pair$da2[unlist(lapply(error_mo[[j]],function(x){which(mo_sis_pair$mo==x)}))]
		
			focal_error_da1_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da1,function(x){which(mo_sis_pair$mo==x)}))]
			focal_error_da2_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da2,function(x){which(mo_sis_pair$mo==x)}))]	
		
			focal_error_da1_da_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da1_da,function(x){which(mo_sis_pair$mo==x)}))]
			focal_error_da2_da_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da2_da,function(x){which(mo_sis_pair$mo==x)}))]	
		
			focal_error_da1_da_da_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da1_da_da,function(x){which(mo_sis_pair$mo==x)}))]
			focal_error_da2_da_da_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da2_da_da,function(x){which(mo_sis_pair$mo==x)}))]	

			focal_error_da1_da_da_da_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da1_da_da_da,function(x){which(mo_sis_pair$mo==x)}))]
			focal_error_da2_da_da_da_da<-mo_sis_pair$da1[unlist(lapply(focal_error_da2_da_da_da,function(x){which(mo_sis_pair$mo==x)}))]	
			
			focal_error_da1_da_da_da_da_da_ind<-unlist(lapply(focal_error_da1_da_da_da_da,function(x){
				temp<-which(mo_sis_pair$mo==x)
				if(length(temp)==0){
					temp<-c(NA,NA)
				}
				return(temp)
			}))
			focal_error_da2_da_da_da_da_da_ind<-unlist(lapply(focal_error_da2_da_da_da_da,function(x){
				temp<-which(mo_sis_pair$mo==x)
				if(length(temp)==0){
					temp<-c(NA,NA)
				}
				return(temp)
			}))
		
		

			focal_error_da1_ind<-c(match(c(focal_error_da1,focal_error_da1_da,focal_error_da1_da_da,focal_error_da1_da_da_da,focal_error_da1_da_da_da_da),feature_names),focal_error_da1_da_da_da_da_da_ind)
			focal_error_da2_ind<-c(match(c(focal_error_da2,focal_error_da2_da,focal_error_da2_da_da,focal_error_da1_da_da_da,focal_error_da2_da_da_da_da),feature_names),focal_error_da2_da_da_da_da_da_ind)
		
			focal_error_da1_da2_nonNAInd<-which(!is.na(focal_error_da1_ind)&!is.na(focal_error_da2_ind))
			focal_error_da1_da2_NAInd<-which(is.na(focal_error_da1_ind)|is.na(focal_error_da2_ind))
		
			temp<-focal_x[j,focal_error_da1_ind[focal_error_da1_da2_nonNAInd]]
			focal_x[j,focal_error_da1_ind[focal_error_da1_da2_nonNAInd]]<-focal_x[j,focal_error_da2_ind[focal_error_da1_da2_nonNAInd]]
			focal_x[j,focal_error_da2_ind[focal_error_da1_da2_nonNAInd]]<-temp
		
			focal_error_NA_ind<-na.omit(c(focal_error_da1_ind[focal_error_da1_da2_NAInd],focal_error_da2_ind[focal_error_da1_da2_NAInd]))
		
			focal_x[j,focal_error_NA_ind]<-NA
		}
	}
	
	X_mo<-focal_x[,match(mo_sis_pair$mo,feature_names)]
	X_da1<-focal_x[,match(mo_sis_pair$da1,feature_names)]
	delta_X<-X_da1-X_mo
	
	ccl_model_mother_coef_error<-ccl_model_mother_coef
	ccl_model_mother_coef_error[c(cell_4_mo_ind,cell_5_mo_ind,cell_6_mo_ind,cell_7_mo_ind,cell_8_mo_ind)]<-sapply(c(cell_4_mo_ind,cell_5_mo_ind,cell_6_mo_ind,cell_7_mo_ind,cell_8_mo_ind),function(iii){summary(lm(delta_X[,iii]~X_mo[,iii]))$coefficients[2,1]})
	
	accuracy_4[i]<-cor(ccl_model_mother_coef_error,ccl_model_mother_coef)
	
	plot(accuracy_4)
}


#

cell_9_da<-intersect(na.omit(unique(c(lineage_info[9,]))),mo_sis_pair$da1)
cell_8_da<-intersect(na.omit(unique(c(lineage_info[8,]))),mo_sis_pair$da1)
cell_7_da<-intersect(na.omit(unique(c(lineage_info[7,]))),mo_sis_pair$da1)
cell_6_da<-intersect(na.omit(unique(c(lineage_info[6,]))),mo_sis_pair$da1)
cell_5_da<-intersect(na.omit(unique(c(lineage_info[5,]))),mo_sis_pair$da1)


tracking_error_CCL_01<-data.frame(accuracy=c(accuracy_4,accuracy_5,accuracy_6,accuracy_7,accuracy_8),
	generation=factor(c(rep(4,50),rep(5,50),rep(6,50),rep(7,50),rep(8,50)),level=c(4:8)),
	error_rate=c(c(1:50)*0.1/50*(length(cell_5_da)+length(cell_6_da)+length(cell_7_da)+length(cell_8_da)+length(cell_9_da))/345,
		c(1:50)*0.1/50*(length(cell_6_da)+length(cell_7_da)+length(cell_8_da)+length(cell_9_da))/345,
		c(1:50)*0.1/50*(length(cell_7_da)+length(cell_8_da)+length(cell_9_da))/345,
		c(1:50)*0.1/50*(length(cell_8_da)+length(cell_9_da))/345,
		c(1:50)*0.1/50*(length(cell_9_da))/345))

save(tracking_error_CCL_01,file=paste0(focal_projectFolderPath,"/data/confounding/tracking_error_CCL_01.RData"))
	

#####################################################################################################################3	


# x 

ggplot(tracking_error_CCL_01,aes(x=error_rate,y=accuracy^2,col=generation))+
geom_point(size=0.2,stroke=0.1)+geom_smooth(lwd=0.2)+
scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
scale_x_continuous(breaks=c(0,0.05,0.1),limits=c(0,0.1))+
labs(x="Tracking error rate",y=expression(paste("Accuracy ( ",italic(R)^2," )")),title=expression(italic(k)[CCL]))+
theme(legend.position=c(0.2,0.6))
	
	
ggsave("effect_tracking_error_01_kCCL.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

