

### marginal model 
### shuffled model 
### full model 
### nested model 
### descendant blocking 
### same-cell-cycle: M1->M11


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

library(foreach)
library(glmnet)
library(parallel)
library(doParallel)
library(ggplot2)



focal_projectFolderPath<-r"(E:\canalizationENV)"
#focal_projectFolderPath<-r"(/home/wjg/backup20221105/canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/data/spatial"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q11_cell350_alignment_xyz_noise.RData"
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


##################################################################################################
                                           _x
##################################################################################################

###

X_mo<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo


#

x_residuals<-sapply(1:dim(X_mo)[2],function(i){
	residuals(lm(delta_X[,i]~X_mo[,i]))
})


#

set.seed(1)
X_da1_shuffle<-X_da1[,sample(706)]
delta_X_shuffle<-X_da1_shuffle-X_mo


### marginal model - mother

xyz_x_model_mother_R<-sapply(1:dim(X_mo)[2],function(i){cor(delta_X[,i],X_mo[,i])})
xyz_x_model_mother_coef<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})
xyz_x_model_mo_da1_R<-sapply(1:dim(X_mo)[2],function(i){cor(X_mo[,i],X_da1[,i])})
xyz_x_model_mo_da1_pvalue<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_mo[,i],X_da1[,i])$p.value})

names(xyz_x_model_mother_R)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_x_model_mother_coef)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_x_model_mo_da1_R)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_x_model_mo_da1_pvalue)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_x_model_mother_R,file="xyz_x_model_mother_R.RData")
save(xyz_x_model_mother_coef,file="xyz_x_model_mother_coef.RData")
save(xyz_x_model_mo_da1_R,file="xyz_x_model_mo_da1_R.RData")
save(xyz_x_model_mo_da1_pvalue,file="xyz_x_model_mo_da1_pvalue.RData")


### marginal model - mother - shuffle

xyz_x_model_mother_R_shuffle<-sapply(1:dim(X_mo)[2],function(i){cor(delta_X_shuffle[,i],X_mo[,i])})
xyz_x_model_mother_coef_shuffle<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X_shuffle[,i]~X_mo[,i]))$coefficients[2,1]})
xyz_x_model_mo_da1_R_shuffle<-sapply(1:dim(X_mo)[2],function(i){cor(X_mo[,i],X_da1_shuffle[,i])})
xyz_x_model_mo_da1_pvalue_shuffle<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_mo[,i],X_da1_shuffle[,i])$p.value})

names(xyz_x_model_mother_R_shuffle)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_x_model_mother_coef_shuffle)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_x_model_mo_da1_R_shuffle)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_x_model_mo_da1_pvalue_shuffle)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_x_model_mother_R_shuffle,file="xyz_x_model_mother_R_shuffle.RData")
save(xyz_x_model_mother_coef_shuffle,file="xyz_x_model_mother_coef_shuffle.RData")
save(xyz_x_model_mo_da1_R_shuffle,file="xyz_x_model_mo_da1_R_shuffle.RData")
save(xyz_x_model_mo_da1_pvalue_shuffle,file="xyz_x_model_mo_da1_pvalue_shuffle.RData")


### marginal model - sister

xyz_x_model_sister_R<-sapply(1:dim(X_mo)[2],function(i){cor(delta_X[,i],X_da2[,i])})
xyz_x_model_sister_coef<-sapply(1:dim(X_mo)[2],function(i){
	if(is.na(colnames(X_da2)[i])){
		return(NA)
	}else{
		summary(lm(delta_X[,i]~X_da2[,i]))$coefficients[2,1]	
	}
})

names(xyz_x_model_sister_R)<-paste0(mo_sis_pair$da2,"_",mo_sis_pair$da1)
names(xyz_x_model_sister_coef)<-paste0(mo_sis_pair$da2,"_",mo_sis_pair$da1)

save(xyz_x_model_sister_R,file="xyz_x_model_sister_R.RData")
save(xyz_x_model_sister_coef,file="xyz_x_model_sister_coef.RData")


### nested model - mother and sister

xyz_x_model_mother_sister_coef<-sapply(1:dim(X_mo)[2],function(i){
	print(i)
	sisterPairInd<-which(colnames(X_mo)==colnames(X_da2)[i])
	
	temp<-rep(0,dim(X_mo)[2]+1)
	names(temp)<-c("intercept",colnames(X_mo))
	
	if(length(sisterPairInd)>0){
		tempCoef<-summary(lm(delta_X[,i]~X_mo[,c(i,sisterPairInd[1])]))$coefficients[,1]
		temp[c(1,c(i,sisterPairInd[1])+1)]<-tempCoef
	}
	
	tempCoef<-temp
	return(tempCoef)		

})

colnames(xyz_x_model_mother_sister_coef)<-colnames(X_da1)

xyz_x_model_mother_sister_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_x_model_mother_sister_coef
xyz_x_model_mother_sister_R<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],xyz_x_model_mother_sister_pred[,i])})
names(xyz_x_model_mother_sister_R)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da2,"_",mo_sis_pair$da1)

save(xyz_x_model_mother_sister_coef,file="xyz_x_model_mother_sister_coef.RData")
save(xyz_x_model_mother_sister_R,file="xyz_x_model_mother_sister_R.RData")


### full model

cl<-makeCluster(10,type="FORK")
registerDoParallel(cl)

xyz_x_model_full_coef<-sapply(1:dim(X_mo)[2],function(i){
	print(i)
	motherPairInd<-which(colnames(X_mo)==colnames(X_mo)[i])
	daughterPairInd<-which(colnames(X_mo)==colnames(X_da1)[i])
	dupInd<-which(duplicated(colnames(X_mo)))
	tempInd<-setdiff(union(dupInd,union(motherPairInd,daughterPairInd)),i)
	
	if(length(tempInd)==0){
		tempCoef<-coefficients(cv.glmnet(X_mo, delta_X[,i], nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
	}else{
		tempCoef<-coefficients(cv.glmnet(X_mo[,-tempInd], delta_X[,i],nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
		temp<-rep(0,dim(X_mo)[2]+1)
		names(temp)<-c("intercept",colnames(X_mo))
		temp[c(1,setdiff(1:dim(X_mo)[2],tempInd)+1)]<-tempCoef
		tempCoef<-temp
	}
	return(tempCoef)
})

colnames(xyz_x_model_full_coef)<-colnames(X_da1)

save(xyz_x_model_full_coef,file="xyz_x_model_full_coef.RData")


delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_x_model_full_coef

xyz_x_model_full_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

save(xyz_x_model_full_R2,file="xyz_x_model_full_R2.RData")



### nested model without mother

xyz_x_model_nested_without_mother_coef<-sapply(1:dim(X_mo)[2],function(i){
	print(i)
	motherPairInd<-which(colnames(X_mo)==colnames(X_mo)[i])
	daughterPairInd<-which(colnames(X_mo)==colnames(X_da1)[i])
	dupInd<-which(duplicated(colnames(X_mo)))
	tempInd<-union(dupInd,union(motherPairInd,daughterPairInd))
	
	if(length(tempInd)==0){
		tempCoef<-coefficients(cv.glmnet(X_mo, delta_X[,i], nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
	}else{
		tempCoef<-coefficients(cv.glmnet(X_mo[,-tempInd], delta_X[,i],nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
		temp<-rep(0,dim(X_mo)[2]+1)
		names(temp)<-c("intercept",colnames(X_mo))
		temp[c(1,setdiff(1:dim(X_mo)[2],tempInd)+1)]<-tempCoef
		tempCoef<-temp
	}
	return(tempCoef)
})

colnames(xyz_x_model_nested_without_mother_coef)<-colnames(X_da1)

save(xyz_x_model_nested_without_mother_coef,file="xyz_x_model_nested_without_mother_coef.RData")


delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_x_model_nested_without_mother_coef

xyz_x_model_nested_without_mother_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

save(xyz_x_model_nested_without_mother_R2,file="xyz_x_model_nested_without_mother_R2.RData")



### nested model without sister 

xyz_x_model_nested_without_sister_coef<-sapply(1:dim(X_mo)[2],function(i){
	print(i)
	motherPairInd<-which(colnames(X_mo)==colnames(X_mo)[i])
	daughterPairInd<-which(colnames(X_mo)==colnames(X_da1)[i])
	sisterPairInd<-which(colnames(X_mo)==colnames(X_da2)[i])
	dupInd<-which(duplicated(colnames(X_mo)))
	tempInd<-setdiff(union(union(dupInd,union(motherPairInd,daughterPairInd)),sisterPairInd),i)

	if(length(tempInd)==0){
		tempCoef<-coefficients(cv.glmnet(X_mo, delta_X[,i], nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
	}else{
		tempCoef<-coefficients(cv.glmnet(X_mo[,-tempInd], delta_X[,i],nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
		temp<-rep(0,dim(delta_X)[2]+1)
		names(temp)<-c("intercept",colnames(X_mo))
		temp[c(1,setdiff(1:dim(delta_X)[2],tempInd)+1)]<-tempCoef
		tempCoef<-temp
	}
	return(tempCoef)
})

colnames(xyz_x_model_nested_without_sister_coef)<-colnames(X_da1)

save(xyz_x_model_nested_without_sister_coef,file="xyz_x_model_nested_without_sister_coef.RData")


delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_x_model_nested_without_sister_coef

xyz_x_model_nested_without_sister_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

save(xyz_x_model_nested_without_sister_R2,file="xyz_x_model_nested_without_sister_R2.RData")


### full model : block descendant (x)

cl<-makeCluster(10,type="FORK")
registerDoParallel(cl)

xyz_x_model_block_descendant_coef<-sapply(1:dim(X_mo)[2],function(i){
	print(i)
	motherPairInd<-which(colnames(X_mo)==colnames(X_mo)[i])
	daughterPairInd<-which(colnames(X_mo)==colnames(X_da1)[i])
	dupInd<-which(duplicated(colnames(X_mo)))

	focal_lineage<-do.call(rbind,lapply(which(lineage_info==unlist(strsplit(colnames(X_da1)[i],"_x"))),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2]))
	#including itself
	descendant<-unique(unlist(lapply(1:dim(focal_lineage)[1],function(j){lineage_info[focal_lineage[j,1]:10,focal_lineage[j,2]]})))
	descendantInd<-unlist(lapply(paste0(descendant,"_x"),function(x){which(colnames(X_mo)==x)}))
	
	tempInd<-setdiff(union(descendantInd,union(dupInd,union(motherPairInd,daughterPairInd))),i)

	
	if(length(tempInd)==0){
		tempCoef<-coefficients(cv.glmnet(X_mo, delta_X[,i], nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
	}else{
		tempCoef<-coefficients(cv.glmnet(X_mo[,-tempInd], delta_X[,i],nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
		temp<-rep(0,dim(X_mo)[2]+1)
		names(temp)<-c("intercept",colnames(X_mo))
		temp[c(1,setdiff(1:dim(X_mo)[2],tempInd)+1)]<-tempCoef
		tempCoef<-temp
	}
	return(tempCoef)
})


colnames(xyz_x_model_block_descendant_coef)<-colnames(X_da1)

save(xyz_x_model_block_descendant_coef,file="xyz_x_model_block_descendant_coef.RData")


delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_x_model_block_descendant_coef

xyz_x_model_block_descendant_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

save(xyz_x_model_block_descendant_R2,file="xyz_x_model_block_descendant_R2.RData")


##################################################################################################
                                           _y
##################################################################################################

#

X_mo<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo

#

set.seed(1)
X_da1_shuffle<-X_da1[,sample(706)]
delta_X_shuffle<-X_da1_shuffle-X_mo


### marginal model - mother

xyz_y_model_mother_R<-sapply(1:dim(X_mo)[2],function(i){cor(delta_X[,i],X_mo[,i])})
xyz_y_model_mother_coef<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})
xyz_y_model_mo_da1_R<-sapply(1:dim(X_mo)[2],function(i){cor(X_mo[,i],X_da1[,i])})
xyz_y_model_mo_da1_pvalue<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_mo[,i],X_da1[,i])$p.value})

names(xyz_y_model_mother_R)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_y_model_mother_coef)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_y_model_mo_da1_R)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_y_model_mo_da1_pvalue)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_y_model_mother_R,file="xyz_y_model_mother_R.RData")
save(xyz_y_model_mother_coef,file="xyz_y_model_mother_coef.RData")
save(xyz_y_model_mo_da1_R,file="xyz_y_model_mo_da1_R.RData")
save(xyz_y_model_mo_da1_pvalue,file="xyz_y_model_mo_da1_pvalue.RData")


### marginal model - mother - shuffle

xyz_y_model_mother_R_shuffle<-sapply(1:dim(X_mo)[2],function(i){cor(delta_X_shuffle[,i],X_mo[,i])})
xyz_y_model_mother_coef_shuffle<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X_shuffle[,i]~X_mo[,i]))$coefficients[2,1]})
xyz_y_model_mo_da1_R_shuffle<-sapply(1:dim(X_mo)[2],function(i){cor(X_mo[,i],X_da1_shuffle[,i])})
xyz_y_model_mo_da1_pvalue_shuffle<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_mo[,i],X_da1_shuffle[,i])$p.value})

names(xyz_y_model_mother_R_shuffle)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_y_model_mother_coef_shuffle)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_y_model_mo_da1_R_shuffle)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_y_model_mo_da1_pvalue_shuffle)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_y_model_mother_R_shuffle,file="xyz_y_model_mother_R_shuffle.RData")
save(xyz_y_model_mother_coef_shuffle,file="xyz_y_model_mother_coef_shuffle.RData")
save(xyz_y_model_mo_da1_R_shuffle,file="xyz_y_model_mo_da1_R_shuffle.RData")
save(xyz_y_model_mo_da1_pvalue_shuffle,file="xyz_y_model_mo_da1_pvalue_shuffle.RData")


### marginal model - sister

xyz_y_model_sister_R<-sapply(1:dim(X_mo)[2],function(i){cor(delta_X[,i],X_da2[,i])})
xyz_y_model_sister_coef<-sapply(1:dim(X_mo)[2],function(i){
	if(is.na(colnames(X_da2)[i])){
		return(NA)
	}else{
		summary(lm(delta_X[,i]~X_da2[,i]))$coefficients[2,1]	
	}
})

names(xyz_y_model_sister_R)<-paste0(mo_sis_pair$da2,"_",mo_sis_pair$da1)
names(xyz_y_model_sister_coef)<-paste0(mo_sis_pair$da2,"_",mo_sis_pair$da1)

save(xyz_y_model_sister_R,file="xyz_y_model_sister_R.RData")
save(xyz_y_model_sister_coef,file="xyz_y_model_sister_coef.RData")


### nested model - mother and sister

xyz_y_model_mother_sister_coef<-sapply(1:dim(X_mo)[2],function(i){
	print(i)
	sisterPairInd<-which(colnames(X_mo)==colnames(X_da2)[i])
	
	temp<-rep(0,dim(X_mo)[2]+1)
	names(temp)<-c("intercept",colnames(X_mo))
	
	if(length(sisterPairInd)>0){
		tempCoef<-summary(lm(delta_X[,i]~X_mo[,c(i,sisterPairInd[1])]))$coefficients[,1]
		temp[c(1,c(i,sisterPairInd[1])+1)]<-tempCoef
	}
	
	tempCoef<-temp
	return(tempCoef)		

})

colnames(xyz_y_model_mother_sister_coef)<-colnames(X_da1)

xyz_y_model_mother_sister_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_y_model_mother_sister_coef
xyz_y_model_mother_sister_R<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],xyz_y_model_mother_sister_pred[,i])})
names(xyz_y_model_mother_sister_R)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da2,"_",mo_sis_pair$da1)

save(xyz_y_model_mother_sister_coef,file="xyz_y_model_mother_sister_coef.RData")
save(xyz_y_model_mother_sister_R,file="xyz_y_model_mother_sister_R.RData")


### full model

cl<-makeCluster(10,type="FORK")
registerDoParallel(cl)

xyz_y_model_full_coef<-sapply(1:dim(X_mo)[2],function(i){
	print(i)
	motherPairInd<-which(colnames(X_mo)==colnames(X_mo)[i])
	daughterPairInd<-which(colnames(X_mo)==colnames(X_da1)[i])
	dupInd<-which(duplicated(colnames(X_mo)))
	tempInd<-setdiff(union(dupInd,union(motherPairInd,daughterPairInd)),i)
	
	if(length(tempInd)==0){
		tempCoef<-coefficients(cv.glmnet(X_mo, delta_X[,i], nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
	}else{
		tempCoef<-coefficients(cv.glmnet(X_mo[,-tempInd], delta_X[,i],nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
		temp<-rep(0,dim(X_mo)[2]+1)
		names(temp)<-c("intercept",colnames(X_mo))
		temp[c(1,setdiff(1:dim(X_mo)[2],tempInd)+1)]<-tempCoef
		tempCoef<-temp
	}
	return(tempCoef)
})

colnames(xyz_y_model_full_coef)<-colnames(X_da1)

save(xyz_y_model_full_coef,file="xyz_y_model_full_coef.RData")


delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_y_model_full_coef

xyz_y_model_full_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

save(xyz_y_model_full_R2,file="xyz_y_model_full_R2.RData")


### nested model without mother

xyz_y_model_nested_without_mother_coef<-sapply(1:dim(X_mo)[2],function(i){
	print(i)
	motherPairInd<-which(colnames(X_mo)==colnames(X_mo)[i])
	daughterPairInd<-which(colnames(X_mo)==colnames(X_da1)[i])
	dupInd<-which(duplicated(colnames(X_mo)))
	tempInd<-union(dupInd,union(motherPairInd,daughterPairInd))
	
	if(length(tempInd)==0){
		tempCoef<-coefficients(cv.glmnet(X_mo, delta_X[,i], nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
	}else{
		tempCoef<-coefficients(cv.glmnet(X_mo[,-tempInd], delta_X[,i],nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
		temp<-rep(0,dim(X_mo)[2]+1)
		names(temp)<-c("intercept",colnames(X_mo))
		temp[c(1,setdiff(1:dim(X_mo)[2],tempInd)+1)]<-tempCoef
		tempCoef<-temp
	}
	return(tempCoef)
})

colnames(xyz_y_model_nested_without_mother_coef)<-colnames(X_da1)

save(xyz_y_model_nested_without_mother_coef,file="xyz_y_model_nested_without_mother_coef.RData")


delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_y_model_nested_without_mother_coef

xyz_y_model_nested_without_mother_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

save(xyz_y_model_nested_without_mother_R2,file="xyz_y_model_nested_without_mother_R2.RData")


### nested model without sister 

xyz_y_model_nested_without_sister_coef<-sapply(1:dim(X_mo)[2],function(i){
	print(i)
	motherPairInd<-which(colnames(X_mo)==colnames(X_mo)[i])
	daughterPairInd<-which(colnames(X_mo)==colnames(X_da1)[i])
	sisterPairInd<-which(colnames(X_mo)==colnames(X_da2)[i])
	dupInd<-which(duplicated(colnames(X_mo)))
	tempInd<-setdiff(union(union(dupInd,union(motherPairInd,daughterPairInd)),sisterPairInd),i)

	if(length(tempInd)==0){
		tempCoef<-coefficients(cv.glmnet(X_mo, delta_X[,i], nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
	}else{
		tempCoef<-coefficients(cv.glmnet(X_mo[,-tempInd], delta_X[,i],nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
		temp<-rep(0,dim(delta_X)[2]+1)
		names(temp)<-c("intercept",colnames(X_mo))
		temp[c(1,setdiff(1:dim(delta_X)[2],tempInd)+1)]<-tempCoef
		tempCoef<-temp
	}
	return(tempCoef)
})

colnames(xyz_y_model_nested_without_sister_coef)<-colnames(X_da1)

save(xyz_y_model_nested_without_sister_coef,file="xyz_y_model_nested_without_sister_coef.RData")


delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_y_model_nested_without_sister_coef

xyz_y_model_nested_without_sister_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

save(xyz_y_model_nested_without_sister_R2,file="xyz_y_model_nested_without_sister_R2.RData")



### full model : block descendant (y)

cl<-makeCluster(10,type="FORK")
registerDoParallel(cl)

xyz_y_model_block_descendant_coef<-sapply(1:dim(X_mo)[2],function(i){
	print(i)
	motherPairInd<-which(colnames(X_mo)==colnames(X_mo)[i])
	daughterPairInd<-which(colnames(X_mo)==colnames(X_da1)[i])
	dupInd<-which(duplicated(colnames(X_mo)))
	
	focal_lineage<-do.call(rbind,lapply(which(lineage_info==unlist(strsplit(colnames(X_da1)[i],"_y"))),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2]))
	#including itself
	descendant<-unique(unlist(lapply(1:dim(focal_lineage)[1],function(j){lineage_info[focal_lineage[j,1]:10,focal_lineage[j,2]]})))
	descendantInd<-unlist(lapply(paste0(descendant,"_y"),function(x){which(colnames(X_mo)==x)}))
	
	tempInd<-setdiff(union(descendantInd,union(dupInd,union(motherPairInd,daughterPairInd))),i)
	
	if(length(tempInd)==0){
		tempCoef<-coefficients(cv.glmnet(X_mo, delta_X[,i], nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
	}else{
		tempCoef<-coefficients(cv.glmnet(X_mo[,-tempInd], delta_X[,i],nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
		temp<-rep(0,dim(X_mo)[2]+1)
		names(temp)<-c("intercept",colnames(X_mo))
		temp[c(1,setdiff(1:dim(X_mo)[2],tempInd)+1)]<-tempCoef
		tempCoef<-temp
	}
	return(tempCoef)
})

colnames(xyz_y_model_block_descendant_coef)<-colnames(X_da1)

save(xyz_y_model_block_descendant_coef,file="xyz_y_model_block_descendant_coef.RData")


delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_y_model_block_descendant_coef

xyz_y_model_block_descendant_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

save(xyz_y_model_block_descendant_R2,file="xyz_y_model_block_descendant_R2.RData")


##################################################################################################
                                           _z
##################################################################################################

#

X_mo<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo

#
set.seed(1)
X_da1_shuffle<-X_da1[,sample(706)]
delta_X_shuffle<-X_da1_shuffle-X_mo


### marginal model - mother

xyz_z_model_mother_R<-sapply(1:dim(X_mo)[2],function(i){cor(delta_X[,i],X_mo[,i])})
xyz_z_model_mother_coef<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})
xyz_z_model_mo_da1_R<-sapply(1:dim(X_mo)[2],function(i){cor(X_mo[,i],X_da1[,i])})
xyz_z_model_mo_da1_pvalue<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_mo[,i],X_da1[,i])$p.value})

names(xyz_z_model_mother_R)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_z_model_mother_coef)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_z_model_mo_da1_R)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_z_model_mo_da1_pvalue)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_z_model_mother_R,file="xyz_z_model_mother_R.RData")
save(xyz_z_model_mother_coef,file="xyz_z_model_mother_coef.RData")
save(xyz_z_model_mo_da1_R,file="xyz_z_model_mo_da1_R.RData")
save(xyz_z_model_mo_da1_pvalue,file="xyz_z_model_mo_da1_pvalue.RData")


### marginal model - mother - shuffle

xyz_z_model_mother_R_shuffle<-sapply(1:dim(X_mo)[2],function(i){cor(delta_X_shuffle[,i],X_mo[,i])})
xyz_z_model_mother_coef_shuffle<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X_shuffle[,i]~X_mo[,i]))$coefficients[2,1]})
xyz_z_model_mo_da1_R_shuffle<-sapply(1:dim(X_mo)[2],function(i){cor(X_mo[,i],X_da1_shuffle[,i])})
xyz_z_model_mo_da1_pvalue_shuffle<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_mo[,i],X_da1_shuffle[,i])$p.value})

names(xyz_z_model_mother_R_shuffle)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_z_model_mother_coef_shuffle)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_z_model_mo_da1_R_shuffle)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(xyz_z_model_mo_da1_pvalue_shuffle)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(xyz_z_model_mother_R_shuffle,file="xyz_z_model_mother_R_shuffle.RData")
save(xyz_z_model_mother_coef_shuffle,file="xyz_z_model_mother_coef_shuffle.RData")
save(xyz_z_model_mo_da1_R_shuffle,file="xyz_z_model_mo_da1_R_shuffle.RData")
save(xyz_z_model_mo_da1_pvalue_shuffle,file="xyz_z_model_mo_da1_pvalue_shuffle.RData")


### marginal model - sister

xyz_z_model_sister_R<-sapply(1:dim(X_mo)[2],function(i){cor(delta_X[,i],X_da2[,i])})
xyz_z_model_sister_coef<-sapply(1:dim(X_mo)[2],function(i){
	if(is.na(colnames(X_da2)[i])){
		return(NA)
	}else{
		summary(lm(delta_X[,i]~X_da2[,i]))$coefficients[2,1]	
	}
})

names(xyz_z_model_sister_R)<-paste0(mo_sis_pair$da2,"_",mo_sis_pair$da1)
names(xyz_z_model_sister_coef)<-paste0(mo_sis_pair$da2,"_",mo_sis_pair$da1)

save(xyz_z_model_sister_R,file="xyz_z_model_sister_R.RData")
save(xyz_z_model_sister_coef,file="xyz_z_model_sister_coef.RData")


### nested model - mother and sister

xyz_z_model_mother_sister_coef<-sapply(1:dim(X_mo)[2],function(i){
	print(i)
	sisterPairInd<-which(colnames(X_mo)==colnames(X_da2)[i])
	
	temp<-rep(0,dim(X_mo)[2]+1)
	names(temp)<-c("intercept",colnames(X_mo))
	
	if(length(sisterPairInd)>0){
		tempCoef<-summary(lm(delta_X[,i]~X_mo[,c(i,sisterPairInd[1])]))$coefficients[,1]
		temp[c(1,c(i,sisterPairInd[1])+1)]<-tempCoef
	}
	
	tempCoef<-temp
	return(tempCoef)		

})

colnames(xyz_z_model_mother_sister_coef)<-colnames(X_da1)

xyz_z_model_mother_sister_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_z_model_mother_sister_coef
xyz_z_model_mother_sister_R<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],xyz_z_model_mother_sister_pred[,i])})
names(xyz_z_model_mother_sister_R)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da2,"_",mo_sis_pair$da1)

save(xyz_z_model_mother_sister_coef,file="xyz_z_model_mother_sister_coef.RData")
save(xyz_z_model_mother_sister_R,file="xyz_z_model_mother_sister_R.RData")


### full model

cl<-makeCluster(10,type="FORK")
registerDoParallel(cl)

xyz_z_model_full_coef<-sapply(1:dim(X_mo)[2],function(i){
	print(i)
	motherPairInd<-which(colnames(X_mo)==colnames(X_mo)[i])
	daughterPairInd<-which(colnames(X_mo)==colnames(X_da1)[i])
	dupInd<-which(duplicated(colnames(X_mo)))
	tempInd<-setdiff(union(dupInd,union(motherPairInd,daughterPairInd)),i)
	
	if(length(tempInd)==0){
		tempCoef<-coefficients(cv.glmnet(X_mo, delta_X[,i], nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
	}else{
		tempCoef<-coefficients(cv.glmnet(X_mo[,-tempInd], delta_X[,i],nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
		temp<-rep(0,dim(X_mo)[2]+1)
		names(temp)<-c("intercept",colnames(X_mo))
		temp[c(1,setdiff(1:dim(X_mo)[2],tempInd)+1)]<-tempCoef
		tempCoef<-temp
	}
	return(tempCoef)
})

colnames(xyz_z_model_full_coef)<-colnames(X_da1)

save(xyz_z_model_full_coef,file="xyz_z_model_full_coef.RData")


delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_z_model_full_coef

xyz_z_model_full_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

save(xyz_z_model_full_R2,file="xyz_z_model_full_R2.RData")


### nested model without mother

xyz_z_model_nested_without_mother_coef<-sapply(1:dim(X_mo)[2],function(i){
	print(i)
	motherPairInd<-which(colnames(X_mo)==colnames(X_mo)[i])
	daughterPairInd<-which(colnames(X_mo)==colnames(X_da1)[i])
	dupInd<-which(duplicated(colnames(X_mo)))
	tempInd<-union(dupInd,union(motherPairInd,daughterPairInd))
	
	if(length(tempInd)==0){
		tempCoef<-coefficients(cv.glmnet(X_mo, delta_X[,i], nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
	}else{
		tempCoef<-coefficients(cv.glmnet(X_mo[,-tempInd], delta_X[,i],nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
		temp<-rep(0,dim(X_mo)[2]+1)
		names(temp)<-c("intercept",colnames(X_mo))
		temp[c(1,setdiff(1:dim(X_mo)[2],tempInd)+1)]<-tempCoef
		tempCoef<-temp
	}
	return(tempCoef)
})

colnames(xyz_z_model_nested_without_mother_coef)<-colnames(X_da1)

save(xyz_z_model_nested_without_mother_coef,file="xyz_z_model_nested_without_mother_coef.RData")


delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_z_model_nested_without_mother_coef

xyz_z_model_nested_without_mother_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

save(xyz_z_model_nested_without_mother_R2,file="xyz_z_model_nested_without_mother_R2.RData")


### nested model without sister 

xyz_z_model_nested_without_sister_coef<-sapply(1:dim(X_mo)[2],function(i){
	print(i)
	motherPairInd<-which(colnames(X_mo)==colnames(X_mo)[i])
	daughterPairInd<-which(colnames(X_mo)==colnames(X_da1)[i])
	sisterPairInd<-which(colnames(X_mo)==colnames(X_da2)[i])
	dupInd<-which(duplicated(colnames(X_mo)))
	tempInd<-setdiff(union(union(dupInd,union(motherPairInd,daughterPairInd)),sisterPairInd),i)

	if(length(tempInd)==0){
		tempCoef<-coefficients(cv.glmnet(X_mo, delta_X[,i], nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
	}else{
		tempCoef<-coefficients(cv.glmnet(X_mo[,-tempInd], delta_X[,i],nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
		temp<-rep(0,dim(delta_X)[2]+1)
		names(temp)<-c("intercept",colnames(X_mo))
		temp[c(1,setdiff(1:dim(delta_X)[2],tempInd)+1)]<-tempCoef
		tempCoef<-temp
	}
	return(tempCoef)
})

colnames(xyz_z_model_nested_without_sister_coef)<-colnames(X_da1)

save(xyz_z_model_nested_without_sister_coef,file="xyz_z_model_nested_without_sister_coef.RData")


delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_z_model_nested_without_sister_coef

xyz_z_model_nested_without_sister_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

save(xyz_z_model_nested_without_sister_R2,file="xyz_z_model_nested_without_sister_R2.RData")



### full model : block descendant (z)

cl<-makeCluster(10,type="FORK")
registerDoParallel(cl)

xyz_z_model_block_descendant_coef<-sapply(1:dim(X_mo)[2],function(i){
	print(i)
	motherPairInd<-which(colnames(X_mo)==colnames(X_mo)[i])
	daughterPairInd<-which(colnames(X_mo)==colnames(X_da1)[i])
	dupInd<-which(duplicated(colnames(X_mo)))
	
	focal_lineage<-do.call(rbind,lapply(which(lineage_info==unlist(strsplit(colnames(X_da1)[i],"_z"))),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2]))
	#including itself
	descendant<-unique(unlist(lapply(1:dim(focal_lineage)[1],function(j){lineage_info[focal_lineage[j,1]:10,focal_lineage[j,2]]})))
	descendantInd<-unlist(lapply(paste0(descendant,"_z"),function(x){which(colnames(X_mo)==x)}))

	tempInd<-setdiff(union(descendantInd,union(dupInd,union(motherPairInd,daughterPairInd))),i)
	
	if(length(tempInd)==0){
		tempCoef<-coefficients(cv.glmnet(X_mo, delta_X[,i], nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
	}else{
		tempCoef<-coefficients(cv.glmnet(X_mo[,-tempInd], delta_X[,i],nlambda = 200,nfolds = 10,parallel=TRUE),s="lambda.min")[,1]
		temp<-rep(0,dim(X_mo)[2]+1)
		names(temp)<-c("intercept",colnames(X_mo))
		temp[c(1,setdiff(1:dim(X_mo)[2],tempInd)+1)]<-tempCoef
		tempCoef<-temp
	}
	return(tempCoef)
})


colnames(xyz_z_model_block_descendant_coef)<-colnames(X_da1)


save(xyz_z_model_block_descendant_coef,file="xyz_z_model_block_descendant_coef.RData")


delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_z_model_block_descendant_coef

xyz_z_model_block_descendant_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

save(xyz_z_model_block_descendant_R2,file="xyz_z_model_block_descendant_R2.RData")


##################################################################################################

##### remove last generation 

cell_gen<-sapply(match(feature_names,lineage_info),function(x){rowColRec(x,dim(lineage_info)[1],dim(lineage_info)[2])[1]})
g9Ind<-which(cell_gen<10)

### M11-M1 (x) : g9

delta_X=q11_cell350_alignment_xyz_noise[,g9Ind]-q1_cell350_alignment_xyz_noise[,g9Ind]
	
X_mo<-q1_cell350_alignment_xyz_noise[,g9Ind]


cl<-makeCluster(35,type="FORK")
registerDoParallel(cl)

xyz_x_M1M11g9_full_coef<-foreach(i=1:dim(delta_X)[2],.combine=cbind)%dopar%{
	print(i)
	tempCoef<-coefficients(cv.glmnet(X_mo, delta_X[,i],nlambda = 200,nfolds = 10,parallel=FALSE),s="lambda.min")[,1]
	return(tempCoef)
}

stopCluster(cl)

colnames(xyz_x_M1M11g9_full_coef)<-rownames(xyz_x_M1M11g9_full_coef)[-1]

delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_x_M1M11g9_full_coef

xyz_x_M1M11g9_full_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

names(xyz_x_M1M11g9_full_R2)<-rownames(xyz_x_M1M11g9_full_coef)[-1]

quantile(Re(eigen(t(xyz_x_M1M11g9_full_coef[-1,]))$values))

quantile(xyz_x_M1M11g9_full_R2)


save(xyz_x_M1M11g9_full_coef,file="xyz_x_M1M11g9_full_coef.RData")
save(xyz_x_M1M11g9_full_R2,file="xyz_x_M1M11g9_full_R2.RData")




### M11-M1 (y) : g9

delta_X=q11_cell350_alignment_xyz_noise[,g9Ind+714]-q1_cell350_alignment_xyz_noise[,g9Ind+714]
	
X_mo<-q1_cell350_alignment_xyz_noise[,g9Ind+714]


cl<-makeCluster(35,type="FORK")
registerDoParallel(cl)

xyz_y_M1M11g9_full_coef<-foreach(i=1:dim(delta_X)[2],.combine=cbind)%dopar%{
	print(i)
	tempCoef<-coefficients(cv.glmnet(X_mo, delta_X[,i],nlambda = 200,nfolds = 10,parallel=FALSE),s="lambda.min")[,1]
	return(tempCoef)
}

stopCluster(cl)


colnames(xyz_y_M1M11g9_full_coef)<-rownames(xyz_y_M1M11g9_full_coef)[-1]

delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_y_M1M11g9_full_coef

xyz_y_M1M11g9_full_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

names(xyz_y_M1M11g9_full_R2)<-rownames(xyz_y_M1M11g9_full_coef)[-1]

quantile(Re(eigen(t(xyz_y_M1M11g9_full_coef[-1,]))$values))

quantile(xyz_y_M1M11g9_full_R2)


save(xyz_y_M1M11g9_full_coef,file="xyz_y_M1M11g9_full_coef.RData")
save(xyz_y_M1M11g9_full_R2,file="xyz_y_M1M11g9_full_R2.RData")


### M11-M1 (z) : g9

delta_X=q11_cell350_alignment_xyz_noise[,g9Ind+(714*2)]-q1_cell350_alignment_xyz_noise[,g9Ind+(714*2)]
	
X_mo<-q1_cell350_alignment_xyz_noise[,g9Ind+(714*2)]


cl<-makeCluster(40,type="FORK")
registerDoParallel(cl)

xyz_z_M1M11g9_full_coef<-foreach(i=1:dim(delta_X)[2],.combine=cbind)%dopar%{
	print(i)
	tempCoef<-coefficients(cv.glmnet(X_mo, delta_X[,i],nlambda = 200,nfolds = 10,parallel=FALSE),s="lambda.min")[,1]
	return(tempCoef)
}

stopCluster(cl)

colnames(xyz_z_M1M11g9_full_coef)<-rownames(xyz_z_M1M11g9_full_coef)[-1]

delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%xyz_z_M1M11g9_full_coef

xyz_z_M1M11g9_full_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

names(xyz_z_M1M11g9_full_R2)<-rownames(xyz_z_M1M11g9_full_coef)[-1]


quantile(Re(eigen(t(xyz_z_M1M11g9_full_coef[-1,]))$values))

quantile(xyz_z_M1M11g9_full_R2)


save(xyz_z_M1M11g9_full_R2,file="xyz_z_M1M11g9_full_R2.RData")
save(xyz_z_M1M11g9_full_coef,file="xyz_z_M1M11g9_full_coef.RData")




























