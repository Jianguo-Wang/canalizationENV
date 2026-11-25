

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

#

kTest<-function(X,Y,refInd,perturbInd){

	kTest_result<-data.frame(t(sapply(1:dim(X)[2],function(j){
		model_ref<-lm(Y[refInd,j]~X[refInd,j])
		model_perturb<-lm(Y[perturbInd,j]~X[perturbInd,j])
		
		k_ref<-coefficients(summary(model_ref))[2,1]
		k_ref_se<-coefficients(summary(model_ref))[2,2]
		k_perturb<-coefficients(summary(model_perturb))[2,1]
		k_perturb_se<-coefficients(summary(model_perturb))[2,2]	
		
		k_diff<-k_perturb-k_ref
		k_diff_se<-sqrt(k_ref_se^2+k_perturb_se^2)
		k_diff_p<-2*min(pt(k_diff/k_diff_se,length(c(refInd,perturbInd))-2),
			1-pt(k_diff/k_diff_se,length(c(refInd,perturbInd))-2))
		return(c(k_diff,k_diff_se,k_diff_p,k_ref,k_perturb))
	})))
	rownames(kTest_result)<-paste(colnames(X),colnames(Y),sep="_")
	names(kTest_result)<-c("k_diff","k_diff_se","k_diff_p","k_ref","k_perturb")
	
	return(kTest_result)

}


#######################################################################

#######################################################################

library(foreach)
library(glmnet)
library(parallel)
library(doParallel)



focal_projectFolderPath<-r"(F:\wangj\Desktop\temp\canalizationENV)"
#focal_projectFolderPath<-r"(/home/wjg/backup20221105/canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/data/temporal"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling.RData"
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling_noise.RData"
	),projectFolderPath=focal_projectFolderPath)



###

feature_names<-unlist(strsplit(colnames(ccl_raw_nonNA_scaling_noise),"_ccl"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]

##

X_mo<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$mo,feature_names)]
X_da1<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$da1,feature_names)]
X_da2<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo

#

set.seed(1)
X_da1_shuffle<-X_da1[,sample(345)]
delta_X_shuffle<-X_da1_shuffle-X_mo


#######################################################################


### marginal model - mother

ccl_model_mother_R<-sapply(1:dim(X_mo)[2],function(i){cor(delta_X[,i],X_mo[,i])})
ccl_model_mother_coef<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[2,1]})
ccl_model_mo_da1_R<-sapply(1:dim(X_mo)[2],function(i){cor(X_mo[,i],X_da1[,i])})
ccl_model_mo_da1_pvalue<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_mo[,i],X_da1[,i])$p.value})

names(ccl_model_mother_R)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(ccl_model_mother_coef)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(ccl_model_mo_da1_R)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(ccl_model_mo_da1_pvalue)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)


save(ccl_model_mother_R,file="ccl_model_mother_R.RData")
save(ccl_model_mother_coef,file="ccl_model_mother_coef.RData")
save(ccl_model_mo_da1_R,file="ccl_model_mo_da1_R.RData")
save(ccl_model_mo_da1_pvalue,file="ccl_model_mo_da1_pvalue.RData")


### marginal model - mother - shuffle

ccl_model_mother_R_shuffle<-sapply(1:dim(X_mo)[2],function(i){cor(delta_X_shuffle[,i],X_mo[,i])})
ccl_model_mother_coef_shuffle<-sapply(1:dim(X_mo)[2],function(i){summary(lm(delta_X_shuffle[,i]~X_mo[,i]))$coefficients[2,1]})
ccl_model_mo_da1_R_shuffle<-sapply(1:dim(X_mo)[2],function(i){cor(X_mo[,i],X_da1_shuffle[,i])})
ccl_model_mo_da1_pvalue_shuffle<-sapply(1:dim(X_mo)[2],function(i){cor.test(X_mo[,i],X_da1_shuffle[,i])$p.value})

names(ccl_model_mother_R_shuffle)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(ccl_model_mother_coef_shuffle)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(ccl_model_mo_da1_R_shuffle)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)
names(ccl_model_mo_da1_pvalue_shuffle)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da1)

save(ccl_model_mother_R_shuffle,file="ccl_model_mother_R_shuffle.RData")
save(ccl_model_mother_coef_shuffle,file="ccl_model_mother_coef_shuffle.RData")
save(ccl_model_mo_da1_R_shuffle,file="ccl_model_mo_da1_R_shuffle.RData")
save(ccl_model_mo_da1_pvalue_shuffle,file="ccl_model_mo_da1_pvalue_shuffle.RData")



### marginal model - sister

ccl_model_sister_R<-sapply(1:dim(X_mo)[2],function(i){cor(delta_X[,i],X_da2[,i])})
ccl_model_sister_coef<-sapply(1:dim(X_mo)[2],function(i){
	if(is.na(colnames(X_da2)[i])){
		return(NA)
	}else{
		summary(lm(delta_X[,i]~X_da2[,i]))$coefficients[2,1]	
	}
})

names(ccl_model_sister_R)<-paste0(mo_sis_pair$da2,"_",mo_sis_pair$da1)
names(ccl_model_sister_coef)<-paste0(mo_sis_pair$da2,"_",mo_sis_pair$da1)

save(ccl_model_sister_R,file="ccl_model_sister_R.RData")
save(ccl_model_sister_coef,file="ccl_model_sister_coef.RData")



### mother and sister


ccl_model_mother_sister_coef<-sapply(1:dim(X_mo)[2],function(i){
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

colnames(ccl_model_mother_sister_coef)<-colnames(X_da1)

save(ccl_model_mother_sister_coef,file="ccl_model_mother_sister_coef.RData")

ccl_model_mother_sister_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%ccl_model_mother_sister_coef
ccl_model_mother_sister_R<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],ccl_model_mother_sister_pred[,i])})
names(ccl_model_mother_sister_R)<-paste0(mo_sis_pair$mo,"_",mo_sis_pair$da2,"_",mo_sis_pair$da1)

save(ccl_model_mother_sister_R,file="ccl_model_mother_sister_R.RData")


### full model

cl<-makeCluster(10,type="FORK")
registerDoParallel(cl)

ccl_model_full_coef<-sapply(1:dim(X_mo)[2],function(i){
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

colnames(ccl_model_full_coef)<-colnames(X_da1)

save(ccl_model_full_coef,file="ccl_model_full_coef.RData")


delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%ccl_model_full_coef

ccl_model_full_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

save(ccl_model_full_R2,file="ccl_model_full_R2.RData")


### nested model without mother

cl<-makeCluster(10,type="FORK")
registerDoParallel(cl)


ccl_model_nested_without_mother_coef<-sapply(1:dim(X_mo)[2],function(i){
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

colnames(ccl_model_nested_without_mother_coef)<-colnames(X_da1)

save(ccl_model_nested_without_mother_coef,file="ccl_model_nested_without_mother_coef.RData")


delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%ccl_model_nested_without_mother_coef

ccl_model_nested_without_mother_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

save(ccl_model_nested_without_mother_R2,file="ccl_model_nested_without_mother_R2.RData")


### nested model without sister 

cl<-makeCluster(10,type="FORK")
registerDoParallel(cl)

ccl_model_nested_without_sister_coef<-sapply(1:dim(X_mo)[2],function(i){
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

colnames(ccl_model_nested_without_sister_coef)<-colnames(X_da1)

save(ccl_model_nested_without_sister_coef,file="ccl_model_nested_without_sister_coef.RData")


delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%ccl_model_nested_without_sister_coef

ccl_model_nested_without_sister_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

save(ccl_model_nested_without_sister_R2,file="ccl_model_nested_without_sister_R2.RData")


### block descent 



### full model : block descendant (z)

cl<-makeCluster(10,type="FORK")
registerDoParallel(cl)

ccl_model_block_descendant_coef<-sapply(1:dim(X_mo)[2],function(i){
	print(i)
	motherPairInd<-which(colnames(X_mo)==colnames(X_mo)[i])
	daughterPairInd<-which(colnames(X_mo)==colnames(X_da1)[i])
	dupInd<-which(duplicated(colnames(X_mo)))
	
	focal_lineage<-do.call(rbind,lapply(which(lineage_info==unlist(strsplit(colnames(X_da1)[i],"_ccl"))),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2]))
	#including itself
	descendant<-unique(unlist(lapply(1:dim(focal_lineage)[1],function(j){lineage_info[focal_lineage[j,1]:10,focal_lineage[j,2]]})))
	descendantInd<-unlist(lapply(paste0(descendant,"_ccl"),function(x){which(colnames(X_mo)==x)}))

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


colnames(ccl_model_block_descendant_coef)<-colnames(X_da1)


save(ccl_model_block_descendant_coef,file="ccl_model_block_descendant_coef.RData")


delta_X_pred<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%ccl_model_block_descendant_coef

ccl_model_block_descendant_R2<-sapply(1:dim(delta_X)[2],function(i){cor(delta_X[,i],delta_X_pred[,i])^2})

save(ccl_model_block_descendant_R2,file="ccl_model_block_descendant_R2.RData")


#######################################################













