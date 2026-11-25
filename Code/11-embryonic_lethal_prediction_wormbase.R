
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


#####################################################################################

#####################################################################################

library(ggplot2)
library(glmnet)



focal_projectFolderPath<-r"(E:\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/data/prediction"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling_noise.RData"
	),projectFolderPath=focal_projectFolderPath)





#############################################################################################################################
							wormbase pheno
#############################################################################################################################

#wormbase_pheno_colname<-sapply(read.table("E:/papers/canalization/data/wormbase_pheno/caenorhabditis_elegans.PRJNA13758.WBPS17.phenotypes.gaf",header=F,sep="\t")[7:19,1],
#	function(x){unlist(strsplit(x,") "))[2]})


wormbase_pheno_colname<-c("DB","Gene_ID","Gene_Symbol","Qualifier","Phenotype_ID","DB_Reference","Evidence_Code","Aspect","Object_Type","Taxon","Release_Date","Assigned_By","Annotation_Extension")

wormbase_pheno<-read.table(paste0(focal_projectFolderPath,"/wormbase/phenotype/caenorhabditis_elegans.PRJNA13758.WBPS18.phenotypes.wb"),header=F,skip=19,sep="\t")
colnames(wormbase_pheno)<-wormbase_pheno_colname

#wormbase_pheno<-wormbase_pheno[wormbase_pheno$Qualifier!="NOT",]

#

xyz_rowGene<-tolower(do.call(rbind,strsplit(rownames(q1_cell350_alignment_xyz_noise),"_"))[,1])

#

wormbase_pheno_mat<-matrix(nrow=dim(q1_cell350_alignment_xyz_noise)[1],ncol=length(table(wormbase_pheno$Phenotype_ID)))

for(i in 1:dim(wormbase_pheno_mat)[2]){
	print(i)
	focal_pheno<-names(table(wormbase_pheno$Phenotype_ID))[i]
	focal_pheno_1_ind<-sort(unlist(lapply(unique(wormbase_pheno$Gene_Symbol[wormbase_pheno$Phenotype_ID==focal_pheno&wormbase_pheno$Qualifier!="NOT"]),function(x){which(xyz_rowGene==x)})))
	wormbase_pheno_mat[focal_pheno_1_ind,i]<-"1"
	focal_pheno_1_ind<-sort(unlist(lapply(unique(wormbase_pheno$Gene_Symbol[wormbase_pheno$Phenotype_ID==focal_pheno&wormbase_pheno$Qualifier=="NOT"]),function(x){which(xyz_rowGene==x)})))
	wormbase_pheno_mat[focal_pheno_1_ind,i]<-"0"	
}

colnames(wormbase_pheno_mat)<-names(table(wormbase_pheno$Phenotype_ID))
rownames(wormbase_pheno_mat)<-rownames(q1_cell350_alignment_xyz_noise)

wormbase_pheno_mat[grep("ctr_",rownames(q1_cell350_alignment_xyz_noise)),]<-"0"

plot(apply(wormbase_pheno_mat,2,function(x){length(which(x=="1"))}),apply(wormbase_pheno_mat,2,function(x){length(which(x=="0"))}))

save(wormbase_pheno_mat,file=paste0(focal_projectFolderPath,"/data/prediction/wormbase_pheno_mat.RData"))

#

wormbase_pheno_mat_select<-wormbase_pheno_mat[,apply(wormbase_pheno_mat,2,function(x){length(which(x=="1"))})>=80]

save(wormbase_pheno_mat_select,file=paste0(focal_projectFolderPath,"/data/prediction/wormbase_pheno_mat_select.RData"))



#######################################################################
					featuer matrix 
#######################################################################


### xyz

##

q1_cell350_alignment_xyz_noise_x<-q1_cell350_alignment_xyz_noise[,1:714]
q1_cell350_alignment_xyz_noise_y<-q1_cell350_alignment_xyz_noise[,(1+714):(714*2)]
q1_cell350_alignment_xyz_noise_z<-q1_cell350_alignment_xyz_noise[,(1+714*2):(714*3)]

#

feature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise[,1:714]),"_x"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]

##x 

X_mo<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo

x_deviation<-sapply(1:dim(delta_X)[2],function(i){return(residuals(lm(delta_X[,i]~X_mo[,i])))})
colnames(x_deviation)<-colnames(delta_X)

##y 

X_mo<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo

y_deviation<-sapply(1:dim(delta_X)[2],function(i){return(residuals(lm(delta_X[,i]~X_mo[,i])))})
colnames(y_deviation)<-colnames(delta_X)

##z 

X_mo<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo

z_deviation<-sapply(1:dim(delta_X)[2],function(i){return(residuals(lm(delta_X[,i]~X_mo[,i])))})
colnames(z_deviation)<-colnames(delta_X)



##################

feature_names_CCL<-unlist(strsplit(colnames(ccl_raw_nonNA_scaling_noise),"_ccl"))

temp1<-unlist(lapply(feature_names_CCL,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names_CCL,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names_CCL,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair_CCL<-cell_mo_da1_da2[intersect(temp1,temp2),]

## feature : CCL

X_mo<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair_CCL$mo,feature_names_CCL)]
X_da1<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair_CCL$da1,feature_names_CCL)]
X_da2<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair_CCL$da2,feature_names_CCL)]

delta_X<-X_da1-X_mo

CCL_deviation<-sapply(1:dim(delta_X)[2],function(i){return(residuals(lm(delta_X[,i]~X_mo[,i])))})
colnames(CCL_deviation)<-colnames(delta_X)

temp<-ccl_raw_nonNA_scaling_noise
temp[,colnames(CCL_deviation)]<-CCL_deviation
CCL_deviation<-temp


###

residual_feature<-cbind(x_deviation,y_deviation,z_deviation)

save(residual_feature,file="residual_feature.RData")

residual_feature_comb<-cbind(residual_feature,CCL_deviation)

save(residual_feature_comb,file="residual_feature_comb.RData")

#################################################################################
						prediction 
#################################################################################


#

cvglm_wormbase<-list()

for(i in 1:dim(wormbase_pheno_mat_select)[2]){
	print(i)
	
	nonNAInd<-which(!is.na(wormbase_pheno_mat_select[,i]))
	
	focalX<-residual_feature_comb[nonNAInd,]
	focalY<-wormbase_pheno_mat_select[nonNAInd,i]
	
	set.seed(1)
	cvglm_wormbase[[i]]<-cv.glmnet(x=focalX,y=focalY,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE)
	
	print(cvglm_wormbase[[i]]$cvm[cvglm_wormbase[[i]]$lambda==cvglm_wormbase[[i]]$lambda.min])
}

names(cvglm_wormbase)<-colnames(wormbase_pheno_mat_select)

save(cvglm_wormbase,file="cvglm_wormbase.RData")



cvglm_wormbase_auc<-sapply(cvglm_wormbase,function(x){x$cvm[x$lambda==x$lambda.min]})

cvglm_wormbase_rocs<-lapply(1:length(cvglm_wormbase),function(i){roc.glmnet(cvglm_wormbase[[i]]$fit.preval, newy = wormbase_pheno_mat_select[which(!is.na(wormbase_pheno_mat_select[,i])),i])[[cvglm_wormbase[[i]]$index["min",]]]})

plot(cvglm_wormbase_rocs[[1]], type = "l",col="red")
sapply(cvglm_wormbase_rocs, lines, col="grey")
abline(0,1,lty=2)


plot(apply(wormbase_pheno_mat_select,2,function(x){length(which(x=="1"))})/apply(wormbase_pheno_mat_select,2,function(x){length(which(!is.na(x)))}),cvglm_wormbase_auc)




wormbase_summary<-data.frame(phenoID=c("WBPhenotype:0000050",names(table(wormbase_pheno$Phenotype_ID))[apply(wormbase_pheno_mat,2,function(x){length(which(x=="1"))})>=80]),
	auc=c(cvglm$cvm[cvglm$lambda==cvglm$lambda.min],cvglm_wormbase_auc),
	num_of_RNAi_emb=apply(cbind(lethal_pheno,wormbase_pheno_mat_select),2,function(x){length(which(x=="1"))}))
wormbase_summary$annotation=wormbase_pheno$Annotation_Extension[match(wormbase_summary$phenoID,wormbase_pheno$Phenotype_ID)]


###


#

cvglm_wormbase_withoutCCL<-list()

for(i in 1:dim(wormbase_pheno_mat_select)[2]){
	print(i)
	
	nonNAInd<-which(!is.na(wormbase_pheno_mat_select[,i]))
	
	focalX<-residual_feature[nonNAInd,]
	focalY<-wormbase_pheno_mat_select[nonNAInd,i]
	
	set.seed(1)
	cvglm_wormbase_withoutCCL[[i]]<-cv.glmnet(x=focalX,y=focalY,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE)
	
	print(cvglm_wormbase_withoutCCL[[i]]$cvm[cvglm_wormbase_withoutCCL[[i]]$lambda==cvglm_wormbase_withoutCCL[[i]]$lambda.min])
}

names(cvglm_wormbase_withoutCCL)<-colnames(wormbase_pheno_mat_select)

save(cvglm_wormbase_withoutCCL,file="cvglm_wormbase_withoutCCL.RData")




cvglm_wormbase_withoutCCL_auc<-sapply(cvglm_wormbase_withoutCCL,function(x){x$cvm[x$lambda==x$lambda.min]})

cvglm_wormbase_withoutCCL_rocs<-lapply(1:length(cvglm_wormbase_withoutCCL),function(i){roc.glmnet(cvglm_wormbase_withoutCCL[[i]]$fit.preval, newy = wormbase_pheno_mat_select[which(!is.na(wormbase_pheno_mat_select[,i])),i])[[cvglm_wormbase_withoutCCL[[i]]$index["min",]]]})

plot(cvglm_wormbase_withoutCCL_rocs[[1]], type = "l",col="red")
sapply(cvglm_wormbase_withoutCCL_rocs, lines, col="grey")
abline(0,1,lty=2)


plot(apply(wormbase_pheno_mat_select,2,function(x){length(which(x=="1"))})/apply(wormbase_pheno_mat_select,2,function(x){length(which(!is.na(x)))}),cvglm_wormbase_withoutCCL_auc)


plot(cvglm_wormbase_auc,cvglm_wormbase_withoutCCL_auc)
abline(0,1)





