


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

train_partition<-function(embryo_gene,trainProp=2/3){
	train_1<-which(embryo_gene=="ctr")
	embryo_gene_sub<-setdiff(embryo_gene,"ctr")
	embryo_gene_sub_type<-names(table(embryo_gene_sub))
	
	set.seed(1)
	
	train_2_type<-sample(embryo_gene_sub_type,round(length(embryo_gene_sub_type)*trainProp))
	train_2<-unlist(lapply(train_2_type,function(x){which(embryo_gene==x)}))
	
	return(sort(c(train_1,train_2)))
}


#####################################################################################

#####################################################################################

library(ggplot2)
library(glmnet)



focal_projectFolderPath<-r"(F:\wangj\Desktop\temp\canalizationENV)"

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



#######################################################################
					lethal phenotype 
#######################################################################


### embryo lethal phenotype 

embryo_lethal_RNAi<-read.csv(paste0(focal_projectFolderPath,"/data/prediction/_download_Supplemental tables_Table_S3_1.csv"))
embryo_lethal_Env<-read.table(paste0(focal_projectFolderPath,"/data/prediction/embryo_state_Env.txt"),header=T)

embryo_lethal<-rep("0",dim(q1_cell350_alignment_xyz_noise)[1])

embryo_lethal[na.omit(match(embryo_lethal_RNAi$Embryo.ID,rownames(q1_cell350_alignment_xyz_noise)))]<-embryo_lethal_RNAi$Organismal.phenotype[which(!is.na(match(embryo_lethal_RNAi$Embryo.ID,rownames(q1_cell350_alignment_xyz_noise))))]

embryo_lethal[na.omit(match(embryo_lethal_Env$embryo,rownames(q1_cell350_alignment_xyz_noise)))]<-embryo_lethal_Env$state[which(!is.na(match(embryo_lethal_Env$embryo,rownames(q1_cell350_alignment_xyz_noise))))]

embryo_lethal[embryo_lethal=="alive"|embryo_lethal=="hatched"]<-"0"
embryo_lethal[embryo_lethal!="0"]<-"1"

save(embryo_lethal,file=paste0(focal_projectFolderPath,"/data/prediction/embryo_lethal.RData"))

###

embryo_lethal<-as.numeric(embryo_lethal)

embryoGene<-c(do.call(rbind,strsplit(rownames(q1_cell350_alignment_xyz_noise)[1:1928],"_"))[,1],do.call(rbind,strsplit(rownames(q1_cell350_alignment_xyz_noise)[1929:2039],"_"))[,3])

trainIndex<-train_partition(embryoGene)

embryo_lethal_train<-embryo_lethal[trainIndex]

### shuffle 

set.seed(1)
embryo_lethal_shuffle<-embryo_lethal[sample(length(embryo_lethal))]

###########################################################################
					phenotype prediction by xyz : residual feature matrix 
###########################################################################


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

## feature : x 

X_mo<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo

x_deviation<-sapply(1:dim(delta_X)[2],function(i){return(residuals(lm(delta_X[,i]~X_mo[,i])))})
colnames(x_deviation)<-colnames(delta_X)

temp<-q1_cell350_alignment_xyz_noise_x
temp[,colnames(x_deviation)]<-x_deviation
x_deviation<-temp

## feature : y 

X_mo<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo

y_deviation<-sapply(1:dim(delta_X)[2],function(i){return(residuals(lm(delta_X[,i]~X_mo[,i])))})
colnames(y_deviation)<-colnames(delta_X)

temp<-q1_cell350_alignment_xyz_noise_y
temp[,colnames(y_deviation)]<-y_deviation
y_deviation<-temp

## feature : z 

X_mo<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo

z_deviation<-sapply(1:dim(delta_X)[2],function(i){return(residuals(lm(delta_X[,i]~X_mo[,i])))})
colnames(z_deviation)<-colnames(delta_X)

temp<-q1_cell350_alignment_xyz_noise_z
temp[,colnames(z_deviation)]<-z_deviation
z_deviation<-temp

###

residual_feature<-cbind(x_deviation,y_deviation,z_deviation)

save(residual_feature,file="residual_feature.RData")

### training-testing 
residual_feature_train<-residual_feature[trainIndex,]


#######################################################################
					phenotype prediction by xyz : prediction 
#######################################################################

### total

set.seed(1)
cvglm_residual<-cv.glmnet(x=residual_feature,y=embryo_lethal,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE)

save(cvglm_residual,file="cvglm_residual.RData")

cvglm_residual_foldid<-cvglm_residual$foldid 
save(cvglm_residual_foldid,file="cvglm_residual_foldid.RData")

### training-testing 

set.seed(1)
cvglm_residual_training<-cv.glmnet(x=residual_feature_train,y=embryo_lethal_train,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE)

save(cvglm_residual_training,file="cvglm_residual_training.RData")


### generation constrained

cell_name<-do.call(rbind,strsplit(colnames(q1_cell350_alignment_xyz_noise)[1:714],"_"))[,1]
cell_gen<-sapply(match(cell_name,lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[1,]


#gen=9

gen=9
selectInd<-c(which(cell_gen<=gen),which(cell_gen<=gen)+714,which(cell_gen<=gen)+714*2)

set.seed(1)
cvglm_residual_g9<-cv.glmnet(x=residual_feature[,selectInd],y=embryo_lethal,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE,foldid=cvglm_residual$foldid)

#gen=8

gen=8
selectInd<-c(which(cell_gen<=gen),which(cell_gen<=gen)+714,which(cell_gen<=gen)+714*2)

set.seed(1)
cvglm_residual_g8<-cv.glmnet(x=residual_feature[,selectInd],y=embryo_lethal,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE,foldid=cvglm_residual$foldid)

#gen=7

gen=7
selectInd<-c(which(cell_gen<=gen),which(cell_gen<=gen)+714,which(cell_gen<=gen)+714*2)

set.seed(1)
cvglm_residual_g7<-cv.glmnet(x=residual_feature[,selectInd],y=embryo_lethal,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE,foldid=cvglm_residual$foldid)

#gen=6

gen=6
selectInd<-c(which(cell_gen<=gen),which(cell_gen<=gen)+714,which(cell_gen<=gen)+714*2)

set.seed(1)
cvglm_residual_g6<-cv.glmnet(x=residual_feature[,selectInd],y=embryo_lethal,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE,foldid=cvglm_residual$foldid)

#gen=5

gen=5
selectInd<-c(which(cell_gen<=gen),which(cell_gen<=gen)+714,which(cell_gen<=gen)+714*2)

set.seed(1)
cvglm_residual_g5<-cv.glmnet(x=residual_feature[,selectInd],y=embryo_lethal,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE,foldid=cvglm_residual$foldid)

#gen=4

gen=4
selectInd<-c(which(cell_gen<=gen),which(cell_gen<=gen)+714,which(cell_gen<=gen)+714*2)

set.seed(1)
cvglm_residual_g4<-cv.glmnet(x=residual_feature[,selectInd],y=embryo_lethal,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE,foldid=cvglm_residual$foldid)

save(cvglm_residual_g9,file="cvglm_residual_g9.RData")
save(cvglm_residual_g8,file="cvglm_residual_g8.RData")
save(cvglm_residual_g7,file="cvglm_residual_g7.RData")
save(cvglm_residual_g6,file="cvglm_residual_g6.RData")
save(cvglm_residual_g5,file="cvglm_residual_g5.RData")
save(cvglm_residual_g4,file="cvglm_residual_g4.RData")

### shuflle 


set.seed(1)
cvglm_residual_shuffle<-cv.glmnet(x=residual_feature,y=embryo_lethal_shuffle,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE,foldid=cvglm_residual$foldid)

save(cvglm_residual_shuffle,file="cvglm_residual_shuffle.RData")




#######################################################################
					phenotype prediction by CCL
#######################################################################


#

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

residual_feature_comb<-cbind(residual_feature,CCL_deviation)

### training-testing 

residual_feature_comb_train<-residual_feature_comb[trainIndex,]
CCL_deviation_train<-CCL_deviation[trainIndex,]

#######################################################################
					phenotype prediction by CCL : prediction 
#######################################################################

### total

penaltyFactor<-rep(1,dim(residual_feature_comb)[2])
penaltyFactor[which(coef(cvglm_residual)[-1]!=0)]<-0

set.seed(1)
cvglm_residual_comb<-cv.glmnet(x=residual_feature_comb,y=embryo_lethal,penalty.factor=penaltyFactor,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE,foldid=cvglm_residual$foldid)

save(cvglm_residual_comb,file="cvglm_residual_comb.RData")


### traing-testing 

set.seed(1)
cvglm_residual_comb_training<-cv.glmnet(x=residual_feature_comb_train,y=embryo_lethal_train,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE)

save(cvglm_residual_comb_training,file="cvglm_residual_comb_training.RData")


### only CCL 

set.seed(1)
cvglm_CCL<-cv.glmnet(x=CCL_deviation,y=embryo_lethal,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE,foldid=cvglm_residual$foldid)

save(cvglm_CCL,file="cvglm_CCL.RData")


### only CCL : shuffled

set.seed(4)
cvglm_CCL_shuffle<-cv.glmnet(x=CCL_deviation,y=embryo_lethal_shuffle,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE)

cvglm_CCL_shuffle$cvm[which(cvglm_CCL_shuffle$lambda==cvglm_CCL_shuffle$lambda.min)]
cvglm_CCL_shuffle$cvsd[which(cvglm_CCL_shuffle$lambda==cvglm_CCL_shuffle$lambda.min)]
plot(roc.glmnet(cvglm_CCL_shuffle$fit.preval, newy = embryo_lethal_shuffle,family="binomial")[[which(cvglm_CCL_shuffle$lambda==cvglm_CCL_shuffle$lambda.min)]])
abline(0,1)


save(cvglm_CCL_shuffle,file="cvglm_CCL_shuffle.RData")


### only CCL : training 

set.seed(1)
cvglm_CCL_train<-cv.glmnet(x=CCL_deviation_train,y=embryo_lethal_train,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE)

save(cvglm_CCL_train,file="cvglm_CCL_train.RData")



####################################################################################################
							End
####################################################################################################

#cvglm_residual_abs<-cv.glmnet(x=abs(residual_feature),y=embryo_lethal,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE)
#cvglm_noise<-cv.glmnet(x=q1_cell350_alignment_xyz_noise,y=embryo_lethal,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE)
#cvglm_noise_abs<-cv.glmnet(x=abs(q1_cell350_alignment_xyz_noise),y=embryo_lethal,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE)


cvglm<-cvglm_residual
cvglm$cvm[cvglm$lambda==cvglm$lambda.min]

cvglm<-cvglm_residual_abs
cvglm$cvm[cvglm$lambda==cvglm$lambda.min]

cvglm<-cvglm_noise
cvglm$cvm[cvglm$lambda==cvglm$lambda.min]

cvglm<-cvglm_noise_abs
cvglm$cvm[cvglm$lambda==cvglm$lambda.min]


cvglm_auc<-cvglm$cvm[cvglm$lambda==cvglm$lambda.min]

cvglm_rocs<-roc.glmnet(cvglm$fit.preval, newy = embryo_lethal)[[cvglm$lambda==cvglm$lambda.min]]

plot(cvglm_rocs, type = "l",col="red")
abline(0,1,lty=2)


### prediction

embryo_names<-do.call(rbind,strsplit(do.call(rbind,strsplit(rownames(residual_feature),"_emb"))[,1],"_s"))[,1]
asso_coef_name<-names(coefficients(cvglm)[which(coefficients(cvglm)!=0),])

prediction<-cvglm$fit.preval[,which(cvglm$lambda==cvglm$lambda.min)]

predion_summary<-data.frame(embryo=embryo_names,lethal=embryo_lethal,prediction)

lethal_alive_diff<-sapply(names(table(predion_summary[,1])),function(x){
	focal_mut<-predion_summary[which(predion_summary[,1]==x),]
	focal_lethal_ind<-which(focal_mut[,2]==1)
	focal_alive_ind<-which(focal_mut[,2]==0)
	if(length(focal_lethal_ind)>0&length(focal_alive_ind)>0){
		return(max(focal_mut[focal_lethal_ind,3])-min(focal_mut[focal_alive_ind,3]))
	}else{
		return(NA)
	}
	
})

sort(lethal_alive_diff)

#

predion_summary[predion_summary[,1]=="SEC-12",]

focal_mut_ind<-which(embryo_names=="SEC-12")


dataTemp2<-data.frame(NO=rep(c(1:length(asso_coef_name[-1])),3),
	cell_feature_name=factor(rep(asso_coef_name[-1],3),level=asso_coef_name[-1]),
	cell_feature=as.numeric(t(q1_cell350_alignment_xyz_noise[which(embryo_names=="SEC-12"),match(asso_coef_name[-1],colnames(q1_cell350_alignment_xyz_noise))])),
	embryo=as.character(sapply(paste0("embryo_",1:length(focal_mut_ind)),rep,length(asso_coef_name[-1]))),
	generation=rep(sapply(match(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,1],lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[1,],3),
	lineage_type=rep(lineage_info[3,sapply(match(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,1],lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[2,]],3),
	lineage=rep(sapply(match(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,1],lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[2,],3))


ggplot(dataTemp2,aes(x=cell_feature_name,y=cell_feature,col=embryo))+
geom_point()+
geom_line(aes(x=NO,y=cell_feature,col=embryo))+
geom_hline(yintercept=0,lty=2)+
geom_vline(xintercept=c(cumsum(table(sapply(match(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,1],lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[1,]))),lty=2)
#facet_wrap(vars(lineage_type),nrow=1,scales="fixed",strip.position = "top")


focal_generation_1<-sapply(match(unique(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,1]),lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[1,]
focal_lineage_1<-unique(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,1])[which(focal_generation_1==10)]

temp<-setdiff(unique(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,1]),as.character(lineage_info[,sapply(match(focal_lineage_1,lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[2,]]))
focal_generation_2<-sapply(match(temp,lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[1,]
focal_lineage_2<-temp[which(focal_generation_2==9)]

temp<-setdiff(temp,as.character(lineage_info[,sapply(match(focal_lineage_2,lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[2,]]))
focal_generation_3<-sapply(match(temp,lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[1,]
focal_lineage_3<-temp[which(focal_generation_3==8)]

temp<-setdiff(temp,as.character(lineage_info[,sapply(match(focal_lineage_3,lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[2,]]))
focal_generation_4<-sapply(match(temp,lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[1,]
focal_lineage_4<-temp[which(focal_generation_4==7)]

temp<-setdiff(temp,as.character(lineage_info[,sapply(match(focal_lineage_4,lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[2,]]))


focal_lineage<-colnames(lineage_info)[sapply(match(c(focal_lineage_1,focal_lineage_2,focal_lineage_3,focal_lineage_4),lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[2,]]


dataTemp3<-data.frame(cell=as.character(lineage_info[4:10,match(focal_lineage,colnames(lineage_info))]))
dataTemp3$noise_x<-q1_cell350_alignment_xyz_noise_x[focal_mut_ind[2],match(dataTemp3$cell,do.call(rbind,strsplit(colnames(q1_cell350_alignment_xyz_noise_x),"_"))[,1])]
dataTemp3$noise_y<-q1_cell350_alignment_xyz_noise_y[focal_mut_ind[2],match(dataTemp3$cell,do.call(rbind,strsplit(colnames(q1_cell350_alignment_xyz_noise_y),"_"))[,1])]
dataTemp3$noise_z<-q1_cell350_alignment_xyz_noise_z[focal_mut_ind[2],match(dataTemp3$cell,do.call(rbind,strsplit(colnames(q1_cell350_alignment_xyz_noise_z),"_"))[,1])]

dataTemp3$generation<-sapply(match(dataTemp3$cell,lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[1,]
dataTemp3$lineage<-as.numeric(sapply(match(focal_lineage,colnames(lineage_info)),rep,7))
dataTemp3<-na.omit(dataTemp3)

dataTemp3<-dataTemp3[order(dataTemp3$lineage,dataTemp3$generation),]
dataTemp3$cell<-factor(dataTemp3$cell,level=unique(dataTemp3$cell))
dataTemp3$lineage<-factor(dataTemp3$lineage,level=unique(dataTemp3$lineage))

dataTemp3$sig_x<-rep(1,dim(dataTemp3)[1])
dataTemp3$sig_x[match(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,2]=="x",1],dataTemp3$cell)]<-2
dataTemp3$sig_y<-rep(1,dim(dataTemp3)[1])
dataTemp3$sig_y[match(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,2]=="y",1],dataTemp3$cell)]<-2
dataTemp3$sig_z<-rep(1,dim(dataTemp3)[1])
dataTemp3$sig_z[match(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,2]=="z",1],dataTemp3$cell)]<-2



ggplot(dataTemp3,aes(x=generation,y=noise_x,col=lineage))+geom_point(show.legend=F,size=dataTemp3$sig_x)+geom_line(show.legend=F)+geom_hline(yintercept=0,lty=2)
ggplot(dataTemp3,aes(x=generation,y=noise_y,col=lineage))+geom_point(show.legend=F)+geom_line(show.legend=F)+geom_hline(yintercept=0,lty=2)
ggplot(dataTemp3,aes(x=generation,y=noise_z,col=lineage))+geom_point(show.legend=F)+geom_line(show.legend=F)+geom_hline(yintercept=0,lty=2)



	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
#############################################################################################################################
							wormbase pheno
#############################################################################################################################

#wormbase_pheno_colname<-sapply(read.table("E:/papers/canalization/data/wormbase_pheno/caenorhabditis_elegans.PRJNA13758.WBPS17.phenotypes.gaf",header=F,sep="\t")[7:19,1],
#	function(x){unlist(strsplit(x,") "))[2]})


wormbase_pheno_colname<-c("DB","Gene_ID","Gene_Symbol","Qualifier","Phenotype_ID","DB_Reference","Evidence_Code","Aspect","Object_Type","Taxon","Release_Date","Assigned_By","Annotation_Extension")

wormbase_pheno<-read.table(paste0(focal_projectFolderPath,"/data/prediction/caenorhabditis_elegans.PRJNA13758.WBPS17.phenotypes.gaf"),header=F,skip=19,sep="\t")
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
					phenotype prediction 
#######################################################################


###

residual_feature<-cbind(x_deviation,y_deviation,z_deviation)

celegans_pheno<-cbind(embryo_lethal,wormbase_pheno_mat_select)




cvglm_wormbase<-list()

for(i in 1:dim(celegans_pheno)[2]){
	print(i)
	
	nonNAInd<-which(!is.na(celegans_pheno[,i]))
	
	focal_y<-celegans_pheno[nonNAInd,i]

	cvglm_wormbase[[i]]<-cv.glmnet(x=residual_feature[nonNAInd,],y=focal_y,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE)
	
	print(cvglm_wormbase[[i]]$cvm[cvglm_wormbase[[i]]$lambda==cvglm_wormbase[[i]]$lambda.min])
}


cvglm_wormbase_auc<-sapply(cvglm_wormbase,function(x){x$cvm[x$lambda==x$lambda.min]})

cvglm_wormbase_rocs<-lapply(1:length(cvglm_wormbase),function(i){roc.glmnet(cvglm_wormbase[[i]]$fit.preval, newy = celegans_pheno[which(!is.na(celegans_pheno[,i])),i])[[cvglm_wormbase[[i]]$index["min",]]]})

plot(cvglm_wormbase_rocs[[1]], type = "l",col="red")
sapply(cvglm_wormbase_rocs, lines, col="grey")
abline(0,1,lty=2)
text(0.7,0.3,"red:embryonic lethal (Du Zhuo)")
text(0.7,0.2,"grey:wormbase traits (#RNAi_emb>80)")


plot(apply(celegans_pheno,2,function(x){length(which(x=="1"))})/apply(celegans_pheno,2,function(x){length(which(!is.na(x)))}),cvglm_wormbase_auc)

### prediction

cvglm_wormbase[[1]]$fit.preval[,which(cvglm_wormbase[[1]]$lambda==cvglm_wormbase[[1]]$lambda.min)]

dataTemp<-data.frame(embryo=do.call(rbind,strsplit(do.call(rbind,strsplit(rownames(residual_feature),"_emb"))[,1],"_s"))[,1],
lethal=as.numeric(celegans_pheno[,1]),prediction=cvglm_wormbase[[1]]$fit.preval[,which(cvglm_wormbase[[1]]$lambda==cvglm_wormbase[[1]]$lambda.min)])

lethal_alive_diff<-sapply(names(table(dataTemp[,1])),function(x){
	focal_mut<-dataTemp[which(dataTemp[,1]==x),]
	focal_lethal_ind<-which(focal_mut[,2]==1)
	focal_alive_ind<-which(focal_mut[,2]==0)
	if(length(focal_lethal_ind)>0&length(focal_alive_ind)>0){
		return(max(focal_mut[focal_lethal_ind,3])-min(focal_mut[focal_alive_ind,3]))
	}else{
		return(NA)
	}
	
})

sort(lethal_alive_diff)

dataTemp[dataTemp[,1]=="SEC-12",]

names(coefficients(cvglm_wormbase[[1]])[which(coefficients(cvglm_wormbase[[1]])!=0),])

embryo_names<-do.call(rbind,strsplit(do.call(rbind,strsplit(rownames(residual_feature),"_emb"))[,1],"_s"))[,1]
asso_coef_name<-names(coefficients(cvglm_wormbase[[1]])[which(coefficients(cvglm_wormbase[[1]])!=0),])
focal_mut_ind<-which(do.call(rbind,strsplit(do.call(rbind,strsplit(rownames(residual_feature),"_emb"))[,1],"_s"))[,1]=="SEC-12")


dataTemp2<-data.frame(NO=rep(c(1:length(asso_coef_name[-1])),3),
	cell_feature_name=factor(rep(asso_coef_name[-1],3),level=asso_coef_name[-1]),
	cell_feature=as.numeric(t(q1_cell350_alignment_xyz_noise[which(embryo_names=="SEC-12"),match(asso_coef_name[-1],colnames(q1_cell350_alignment_xyz_noise))])),
	embryo=as.character(sapply(paste0("embryo_",1:length(focal_mut_ind)),rep,length(asso_coef_name[-1]))),
	generation=rep(sapply(match(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,1],lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[1,],3),
	lineage_type=rep(lineage_info[3,sapply(match(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,1],lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[2,]],3),
	lineage=rep(sapply(match(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,1],lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[2,],3))


ggplot(dataTemp2,aes(x=cell_feature_name,y=cell_feature,col=embryo))+
geom_point()+
geom_line(aes(x=NO,y=cell_feature,col=embryo))+
geom_hline(yintercept=0,lty=2)+
geom_vline(xintercept=c(cumsum(table(sapply(match(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,1],lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[1,]))),lty=2)
#facet_wrap(vars(lineage_type),nrow=1,scales="fixed",strip.position = "top")


focal_generation_1<-sapply(match(unique(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,1]),lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[1,]
focal_lineage_1<-unique(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,1])[which(focal_generation_1==10)]

temp<-setdiff(unique(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,1]),as.character(lineage_info[,sapply(match(focal_lineage_1,lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[2,]]))
focal_generation_2<-sapply(match(temp,lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[1,]
focal_lineage_2<-temp[which(focal_generation_2==9)]

temp<-setdiff(temp,as.character(lineage_info[,sapply(match(focal_lineage_2,lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[2,]]))
focal_generation_3<-sapply(match(temp,lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[1,]
focal_lineage_3<-temp[which(focal_generation_3==8)]

temp<-setdiff(temp,as.character(lineage_info[,sapply(match(focal_lineage_3,lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[2,]]))
focal_generation_4<-sapply(match(temp,lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[1,]
focal_lineage_4<-temp[which(focal_generation_4==7)]

temp<-setdiff(temp,as.character(lineage_info[,sapply(match(focal_lineage_4,lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[2,]]))


focal_lineage<-colnames(lineage_info)[sapply(match(c(focal_lineage_1,focal_lineage_2,focal_lineage_3,focal_lineage_4),lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[2,]]


dataTemp3<-data.frame(cell=as.character(lineage_info[4:10,match(focal_lineage,colnames(lineage_info))]))
dataTemp3$noise_x<-q1_cell350_alignment_xyz_noise_x[focal_mut_ind[2],match(dataTemp3$cell,do.call(rbind,strsplit(colnames(q1_cell350_alignment_xyz_noise_x),"_"))[,1])]
dataTemp3$noise_y<-q1_cell350_alignment_xyz_noise_y[focal_mut_ind[2],match(dataTemp3$cell,do.call(rbind,strsplit(colnames(q1_cell350_alignment_xyz_noise_y),"_"))[,1])]
dataTemp3$noise_z<-q1_cell350_alignment_xyz_noise_z[focal_mut_ind[2],match(dataTemp3$cell,do.call(rbind,strsplit(colnames(q1_cell350_alignment_xyz_noise_z),"_"))[,1])]

dataTemp3$generation<-sapply(match(dataTemp3$cell,lineage_info),rowColRec,dim(lineage_info)[1],dim(lineage_info)[2])[1,]
dataTemp3$lineage<-as.numeric(sapply(match(focal_lineage,colnames(lineage_info)),rep,7))
dataTemp3<-na.omit(dataTemp3)

dataTemp3<-dataTemp3[order(dataTemp3$lineage,dataTemp3$generation),]
dataTemp3$cell<-factor(dataTemp3$cell,level=unique(dataTemp3$cell))
dataTemp3$lineage<-factor(dataTemp3$lineage,level=unique(dataTemp3$lineage))

dataTemp3$sig_x<-rep(1,dim(dataTemp3)[1])
dataTemp3$sig_x[match(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,2]=="x",1],dataTemp3$cell)]<-2
dataTemp3$sig_y<-rep(1,dim(dataTemp3)[1])
dataTemp3$sig_y[match(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,2]=="y",1],dataTemp3$cell)]<-2
dataTemp3$sig_z<-rep(1,dim(dataTemp3)[1])
dataTemp3$sig_z[match(do.call(rbind,strsplit(asso_coef_name[-1],"_"))[do.call(rbind,strsplit(asso_coef_name[-1],"_"))[,2]=="z",1],dataTemp3$cell)]<-2



ggplot(dataTemp3,aes(x=generation,y=noise_x,col=lineage))+geom_point(show.legend=F,size=dataTemp3$sig_x)+geom_line(show.legend=F)+geom_hline(yintercept=0,lty=2)
ggplot(dataTemp3,aes(x=generation,y=noise_y,col=lineage))+geom_point(show.legend=F)+geom_line(show.legend=F)+geom_hline(yintercept=0,lty=2)
ggplot(dataTemp3,aes(x=generation,y=noise_z,col=lineage))+geom_point(show.legend=F)+geom_line(show.legend=F)+geom_hline(yintercept=0,lty=2)





prauc<-c()
for(i in 1:21){
logODDs<-cvglm_wormbase[[i]]$fit.preval[,which(cvglm_wormbase[[i]]$lambda==cvglm_wormbase[[i]]$lambda.min)]

preidctProb<-exp(logODDs)/(1+exp(logODDs))

truth<-celegans_pheno[which(!is.na(celegans_pheno[,i])),i]
truth[which(truth=="0")]="Class2"
truth[which(truth=="1")]="Class1"

truth<-factor(truth,level=c("Class1","Class2"))

dataTemp<-data.frame(truth,Class1=preidctProb)

pr_curve(dataTemp, truth, Class1) %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_path() +
  coord_equal() +
  theme_bw()

abline(h=length(which(celegans_pheno[,i]=="1"))/length(which(!is.na(celegans_pheno[,i]))),lty=2,col="red")

prauc[i]<-as.numeric(pr_auc(dataTemp, truth, Class1)[,3])

}

dataTemp<-data.frame(Area_auc=cvglm_wormbase_auc,Area_pr=prauc-sapply(1:21,function(i){length(which(celegans_pheno[,i]=="1"))/length(which(!is.na(celegans_pheno[,i])))}))
plot(dataTemp[,1]-0.5,dataTemp[,2],xlab="area of AUC curve-0.5",ylab="area of PR curve")
abline(0,1,col="red")

plot(sapply(1:21,function(i){length(which(celegans_pheno[,i]=="1"))/length(which(!is.na(celegans_pheno[,i])))}),dataTemp[,2])


abline(lm(dataTemp[,2]~dataTemp[,1]),col="red")

### weight

prauc_weight<-c()

for(j in 1:10){
i=1
	print(j)
	
	nonNAInd<-which(!is.na(celegans_pheno[,i]))
	
	weight1<-rep(1,length(celegans_pheno[,i]))
	weight1[which(celegans_pheno[,i]=="0")]<-0.1*j
	weight1<-weight1[nonNAInd]
	
	focal_y<-celegans_pheno[nonNAInd,i]
	

	cvglm_wormbase_weight<-cv.glmnet(x=residual_feature[nonNAInd,],y=focal_y,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE,weight=weight1)


logODDs<-cvglm_wormbase_weight$fit.preval[,which(cvglm_wormbase_weight$lambda==cvglm_wormbase_weight$lambda.min)]

preidctProb<-exp(logODDs)/(1+exp(logODDs))

truth<-celegans_pheno[which(!is.na(celegans_pheno[,i])),i]
truth[which(truth=="0")]="Class2"
truth[which(truth=="1")]="Class1"

truth<-factor(truth,level=c("Class1","Class2"))

dataTemp<-data.frame(truth,Class1=preidctProb)


prauc_weight[j]<-as.numeric(pr_auc(dataTemp, truth, Class1)[,3])

}



















wormbase_summary<-data.frame(phenoID=c("WBPhenotype:0000050",names(table(wormbase_pheno$Phenotype_ID))[apply(wormbase_pheno_mat,2,function(x){length(which(x=="1"))})>=80]),
	auc=c(cvglm$cvm[cvglm$lambda==cvglm$lambda.min],cvglm_wormbase_auc),
	num_of_RNAi_emb=apply(cbind(lethal_pheno,wormbase_pheno_mat_select),2,function(x){length(which(x=="1"))}))
wormbase_summary$annotation=wormbase_pheno$Annotation_Extension[match(wormbase_summary$phenoID,wormbase_pheno$Phenotype_ID)]
wormbase_summary$annotation[1]<-paste0(wormbase_summary$annotation[1],"(Do Zhuo)")






#################################################################################################

###


load("E:/papers/canalization/data/basic/cell_stage_df.RData")





######


cvglm<-cv.glmnet(x=residual_feature,y=lethal_pheno,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE)

plot(cvglm)

rocs <- roc.glmnet(cvglm$fit.preval, newy = lethal_pheno)
best <- cvglm$index["min",]
plot(rocs[[best]], type = "l")
invisible(sapply(rocs, lines, col="grey"))
lines(rocs[[best]], lwd = 2,col = "red")

dataTemp<-data.frame(coeff=coefficients(cvglm,s="lambda.min")[which(coefficients(cvglm,s="lambda.min")[,1]!=0),][-1],
	cell=sapply(names(coefficients(cvglm,s="lambda.min")[which(coefficients(cvglm,s="lambda.min")[,1]!=0),]),function(x){unlist(strsplit(x,"_"))[1]})[-1],
	coo=sapply(names(coefficients(cvglm,s="lambda.min")[which(coefficients(cvglm,s="lambda.min")[,1]!=0),]),function(x){unlist(strsplit(x,"_"))[2]})[-1])
	
dataTemp$gen=sapply(match(dataTemp$cell,lineage_info),function(x){rowColRec(x,dim(lineage_info)[1],dim(lineage_info)[2])[1]})

ggplot(dataTemp,aes(x=gen,y=coeff,col=coo))+geom_jitter(width=0.25)+
geom_hline(yintercept=0)


cvglm_wormbase<-list()

for(i in 1:dim(wormbase_pheno_mat_select)[2]){
	print(i)
	
	focal_y<-wormbase_pheno_mat_select[,i]

	cvglm_wormbase[[i]]<-cv.glmnet(x=residual_feature,y=focal_y,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE)
	
}

cvglm_wormbase_auc<-sapply(cvglm_wormbase,function(x){x$cvm[x$lambda==x$lambda.min]})

cvglm_wormbase_rocs<-lapply(1:dim(wormbase_pheno_mat_select)[2],function(i){roc.glmnet(cvglm_wormbase[[i]]$fit.preval, newy = wormbase_pheno_mat_select[,i])[[cvglm_wormbase[[i]]$index["min",]]]})

plot(cvglm_wormbase_rocs[[1]])
plot(rocs[[best]], type = "l",col="red")
sapply(cvglm_wormbase_rocs, lines, col="grey")
abline(0,1,lty=2)
text(0.7,0.3,"red:embryonic lethal (Du Zhuo)")
text(0.7,0.2,"grey:wormbase traits (#RNAi_emb>80)")

wormbase_summary<-data.frame(phenoID=c("WBPhenotype:0000050",names(table(wormbase_pheno$Phenotype_ID))[apply(wormbase_pheno_mat,2,function(x){length(which(x=="1"))})>=80]),
	auc=c(cvglm$cvm[cvglm$lambda==cvglm$lambda.min],cvglm_wormbase_auc),
	num_of_RNAi_emb=apply(cbind(lethal_pheno,wormbase_pheno_mat_select),2,function(x){length(which(x=="1"))}))
wormbase_summary$annotation=wormbase_pheno$Annotation_Extension[match(wormbase_summary$phenoID,wormbase_pheno$Phenotype_ID)]
wormbase_summary$annotation[1]<-paste0(wormbase_summary$annotation[1],"(Do Zhuo)")


### testing set 

set.seed(1)
testInd<-sample(1:1926,642)
trainInd<-setdiff(1:1926,testInd)

trainX<-residual_feature[trainInd,]
testX<-residual_feature[testInd,]

trainY<-wormbase_pheno_mat_select[trainInd,1]
testY<-wormbase_pheno_mat_select[testInd,1]

cvglm_wormbase_train_1<-cv.glmnet(x=trainX,y=trainY,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE)

data.frame(testY,predict(cvglm_wormbase_train_1,testX,s="lambda.min"))


####
traitMat<-apply(cbind(lethal_pheno,wormbase_pheno_mat),2,as.numeric)
colnames(traitMat)<-c("lethal_Du",names(table(wormbase_pheno$Phenotype_ID)))

sumstat_list<-list()
for(i in 1:dim(traitMat)[2]){
	print(i)
	sumstat_list[[i]]<-t(apply(residual_feature,2,function(x){summary(glm(traitMat[,i]~x,family="binomial"))$coefficients[2,c(1,2,4)]}))
}

plot(-log10(sumstat_list[[1]][,3]),xlab="xyz features",ylab="-log10(P)")
abline(h=-log10(0.05/dim(sumstat_list[[1]])[1]),col="red")



dataTemp<-data.frame(x=rep(cell_stage_df$x[match(unlist(strsplit(colnames(delta_X),"_x")),rownames(cell_stage_df))],3),
	y=rep(cell_stage_df$y[match(unlist(strsplit(colnames(delta_X),"_x")),rownames(cell_stage_df))],3),
	z=rep(cell_stage_df$z[match(unlist(strsplit(colnames(delta_X),"_x")),rownames(cell_stage_df))],3),
	sumstat_list[[1]],
	coo=c(rep("x",706),rep("y",706),rep("z",706)),
	gen=)

colnames(dataTemp)<-c("x","y","z","logOR","SE","P","coo")

coo_col<-c(rep("#F8766D",706),rep("#00ba38",706),rep("#619CFF",706))
coo_col[which(dataTemp$P>0.05/dim(sumstat_list[[1]])[1])]<-"grey"

ggplot(dataTemp,aes(x=x,y=y,size=-log10(P)))+
geom_point(alpha=0.25,width=0.2,height=0.2,col=coo_col)+




############################################################################################

#######################################################################################################

library(plotROC)
rocdata <- data.frame(D = as.numeric(testY),
                     M = predict(cvglm_wormbase_train_1,testX,s="lambda.min")[,1])

ggroc <- ggplot(rocdata, aes(m = M, d = D)) + geom_roc()
calc_auc(ggroc)


cvglm_result<-data.frame(obs=lethal_pheno,pred=predict(cvglm,residual_feature)[,1])
ggplot(cvglm_result, aes(d = obs, m = pred)) + geom_roc()

library(pROC)

auc(cvglm)











