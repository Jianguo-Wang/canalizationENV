




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


#######################################################################

#######################################################################

library(ggplot2)
library(reshape)
library(gplots)
library(glmnet)


focal_projectFolderPath<-r"(E:\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/figs"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"E:/papers/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/temporal/ccl_model_mother_R.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/temporal/ccl_model_mother_coef.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/temporal/ccl_model_mo_da1_R.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/temporal/ccl_model_mother_coef_shuffle.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/temporal/ccl_model_mo_da1_R_shuffle.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/temporal/ccl_model_mo_da1_pvalue.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/temporal/ccl_model_sister_R.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/temporal/ccl_model_sister_coef.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/temporal/ccl_model_mother_sister_R.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/temporal/ccl_model_mother_sister_coef.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/temporal/ccl_model_full_coef.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/temporal/ccl_model_full_R2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/temporal/ccl_model_block_descendant_coef.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/temporal/ccl_model_block_descendant_R2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/temporal/ccl_model_nested_without_mother_coef.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/temporal/ccl_model_nested_without_sister_coef.RData",			
	"F:/wangj/Desktop/temp/canalizationENV/data/spatial/xyz_x_model_mother_coef.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/spatial/xyz_y_model_mother_coef.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/spatial/xyz_z_model_mother_coef.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/embryo_lethal.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/depth/depth_summary.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/depth/ccl_depth_abs.RData"
	),projectFolderPath=focal_projectFolderPath)


	
###
	
embryo_lethal<-as.numeric(embryo_lethal)

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

temp<-q1_cell350_alignment_xyz_noise_x
temp[,colnames(x_deviation)]<-x_deviation
x_deviation<-temp

##y 

X_mo<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo

y_deviation<-sapply(1:dim(delta_X)[2],function(i){return(residuals(lm(delta_X[,i]~X_mo[,i])))})
colnames(y_deviation)<-colnames(delta_X)

temp<-q1_cell350_alignment_xyz_noise_y
temp[,colnames(y_deviation)]<-y_deviation
y_deviation<-temp

##z 

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

feature_names<-unlist(strsplit(colnames(ccl_raw_nonNA_scaling_noise),"_ccl"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]


X_mo<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$mo,feature_names)]
X_da1<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$da1,feature_names)]

delta_X<-X_da1-X_mo


###

ccl_deviation<-sapply(1:dim(delta_X)[2],function(i){return(residuals(lm(delta_X[,i]~X_mo[,i])))})
colnames(ccl_deviation)<-colnames(delta_X)

temp<-ccl_raw_nonNA_scaling_noise
temp[,colnames(ccl_deviation)]<-ccl_deviation
ccl_deviation<-temp

### only CCL 

residual_feature<-ccl_deviation

# total

set.seed(1)
cvglm_residual_ccl<-cv.glmnet(x=residual_feature,y=embryo_lethal,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE)

#save(cvglm_residual_ccl,file="cvglm_residual_ccl.RData")


# training-testing 

embryoGene<-c(do.call(rbind,strsplit(rownames(ccl_raw_nonNA_scaling_noise)[1:1928],"_"))[,1],do.call(rbind,strsplit(rownames(ccl_raw_nonNA_scaling_noise)[1929:2039],"_"))[,3])

trainIndex<-train_partition(embryoGene)

residual_feature_train<-residual_feature[trainIndex,]
embryo_lethal_train<-embryo_lethal[trainIndex]

residual_feature_test<-residual_feature[setdiff(1:dim(residual_feature)[1],trainIndex),]
embryo_lethal_test<-embryo_lethal[setdiff(1:dim(residual_feature)[1],trainIndex)]

#

set.seed(1)
cvglm_residual_ccl_training<-cv.glmnet(x=residual_feature_train,y=embryo_lethal_train,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE)

#save(cvglm_residual_ccl_training,file="cvglm_residual_ccl_training.RData")

### combine ccl and xyz 

residual_feature_comb<-cbind(ccl_deviation,x_deviation,y_deviation,z_deviation)

#total 
set.seed(1)
cvglm_residual_comb<-cv.glmnet(x=residual_feature_comb,y=embryo_lethal,nlambda=200,type.measure="auc",nfolds = 10,family="binomial",keep=TRUE)

cvglm_residual_comb$cvm[cvglm_residual_comb$lambda==cvglm_residual_comb$lambda.min]


### coefficients of relatives

full_coef_mother<-sapply(1:dim(mo_sis_pair)[1],function(i){
	return(setdiff(ccl_model_full_coef[which(rownames(ccl_model_full_coef)==paste0(mo_sis_pair$mo[i],"_ccl")),i],0))
})


full_coef_sister<-sapply(1:dim(mo_sis_pair)[1],function(i){
	temp<-setdiff(ccl_model_full_coef[which(rownames(ccl_model_full_coef)==paste0(mo_sis_pair$da2[i],"_ccl")),i],0)
	return(ifelse(is.null(temp),NA,temp))
})


full_coef_aunt<-sapply(1:dim(mo_sis_pair)[1],function(i){
	focal_aunt<-mo_sis_pair$da2[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	focal_aunt_ind<-which(rownames(ccl_model_full_coef)==paste0(focal_aunt,"_ccl"))
	if(length(focal_aunt_ind)==0){
		return(NA)
	}else{
		temp<-setdiff(ccl_model_full_coef[focal_aunt_ind,i],0)
		return(ifelse(is.null(temp),NA,temp))	
	}
})


full_coef_grandmother<-sapply(1:dim(mo_sis_pair)[1],function(i){
	focal_grandmother<-mo_sis_pair$mo[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	focal_grandmother_ind<-which(rownames(ccl_model_full_coef)==paste0(focal_grandmother,"_ccl"))
	if(length(focal_grandmother_ind)==0){
		return(NA)
	}else{
		temp<-setdiff(ccl_model_full_coef[focal_grandmother_ind,i],0)
		return(ifelse(is.null(temp),NA,temp))	
	}
})


full_coef_cousin<-unlist(lapply(1:dim(mo_sis_pair)[1],function(i){
	focal_aunt<-mo_sis_pair$da2[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	if(is.na(focal_aunt)){
		return(NA)
	}else{
		tempInd<-which(mo_sis_pair$mo==focal_aunt)
		if(length(tempInd)==0){
			return(NA)
		}else{
			focal_cousin<-unique(c(mo_sis_pair$da1[tempInd],mo_sis_pair$da2[tempInd]))
			focal_cousin_ind<-match(paste0(focal_cousin,"_ccl"),rownames(ccl_model_full_coef))
			if(length(focal_cousin_ind)==0){
				return(NA)
			}else{
				temp<-setdiff(ccl_model_full_coef[focal_cousin_ind,i],0)
				if(is.null(temp)){
					return(NA)
				}else{
					return(temp)
				}
			}			
		}
	}
}))



full_coef_other<-unlist(lapply(1:dim(mo_sis_pair)[1],function(i){

	ind_mo<-which(rownames(ccl_model_full_coef)==paste0(mo_sis_pair$mo[i],"_ccl"))
	ind_sis<-which(rownames(ccl_model_full_coef)==paste0(mo_sis_pair$da2[i],"_ccl"))
	
	focal_aunt<-mo_sis_pair$da2[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	ind_aunt<-which(rownames(ccl_model_full_coef)==paste0(focal_aunt,"_ccl"))
	
	focal_grandmother<-mo_sis_pair$mo[match(mo_sis_pair$mo[i],mo_sis_pair$da1)]
	ind_grandmother<-which(rownames(ccl_model_full_coef)==paste0(focal_grandmother,"_ccl"))	
	
	ind_cousin<-NA
	
	if(!is.na(focal_aunt)){
		tempInd<-which(mo_sis_pair$mo==focal_aunt)
		if(length(tempInd)!=0){
			focal_cousin<-unique(c(mo_sis_pair$da1[tempInd],mo_sis_pair$da2[tempInd]))
			ind_cousin<-unlist(lapply(paste0(focal_cousin,"_ccl"),function(x){which(rownames(ccl_model_full_coef)==x)}))
		}
	}
	
	ind_non_relative<-setdiff(2:346,c(ind_mo,ind_sis,ind_aunt,ind_grandmother,ind_cousin))
	temp<-setdiff(ccl_model_full_coef[ind_non_relative,i],0)
	
	return(temp)

}))

###

coef_df<-data.frame(mo_sis_pair,mo_CCL=diag(ccl_model_full_coef[-1,]),
	sis_CCL=sapply(1:dim(mo_sis_pair)[1],function(i){ccl_model_full_coef[match(paste0(mo_sis_pair[i,3],"_ccl"),rownames(ccl_model_full_coef)),i]})
	)




###########################################################################################



### distribution of kCCL


dataTemp<-data.frame(Coef=ccl_model_mother_coef)


ggplot(dataTemp,aes(x=Coef))+
geom_histogram(colour="grey",fill="grey39",size=0.2)+
geom_vline(xintercept=c(-1,0),lty=2,col="grey",lwd=0.25)+
scale_x_continuous(limits=c(-1.2,0),breaks=c(-1,-0.5,0))+
scale_y_continuous(breaks=c(0,30,60))+
labs(x=expression(italic(k)[CCL]),y="Count",title="")


ggsave("temporal_kCCL.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))




### distribution of p-value for kCCL 

dataTemp<-data.frame(k=ccl_model_mother_coef,
	pvalue=ccl_model_mo_da1_pvalue)

ggplot(dataTemp,aes(x=k,y=-log10(pvalue)))+geom_point(size=0.2,stroke=0.2)+
geom_hline(yintercept=-log10(0.05/345),lty=2,lwd=0.3)+
labs(x=expression(italic(k)[CCL]),y=expression(paste("-",log[10],"(",italic(p),")")),title=expression(paste("78% of ",italic(k)[CCL]," with ",italic(p),"<0.05/345")))+
theme(legend.position=c(0.2,0.8))

length(which(dataTemp$pvalue<0.05/345))/345

ggsave("temporal_k_pvalue.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



### comparison between kCCL and kx, ky and kz 

overlap_pair<-intersect(names(ccl_model_mother_coef),names(xyz_x_model_mother_coef))
ind1<-match(overlap_pair,names(ccl_model_mother_coef))
ind2<-match(overlap_pair,names(xyz_x_model_mother_coef))


dataTemp<-data.frame(kCCL=rep(ccl_model_mother_coef[ind1],3),
	kxyz=c(xyz_x_model_mother_coef[ind2],xyz_y_model_mother_coef[ind2],xyz_z_model_mother_coef[ind2]),
	coo=c(rep("x",length(ind1)),rep("y",length(ind1)),rep("z",length(ind1))))

cor.test(dataTemp[dataTemp$coo=="x",1],dataTemp[dataTemp$coo=="x",2])
cor.test(dataTemp[dataTemp$coo=="y",1],dataTemp[dataTemp$coo=="y",2])
cor.test(dataTemp[dataTemp$coo=="z",1],dataTemp[dataTemp$coo=="z",2])

ggplot(dataTemp,aes(x=kCCL,y=kxyz,col=coo))+
geom_point(size=0.1,stroke=0.1,alpha=0.5,show.legend=FALSE)+
geom_abline(intercept=0,slope=1,lty=2,size=0.3)+
scale_x_continuous(breaks=c(-1,-0.5,0),limits=c(min(dataTemp[,1:2]),max(max(dataTemp[,1:2]),0)))+
scale_y_continuous(breaks=c(-1,-0.5,0),limits=c(min(dataTemp[,1:2]),max(max(dataTemp[,1:2]),0)))+
labs(x=expression(italic(k)[CCL]),y=expression(paste(italic(k)[x]," / ",italic(k)[y]," / ",italic(k)[z])),title="")

ggsave("temporal_kCCL_compare_with_kx_ky_kz.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))




### distribution of depth

dataTemp<-data.frame(depth=ccl_depth_abs[,1],pvalue=ccl_depth_abs[,4])

ggplot(dataTemp,aes(x=depth*100))+
geom_histogram(colour="grey",fill="grey39",size=0.2)+
labs(x=expression(paste(italic(d),"(%)")),y="Count",title="")

ggsave("temporal_depth_distribution.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


# depth and pvalue 

ggplot(dataTemp,aes(x=depth*100,y=-log10(p.adjust(pvalue,method="BH"))))+
geom_point(size=0.2,stroke=0.2,colour="grey39")+
geom_hline(yintercept=-log10(0.05),lwd=0.3,lty=2)+
#scale_x_continuous(breaks=c(-0.12,-0.06,0))+
#scale_y_continuous(breaks=c(0,10,20))+
labs(x=expression(paste(italic(d),"(%)")),y=expression(paste("-",log[10],"(",italic(p),")")),title="")

ggsave("temporal_depth_pvalues.pdf",width=unit(2.1,"in"),height=unit(2.2,"in"))



# compare between depth_abs and slope k 

dataTemp<-data.frame(depth_ccl=ccl_depth_abs[match(mo_sis_pair$da1,feature_names),1]*100,
	k=ccl_model_mother_coef)

ggplot(dataTemp,aes(x=k,y=depth_ccl))+
geom_point(size=0.2,stroke=0.2,col="grey39")+
geom_smooth(method="lm",lwd=0.3)+
scale_x_continuous(limits=c(min(ccl_model_mother_coef),0),breaks=c(-1,-0.5,0))+
scale_y_continuous(breaks=c(10,20,30))+
labs(x=expression(italic(k)),y=expression(paste(italic(d)," (%)")),title="")

ggsave("temporal_depth_Vs_k.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


cor.test(dataTemp[,1],dataTemp[,2])

outlierInd<-union(match(boxplot.stats(dataTemp[,1])$out,dataTemp[,1]),
match(boxplot.stats(dataTemp[,2])$out,dataTemp[,2]))

cor.test(dataTemp[-outlierInd,1],dataTemp[-outlierInd,2])


### comparison between k_Contr, k_RNAi, k_Envir for CCL (main)

k_wt_CCL<-sapply(1:dim(X_mo)[2],function(i){lm(delta_X[1:105,i]~X_mo[1:105,i])$coefficients[2]})
k_mt_CCL<-sapply(1:dim(X_mo)[2],function(i){lm(delta_X[106:1928,i]~X_mo[106:1928,i])$coefficients[2]})
k_env_CCL<-sapply(1:dim(X_mo)[2],function(i){lm(delta_X[1929:2039,i]~X_mo[1929:2039,i])$coefficients[2]})


cor.test(k_wt_CCL,k_mt_CCL)
cor.test(k_wt_CCL,k_env_CCL)

dataTemp<-data.frame(k_wt_CCL,k_mt_env=c(k_mt_CCL,k_env_CCL),mt_env=c(rep("RNAi",345),rep("Env",345)))

ggplot(dataTemp,aes(x=k_wt_CCL,y=k_mt_env,col=mt_env))+
geom_point(size=0.1,stroke=0.1)+
geom_abline(intercept=0,slope=1,lwd=0.3,lty=2)+
scale_x_continuous(limits=c(-1.5,0),breaks=c(-1,-0.5,0))+
scale_y_continuous(limits=c(-1.5,0),breaks=c(-1,-0.5,0))+
labs(x=expression(italic(k)[Control]),y=expression(paste(italic(k)[RNAi]," or ",italic(k)[Environment])),title="")+
theme(legend.position=c(0.1,0.9))

ggsave("temporal_kContr_kRNAi_kENV.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### M-D inheritance recovering under specific environment 




### shuffled cell pairs for CCL 


#R comparison

dataTemp<-data.frame(R=c(ccl_model_mo_da1_R,ccl_model_mo_da1_R_shuffle),type=c(rep("Observed",345),rep("Shuffled",345)))
ggplot(dataTemp,aes(x=R,col=type,fill=type))+
geom_density(alpha=0.5,lwd=0.1)+
geom_vline(xintercept=0,lty=2,lwd=0.25,col="grey")+
scale_y_continuous(limits=c(0,4),breaks=c(0,3))+
labs(x=expression(paste("Pearson's ",italic(R))),y="Density",title="Mother-daughter correlation (CCL)")+
theme(legend.position=c(0.35,0.9))

ggsave("temporal_shuffle_observe_R.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


#k comparison 

dataTemp<-data.frame(kCCL=c(ccl_model_mother_coef,ccl_model_mother_coef_shuffle),
type=c(rep("Observed",345),rep("Shuffled",345)))

ggplot(dataTemp,aes(x=kCCL,col=type,fill=type))+
geom_density(alpha=0.5,lwd=0.1)+
geom_vline(xintercept=c(-1,0),lty=2,lwd=0.25,col=c("#619CFF","#619CFF"))+
scale_y_continuous()+
labs(x=expression(italic(k)[x]),y="Density", title="")+
theme(legend.position=c(0.15,0.9))

ggsave("temporal_shuffle_observe_k.pdf",width=unit(3.2,"in"),height=unit(1.4,"in"))





### compare sister 

uniqueInd<-match(names(table(mo_sis_pair[,1]))[table(mo_sis_pair[,1])==2],mo_sis_pair[,1])

coef_mo_da1<-ccl_model_mother_coef[match(paste0(mo_sis_pair[,1],"_",mo_sis_pair[,2]),names(ccl_model_mother_coef))]
coef_mo_da2<-ccl_model_mother_coef[match(paste0(mo_sis_pair[,1],"_",mo_sis_pair[,3]),names(ccl_model_mother_coef))]

dataTemp<-data.frame(coef_mo_da1=coef_mo_da1[uniqueInd],coef_mo_da2=coef_mo_da2[uniqueInd])


ggplot(dataTemp,aes(x=coef_mo_da1[uniqueInd],y=coef_mo_da2[uniqueInd]))+
geom_point(size=0.5,stroke=0.2,col="grey39")+
geom_abline(slope=1,intercept=0,lwd=0.3,lty=2)+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])),breaks=c(-1,0))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])),breaks=c(-1,0))+
labs(x=expression(italic(k)[Daughter[1]]),y=expression(italic(k)[Daughter[2]]),title="")

ggsave("temporal_kCCL_compare_between_sister.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



### cell noise level 


ccl_wt_mean<-apply(ccl_raw_nonNA_scaling[1:105,],2,mean)

# noise level comparison 

sd_hatch<-sapply(1:dim(ccl_raw_nonNA_scaling[embryo_lethal==0,])[2],function(i){sqrt(sum((ccl_raw_nonNA_scaling[embryo_lethal==0,i]-ccl_wt_mean[i])^2)/length(ccl_raw_nonNA_scaling[embryo_lethal==0,i]))})
sd_lethal<-sapply(1:dim(ccl_raw_nonNA_scaling[embryo_lethal==1,])[2],function(i){sqrt(sum((ccl_raw_nonNA_scaling[embryo_lethal==1,i]-ccl_wt_mean[i])^2)/length(ccl_raw_nonNA_scaling[embryo_lethal==1,i]))})

dataTemp<-data.frame(sd_hatch,sd_lethal)

ggplot(dataTemp,aes(x=sd_hatch,y=sd_lethal))+
geom_point(size=0.1,stroke=0.1,col="grey39")+
geom_abline(slope=1,intercept=0,lwd=0.3,lty=2)+
scale_x_continuous(limits=range(dataTemp),breaks=c(4,8,12))+
scale_y_continuous(limits=range(dataTemp),breaks=c(4,8,12))+
labs(x="Hatched embryos",y="Lethal embryos",title="SD of cell noise")

ggsave("temporal_prediction_compare_SD.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


###

#feature_names2<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise[,1:714]),"_x"))

#temp1<-unlist(lapply(feature_names2,function(x){which(cell_mo_da1_da2$mo==x)}))
#temp2<-unlist(lapply(feature_names2,function(x){which(cell_mo_da1_da2$da1==x)}))
#temp3<-unlist(lapply(feature_names2,function(x){which(cell_mo_da1_da2$da2==x)}))

#mo_sis_pair2<-cell_mo_da1_da2[intersect(temp1,temp2),]


#dataTemp2<-data.frame(depth_x=depth_summary$depth_abs[match(mo_sis_pair2$da1,feature_names2)],
#	depth_y=depth_summary$depth_abs[match(mo_sis_pair2$da1,feature_names2)+714],
#	depth_z=depth_summary$depth_abs[match(mo_sis_pair2$da1,feature_names2)+714*2],
#	kx=xyz_x_model_mother_coef,
#	ky=xyz_y_model_mother_coef,
#	kz=xyz_z_model_mother_coef)

#dataTemp3<-data.frame(dataTemp,dataTemp2[match(names(ccl_model_mother_coef),names(xyz_x_model_mother_coef)),])	

###

# ROC curve and AUC by logistic regression under machine learning framework 

auc_summary<-data.frame(Subset=c("Total","Training","Testing"),AUC=c(
cvglm_residual_ccl$cvm[cvglm_residual_ccl$lambda==cvglm_residual_ccl$lambda.min],
cvglm_residual_ccl_training$cvm[cvglm_residual_ccl_training$lambda==cvglm_residual_ccl_training$lambda.min],
assess.glmnet(predict(cvglm_residual_ccl_training,residual_feature_test,s="lambda.min"), newy = embryo_lethal_test,family="binomial")$auc
))


roc_list<-list()
roc_list[[1]]<-roc.glmnet(cvglm_residual_ccl$fit.preval, newy = embryo_lethal)[[which(cvglm_residual_ccl$lambda==cvglm_residual_ccl$lambda.min)]]
roc_list[[2]]<-roc.glmnet(cvglm_residual_ccl_training$fit.preval, newy = embryo_lethal_train)[[which(cvglm_residual_ccl_training$lambda==cvglm_residual_ccl_training$lambda.min)]]
roc_list[[3]]<-roc.glmnet(predict(cvglm_residual_ccl_training,residual_feature_test,s="lambda.min"), newy = embryo_lethal_test)
names(roc_list)<-c("Total","Training","Testing")


dataTemp<-data.frame(do.call(rbind,roc_list),Subset=factor(c(rep("Total",dim(roc_list[[1]])[1]),rep("Training",dim(roc_list[[2]])[1]),rep("Testing",dim(roc_list[[3]])[1])),level=c("Training","Testing","Total")))

ggplot(dataTemp,aes(x=FPR,y=TPR,col=Subset))+
geom_line()+
geom_abline(slope=1,intercept=0,lty=2)+
theme(legend.position=c(0.5,0.2))+
labs(title="Prediction of embryonic lethality")

#ggsave("temporal_Prediction_of_embryonic_lethality_ROC.pdf",width=unit(2,"in"),height=unit(2.1,"in"))




#########################################################################################
						high-dimensional stability 
#########################################################################################


### compare coefficients of mother and sister 

dataTemp<-data.frame(coef_mo=coef_df$mo_CCL,
	coef_sis=coef_df$sis_CCL)

#x 
	
ggplot(na.omit(dataTemp),aes(x=coef_mo,y=coef_sis))+
geom_point(size=0.5,stroke=0.5,show.legend=FALSE,col=colors()[182])+
geom_smooth(method="lm",show.legend=FALSE)+
labs(x="Coefficients of mother",y="Coefficients of sister",title="CCL")

ggsave("temporal_full_model_coef_mo_sis_comparison.pdf",width=unit(2,"in"),height=unit(2.1,"in"))

### coefficient matrix 

da_cell<-unlist(strsplit(rownames(t(ccl_model_full_coef[-1,])),"_ccl"))
orderInd<-na.omit(match(unlist(lapply(1:9,function(i){names(table(lineage_info[i,]))})),da_cell))

temp<-t(ccl_model_full_coef[-1,])
temp<-temp[orderInd,orderInd]

colnames(temp)<-unlist(strsplit(colnames(temp),"_ccl"))
rownames(temp)<-unlist(strsplit(rownames(temp),"_ccl"))
rownames(temp)<-c(1:345)
colnames(temp)<-c(1:345)


temp_melt<-melt(temp,c("daughter","mother"))

temp_melt<-temp_melt[which(temp_melt$value!=0),]
temp_melt[which(temp_melt[,3]>0),3]<-temp_melt[which(temp_melt[,3]>0),3]*2

ggplot(temp_melt,aes(x=mother,y=daughter,col=value))+geom_point(size=0.2,stroke=0.05)+
scale_colour_gradient2(low="blue",high="red",mid="white",midpoint=0)+
scale_x_continuous(breaks=c(1,345))+
scale_y_continuous(breaks=c(1,345))+
labs(x=expression(italic(L)[M]),y=expression(paste(italic(Delta),italic(L))),title=expression(italic(K)[CCL]))


ggsave("temporal_full_model_coef.pdf",width=unit(2.6,"in"),height=unit(2.1,"in"))



# coefficients of relatives 

dataTemp<-data.frame(
	Coef=c(full_coef_mother,full_coef_sister,full_coef_aunt,full_coef_grandmother,full_coef_cousin,full_coef_other),
	Relative=factor(c(rep("Mother",length(full_coef_mother)),
		rep("Sister",length(full_coef_sister)),
		rep("Aunt",length(full_coef_aunt)),
		rep("Grandmother",length(full_coef_grandmother)),
		rep("Cousin",length(full_coef_cousin)),
		rep("Other",length(full_coef_other))),level=c("Other","Grandmother","Aunt","Cousin","Sister","Mother"))
)


ggplot(dataTemp,aes(y=Relative,x=Coef))+
geom_vline(xintercept=0,col="red",lwd=0.2)+
geom_boxplot(col=colors()[180],fill="grey",show.legend=FALSE,size=0.2,outlier.size=0.2,outlier.stroke=0.1,width=0.5)+
scale_x_continuous(breaks=c(-1,0,0.5))+
labs(x="Coefficients",y="",title="")+
theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.5))


ggsave("temporal_ODE_relative_coefficients.pdf",width=unit(1.7,"in"),height=unit(1.6,"in"))


### R2 comparison - mother and sister

pred_full<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%ccl_model_full_coef
R2_full<-sapply(1:345,function(i){cor(delta_X[,i],pred_full[,i])^2})

R2_mother<-ccl_model_mother_R^2

R2_sister<-ccl_model_sister_R^2

pred_mother_sister<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%ccl_model_mother_sister_coef
R2_mother_sister<-sapply(1:345,function(i){cor(delta_X[,i],pred_mother_sister[,i])^2})


with_sister_ind<-which(!is.na(match(mo_sis_pair$da2[match(unlist(strsplit(colnames(delta_X),"_ccl")),mo_sis_pair$da1)],unlist(strsplit(colnames(X_mo),"_ccl")))))

dataTemp<-data.frame(R2=c(R2_full[with_sister_ind],R2_mother[with_sister_ind],R2_sister[with_sister_ind],R2_mother_sister[with_sister_ind]),
	Model=factor(c(rep("Full",length(with_sister_ind)),
				   rep("Mother",length(with_sister_ind)),
				   rep("Sister",length(with_sister_ind)),
				   rep("Mother and sister",length(with_sister_ind))),
				level=c("Mother and sister","Sister","Mother","Full")))
	
ggplot(dataTemp,aes(y=Model,x=R2))+
geom_boxplot(show.legend=FALSE,width=0.4,col=colors()[180],fill="grey",outlier.shape=NA,size=0.2)+
geom_vline(xintercept=median(R2_full[with_sister_ind]),lty=2,lwd=0.2)+
scale_x_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.5))+
labs(x=expression(italic(R)^2),y="",title="")


ggsave("temporal_ODE_R2_mother_sister.pdf",width=unit(1.7,"in"),height=unit(1.6,"in"))




### full model

A<-t(ccl_model_full_coef[-1,])

plot(eigen(A)$values)

dataTemp<-data.frame(Re=Re(eigen(A)$values),Im=Im(eigen(A)$values))

ggplot(dataTemp,aes(x=Re,y=Im))+geom_point(col="#F8766D",size=0.5,stroke=0.2)+
geom_vline(xintercept=0,lwd=0.2)+geom_hline(yintercept=0,lwd=0.2)+geom_vline(xintercept=max(dataTemp$Re),lty=2,lwd=0.2)+
scale_x_continuous(limits=c(-2.7,2.7),breaks=c(-2,0,2))+
scale_y_continuous(limits=c(-0.38,0.38),breaks=c(-0.3,0,0.3))+
labs(x="Real",y="Imaginary",title=expression(paste("Eigenvalues of ",italic(K)[CCL])))

ggsave("temporal_ODE_eigen_full.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


### nested model without mother

A<-t(ccl_model_nested_without_mother_coef[-1,])

plot(eigen(A)$values)

dataTemp<-data.frame(Re=Re(eigen(A)$values),Im=Im(eigen(A)$values))

ggplot(dataTemp,aes(x=Re,y=Im))+geom_point(col="#F8766D",size=0.5,stroke=0.2)+
geom_vline(xintercept=0,lwd=0.2)+geom_hline(yintercept=0,lwd=0.2)+geom_vline(xintercept=max(dataTemp$Re),lty=2,lwd=0.2)+
scale_x_continuous(limits=c(-2.7,2.7),breaks=c(-2,0,2))+
scale_y_continuous(limits=c(-0.38,0.38),breaks=c(-0.3,0,0.3))+
labs(x="Real",y="Imaginary",title="No mother")

ggsave("temporal_ODE_eigen_nested_without_mother.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


### nested model without sister

A<-t(ccl_model_nested_without_sister_coef[-1,])

plot(eigen(A)$values)

dataTemp<-data.frame(Re=Re(eigen(A)$values),Im=Im(eigen(A)$values))

ggplot(dataTemp,aes(x=Re,y=Im))+geom_point(col="#F8766D",size=0.5,stroke=0.2)+
geom_vline(xintercept=0,lwd=0.2)+geom_hline(yintercept=0,lwd=0.2)+geom_vline(xintercept=max(dataTemp$Re),lty=2,lwd=0.2)+
scale_x_continuous(limits=c(-2.7,2.7),breaks=c(-2,0,2))+
scale_y_continuous(limits=c(-0.38,0.38),breaks=c(-0.3,0,0.3))+
labs(x="Real",y="Imaginary",title="No sister")

ggsave("temporal_ODE_eigen_nested_without_sister.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))



### R2 comparison - without mother or sister

pred_full<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%ccl_model_full_coef
R2_full<-sapply(1:345,function(i){cor(delta_X[,i],pred_full[,i])^2})

pred_without_mother<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%ccl_model_nested_without_mother_coef
R2_without_mother<-sapply(1:345,function(i){cor(delta_X[,i],pred_without_mother[,i])^2})

pred_without_sister<-cbind(rep(1,dim(X_mo)[1]),X_mo)%*%ccl_model_nested_without_sister_coef
R2_without_sister<-sapply(1:345,function(i){cor(delta_X[,i],pred_without_sister[,i])^2})

with_sister_ind<-which(!is.na(match(mo_sis_pair$da2[match(unlist(strsplit(colnames(delta_X),"_ccl")),mo_sis_pair$da1)],unlist(strsplit(colnames(X_mo),"_ccl")))))

dataTemp<-data.frame(R2=c(R2_full[with_sister_ind],R2_without_mother[with_sister_ind],R2_without_sister[with_sister_ind]),
	Model=factor(c(rep("Full",length(with_sister_ind)),
				   rep("No mother",length(with_sister_ind)),
				   rep("No sister",length(with_sister_ind))),
				level=c("No sister","No mother","Full")))
	
ggplot(dataTemp,aes(y=Model,x=R2))+
geom_boxplot(show.legend=FALSE,width=0.4,col=colors()[180],fill="grey",outlier.shape=NA,size=0.2)+
geom_vline(xintercept=median(R2_full[with_sister_ind]),lty=2,lwd=0.2)+
scale_x_continuous(limits=c(0,1),breaks=c(0,0.5,1))+
theme(axis.text.y=element_text(angle=0,hjust=1,vjust=0.5))+
labs(x=expression(paste("Model performance (",italic(R)^2,")")),y="",title="")


ggsave("temporal_ODE_R2_without_mother_sister.pdf",width=unit(1.9,"in"),height=unit(1.6,"in"))


######################################################################################
						cell contact 
######################################################################################



### cell contact 

cell_contact<-read.csv(paste0(focal_projectFolderPath,"/data/annotation/effective_connect_list.csv"))

cell_contact_pair<-do.call(rbind,strsplit(cell_contact[,1],"->"))

sort(table(c(apply(cell_contact_pair,1,function(x){paste0(x[1],"_",x[2])}),apply(cell_contact_pair,1,function(x){paste0(x[2],"_",x[1])}))))


cell_contact_pair2<-do.call(c,lapply(unique(c(cell_contact_pair)),function(x){
paste0(x,"_",c(cell_contact_pair[which(cell_contact_pair[,2]==x),1],cell_contact_pair[which(cell_contact_pair[,1]==x),2]))}))



#ccl 

ccl_full_coef_df<-do.call(rbind,lapply(1:dim(ccl_model_full_coef)[2],function(i){
	da_name<-colnames(ccl_model_full_coef)[i]
	mo_name<-rownames(ccl_model_full_coef)[i+1]
	focalInd<-setdiff(2:dim(ccl_model_full_coef)[1],setdiff(union(which(duplicated(rownames(ccl_model_full_coef))),
		union(which(rownames(ccl_model_full_coef)==da_name),which(rownames(ccl_model_full_coef)==mo_name))),i+1))
		
	focalCell<-unlist(strsplit(rownames(ccl_model_full_coef)[focalInd],"_ccl"))
	dataTemp<-data.frame(cell_pair=paste0(unlist(strsplit(colnames(ccl_model_full_coef)[i],"_ccl")),"_",focalCell),
		coef=ccl_model_full_coef[focalInd,i])
	
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

temp=ccl_full_coef_df[ccl_full_coef_df$mo_sis=="other",]
c1=length(which(temp$coef!=0&temp$contact==0))/length(which(temp$contact==0))
c2=length(which(temp$coef!=0&temp$contact==1))/length(which(temp$contact==1))

c1_se=sqrt(c1*(1-c1)/length(which(temp$contact==0)))
c2_se=sqrt(c2*(1-c2)/length(which(temp$contact==1)))


dataTemp<-data.frame(proportion=c(c2,c1),
	proportion_se=c(c2_se,c1_se)*2,
	type=c("contact","no contact"))
	
ggplot(dataTemp,aes(x=type,y=proportion,ymin=proportion-proportion_se,ymax=proportion+proportion_se))+
geom_bar(position=position_dodge2(padding=0.3),stat="identity",width=0.65,fill=colors()[182])+
geom_errorbar(position=position_dodge2(padding=0.3),stat="identity",width=0.25,lwd=0.2,show.legend=FALSE)+
scale_y_continuous(breaks=c(0,0.3,0.6))+
labs(x="Cell type",y="Proportion of\nnon-zero coefficients",title="")

ggsave("temporal_full_contact_non_zero_proportion.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


# effect size of coefficients 

dataTemp<-ccl_full_coef_df[ccl_full_coef_df$mo_sis=="other"&ccl_full_coef_df$coef!=0,]

dataTemp$contact<-factor(ifelse(dataTemp$contact==1,"contact","no contact"),level=c("contact","no contact"))
	
ggplot(dataTemp,aes(x=abs(coef),lty=contact))+
geom_density(show.legend=FALSE)+
scale_x_log10(limits=c(1e-5,1),breaks=c(1e-4,1e-2,1))+
scale_y_continuous(breaks=c(0,0.4,0.8))+
labs(x=expression(paste("Absolute value of\nnon-zero coefficients in ",italic(K))),y="Density",title="")
	
ggsave("temporal_full_contact_non_zero_distribution.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))

	
#########################################################################################################3


### block descendant 


## R2 compare with full model 

dataTemp<-data.frame(ccl_model_full_R2,ccl_model_block_descendant_R2)

#R2 comparison  

ggplot(dataTemp,aes(x=ccl_model_full_R2,y=ccl_model_block_descendant_R2))+
geom_point(size=0.3,col=colors()[182])+
geom_abline(slope=1,intercept=0,lwd=0.3)+
labs(x="full model",y="descendant blocked",title=expression(italic(R)[CCL]^2))

ggsave("temporal_descendant_block_R2_compare.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))

# eigenvalue

A<-t(ccl_model_block_descendant_coef[-1,])

plot(eigen(A)$values)

dataTemp<-data.frame(Re=Re(eigen(A)$values),Im=Im(eigen(A)$values))

ggplot(dataTemp,aes(x=Re,y=Im))+geom_point(col="#F8766D",size=0.5,stroke=0.2)+
geom_vline(xintercept=0,lwd=0.2)+geom_hline(yintercept=0,lwd=0.2)+geom_vline(xintercept=max(dataTemp$Re),lty=2,lwd=0.2)+
scale_x_continuous(limits=c(-2.7,2.7),breaks=c(-2,0,2))+
scale_y_continuous(limits=c(-0.38,0.38),breaks=c(-0.3,0,0.3))+
labs(x="Real",y="Imaginary",title=expression(paste(italic(K)[CCL]," of descendant blocked")))

ggsave("temporal_ODE_eigen_descendant_blocked.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


########################################################################

### compare kContr, kRNAi, kEnv , susceptability 

#

feature_names<-unlist(strsplit(colnames(ccl_raw_nonNA_scaling_noise),"_ccl"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]

#

X_mo_CCL<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$mo,feature_names)]
X_da1_CCL<-ccl_raw_nonNA_scaling_noise[,match(mo_sis_pair$da1,feature_names)]
delta_X_CCL<-X_da1_CCL-X_mo_CCL

kTest_kWT_kRNAi_CCL<-kTest(X_mo_CCL,delta_X_CCL,refInd=c(1:105),perturbInd=c(106:1928))

kTest_kWT_kEnv_CCL<-kTest(X_mo_CCL,delta_X_CCL,refInd=c(1:105),perturbInd=c(1929:2039))


#

dataTemp<-data.frame(kWT=kTest_kWT_kEnv_CCL[,4],kEnv=kTest_kWT_kEnv_CCL[,5],kRNAi=kTest_kWT_kRNAi_CCL[,5],kEnv_sus=kTest_kWT_kEnv_CCL[,1],kRNAi_sus=kTest_kWT_kRNAi_CCL[,1])

ggplot(dataTemp,aes(x=kWT,y=kEnv))+
geom_point(size=0.2,colour=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,c(1,2)]),breaks=c(-1,-0.5,0))+
scale_y_continuous(limits=range(dataTemp[,c(1,2)]),breaks=c(-1,-0.5,0))+
labs(x=expression(italic(k)[control]),y=expression(italic(k)[Env]),title="CCL")

ggsave("temporal_kContr_kENV.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


ggplot(dataTemp,aes(x=kWT,y=kRNAi))+
geom_point(size=0.2,colour=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,c(1,3)]),breaks=c(-1,-0.5,0))+
scale_y_continuous(limits=range(dataTemp[,c(1,3)]),breaks=c(-1,-0.5,0))+
labs(x=expression(italic(k)[control]),y=expression(italic(k)[RNAi]),title="CCL")

ggsave("temporal_kContr_kRNAi.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


ggplot(dataTemp,aes(x=kEnv_sus,y=kRNAi_sus))+
geom_point(size=0.2,stroke=0.2,colour=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
geom_hline(yintercept=0,lty=2,lwd=0.3)+
geom_vline(xintercept=0,lty=2,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,c(4,5)]),breaks=c(-0.5,0,0.5))+
scale_y_continuous(limits=range(dataTemp[,c(4,5)]),breaks=c(-0.5,0,0.5))+
labs(x=expression(paste(italic(k)[Env],"-",italic(k)[control])),y=expression(paste(italic(k)[RNAi],"-",italic(k)[control])),title="CCL")

ggsave("temporal_susceptibility_kEnv_kRNAi.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


###########################################################################33
