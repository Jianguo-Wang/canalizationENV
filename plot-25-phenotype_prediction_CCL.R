




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
library(yardstick)


focal_projectFolderPath<-r"(F:\wangj\Desktop\temp\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/figs"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/plot_triangle_tree.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/mo_da_link_cut.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/functions/triangle_tree_cut.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/embryo_lethal.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/cvglm_residual.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/cvglm_residual_comb.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/cvglm_CCL.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/cvglm_CCL_shuffle.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/residual_feature_comb.RData"
	),projectFolderPath=focal_projectFolderPath)



###

embryo_lethal<-as.numeric(embryo_lethal)

### shuffle 

set.seed(1)
embryo_lethal_shuffle<-embryo_lethal[sample(length(embryo_lethal))]


###############################################################################


# ROC curve and AUC by logistic regression under machine learning framework 

auc_summary<-data.frame(Subset=c("xyz","CCL","xyz+CCL","CCL(shuffled)"),
	AUC_cvm=c(cvglm_residual$cvm[which(cvglm_residual$lambda==cvglm_residual$lambda.min)],
		cvglm_CCL$cvm[which(cvglm_CCL$lambda==cvglm_CCL$lambda.min)],
		cvglm_residual_comb$cvm[which(cvglm_residual_comb$lambda==cvglm_residual_comb$lambda.min)],
		cvglm_CCL_shuffle$cvm[which(cvglm_CCL_shuffle$lambda==cvglm_CCL_shuffle$lambda.min)]),
	AUC_cvsd=c(cvglm_residual$cvsd[which(cvglm_residual$lambda==cvglm_residual$lambda.min)],
		cvglm_CCL$cvsd[which(cvglm_CCL$lambda==cvglm_CCL$lambda.min)],
		cvglm_residual_comb$cvsd[which(cvglm_residual_comb$lambda==cvglm_residual_comb$lambda.min)],
		cvglm_CCL_shuffle$cvsd[which(cvglm_CCL_shuffle$lambda==cvglm_CCL_shuffle$lambda.min)])
)


roc_list<-list(
	roc.glmnet(cvglm_residual$fit.preval, newy = embryo_lethal,family="binomial")[[which(cvglm_residual$lambda==cvglm_residual$lambda.min)]],
	roc.glmnet(cvglm_CCL$fit.preval, newy = embryo_lethal,family="binomial")[[which(cvglm_CCL$lambda==cvglm_CCL$lambda.min)]],
	roc.glmnet(cvglm_residual_comb$fit.preval, newy = embryo_lethal,family="binomial")[[which(cvglm_residual_comb$lambda==cvglm_residual_comb$lambda.min)]],
	roc.glmnet(cvglm_CCL_shuffle$fit.preval, newy = embryo_lethal_shuffle,family="binomial")[[which(cvglm_CCL_shuffle$lambda==cvglm_CCL_shuffle$lambda.min)]]
)

names(roc_list)<-c("xyz","CCL","xyz+CCL","CCL(shuffled)")


dataTemp<-data.frame(do.call(rbind,roc_list),
	Model=factor(c(rep("xyz",dim(roc_list[[1]])[1]),rep("CCL",dim(roc_list[[2]])[1]),rep("xyz+CCL",dim(roc_list[[3]])[1]),rep("CCL(shuffled)",dim(roc_list[[4]])[1])),
		level=c("xyz+CCL","xyz","CCL","CCL(shuffled)")))

ggplot(dataTemp,aes(x=FPR,y=TPR,col=Model))+
geom_line()+
geom_abline(slope=1,intercept=0,lty=2)+
scale_x_continuous(breaks=c(0,0.5,1))+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(title="Prediction of lethality")+
theme(legend.position=c(0.7,0.3))


ggsave("Prediction_embryonic_lethal_CCL.pdf",width=unit(2,"in"),height=unit(2.1,"in"))


paste("(AUC=",paste(round(auc_summary$AUC_cvm[order(auc_summary$AUC_cvm,decreasing=T)],2),round(auc_summary$AUC_cvsd[order(auc_summary$AUC_cvm,decreasing=T)],2),sep="±"),")",sep="")


# lineage with effectCells highlighted
effectCells<-colnames(q1_cell350_alignment_xyz_noise)[which(coef(cvglm_residual)!=0)[-1]]
effectCells_CCL<-colnames(ccl_raw_nonNA_scaling_noise)[which(coef(cvglm_CCL)!=0)[-1]]

effectCells_name<-do.call(rbind,strsplit(effectCells,"_"))[,1]
effectCells_CCL_name<-do.call(rbind,strsplit(effectCells_CCL,"_"))[,1]

plot_triangle_tree(triangle_tree_cut,focal_cells=c(effectCells_name,effectCells_CCL_name),
	focal_cells_col=c(rep("#F8766D",length(effectCells)),rep("#00BFC4",length(effectCells_CCL))),
	focal_cells_size=c(rep(0.2,length(effectCells)),rep(0.5,length(effectCells_CCL))),
	back_col="grey")


ggsave("effectCells_lineage_distribution_xyz_and_CCL.pdf",width=unit(4.2,"in"),height=unit(2.3,"in"))








