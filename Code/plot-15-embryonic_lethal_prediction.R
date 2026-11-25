



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
library(yardstick)


focal_projectFolderPath<-r"(E:\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/figs"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/plot_triangle_tree.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/mo_da_link_cut.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/functions/triangle_tree_cut.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info_short.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/embryo_lethal.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/spatial/xyz_x_model_mother_coef.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/spatial/xyz_y_model_mother_coef.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/spatial/xyz_z_model_mother_coef.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/cvglm_residual_training.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/cvglm_residual.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/cvglm_residual_shuffle.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/cvglm_residual_g9.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/cvglm_residual_g8.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/cvglm_residual_g7.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/cvglm_residual_g6.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/cvglm_residual_g5.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/cvglm_residual_g4.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/prediction/residual_feature.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/depth/depth_summary.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData"
	),projectFolderPath=focal_projectFolderPath)




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

###

embryo_lethal<-as.numeric(embryo_lethal)

###

embryoGene<-c(do.call(rbind,strsplit(rownames(q1_cell350_alignment_xyz_noise)[1:1928],"_"))[,1],do.call(rbind,strsplit(rownames(q1_cell350_alignment_xyz_noise)[1929:2039],"_"))[,3])

trainIndex<-train_partition(embryoGene)

residual_feature_train<-residual_feature[trainIndex,]
embryo_lethal_train<-embryo_lethal[trainIndex]

residual_feature_test<-residual_feature[setdiff(1:dim(residual_feature)[1],trainIndex),]
embryo_lethal_test<-embryo_lethal[setdiff(1:dim(residual_feature)[1],trainIndex)]

coo_x<-apply(q1_cell350_alignment_xyz[1:105,1:714],2,mean)
coo_y<-apply(q1_cell350_alignment_xyz[1:105,c(1:714)+714],2,mean)
coo_z<-apply(q1_cell350_alignment_xyz[1:105,c(1:714)+714*2],2,mean)


##########################################################################################
							plot : depth
##########################################################################################

# definition

logistic_func<-function(x){exp(x)/(1+exp(x))}

ggplot()+geom_function(fun=logistic_func,xlim=c(-5,5))+
geom_vline(xintercept=0,lty=2,lwd=0.3)+
geom_hline(yintercept=c(0,0.5,1),lty=2,lwd=0.3)+
scale_x_continuous(breaks=c(-5,0,5))+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(x="Noise level",y="Probability of lethality",title=expression(paste("Noise tolerance (",italic(d),")")))

ggsave("spatial_depth_definition.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



# distribution of depth 


dataTemp<-data.frame(depth_summary,coo=c(rep("x",714),rep("y",714),rep("z",714)))

ggplot(dataTemp,aes(x=depth_abs,col=coo))+
#geom_histogram(colour="grey",fill="grey39",size=0.2)+
geom_density(lwd=0.5,fill=NA)+
#scale_x_continuous(breaks=round(c(range(depth_summary$depth_abs)[1],median(depth_summary$depth_abs),range(depth_summary$depth_abs)[2]),0))+
theme(legend.position=c(0.8,0.8))+
labs(x=expression(italic(d)),y="Count",title="")

ggsave("spatial_depth_distribution.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


# compare between depth_abs and slope k 

dataTemp<-data.frame(depth=depth_summary$depth_abs[c(match(mo_sis_pair$da1,feature_names),match(mo_sis_pair$da1,feature_names)+714,match(mo_sis_pair$da1,feature_names)+714*2)],
	k=c(xyz_x_model_mother_coef,xyz_y_model_mother_coef,xyz_z_model_mother_coef),
	coo=c(rep("x",706),rep("y",706),rep("z",706)))


ggplot(dataTemp,aes(x=k,y=depth,col=coo))+
geom_point(size=0.15,stroke=0.1,show.legend=FALSE)+
geom_smooth(method="lm",lwd=0.3,show.legend=FALSE)+
scale_x_continuous(breaks=c(-1,-0.5,0))+
scale_y_continuous(breaks=c(25,50,75))+
labs(x=expression(italic(k)),y=expression(italic(d)),title="")

ggsave("spatial_depth_Vs_k.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



cor.test(dataTemp[dataTemp$coo=="x",1],dataTemp[dataTemp$coo=="x",2])
cor.test(dataTemp[dataTemp$coo=="y",1],dataTemp[dataTemp$coo=="y",2])
cor.test(dataTemp[dataTemp$coo=="z",1],dataTemp[dataTemp$coo=="z",2])

#
noise_level<-apply(q1_cell350_alignment_xyz,2,sd)

dataTemp<-data.frame(depth=depth_summary$depth_abs[c(match(mo_sis_pair$da1,feature_names),match(mo_sis_pair$da1,feature_names)+714,match(mo_sis_pair$da1,feature_names)+714*2)],
	k=c(xyz_x_model_mother_coef,xyz_y_model_mother_coef,xyz_z_model_mother_coef),
	noise_level=noise_level[c(match(mo_sis_pair$da1,feature_names),match(mo_sis_pair$da1,feature_names)+714,match(mo_sis_pair$da1,feature_names)+714*2)],
	coo=c(rep("x",706),rep("y",706),rep("z",706)))

plot(dataTemp$k,dataTemp$noise_level)
plot(dataTemp$depth,dataTemp$noise_level)

cor.test(residuals(lm(dataTemp$depth~dataTemp$noise_level+dataTemp$coo)),residuals(lm(dataTemp$k~dataTemp$noise_level+dataTemp$coo)))

dataTemp2<-data.frame(depth_res=residuals(lm(dataTemp$depth~dataTemp$noise_level+dataTemp$coo)),
	k_res=residuals(lm(dataTemp$k~dataTemp$noise_level+dataTemp$coo)),
	coo=dataTemp$coo)

ggplot(dataTemp2,aes(x=k_res,y=depth_res,col=coo))+
geom_point(size=0.2,stroke=0.2,show.legend=FALSE)+
geom_smooth(method="lm",lwd=0.3,show.legend=FALSE)



# depth and pvalue 

ggplot(depth_summary,aes(x=depth_abs,y=-log10(p.adjust(p_abs,method="BH"))))+
geom_point(size=0.2,stroke=0.2,colour="grey39")+
geom_hline(yintercept=-log10(0.05),lwd=0.3,lty=2)+
scale_y_continuous(breaks=c(0,10,20))+
labs(x=expression(italic(d)),y=expression(paste("-",log[10],"(",italic(p),")")),title="")

ggsave("spatial_depth_pvalues.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


# compare between depth_pos and depth_neg 

col1<-rep("grey39",dim(depth_summary)[1])
col1[which(p.adjust(depth_summary$p_diff,method="BH")<0.05)]="#F8766D"

size1<-rep(0.2,dim(depth_summary)[1])
size1[which(p.adjust(depth_summary$p_diff,method="BH")<0.05)]=1

alpha1<-rep(0.75,dim(depth_summary)[1])
alpha1[which(p.adjust(depth_summary$p_diff,method="BH")<0.05)]=1


ggplot(depth_summary,aes(x=depth_pos,y=depth_neg))+
geom_point(size=size1,stroke=0.2,col=col1,alpha=alpha1)+
geom_abline(slope=1,intercept=0,lwd=0.3)+
scale_x_continuous(limits=c(-0.14,0),breaks=c(-0.14,-0.07,0))+
scale_y_continuous(limits=c(-0.14,0),breaks=c(-0.14,-0.07,0))+
labs(x="Depth of canals (+)",y="Depth of canals (-)",title="")


# compare between depth_abs and depth_pos 

ggplot(depth_summary,aes(x=depth_abs,y=depth_pos))+
geom_point(size=0.2,stroke=0.2,col="grey39")+
geom_abline(slope=1,intercept=0,lwd=0.3)+
scale_x_continuous(limits=c(-0.14,0),breaks=c(-0.14,-0.07,0))+
scale_y_continuous(limits=c(-0.14,0),breaks=c(-0.14,-0.07,0))+
labs(x="Depth of canals",y="Depth of canals (+)",title="")


# compare between depth_abs and depth_neg  

ggplot(depth_summary,aes(x=depth_abs,y=depth_neg))+
geom_point(size=0.2,stroke=0.2,col="grey39")+
geom_abline(slope=1,intercept=0,lwd=0.3)+
scale_x_continuous(limits=c(-0.14,0),breaks=c(-0.14,-0.07,0))+
scale_y_continuous(limits=c(-0.14,0),breaks=c(-0.14,-0.07,0))+
labs(x="Depth of canals",y="Depth of canals (-)",title="")


### top depth enrichment 

cells_d_large<-feature_names[which(depth_summary[c(1:714),1]>quantile(depth_summary[c(1:714),1],0.95)|depth_summary[c(1:714)+714,1]>quantile(depth_summary[c(1:714)+714,1],0.95)|depth_summary[c(1:714)+714*2,1]>quantile(depth_summary[c(1:714)+714*2,1],0.95))]
cells_d_small<-feature_names[which(depth_summary[c(1:714),1]<quantile(depth_summary[c(1:714),1],0.05)|depth_summary[c(1:714)+714,1]<quantile(depth_summary[c(1:714)+714,1],0.05)|depth_summary[c(1:714)+714*2,1]<quantile(depth_summary[c(1:714)+714*2,1],0.05))]


TissueFC_group<-read.csv(paste0(focal_projectFolderPath,"/data/annotation/TissueFC_group.csv"))
TissueFC_group<-TissueFC_group[TissueFC_group[,2]!="Germline"&TissueFC_group[,2]!="Death",]


get_chisq_levelP<-function(cell_vec,Group){
	
	if(length(cell_vec)>0&length(na.omit(match(cell_vec,Group[,1])))>0){
		focal_index<-na.omit(match(cell_vec,Group[,1]))
		focal_list<-Group[focal_index,2]
		Expect<-table(Group[match(feature_names,Group[,1]),2])
		Observ<-sapply(names(Expect),function(x){length(which(focal_list==x))})
	
		temp<-sapply(names(Expect),function(x){
			Expect_focalLevel<-c(Expect[names(Expect)==x],Expect[names(Expect)!=x])
			Observ_focalLevel<-c(Observ[names(Expect)==x],Observ[names(Expect)!=x])
			return(binom.test(Observ_focalLevel[1],sum(Observ_focalLevel),p=Expect_focalLevel[1]/sum(Expect_focalLevel),alternative="greater")$p.value)

		})
	
		names(temp)<-names(Expect)		
	}else{
		temp<-rep(NA,length(table(Group[,2])))
	}
	
	
	return(temp)

}



table(TissueFC_group[match(cells_d_large,TissueFC_group[,1]),2])
table(TissueFC_group[match(cells_d_small,TissueFC_group[,1]),2])

dataTemp<-data.frame(cell=rep(feature_names,3),
	type=rep(TissueFC_group[match(feature_names,TissueFC_group[,1]),2],3),
	d=depth_summary[,1],
	k=c(xyz_x_model_mother_coef[match(feature_names,mo_sis_pair$da1)],xyz_y_model_mother_coef[match(feature_names,mo_sis_pair$da1)],xyz_z_model_mother_coef[match(feature_names,mo_sis_pair$da1)]),
	coo=c(rep("x",714),rep("y",714),rep("z",714)),
	coo_x=rep(coo_x,3),
	coo_y=rep(coo_y,3),
	coo_z=rep(coo_z,3),
	generation=sapply(match(feature_names,lineage_info),function(x){rowColRec(x,dim(lineage_info)[1],dim(lineage_info)[2])[1]})
	)

dataTemp2<-data.frame(tissue=names(get_chisq_levelP(cells_d_large,TissueFC_group)),
	p_large=get_chisq_levelP(cells_d_large,TissueFC_group),
	p_small=get_chisq_levelP(cells_d_small,TissueFC_group))

ggplot(dataTemp,aes(x=type,y=d,col=coo))+
geom_boxplot()
geom_jitter(height=0,width=0.15)

ggplot(dataTemp,aes(x=generation,y=d,col=coo))+
#geom_jitter(height=0,alpha=0.25)+
geom_smooth()+
labs(x="Generation",y=expression(italic(d)),title="")

ggsave("spatial_d_generation_distribution.pdf",width=unit(3,"in"),height=unit(1.7,"in"))


ggplot(dataTemp,aes(x=generation,y=k,col=coo))+
#geom_jitter(height=0,alpha=0.25)+
geom_smooth()+
labs(x="Generation",y=expression(italic(k)),title="")

ggsave("spatial_k_generation_distribution.pdf",width=unit(3,"in"),height=unit(1.7,"in"))


dataTemp4<-data.frame(coo_x,coo_y,coo_z,d_x=depth_summary[1:714,1],d_y=depth_summary[c(1:714)+714,1],d_x=depth_summary[c(1:714)+714*2,1])

dataTemp5<-data.frame(d=rep(depth_summary[,1],3),
	k=rep(c(xyz_x_model_mother_coef[match(feature_names,mo_sis_pair$da1)],xyz_y_model_mother_coef[match(feature_names,mo_sis_pair$da1)],xyz_z_model_mother_coef[match(feature_names,mo_sis_pair$da1)]),3),
	coo=c(rep(coo_x,3),rep(coo_y,3),rep(coo_z,3)),
	d_type=rep(c(rep("d_x",714),rep("d_y",714),rep("d_z",714)),3),
	k_type=rep(c(rep("k_x",714),rep("k_y",714),rep("k_z",714)),3),
	coo_type=c(rep("coo_x",714*3),rep("coo_y",714*3),rep("coo_z",714*3)))

ggplot(dataTemp5,aes(x=coo,y=d))+
geom_point(size=0.1,col="grey",alpha=0.25)+
geom_smooth(lwd=0.3)+
facet_grid(vars(d_type),vars(coo_type),scales="free")+
labs(x="Coordinate",y=expression(italic(d)),title="")

ggsave("spatial_d_spatial_distribution.pdf",width=unit(3,"in"),height=unit(4,"in"))


ggplot(dataTemp5,aes(x=coo,y=k))+
geom_point(size=0.1,col="grey",alpha=0.25)+
geom_smooth(lwd=0.3)+
facet_grid(vars(k_type),vars(coo_type),scales="free")+
labs(x="Coordinate",y=expression(italic(k)),title="")

ggsave("spatial_k_spatial_distribution.pdf",width=unit(3,"in"),height=unit(4,"in"))


focal_cells<-lineage_info[4:10,match(c("ABalaaaala","ABaraaaaaa","ABplaaaaaa","ABpraaaaaa","MSaaaaaa","Ealaa","Caaaaa","Daaa"),colnames(lineage_info))]

dataTemp[na.omit(match(focal_cells[,8],feature_names)),]

dataTemp3<-dataTemp[,]


ggplot(dataTemp2,aes(x=-log10(p_large),y=-log10(p_small)))+
geom_point(col=colors()[182])+
geom_vline(xintercept=1.3,lty=2,lwd=0.3,col="grey")+
geom_hline(yintercept=1.3,lty=2,lwd=0.3,col="grey")+
labs(x=expression(paste("-",log[10],"(",italic(p),") for large-",italic(d)," cells")),
	y=expression(paste("-",log[10],"(",italic(p),") for small-",italic(d)," cells")),title="Tissue enrichment")

ggsave("spatial_depth_tissue_enrich.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### correlation between three coordinates
temp<-data.frame(depth_summary[1:714,1],depth_summary[c(1:714)+714,1],depth_summary[c(1:714)+714*2,1])
cor(temp)
par(mfrow=c(2,2))
plot(temp[,1],temp[,2])
plot(temp[,1],temp[,3])
plot(temp[,2],temp[,3])


############################################################################### 
						plot : prediction 
###############################################################################

###

effectCells<-colnames(q1_cell350_alignment_xyz_noise)[which(coef(cvglm_residual)!=0)[-1]]

unique(do.call(rbind,strsplit(effectCells,"_"))[,1])

coo_x<-apply(q1_cell350_alignment_xyz[1:105,1:714],2,mean)
coo_y<-apply(q1_cell350_alignment_xyz[1:105,c(1:714)+714],2,mean)
coo_z<-apply(q1_cell350_alignment_xyz[1:105,c(1:714)+714*2],2,mean)

coo_wt<-apply(q1_cell350_alignment_xyz[1:105,],2,mean)

# noise level comparison 

sd_hatch<-sapply(1:dim(q1_cell350_alignment_xyz[embryo_lethal==0,])[2],function(i){sqrt(sum((q1_cell350_alignment_xyz[embryo_lethal==0,i]-coo_wt[i])^2)/length(q1_cell350_alignment_xyz[embryo_lethal==0,i]))})
sd_lethal<-sapply(1:dim(q1_cell350_alignment_xyz[embryo_lethal==1,])[2],function(i){sqrt(sum((q1_cell350_alignment_xyz[embryo_lethal==1,i]-coo_wt[i])^2)/length(q1_cell350_alignment_xyz[embryo_lethal==1,i]))})

dataTemp<-data.frame(sd_hatch,sd_lethal)

ggplot(dataTemp,aes(x=sd_hatch,y=sd_lethal))+
geom_point(size=0.1,stroke=0.1,col="grey39")+
geom_abline(slope=1,intercept=0,lwd=0.3,lty=2)+
scale_x_continuous(limits=range(dataTemp),breaks=c(10,25,40))+
scale_y_continuous(limits=range(dataTemp),breaks=c(10,25,40))+
labs(x="Hatched embryos",y="Lethal embryos",title="SD of cell noise")

ggsave("Prediction_compare_SD.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


# proportion of lethal embryos 

embryo_lethal_cont<-data.frame(Total=c(105,1823,111),
	Lethal=c(length(which(embryo_lethal[1:105]==1)),length(which(embryo_lethal[106:1928]==1)),length(which(embryo_lethal[1929:2039]==1))))
embryo_lethal_cont$Alive<-embryo_lethal_cont$Total-embryo_lethal_cont$Lethal

dataTemp<-data.frame(Num=c(embryo_lethal_cont$Lethal,embryo_lethal_cont$Alive),state=c(rep("lethal",3),rep("hatched",3)),type=rep(c("WT","RNAi","Env"),2))

ggplot(dataTemp,aes(x=type,y=Num,fill=state))+
geom_bar(position="fill",stat="identity")+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(x="Type of embryos",y="Proportion of embryos",title="")

ggsave("Prediction_proportion_of_lethal_embryos.pdf",width=unit(2.5,"in"),height=unit(2.1,"in"))


# ROC curve and AUC by logistic regression under machine learning framework 

auc_summary<-data.frame(Subset=c("Total","Training","Testing","Shuffle"),
	AUC_cvm=c(cvglm_residual$cvm[cvglm_residual$lambda==cvglm_residual$lambda.min],
		cvglm_residual_training$cvm[cvglm_residual_training$lambda==cvglm_residual_training$lambda.min],
		assess.glmnet(predict(cvglm_residual_training,residual_feature_test,s="lambda.min"), newy = embryo_lethal_test,family="binomial")$auc,
		cvglm_residual_shuffle$cvm[cvglm_residual_shuffle$lambda==cvglm_residual_shuffle$lambda.min]),
	AUC_cvsd=c(cvglm_residual$cvsd[cvglm_residual$lambda==cvglm_residual$lambda.min],
		cvglm_residual_training$cvsd[cvglm_residual_training$lambda==cvglm_residual_training$lambda.min],
		NA,
		cvglm_residual_shuffle$cvsd[cvglm_residual_shuffle$lambda==cvglm_residual_shuffle$lambda.min]))

cal_AUC<-function(df){
	auc_score <- sum(diff(df$FPR) * (df$TPR[-1] + df$TPR[-length(df$TPR)]) / 2)
	return(auc_score)
}

#total
total_roc<-lapply(1:10,function(i){roc.glmnet(cvglm_residual$fit.preval[cvglm_residual$foldid!=i,cvglm_residual$lambda==cvglm_residual$lambda.min],newy=embryo_lethal[cvglm_residual$foldid!=i])})
total_roc_refine<-lapply(total_roc,function(df){temp<-approx(df,xout=total_roc[[1]][,1]);return(data.frame(FPR=temp$x,TPR=temp$y))})
total_roc_average<-data.frame(FPR=total_roc_refine[[1]][,1],TPR=rowMeans(sapply(total_roc_refine,function(df){df[,2]})))
cal_AUC(total_roc_average)

#train
train_roc<-lapply(1:10,function(i){roc.glmnet(cvglm_residual_training$fit.preval[cvglm_residual_training$foldid!=i,cvglm_residual_training$lambda==cvglm_residual_training$lambda.min],newy=embryo_lethal_train[cvglm_residual_training$foldid!=i])})
train_roc_refine<-lapply(train_roc,function(df){temp<-approx(df,xout=train_roc[[1]][,1]);return(data.frame(FPR=temp$x,TPR=temp$y))})
train_roc_average<-data.frame(FPR=train_roc_refine[[1]][,1],TPR=rowMeans(sapply(train_roc_refine,function(df){df[,2]})))
cal_AUC(train_roc_average)

#test 
train_roc<-roc.glmnet(predict(cvglm_residual_training,residual_feature_test,s="lambda.min"), newy = embryo_lethal_test)

#shuffle
set.seed(1)
embryo_lethal_shuffle<-embryo_lethal[sample(length(embryo_lethal))]

shuffle_roc<-lapply(1:10,function(i){roc.glmnet(cvglm_residual_shuffle$fit.preval[cvglm_residual_shuffle$foldid!=i,cvglm_residual_shuffle$lambda==cvglm_residual_shuffle$lambda.min],newy=embryo_lethal_shuffle[cvglm_residual_shuffle$foldid!=i])})
shuffle_roc_refine<-lapply(shuffle_roc,function(df){temp<-approx(df,xout=shuffle_roc[[1]][,1]);return(data.frame(FPR=temp$x,TPR=temp$y))})
shuffle_roc_average<-data.frame(FPR=shuffle_roc_refine[[1]][,1],TPR=rowMeans(sapply(shuffle_roc_refine,function(df){df[,2]}),na.rm=T))
cal_AUC(shuffle_roc_average)


roc_list<-list(
roc.glmnet(cvglm_residual$fit.preval, newy = embryo_lethal)[[which(cvglm_residual$lambda==cvglm_residual$lambda.min)]],
roc.glmnet(cvglm_residual_training$fit.preval, newy = embryo_lethal_train)[[which(cvglm_residual_training$lambda==cvglm_residual_training$lambda.min)]],
roc.glmnet(predict(cvglm_residual_training,residual_feature_test,s="lambda.min"), newy = embryo_lethal_test),
roc.glmnet(cvglm_residual_shuffle$fit.preval, newy = embryo_lethal_shuffle)[[which(cvglm_residual_shuffle$lambda==cvglm_residual_shuffle$lambda.min)]]
)

names(roc_list)<-c("Total","Training","Testing","Shuffle")


dataTemp<-data.frame(do.call(rbind,roc_list),
	Subset=factor(c(rep("Total",dim(roc_list[[1]])[1]),rep("Training",dim(roc_list[[2]])[1]),rep("Testing",dim(roc_list[[3]])[1]),rep("Shuffle",dim(roc_list[[4]])[1])),
		level=c("Total","Training","Testing","Shuffle")))

ggplot(dataTemp,aes(x=FPR,y=TPR,col=Subset))+
geom_line(lwd=0.5)+
geom_abline(slope=1,intercept=0,lty=2)+
theme(legend.position=c(0.5,0.2))+
scale_x_continuous(breaks=c(0,0.5,1))+
scale_y_continuous(breaks=c(0,0.5,1))+
labs(title="Prediction of lethality")

ggsave("Prediction_of_embryonic_lethality_ROC.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))

paste("(AUC=",paste(round(auc_summary$AUC_cvm[order(auc_summary$AUC_cvm,decreasing=T)],2),round(auc_summary$AUC_cvsd[order(auc_summary$AUC_cvm,decreasing=T)],2),sep="±"),")",sep="")


# lineage with effectCells highlighted

plot_triangle_tree(triangle_tree_cut,focal_cells=do.call(rbind,strsplit(effectCells,"_"))[,1],focal_cells_col="#F8766D",focal_cells_size=0.2,back_col="grey")

ggsave("effectCells_lineage_distribution.pdf",width=unit(4,"in"),height=unit(2.3,"in"))


# AUC ~ generation 

dataTemp<-data.frame(generation=factor(c(10:4),level=c(4:10)),
AUC=c(cvglm_residual$cvm[cvglm_residual$lambda==cvglm_residual$lambda.min],
cvglm_residual_g9$cvm[cvglm_residual_g9$lambda==cvglm_residual_g9$lambda.min],
cvglm_residual_g8$cvm[cvglm_residual_g8$lambda==cvglm_residual_g8$lambda.min],
cvglm_residual_g7$cvm[cvglm_residual_g7$lambda==cvglm_residual_g7$lambda.min],
cvglm_residual_g6$cvm[cvglm_residual_g6$lambda==cvglm_residual_g6$lambda.min],
cvglm_residual_g5$cvm[cvglm_residual_g5$lambda==cvglm_residual_g5$lambda.min],
cvglm_residual_g4$cvm[cvglm_residual_g4$lambda==cvglm_residual_g4$lambda.min]),
AUC_SE=c(cvglm_residual$cvsd[cvglm_residual$lambda==cvglm_residual$lambda.min],
cvglm_residual_g9$cvsd[cvglm_residual_g9$lambda==cvglm_residual_g9$lambda.min],
cvglm_residual_g8$cvsd[cvglm_residual_g8$lambda==cvglm_residual_g8$lambda.min],
cvglm_residual_g7$cvsd[cvglm_residual_g7$lambda==cvglm_residual_g7$lambda.min],
cvglm_residual_g6$cvsd[cvglm_residual_g6$lambda==cvglm_residual_g6$lambda.min],
cvglm_residual_g5$cvsd[cvglm_residual_g5$lambda==cvglm_residual_g5$lambda.min],
cvglm_residual_g4$cvsd[cvglm_residual_g4$lambda==cvglm_residual_g4$lambda.min]))

ggplot(dataTemp,aes(x=generation,y=AUC,ymin=AUC-AUC_SE,ymax=AUC+AUC_SE))+
geom_point(size=0.5)+
geom_line(aes(x=generation,y=AUC,group=NA))+
geom_hline(yintercept=1.0,lty=2)+
geom_errorbar(width=0.2,lwd=0.3)+
scale_y_continuous(limits=c(0.5,1.0),breaks=c(0.5,1.0))+
labs(y="AUC",x="Generation",title="")

ggsave("Prediction_of_embryonic_lethality_ROC_generation.pdf",width=unit(1.6,"in"),height=unit(1.3,"in"))


###






# case showing alive and lethal fates of a same RNAi 

dataTemp<-data.frame(embryo_lethal,predicted=ifelse(predict(cvglm_residual,residual_feature,s="lambda.min",family="binomial")>0,1,0))

do.call(rbind,lapply(unique(embryoGene[which(dataTemp[,1]==dataTemp[,2]&dataTemp[,1]==1)]),function(x){dataTemp[embryoGene==x,]}))

needCheckEmbID<-rownames(q1_cell350_alignment_xyz)[which(dataTemp[,1]==dataTemp[,2]&dataTemp[,1]==1)]

needCheckEmbID_cell<-lapply(needCheckEmbID,function(x){
	na.omit(sapply(1:length(effectCells),function(i){
		temp<-q1_cell350_alignment_xyz_noise[rownames(q1_cell350_alignment_xyz)==x,effectCells[i]]
		temp1<-q1_cell350_alignment_xyz_noise[,effectCells[i]]
		
		temp2<-ifelse(temp>quantile(temp1,0.95)|temp<quantile(temp1,0.05),effectCells[i],NA)
		
		if(temp>quantile(temp1,0.95)){names(temp2)<-"0.95"}
		if(temp<quantile(temp1,0.05)){names(temp2)<-"0.05"}
		
		
		return(temp2)
			
	}))

})

names(needCheckEmbID_cell)<-needCheckEmbID


needCheckEmbID_cell_df<-data.frame(embryo=unlist(lapply(1:length(needCheckEmbID),function(i){rep(needCheckEmbID[i],length(needCheckEmbID_cell[[i]]))})),cell=do.call(c,needCheckEmbID_cell),outlier=unlist(lapply(needCheckEmbID_cell,names)))
write.csv(needCheckEmbID_cell_df,paste0(focal_projectFolderPath,"/data/prediction/needCheckEmbID_cell_df.csv"),row.names=F)


#

focalGeneVec<-c("CSNK-1","DHS-4","DXBP-1","F43G9.12","FKB-5","PRP-4","RPL-31","SEC-12","SKR-2","SMN-1")


focalGene<-focalGeneVec[6]



which(embryoGene==focalGene)

dataTemp<-data.frame(effectCells=factor(rep(effectCells,length(which(embryoGene==focalGene))),level=effectCells),
	Noise=c(t(q1_cell350_alignment_xyz_noise[which(embryoGene==focalGene),which(coef(cvglm_residual)!=0)[-1]])),
	embryoID=c(sapply(paste(paste("PRP-4 embryo",c(1:length(which(embryoGene==focalGene))),sep="_"),ifelse(embryo_lethal[which(embryoGene==focalGene)]==1,"lethal","hatched"),sep=" : "),rep,length(effectCells))))

embryo_lethal[which(embryoGene==focalGene)]

percent001<-data.frame(t(apply(q1_cell350_alignment_xyz_noise[,effectCells],2,quantile,c(0.005,0.995))))
names(percent001)<-c("p0005","p0095")

percent01<-data.frame(t(apply(q1_cell350_alignment_xyz_noise[,effectCells],2,quantile,c(0.05,0.95))))
names(percent01)<-c("p005","p095")

pPoly<-data.frame(id=c(sapply(1:(length(effectCells)-1),rep,4)),
	value=rep("grey",(length(effectCells)-1)*4),
	y1=c(sapply(1:(length(effectCells)-1),function(i){c(percent001[i:(i+1),1],rev(percent001[i:(i+1),2]))})),
	x1=c(sapply(1:(length(effectCells)-1),function(i){c(i:(i+1),(i+1):i)})))

pPoly2<-data.frame(id=c(sapply(1:(length(effectCells)-1),rep,4)),
	value=rep("grey",(length(effectCells)-1)*4),
	y1=c(sapply(1:(length(effectCells)-1),function(i){c(percent01[i:(i+1),1],rev(percent01[i:(i+1),2]))})),
	x1=c(sapply(1:(length(effectCells)-1),function(i){c(i:(i+1),(i+1):i)})))

#c(sapply(1:length(which(embryoGene==focalGene)),function(i){effectCellsOrder+length(effectCells)*(i-1)}))

ggplot(dataTemp,aes(x=effectCells,y=Noise,col=embryoID,group=embryoID))+
geom_point(size=0.5)+
geom_polygon(data=pPoly,aes(x=x1,y=y1,group=id),color="grey",fill="lightGrey",alpha=0.15,lwd=0.2)+
#geom_polygon(data=pPoly2,aes(x=x1,y=y1,group=id),color="grey",fill="Grey1",alpha=0.1)+
scale_y_continuous(breaks=c(-40,0,40))+
labs(x="Effect cells",y="Noise level",title="")+
theme(axis.text.x=element_text(size=4,angle=90,hjust=1,vjust=0.25),legend.position=c(0.7,0.95))


ggsave("effectCells_PRP-4_noise_level.pdf",width=unit(4.5,"in"),height=unit(2.56,"in"))


# emb : 2 

focalEmbIndex=which(embryoGene==focalGeneVec[6])[2]
focalCells<-effectCells[which(abs(q1_cell350_alignment_xyz_noise[focalEmbIndex,which(coef(cvglm_residual)!=0)[-1]])>30)]
focalCells 

focalLineage<-lineage_info[4:10,sapply(do.call(rbind,strsplit(focalCells,"_"))[,1],function(x){rowColRec(which(lineage_info==x)[1],dim(lineage_info)[1],dim(lineage_info)[2])[2]})]
focal_lineage_x<-unlist(lapply(1:dim(focalLineage)[2],function(i){q1_cell350_alignment_xyz[focalEmbIndex,na.omit(match(paste0(focalLineage[,i],"_x"),colnames(q1_cell350_alignment_xyz)))]}))
contr_lineage_x<-unlist(lapply(1:dim(focalLineage)[2],function(i){coo_x[na.omit(match(paste0(focalLineage[,i],"_x"),colnames(q1_cell350_alignment_xyz)))]}))
focal_lineage_y<-unlist(lapply(1:dim(focalLineage)[2],function(i){q1_cell350_alignment_xyz[focalEmbIndex,na.omit(match(paste0(focalLineage[,i],"_y"),colnames(q1_cell350_alignment_xyz)))]}))
contr_lineage_y<-unlist(lapply(1:dim(focalLineage)[2],function(i){coo_y[na.omit(match(paste0(focalLineage[,i],"_x"),colnames(q1_cell350_alignment_xyz)))]}))

dataTemp<-data.frame(coo_x=c(coo_x,focal_lineage_x,contr_lineage_x),coo_y=c(coo_y,focal_lineage_y,contr_lineage_y),
lineage=c(rep("WT",length(coo_x)),unlist(lapply(1:length(focalCells),function(i){rep(focalCells[i],length(na.omit(match(paste0(focalLineage[,i],"_x"),colnames(q1_cell350_alignment_xyz)))))})),
	unlist(lapply(1:length(focalCells),function(i){rep(paste0(focalCells[i],"(WT)"),length(na.omit(match(paste0(focalLineage[,i],"_x"),colnames(q1_cell350_alignment_xyz)))))}))))

ggplot(data=dataTemp[dataTemp$lineage=="WT",],aes(x=coo_x,y=coo_y))+
geom_point(size=0.2,col="lightGrey")+
geom_point(data=dataTemp[dataTemp$lineage==focalCells[1]|dataTemp$lineage==paste0(focalCells[1],"(WT)"),],aes(x=coo_x,y=coo_y,col=lineage),size=0.75,alpha=0.7,show.legend=F)+
geom_path(data=dataTemp[dataTemp$lineage==focalCells[1]|dataTemp$lineage==paste0(focalCells[1],"(WT)"),],aes(x=coo_x,y=coo_y,col=lineage),lwd=0.5,alpha=0.7,show.legend=F)+
labs(x="x-coordinate (embryo : 2)",y="y-coordinate",title=paste0("Lineage: ABpl —> ",strsplit(focalCells[1],"_")[[1]][1]))

#geom_point(data=dataTemp[dataTemp$lineage==focalCells[1],],aes(x=coo_x,y=coo_y),size=1,col="#F8766D")+
#geom_path(data=dataTemp[dataTemp$lineage==focalCells[1],],aes(x=coo_x,y=coo_y),lwd=0.5,col="#F8766D")+
#geom_point(data=dataTemp[dataTemp$lineage==paste0(focalCells[1],"(WT)"),],aes(x=coo_x,y=coo_y),size=1,col="#00BFC4")+
#geom_path(data=dataTemp[dataTemp$lineage==paste0(focalCells[1],"(WT)"),],aes(x=coo_x,y=coo_y),lwd=0.5,col="#00BFC4")
#geom_point(data=dataTemp[dataTemp$lineage==focalCells[2],],aes(x=coo_x,y=coo_y),size=1,col="#00BFC4")+
#geom_path(data=dataTemp[dataTemp$lineage==focalCells[2],],aes(x=coo_x,y=coo_y),lwd=0.5,col="#00BFC4")+
#geom_point(data=dataTemp[dataTemp$lineage==paste0(focalCells[2],"(WT)"),],aes(x=coo_x,y=coo_y),size=1,col="#00BFC4")+
#geom_path(data=dataTemp[dataTemp$lineage==paste0(focalCells[2],"(WT)"),],aes(x=coo_x,y=coo_y),lwd=0.5,col="#00BFC4",lty=2)

ggsave("effectCells_PRP-4_spatial_distribution_emb2.pdf",width=unit(3.2,"in"),height=unit(2,"in"))



# emb : 3

focalEmbIndex=which(embryoGene==focalGeneVec[6])[3]
focalCells<-effectCells[which(abs(q1_cell350_alignment_xyz_noise[focalEmbIndex,which(coef(cvglm_residual)!=0)[-1]])>30)]
focalCells 

focalLineage<-lineage_info[4:10,sapply(do.call(rbind,strsplit(focalCells,"_"))[,1],function(x){rowColRec(which(lineage_info==x)[1],dim(lineage_info)[1],dim(lineage_info)[2])[2]})]
focal_lineage_x<-unlist(lapply(1:dim(focalLineage)[2],function(i){q1_cell350_alignment_xyz[focalEmbIndex,na.omit(match(paste0(focalLineage[,i],"_x"),colnames(q1_cell350_alignment_xyz)))]}))
contr_lineage_x<-unlist(lapply(1:dim(focalLineage)[2],function(i){coo_x[na.omit(match(paste0(focalLineage[,i],"_x"),colnames(q1_cell350_alignment_xyz)))]}))
focal_lineage_y<-unlist(lapply(1:dim(focalLineage)[2],function(i){q1_cell350_alignment_xyz[focalEmbIndex,na.omit(match(paste0(focalLineage[,i],"_y"),colnames(q1_cell350_alignment_xyz)))]}))
contr_lineage_y<-unlist(lapply(1:dim(focalLineage)[2],function(i){coo_y[na.omit(match(paste0(focalLineage[,i],"_x"),colnames(q1_cell350_alignment_xyz)))]}))

dataTemp<-data.frame(coo_x=c(coo_x,focal_lineage_x,contr_lineage_x),coo_y=c(coo_y,focal_lineage_y,contr_lineage_y),
lineage=c(rep("WT",length(coo_x)),unlist(lapply(1:length(focalCells),function(i){rep(focalCells[i],length(na.omit(match(paste0(focalLineage[,i],"_x"),colnames(q1_cell350_alignment_xyz)))))})),
	unlist(lapply(1:length(focalCells),function(i){rep(paste0(focalCells[i],"(WT)"),length(na.omit(match(paste0(focalLineage[,i],"_x"),colnames(q1_cell350_alignment_xyz)))))}))))


ggplot(data=dataTemp[dataTemp$lineage=="WT",],aes(x=coo_x,y=coo_y))+
geom_point(size=0.2,col="lightGrey")+
geom_point(data=dataTemp[dataTemp$lineage==focalCells[2]|dataTemp$lineage==paste0(focalCells[2],"(WT)"),],aes(x=coo_x,y=coo_y,col=lineage),size=0.75,alpha=0.7,show.legend=F)+
geom_path(data=dataTemp[dataTemp$lineage==focalCells[2]|dataTemp$lineage==paste0(focalCells[2],"(WT)"),],aes(x=coo_x,y=coo_y,col=lineage),lwd=0.5,alpha=0.7,show.legend=F)+
labs(x="x-coordinate (embryo : 3)",y="y-coordinate",title=paste0("Lineage: MS —> ",strsplit(focalCells[2],"_")[[1]][1]))

ggsave("effectCells_PRP-4_spatial_distribution_emb3.pdf",width=unit(3.2,"in"),height=unit(2,"in"))


## 32-cell embryo 

earlyCell<-unique(as.character(lineage_info[4:5,]))
dataTemp<-data.frame(coo_x=c(coo_x,coo_x[paste0(earlyCell,"_x")]),coo_y=c(coo_y,coo_y[paste0(earlyCell,"_y")]),
	lineage=c(rep("all",714),lineage_info[3,sapply(match(earlyCell,lineage_info),function(x){rowColRec(x,dim(lineage_info)[1],dim(lineage_info)[2])[2]})]))

ggplot(data=dataTemp[dataTemp$lineage=="all",],aes(x=coo_x,y=coo_y))+
geom_point(size=0.1,col="lightGrey")+
geom_point(data=dataTemp[dataTemp$lineage!="all",],aes(x=coo_x,y=coo_y,col=lineage),size=0.4,alpha=0.7,show.legend=F)+
labs(x="x-coordinate",y="y-coordinate",title=expression(paste(4^th," to ",5^th)))

ggsave("earlyCell_4_to_5.pdf",width=unit(1.5,"in"),height=unit(1.2,"in"))

	
	
earlyCell<-unique(as.character(lineage_info[4:6,]))
dataTemp<-data.frame(coo_x=c(coo_x,coo_x[paste0(earlyCell,"_x")]),coo_y=c(coo_y,coo_y[paste0(earlyCell,"_y")]),
	lineage=c(rep("all",714),lineage_info[3,sapply(match(earlyCell,lineage_info),function(x){rowColRec(x,dim(lineage_info)[1],dim(lineage_info)[2])[2]})]))

ggplot(data=dataTemp[dataTemp$lineage=="all",],aes(x=coo_x,y=coo_y))+
geom_point(size=0.1,col="lightGrey")+
geom_point(data=dataTemp[dataTemp$lineage!="all",],aes(x=coo_x,y=coo_y,col=lineage),size=0.4,alpha=0.7,show.legend=F)+
labs(x="x-coordinate",y="y-coordinate",title=expression(paste(4^th," to ",6^th)))

ggsave("earlyCell_4_to_6.pdf",width=unit(1.5,"in"),height=unit(1.2,"in"))




##################################################################################################
								End
##################################################################################################







### spatial distribution of effectCells for case 

#check for embryo 2 

focalEmbIndex=which(embryoGene==focalGeneVec[4])[2]
focalCells<-effectCells[which(abs(q1_cell350_alignment_xyz_noise[focalEmbIndex,which(coef(cvglm_residual)!=0)[-1]])>30)]


{

par(mfrow=c(3,2))
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpl")]))
plot(coo_x[paste0(tempCells,"_x")],coo_y[paste0(tempCells,"_y")])
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpr")]))
plot(coo_x[paste0(tempCells,"_x")],coo_y[paste0(tempCells,"_y")])

focal_coo_x<-q1_cell350_alignment_xyz[focalEmbIndex,1:714]
focal_coo_y<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714]
focal_coo_z<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714*2]

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpl")]))
plot(focal_coo_x[paste0(tempCells,"_x")],focal_coo_y[paste0(tempCells,"_y")])
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpr")]))
plot(focal_coo_x[paste0(tempCells,"_x")],focal_coo_y[paste0(tempCells,"_y")])

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpl")]))
plot(coo_x[paste0(tempCells,"_x")],focal_coo_x[paste0(tempCells,"_x")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpr")]))
plot(coo_x[paste0(tempCells,"_x")],focal_coo_x[paste0(tempCells,"_x")])
abline(0,1)

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpl")]))
plot(coo_y[paste0(tempCells,"_y")],focal_coo_y[paste0(tempCells,"_y")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpr")]))
plot(coo_y[paste0(tempCells,"_y")],focal_coo_y[paste0(tempCells,"_y")])
abline(0,1)

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpl")]))
plot(coo_z[paste0(tempCells,"_z")],focal_coo_z[paste0(tempCells,"_z")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpr")]))
plot(coo_z[paste0(tempCells,"_z")],focal_coo_z[paste0(tempCells,"_z")])
abline(0,1)



par(mfrow=c(3,2))
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABal")]))
plot(coo_x[paste0(tempCells,"_x")],coo_y[paste0(tempCells,"_y")])
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABar")]))
plot(coo_x[paste0(tempCells,"_x")],coo_y[paste0(tempCells,"_y")])

focal_coo_x<-q1_cell350_alignment_xyz[focalEmbIndex,1:714]
focal_coo_y<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714]
focal_coo_z<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714*2]

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABal")]))
plot(focal_coo_x[paste0(tempCells,"_x")],focal_coo_y[paste0(tempCells,"_y")])
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABar")]))
plot(focal_coo_x[paste0(tempCells,"_x")],focal_coo_y[paste0(tempCells,"_y")])

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABal")]))
plot(coo_x[paste0(tempCells,"_x")],focal_coo_x[paste0(tempCells,"_x")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABar")]))
plot(coo_x[paste0(tempCells,"_x")],focal_coo_x[paste0(tempCells,"_x")])
abline(0,1)

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABal")]))
plot(coo_y[paste0(tempCells,"_y")],focal_coo_y[paste0(tempCells,"_y")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABar")]))
plot(coo_y[paste0(tempCells,"_y")],focal_coo_y[paste0(tempCells,"_y")])
abline(0,1)

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABal")]))
plot(coo_z[paste0(tempCells,"_z")],focal_coo_z[paste0(tempCells,"_z")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABar")]))
plot(coo_z[paste0(tempCells,"_z")],focal_coo_z[paste0(tempCells,"_z")])
abline(0,1)




par(mfrow=c(3,2))
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="MS")]))
plot(coo_x[paste0(tempCells,"_x")],coo_y[paste0(tempCells,"_y")])
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="E")]))
plot(coo_x[paste0(tempCells,"_x")],coo_y[paste0(tempCells,"_y")])

focal_coo_x<-q1_cell350_alignment_xyz[focalEmbIndex,1:714]
focal_coo_y<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714]
focal_coo_z<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714*2]

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="MS")]))
plot(focal_coo_x[paste0(tempCells,"_x")],focal_coo_y[paste0(tempCells,"_y")])
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="E")]))
plot(focal_coo_x[paste0(tempCells,"_x")],focal_coo_y[paste0(tempCells,"_y")])

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="MS")]))
plot(coo_x[paste0(tempCells,"_x")],focal_coo_x[paste0(tempCells,"_x")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="E")]))
plot(coo_x[paste0(tempCells,"_x")],focal_coo_x[paste0(tempCells,"_x")])
abline(0,1)

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="MS")]))
plot(coo_y[paste0(tempCells,"_y")],focal_coo_y[paste0(tempCells,"_y")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="E")]))
plot(coo_y[paste0(tempCells,"_y")],focal_coo_y[paste0(tempCells,"_y")])
abline(0,1)

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="MS")]))
plot(coo_z[paste0(tempCells,"_z")],focal_coo_z[paste0(tempCells,"_z")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="E")]))
plot(coo_z[paste0(tempCells,"_z")],focal_coo_z[paste0(tempCells,"_z")])
abline(0,1)




par(mfrow=c(3,2))
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="C")]))
plot(coo_x[paste0(tempCells,"_x")],coo_y[paste0(tempCells,"_y")])
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="P3")]))
plot(coo_x[paste0(tempCells,"_x")],coo_y[paste0(tempCells,"_y")])

focal_coo_x<-q1_cell350_alignment_xyz[focalEmbIndex,1:714]
focal_coo_y<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714]
focal_coo_z<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714*2]

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="C")]))
plot(focal_coo_x[paste0(tempCells,"_x")],focal_coo_y[paste0(tempCells,"_y")])
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="P3")]))
plot(focal_coo_x[paste0(tempCells,"_x")],focal_coo_y[paste0(tempCells,"_y")])

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="C")]))
plot(coo_x[paste0(tempCells,"_x")],focal_coo_x[paste0(tempCells,"_x")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="P3")]))
plot(coo_x[paste0(tempCells,"_x")],focal_coo_x[paste0(tempCells,"_x")])
abline(0,1)

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="C")]))
plot(coo_y[paste0(tempCells,"_y")],focal_coo_y[paste0(tempCells,"_y")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="P3")]))
plot(coo_y[paste0(tempCells,"_y")],focal_coo_y[paste0(tempCells,"_y")])
abline(0,1)

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="C")]))
plot(coo_z[paste0(tempCells,"_z")],focal_coo_z[paste0(tempCells,"_z")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="P3")]))
plot(coo_z[paste0(tempCells,"_z")],focal_coo_z[paste0(tempCells,"_z")])
abline(0,1)


}


#check for embryo 3 

focalEmbIndex=which(embryoGene==focalGeneVec[4])[3]
focalCells<-effectCells[which(abs(q1_cell350_alignment_xyz_noise[focalEmbIndex,which(coef(cvglm_residual)!=0)[-1]])>30)]



{

par(mfrow=c(3,2))
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpl")]))
plot(coo_x[paste0(tempCells,"_x")],coo_y[paste0(tempCells,"_y")])
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpr")]))
plot(coo_x[paste0(tempCells,"_x")],coo_y[paste0(tempCells,"_y")])

focal_coo_x<-q1_cell350_alignment_xyz[focalEmbIndex,1:714]
focal_coo_y<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714]
focal_coo_z<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714*2]

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpl")]))
plot(focal_coo_x[paste0(tempCells,"_x")],focal_coo_y[paste0(tempCells,"_y")])
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpr")]))
plot(focal_coo_x[paste0(tempCells,"_x")],focal_coo_y[paste0(tempCells,"_y")])

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpl")]))
plot(coo_x[paste0(tempCells,"_x")],focal_coo_x[paste0(tempCells,"_x")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpr")]))
plot(coo_x[paste0(tempCells,"_x")],focal_coo_x[paste0(tempCells,"_x")])
abline(0,1)

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpl")]))
plot(coo_y[paste0(tempCells,"_y")],focal_coo_y[paste0(tempCells,"_y")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpr")]))
plot(coo_y[paste0(tempCells,"_y")],focal_coo_y[paste0(tempCells,"_y")])
abline(0,1)

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpl")]))
plot(coo_z[paste0(tempCells,"_z")],focal_coo_z[paste0(tempCells,"_z")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABpr")]))
plot(coo_z[paste0(tempCells,"_z")],focal_coo_z[paste0(tempCells,"_z")])
abline(0,1)



par(mfrow=c(3,2))
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABal")]))
plot(coo_x[paste0(tempCells,"_x")],coo_y[paste0(tempCells,"_y")])
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABar")]))
plot(coo_x[paste0(tempCells,"_x")],coo_y[paste0(tempCells,"_y")])

focal_coo_x<-q1_cell350_alignment_xyz[focalEmbIndex,1:714]
focal_coo_y<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714]
focal_coo_z<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714*2]

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABal")]))
plot(focal_coo_x[paste0(tempCells,"_x")],focal_coo_y[paste0(tempCells,"_y")])
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABar")]))
plot(focal_coo_x[paste0(tempCells,"_x")],focal_coo_y[paste0(tempCells,"_y")])

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABal")]))
plot(coo_x[paste0(tempCells,"_x")],focal_coo_x[paste0(tempCells,"_x")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABar")]))
plot(coo_x[paste0(tempCells,"_x")],focal_coo_x[paste0(tempCells,"_x")])
abline(0,1)

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABal")]))
plot(coo_y[paste0(tempCells,"_y")],focal_coo_y[paste0(tempCells,"_y")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABar")]))
plot(coo_y[paste0(tempCells,"_y")],focal_coo_y[paste0(tempCells,"_y")])
abline(0,1)

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABal")]))
plot(coo_z[paste0(tempCells,"_z")],focal_coo_z[paste0(tempCells,"_z")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="ABar")]))
plot(coo_z[paste0(tempCells,"_z")],focal_coo_z[paste0(tempCells,"_z")])
abline(0,1)




par(mfrow=c(3,2))
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="MS")]))
plot(coo_x[paste0(tempCells,"_x")],coo_y[paste0(tempCells,"_y")])
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="E")]))
plot(coo_x[paste0(tempCells,"_x")],coo_y[paste0(tempCells,"_y")])

focal_coo_x<-q1_cell350_alignment_xyz[focalEmbIndex,1:714]
focal_coo_y<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714]
focal_coo_z<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714*2]

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="MS")]))
plot(focal_coo_x[paste0(tempCells,"_x")],focal_coo_y[paste0(tempCells,"_y")])
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="E")]))
plot(focal_coo_x[paste0(tempCells,"_x")],focal_coo_y[paste0(tempCells,"_y")])

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="MS")]))
plot(coo_x[paste0(tempCells,"_x")],focal_coo_x[paste0(tempCells,"_x")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="E")]))
plot(coo_x[paste0(tempCells,"_x")],focal_coo_x[paste0(tempCells,"_x")])
abline(0,1)

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="MS")]))
plot(coo_y[paste0(tempCells,"_y")],focal_coo_y[paste0(tempCells,"_y")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="E")]))
plot(coo_y[paste0(tempCells,"_y")],focal_coo_y[paste0(tempCells,"_y")])
abline(0,1)

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="MS")]))
plot(coo_z[paste0(tempCells,"_z")],focal_coo_z[paste0(tempCells,"_z")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="E")]))
plot(coo_z[paste0(tempCells,"_z")],focal_coo_z[paste0(tempCells,"_z")])
abline(0,1)




par(mfrow=c(3,2))
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="C")]))
plot(coo_x[paste0(tempCells,"_x")],coo_y[paste0(tempCells,"_y")])
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="P3")]))
plot(coo_x[paste0(tempCells,"_x")],coo_y[paste0(tempCells,"_y")])

focal_coo_x<-q1_cell350_alignment_xyz[focalEmbIndex,1:714]
focal_coo_y<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714]
focal_coo_z<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714*2]

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="C")]))
plot(focal_coo_x[paste0(tempCells,"_x")],focal_coo_y[paste0(tempCells,"_y")])
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="P3")]))
plot(focal_coo_x[paste0(tempCells,"_x")],focal_coo_y[paste0(tempCells,"_y")])

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="C")]))
plot(coo_x[paste0(tempCells,"_x")],focal_coo_x[paste0(tempCells,"_x")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="P3")]))
plot(coo_x[paste0(tempCells,"_x")],focal_coo_x[paste0(tempCells,"_x")])
abline(0,1)

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="C")]))
plot(coo_y[paste0(tempCells,"_y")],focal_coo_y[paste0(tempCells,"_y")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="P3")]))
plot(coo_y[paste0(tempCells,"_y")],focal_coo_y[paste0(tempCells,"_y")])
abline(0,1)

tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="C")]))
plot(coo_z[paste0(tempCells,"_z")],focal_coo_z[paste0(tempCells,"_z")])
abline(0,1)
tempCells<-unique(c(lineage_info[5:10,which(lineage_info[4,]=="P3")]))
plot(coo_z[paste0(tempCells,"_z")],focal_coo_z[paste0(tempCells,"_z")])
abline(0,1)


}


































focalCellIndex<-match(c("ABaraaapaa_z","ABarapapaa_z","ABarpaapa_z","ABarppapp_z"),colnames(q1_cell350_alignment_xyz_noise_z))

match(c("ABaraaapaa_z","ABarapapaa_z","ABarpaapa_z","ABarppapp_z"),colnames(q1_cell350_alignment_xyz_noise_z))
focalSisterIndex<-match(c("ABaraaapap_z","ABarapapap_z","ABarpaapp_z","ABarppapa_z"),colnames(q1_cell350_alignment_xyz_noise_z))
focalMoIndex<-match(c("ABaraaapa_z","ABarapapa_z","ABarpaap_z","ABarppap_z"),colnames(q1_cell350_alignment_xyz_noise_z))


focal_coo_x<-q1_cell350_alignment_xyz[focalEmbIndex,1:714]
focal_coo_y<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714]
focal_coo_z<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714*2]


dataTemp<-data.frame(coo_x=c(focal_coo_x[-focalCellIndex],focal_coo_x[focalCellIndex],coo_x[focalCellIndex]),
	coo_z=c(focal_coo_z[-focalCellIndex],focal_coo_z[focalCellIndex],coo_z[focalCellIndex]),
	type=factor(c(rep("other cells",length(focal_coo_x[-focalCellIndex])),rep("C1:C4 (expected)",length(focalCellIndex)),rep("C1:C4 (observed)",length(focalCellIndex))),level=c("other cells","C1:C4 (expected)","C1:C4 (observed)")))

dataTemp_C1<-data.frame(coo_x2=coo_x[focalCellIndex],focal_coo_x2=focal_coo_x[focalCellIndex],coo_z2=coo_z[focalCellIndex],focal_coo_z2=focal_coo_z[focalCellIndex])[1,]
dataTemp_C2<-data.frame(coo_x2=coo_x[focalCellIndex],focal_coo_x2=focal_coo_x[focalCellIndex],coo_z2=coo_z[focalCellIndex],focal_coo_z2=focal_coo_z[focalCellIndex])[2,]
dataTemp_C3<-data.frame(coo_x2=coo_x[focalCellIndex],focal_coo_x2=focal_coo_x[focalCellIndex],coo_z2=coo_z[focalCellIndex],focal_coo_z2=focal_coo_z[focalCellIndex])[3,]
dataTemp_C4<-data.frame(coo_x2=coo_x[focalCellIndex],focal_coo_x2=focal_coo_x[focalCellIndex],coo_z2=coo_z[focalCellIndex],focal_coo_z2=focal_coo_z[focalCellIndex])[4,]

dataTemp_S1<-data.frame(coo_x2=coo_x[focalSisterIndex],focal_coo_x2=focal_coo_x[focalSisterIndex],coo_z2=coo_z[focalSisterIndex],focal_coo_z2=focal_coo_z[focalSisterIndex])[1,]
dataTemp_S2<-data.frame(coo_x2=coo_x[focalSisterIndex],focal_coo_x2=focal_coo_x[focalSisterIndex],coo_z2=coo_z[focalSisterIndex],focal_coo_z2=focal_coo_z[focalSisterIndex])[2,]
dataTemp_S3<-data.frame(coo_x2=coo_x[focalSisterIndex],focal_coo_x2=focal_coo_x[focalSisterIndex],coo_z2=coo_z[focalSisterIndex],focal_coo_z2=focal_coo_z[focalSisterIndex])[3,]
dataTemp_S4<-data.frame(coo_x2=coo_x[focalSisterIndex],focal_coo_x2=focal_coo_x[focalSisterIndex],coo_z2=coo_z[focalSisterIndex],focal_coo_z2=focal_coo_z[focalSisterIndex])[4,]

dataTemp_M1<-data.frame(coo_x2=coo_x[focalMoIndex],focal_coo_x2=focal_coo_x[focalMoIndex],coo_z2=coo_z[focalMoIndex],focal_coo_z2=focal_coo_z[focalMoIndex])[1,]
dataTemp_M2<-data.frame(coo_x2=coo_x[focalMoIndex],focal_coo_x2=focal_coo_x[focalMoIndex],coo_z2=coo_z[focalMoIndex],focal_coo_z2=focal_coo_z[focalMoIndex])[2,]
dataTemp_M3<-data.frame(coo_x2=coo_x[focalMoIndex],focal_coo_x2=focal_coo_x[focalMoIndex],coo_z2=coo_z[focalMoIndex],focal_coo_z2=focal_coo_z[focalMoIndex])[3,]
dataTemp_M4<-data.frame(coo_x2=coo_x[focalMoIndex],focal_coo_x2=focal_coo_x[focalMoIndex],coo_z2=coo_z[focalMoIndex],focal_coo_z2=focal_coo_z[focalMoIndex])[4,]

ggplot(data=dataTemp,aes(x=coo_x,y=coo_z,col=type))+
geom_point(alpha=1,show.legend=TRUE)+
scale_y_continuous(breaks=c(-40,0,40))+
geom_segment(data=dataTemp_C1,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_C2,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_C3,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_C4,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_S1,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment2"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_S2,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment2"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_S3,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment2"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_S4,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment2"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_M1,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment3"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_M2,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment3"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_M3,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment3"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_M4,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment3"),arrow=arrow(length=unit(0.1,"in")))+
scale_color_manual(values=c("grey39","#F8766D","lightGrey","#F8766D","blue","orange"))+
labs(x="x-coordinate",y="z-coordinate",title="SEC-12 : embryo : 2")
























# spatial distribution of effectCells for case 

focalEmbIndex=which(embryoGene==focalGene)[2]
focalCellIndex<-match(c("ABaraaapaa_z","ABarapapaa_z","ABarpaapa_z","ABarppapp_z"),colnames(q1_cell350_alignment_xyz_noise_z))

match(c("ABaraaapaa_z","ABarapapaa_z","ABarpaapa_z","ABarppapp_z"),colnames(q1_cell350_alignment_xyz_noise_z))
focalSisterIndex<-match(c("ABaraaapap_z","ABarapapap_z","ABarpaapp_z","ABarppapa_z"),colnames(q1_cell350_alignment_xyz_noise_z))
focalMoIndex<-match(c("ABaraaapa_z","ABarapapa_z","ABarpaap_z","ABarppap_z"),colnames(q1_cell350_alignment_xyz_noise_z))

coo_x<-apply(q1_cell350_alignment_xyz[,1:714],2,mean)
coo_y<-apply(q1_cell350_alignment_xyz[,c(1:714)+714],2,mean)
coo_z<-apply(q1_cell350_alignment_xyz[,c(1:714)+714*2],2,mean)

focal_coo_x<-q1_cell350_alignment_xyz[focalEmbIndex,1:714]
focal_coo_y<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714]
focal_coo_z<-q1_cell350_alignment_xyz[focalEmbIndex,c(1:714)+714*2]


dataTemp<-data.frame(coo_x=c(focal_coo_x[-focalCellIndex],focal_coo_x[focalCellIndex],coo_x[focalCellIndex]),
	coo_z=c(focal_coo_z[-focalCellIndex],focal_coo_z[focalCellIndex],coo_z[focalCellIndex]),
	type=factor(c(rep("other cells",length(focal_coo_x[-focalCellIndex])),rep("C1:C4 (expected)",length(focalCellIndex)),rep("C1:C4 (observed)",length(focalCellIndex))),level=c("other cells","C1:C4 (expected)","C1:C4 (observed)")))

dataTemp_C1<-data.frame(coo_x2=coo_x[focalCellIndex],focal_coo_x2=focal_coo_x[focalCellIndex],coo_z2=coo_z[focalCellIndex],focal_coo_z2=focal_coo_z[focalCellIndex])[1,]
dataTemp_C2<-data.frame(coo_x2=coo_x[focalCellIndex],focal_coo_x2=focal_coo_x[focalCellIndex],coo_z2=coo_z[focalCellIndex],focal_coo_z2=focal_coo_z[focalCellIndex])[2,]
dataTemp_C3<-data.frame(coo_x2=coo_x[focalCellIndex],focal_coo_x2=focal_coo_x[focalCellIndex],coo_z2=coo_z[focalCellIndex],focal_coo_z2=focal_coo_z[focalCellIndex])[3,]
dataTemp_C4<-data.frame(coo_x2=coo_x[focalCellIndex],focal_coo_x2=focal_coo_x[focalCellIndex],coo_z2=coo_z[focalCellIndex],focal_coo_z2=focal_coo_z[focalCellIndex])[4,]

dataTemp_S1<-data.frame(coo_x2=coo_x[focalSisterIndex],focal_coo_x2=focal_coo_x[focalSisterIndex],coo_z2=coo_z[focalSisterIndex],focal_coo_z2=focal_coo_z[focalSisterIndex])[1,]
dataTemp_S2<-data.frame(coo_x2=coo_x[focalSisterIndex],focal_coo_x2=focal_coo_x[focalSisterIndex],coo_z2=coo_z[focalSisterIndex],focal_coo_z2=focal_coo_z[focalSisterIndex])[2,]
dataTemp_S3<-data.frame(coo_x2=coo_x[focalSisterIndex],focal_coo_x2=focal_coo_x[focalSisterIndex],coo_z2=coo_z[focalSisterIndex],focal_coo_z2=focal_coo_z[focalSisterIndex])[3,]
dataTemp_S4<-data.frame(coo_x2=coo_x[focalSisterIndex],focal_coo_x2=focal_coo_x[focalSisterIndex],coo_z2=coo_z[focalSisterIndex],focal_coo_z2=focal_coo_z[focalSisterIndex])[4,]

dataTemp_M1<-data.frame(coo_x2=coo_x[focalMoIndex],focal_coo_x2=focal_coo_x[focalMoIndex],coo_z2=coo_z[focalMoIndex],focal_coo_z2=focal_coo_z[focalMoIndex])[1,]
dataTemp_M2<-data.frame(coo_x2=coo_x[focalMoIndex],focal_coo_x2=focal_coo_x[focalMoIndex],coo_z2=coo_z[focalMoIndex],focal_coo_z2=focal_coo_z[focalMoIndex])[2,]
dataTemp_M3<-data.frame(coo_x2=coo_x[focalMoIndex],focal_coo_x2=focal_coo_x[focalMoIndex],coo_z2=coo_z[focalMoIndex],focal_coo_z2=focal_coo_z[focalMoIndex])[3,]
dataTemp_M4<-data.frame(coo_x2=coo_x[focalMoIndex],focal_coo_x2=focal_coo_x[focalMoIndex],coo_z2=coo_z[focalMoIndex],focal_coo_z2=focal_coo_z[focalMoIndex])[4,]

ggplot(data=dataTemp,aes(x=coo_x,y=coo_z,col=type))+
geom_point(alpha=1,show.legend=TRUE)+
scale_y_continuous(breaks=c(-40,0,40))+
geom_segment(data=dataTemp_C1,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_C2,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_C3,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_C4,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_S1,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment2"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_S2,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment2"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_S3,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment2"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_S4,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment2"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_M1,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment3"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_M2,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment3"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_M3,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment3"),arrow=arrow(length=unit(0.1,"in")))+
geom_segment(data=dataTemp_M4,aes(xend=coo_x2,yend=coo_z2,x=focal_coo_x2,y=focal_coo_z2,colour="segment3"),arrow=arrow(length=unit(0.1,"in")))+
scale_color_manual(values=c("grey39","#F8766D","lightGrey","#F8766D","blue","orange"))+
labs(x="x-coordinate",y="z-coordinate",title="SEC-12 : embryo : 2")


ggsave(".pdf",width=unit(3.6,"in"),height=unit(1.5,"in"))

###

rowColRec(match("ABarpaapa",lineage_info),dim(lineage_info)[1],dim(lineage_info)[2])

paste0(lineage_info[4:10,121],"_z")

focal_coo<-data.frame(x=q1_cell350_alignment_xyz[focalEmbIndex,paste0(lineage_info[4:10,121],"_x")],
	y=q1_cell350_alignment_xyz[focalEmbIndex,paste0(lineage_info[4:10,121],"_y")],
	z=q1_cell350_alignment_xyz[focalEmbIndex,paste0(lineage_info[4:10,121],"_z")])

focal_coo_sister<-data.frame(x=q1_cell350_alignment_xyz[focalEmbIndex,paste0(lineage_info[4:10,102],"_x")],
	y=q1_cell350_alignment_xyz[focalEmbIndex,paste0(lineage_info[4:10,102],"_y")],
	z=q1_cell350_alignment_xyz[focalEmbIndex,paste0(lineage_info[4:10,102],"_z")])

focal_coo_contr<-data.frame(x=q1_cell350_alignment_xyz[1,paste0(lineage_info[4:10,121],"_x")],
	y=q1_cell350_alignment_xyz[1,paste0(lineage_info[4:10,121],"_y")],
	z=q1_cell350_alignment_xyz[1,paste0(lineage_info[4:10,121],"_z")])

focal_coo_contr2<-data.frame(x=q1_cell350_alignment_xyz[1,paste0(lineage_info[4:10,102],"_x")],
	y=q1_cell350_alignment_xyz[1,paste0(lineage_info[4:10,102],"_y")],
	z=q1_cell350_alignment_xyz[1,paste0(lineage_info[4:10,102],"_z")])

	
#

focal_coo_abarapapaa<-data.frame(x=q1_cell350_alignment_xyz[focalEmbIndex,paste0(lineage_info[4:10,39],"_x")],
	y=q1_cell350_alignment_xyz[focalEmbIndex,paste0(lineage_info[4:10,39],"_y")],
	z=q1_cell350_alignment_xyz[focalEmbIndex,paste0(lineage_info[4:10,39],"_z")])

focal_coo_contr_abarapapaa<-data.frame(x=q1_cell350_alignment_xyz[1,paste0(lineage_info[4:10,39],"_x")],
	y=q1_cell350_alignment_xyz[1,paste0(lineage_info[4:10,39],"_y")],
	z=q1_cell350_alignment_xyz[1,paste0(lineage_info[4:10,39],"_z")])

#
focal_coo_ABarpaapa<-data.frame(x=q1_cell350_alignment_xyz[focalEmbIndex,paste0(lineage_info[4:10,428],"_x")],
	y=q1_cell350_alignment_xyz[focalEmbIndex,paste0(lineage_info[4:10,428],"_y")],
	z=q1_cell350_alignment_xyz[focalEmbIndex,paste0(lineage_info[4:10,428],"_z")])

focal_coo_contr_ABarpaapa<-data.frame(x=q1_cell350_alignment_xyz[1,paste0(lineage_info[4:10,428],"_x")],
	y=q1_cell350_alignment_xyz[1,paste0(lineage_info[4:10,428],"_y")],
	z=q1_cell350_alignment_xyz[1,paste0(lineage_info[4:10,428],"_z")])

##

par(mfrow=c(2,3))

plot(q1_cell350_alignment_xyz[1,match(c(paste0(lineage_info[10,lineage_info[5,]=="ABara"],"_x"),paste0(lineage_info[10,lineage_info[5,]=="ABarp"],"_x")),colnames(q1_cell350_alignment_xyz))],
q1_cell350_alignment_xyz[focalEmbIndex,match(c(paste0(lineage_info[10,lineage_info[5,]=="ABara"],"_x"),paste0(lineage_info[10,lineage_info[5,]=="ABarp"],"_x")),colnames(q1_cell350_alignment_xyz))])

plot(q1_cell350_alignment_xyz[1,match(c(paste0(lineage_info[10,lineage_info[5,]=="ABara"],"_y"),paste0(lineage_info[10,lineage_info[5,]=="ABarp"],"_y")),colnames(q1_cell350_alignment_xyz))],
q1_cell350_alignment_xyz[focalEmbIndex,match(c(paste0(lineage_info[10,lineage_info[5,]=="ABara"],"_y"),paste0(lineage_info[10,lineage_info[5,]=="ABarp"],"_y")),colnames(q1_cell350_alignment_xyz))])

plot(q1_cell350_alignment_xyz[1,match(c(paste0(lineage_info[10,lineage_info[5,]=="ABara"],"_z"),paste0(lineage_info[10,lineage_info[5,]=="ABarp"],"_z")),colnames(q1_cell350_alignment_xyz))],
q1_cell350_alignment_xyz[focalEmbIndex,match(c(paste0(lineage_info[10,lineage_info[5,]=="ABara"],"_z"),paste0(lineage_info[10,lineage_info[5,]=="ABarp"],"_z")),colnames(q1_cell350_alignment_xyz))])


plot(q1_cell350_alignment_xyz[1,match(c(paste0(lineage_info[10,lineage_info[5,]=="ABala"],"_x"),paste0(lineage_info[10,lineage_info[5,]=="ABalp"],"_x")),colnames(q1_cell350_alignment_xyz))],
q1_cell350_alignment_xyz[focalEmbIndex,match(c(paste0(lineage_info[10,lineage_info[5,]=="ABala"],"_x"),paste0(lineage_info[10,lineage_info[5,]=="ABalp"],"_x")),colnames(q1_cell350_alignment_xyz))])

plot(q1_cell350_alignment_xyz[1,match(c(paste0(lineage_info[10,lineage_info[5,]=="ABala"],"_y"),paste0(lineage_info[10,lineage_info[5,]=="ABalp"],"_y")),colnames(q1_cell350_alignment_xyz))],
q1_cell350_alignment_xyz[focalEmbIndex,match(c(paste0(lineage_info[10,lineage_info[5,]=="ABala"],"_y"),paste0(lineage_info[10,lineage_info[5,]=="ABalp"],"_y")),colnames(q1_cell350_alignment_xyz))])

plot(q1_cell350_alignment_xyz[1,match(c(paste0(lineage_info[10,lineage_info[5,]=="ABala"],"_z"),paste0(lineage_info[10,lineage_info[5,]=="ABalp"],"_z")),colnames(q1_cell350_alignment_xyz))],
q1_cell350_alignment_xyz[focalEmbIndex,match(c(paste0(lineage_info[10,lineage_info[5,]=="ABala"],"_z"),paste0(lineage_info[10,lineage_info[5,]=="ABalp"],"_z")),colnames(q1_cell350_alignment_xyz))])




plot(q1_cell350_alignment_xyz[which(embryoGene==focalGeneVec[4])[3],],q1_cell350_alignment_xyz[1,])


plot(q1_cell350_alignment_xyz[which(embryoGene=="SKR-2")[3],1:714],q1_cell350_alignment_xyz[1,1:714])
plot(q1_cell350_alignment_xyz[which(embryoGene=="SKR-2")[3],c(1:714)+714],q1_cell350_alignment_xyz[1,c(1:714)+714])
plot(q1_cell350_alignment_xyz[which(embryoGene=="SKR-2")[3],c(1:714)+714*2],q1_cell350_alignment_xyz[1,c(1:714)+714*2])


plot(q1_cell350_alignment_xyz[which(embryoGene=="SEC-12")[2],c(1:714)+714*2],q1_cell350_alignment_xyz[1,c(1:714)+714*2])












cvglm_residual_rocs<-
cvglm_residual_training_rocs<-
cvglm_residual_test_rocs<-


predDF<-data.frame(Class1=predict(cvglm_residual_training,residual_feature_test,s="lambda.min")[,1])
predDF$truth<-rep("Class1",length(embryo_lethal_test))
predDF$truth[which(embryo_lethal_test==0)]="Class2"
predDF$truth<-as.factor(predDF$truth)

roc_auc(predDF,truth,Class1)[,3]




logODDs<-cvglm_wormbase[[i]]$fit.preval[,which(cvglm_wormbase[[i]]$lambda==cvglm_wormbase[[i]]$lambda.min)]

preidctProb<-exp(logODDs)/(1+exp(logODDs))

truth<-celegans_pheno[which(!is.na(celegans_pheno[,i])),i]
truth[which(truth=="0")]="Class2"
truth[which(truth=="1")]="Class1"
truth<-factor(truth,level=c("Class1","Class2"))

dataTemp<-data.frame(truth,Class1=preidctProb)

prauc[i]<-as.numeric(pr_auc(dataTemp, truth, Class1)[,3])







