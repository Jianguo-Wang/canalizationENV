

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


library(mnormt)
library(ggplot2)



focal_projectFolderPath<-r"(F:\wangj\Desktop\temp\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/data/canalGene"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_size.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData"
	),projectFolderPath=focal_projectFolderPath)



###

q1_cell350_alignment_xyz_noise_x<-q1_cell350_alignment_xyz_noise[1:1928,1:714]
q1_cell350_alignment_xyz_noise_y<-q1_cell350_alignment_xyz_noise[1:1928,(1+714):(714*2)]
q1_cell350_alignment_xyz_noise_z<-q1_cell350_alignment_xyz_noise[1:1928,(1+714*2):(714*3)]

#

feature_names<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise[1:1928,1:714]),"_x"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]


#####################################################################################
                           x 
#####################################################################################

###

X_mo<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo


### mutants with deviation from population mean of WT

embryo_geneName<-do.call(rbind,strsplit(rownames(delta_X),"_"))[,1]


#

mean_WT_x<-apply(q1_cell350_alignment_xyz[1:105,1:714],2,mean)
dist_with_WT_mean_x<-sqrt(sapply(1:1928,function(i){sum((mean_WT_x-q1_cell350_alignment_xyz[i,1:714])^2)}))


cutoff_1_x<-qnorm(1-0.05/(1928-105),mean(dist_with_WT_mean_x[1:105]),sd(dist_with_WT_mean_x[1:105]))
sigEmbryoInd_dist_with_WT_mean_x<-which(dist_with_WT_mean_x>cutoff_1_x)


##

sigEmbryoInd_dist_with_WT_mean_summary_x=data.frame(
	embryoInd=sigEmbryoInd_dist_with_WT_mean_x,
	embryoID=rownames(delta_X)[sigEmbryoInd_dist_with_WT_mean_x],
	geneName=embryo_geneName[sigEmbryoInd_dist_with_WT_mean_x],
	totalEmbryoNum=as.numeric(table(embryo_geneName)[match(embryo_geneName[sigEmbryoInd_dist_with_WT_mean_x],names(table(embryo_geneName)))]),
	sigEmbryoNum=as.numeric(table(embryo_geneName[sigEmbryoInd_dist_with_WT_mean_x])[match(embryo_geneName[sigEmbryoInd_dist_with_WT_mean_x],names(table(embryo_geneName[sigEmbryoInd_dist_with_WT_mean_x])))]),
	dist_with_WT_mean=dist_with_WT_mean_x[sigEmbryoInd_dist_with_WT_mean_x]
)

sigEmbryoInd_dist_with_WT_mean_summary_x$ratio<-round(sigEmbryoInd_dist_with_WT_mean_summary_x$sigEmbryoNum/sigEmbryoInd_dist_with_WT_mean_summary_x$totalEmbryoNum,2)

sigEmbryoInd_dist_with_WT_mean_summary_x_replicate<-sigEmbryoInd_dist_with_WT_mean_summary_x[sigEmbryoInd_dist_with_WT_mean_summary_x$ratio>0.5,]
sigEmbryoInd_dist_with_WT_mean_summary_x_replicate<-sigEmbryoInd_dist_with_WT_mean_summary_x_replicate[order(sigEmbryoInd_dist_with_WT_mean_summary_x_replicate$ratio,sigEmbryoInd_dist_with_WT_mean_summary_x_replicate$totalEmbryoNum,decreasing=T),]

sigGeneName_dist_with_WT_mean_x<-unique(sigEmbryoInd_dist_with_WT_mean_summary_x_replicate$geneName)


### mother model 

delta_X_pred_mo<-sapply(1:dim(delta_X)[2],function(i){
	return(cbind(rep(1,dim(delta_X)[1]),X_mo[,i])%*%
		summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[,1])
})

Res2_mo_x<-sapply(1:1928,function(i){sum((delta_X[i,]-delta_X_pred_mo[i,])^2)/sum(delta_X[i,]^2)})
sigEmbryoInd_dist_with_WT_mean_summary_x_replicate$Res2_mo_x<-round(Res2_mo_x[sigEmbryoInd_dist_with_WT_mean_summary_x_replicate$embryoInd],2)

cutoff_2_x<-qnorm(1-0.05/length(sigEmbryoInd_dist_with_WT_mean_summary_x_replicate$embryoInd),mean(Res2_mo_x[1:105]),sd(Res2_mo_x[1:105]))

sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_x<-sigEmbryoInd_dist_with_WT_mean_summary_x_replicate[sigEmbryoInd_dist_with_WT_mean_summary_x_replicate$Res2_mo_x>cutoff_2_x,]


sigGeneName_mother_defect_x<-unique(sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_x$geneName)


write.csv(sigGeneName_mother_defect_x,"sigGeneName_mother_defect_x.csv")


#####################################################################################
                           y 
#####################################################################################

###y

X_mo<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo

### mutants with deviation from population mean of WT

embryo_geneName<-do.call(rbind,strsplit(rownames(delta_X),"_"))[,1]


#

mean_WT_y<-apply(q1_cell350_alignment_xyz[1:105,c((1+714):(714*2))],2,mean)
dist_with_WT_mean_y<-sqrt(sapply(1:1928,function(i){sum((mean_WT_y-q1_cell350_alignment_xyz[i,c((1+714):(714*2))])^2)}))

#1-pnorm(max(dist_with_WT_mean[1:105]),mean(dist_with_WT_mean[1:105]),sd(dist_with_WT_mean[1:105]))=4.032075e-07<0.05/(1928-105)
#FDR=length(which(dist_with_WT_mean[1:105]>cutoff_1))/105<1%

cutoff_1_y<-qnorm(1-0.05/(1928-105),mean(dist_with_WT_mean_y[1:105]),sd(dist_with_WT_mean_y[1:105]))
sigEmbryoInd_dist_with_WT_mean_y<-which(dist_with_WT_mean_y>cutoff_1_y)


##

sigEmbryoInd_dist_with_WT_mean_summary_y=data.frame(
	embryoInd=sigEmbryoInd_dist_with_WT_mean_y,
	embryoID=rownames(delta_X)[sigEmbryoInd_dist_with_WT_mean_y],
	geneName=embryo_geneName[sigEmbryoInd_dist_with_WT_mean_y],
	totalEmbryoNum=as.numeric(table(embryo_geneName)[match(embryo_geneName[sigEmbryoInd_dist_with_WT_mean_y],names(table(embryo_geneName)))]),
	sigEmbryoNum=as.numeric(table(embryo_geneName[sigEmbryoInd_dist_with_WT_mean_y])[match(embryo_geneName[sigEmbryoInd_dist_with_WT_mean_y],names(table(embryo_geneName[sigEmbryoInd_dist_with_WT_mean_y])))]),
	dist_with_WT_mean=dist_with_WT_mean_y[sigEmbryoInd_dist_with_WT_mean_y]
)

sigEmbryoInd_dist_with_WT_mean_summary_y$ratio<-round(sigEmbryoInd_dist_with_WT_mean_summary_y$sigEmbryoNum/sigEmbryoInd_dist_with_WT_mean_summary_y$totalEmbryoNum,2)

sigEmbryoInd_dist_with_WT_mean_summary_y_replicate<-sigEmbryoInd_dist_with_WT_mean_summary_y[sigEmbryoInd_dist_with_WT_mean_summary_y$ratio>0.5,]
sigEmbryoInd_dist_with_WT_mean_summary_y_replicate<-sigEmbryoInd_dist_with_WT_mean_summary_y_replicate[order(sigEmbryoInd_dist_with_WT_mean_summary_y_replicate$ratio,sigEmbryoInd_dist_with_WT_mean_summary_y_replicate$totalEmbryoNum,decreasing=T),]

sigGeneName_dist_with_WT_mean_y<-unique(sigEmbryoInd_dist_with_WT_mean_summary_y_replicate$geneName)


### mother model 

delta_X_pred_mo<-sapply(1:dim(delta_X)[2],function(i){
	return(cbind(rep(1,dim(delta_X)[1]),X_mo[,i])%*%
		summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[,1])
})

Res2_mo_y<-sapply(1:1928,function(i){sum((delta_X[i,]-delta_X_pred_mo[i,])^2)/sum(delta_X[i,]^2)})
sigEmbryoInd_dist_with_WT_mean_summary_y_replicate$Res2_mo_y<-round(Res2_mo_y[sigEmbryoInd_dist_with_WT_mean_summary_y_replicate$embryoInd],2)

cutoff_2_y<-qnorm(1-0.05/length(sigEmbryoInd_dist_with_WT_mean_summary_y_replicate$embryoInd),mean(Res2_mo_y[1:105]),sd(Res2_mo_y[1:105]))

sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_y<-sigEmbryoInd_dist_with_WT_mean_summary_y_replicate[sigEmbryoInd_dist_with_WT_mean_summary_y_replicate$Res2_mo_y>cutoff_2_y,]


sigGeneName_mother_defect_y<-unique(sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_y$geneName)


write.csv(sigGeneName_mother_defect_y,"sigGeneName_mother_defect_y.csv")


#####################################################################################
                           z 
#####################################################################################

###

X_mo<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo

### mutants with deviation from population mean of WT

embryo_geneName<-do.call(rbind,strsplit(rownames(delta_X),"_"))[,1]


#

mean_WT_z<-apply(q1_cell350_alignment_xyz[1:105,c((1+714*2):(714*3))],2,mean)
dist_with_WT_mean_z<-sqrt(sapply(1:1928,function(i){sum((mean_WT_z-q1_cell350_alignment_xyz[i,c((1+714*2):(714*3))])^2)}))

#1-pnorm(max(dist_with_WT_mean[1:105]),mean(dist_with_WT_mean[1:105]),sd(dist_with_WT_mean[1:105]))=4.032075e-07<0.05/(1928-105)
#FDR=length(which(dist_with_WT_mean[1:105]>cutoff_1))/105<1%

cutoff_1_z<-qnorm(1-0.05/(1928-105),mean(dist_with_WT_mean_z[1:105]),sd(dist_with_WT_mean_z[1:105]))
sigEmbryoInd_dist_with_WT_mean_z<-which(dist_with_WT_mean_z>cutoff_1_z)


##

sigEmbryoInd_dist_with_WT_mean_summary_z=data.frame(
	embryoInd=sigEmbryoInd_dist_with_WT_mean_z,
	embryoID=rownames(delta_X)[sigEmbryoInd_dist_with_WT_mean_z],
	geneName=embryo_geneName[sigEmbryoInd_dist_with_WT_mean_z],
	totalEmbryoNum=as.numeric(table(embryo_geneName)[match(embryo_geneName[sigEmbryoInd_dist_with_WT_mean_z],names(table(embryo_geneName)))]),
	sigEmbryoNum=as.numeric(table(embryo_geneName[sigEmbryoInd_dist_with_WT_mean_z])[match(embryo_geneName[sigEmbryoInd_dist_with_WT_mean_z],names(table(embryo_geneName[sigEmbryoInd_dist_with_WT_mean_z])))]),
	dist_with_WT_mean=dist_with_WT_mean_z[sigEmbryoInd_dist_with_WT_mean_z]
)

sigEmbryoInd_dist_with_WT_mean_summary_z$ratio<-round(sigEmbryoInd_dist_with_WT_mean_summary_z$sigEmbryoNum/sigEmbryoInd_dist_with_WT_mean_summary_z$totalEmbryoNum,2)

sigEmbryoInd_dist_with_WT_mean_summary_z_replicate<-sigEmbryoInd_dist_with_WT_mean_summary_z[sigEmbryoInd_dist_with_WT_mean_summary_z$ratio>0.5,]
sigEmbryoInd_dist_with_WT_mean_summary_z_replicate<-sigEmbryoInd_dist_with_WT_mean_summary_z_replicate[order(sigEmbryoInd_dist_with_WT_mean_summary_z_replicate$ratio,sigEmbryoInd_dist_with_WT_mean_summary_z_replicate$totalEmbryoNum,decreasing=T),]

sigGeneName_dist_with_WT_mean_z<-unique(sigEmbryoInd_dist_with_WT_mean_summary_z_replicate$geneName)


### mother model 

delta_X_pred_mo<-sapply(1:dim(delta_X)[2],function(i){
	return(cbind(rep(1,dim(delta_X)[1]),X_mo[,i])%*%
		summary(lm(delta_X[,i]~X_mo[,i]))$coefficients[,1])
})

Res2_mo_z<-sapply(1:1928,function(i){sum((delta_X[i,]-delta_X_pred_mo[i,])^2)/sum(delta_X[i,]^2)})
sigEmbryoInd_dist_with_WT_mean_summary_z_replicate$Res2_mo_z<-round(Res2_mo_z[sigEmbryoInd_dist_with_WT_mean_summary_z_replicate$embryoInd],2)

cutoff_2_z<-qnorm(1-0.05/length(sigEmbryoInd_dist_with_WT_mean_summary_z_replicate$embryoInd),mean(Res2_mo_z[1:105]),sd(Res2_mo_z[1:105]))

sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_z<-sigEmbryoInd_dist_with_WT_mean_summary_z_replicate[sigEmbryoInd_dist_with_WT_mean_summary_z_replicate$Res2_mo_z>cutoff_2_z,]


sigGeneName_mother_defect_z<-unique(sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_z$geneName)


write.csv(sigGeneName_mother_defect_z,"sigGeneName_mother_defect_z.csv")




#################################################################################################

#################################################################################################



### spatial_mutant_res_x

spatial_mutant_res_x<-data.frame(
	EmbryoID=sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_x$embryoID,
	Gene_name=sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_x$geneName,
	Distance_from_average=sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_x$dist_with_WT_mean,
	Distance_from_average_pvalue=1-pnorm(sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_x$dist_with_WT_mean,mean(dist_with_WT_mean_x[1:105]),sd(dist_with_WT_mean_x[1:105])),
	SSE_SST=sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_x$Res2_mo_x,
	SSE_SST_pvalue=1-pnorm(sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_x$Res2_mo_x,mean(Res2_mo_x[1:105]),sd(Res2_mo_x[1:105]))
)


write.csv(spatial_mutant_res_x,"Tab.S2-spatial_mutant_res_x.csv",row.names=FALSE)


### spatial_mutant_res_y

spatial_mutant_res_y<-data.frame(
	EmbryoID=sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_y$embryoID,
	Gene_name=sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_y$geneName,
	Distance_from_average=sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_y$dist_with_WT_mean,
	Distance_from_average_pvalue=1-pnorm(sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_y$dist_with_WT_mean,mean(dist_with_WT_mean_y[1:105]),sd(dist_with_WT_mean_y[1:105])),
	SSE_SST=sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_y$Res2_mo_y,
	SSE_SST_pvalue=1-pnorm(sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_y$Res2_mo_y,mean(Res2_mo_y[1:105]),sd(Res2_mo_y[1:105]))
)


write.csv(spatial_mutant_res_y,"Tab.S2-spatial_mutant_res_y.csv",row.names=FALSE)



### spatial_mutant_res_z

spatial_mutant_res_z<-data.frame(
	EmbryoID=sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_z$embryoID,
	Gene_name=sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_z$geneName,
	Distance_from_average=sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_z$dist_with_WT_mean,
	Distance_from_average_pvalue=1-pnorm(sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_z$dist_with_WT_mean,mean(dist_with_WT_mean_z[1:105]),sd(dist_with_WT_mean_z[1:105])),
	SSE_SST=sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_z$Res2_mo_z,
	SSE_SST_pvalue=1-pnorm(sigEmbryoInd_dist_with_WT_mean_summary_replicate_mother_defect_z$Res2_mo_z,mean(Res2_mo_z[1:105]),sd(Res2_mo_z[1:105]))
)


write.csv(spatial_mutant_res_z,"Tab.S2-spatial_mutant_res_z.csv",row.names=FALSE)






