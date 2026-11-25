


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
library(ggVennDiagram)



focal_projectFolderPath<-r"(E:\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/figs"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_duration_prcomp.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_raw_nonNA_scaling_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData"
	),projectFolderPath=focal_projectFolderPath)


temporalGene<-read.csv(paste0(focal_projectFolderPath,"/data/canalGene/sigGeneName_temporal.csv"))
spatialGene_x<-read.csv(paste0(focal_projectFolderPath,"/data/canalGene/Tab.S2-spatial_mutant_res_x.csv"))
spatialGene_y<-read.csv(paste0(focal_projectFolderPath,"/data/canalGene/Tab.S2-spatial_mutant_res_y.csv"))
spatialGene_z<-read.csv(paste0(focal_projectFolderPath,"/data/canalGene/Tab.S2-spatial_mutant_res_z.csv"))

temporal_mutant_res<-read.csv(paste0(focal_projectFolderPath,"/data/canalGene/Tab.S1-temporal_mutant_res.csv"))
spatial_mutant_res_x<-read.csv(paste0(focal_projectFolderPath,"/data/canalGene/Tab.S2-spatial_mutant_res_x.csv"))
spatial_mutant_res_y<-read.csv(paste0(focal_projectFolderPath,"/data/canalGene/Tab.S2-spatial_mutant_res_y.csv"))
spatial_mutant_res_z<-read.csv(paste0(focal_projectFolderPath,"/data/canalGene/Tab.S2-spatial_mutant_res_z.csv"))

final_mutant_ind_CCL<-match(temporal_mutant_res$EmbryoID[which(!duplicated(temporal_mutant_res$Gene_name))],rownames(ccl_raw_nonNA_scaling))
final_mutant_ind_x<-match(spatial_mutant_res_x$EmbryoID[which(!duplicated(spatial_mutant_res_x$Gene_name))],rownames(q1_cell350_alignment_xyz))
final_mutant_ind_y<-match(spatial_mutant_res_y$EmbryoID[which(!duplicated(spatial_mutant_res_y$Gene_name))],rownames(q1_cell350_alignment_xyz))
final_mutant_ind_z<-match(spatial_mutant_res_z$EmbryoID[which(!duplicated(spatial_mutant_res_z$Gene_name))],rownames(q1_cell350_alignment_xyz))


###############################################################


### gene overlap

mutantGene<- list(CCL=unique(temporalGene[,2]),x=unique(spatialGene_x[,2]),y=unique(spatialGene_y[,2]),z=unique(spatialGene_z[,2]))

ggVennDiagram(mutantGene,label="count",label_alpha=0,category.names=c("CCL (24)","x (17)","y (10)","z (13)"),edge_size=0.1,lwd=0.1,set_size=3,label_size=8/3)+
scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
theme(legend.position = "none")

ggsave("canalGene_overlap.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


#CCL gene
#num_m<-length(unique(temporalGene[,2]))
#non CCL gene
#num_n<-length(unique(do.call(rbind,strsplit(rownames(ccl_raw_nonNA_scaling),"_"))[,1]))-1-num_m
#coo_x gene
#num_k<-length(unique(spatialGene_x[,2]))
#overlap between CCL genes and coo_x genes 
#num_overlap<-length(intersect(unique(temporalGene[,2]),unique(spatialGene_x[,2])))
#p<e-7
#overlap_pvalue<-binom.test(num_overlap,num_k,num_m/(num_m+num_n))




### GO anotation

write.table(unique(unlist(mutantGene)),paste0(focal_projectFolderPath,"/data/canalGene/sigGene_union_CCL_xyz.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE)

write.table(unique(do.call(rbind,strsplit(rownames(ccl_raw_nonNA_scaling)[1:1928],"_"))[,1])[-1],paste0(focal_projectFolderPath,"/data/canalGene/totalMutants.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE)


## PEA in wormbase 

PEA_res<-read.table(paste0(focal_projectFolderPath,"/data/canalGene/wormbase_PEA_results.txt"),header=T,sep="\t")

dataTemp<-data.frame(Pheno=do.call(rbind,strsplit(PEA_res$Term," WBPhenotype:"))[,1],
	qvalue=-log10(PEA_res$Q.value),
	Count=PEA_res$Observed,
	FC=PEA_res$Enrichment.Fold.Change)

ggplot(dataTemp,aes(x=qvalue,y=Pheno,size=FC))+
geom_point(col="#F8766D",stroke=0.1)+
scale_size(range=c(1,3))+
theme(axis.text.y=element_text(angle=20,hjust=1,vjust=1))+
scale_x_continuous(limits=c(1,3),breaks=c(1,2,3))+
labs(x=expression(paste("-",log[10],"(",italic(q),")")),y="Phenotype",title=" ")


ggsave("wormbase_PEA_results.pdf",width=unit(3.7,"in"),height=unit(1.7,"in"))


## GEA in wormbase 


GEA_res<-read.table(paste0(focal_projectFolderPath,"/data/canalGene/wormbase_GEA_results.txt"),header=T,sep="\t")

dataTemp<-data.frame(Pheno=do.call(rbind,strsplit(GEA_res$Term," GO:"))[,1],
	qvalue=-log10(GEA_res$Q.value),
	Count=GEA_res$Observed,
	FC=GEA_res$Enrichment.Fold.Change)

ggplot(dataTemp,aes(x=qvalue,y=Pheno,size=FC))+
geom_point(col="#F8766D",stroke=0.1)+
scale_size(range=c(1,3.5))+
theme(axis.text.y=element_text(angle=30,hjust=1,vjust=1))+
scale_x_continuous(limits=c(1,3.5),breaks=c(1,2,3))+
labs(x=expression(paste("-",log[10],"(",italic(q),")")),y="GO term",title=" ")


ggsave("wormbase_GEA_results.pdf",width=unit(4,"in"),height=unit(2.7,"in"))



###############################################################
				xyz
###############################################################

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


### x

X_mo<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo


### mother defect

slope_WT_x<-sapply(1:dim(delta_X)[2],function(i){summary(lm(delta_X[1:105,i]~X_mo[1:105,i]))$coefficients[2,1]})

L_delta_L_combine_x<-data.frame(Mo=as.numeric(t(sapply(1:105,function(i){X_mo[i,]*abs(slope_WT_x)}))),Delta=as.numeric(delta_X[1:105,]))

plot_index_x<-final_mutant_ind_x

L_delta_L_combine_case_x<-data.frame(
	Mo=as.numeric(t(sapply(plot_index_x,function(i){X_mo[i,]*abs(slope_WT_x)}))),
	Delta=as.numeric(delta_X[plot_index_x,]),
	Gene=factor(as.character(t(sapply(do.call(rbind,strsplit(rownames(X_mo)[plot_index_x],"_"))[,1],rep,dim(X_mo)[2]))),level=do.call(rbind,strsplit(rownames(X_mo)[plot_index_x],"_"))[,1])
)


###

i=1

L_delta_L_combine_case=data.frame(Mo=as.numeric(X_mo[i,]*abs(slope_WT_x)),Delta=as.numeric(delta_X[i,]))

ggplot(L_delta_L_combine_x,aes(x=Mo,y=Delta))+geom_density2d_filled(contour_var = "ndensity",alpha=0.65,show.legend=FALSE)+scale_fill_brewer(direction=-1,palette =9)+
geom_point(data=L_delta_L_combine_case,aes(x=Mo,y=Delta),col="#F8766D",size=0.3,stroke=0.2,alpha=0.75)+
geom_hline(yintercept=0,lty=2,lwd=0.2)+
geom_vline(xintercept=0,lty=2,lwd=0.2)+
labs(x=expression(paste("-",italic(kX)[M])),y=expression(paste("Δ",italic(X))),title=rownames(X_mo)[i])

ggsave(paste0("spatial_x_mutant_SSE_SST_mother_defect_case_wt_",i,".pdf"),width=unit(1.6,"in"),height=unit(1.7,"in"))


i=3

L_delta_L_combine_case<-L_delta_L_combine_case_x[L_delta_L_combine_case_x$Gene==L_delta_L_combine_case_x$Gene[i],]

ggplot(L_delta_L_combine_x,aes(x=Mo,y=Delta))+geom_density2d_filled(contour_var = "ndensity",alpha=0.65,show.legend=FALSE)+scale_fill_brewer(direction=-1,palette =9)+
geom_point(data=L_delta_L_combine_case,aes(x=Mo,y=Delta),col="#F8766D",size=0.3,stroke=0.2,alpha=0.75)+
geom_hline(yintercept=0,lty=2,lwd=0.2)+
geom_vline(xintercept=0,lty=2,lwd=0.2)+
labs(x=expression(paste("-",italic(kX)[M])),y=expression(paste("Δ",italic(X))),title=L_delta_L_combine_case_x$Gene[i])

ggsave(paste0("spatial_x_mutant_SSE_SST_mother_defect_case_mt_",rownames(X_mo)[plot_index_x[i]],".pdf"),width=unit(1.6,"in"),height=unit(1.7,"in"))


i=4

L_delta_L_combine_case<-L_delta_L_combine_case_x[L_delta_L_combine_case_x$Gene==L_delta_L_combine_case_x$Gene[i],]

ggplot(L_delta_L_combine_x,aes(x=Mo,y=Delta))+geom_density2d_filled(contour_var = "ndensity",alpha=0.65,show.legend=FALSE)+scale_fill_brewer(direction=-1,palette =9)+
geom_point(data=L_delta_L_combine_case,aes(x=Mo,y=Delta),col="#F8766D",size=0.3,stroke=0.2,alpha=0.75)+
geom_hline(yintercept=0,lty=2,lwd=0.2)+
geom_vline(xintercept=0,lty=2,lwd=0.2)+
labs(x=expression(paste("-",italic(kX)[M])),y=expression(paste("Δ",italic(X))),title=L_delta_L_combine_case_x$Gene[i])

ggsave(paste0("spatial_x_mutant_SSE_SST_mother_defect_case_mt_",rownames(X_mo)[plot_index_x[i]],".pdf"),width=unit(1.6,"in"),height=unit(1.7,"in"))




###y

X_mo<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo

slope_WT_y<-sapply(1:dim(delta_X)[2],function(i){summary(lm(delta_X[1:105,i]~X_mo[1:105,i]))$coefficients[2,1]})

L_delta_L_combine_y<-data.frame(Mo=as.numeric(t(sapply(1:105,function(i){X_mo[i,]*abs(slope_WT_y)}))),Delta=as.numeric(delta_X[1:105,]))

plot_index_y<-final_mutant_ind_y

L_delta_L_combine_case_y<-data.frame(
	Mo=as.numeric(t(sapply(plot_index_y,function(i){X_mo[i,]*abs(slope_WT_y)}))),
	Delta=as.numeric(delta_X[plot_index_y,]),
	Gene=factor(as.character(t(sapply(do.call(rbind,strsplit(rownames(X_mo)[plot_index_y],"_"))[,1],rep,dim(X_mo)[2]))),level=do.call(rbind,strsplit(rownames(X_mo)[plot_index_y],"_"))[,1])
)


###

i=1

L_delta_L_combine_case=data.frame(Mo=as.numeric(X_mo[i,]*abs(slope_WT_y)),Delta=as.numeric(delta_X[i,]))

ggplot(L_delta_L_combine_y,aes(x=Mo,y=Delta))+geom_density2d_filled(contour_var = "ndensity",alpha=0.65,show.legend=FALSE)+scale_fill_brewer(direction=-1,palette =9)+
geom_point(data=L_delta_L_combine_case,aes(x=Mo,y=Delta),col="#F8766D",size=0.3,stroke=0.2,alpha=0.75)+
geom_hline(yintercept=0,lty=2,lwd=0.2)+
geom_vline(xintercept=0,lty=2,lwd=0.2)+
labs(x=expression(paste("-",italic(kY)[M])),y=expression(paste("Δ",italic(Y))),title=rownames(X_mo)[i])

ggsave(paste0("spatial_y_mutant_SSE_SST_mother_defect_case_wt_",i,".pdf"),width=unit(1.6,"in"),height=unit(1.7,"in"))


i=1

L_delta_L_combine_case<-L_delta_L_combine_case_y[L_delta_L_combine_case_y$Gene==L_delta_L_combine_case_y$Gene[i],]

ggplot(L_delta_L_combine_y,aes(x=Mo,y=Delta))+geom_density2d_filled(contour_var = "ndensity",alpha=0.65,show.legend=FALSE)+scale_fill_brewer(direction=-1,palette =9)+
geom_point(data=L_delta_L_combine_case,aes(x=Mo,y=Delta),col="#F8766D",size=0.3,stroke=0.2,alpha=0.75)+
geom_hline(yintercept=0,lty=2,lwd=0.2)+
geom_vline(xintercept=0,lty=2,lwd=0.2)+
labs(x=expression(paste("-",italic(kY)[M])),y=expression(paste("Δ",italic(Y))),title=L_delta_L_combine_case_y$Gene[i])

ggsave(paste0("spatial_y_mutant_SSE_SST_mother_defect_case_mt_",rownames(X_mo)[plot_index_y[i]],".pdf"),width=unit(1.6,"in"),height=unit(1.7,"in"))


###z

X_mo<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$mo,feature_names)]
X_da1<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da1,feature_names)]
X_da2<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo

slope_WT_z<-sapply(1:dim(delta_X)[2],function(i){summary(lm(delta_X[1:105,i]~X_mo[1:105,i]))$coefficients[2,1]})

L_delta_L_combine_z<-data.frame(Mo=as.numeric(t(sapply(1:105,function(i){X_mo[i,]*abs(slope_WT_z)}))),Delta=as.numeric(delta_X[1:105,]))

plot_index_z<-final_mutant_ind_z

L_delta_L_combine_case_z<-data.frame(
	Mo=as.numeric(t(sapply(plot_index_z,function(i){X_mo[i,]*abs(slope_WT_z)}))),
	Delta=as.numeric(delta_X[plot_index_z,]),
	Gene=factor(as.character(t(sapply(do.call(rbind,strsplit(rownames(X_mo)[plot_index_z],"_"))[,1],rep,dim(X_mo)[2]))),level=do.call(rbind,strsplit(rownames(X_mo)[plot_index_z],"_"))[,1])
)



###

i=1

L_delta_L_combine_case=data.frame(Mo=as.numeric(X_mo[i,]*abs(slope_WT_z)),Delta=as.numeric(delta_X[i,]))

ggplot(L_delta_L_combine_z,aes(x=Mo,y=Delta))+geom_density2d_filled(contour_var = "ndensity",alpha=0.65,show.legend=FALSE)+scale_fill_brewer(direction=-1,palette =9)+
geom_point(data=L_delta_L_combine_case,aes(x=Mo,y=Delta),col="#F8766D",size=0.3,stroke=0.2,alpha=0.75)+
geom_hline(yintercept=0,lty=2,lwd=0.2)+
geom_vline(xintercept=0,lty=2,lwd=0.2)+
labs(x=expression(paste("-",italic(kZ)[M])),y=expression(paste("Δ",italic(Z))),title=rownames(X_mo)[i])

ggsave(paste0("spatial_z_mutant_SSE_SST_mother_defect_case_wt_",i,".pdf"),width=unit(1.6,"in"),height=unit(1.7,"in"))


i=1

L_delta_L_combine_case<-L_delta_L_combine_case_z[L_delta_L_combine_case_z$Gene==L_delta_L_combine_case_z$Gene[i],]

ggplot(L_delta_L_combine_z,aes(x=Mo,y=Delta))+geom_density2d_filled(contour_var = "ndensity",alpha=0.65,show.legend=FALSE)+scale_fill_brewer(direction=-1,palette =9)+
geom_point(data=L_delta_L_combine_case,aes(x=Mo,y=Delta),col="#F8766D",size=0.3,stroke=0.2,alpha=0.75)+
geom_hline(yintercept=0,lty=2,lwd=0.2)+
geom_vline(xintercept=0,lty=2,lwd=0.2)+
labs(x=expression(paste("-",italic(kZ)[M])),y=expression(paste("Δ",italic(Z))),title=L_delta_L_combine_case_z$Gene[i])

ggsave(paste0("spatial_z_mutant_SSE_SST_mother_defect_case_mt_",rownames(X_mo)[plot_index_z[i]],".pdf"),width=unit(1.6,"in"),height=unit(1.7,"in"))


####################################################################

###################################################################

###CCL 

feature_names<-unlist(strsplit(colnames(ccl_raw_nonNA_scaling_noise),"_ccl"))

temp1<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da1==x)}))
temp3<-unlist(lapply(feature_names,function(x){which(cell_mo_da1_da2$da2==x)}))

mo_sis_pair<-cell_mo_da1_da2[intersect(temp1,temp2),]

#

X_mo<-ccl_raw_nonNA_scaling_noise[1:1928,match(mo_sis_pair$mo,feature_names)]
X_da1<-ccl_raw_nonNA_scaling_noise[1:1928,match(mo_sis_pair$da1,feature_names)]
X_da2<-ccl_raw_nonNA_scaling_noise[1:1928,match(mo_sis_pair$da2,feature_names)]

delta_X<-X_da1-X_mo


#### representative cases

slope_WT_CCL<-sapply(1:dim(delta_X)[2],function(i){summary(lm(delta_X[1:105,i]~X_mo[1:105,i]))$coefficients[2,1]})

L_delta_L_combine_CCL<-data.frame(Mo=as.numeric(t(sapply(1:105,function(i){X_mo[i,]*abs(slope_WT_CCL)}))),Delta=as.numeric(delta_X[1:105,]))

plot_index_CCL<-final_mutant_ind_CCL

L_delta_L_combine_case_CCL<-data.frame(
	Mo=as.numeric(t(sapply(plot_index_CCL,function(i){X_mo[i,]*abs(slope_WT_CCL)}))),
	Delta=as.numeric(delta_X[plot_index_CCL,]),
	Gene=factor(as.character(t(sapply(do.call(rbind,strsplit(rownames(X_mo)[plot_index_CCL],"_"))[,1],rep,dim(X_mo)[2]))),level=do.call(rbind,strsplit(rownames(X_mo)[plot_index_CCL],"_"))[,1])
)



###

i=1

L_delta_L_combine_case=data.frame(Mo=as.numeric(X_mo[i,]*abs(slope_WT_CCL)),Delta=as.numeric(delta_X[i,]))

ggplot(L_delta_L_combine_CCL,aes(x=Mo,y=Delta))+geom_density2d_filled(contour_var = "ndensity",alpha=0.65,show.legend=FALSE)+scale_fill_brewer(direction=-1,palette =9)+
geom_point(data=L_delta_L_combine_case,aes(x=Mo,y=Delta),col="#F8766D",size=0.3,stroke=0.2,alpha=0.75)+
geom_hline(yintercept=0,lty=2,lwd=0.2)+
geom_vline(xintercept=0,lty=2,lwd=0.2)+
labs(x=expression(paste("-",italic(kL)[M])),y=expression(paste("Δ",italic(L))),title=rownames(X_mo)[i])

ggsave(paste0("temporal_mutant_SSE_SST_mother_defect_case_wt_",i,".pdf"),width=unit(1.6,"in"),height=unit(1.7,"in"))


i=3

L_delta_L_combine_case<-L_delta_L_combine_case_CCL[L_delta_L_combine_case_CCL$Gene==L_delta_L_combine_case_CCL$Gene[i],]

ggplot(L_delta_L_combine_CCL,aes(x=Mo,y=Delta))+geom_density2d_filled(contour_var = "ndensity",alpha=0.65,show.legend=FALSE)+scale_fill_brewer(direction=-1,palette =9)+
geom_point(data=L_delta_L_combine_case,aes(x=Mo,y=Delta),col="#F8766D",size=0.3,stroke=0.2,alpha=0.75)+
geom_hline(yintercept=0,lty=2,lwd=0.2)+
geom_vline(xintercept=0,lty=2,lwd=0.2)+
labs(x=expression(paste("-",italic(kL)[M])),y=expression(paste("Δ",italic(L))),title=L_delta_L_combine_case_CCL$Gene[i])


ggsave(paste0("temporal_mutant_SSE_SST_mother_defect_case_mt_",rownames(X_mo)[plot_index_CCL[i]],".pdf"),width=unit(1.6,"in"),height=unit(1.7,"in"))


i=2

L_delta_L_combine_case<-L_delta_L_combine_case_CCL[L_delta_L_combine_case_CCL$Gene==L_delta_L_combine_case_CCL$Gene[i],]

ggplot(L_delta_L_combine_CCL,aes(x=Mo,y=Delta))+geom_density2d_filled(contour_var = "ndensity",alpha=0.65,show.legend=FALSE)+scale_fill_brewer(direction=-1,palette =9)+
geom_point(data=L_delta_L_combine_case,aes(x=Mo,y=Delta),col="#F8766D",size=0.3,stroke=0.2,alpha=0.75)+
geom_hline(yintercept=0,lty=2,lwd=0.2)+
geom_vline(xintercept=0,lty=2,lwd=0.2)+
labs(x=expression(paste("-",italic(kL)[M])),y=expression(paste("Δ",italic(L))),title=L_delta_L_combine_case_CCL$Gene[i])


ggsave(paste0("temporal_mutant_SSE_SST_mother_defect_case_mt_",rownames(X_mo)[plot_index_CCL[i]],".pdf"),width=unit(1.6,"in"),height=unit(1.7,"in"))


























