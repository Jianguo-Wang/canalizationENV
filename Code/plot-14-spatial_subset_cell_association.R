



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

#

get_kTest_p<-function(kTest_result){
	kTest_p<-sapply(kTest_result,function(x){x[,3]})
	colnames(kTest_p)<-names(kTest_result)
	rownames(kTest_p)<-rownames(kTest_result[[1]])
	
	return(kTest_p)
}

#

pMat.adjust<-function(p_now,method="BH",rowCol=0){
	if(rowCol==1){
		return(matrix(t(apply(p_now,rowCol,p.adjust,method=method)),
			nrow=dim(p_now)[1],ncol=dim(p_now)[2],
			dimnames=list(rownames(p_now),colnames(p_now))))			
	}
	if(rowCol==2){
		return(matrix(apply(p_now,rowCol,p.adjust,method=method),
			nrow=dim(p_now)[1],ncol=dim(p_now)[2],
			dimnames=list(rownames(p_now),colnames(p_now))))		
	}
	if(rowCol==0){
		return(matrix(p.adjust(p_now,method=method),
			nrow=dim(p_now)[1],ncol=dim(p_now)[2],
			dimnames=list(rownames(p_now),colnames(p_now))))	
	}

}


#

get_chisq_pvalue<-function(cell_vec,Group){
	
	focal_index<-na.omit(match(cell_vec,Group[,1]))
	if(length(focal_index)==0){
		return(NA)
	}else{
		focal_list<-Group[focal_index,2]
		Expect<-table(Group[na.omit(match(feature_names,Group[,1])),2])
		Observ<-sapply(names(Expect),function(x){length(which(focal_list==x))})
		return(chisq.test(Observ,p=Expect/sum(Expect))$p.value)	
	}
}


#

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

#

get_chisq_foldChange<-function(cell_vec,Group){
	
	if(length(cell_vec)>0&length(na.omit(match(cell_vec,Group[,1])))>0){
		focal_index<-na.omit(match(cell_vec,Group[,1]))
		focal_list<-Group[focal_index,2]
		Expect<-table(Group[match(feature_names,Group[,1]),2])
		Observ<-sapply(names(Expect),function(x){length(which(focal_list==x))})
	
		temp<-sapply(names(Expect),function(x){
			Expect_focalLevel<-c(Expect[names(Expect)==x],Expect[names(Expect)!=x])
			Observ_focalLevel<-c(Observ[names(Expect)==x],Observ[names(Expect)!=x])
			return((Observ_focalLevel[1]/sum(Observ_focalLevel))/(Expect_focalLevel[1]/sum(Expect_focalLevel)))

		})
	
		names(temp)<-names(Expect)		
	}else{
		temp<-rep(NA,length(table(Group[,2])))
	}
	
	
	return(temp)

}


#######################################################################

#######################################################################

library(elitism)
library(ggplot2)
library(pheatmap)


focal_projectFolderPath<-r"(E:\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/figs"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info_short.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/onto_ana.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/onto_phe.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/anatomyInd_list.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/envInd_list.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/phenotypeInd_list.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/anatomy_subset.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/env_subset.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/phenotype_subset.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_env_x.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_env_x_adjP.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_env_y.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_env_y_adjP.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_env_z.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_env_z_adjP.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_anatomy_x.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_anatomy_x_adjP.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_anatomy_y.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_anatomy_y_adjP.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_anatomy_z.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_anatomy_z_adjP.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_phenotype_x.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_phenotype_x_adjP.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_phenotype_y.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_phenotype_y_adjP.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_phenotype_z.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/kTest_phenotype_z_adjP.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/TissueFC_res_env.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/TissueFC_res_anatomy.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/subsetDiff/TissueFC_res_phenotype.RData"
	),projectFolderPath=focal_projectFolderPath)


#

TissueFC_group<-read.csv(paste0(focal_projectFolderPath,"/data/annotation/TissueFC_group.csv"))

TissueFC_group<-TissueFC_group[TissueFC_group[,2]!="Germline"&TissueFC_group[,2]!="Death",]


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

#

X_mo_x<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$mo,feature_names)]
X_da1_x<-q1_cell350_alignment_xyz_noise_x[,match(mo_sis_pair$da1,feature_names)]
delta_X_x<-X_da1_x-X_mo_x

X_mo_y<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$mo,feature_names)]
X_da1_y<-q1_cell350_alignment_xyz_noise_y[,match(mo_sis_pair$da1,feature_names)]
delta_X_y<-X_da1_y-X_mo_y

X_mo_z<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$mo,feature_names)]
X_da1_z<-q1_cell350_alignment_xyz_noise_z[,match(mo_sis_pair$da1,feature_names)]
delta_X_z<-X_da1_z-X_mo_z

#
embryoGene<-tolower(do.call(rbind,strsplit(rownames(q1_cell350_alignment_xyz_noise)[106:1928],"_"))[,1])

embryoGene_uniq<-unique(embryoGene)


#

TissueFC_group[na.omit(match(env_subset[[3]],TissueFC_group[,1])),]
TissueFC_group[na.omit(match(anatomy_subset[[1]],TissueFC_group[,1])),]
TissueFC_group[na.omit(match(anatomy_subset[[3]],TissueFC_group[,1])),]
TissueFC_group[na.omit(match(phenotype_subset[[1]],TissueFC_group[,1])),]
TissueFC_group[na.omit(match(phenotype_subset[[2]],TissueFC_group[,1])),]
TissueFC_group[na.omit(match(phenotype_subset[[3]],TissueFC_group[,1])),]
TissueFC_group[na.omit(match(phenotype_subset[[4]],TissueFC_group[,1])),]

#########################################################################################
						plot - susceptibility
#########################################################################################

# WT 1:105
# RNAi 106:1928
# ENV 1929:2039

# compare between kWT and kRNAi


kTest_kWT_kRNAi_x<-kTest(X_mo_x,delta_X_x,refInd=c(1:105),perturbInd=c(106:1928))
kTest_kWT_kRNAi_y<-kTest(X_mo_y,delta_X_y,refInd=c(1:105),perturbInd=c(106:1928))
kTest_kWT_kRNAi_z<-kTest(X_mo_z,delta_X_z,refInd=c(1:105),perturbInd=c(106:1928))

kTest_kWT_kRNAi_x$adjP<-p.adjust(kTest_kWT_kRNAi_x$k_diff_p,method="BH")
kTest_kWT_kRNAi_y$adjP<-p.adjust(kTest_kWT_kRNAi_y$k_diff_p,method="BH")
kTest_kWT_kRNAi_z$adjP<-p.adjust(kTest_kWT_kRNAi_z$k_diff_p,method="BH")

# compare between kWT and kEnv

kTest_kWT_kEnv_x<-kTest(X_mo_x,delta_X_x,refInd=c(1:105),perturbInd=c(1929:2039))
kTest_kWT_kEnv_y<-kTest(X_mo_y,delta_X_y,refInd=c(1:105),perturbInd=c(1929:2039))
kTest_kWT_kEnv_z<-kTest(X_mo_z,delta_X_z,refInd=c(1:105),perturbInd=c(1929:2039))

kTest_kWT_kEnv_x$adjP<-p.adjust(kTest_kWT_kEnv_x$k_diff_p,method="BH")
kTest_kWT_kEnv_y$adjP<-p.adjust(kTest_kWT_kEnv_y$k_diff_p,method="BH")
kTest_kWT_kEnv_z$adjP<-p.adjust(kTest_kWT_kEnv_z$k_diff_p,method="BH")


### x 


color1<-rep("grey",706)
color1[intersect(which(kTest_kWT_kRNAi_x$k_diff_p<0.05),which(kTest_kWT_kEnv_x$k_diff_p<0.05))]<-"#F8766D"
dataTemp<-data.frame(kWT=kTest_kWT_kEnv_x[,4],kEnv=kTest_kWT_kEnv_x[,5],kRNAi=kTest_kWT_kRNAi_x[,5],kEnv_sus=kTest_kWT_kEnv_x[,1],kRNAi_sus=kTest_kWT_kRNAi_x[,1],signi=color1)

ggplot(dataTemp,aes(x=kWT,y=kEnv))+
geom_point(size=0.2,colour=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,c(1,2)]),breaks=c(-1,-0.5,0))+
scale_y_continuous(limits=range(dataTemp[,c(1,2)]),breaks=c(-1,-0.5,0))+
labs(x=expression(italic(k)[control]),y=expression(italic(k)[Env]),title="")

ggsave("spatial_kContr_kENV_x.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


ggplot(dataTemp,aes(x=kWT,y=kRNAi))+
geom_point(size=0.2,colour=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,c(1,3)]),breaks=c(-1,-0.5,0))+
scale_y_continuous(limits=range(dataTemp[,c(1,3)]),breaks=c(-1,-0.5,0))+
labs(x=expression(italic(k)[control]),y=expression(italic(k)[RNAi]),title="")

ggsave("spatial_kContr_kRNAi_x.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



ggplot(dataTemp,aes(x=kEnv_sus,y=kRNAi_sus))+
geom_point(size=0.2,stroke=0.2,col=dataTemp$signi,alpha=ifelse(dataTemp$signi=="grey",0.5,1))+
geom_abline(intercept=0,slope=1,lwd=0.3)+
geom_hline(yintercept=0,lty=2,lwd=0.3)+
geom_vline(xintercept=0,lty=2,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,c(4,5)]),breaks=c(-0.5,0,0.5))+
scale_y_continuous(limits=range(dataTemp[,c(4,5)]),breaks=c(-0.5,0,0.5))+
labs(x=expression(paste(italic(k)[Env],"-",italic(k)[control])),y=expression(paste(italic(k)[RNAi],"-",italic(k)[control])),title="")

ggsave("spatial_susceptibility_kEnv_kRNAi_x.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


### y

color1<-rep("grey",706)
color1[intersect(which(kTest_kWT_kRNAi_y$k_diff_p<0.05),which(kTest_kWT_kRNAi_y$k_diff_p<0.05))]<-"#F8766D"
dataTemp<-data.frame(kWT=kTest_kWT_kEnv_y[,4],kEnv=kTest_kWT_kEnv_y[,5],kRNAi=kTest_kWT_kRNAi_y[,5],kEnv_sus=kTest_kWT_kEnv_y[,1],kRNAi_sus=kTest_kWT_kRNAi_y[,1],signi=color1)


ggplot(dataTemp,aes(x=kEnv_sus,y=kRNAi_sus))+
geom_point(size=0.2,stroke=0.2,col=dataTemp$signi,alpha=ifelse(dataTemp$signi=="grey",0.5,1))+
geom_abline(intercept=0,slope=1,lwd=0.3)+
geom_hline(yintercept=0,lty=2,lwd=0.3)+
geom_vline(xintercept=0,lty=2,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,c(4,5)]),breaks=c(-0.5,0,0.5))+
scale_y_continuous(limits=range(dataTemp[,c(4,5)]),breaks=c(-0.5,0,0.5))+
labs(x=expression(paste(italic(k)[Env],"-",italic(k)[control])),y=expression(paste(italic(k)[RNAi],"-",italic(k)[control])),title="")

ggsave("spatial_susceptibility_kEnv_kRNAi_y.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



### z

color1<-rep("grey",706)
color1[intersect(which(kTest_kWT_kRNAi_z$k_diff_p<0.05),which(kTest_kWT_kRNAi_z$k_diff_p<0.05))]<-"#F8766D"
dataTemp<-data.frame(kWT=kTest_kWT_kEnv_z[,4],kEnv=kTest_kWT_kEnv_z[,5],kRNAi=kTest_kWT_kRNAi_z[,5],kEnv_sus=kTest_kWT_kEnv_z[,1],kRNAi_sus=kTest_kWT_kRNAi_z[,1],signi=color1)


ggplot(dataTemp,aes(x=kEnv_sus,y=kRNAi_sus))+
geom_point(size=0.2,stroke=0.2,col=dataTemp$signi,alpha=ifelse(dataTemp$signi=="grey",0.5,1))+
geom_abline(intercept=0,slope=1,lwd=0.3)+
geom_hline(yintercept=0,lty=2,lwd=0.3)+
geom_vline(xintercept=0,lty=2,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,c(4,5)]),breaks=c(-0.5,0,0.5))+
scale_y_continuous(limits=range(dataTemp[,c(4,5)]),breaks=c(-0.5,0,0.5))+
labs(x=expression(paste(italic(k)[Env],"-",italic(k)[control])),y=expression(paste(italic(k)[RNAi],"-",italic(k)[control])),title="")

ggsave("spatial_susceptibility_kEnv_kRNAi_z.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



### enrichment of susceptibility cells 

temp_x<-intersect(which(kTest_kWT_kRNAi_x$k_diff_p<0.05),which(kTest_kWT_kEnv_x$k_diff_p<0.05))
temp_y<-intersect(which(kTest_kWT_kRNAi_y$k_diff_p<0.05),which(kTest_kWT_kEnv_y$k_diff_p<0.05))
temp_z<-intersect(which(kTest_kWT_kRNAi_z$k_diff_p<0.05),which(kTest_kWT_kEnv_z$k_diff_p<0.05))

length(which(kTest_kWT_kRNAi_x$k_diff[temp_x]>0&kTest_kWT_kEnv_x$k_diff[temp_x]>0))/length(temp_x)
length(which(kTest_kWT_kRNAi_y$k_diff[temp_x]>0&kTest_kWT_kEnv_y$k_diff[temp_x]>0))/length(temp_y)
length(which(kTest_kWT_kRNAi_z$k_diff[temp_x]>0&kTest_kWT_kEnv_z$k_diff[temp_x]>0))/length(temp_z)

dataTemp<-data.frame(kdiff_RNAi=c(kTest_kWT_kRNAi_x[temp_x,1],kTest_kWT_kRNAi_y[temp_y,1],kTest_kWT_kRNAi_z[temp_z,1]),
	kdiff_Env=c(kTest_kWT_kEnv_x[temp_x,1],kTest_kWT_kEnv_y[temp_y,1],kTest_kWT_kEnv_z[temp_z,1]),
	coo=c(rep("x",length(temp_x)),rep("y",length(temp_y)),rep("z",length(temp_z))))

ggplot(dataTemp,aes(x=kdiff_Env,y=kdiff_RNAi,col=coo))+
geom_point(size=0.2,stroke=0.2)+
geom_hline(yintercept=0,lty=2,lwd=0.3)+
geom_vline(xintercept=0,lty=2,lwd=0.3)+
scale_x_continuous(breaks=c(0,1.0))+
scale_y_continuous(breaks=c(0,0.5))+
labs(x=expression(paste(italic(k)[Env],"-",italic(k)[control])),y=expression(paste(italic(k)[RNAi],"-",italic(k)[control])),title="significant cells")

ggsave("spatial_susceptibility_kEnv_kRNAi_significant_cells.pdf",width=unit(1.75,"in"),height=unit(1.3,"in"))



#ignore
sigCells<-unique(do.call(rbind,strsplit(c(rownames(kTest_kWT_kRNAi_x)[temp_x],rownames(kTest_kWT_kRNAi_y)[temp_y],rownames(kTest_kWT_kRNAi_z)[temp_z]),"_"))[,3])

sigCells_levelP<-get_chisq_levelP(sigCells,TissueFC_group)
sigCells_foldChange<-get_chisq_foldChange(sigCells,TissueFC_group)


dataTemp<-data.frame(p=sigCells_levelP,tissue=names(sigCells_levelP),FC=sigCells_foldChange)

ggplot(dataTemp,aes(y=tissue,x=-log10(p)))+
geom_bar(position="stack",stat="identity",fill=colors()[182],width=0.7)+
geom_vline(xintercept=1.3,lty=2,lwd=0.3,col="red")+
scale_x_continuous(breaks=c(0,2))+
labs(y="",x=expression(paste("-",log[10],"(",italic(p),")")),title="")+
theme(axis.text.y=element_text(angle=0,vjust=0.5,hjust=1))

ggsave("spatial_susceptibility_kEnv_kRNAi_significant_cells_tissue_enrich.pdf",width=unit(1.7,"in"),height=unit(1.7,"in"))

#
kTest_kWT_kRNAi_xyz_signif<-rbind(kTest_kWT_kRNAi_x[temp_x,],kTest_kWT_kRNAi_y[temp_y,],kTest_kWT_kRNAi_z[temp_z,])
kTest_kWT_kRNAi_xyz_signif_pos<-kTest_kWT_kRNAi_xyz_signif[kTest_kWT_kRNAi_xyz_signif[,1]>0,]
kTest_kWT_kRNAi_xyz_signif_neg<-kTest_kWT_kRNAi_xyz_signif[kTest_kWT_kRNAi_xyz_signif[,1]<0,]

sigCells_pos<-do.call(rbind,strsplit(rownames(kTest_kWT_kRNAi_xyz_signif_pos),"_"))[,3]
sigCells_pos_levelP<-get_chisq_levelP(sigCells_pos,TissueFC_group)
sigCells_pos_foldChange<-get_chisq_foldChange(sigCells_pos,TissueFC_group)

sigCells_neg<-do.call(rbind,strsplit(rownames(kTest_kWT_kRNAi_xyz_signif_neg),"_"))[,3]
sigCells_neg_levelP<-get_chisq_levelP(sigCells_neg,TissueFC_group)
sigCells_neg_foldChange<-get_chisq_foldChange(sigCells_neg,TissueFC_group)

dataTemp<-data.frame(p_pos=sigCells_pos_levelP,p_neg=sigCells_neg_levelP)

ggplot(dataTemp,aes(x=-log10(p_pos),y=-log10(p_neg)))+
geom_point(col=colors()[182],size=0.5)+
geom_vline(xintercept=1.3,lty=2,lwd=0.3,col="grey")+
geom_hline(yintercept=1.3,lty=2,lwd=0.3,col="grey")+
labs(x=expression(paste("-",log[10],"(",italic(p),") for increased-",italic(k)," cells")),
	y=expression(paste("-",log[10],"(",italic(p),") for decreased-",italic(k)," cells")),title="Tissue enrichment")


ggsave("spatial_susceptibility_kEnv_kRNAi_significant_cells_tissue_enrich_bidirection.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))



### y

dataTemp<-data.frame(kWT=kTest_kWT_kEnv_y[,4],kEnv=kTest_kWT_kEnv_y[,5],kRNAi=kTest_kWT_kRNAi_y[,5],kEnv_sus=kTest_kWT_kEnv_y[,1],kRNAi_sus=kTest_kWT_kRNAi_y[,1])

ggplot(dataTemp,aes(x=kWT,y=kEnv))+
geom_point(size=0.2,colour=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,c(1,2)]),breaks=c(-1,-0.5,0))+
scale_y_continuous(limits=range(dataTemp[,c(1,2)]),breaks=c(-1,-0.5,0))+
labs(x=expression(italic(k)[control]),y=expression(italic(k)[Env]),title="y")

ggsave("spatial_kContr_kENV_y.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


ggplot(dataTemp,aes(x=kWT,y=kRNAi))+
geom_point(size=0.2,colour=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,c(1,3)]),breaks=c(-1,-0.5,0))+
scale_y_continuous(limits=range(dataTemp[,c(1,3)]),breaks=c(-1,-0.5,0))+
labs(x=expression(italic(k)[control]),y=expression(italic(k)[RNAi]),title="y")

ggsave("spatial_kContr_kRNAi_y.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


ggplot(dataTemp,aes(x=kEnv_sus,y=kRNAi_sus))+
geom_point(size=0.2,stroke=0.2,colour=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
geom_hline(yintercept=0,lty=2,lwd=0.3)+
geom_vline(xintercept=0,lty=2,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,c(4,5)]),breaks=c(-0.5,0,0.5))+
scale_y_continuous(limits=range(dataTemp[,c(4,5)]),breaks=c(-0.5,0,0.5))+
labs(x=expression(paste(italic(k)[Env],"-",italic(k)[control])),y=expression(paste(italic(k)[RNAi],"-",italic(k)[control])),title="y")

ggsave("spatial_susceptibility_kEnv_kRNAi_y.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


### z

dataTemp<-data.frame(kWT=kTest_kWT_kEnv_z[,4],kEnv=kTest_kWT_kEnv_z[,5],kRNAi=kTest_kWT_kRNAi_z[,5],kEnv_sus=kTest_kWT_kEnv_z[,1],kRNAi_sus=kTest_kWT_kRNAi_z[,1])

ggplot(dataTemp,aes(x=kWT,y=kEnv))+
geom_point(size=0.2,colour=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,c(1,2)]),breaks=c(-1,-0.5,0))+
scale_y_continuous(limits=range(dataTemp[,c(1,2)]),breaks=c(-1,-0.5,0))+
labs(x=expression(italic(k)[control]),y=expression(italic(k)[Env]),title="z")

ggsave("spatial_kContr_kENV_z.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


ggplot(dataTemp,aes(x=kWT,y=kRNAi))+
geom_point(size=0.2,colour=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,c(1,3)]),breaks=c(-1,-0.5,0))+
scale_y_continuous(limits=range(dataTemp[,c(1,3)]),breaks=c(-1,-0.5,0))+
labs(x=expression(italic(k)[control]),y=expression(italic(k)[RNAi]),title="z")

ggsave("spatial_kContr_kRNAi_z.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


ggplot(dataTemp,aes(x=kEnv_sus,y=kRNAi_sus))+
geom_point(size=0.2,stroke=0.2,colour=colors()[182])+
geom_abline(intercept=0,slope=1,lwd=0.3)+
geom_hline(yintercept=0,lty=2,lwd=0.3)+
geom_vline(xintercept=0,lty=2,lwd=0.3)+
scale_x_continuous(limits=range(dataTemp[,c(4,5)]),breaks=c(-0.5,0,0.5))+
scale_y_continuous(limits=range(dataTemp[,c(4,5)]),breaks=c(-0.5,0,0.5))+
labs(x=expression(paste(italic(k)[Env],"-",italic(k)[control])),y=expression(paste(italic(k)[RNAi],"-",italic(k)[control])),title="z")

ggsave("spatial_susceptibility_kEnv_kRNAi_z.pdf",width=unit(1.2,"in"),height=unit(1.3,"in"))


###

dataTemp<-data.frame(Epra_x=X_mo_x[,605],Epraa_x=X_da1_x[,605],type=factor(c(rep("control",105),rep("RNAi",1823),rep("Env",111)),level=c("control","Env","RNAi")))

lm(delta_X_x[,605]~X_mo_x[,605])
lm(X_da1_x[1:105,605]~X_mo_x[1:105,605])
lm(X_da1_x[106:1928,605]~X_mo_x[106:1928,605])
lm(X_da1_x[1929:2039,605]~X_mo_x[1929:2039,605])

ggplot(dataTemp,aes(x=Epra_x,y=Epraa_x))+
geom_point(size=0.2,stroke=0.2)+
geom_smooth(method="lm",lwd=0.3)+
scale_x_continuous(limits=c(-25,25),breaks=c(-20,0,20))+
scale_y_continuous(limits=c(-25,25),breaks=c(-20,0,20))+
facet_wrap(vars(type))+
labs(x="Mother noise (Epra_x)",y="Daughter noise (Epraa_x)")

ggsave("spatial_susceptibility_kEnv_kRNAi_x_case_Epra_Epraa.pdf",width=unit(2.7,"in"),height=unit(1.3,"in"))


###

dataTemp<-cbind(kTest_env_x_adjP,kTest_env_y_adjP,kTest_env_z_adjP)
colnames(dataTemp)<-c(paste0(colnames(kTest_env_x_adjP),"(x)"),paste0(colnames(kTest_env_x_adjP),"(y)"),paste0(colnames(kTest_env_x_adjP),"(z)"))
dataTemp[which(dataTemp==0)]<-min(setdiff(as.numeric(dataTemp),0))

dataTemp2<-dataTemp[,c(sapply(1:8,function(x){x+c(0,8,16)}))]

pdf("env_cell_association.pdf",height=2.3,width=3.5)

pheatmap(t(sqrt(-log10(dataTemp2))),border_color=NA,labels_col=FALSE,cluster_cols=F,cluster_rows=F,fontsize_row=5,legend_breaks=c(1,2,3))

dev.off()


#


dataTemp<-cbind(kTest_env_x_adjP,kTest_env_y_adjP,kTest_env_z_adjP)
colnames(dataTemp)<-c(paste0(colnames(kTest_env_x_adjP),"(x)"),paste0(colnames(kTest_env_x_adjP),"(y)"),paste0(colnames(kTest_env_x_adjP),"(z)"))
dataTemp[which(dataTemp==0)]<-min(setdiff(as.numeric(dataTemp),0))

temp<-apply(dataTemp,2,function(x){which(x<0.05)})

dataTemp3<-data.frame(Count=sapply(1:8,function(i){length(unique(c(temp[[i]],temp[[i+8]],temp[[i+16]])))}),
	Env=factor(c("Temperature","Heat shock","Bacterial food","Stearic acid","Glucose","NaCl","FeCl3","Paraquat"),
	level=c("Paraquat","FeCl3","NaCl","Glucose","Stearic acid","Bacterial food","Heat shock","Temperature")))

ggplot(dataTemp3,aes(x=Env,y=Count))+
geom_bar(position="stack",stat="identity",col=colors()[182],width=0.6)+
scale_y_continuous(limits=c(0,90),breaks=c(0,45,90))+
labs(x="Environment",y="Number of Cells",title="")+
theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))

ggsave("env_cell_association_summary.pdf",width=unit(1.6,"in"),height=unit(1.7,"in"))


###

pdf("env_tissue_association.pdf",height=2.3,width=2.1)

dataTemp3<--log10(TissueFC_res_env[,4:8])
rownames(dataTemp3)<-TissueFC_res_env[,"description"]
pheatmap(dataTemp3,angle_col=90,border_color=NA,cellwidth=8,cellheight=8,
fontsize_row=7,fontsize_col=7,main="Tissue type",cluster_rows=F,cluster_cols=F,legend_breaks=c(0:3))


dev.off()


###

col1<-rep("grey",706)
col1[which(kTest_env_y_adjP[,3]<0.05)]<-"#F8766D"

size1<-rep(0.2,706)
size1[which(kTest_env_y_adjP[,3]<0.05)]<-0.5

alpha1<-rep(0.75,706)
alpha1[which(kTest_env_y_adjP[,3]<0.05)]<-1


ggplot(kTest_env_y[[3]],aes(x=k_ref,y=k_perturb))+
geom_point(size=size1,stroke=0.2,colour=col1,alpha=alpha1)+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=c(-2,1),breaks=c(-1,0))+
scale_y_continuous(limits=c(-2,1),breaks=c(-1,0))+
labs(x=expression(italic(k)[control]),y=expression(italic(k)[Bacterial]),title="Bacterial food\n(y-coordinate)")

ggsave("spatial_kContr_kBacteria_y.pdf",width=unit(1.8,"in"),height=unit(1.9,"in"))


###


dataTemp<-data.frame(Mother=c(c(X_mo_y[envInd_list[[3]],c(585,586,593,594)]),c(X_mo_y[1:105,c(585,586,593,594)])),
	Daughter=c(c(X_da1_y[envInd_list[[3]],c(585,586,593,594)]),c(X_da1_y[1:105,c(585,586,593,594)])),
	Embryo_group=factor(c(rep("Bacterial",length(envInd_list[[3]])*4),rep("control",105*4)),level=c("control","Bacterial")),
	Cell_pair=c(c(sapply(do.call(rbind,strsplit(colnames(X_da1_y)[c(585,586,593,594)],"_"))[,1],rep,length(envInd_list[[3]]))),
		c(sapply(do.call(rbind,strsplit(colnames(X_da1_y)[c(585,586,593,594)],"_"))[,1],rep,105))))

ggplot(dataTemp,aes(x=Mother,y=Daughter))+geom_point(colour=colors()[182],size=0.3,stroke=0.2)+geom_smooth(method="lm",lwd=0.5)+
facet_grid(vars(Embryo_group),vars(Cell_pair))+
scale_y_continuous(breaks=c(-20,0,20))+
labs(x="Mother noise (y)",y="Daughter noise (y)",title="")+
theme(strip.text=element_text(size=7))

ggsave("kTest_Bacterial_food_kdiff_y.pdf",width=unit(2.87,"in"),height=unit(1.9,"in"))


###### anatomy 
###### dorsal nerve cord 

###

col1<-rep("grey",706)
col1[which(kTest_anatomy_z_adjP[,1]<0.05)]<-"#F8766D"

size1<-rep(0.2,706)
size1[which(kTest_anatomy_z_adjP[,1]<0.05)]<-0.5

alpha1<-rep(0.75,706)
alpha1[which(kTest_anatomy_z_adjP[,1]<0.05)]<-1


ggplot(kTest_anatomy_z[[1]],aes(x=k_ref,y=k_perturb))+
geom_point(size=size1,stroke=0.2,colour=col1,alpha=alpha1)+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(kTest_anatomy_z[[1]][,4:5]),breaks=c(-1,0))+
scale_y_continuous(limits=range(kTest_anatomy_z[[1]][,4:5]),breaks=c(-1,0))+
labs(x=expression(italic(k)[control]),y=expression(italic(k)[dorsal]),title="Dorsal nerve cord\n(z-coordinate)")

ggsave("spatial_kContr_kDorsalNerveCord_z.pdf",width=unit(1.8,"in"),height=unit(1.9,"in"))


###


dataTemp<-data.frame(Mother=c(c(X_mo_z[anatomyInd_list[[1]],c(321,322,351,352)]),c(X_mo_z[1:105,c(321,322,351,352)])),
	Daughter=c(c(X_da1_z[anatomyInd_list[[1]],c(321,322,351,352)]),c(X_da1_z[1:105,c(321,322,351,352)])),
	Embryo_group=factor(c(rep("Dorsal",length(anatomyInd_list[[1]])*4),rep("control",105*4)),level=c("control","Dorsal")),
	Cell_pair=c(c(sapply(do.call(rbind,strsplit(colnames(X_da1_z)[c(321,322,351,352)],"_"))[,1],rep,length(anatomyInd_list[[1]]))),
		c(sapply(do.call(rbind,strsplit(colnames(X_da1_z)[c(321,322,351,352)],"_"))[,1],rep,105))))

ggplot(dataTemp,aes(x=Mother,y=Daughter))+
geom_point(colour=colors()[182],size=0.3,stroke=0.2)+
geom_smooth(method="lm",lwd=0.5)+
scale_x_continuous(breaks=c(-40,0))+
scale_y_continuous(breaks=c(-40,0))+
facet_grid(vars(Embryo_group),vars(Cell_pair))+
labs(x="Mother noise (z)",y="Daughter noise",title="")+
theme(strip.text=element_text(size=7))

ggsave("kTest_DorsalNerveCord_kdiff.pdf",width=unit(2.87,"in"),height=unit(1.9,"in"))

###### neurons controlling muscle

###

col1<-rep("grey",706)
col1[which(kTest_anatomy_y_adjP[,3]<0.05)]<-"#F8766D"

size1<-rep(0.2,706)
size1[which(kTest_anatomy_y_adjP[,3]<0.05)]<-0.5

alpha1<-rep(0.75,706)
alpha1[which(kTest_anatomy_y_adjP[,3]<0.05)]<-1


ggplot(kTest_anatomy_y[[3]],aes(x=k_ref,y=k_perturb))+
geom_point(size=size1,stroke=0.2,colour=col1,alpha=alpha1)+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(kTest_anatomy_y[[3]][,4:5]),breaks=c(-1,0))+
scale_y_continuous(limits=range(kTest_anatomy_y[[3]][,4:5]),breaks=c(-1,0))+
labs(x=expression(italic(k)[control]),y=expression(italic(k)[RNAi]),title="Neurons controlling muscle")

ggsave("spatial_kContr_NeuronsControllingMuscle_x.pdf",width=unit(1.8,"in"),height=unit(1.9,"in"))


###


dataTemp<-data.frame(Mother=c(c(X_mo_y[anatomyInd_list[[3]],c(525:528,531:533,561)]),c(X_mo_y[1:105,c(525:528,531:533,561)])),
	Daughter=c(c(X_da1_y[anatomyInd_list[[3]],c(525:528,531:533,561)]),c(X_da1_y[1:105,c(525:528,531:533,561)])),
	Embryo_group=factor(c(rep("N2M",length(anatomyInd_list[[3]])*8),rep("control",105*8)),level=c("control","N2M")),
	Cell_pair=c(c(sapply(colnames(X_da1_y)[c(525:528,531:533,561)],rep,length(anatomyInd_list[[3]]))),
		c(sapply(colnames(X_da1_y)[c(525:528,531:533,561)],rep,105))))

ggplot(dataTemp,aes(x=Mother,y=Daughter))+
geom_point(colour=colors()[182],size=0.3,stroke=0.2)+
geom_smooth(method="lm",lwd=0.5)+
scale_x_continuous(breaks=c(-40,0))+
scale_y_continuous(breaks=c(-40,0))+
facet_grid(vars(Embryo_group),vars(Cell_pair))+
labs(x="Mother noise (y)",y="Daughter noise",title="")+
theme(strip.text=element_text(size=7))

ggsave("kTest_Neuromuscular_kdiff.pdf",width=unit(5.5,"in"),height=unit(1.9,"in"))



###### phenotype 
###### organism oxidative stress resistant 

###

col1<-rep("grey",706)
col1[which(kTest_phenotype_y_adjP[,1]<0.05)]<-"#F8766D"

size1<-rep(0.2,706)
size1[which(kTest_phenotype_y_adjP[,1]<0.05)]<-0.5

alpha1<-rep(0.75,706)
alpha1[which(kTest_phenotype_y_adjP[,1]<0.05)]<-1


ggplot(kTest_phenotype_y[[1]],aes(x=k_ref,y=k_perturb))+
geom_point(size=size1,stroke=0.2,colour=col1,alpha=alpha1)+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(kTest_phenotype_y[[1]][,4:5]),breaks=c(-1,0))+
scale_y_continuous(limits=range(kTest_phenotype_y[[1]][,4:5]),breaks=c(-1,0))+
labs(x=expression(italic(k)[control]),y=expression(italic(k)[RNAi]),title="Oxidative stress resistant")

ggsave("spatial_kContr_oxidativeStreeResistant_x.pdf",width=unit(1.8,"in"),height=unit(1.9,"in"))


###



dataTemp<-data.frame(Mother=c(c(X_mo_y[phenotypeInd_list[[1]],c(585:588)]),c(X_mo_y[1:105,c(585:588)])),
	Daughter=c(c(X_da1_y[phenotypeInd_list[[1]],c(585:588)]),c(X_da1_y[1:105,c(585:588)])),
	Embryo_group=factor(c(rep("Oxidative",length(phenotypeInd_list[[1]])*4),rep("control",105*4)),level=c("control","Oxidative")),
	Cell_pair=c(c(sapply(do.call(rbind,strsplit(colnames(X_da1_y)[c(585:588)],"_"))[,1],rep,length(phenotypeInd_list[[1]]))),
		c(sapply(do.call(rbind,strsplit(colnames(X_da1_y)[c(585:588)],"_"))[,1],rep,105))))

ggplot(dataTemp,aes(x=Mother,y=Daughter))+
geom_point(colour=colors()[182],size=0.3,stroke=0.2)+
geom_smooth(method="lm",lwd=0.5)+
facet_grid(vars(Embryo_group),vars(Cell_pair))+
scale_x_continuous(breaks=c(-30,0))+
scale_y_continuous(breaks=c(-30,0))+
labs(x="Mother noise (y)",y="Daughter noise",title="")+
theme(strip.text=element_text(size=7))

ggsave("kTest_oxidativeStreeResistant_kdiff.pdf",width=unit(2.87,"in"),height=unit(1.9,"in"))


###### aldicarb hypersensitive 

###

col1<-rep("grey",706)
col1[which(kTest_phenotype_y_adjP[,2]<0.05)]<-"#F8766D"

size1<-rep(0.2,706)
size1[which(kTest_phenotype_y_adjP[,1]<0.05)]<-0.5

alpha1<-rep(0.75,706)
alpha1[which(kTest_phenotype_y_adjP[,1]<0.05)]<-1


ggplot(kTest_phenotype_y[[1]],aes(x=k_ref,y=k_perturb))+
geom_point(size=size1,stroke=0.2,colour=col1,alpha=alpha1)+
geom_abline(intercept=0,slope=1,lwd=0.3)+
scale_x_continuous(limits=range(kTest_phenotype_y[[1]][,4:5]),breaks=c(-1,0))+
scale_y_continuous(limits=range(kTest_phenotype_y[[1]][,4:5]),breaks=c(-1,0))+
labs(x=expression(italic(k)[control]),y=expression(italic(k)[RNAi]),title="oxidative stress resistant")

ggsave("spatial_kContr_oxidativeStreeResistant_x.pdf",width=unit(1.8,"in"),height=unit(1.9,"in"))


### error



dataTemp<-data.frame(Mother=c(c(X_mo_y[phenotypeInd_list[[1]],c(585:588)]),c(X_mo_y[1:105,c(585:588)])),
	Daughter=c(c(X_da1_y[phenotypeInd_list[[1]],c(585:588)]),c(X_da1_y[1:105,c(585:588)])),
	Embryo_group=factor(c(rep("Oxidative",length(phenotypeInd_list[[1]])*4),rep("control",105*4)),level=c("control","Oxidative")),
	Cell_pair=c(c(sapply(colnames(X_da1_y)[c(585:588)],rep,length(phenotypeInd_list[[1]]))),
		c(sapply(colnames(X_da1_y)[c(585:588)],rep,105))))

ggplot(dataTemp,aes(x=Mother,y=Daughter))+
geom_point(colour=colors()[182],size=0.3,stroke=0.2)+
geom_smooth(method="lm",lwd=0.5)+
facet_grid(vars(Embryo_group),vars(Cell_pair))+
labs(x="Mother noise (y)",y="Daughter noise",title="")+
theme(strip.text=element_text(size=7))

ggsave("kTest_oxidativeStreeResistant_kdiff.pdf",width=unit(2.87,"in"),height=unit(1.9,"in"))



##########################################################################################################

### spatial distribution 
color_label=c("Insignificant","other","Intestine","Muscle","Neuronal system","Pharynx","Skin")
color_name=c("lightGrey",colors()[192],"#F8766D","#A3A500","#00B0F6","#79DEEC","#E76BF3")

color_set<-function(tiss_vec){
	temp<-rep(colors()[192],length(tiss_vec))
	temp[tiss_vec=="Intestine"]<-"#F8766D"
	temp[tiss_vec=="Muscle"]<-"#A3A500"
	temp[tiss_vec=="Neuronal system"]<-"#79DEEC"
	temp[tiss_vec=="Pharynx"]<-"#00B0F6"
	temp[tiss_vec=="Skin"]<-"#E76BF3"
	
	return(temp)
}

coo_x<-apply(q1_cell350_alignment_xyz[1:105,1:714],2,mean)
coo_y<-apply(q1_cell350_alignment_xyz[1:105,c(1:714)+714],2,mean)
coo_z<-apply(q1_cell350_alignment_xyz[1:105,c(1:714)+714*2],2,mean)

### bacterial food 

focalCell<-env_subset[[3]]
focalCellIndex<-match(focalCell,feature_names)
focalCellTissue<-TissueFC_group[match(focalCell,TissueFC_group[,1]),2]
focalCellTissue[which(is.na(focalCellTissue))]<-"other"
focalCellColor<-color_set(focalCellTissue)

tissueType=rep("Insignificant",714)
tissueType[focalCellIndex]<-focalCellTissue

cellSize1=rep(0.1,714)
cellSize1[focalCellIndex]<-1

alpha1=rep(0.75,714)
alpha1[focalCellIndex]<-1

dataTemp<-data.frame(x=coo_x,y=coo_y,z=coo_z,tissue=factor(tissueType,level=color_label),size1=ifelse(cellSize1==1,1,0.1))

# x-y 

ggplot(dataTemp,aes(x=x,y=y,col=tissue))+
geom_point(size=cellSize1,alpha=alpha1)+
scale_color_manual(values=color_name)+
labs(x="x-coordinate",y="y-coordinate",title="Bacterial food")

ggsave("spatial_distribution_bacterial_food_xy.pdf",width=unit(3.5,"in"),height=unit(1.55,"in"))


## non-random distribution 
dataTemp_signif<-dataTemp[dataTemp$tissue!="Insignificant",]
dataTemp_insignif<-dataTemp[dataTemp$tissue=="Insignificant",]

between_signif<-matrix(nrow=dim(dataTemp_signif)[1],ncol=dim(dataTemp_signif)[1])
between_sig_insig<-matrix(nrow=dim(dataTemp_signif)[1],ncol=dim(dataTemp_insignif)[1])

for(i in 1:(dim(dataTemp_signif)[1]-1)){
	for(j in (i+1):dim(dataTemp_signif)[1]){
		between_signif[i,j]<-sqrt(sum((dataTemp_signif[i,1]-dataTemp_signif[j,1])^2+(dataTemp_signif[i,2]-dataTemp_signif[j,2])^2+(dataTemp_signif[i,3]-dataTemp_signif[j,3])^2))	
	}

}

for(i in 1:dim(dataTemp_signif)[1]){
	for(j in 1:dim(dataTemp_insignif)[1]){
		between_sig_insig[i,j]<-sqrt(sum((dataTemp_signif[i,1]-dataTemp_insignif[j,1])^2+(dataTemp_signif[i,2]-dataTemp_insignif[j,2])^2+(dataTemp_signif[i,3]-dataTemp_insignif[j,3])^2))	
	}

}

par(mfrow=c(1,2))
hist(between_signif[upper.tri(between_signif)])
hist(between_sig_insig)

temp<-data.frame(distance=c(c(between_signif[upper.tri(between_signif)]),c(between_sig_insig)),
type=c(rep("within",length(c(between_signif[upper.tri(between_signif)]))),rep("between",length(c(between_sig_insig)))))

ggplot(temp,aes(x=distance,col=type))+geom_density()

t.test(c(between_signif[upper.tri(between_signif)]),c(between_sig_insig),alternative="less")
wilcox.test(c(between_signif[upper.tri(between_signif)]),c(between_sig_insig),alternative="less")



# x-z 

ggplot(dataTemp,aes(x=x,y=z,col=tissue))+
geom_point(size=cellSize1,alpha=alpha1)+
scale_color_manual(values=color_name)+
labs(x="x-coordinate",y="z-coordinate",title="Bacterial food")

ggsave("spatial_distribution_bacterial_food_xz.pdf",width=unit(3.5,"in"),height=unit(1.25,"in"))

# y-z 

ggplot(dataTemp,aes(x=y,y=z,col=tissue))+
geom_point(size=cellSize1,alpha=alpha1)+
scale_color_manual(values=color_name)+
labs(x="y-coordinate",y="z-coordinate",title="Bacterial food")

ggsave("spatial_distribution_bacterial_food_yz.pdf",width=unit(2.7,"in"),height=unit(1.25,"in"))



### dorsal nerve cord

focalCell<-anatomy_subset[[1]]
focalCellIndex<-match(focalCell,feature_names)
focalCellTissue<-TissueFC_group[match(focalCell,TissueFC_group[,1]),2]
focalCellTissue[which(is.na(focalCellTissue))]<-"other"
focalCellColor<-color_set(focalCellTissue)

tissueType=rep("Insignificant",714)
tissueType[focalCellIndex]<-focalCellTissue

cellSize2=rep(0.1,714)
cellSize2[focalCellIndex]<-1

alpha1=rep(0.75,714)
alpha1[focalCellIndex]<-1

dataTemp2<-data.frame(x=coo_x,y=coo_y,z=coo_z,tissue=factor(tissueType,level=intersect(color_label,tissueType)),size1=ifelse(cellSize2==1,1,0.1))

# dorsal nerve cord (x-y)

ggplot(dataTemp2,aes(x=x,y=y,col=tissue))+
geom_point(size=cellSize2,alpha=alpha1)+
scale_color_manual(values=color_name[match(intersect(color_label,tissueType),color_label)])+
labs(x="x-coordinate",y="y-coordinate",title="Dorsal nerve cord")

ggsave("spatial_distribution_dorsal_nerve_cord_xy.pdf",width=unit(3.5,"in"),height=unit(1.55,"in"))


# dorsal nerve cord (x-z)

ggplot(dataTemp2,aes(x=x,y=z,col=tissue))+
geom_point(size=cellSize2,alpha=alpha1)+
scale_color_manual(values=color_name[match(intersect(color_label,tissueType),color_label)])+
labs(x="x-coordinate",y="z-coordinate",title="Dorsal nerve cord")

ggsave("spatial_distribution_dorsal_nerve_cord_xz.pdf",width=unit(3.5,"in"),height=unit(1.25,"in"))


# dorsal nerve cord (y-z)

ggplot(dataTemp2,aes(x=y,y=z,col=tissue))+
geom_point(size=cellSize2,alpha=alpha1)+
scale_color_manual(values=color_name[match(intersect(color_label,tissueType),color_label)])+
labs(x="y-coordinate",y="z-coordinate",title="Dorsal nerve cord")

ggsave("spatial_distribution_dorsal_nerve_cord_yz.pdf",width=unit(2.7,"in"),height=unit(1.25,"in"))


### comb 

dataTemp3<-rbind(dataTemp,dataTemp2)
dataTemp3$trait=c(rep("bacterial food",714),rep("Dorsal nerve cord",714))
dataTemp3$tissue<-factor(dataTemp3$tissue,level=color_label)

ggplot(dataTemp3,aes(x=x,y=y,col=tissue))+
geom_point(size=c(cellSize1,cellSize2))+
scale_color_manual(values=color_name)+
labs(x="x-coordinate",y="y-coordinate",title="")+
facet_wrap(vars(trait))

ggsave("spatial_distribution_bacterial_food_dorsal_nerve_cord.pdf",width=unit(6.4,"in"),height=unit(2.1,"in"))



