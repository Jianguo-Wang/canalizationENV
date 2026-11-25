


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

extract_onto<-function(focal_onto){
	onto_temp<-data.frame(id=do.call(rbind,strsplit(focal_onto[which(focal_onto=="[Term]")+1],"id: "))[,2],
		description=do.call(rbind,strsplit(focal_onto[which(focal_onto=="[Term]")+2],"name: "))[,2])
	if(length(which(!grepl("name: ",focal_onto[which(focal_onto=="[Term]")+2])))>0){
		onto_temp<-onto_temp[-which(!grepl("name: ",focal_onto[which(focal_onto=="[Term]")+2])),]
	}
	return(onto_temp)
}

#

annot_onto<-function(focal_df,focal_var,focal_onto){
	temp<-focal_df[,focal_var]
	temp[which(!is.na(match(focal_df[,focal_var],focal_onto$id)))]<-focal_onto[na.omit(match(focal_df[,focal_var],focal_onto$id)),"description"]
	return(temp)
}


# 

plot_embryo_TissueFC<-function(cell_vec,
	color_label=c("Insignificant","Significant (Intestine)","Significant (Muscle)","Significant (Neuronal system)","Significant (Pharynx)","Significant (Skin)"),
	color_name=c("grey","#F8766D","#A3A500","#79DEEC","#00B0F6","#E76BF3"),
	title1,xlab1="x-coordinate",ylab1="y-coordinate",
	showLegend=T){

	coo_x<-apply(q1_cell350_alignment_xyz[,1:714],2,mean)
	coo_y<-apply(q1_cell350_alignment_xyz[,c(1:714)+714],2,mean)
	coo_z<-apply(q1_cell350_alignment_xyz[,c(1:714)+714*2],2,mean)
	
	col1=rep(color_label[1],714)
	col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Intestine"),1],cell_vec),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[2]
	col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Muscle"),1],cell_vec),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[3]
	col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Neuronal system"),1],cell_vec),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[4]
	col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Pharynx"),1],cell_vec),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[5]
	col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Skin"),1],cell_vec),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[6]

	alpha1<-rep(0.25,length(col1))
	alpha1[col1!=color_label[1]]<-1
	
	dataTemp<-data.frame(coo_x,coo_y,coo_z,col1=factor(col1,level=color_label[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})]))
	ggplot(dataTemp,aes(x=coo_x,y=coo_y,col=col1))+
	geom_point(alpha=alpha1,show.legend=showLegend)+
	scale_color_manual(values=color_name[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})])+
	labs(x=xlab1,y=ylab1,title=title1)


}



#######################################################################

#######################################################################

library(elitism)
library(ggplot2)
library(pheatmap)


focal_projectFolderPath<-r"(F:\wangj\Desktop\temp\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/data/subsetDiff"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info_short.RData"),projectFolderPath=focal_projectFolderPath)



##############################################################

##############################################################

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


#################################################################################################
                                    extract ontoloy
#################################################################################################


anatomyOntology<-readLines(paste0(focal_projectFolderPath,"/wormbase/ontology/anatomy_ontology.WS287.obo"))
onto_ana<-extract_onto(anatomyOntology)

#

developmentOntology<-readLines(paste0(focal_projectFolderPath,"/wormbase/ontology/development_ontology.WS287.obo"))
onto_dev<-extract_onto(developmentOntology)

#

geneOntology<-readLines(paste0(focal_projectFolderPath,"/wormbase/ontology/gene_ontology.WS287.obo"))
onto_GO<-extract_onto(geneOntology)

#

phenotypeOntology<-readLines(paste0(focal_projectFolderPath,"/wormbase/ontology/phenotype_ontology.WS287.obo"))
onto_phe<-extract_onto(phenotypeOntology)


save(onto_ana,file="onto_ana.RData")
save(onto_dev,file="onto_dev.RData")
save(onto_GO,file="onto_GO.RData")
save(onto_phe,file="onto_phe.RData")



##################################################################################################
                                           # _env
##################################################################################################


### embryo index by environment 

envInd_list<-list(C=c(1929:1943),
	HS=c(1983:1997),
	Bacteria=c(1980:1982,2015:2016,1944:1946,1998:2000,1947:1949),	
	SA=c(2025:2039),
	GLU=c(1965:1979),
	NaCl=c(2001:2014),
	Fe=c(1950:1964),
	PQ=c(2017:2024))

save(envInd_list,file="envInd_list.RData")


### test the difference of k-wt and k-env

#x

kTest_env_x<-lapply(1:8,function(i){
	print(i)
	return(kTest(X=X_mo_x,Y=delta_X_x,refInd=c(1:105),perturbInd=envInd_list[[i]]))

})
names(kTest_env_x)<-names(envInd_list)

kTest_env_x_p<-get_kTest_p(kTest_env_x)

kTest_env_x_adjP<-pMat.adjust(kTest_env_x_p)

lapply(1:dim(kTest_env_x_adjP)[2],function(i){rownames(kTest_env_x_adjP)[which(kTest_env_x_adjP[,i]<0.05)]})

save(kTest_env_x,file="kTest_env_x.RData")
save(kTest_env_x_adjP,file="kTest_env_x_adjP.RData")


#y

kTest_env_y<-lapply(1:8,function(i){
	print(i)
	return(kTest(X=X_mo_y,Y=delta_X_y,refInd=c(1:105),perturbInd=envInd_list[[i]]))

})
names(kTest_env_y)<-names(envInd_list)

kTest_env_y_p<-get_kTest_p(kTest_env_y)

kTest_env_y_adjP<-pMat.adjust(kTest_env_y_p)

lapply(1:dim(kTest_env_y_adjP)[2],function(i){rownames(kTest_env_y_adjP)[which(kTest_env_y_adjP[,i]<0.05)]})

save(kTest_env_y,file="kTest_env_y.RData")
save(kTest_env_y_adjP,file="kTest_env_y_adjP.RData")

#z

kTest_env_z<-lapply(1:8,function(i){
	print(i)
	return(kTest(X=X_mo_z,Y=delta_X_z,refInd=c(1:105),perturbInd=envInd_list[[i]]))

})
names(kTest_env_z)<-names(envInd_list)

kTest_env_z_p<-get_kTest_p(kTest_env_z)

kTest_env_z_adjP<-pMat.adjust(kTest_env_z_p)

lapply(1:dim(kTest_env_z_adjP)[2],function(i){rownames(kTest_env_z_adjP)[which(kTest_env_z_adjP[,i]<0.05)]})

save(kTest_env_z,file="kTest_env_z.RData")
save(kTest_env_z_adjP,file="kTest_env_z_adjP.RData")


###

env_subset<-lapply(1:length(envInd_list),function(i){
	return(unique(do.call(rbind,strsplit(c(rownames(kTest_env_x_adjP)[which(kTest_env_x_adjP[,i]<0.05)],
	rownames(kTest_env_y_adjP)[which(kTest_env_y_adjP[,i]<0.05)],
	rownames(kTest_env_z_adjP)[which(kTest_env_z_adjP[,i]<0.05)]),"_"))[,3]))
})

names(env_subset)<-names(envInd_list)

save(env_subset,file="env_subset.RData")

# TissueFC

TissueFC_pvalue<-sapply(env_subset,get_chisq_pvalue,TissueFC_group)

TissueFC_res<-data.frame(
	embSubset=names(env_subset),
	chisq_p=signif(TissueFC_pvalue,3),
	adjP=signif(p.adjust(TissueFC_pvalue,method="BH"),3),
	Group=rep("TissueFC",8))

TissueFC_res_levelP<-t(sapply(1:dim(TissueFC_res)[1],function(i){
	print(i)
	return(get_chisq_levelP(env_subset[[TissueFC_res[i,1]]],TissueFC_group))

}))

TissueFC_res<-data.frame(TissueFC_res[,1:3],signif(TissueFC_res_levelP,3))
names(TissueFC_res)<-c(names(TissueFC_res)[1:3],names(table(TissueFC_group[,2])))

TissueFC_res$description<-c("Temperature","Heat shock","Bacterial food","Stearic acid","Glucose","NaCl","FeCl3","Paraquat")

TissueFC_res_env<-TissueFC_res
save(TissueFC_res_env,file="TissueFC_res_env.RData")



##################################################################################################
                                           # _anatomy
##################################################################################################


# anatomy ontology curation

assoc_anatomy<-read.table(paste0(focal_projectFolderPath,"/wormbase/ontology/anatomy_association.WS287.wb"),sep="\t",comment.char ="!")

assoc_anatomy_1<-assoc_anatomy[assoc_anatomy[,4]=="Certain",]

overlapInd<-sort(unlist(lapply(embryoGene_uniq,function(x){which(assoc_anatomy_1[,3]==x)})))

assoc_anatomy_2<-assoc_anatomy_1[overlapInd,]

#curatedOntology_anatomy<-names(table(assoc_anatomy_2[,5]))[table(assoc_anatomy_2[,5])>=10]
curatedOntology_anatomy<-c("WBbt:0006750","WBbt:0005784")

anatomyInd_list<-lapply(curatedOntology_anatomy,function(x){
	focal_geneSet<-assoc_anatomy_2[assoc_anatomy_2[,5]==x,3]
	return(sort(unique(unlist(lapply(focal_geneSet,function(y){
		which(embryoGene==y)+105
	})))))
})

names(anatomyInd_list)<-curatedOntology_anatomy

###

neuron_list<-read.csv("F:/wangj/Desktop/temp/canalizationENV/data/annotation/neuron_list.csv")

neuron_list<-neuron_list[grep("muscle",neuron_list[,4]),1]

neuron_id<-onto_ana[na.omit(match(neuron_list,onto_ana[,2])),1]

neuron_gene<-unique(unlist(lapply(neuron_id,function(x){assoc_anatomy_2[assoc_anatomy_2[,5]==x,3]})))

sapply(lapply(neuron_id,function(x){assoc_anatomy_2[assoc_anatomy_2[,5]==x,3]}),length)

#

neuron_embInd<-sort(unique(unlist(lapply(neuron_gene,function(x){which(embryoGene==x)+105}))))

anatomyInd_list<-c(anatomyInd_list,list(neuron_embInd))
names(anatomyInd_list)<-c(curatedOntology_anatomy,"individual neurons controlling muscle")


sort(sapply(anatomyInd_list,length))

save(anatomyInd_list,file="anatomyInd_list.RData")


### for anatomy

rowCol1=2
#x

kTest_anatomy_x<-lapply(1:length(anatomyInd_list),function(i){
	print(i)
	return(kTest(X=X_mo_x,Y=delta_X_x,refInd=c(1:105),perturbInd=anatomyInd_list[[i]]))

})
names(kTest_anatomy_x)<-names(anatomyInd_list)

kTest_anatomy_x_p<-get_kTest_p(kTest_anatomy_x)

kTest_anatomy_x_adjP<-pMat.adjust(kTest_anatomy_x_p,rowCol=rowCol1)

lapply(1:dim(kTest_anatomy_x_adjP)[2],function(i){rownames(kTest_anatomy_x_adjP)[which(kTest_anatomy_x_adjP[,i]<0.05)]})

save(kTest_anatomy_x,file="kTest_anatomy_x.RData")
save(kTest_anatomy_x_adjP,file="kTest_anatomy_x_adjP.RData")


#y

kTest_anatomy_y<-lapply(1:length(anatomyInd_list),function(i){
	print(i)
	return(kTest(X=X_mo_y,Y=delta_X_y,refInd=c(1:105),perturbInd=anatomyInd_list[[i]]))

})
names(kTest_anatomy_y)<-names(anatomyInd_list)

kTest_anatomy_y_p<-get_kTest_p(kTest_anatomy_y)

kTest_anatomy_y_adjP<-pMat.adjust(kTest_anatomy_y_p,rowCol=rowCol1)

lapply(1:dim(kTest_anatomy_y_adjP)[2],function(i){rownames(kTest_anatomy_y_adjP)[which(kTest_anatomy_y_adjP[,i]<0.05)]})

save(kTest_anatomy_y,file="kTest_anatomy_y.RData")
save(kTest_anatomy_y_adjP,file="kTest_anatomy_y_adjP.RData")

#z

kTest_anatomy_z<-lapply(1:length(anatomyInd_list),function(i){
	print(i)
	return(kTest(X=X_mo_z,Y=delta_X_z,refInd=c(1:105),perturbInd=anatomyInd_list[[i]]))

})
names(kTest_anatomy_z)<-names(anatomyInd_list)

kTest_anatomy_z_p<-get_kTest_p(kTest_anatomy_z)

kTest_anatomy_z_adjP<-pMat.adjust(kTest_anatomy_z_p,rowCol=rowCol1)

lapply(1:dim(kTest_anatomy_z_adjP)[2],function(i){rownames(kTest_anatomy_z_adjP)[which(kTest_anatomy_z_adjP[,i]<0.05)]})

save(kTest_anatomy_z,file="kTest_anatomy_z.RData")
save(kTest_anatomy_z_adjP,file="kTest_anatomy_z_adjP.RData")


###

anatomy_subset<-lapply(1:length(anatomyInd_list),function(i){
	return(unique(do.call(rbind,strsplit(c(rownames(kTest_anatomy_x_adjP)[which(kTest_anatomy_x_adjP[,i]<0.05)],
	rownames(kTest_anatomy_y_adjP)[which(kTest_anatomy_y_adjP[,i]<0.05)],
	rownames(kTest_anatomy_z_adjP)[which(kTest_anatomy_z_adjP[,i]<0.05)]),"_"))[,3]))
})

names(anatomy_subset)<-names(anatomyInd_list)

save(anatomy_subset,file="anatomy_subset.RData")

# TissueFC

TissueFC_pvalue<-sapply(anatomy_subset,get_chisq_pvalue,TissueFC_group)

TissueFC_res<-data.frame(
	embSubset=names(anatomy_subset),
	chisq_p=signif(TissueFC_pvalue,3),
	adjP=signif(p.adjust(TissueFC_pvalue,method="BH"),3),
	Group=rep("TissueFC",length(anatomy_subset)))

	
TissueFC_res_levelP<-t(sapply(1:dim(TissueFC_res)[1],function(i){
	print(i)
	return(get_chisq_levelP(anatomy_subset[[TissueFC_res[i,1]]],TissueFC_group))

}))

TissueFC_res<-data.frame(TissueFC_res[,1:3],TissueFC_res_levelP)
names(TissueFC_res)<-c(names(TissueFC_res)[1:3],names(table(TissueFC_group[,2])))

TissueFC_res$description<-annot_onto(TissueFC_res,"embSubset",onto_ana)

TissueFC_res$cellCount<-sapply(anatomy_subset,length)

TissueFC_res[order(TissueFC_res$cellCount),]

TissueFC_res_anatomy<-TissueFC_res
save(TissueFC_res_anatomy,file="TissueFC_res_anatomy.RData")



##################################################################################################
                                           # _phenotype
##################################################################################################

# phenotype ontology curation

assoc_phenotype<-read.table(paste0(focal_projectFolderPath,"/wormbase/ontology/phenotype_association.WS287.wb"),sep="\t",comment.char ="!")

assoc_phenotype_1<-assoc_phenotype[assoc_phenotype[,4]!="NOT",]

overlapInd<-sort(unlist(lapply(embryoGene_uniq,function(x){which(assoc_phenotype_1[,3]==x)})))

assoc_phenotype_2<-assoc_phenotype_1[overlapInd,]

#curatedOntology_phenotype<-names(table(assoc_phenotype_2[,5]))[table(assoc_phenotype_2[,5])>=10]

curatedOntology_phenotype<-c("WBPhenotype:0001663","WBPhenotype:0000016","WBPhenotype:0000469","WBPhenotype:0000017")

phenotypeInd_list<-lapply(curatedOntology_phenotype,function(x){
	focal_geneSet<-assoc_phenotype_2[assoc_phenotype_2[,5]==x,3]
	return(sort(unique(unlist(lapply(focal_geneSet,function(y){
		which(embryoGene==y)+105
	})))))
})

sort(sapply(phenotypeInd_list,length))

names(phenotypeInd_list)<-curatedOntology_phenotype

save(phenotypeInd_list,file="phenotypeInd_list.RData")

### for phenotype
rowCol1=2
#x

kTest_phenotype_x<-lapply(1:length(phenotypeInd_list),function(i){
	print(i)
	return(kTest(X=X_mo_x,Y=delta_X_x,refInd=c(1:105),perturbInd=phenotypeInd_list[[i]]))

})
names(kTest_phenotype_x)<-names(phenotypeInd_list)

kTest_phenotype_x_p<-get_kTest_p(kTest_phenotype_x)

kTest_phenotype_x_adjP<-pMat.adjust(kTest_phenotype_x_p,rowCol=rowCol1)

lapply(1:dim(kTest_phenotype_x_adjP)[2],function(i){rownames(kTest_phenotype_x_adjP)[which(kTest_phenotype_x_adjP[,i]<0.05)]})

save(kTest_phenotype_x,file="kTest_phenotype_x.RData")
save(kTest_phenotype_x_adjP,file="kTest_phenotype_x_adjP.RData")


#y

kTest_phenotype_y<-lapply(1:length(phenotypeInd_list),function(i){
	print(i)
	return(kTest(X=X_mo_y,Y=delta_X_y,refInd=c(1:105),perturbInd=phenotypeInd_list[[i]]))

})
names(kTest_phenotype_y)<-names(phenotypeInd_list)

kTest_phenotype_y_p<-get_kTest_p(kTest_phenotype_y)

kTest_phenotype_y_adjP<-pMat.adjust(kTest_phenotype_y_p,rowCol=rowCol1)

lapply(1:dim(kTest_phenotype_y_adjP)[2],function(i){rownames(kTest_phenotype_y_adjP)[which(kTest_phenotype_y_adjP[,i]<0.05)]})

save(kTest_phenotype_y,file="kTest_phenotype_y.RData")
save(kTest_phenotype_y_adjP,file="kTest_phenotype_y_adjP.RData")

#z

kTest_phenotype_z<-lapply(1:length(phenotypeInd_list),function(i){
	print(i)
	return(kTest(X=X_mo_z,Y=delta_X_z,refInd=c(1:105),perturbInd=phenotypeInd_list[[i]]))

})
names(kTest_phenotype_z)<-names(phenotypeInd_list)

kTest_phenotype_z_p<-get_kTest_p(kTest_phenotype_z)

kTest_phenotype_z_adjP<-pMat.adjust(kTest_phenotype_z_p,rowCol=rowCol1)

lapply(1:dim(kTest_phenotype_z_adjP)[2],function(i){rownames(kTest_phenotype_z_adjP)[which(kTest_phenotype_z_adjP[,i]<0.05)]})

save(kTest_phenotype_z,file="kTest_phenotype_z.RData")
save(kTest_phenotype_z_adjP,file="kTest_phenotype_z_adjP.RData")


### 

phenotype_subset<-lapply(1:length(phenotypeInd_list),function(i){
	return(unique(do.call(rbind,strsplit(c(rownames(kTest_phenotype_x_adjP)[which(kTest_phenotype_x_adjP[,i]<0.05)],
	rownames(kTest_phenotype_y_adjP)[which(kTest_phenotype_y_adjP[,i]<0.05)],
	rownames(kTest_phenotype_z_adjP)[which(kTest_phenotype_z_adjP[,i]<0.05)]),"_"))[,3]))
})

names(phenotype_subset)<-names(phenotypeInd_list)


save(phenotype_subset,file="phenotype_subset.RData")

# TissueFC

TissueFC_pvalue<-sapply(phenotype_subset,get_chisq_pvalue,TissueFC_group)

TissueFC_res<-data.frame(
	embSubset=names(phenotype_subset),
	chisq_p=signif(TissueFC_pvalue,3),
	adjP=signif(p.adjust(TissueFC_pvalue,method="BH"),3),
	Group=rep("TissueFC",length(phenotype_subset)))

TissueFC_res_levelP<-t(sapply(1:dim(TissueFC_res)[1],function(i){
	print(i)
	return(get_chisq_levelP(phenotype_subset[[TissueFC_res[i,1]]],TissueFC_group))

}))

TissueFC_res<-data.frame(TissueFC_res[,1:3],TissueFC_res_levelP)
names(TissueFC_res)<-c(names(TissueFC_res)[1:3],names(table(TissueFC_group[,2])))

TissueFC_res$description<-annot_onto(TissueFC_res,"embSubset",onto_phe)

TissueFC_res$cellCount<-sapply(phenotype_subset,length)


temp<-na.omit(TissueFC_group[match(phenotype_subset[[4]],TissueFC_group[,1]),])
temp[temp[,2]=="Pharynx",1]

temp2<-na.omit(TissueFC_group[match(phenotype_subset[[2]],TissueFC_group[,1]),])
temp2[temp2[,2]=="Pharynx",1]

na.omit(TissueFC_group[match(phenotype_subset[[1]],TissueFC_group[,1]),])

TissueFC_res_phenotype<-TissueFC_res
save(TissueFC_res_phenotype,file="TissueFC_res_phenotype.RData")


#########################################################################################3

#########################################################################################3

#########################################################################################3










### env 

color_labels<-c("Insignificant",names(table(TissueFC_group[,2])))
color_names<-c("grey","#F8766D","#A3A500","#79DEEC","#00B0F6","#E76BF3")
coo_x<-apply(q1_cell350_alignment_xyz[,1:714],2,mean)
coo_y<-apply(q1_cell350_alignment_xyz[,c(1:714)+714],2,mean)
coo_z<-apply(q1_cell350_alignment_xyz[,c(1:714)+714*2],2,mean)
	
dataTemp<-do.call(rbind,lapply(1:length(table(TissueFC_group[,2])),function(i){
	
	col1=rep(color_labels[1],714)
	col1[match(paste0(TissueFC_group[which(TissueFC_group[,2]==names(table(TissueFC_group[,2]))[i]),1],"_x"),colnames(q1_cell350_alignment_xyz))]<-color_labels[i+1]

	alpha1<-rep(0.25,714)
	alpha1[col1!=color_labels[1]]<-1
	
	return(data.frame(coo_x,coo_y,coo_z,col1=factor(col1,level=color_labels),alpha1,Tissue=rep(names(table(TissueFC_group[,2]))[i],714)))

}))


ggplot(dataTemp,aes(x=coo_x,y=coo_y,col=col1))+
geom_point(alpha=dataTemp$alpha1,show.legend=FALSE,size=0.2)+
scale_color_manual(values=color_names)+
labs(x="x-coordinate",y="y-coordinate")+
scale_y_continuous(breaks=c(-50,0,50))+
scale_x_continuous(breaks=c(-100,0,100))+
facet_wrap(vars(Tissue),nrow=1)+
theme(strip.background=element_blank())

ggsave(paste0(focal_projectFolderPath,"/figs/embryo_layout_five_tissues.pdf"),width=unit(6.4,"in"),height=unit(1.3,"in"))



###

dataTemp<-data.frame(Mother=c(c(X_mo_y[envInd_list[[3]],c(585,586,593,594)]),c(X_mo_y[1:105,c(585,586,593,594)])),
	Daughter=c(c(X_da1_y[envInd_list[[3]],c(585,586,593,594)]),c(X_da1_y[1:105,c(585,586,593,594)])),
	Embryo_group=factor(c(rep("Bacterial",length(envInd_list[[3]])*4),rep("Wild-type",105*4)),level=c("Wild-type","Bacterial")),
	Cell_pair=c(c(sapply(colnames(X_da1_y)[c(585,586,593,594)],rep,length(envInd_list[[3]]))),
		c(sapply(colnames(X_da1_y)[c(585,586,593,594)],rep,105))))

ggplot(dataTemp,aes(x=Mother,y=Daughter))+geom_point(colour=colors()[182],size=0.3,stroke=0.2)+geom_smooth(method="lm",lwd=0.5)+
facet_grid(vars(Embryo_group),vars(Cell_pair))+
labs(x="Noise of mother",y="Noise of daughter",title="")+
theme(strip.text=element_text(size=7))

ggsave(paste0(focal_projectFolderPath,"/figs/kTest_Bacterial_food_kdiff.pdf"),width=unit(2.87,"in"),height=unit(1.9,"in"))

###

i=3

color_label=c("Insignificant","Intestine","Muscle","Neuronal system","Pharynx","Skin")
color_name=c("grey","#F8766D","#A3A500","#79DEEC","#00B0F6","#E76BF3")
coo_x<-apply(q1_cell350_alignment_xyz[,1:714],2,mean)
coo_y<-apply(q1_cell350_alignment_xyz[,c(1:714)+714],2,mean)
coo_z<-apply(q1_cell350_alignment_xyz[,c(1:714)+714*2],2,mean)
	
col1=rep(color_label[1],714)
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Intestine"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[2]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Muscle"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[3]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Neuronal system"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[4]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Pharynx"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[5]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Skin"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[6]

alpha1<-rep(0.25,length(col1))
alpha1[col1!=color_label[1]]<-1
	
dataTemp<-data.frame(coo_x,coo_y,coo_z,col1=factor(col1,level=color_label[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})]))
ggplot(dataTemp,aes(x=coo_x,y=coo_y,col=col1))+
geom_point(alpha=alpha1,show.legend=TRUE)+
scale_color_manual(values=color_name[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})])+
labs(x="x-coordinate",y="y-coordinate",title="Environment: bacterial food")

ggsave(paste0(focal_projectFolderPath,"/figs/kTest_Bacterial_food_distribution.pdf"),width=unit(3.85,"in"),height=unit(1.7,"in"))


### 

pdf(paste0(focal_projectFolderPath,"/figs/env_tissue_association.pdf"),height=2.3,width=2.3)

dataTemp<--log10(TissueFC_res[,4:8])
rownames(dataTemp)<-TissueFC_res[,"description"]
pheatmap(t(dataTemp),angle_col=90,border_color=NA,cellwidth=7,cellheight=7,
fontsize_row=7,fontsize_col=7,main="Tissue type",cluster_rows=F,cluster_cols=F,legend_breaks=c(0:3))

sapply(which(p.adjust(as.matrix(TissueFC_res[,4:8]),method="BH")<0.05),rowColRec,dim(TissueFC_res[,4:8])[1],dim(TissueFC_res[,4:8])[2])
sapply(which(TissueFC_res[,4:8]<0.05),rowColRec,dim(TissueFC_res[,4:8])[1],dim(TissueFC_res[,4:8])[2])

dev.off()




###################


temp<-na.omit(TissueFC_group[match(anatomy_subset[[3]],TissueFC_group[,1]),])
temp[order(temp[,2]),]

dataTemp<--log10(TissueFC_res[,4:8])
rownames(dataTemp)<-TissueFC_res[,"description"]
#dataTemp<-dataTemp[which(apply(dataTemp,1,function(x){length(which(x<0.05))})>0),]
dataTemp<-na.omit(dataTemp)
pheatmap(t(dataTemp),angle_col=45,cellwidth=15,cellheight=15,main="TissueFC",cluster_rows=T,cluster_cols=T)

sapply(which(p.adjust(as.matrix(TissueFC_res[,4:8]),method="BH")<0.05),rowColRec,dim(TissueFC_res[,4:8])[1],dim(TissueFC_res[,4:9])[2])

plot_embryo_TissueFC(anatomy_subset[[20]],title1=TissueFC_res$description[20])

###

i=30
color_label=c("Insignificant","Intestine","Muscle","Neuronal system","Pharynx","Skin")
color_name=c("grey","#F8766D","#A3A500","#79DEEC","#00B0F6","#E76BF3")
coo_x<-apply(q1_cell350_alignment_xyz[,1:714],2,mean)
coo_y<-apply(q1_cell350_alignment_xyz[,c(1:714)+714],2,mean)
coo_z<-apply(q1_cell350_alignment_xyz[,c(1:714)+714*2],2,mean)
	
col1=rep(color_label[1],714)
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Intestine"),1],anatomy_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[2]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Muscle"),1],anatomy_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[3]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Neuronal system"),1],anatomy_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[4]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Pharynx"),1],anatomy_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[5]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Skin"),1],anatomy_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[6]

alpha1<-rep(0.25,length(col1))
alpha1[col1!=color_label[1]]<-1
	
dataTemp<-data.frame(coo_x,coo_y,coo_z,col1=factor(col1,level=color_label[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})]))
ggplot(dataTemp,aes(x=coo_x,y=coo_y,col=col1))+
geom_point(alpha=alpha1,show.legend=FALSE,size=0.7)+
scale_color_manual(values=color_name[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})])+
labs(x="x-coordinate",y="y-coordinate",title="Anatomy: intestinal cell")

ggsave(paste0(focal_projectFolderPath,"/figs/kTest_anatomy_intestine_distribution.pdf"),width=unit(2,"in"),height=unit(1.4,"in"))


###

index=30

k_wt_x<-sapply(1:dim(X_mo_x)[2],function(i){summary(lm(delta_X_x[1:105,i]~X_mo_x[1:105,i]))$coefficients[2,1]})
k_wbp_x<-sapply(1:dim(X_mo_x)[2],function(i){summary(lm(delta_X_x[anatomyInd_list[[index]],i]~X_mo_x[anatomyInd_list[[index]],i]))$coefficients[2,1]})
k_wt_y<-sapply(1:dim(X_mo_x)[2],function(i){summary(lm(delta_X_y[1:105,i]~X_mo_y[1:105,i]))$coefficients[2,1]})
k_wbp_y<-sapply(1:dim(X_mo_x)[2],function(i){summary(lm(delta_X_y[anatomyInd_list[[index]],i]~X_mo_y[anatomyInd_list[[index]],i]))$coefficients[2,1]})
k_wt_z<-sapply(1:dim(X_mo_x)[2],function(i){summary(lm(delta_X_z[1:105,i]~X_mo_z[1:105,i]))$coefficients[2,1]})
k_wbp_z<-sapply(1:dim(X_mo_x)[2],function(i){summary(lm(delta_X_z[anatomyInd_list[[index]],i]~X_mo_z[anatomyInd_list[[index]],i]))$coefficients[2,1]})

col1_x<-rep("Insignificant",706)
col1_x[kTest_anatomy_x_adjP[,match(names(anatomy_subset)[index],names(anatomyInd_list))]<0.05]<-"Significant"
col1_y<-rep("Insignificant",706)
col1_y[kTest_anatomy_y_adjP[,match(names(anatomy_subset)[index],names(anatomyInd_list))]<0.05]<-"Significant"
col1_z<-rep("Insignificant",706)
col1_z[kTest_anatomy_z_adjP[,match(names(anatomy_subset)[index],names(anatomyInd_list))]<0.05]<-"Significant"


dataTemp<-data.frame(k_wt=c(k_wt_x,k_wt_y,k_wt_z),k_wbp=c(k_wbp_x,k_wbp_y,k_wbp_z),col1=c(col1_x,col1_y,col1_z))
dataTemp$size1<-rep(0.1,dim(dataTemp)[1])
dataTemp$size1[which(dataTemp$col1=="Significant")]=0.3

ggplot(dataTemp,aes(x=k_wt,y=k_wbp,col=col1))+
geom_point(stroke=0.2,size=dataTemp$size1)+
geom_abline(intercept=0,slope=1)+
scale_color_manual(values=c("grey","#F8766D"))+
scale_x_continuous(limits=range(c(k_wt,k_wbp)))+
scale_y_continuous(limits=range(c(k_wt,k_wbp)))+
labs(x="k of WT",y="k of RNAi",title="Intestinal cell")+
theme(legend.position=c(0.1,0.8))

ggsave(paste0(focal_projectFolderPath,"/figs/kTest_anatomy_intestine_kdiff.pdf"),width=unit(1.3,"in"),height=unit(1.4,"in"))

#

pdf(paste0(focal_projectFolderPath,"/figs/anatomy_tissue_association_intestine.pdf"),width=2.3,height=2.3)

dataTemp<--log10(TissueFC_res[,4:8])
rownames(dataTemp)<-TissueFC_res[,"description"]

dataTemp2<-dataTemp[match("intestinal cell",rownames(dataTemp)),]
pheatmap(t(dataTemp2),angle_col=90,,border_color=NA,cellwidth=7,cellheight=7,
fontsize_row=7,fontsize_col=7,main="Tissue type",cluster_rows=F,cluster_cols=F,legend_breaks=c(1,4))

dev.off()


























































##################################################################################################

######################################## end #####################################################

##################################################################################################



##################################################################################################
                                           # _env
##################################################################################################


### embryo index by environment 

envInd_list<-list(C=c(1929:1943),
	HS=c(1983:1997),
	Bacteria=c(1980:1982,2015:2016,1944:1946,1998:2000,1947:1949),	
	SA=c(2025:2039),
	GLU=c(1965:1979),
	NaCl=c(2001:2014),
	Fe=c(1950:1964),
	PQ=c(2017:2024))

save(envInd_list,file="envInd_list.RData")


### test the difference of k-wt and k-env

#x

kTest_env_x<-lapply(1:8,function(i){
	print(i)
	return(kTest(X=X_mo_x,Y=delta_X_x,refInd=c(1:105),perturbInd=envInd_list[[i]]))

})
names(kTest_env_x)<-names(envInd_list)

kTest_env_x_p<-get_kTest_p(kTest_env_x)

kTest_env_x_adjP<-pMat.adjust(kTest_env_x_p)

lapply(1:dim(kTest_env_x_adjP)[2],function(i){rownames(kTest_env_x_adjP)[which(kTest_env_x_adjP[,i]<0.05)]})

save(kTest_env_x,file="kTest_env_x.RData")
save(kTest_env_x_adjP,file="kTest_env_x_adjP.RData")


#y

kTest_env_y<-lapply(1:8,function(i){
	print(i)
	return(kTest(X=X_mo_y,Y=delta_X_y,refInd=c(1:105),perturbInd=envInd_list[[i]]))

})
names(kTest_env_y)<-names(envInd_list)

kTest_env_y_p<-get_kTest_p(kTest_env_y)

kTest_env_y_adjP<-pMat.adjust(kTest_env_y_p)

lapply(1:dim(kTest_env_y_adjP)[2],function(i){rownames(kTest_env_y_adjP)[which(kTest_env_y_adjP[,i]<0.05)]})

save(kTest_env_y,file="kTest_env_y.RData")
save(kTest_env_y_adjP,file="kTest_env_y_adjP.RData")

#z

kTest_env_z<-lapply(1:8,function(i){
	print(i)
	return(kTest(X=X_mo_z,Y=delta_X_z,refInd=c(1:105),perturbInd=envInd_list[[i]]))

})
names(kTest_env_z)<-names(envInd_list)

kTest_env_z_p<-get_kTest_p(kTest_env_z)

kTest_env_z_adjP<-pMat.adjust(kTest_env_z_p)

lapply(1:dim(kTest_env_z_adjP)[2],function(i){rownames(kTest_env_z_adjP)[which(kTest_env_z_adjP[,i]<0.05)]})

save(kTest_env_z,file="kTest_env_z.RData")
save(kTest_env_z_adjP,file="kTest_env_z_adjP.RData")


###

env_subset<-lapply(1:length(envInd_list),function(i){
	return(unique(do.call(rbind,strsplit(c(rownames(kTest_env_x_adjP)[which(kTest_env_x_adjP[,i]<0.05)],
	rownames(kTest_env_y_adjP)[which(kTest_env_y_adjP[,i]<0.05)],
	rownames(kTest_env_z_adjP)[which(kTest_env_z_adjP[,i]<0.05)]),"_"))[,3]))
})

names(env_subset)<-names(envInd_list)


# TissueFC

TissueFC_pvalue<-sapply(env_subset,get_chisq_pvalue,TissueFC_group)

TissueFC_res<-data.frame(
	embSubset=names(env_subset),
	chisq_p=signif(TissueFC_pvalue,3),
	adjP=signif(p.adjust(TissueFC_pvalue,method="BH"),3),
	Group=rep("TissueFC",8))

TissueFC_res_levelP<-t(sapply(1:dim(TissueFC_res)[1],function(i){
	print(i)
	return(get_chisq_levelP(env_subset[[TissueFC_res[i,1]]],TissueFC_group))

}))

TissueFC_res<-data.frame(TissueFC_res[,1:3],signif(TissueFC_res_levelP,3))
names(TissueFC_res)<-c(names(TissueFC_res)[1:3],names(table(TissueFC_group[,2])))

TissueFC_res$description<-c("Temperature","Heat shock","Bacterial food","Stearic acid","Glucose","NaCl","FeCl3","Paraquat")

#env_subset_sigCell<-env_subset
#save(env_subset,file="env_subset.RData")
#env_TissueFC_res<-TissueFC_res
#save(env_TissueFC_res,file="env_TissueFC_res.RData")


###

color_labels<-c("Insignificant",names(table(TissueFC_group[,2])))
color_names<-c("grey","#F8766D","#A3A500","#79DEEC","#00B0F6","#E76BF3")
coo_x<-apply(q1_cell350_alignment_xyz[,1:714],2,mean)
coo_y<-apply(q1_cell350_alignment_xyz[,c(1:714)+714],2,mean)
coo_z<-apply(q1_cell350_alignment_xyz[,c(1:714)+714*2],2,mean)
	
dataTemp<-do.call(rbind,lapply(1:length(table(TissueFC_group[,2])),function(i){
	
	col1=rep(color_labels[1],714)
	col1[match(paste0(TissueFC_group[which(TissueFC_group[,2]==names(table(TissueFC_group[,2]))[i]),1],"_x"),colnames(q1_cell350_alignment_xyz))]<-color_labels[i+1]

	alpha1<-rep(0.25,714)
	alpha1[col1!=color_labels[1]]<-1
	
	return(data.frame(coo_x,coo_y,coo_z,col1=factor(col1,level=color_labels),alpha1,Tissue=rep(names(table(TissueFC_group[,2]))[i],714)))

}))


ggplot(dataTemp,aes(x=coo_x,y=coo_y,col=col1))+
geom_point(alpha=dataTemp$alpha1,show.legend=FALSE,size=0.2)+
scale_color_manual(values=color_names)+
labs(x="x-coordinate",y="y-coordinate")+
scale_y_continuous(breaks=c(-50,0,50))+
scale_x_continuous(breaks=c(-100,0,100))+
facet_wrap(vars(Tissue),nrow=1)+
theme(strip.background=element_blank())

ggsave(paste0(focal_projectFolderPath,"/figs/embryo_layout_five_tissues.pdf"),width=unit(6.4,"in"),height=unit(1.3,"in"))



###

dataTemp<-data.frame(Mother=c(c(X_mo_y[envInd_list[[3]],c(585,586,593,594)]),c(X_mo_y[1:105,c(585,586,593,594)])),
	Daughter=c(c(X_da1_y[envInd_list[[3]],c(585,586,593,594)]),c(X_da1_y[1:105,c(585,586,593,594)])),
	Embryo_group=factor(c(rep("Bacterial",length(envInd_list[[3]])*4),rep("Wild-type",105*4)),level=c("Wild-type","Bacterial")),
	Cell_pair=c(c(sapply(colnames(X_da1_y)[c(585,586,593,594)],rep,length(envInd_list[[3]]))),
		c(sapply(colnames(X_da1_y)[c(585,586,593,594)],rep,105))))

ggplot(dataTemp,aes(x=Mother,y=Daughter))+geom_point(colour=colors()[182],size=0.3,stroke=0.2)+geom_smooth(method="lm",lwd=0.5)+
facet_grid(vars(Embryo_group),vars(Cell_pair))+
labs(x="Noise of mother",y="Noise of daughter",title="")+
theme(strip.text=element_text(size=7))

ggsave(paste0(focal_projectFolderPath,"/figs/kTest_Bacterial_food_kdiff.pdf"),width=unit(2.87,"in"),height=unit(1.9,"in"))

###

i=3

color_label=c("Insignificant","Intestine","Muscle","Neuronal system","Pharynx","Skin")
color_name=c("grey","#F8766D","#A3A500","#79DEEC","#00B0F6","#E76BF3")
coo_x<-apply(q1_cell350_alignment_xyz[,1:714],2,mean)
coo_y<-apply(q1_cell350_alignment_xyz[,c(1:714)+714],2,mean)
coo_z<-apply(q1_cell350_alignment_xyz[,c(1:714)+714*2],2,mean)
	
col1=rep(color_label[1],714)
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Intestine"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[2]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Muscle"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[3]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Neuronal system"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[4]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Pharynx"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[5]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Skin"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[6]

alpha1<-rep(0.25,length(col1))
alpha1[col1!=color_label[1]]<-1
	
dataTemp<-data.frame(coo_x,coo_y,coo_z,col1=factor(col1,level=color_label[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})]))
ggplot(dataTemp,aes(x=coo_x,y=coo_y,col=col1))+
geom_point(alpha=alpha1,show.legend=TRUE)+
scale_color_manual(values=color_name[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})])+
labs(x="x-coordinate",y="y-coordinate",title="Environment: bacterial food")

ggsave(paste0(focal_projectFolderPath,"/figs/kTest_Bacterial_food_distribution.pdf"),width=unit(3.85,"in"),height=unit(1.7,"in"))


### 

pdf(paste0(focal_projectFolderPath,"/figs/env_tissue_association.pdf"),height=2.3,width=2.3)

dataTemp<--log10(TissueFC_res[,4:8])
rownames(dataTemp)<-TissueFC_res[,"description"]
pheatmap(t(dataTemp),angle_col=90,border_color=NA,cellwidth=7,cellheight=7,
fontsize_row=7,fontsize_col=7,main="Tissue type",cluster_rows=F,cluster_cols=F,legend_breaks=c(0:3))

sapply(which(p.adjust(as.matrix(TissueFC_res[,4:8]),method="BH")<0.05),rowColRec,dim(TissueFC_res[,4:8])[1],dim(TissueFC_res[,4:8])[2])
sapply(which(TissueFC_res[,4:8]<0.05),rowColRec,dim(TissueFC_res[,4:8])[1],dim(TissueFC_res[,4:8])[2])

dev.off()


###2 7

#

i=2

color_label=c("Insignificant","Intestine","Muscle","Neuronal system","Pharynx","Skin")
color_name=c("grey","#F8766D","#A3A500","#79DEEC","#00B0F6","#E76BF3")
coo_x<-apply(q1_cell350_alignment_xyz[,1:714],2,mean)
coo_y<-apply(q1_cell350_alignment_xyz[,c(1:714)+714],2,mean)
coo_z<-apply(q1_cell350_alignment_xyz[,c(1:714)+714*2],2,mean)
	
col1=rep(color_label[1],714)
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Intestine"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[2]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Muscle"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[3]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Neuronal system"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[4]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Pharynx"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[5]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Skin"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[6]

alpha1<-rep(0.25,length(col1))
alpha1[col1!=color_label[1]]<-1
	
dataTemp<-data.frame(coo_x,coo_y,coo_z,col1=factor(col1,level=color_label[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})]))
ggplot(dataTemp,aes(x=coo_x,y=coo_y,col=col1))+
geom_point(alpha=alpha1,show.legend=FALSE,size=0.7)+
scale_color_manual(values=color_name[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})])+
labs(x="x-coordinate",y="y-coordinate",title="Environment: heat shock")

ggsave(paste0(focal_projectFolderPath,"/figs/kTest_Heat_shock_distribution.pdf"),width=unit(2,"in"),height=unit(1.4,"in"))


#

i=7

color_label=c("Insignificant","Intestine","Muscle","Neuronal system","Pharynx","Skin")
color_name=c("grey","#F8766D","#A3A500","#79DEEC","#00B0F6","#E76BF3")
coo_x<-apply(q1_cell350_alignment_xyz[,1:714],2,mean)
coo_y<-apply(q1_cell350_alignment_xyz[,c(1:714)+714],2,mean)
coo_z<-apply(q1_cell350_alignment_xyz[,c(1:714)+714*2],2,mean)
	
col1=rep(color_label[1],714)
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Intestine"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[2]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Muscle"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[3]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Neuronal system"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[4]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Pharynx"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[5]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Skin"),1],env_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[6]

alpha1<-rep(0.25,length(col1))
alpha1[col1!=color_label[1]]<-1
	
dataTemp<-data.frame(coo_x,coo_y,coo_z,col1=factor(col1,level=color_label[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})]))
ggplot(dataTemp,aes(x=coo_x,y=coo_y,col=col1))+
geom_point(alpha=alpha1,show.legend=FALSE,size=0.7)+
scale_color_manual(values=color_name[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})])+
labs(x="x-coordinate",y="y-coordinate",title="Environment: FeCl3")

ggsave(paste0(focal_projectFolderPath,"/figs/kTest_FeCl3_distribution.pdf"),width=unit(2,"in"),height=unit(1.4,"in"))


#

intersect(env_subset[[i]],TissueFC_group[which(TissueFC_group[,2]=="Intestine"),1])
kTest_env_y_adjP[kTest_env_y_adjP[,i]<0.05,i]
kTest_env_y[[i]][kTest_env_y_adjP[,i]<0.05,]
which(kTest_env_y_adjP[,i]<0.05)
colnames(X_mo_y)[c(585, 586, 593, 594)]
colnames(delta_X_y)[c(585, 586, 593, 594)]

summary(lm(delta_X_y[envInd_list[[3]],585]~X_mo_y[envInd_list[[3]],585]))
summary(lm(delta_X_y[1:105,585]~X_mo_y[1:105,585]))

#

temp<-sapply(1:8,function(i){
	print(i)
	sapply(1:5,function(j){
		print(j)
		if(length(intersect(TissueFC_group[TissueFC_group[,2]==names(table(TissueFC_group[,2]))[j],1],env_subset[[i]]))>0){
			temp<-table(sign(coo_y[match(intersect(TissueFC_group[TissueFC_group[,2]==names(table(TissueFC_group[,2]))[j],1],env_subset[[i]]),do.call(rbind,strsplit(names(coo_y),"_y"))[,1])]))
			temp2<-table(sign(coo_y[match(TissueFC_group[TissueFC_group[,2]==names(table(TissueFC_group[,2]))[j],1],do.call(rbind,strsplit(names(coo_y),"_y"))[,1])]))
			p<-temp2[1]/sum(temp2)
			return(binom.test(temp[1],sum(temp),p)$p.value)		
		}else{
			return(NA)
		}

	})
})
rownames(temp)<-names(table(TissueFC_group[,2]))
colnames(temp)<-names(env_subset)
temp


##################################################################################################
                                           # _anatomy
##################################################################################################


# anatomy ontology curation

assoc_anatomy<-read.table(paste0(focal_projectFolderPath,"/wormbase/ontology/anatomy_association.WS287.wb"),sep="\t",comment.char ="!")

assoc_anatomy_1<-assoc_anatomy[assoc_anatomy[,4]=="Certain",]

overlapInd<-sort(unlist(lapply(embryoGene_uniq,function(x){which(assoc_anatomy_1[,3]==x)})))

assoc_anatomy_2<-assoc_anatomy_1[overlapInd,]

curatedOntology_anatomy<-names(table(assoc_anatomy_2[,5]))[table(assoc_anatomy_2[,5])>=10]


anatomyInd_list<-lapply(curatedOntology_anatomy,function(x){
	focal_geneSet<-assoc_anatomy_2[assoc_anatomy_2[,5]==x,3]
	return(sort(unique(unlist(lapply(focal_geneSet,function(y){
		which(embryoGene==y)+105
	})))))
})

sort(sapply(anatomyInd_list,length))

names(anatomyInd_list)<-curatedOntology_anatomy

save(anatomyInd_list,file="anatomyInd_list.RData")

### for anatomy

rowCol1=2
#x

kTest_anatomy_x<-lapply(1:length(anatomyInd_list),function(i){
	print(i)
	return(kTest(X=X_mo_x,Y=delta_X_x,refInd=c(1:105),perturbInd=anatomyInd_list[[i]]))

})
names(kTest_anatomy_x)<-names(anatomyInd_list)

kTest_anatomy_x_p<-get_kTest_p(kTest_anatomy_x)

kTest_anatomy_x_adjP<-pMat.adjust(kTest_anatomy_x_p,rowCol=rowCol1)

lapply(1:dim(kTest_anatomy_x_adjP)[2],function(i){rownames(kTest_anatomy_x_adjP)[which(kTest_anatomy_x_adjP[,i]<0.05)]})

save(kTest_anatomy_x,file="kTest_anatomy_x.RData")
save(kTest_anatomy_x_adjP,file="kTest_anatomy_x_adjP.RData")


#y

kTest_anatomy_y<-lapply(1:length(anatomyInd_list),function(i){
	print(i)
	return(kTest(X=X_mo_y,Y=delta_X_y,refInd=c(1:105),perturbInd=anatomyInd_list[[i]]))

})
names(kTest_anatomy_y)<-names(anatomyInd_list)

kTest_anatomy_y_p<-get_kTest_p(kTest_anatomy_y)

kTest_anatomy_y_adjP<-pMat.adjust(kTest_anatomy_y_p,rowCol=rowCol1)

lapply(1:dim(kTest_anatomy_y_adjP)[2],function(i){rownames(kTest_anatomy_y_adjP)[which(kTest_anatomy_y_adjP[,i]<0.05)]})

save(kTest_anatomy_y,file="kTest_anatomy_y.RData")
save(kTest_anatomy_y_adjP,file="kTest_anatomy_y_adjP.RData")

#z

kTest_anatomy_z<-lapply(1:length(anatomyInd_list),function(i){
	print(i)
	return(kTest(X=X_mo_z,Y=delta_X_z,refInd=c(1:105),perturbInd=anatomyInd_list[[i]]))

})
names(kTest_anatomy_z)<-names(anatomyInd_list)

kTest_anatomy_z_p<-get_kTest_p(kTest_anatomy_z)

kTest_anatomy_z_adjP<-pMat.adjust(kTest_anatomy_z_p,rowCol=rowCol1)

lapply(1:dim(kTest_anatomy_z_adjP)[2],function(i){rownames(kTest_anatomy_z_adjP)[which(kTest_anatomy_z_adjP[,i]<0.05)]})

save(kTest_anatomy_z,file="kTest_anatomy_z.RData")
save(kTest_anatomy_z_adjP,file="kTest_anatomy_z_adjP.RData")


###


anatomy_subset<-lapply(1:length(anatomyInd_list),function(i){
	return(unique(do.call(rbind,strsplit(c(rownames(kTest_anatomy_x_adjP)[which(kTest_anatomy_x_adjP[,i]<0.05)],
	rownames(kTest_anatomy_y_adjP)[which(kTest_anatomy_y_adjP[,i]<0.05)],
	rownames(kTest_anatomy_z_adjP)[which(kTest_anatomy_z_adjP[,i]<0.05)]),"_"))[,3]))
})

names(anatomy_subset)<-names(anatomyInd_list)


# TissueFC

TissueFC_pvalue<-sapply(anatomy_subset,get_chisq_pvalue,TissueFC_group)

TissueFC_res<-data.frame(
	embSubset=names(anatomy_subset),
	chisq_p=signif(TissueFC_pvalue,3),
	adjP=signif(p.adjust(TissueFC_pvalue,method="BH"),3),
	Group=rep("TissueFC",length(anatomy_subset)))

	
TissueFC_res_levelP<-t(sapply(1:dim(TissueFC_res)[1],function(i){
	print(i)
	return(get_chisq_levelP(anatomy_subset[[TissueFC_res[i,1]]],TissueFC_group))

}))

TissueFC_res<-data.frame(TissueFC_res[,1:3],TissueFC_res_levelP)
names(TissueFC_res)<-c(names(TissueFC_res)[1:3],names(table(TissueFC_group[,2])))

TissueFC_res$description<-annot_onto(TissueFC_res,"embSubset",onto_ana)

TissueFC_res$cellCount<-sapply(anatomy_subset,length)

TissueFC_res[order(TissueFC_res$cellCount),]

temp<-na.omit(TissueFC_group[match(anatomy_subset[[3]],TissueFC_group[,1]),])
temp[order(temp[,2]),]

dataTemp<--log10(TissueFC_res[,4:8])
rownames(dataTemp)<-TissueFC_res[,"description"]
#dataTemp<-dataTemp[which(apply(dataTemp,1,function(x){length(which(x<0.05))})>0),]
dataTemp<-na.omit(dataTemp)
pheatmap(t(dataTemp),angle_col=45,cellwidth=15,cellheight=15,main="TissueFC",cluster_rows=T,cluster_cols=T)

sapply(which(p.adjust(as.matrix(TissueFC_res[,4:8]),method="BH")<0.05),rowColRec,dim(TissueFC_res[,4:8])[1],dim(TissueFC_res[,4:9])[2])

plot_embryo_TissueFC(anatomy_subset[[20]],title1=TissueFC_res$description[20])

###

i=30
color_label=c("Insignificant","Intestine","Muscle","Neuronal system","Pharynx","Skin")
color_name=c("grey","#F8766D","#A3A500","#79DEEC","#00B0F6","#E76BF3")
coo_x<-apply(q1_cell350_alignment_xyz[,1:714],2,mean)
coo_y<-apply(q1_cell350_alignment_xyz[,c(1:714)+714],2,mean)
coo_z<-apply(q1_cell350_alignment_xyz[,c(1:714)+714*2],2,mean)
	
col1=rep(color_label[1],714)
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Intestine"),1],anatomy_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[2]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Muscle"),1],anatomy_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[3]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Neuronal system"),1],anatomy_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[4]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Pharynx"),1],anatomy_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[5]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Skin"),1],anatomy_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[6]

alpha1<-rep(0.25,length(col1))
alpha1[col1!=color_label[1]]<-1
	
dataTemp<-data.frame(coo_x,coo_y,coo_z,col1=factor(col1,level=color_label[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})]))
ggplot(dataTemp,aes(x=coo_x,y=coo_y,col=col1))+
geom_point(alpha=alpha1,show.legend=FALSE,size=0.7)+
scale_color_manual(values=color_name[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})])+
labs(x="x-coordinate",y="y-coordinate",title="Anatomy: intestinal cell")

ggsave(paste0(focal_projectFolderPath,"/figs/kTest_anatomy_intestine_distribution.pdf"),width=unit(2,"in"),height=unit(1.4,"in"))


###

index=30

k_wt_x<-sapply(1:dim(X_mo_x)[2],function(i){summary(lm(delta_X_x[1:105,i]~X_mo_x[1:105,i]))$coefficients[2,1]})
k_wbp_x<-sapply(1:dim(X_mo_x)[2],function(i){summary(lm(delta_X_x[anatomyInd_list[[index]],i]~X_mo_x[anatomyInd_list[[index]],i]))$coefficients[2,1]})
k_wt_y<-sapply(1:dim(X_mo_x)[2],function(i){summary(lm(delta_X_y[1:105,i]~X_mo_y[1:105,i]))$coefficients[2,1]})
k_wbp_y<-sapply(1:dim(X_mo_x)[2],function(i){summary(lm(delta_X_y[anatomyInd_list[[index]],i]~X_mo_y[anatomyInd_list[[index]],i]))$coefficients[2,1]})
k_wt_z<-sapply(1:dim(X_mo_x)[2],function(i){summary(lm(delta_X_z[1:105,i]~X_mo_z[1:105,i]))$coefficients[2,1]})
k_wbp_z<-sapply(1:dim(X_mo_x)[2],function(i){summary(lm(delta_X_z[anatomyInd_list[[index]],i]~X_mo_z[anatomyInd_list[[index]],i]))$coefficients[2,1]})

col1_x<-rep("Insignificant",706)
col1_x[kTest_anatomy_x_adjP[,match(names(anatomy_subset)[index],names(anatomyInd_list))]<0.05]<-"Significant"
col1_y<-rep("Insignificant",706)
col1_y[kTest_anatomy_y_adjP[,match(names(anatomy_subset)[index],names(anatomyInd_list))]<0.05]<-"Significant"
col1_z<-rep("Insignificant",706)
col1_z[kTest_anatomy_z_adjP[,match(names(anatomy_subset)[index],names(anatomyInd_list))]<0.05]<-"Significant"


dataTemp<-data.frame(k_wt=c(k_wt_x,k_wt_y,k_wt_z),k_wbp=c(k_wbp_x,k_wbp_y,k_wbp_z),col1=c(col1_x,col1_y,col1_z))
dataTemp$size1<-rep(0.1,dim(dataTemp)[1])
dataTemp$size1[which(dataTemp$col1=="Significant")]=0.3

ggplot(dataTemp,aes(x=k_wt,y=k_wbp,col=col1))+
geom_point(stroke=0.2,size=dataTemp$size1)+
geom_abline(intercept=0,slope=1)+
scale_color_manual(values=c("grey","#F8766D"))+
scale_x_continuous(limits=range(c(k_wt,k_wbp)))+
scale_y_continuous(limits=range(c(k_wt,k_wbp)))+
labs(x="k of WT",y="k of RNAi",title="Intestinal cell")+
theme(legend.position=c(0.1,0.8))

ggsave(paste0(focal_projectFolderPath,"/figs/kTest_anatomy_intestine_kdiff.pdf"),width=unit(1.3,"in"),height=unit(1.4,"in"))

#

pdf(paste0(focal_projectFolderPath,"/figs/anatomy_tissue_association_intestine.pdf"),width=2.3,height=2.3)

dataTemp<--log10(TissueFC_res[,4:8])
rownames(dataTemp)<-TissueFC_res[,"description"]

dataTemp2<-dataTemp[match("intestinal cell",rownames(dataTemp)),]
pheatmap(t(dataTemp2),angle_col=90,,border_color=NA,cellwidth=7,cellheight=7,
fontsize_row=7,fontsize_col=7,main="Tissue type",cluster_rows=F,cluster_cols=F,legend_breaks=c(1,4))

dev.off()


##################################################################################################
                                           # _phenotype
##################################################################################################

# phenotype ontology curation

assoc_phenotype<-read.table(paste0(focal_projectFolderPath,"/wormbase/ontology/phenotype_association.WS287.wb"),sep="\t",comment.char ="!")

assoc_phenotype_1<-assoc_phenotype[assoc_phenotype[,4]!="NOT",]

overlapInd<-sort(unlist(lapply(embryoGene_uniq,function(x){which(assoc_phenotype_1[,3]==x)})))

assoc_phenotype_2<-assoc_phenotype_1[overlapInd,]

curatedOntology_phenotype<-names(table(assoc_phenotype_2[,5]))[table(assoc_phenotype_2[,5])>=10]

phenotypeInd_list<-lapply(curatedOntology_phenotype,function(x){
	focal_geneSet<-assoc_phenotype_2[assoc_phenotype_2[,5]==x,3]
	return(sort(unique(unlist(lapply(focal_geneSet,function(y){
		which(embryoGene==y)+105
	})))))
})

sort(sapply(phenotypeInd_list,length))

names(phenotypeInd_list)<-curatedOntology_phenotype

save(phenotypeInd_list,file="phenotypeInd_list.RData")

### for phenotype
rowCol1=2
#x

kTest_phenotype_x<-lapply(1:length(phenotypeInd_list),function(i){
	print(i)
	return(kTest(X=X_mo_x,Y=delta_X_x,refInd=c(1:105),perturbInd=phenotypeInd_list[[i]]))

})
names(kTest_phenotype_x)<-names(phenotypeInd_list)

kTest_phenotype_x_p<-get_kTest_p(kTest_phenotype_x)

kTest_phenotype_x_adjP<-pMat.adjust(kTest_phenotype_x_p,rowCol=rowCol1)

lapply(1:dim(kTest_phenotype_x_adjP)[2],function(i){rownames(kTest_phenotype_x_adjP)[which(kTest_phenotype_x_adjP[,i]<0.05)]})

save(kTest_phenotype_x,file="kTest_phenotype_x.RData")
save(kTest_phenotype_x_adjP,file="kTest_phenotype_x_adjP.RData")


#y

kTest_phenotype_y<-lapply(1:length(phenotypeInd_list),function(i){
	print(i)
	return(kTest(X=X_mo_y,Y=delta_X_y,refInd=c(1:105),perturbInd=phenotypeInd_list[[i]]))

})
names(kTest_phenotype_y)<-names(phenotypeInd_list)

kTest_phenotype_y_p<-get_kTest_p(kTest_phenotype_y)

kTest_phenotype_y_adjP<-pMat.adjust(kTest_phenotype_y_p,rowCol=rowCol1)

lapply(1:dim(kTest_phenotype_y_adjP)[2],function(i){rownames(kTest_phenotype_y_adjP)[which(kTest_phenotype_y_adjP[,i]<0.05)]})

save(kTest_phenotype_y,file="kTest_phenotype_y.RData")
save(kTest_phenotype_y_adjP,file="kTest_phenotype_y_adjP.RData")

#z

kTest_phenotype_z<-lapply(1:length(phenotypeInd_list),function(i){
	print(i)
	return(kTest(X=X_mo_z,Y=delta_X_z,refInd=c(1:105),perturbInd=phenotypeInd_list[[i]]))

})
names(kTest_phenotype_z)<-names(phenotypeInd_list)

kTest_phenotype_z_p<-get_kTest_p(kTest_phenotype_z)

kTest_phenotype_z_adjP<-pMat.adjust(kTest_phenotype_z_p,rowCol=rowCol1)

lapply(1:dim(kTest_phenotype_z_adjP)[2],function(i){rownames(kTest_phenotype_z_adjP)[which(kTest_phenotype_z_adjP[,i]<0.05)]})

save(kTest_phenotype_z,file="kTest_phenotype_z.RData")
save(kTest_phenotype_z_adjP,file="kTest_phenotype_z_adjP.RData")


### phenotype - spatial PC 

phenotype_subset<-lapply(1:length(phenotypeInd_list),function(i){
	return(unique(do.call(rbind,strsplit(c(rownames(kTest_phenotype_x_adjP)[which(kTest_phenotype_x_adjP[,i]<0.05)],
	rownames(kTest_phenotype_y_adjP)[which(kTest_phenotype_y_adjP[,i]<0.05)],
	rownames(kTest_phenotype_z_adjP)[which(kTest_phenotype_z_adjP[,i]<0.05)]),"_"))[,3]))
})

names(phenotype_subset)<-names(phenotypeInd_list)



# TissueFC

TissueFC_pvalue<-sapply(phenotype_subset,get_chisq_pvalue,TissueFC_group)

TissueFC_res<-data.frame(
	embSubset=names(phenotype_subset),
	chisq_p=signif(TissueFC_pvalue,3),
	adjP=signif(p.adjust(TissueFC_pvalue,method="BH"),3),
	Group=rep("TissueFC",length(phenotype_subset)))

TissueFC_res_levelP<-t(sapply(1:dim(TissueFC_res)[1],function(i){
	print(i)
	return(get_chisq_levelP(phenotype_subset[[TissueFC_res[i,1]]],TissueFC_group))

}))

TissueFC_res<-data.frame(TissueFC_res[,1:3],TissueFC_res_levelP)
names(TissueFC_res)<-c(names(TissueFC_res)[1:3],names(table(TissueFC_group[,2])))

TissueFC_res$description<-annot_onto(TissueFC_res,"embSubset",onto_phe)

TissueFC_res$cellCount<-sapply(phenotype_subset,length)
#

TEMP<-TissueFC_res
rownames(TEMP)<-NULL
TEMP[,c(4:8)]<-signif(TEMP[,c(4:8)],3)
TEMP<-TEMP[,-c(2:3)]
TEMP<-TEMP[apply(TEMP[,2:6],1,function(x){length(which(x<0.05))})>0,]
TEMP[order(TEMP$cellCount),]


dataTemp<--log10(TissueFC_res[,4:8])
rownames(dataTemp)<-TissueFC_res[,"description"]

dataTemp<-dataTemp[which(apply(dataTemp,1,function(x){length(which(x>-log10(0.05)))})>0),]
dataTemp<-na.omit(dataTemp)
pheatmap(t(dataTemp),angle_col=45,cellwidth=10,cellheight=10,main="TissueFC",cluster_rows=T,cluster_cols=T)

#intestine
TissueFC_res[match(c("locomotion reduced","no induction of antimicrobial peptide expression after infection","egg laying variant","chemical hypersensitive"),TissueFC_res$description),]
#intersect(phenotypeInd_list[[match("reduced brood size",TissueFC_res$description)]],phenotypeInd_list[[match("egg laying variant",TissueFC_res$description)]])

#muscle
TissueFC_res[match(c("gonad development variant","exploded through vulva"),TissueFC_res$description),]

#skin 
TissueFC_res[match(c("Bacillus thuringiensis toxin hypersensitive","pathogen susceptibility increased","fewer seam cells","molt defect","gonad vesiculated","avoids bacterial lawn","aldicarb resistant","male tail morphology variant"),TissueFC_res$description),]

#neuronal 
TissueFC_res[match("fat content reduced",TissueFC_res$description),]

#pharynx 
TissueFC_res[match("foraging behavior variant",TissueFC_res$description),]


###

selectPheno<-c("foraging behavior variant",
	"fat content reduced",
	"gonad development variant",
	"exploded through vulva",
	"gonad vesiculated",
	"fewer seam cells",
	"molt defect",
	"male tail morphology variant",	
	"Bacillus thuringiensis toxin hypersensitive",
	"pathogen susceptibility increased",
	"avoids bacterial lawn",
	"aldicarb resistant",
	"no induction of antimicrobial peptide expression after infection",
	"locomotion reduced",
	"egg laying variant",
	"chemical hypersensitive")


dataTemp<--log10(TissueFC_res[,4:8])
rownames(dataTemp)<-TissueFC_res[,"description"]

dataTemp<-dataTemp[match(selectPheno,TissueFC_res$description),]


pdf(paste0(focal_projectFolderPath,"/figs/phenotype_tissue_association.pdf"),width=3,height=3)

pheatmap(t(dataTemp),angle_col=90,,border_color="grey",cellwidth=6,cellheight=6,
fontsize_row=6,fontsize_col=6,main="Tissue type",cluster_rows=F,cluster_cols=F,legend_breaks=c(0,4,8))

dev.off()


###

coo_x<-apply(q1_cell350_alignment_xyz[,1:714],2,mean)
coo_y<-apply(q1_cell350_alignment_xyz[,c(1:714)+714],2,mean)
coo_z<-apply(q1_cell350_alignment_xyz[,c(1:714)+714*2],2,mean)

color_label=c("Insignificant","Intestine","Muscle","Neuronal system","Pharynx","Skin")
color_name=c("grey","#F8766D","#A3A500","#79DEEC","#00B0F6","#E76BF3")


### 1 2 5 9 14

# 1

i=match(selectPheno,TissueFC_res$description)[1]

col1=rep(color_label[1],714)
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Intestine"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[2]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Muscle"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[3]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Neuronal system"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[4]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Pharynx"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[5]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Skin"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[6]

alpha1<-rep(0.25,length(col1))
alpha1[col1!=color_label[1]]<-1
	
dataTemp<-data.frame(coo_x,coo_y,coo_z,col1=factor(col1,level=color_label[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})]))

ggplot(dataTemp,aes(x=coo_x,y=coo_y,col=col1))+
geom_point(alpha=alpha1,show.legend=FALSE,size=0.7)+
scale_color_manual(values=color_name[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})])+
labs(x="x-coordinate",y="y-coordinate",title=TissueFC_res$description[i])


ggsave(paste0(focal_projectFolderPath,"/figs/kTest_phenotype_foraging_distribution.pdf"),width=unit(2,"in"),height=unit(1.4,"in"))

# 2

i=match(selectPheno,TissueFC_res$description)[2]

col1=rep(color_label[1],714)
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Intestine"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[2]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Muscle"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[3]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Neuronal system"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[4]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Pharynx"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[5]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Skin"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[6]

alpha1<-rep(0.25,length(col1))
alpha1[col1!=color_label[1]]<-1
	
dataTemp<-data.frame(coo_x,coo_y,coo_z,col1=factor(col1,level=color_label[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})]))

ggplot(dataTemp,aes(x=coo_x,y=coo_y,col=col1))+
geom_point(alpha=alpha1,show.legend=FALSE,size=0.7)+
scale_color_manual(values=color_name[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})])+
labs(x="x-coordinate",y="y-coordinate",title=TissueFC_res$description[i])


ggsave(paste0(focal_projectFolderPath,"/figs/kTest_phenotype_fat_distribution.pdf"),width=unit(2,"in"),height=unit(1.4,"in"))

# 5

i=match(selectPheno,TissueFC_res$description)[5]

col1=rep(color_label[1],714)
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Intestine"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[2]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Muscle"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[3]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Neuronal system"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[4]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Pharynx"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[5]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Skin"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[6]

alpha1<-rep(0.25,length(col1))
alpha1[col1!=color_label[1]]<-1
	
dataTemp<-data.frame(coo_x,coo_y,coo_z,col1=factor(col1,level=color_label[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})]))

ggplot(dataTemp,aes(x=coo_x,y=coo_y,col=col1))+
geom_point(alpha=alpha1,show.legend=FALSE,size=0.7)+
scale_color_manual(values=color_name[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})])+
labs(x="x-coordinate",y="y-coordinate",title=TissueFC_res$description[i])


ggsave(paste0(focal_projectFolderPath,"/figs/kTest_phenotype_gonad_distribution.pdf"),width=unit(2,"in"),height=unit(1.4,"in"))

# 9

i=match(selectPheno,TissueFC_res$description)[9]

col1=rep(color_label[1],714)
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Intestine"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[2]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Muscle"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[3]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Neuronal system"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[4]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Pharynx"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[5]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Skin"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[6]

alpha1<-rep(0.25,length(col1))
alpha1[col1!=color_label[1]]<-1
	
dataTemp<-data.frame(coo_x,coo_y,coo_z,col1=factor(col1,level=color_label[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})]))

ggplot(dataTemp,aes(x=coo_x,y=coo_y,col=col1))+
geom_point(alpha=alpha1,show.legend=FALSE,size=0.7)+
scale_color_manual(values=color_name[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})])+
labs(x="x-coordinate",y="y-coordinate",title=TissueFC_res$description[i])


ggsave(paste0(focal_projectFolderPath,"/figs/kTest_phenotype_Bacillus_distribution.pdf"),width=unit(2,"in"),height=unit(1.4,"in"))


# 14

i=match(selectPheno,TissueFC_res$description)[14]

col1=rep(color_label[1],714)
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Intestine"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[2]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Muscle"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[3]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Neuronal system"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[4]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Pharynx"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[5]
col1[match(paste0(intersect(TissueFC_group[which(TissueFC_group[,2]=="Skin"),1],phenotype_subset[[i]]),"_x"),colnames(q1_cell350_alignment_xyz))]<-color_label[6]

alpha1<-rep(0.25,length(col1))
alpha1[col1!=color_label[1]]<-1
	
dataTemp<-data.frame(coo_x,coo_y,coo_z,col1=factor(col1,level=color_label[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})]))

ggplot(dataTemp,aes(x=coo_x,y=coo_y,col=col1))+
geom_point(alpha=alpha1,show.legend=FALSE,size=0.7)+
scale_color_manual(values=color_name[sapply(color_label,function(x){ifelse(length(intersect(col1,x))>0,TRUE,FALSE)})])+
labs(x="x-coordinate",y="y-coordinate",title=TissueFC_res$description[i])


ggsave(paste0(focal_projectFolderPath,"/figs/kTest_phenotype_locomotion_distribution.pdf"),width=unit(2,"in"),height=unit(1.4,"in"))







##########################################################################












