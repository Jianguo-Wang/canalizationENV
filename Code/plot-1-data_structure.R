



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


library(ggplot2)
library(ggtree)
library(treeio)
library(scatterplot3d)

library(extrafont)
#font_import()
loadfonts()




focal_projectFolderPath<-r"(E:\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/figs"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"E:/papers/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info_short.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/ccl_wt.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cdt_wt.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/coo_list_q1_cell350_translation_rotation_scaling.RData"
	),projectFolderPath=focal_projectFolderPath)



#############################################################################################################

### lineage tree of a focal embryo, here is the ctr_emb1 


temp1<-unlist(lapply(rownames(ccl_wt),function(x){which(cell_mo_da1_da2$mo==x)}))
temp2<-unlist(lapply(rownames(ccl_wt),function(x){which(cell_mo_da1_da2$da1==x)}))

mo_da<-cell_mo_da1_da2[intersect(temp1,temp2),]


#

totLen<-210


i=10
tab9<-table(lineage_info[i-1,])
tab9_mo<-names(tab9)[tab9==2]
tab9_mo_tree<-sapply(tab9_mo,function(x){
	temp<-lineage_info[i,which(lineage_info[i-1,]==x)]
	tipLen<-sapply(temp,function(y){
		rowCol<-rowColRec(match(y,lineage_info),dim(lineage_info)[1],dim(lineage_info)[2])
		return(totLen-max(cdt_wt[na.omit(match(lineage_info[,rowCol[2]],rownames(cdt_wt))),1],na.rm=T))
	})
	paste0("(",paste0(paste0(temp,":",tipLen),collapse=","),")")
})

i=9
lineage_info_9<-lineage_info[,which(!duplicated(lineage_info[i,]))]
tab8<-table(lineage_info_9[i-1,])
tab8_mo<-names(tab8)[tab8==2]
tab8_mo_tree<-sapply(tab8_mo,function(x){
	temp1<-lineage_info_9[i,which(lineage_info_9[i-1,]==x)]
	temp2<-sapply(temp1,function(y){
		if(length(intersect(y,tab9_mo))==1){
			return(tab9_mo_tree[match(y,tab9_mo)])
		}else{
			return(y)
		}
	})
	tipLen<-sapply(temp1,function(z){
		tempInd<-match(z,rownames(ccl_wt))
		if(is.na(tempInd)){
			rowCol<-rowColRec(match(z,lineage_info),dim(lineage_info)[1],dim(lineage_info)[2])
			return(totLen-max(cdt_wt[na.omit(match(lineage_info[,rowCol[2]],rownames(cdt_wt))),1],na.rm=T))
		}else{
			return(ccl_wt[tempInd,1])
		}
	})
	return(paste0("(",paste0(paste0(temp2,":",tipLen),collapse=","),")"))
})


i=8
lineage_info_8<-lineage_info[,which(!duplicated(lineage_info[i,]))]
tab7<-table(lineage_info_8[i-1,])
tab7_mo<-names(tab7)[tab7==2]
tab7_mo_tree<-sapply(tab7_mo,function(x){
	temp1<-lineage_info_8[i,which(lineage_info_8[i-1,]==x)]
	temp2<-sapply(temp1,function(y){
		if(length(intersect(y,tab8_mo))==1){
			return(tab8_mo_tree[match(y,tab8_mo)])
		}else{
			return(y)
		}
	})
	tipLen<-sapply(temp1,function(z){
		tempInd<-match(z,rownames(ccl_wt))
		if(is.na(tempInd)){
			rowCol<-rowColRec(match(z,lineage_info),dim(lineage_info)[1],dim(lineage_info)[2])
			return(totLen-max(cdt_wt[na.omit(match(lineage_info[,rowCol[2]],rownames(cdt_wt))),1],na.rm=T))
		}else{
			return(ccl_wt[tempInd,1])
		}
	})	
	return(paste0("(",paste0(paste0(temp2,":",tipLen),collapse=","),")"))
})
	

i=7
lineage_info_7<-lineage_info[,which(!duplicated(lineage_info[i,]))]
tab6<-table(lineage_info_7[i-1,])
tab6_mo<-names(tab6)[tab6==2]
tab6_mo_tree<-sapply(tab6_mo,function(x){
	temp1<-lineage_info_7[i,which(lineage_info_7[i-1,]==x)]
	temp2<-sapply(temp1,function(y){
		if(length(intersect(y,tab7_mo))==1){
			return(tab7_mo_tree[match(y,tab7_mo)])
		}else{
			return(y)
		}
	})
	tipLen<-sapply(temp1,function(z){
		tempInd<-match(z,rownames(ccl_wt))
		if(is.na(tempInd)){
			rowCol<-rowColRec(match(z,lineage_info),dim(lineage_info)[1],dim(lineage_info)[2])
			return(totLen-max(cdt_wt[na.omit(match(lineage_info[,rowCol[2]],rownames(cdt_wt))),1],na.rm=T))
		}else{
			return(ccl_wt[tempInd,1])
		}
	})		
	return(paste0("(",paste0(paste0(temp2,":",tipLen),collapse=","),")"))
})

	
i=6
lineage_info_6<-lineage_info[,which(!duplicated(lineage_info[i,]))]
tab5<-table(lineage_info_6[i-1,])
tab5_mo<-names(tab5)[tab5==2]
tab5_mo_tree<-sapply(tab5_mo,function(x){
	temp1<-lineage_info_6[i,which(lineage_info_6[i-1,]==x)]
	temp2<-sapply(temp1,function(y){
		if(length(intersect(y,tab6_mo))==1){
			return(tab6_mo_tree[match(y,tab6_mo)])
		}else{
			return(y)
		}
	})
	tipLen<-sapply(temp1,function(z){
		tempInd<-match(z,rownames(ccl_wt))
		if(is.na(tempInd)){
			rowCol<-rowColRec(match(z,lineage_info),dim(lineage_info)[1],dim(lineage_info)[2])
			return(totLen-max(cdt_wt[na.omit(match(lineage_info[,rowCol[2]],rownames(cdt_wt))),1],na.rm=T))
		}else{
			return(ccl_wt[tempInd,1])
		}
	})	
	return(paste0("(",paste0(paste0(temp2,":",tipLen),collapse=","),")"))
})
	

i=5
lineage_info_5<-lineage_info[,which(!duplicated(lineage_info[i,]))]
tab4<-table(lineage_info_5[i-1,])
tab4_mo<-names(tab4)[tab4==2]
tab4_mo_tree<-sapply(tab4_mo,function(x){
	temp1<-lineage_info_5[i,which(lineage_info_5[i-1,]==x)]
	temp2<-sapply(temp1,function(y){
		if(length(intersect(y,tab5_mo))==1){
			return(tab5_mo_tree[match(y,tab5_mo)])
		}else{
			return(y)
		}
	})
	tipLen<-sapply(temp1,function(z){
		tempInd<-match(z,rownames(ccl_wt))
		return(ccl_wt[tempInd,1])
	})
	return(paste0("(",paste0(paste0(temp2,":",tipLen),collapse=","),")"))
})


i=4
lineage_info_4<-lineage_info[,which(!duplicated(lineage_info[i,]))]
tab3<-table(lineage_info_4[i-1,])
tab3_mo<-names(tab3)[tab3==2]
tab3_mo_tree<-sapply(tab3_mo,function(x){
	temp1<-lineage_info_4[i,which(lineage_info_4[i-1,]==x)]
	temp2<-sapply(temp1,function(y){
		if(length(intersect(y,tab4_mo))==1){
			return(tab4_mo_tree[match(y,tab4_mo)])
		}else{
			return(y)
		}
	})
	tipLen<-sapply(temp1,function(z){
		tempInd<-match(z,rownames(ccl_wt))
		return(ccl_wt[tempInd,1])
	})
	return(paste0("(",paste0(paste0(temp2,":",tipLen),collapse=","),")"))
})

	
i=3
lineage_info_3<-lineage_info[,which(!duplicated(lineage_info[i,]))]
tab2<-table(lineage_info_3[i-1,])
tab2_mo<-names(tab2)[tab2==2]
tab2_mo_tree<-sapply(tab2_mo,function(x){
	temp1<-lineage_info_3[i,which(lineage_info_3[i-1,]==x)]
	temp2<-sapply(temp1,function(y){
		if(length(intersect(y,tab3_mo))==1){
			return(tab3_mo_tree[match(y,tab3_mo)])
		}else{
			return(y)
		}
	})
	tipLen<-c(10,10)+cdt_wt[match(temp1,rownames(cdt_wt)),1]
	return(paste0("(",paste0(paste0(temp2,":",tipLen),collapse=","),")"))
})
	

i=2
lineage_info_2<-lineage_info[,which(!duplicated(lineage_info[i,]))]
tab1<-table(lineage_info_2[i-1,])
tab1_mo<-names(tab1)[tab1==2]
tab1_mo_tree<-sapply(tab1_mo,function(x){
	temp1<-lineage_info_2[i,which(lineage_info_2[i-1,]==x)]
	temp2<-sapply(temp1,function(y){
		if(length(intersect(y,tab2_mo))==1){
			return(tab2_mo_tree[match(y,tab2_mo)])
		}else{
			return(y)
		}
	})
	tipLen<-c(10,10)
	return(paste0("(",paste0(paste0(temp2,":",tipLen),collapse=","),")"))
})


tab0_mo_tree<-paste0(tab1_mo_tree,":10;")

write.table(tab0_mo_tree,paste0(focal_projectFolderPath,"/data/basic/tab0_mo_tree.newick"),row.names=F,col.names=F,quote=F)

#####################################################################################################################


### spatial distribution of a focal embryo for different stages 


cell_names<-rownames(coo_list_q1_cell350_translation_rotation_scaling[[1]])

lineage_ABa<-cell_names[grep("ABa",cell_names)]
lineage_ABp<-cell_names[grep("ABp",cell_names)]
lineage_EMS<-cell_names[c(grep("MS",cell_names),grep("E",cell_names))]
lineage_P2<-setdiff(cell_names,c(lineage_ABa,lineage_ABp,lineage_EMS))

#

timing_stage<-c()
timing_stage[1]<-min(cdt_wt[match(cell_names[intersect(grep("AB",cell_names),which(nchar(cell_names)==5))],rownames(cdt_wt)),1])
timing_stage[2]<-min(cdt_wt[match(cell_names[intersect(grep("AB",cell_names),which(nchar(cell_names)==6))],rownames(cdt_wt)),1])
timing_stage[3]<-min(cdt_wt[match(cell_names[intersect(grep("AB",cell_names),which(nchar(cell_names)==7))],rownames(cdt_wt)),1])
timing_stage[4]<-min(cdt_wt[match(cell_names[intersect(grep("AB",cell_names),which(nchar(cell_names)==8))],rownames(cdt_wt)),1])
timing_stage[5]<-min(cdt_wt[match(cell_names[intersect(grep("AB",cell_names),which(nchar(cell_names)==9))],rownames(cdt_wt)),1])
timing_stage[6]<-max(cdt_wt[match(cell_names[intersect(grep("AB",cell_names),which(nchar(cell_names)==9))],rownames(cdt_wt)),1])


cell_stage<-list()
cell_stage[[1]]<-intersect(cell_names,rownames(cdt_wt)[which(cdt_wt[,1]<timing_stage[1])])

cell_stage[[2]]<-setdiff(intersect(cell_names,rownames(cdt_wt)[which(cdt_wt[,1]<timing_stage[2])]),
				          intersect(cell_names,rownames(cdt_wt)[which(cdt_wt[,1]<timing_stage[1])]))
						  
cell_stage[[3]]<-setdiff(intersect(cell_names,rownames(cdt_wt)[which(cdt_wt[,1]<timing_stage[3])]),
				          intersect(cell_names,rownames(cdt_wt)[which(cdt_wt[,1]<timing_stage[2])]))
						  
cell_stage[[4]]<-setdiff(intersect(cell_names,rownames(cdt_wt)[which(cdt_wt[,1]<timing_stage[4])]),
				          intersect(cell_names,rownames(cdt_wt)[which(cdt_wt[,1]<timing_stage[3])]))
						  
cell_stage[[5]]<-setdiff(intersect(cell_names,rownames(cdt_wt)[which(cdt_wt[,1]<timing_stage[5])]),
				          intersect(cell_names,rownames(cdt_wt)[which(cdt_wt[,1]<timing_stage[4])]))
						  
cell_stage[[6]]<-setdiff(intersect(cell_names,rownames(cdt_wt)[which(cdt_wt[,1]<=timing_stage[6])]),
				          intersect(cell_names,rownames(cdt_wt)[which(cdt_wt[,1]<timing_stage[5])]))
						  
cell_stage[[7]]<-setdiff(cell_names,unlist(cell_stage[1:6]))


###

dataTemp<-data.frame(coo_list_q1_cell350_translation_rotation_scaling[[1]])

dataTemp$lineage=rep(NA,dim(dataTemp)[1])
dataTemp$lineage[match(lineage_ABa,cell_names)]<-"ABa"
dataTemp$lineage[match(lineage_ABp,cell_names)]<-"ABp"
dataTemp$lineage[match(lineage_EMS,cell_names)]<-"EMS"
dataTemp$lineage[match(lineage_P2,cell_names)]<-"P2"

dataTemp$stage<-rep(0,dim(dataTemp)[1])
dataTemp$stage[match(cell_stage[[1]],cell_names)]<-1
dataTemp$stage[match(cell_stage[[2]],cell_names)]<-2
dataTemp$stage[match(cell_stage[[3]],cell_names)]<-3
dataTemp$stage[match(cell_stage[[4]],cell_names)]<-4
dataTemp$stage[match(cell_stage[[5]],cell_names)]<-5
dataTemp$stage[match(cell_stage[[6]],cell_names)]<-6
dataTemp$stage[match(cell_stage[[7]],cell_names)]<-7


cell_stage_df<-dataTemp
save(cell_stage_df,file=paste0(focal_projectFolderPath,"/data/basic/cell_stage_df.RData"))


####################################################################################################################
								plot
####################################################################################################################

### plot: lineage tree of a focal embryo, here is the ctr_emb1 

lineage_tree<-read.newick(paste0(focal_projectFolderPath,"/data/basic/tab0_mo_tree.newick"))

ggtree(lineage_tree,lwd=0.1,root.position =-20)+theme_tree2()+geom_rootedge(lwd=0.1)+scale_y_reverse()

ggsave("lineage_tree_of_ctr_emb1.pdf",width=unit(3.4,"in"),height=unit(3.8,"in"))

plot(lineage_tree)

### plot: spatial distribution of a focal embryo 

ABa_cell<-intersect(rownames(coo_list_q1_cell350_translation_rotation_scaling[[1]]),unique(c(lineage_info[,lineage_info[3,]=="ABa"])))
ABp_cell<-intersect(rownames(coo_list_q1_cell350_translation_rotation_scaling[[1]]),unique(c(lineage_info[,lineage_info[3,]=="ABp"])))
EMS_cell<-intersect(rownames(coo_list_q1_cell350_translation_rotation_scaling[[1]]),unique(c(lineage_info[,lineage_info[3,]=="EMS"])))
P2_cell<-intersect(rownames(coo_list_q1_cell350_translation_rotation_scaling[[1]]),unique(c(lineage_info[,lineage_info[3,]=="P2"])))

color1<-rep(colors()[192],714)

color1[match(ABa_cell,rownames(coo_list_q1_cell350_translation_rotation_scaling[[1]]))]<-"#F8766D" 
color1[match(ABp_cell,rownames(coo_list_q1_cell350_translation_rotation_scaling[[1]]))]<-"#7CAE00" 
color1[match(EMS_cell,rownames(coo_list_q1_cell350_translation_rotation_scaling[[1]]))]<-"#00BFC4" 
color1[match(P2_cell,rownames(coo_list_q1_cell350_translation_rotation_scaling[[1]]))]<-"#C77CFF"


pdf("spatial_distribution_of_ctr_emb1.pdf")

s3d<-scatterplot3d(coo_list_q1_cell350_translation_rotation_scaling[[1]],xlim=c(-160,160),ylim=c(-160,160),zlim=c(-160,160),cex.symbols=0.75,color=color1,pch=16)

dev.off()



### plot: spatial distribution of a focal embryo for different stages 

dataTemp<-cell_stage_df

i=1
ggplot(dataTemp[dataTemp$stage==i,],aes(x=X,y=Y,col=lineage))+
geom_point(show.legend=FALSE,size=1.5,stroke=1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("spatial_distribution_of_ctr_emb1_stage_1.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=2
ggplot(dataTemp[dataTemp$stage==i,],aes(x=X,y=Y,col=lineage))+
geom_point(show.legend=FALSE,size=1,stroke=1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("spatial_distribution_of_ctr_emb1_stage_2.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=3
ggplot(dataTemp[dataTemp$stage==i,],aes(x=X,y=Y,col=lineage))+
geom_point(show.legend=FALSE,size=0.8,stroke=1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("spatial_distribution_of_ctr_emb1_stage_3.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=4
ggplot(dataTemp[dataTemp$stage==i,],aes(x=X,y=Y,col=lineage))+
geom_point(show.legend=FALSE,size=0.6,stroke=0.8)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("spatial_distribution_of_ctr_emb1_stage_4.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=5
ggplot(dataTemp[dataTemp$stage==i,],aes(x=X,y=Y,col=lineage))+
geom_point(show.legend=FALSE,size=0.5,stroke=0.5)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("spatial_distribution_of_ctr_emb1_stage_5.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=6
ggplot(dataTemp[dataTemp$stage==i,],aes(x=X,y=Y,col=lineage))+
geom_point(show.legend=FALSE,size=0.5,stroke=0.3)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("spatial_distribution_of_ctr_emb1_stage_6.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=7
ggplot(dataTemp[dataTemp$stage==i,],aes(x=X,y=Y,col=lineage))+
geom_point(show.legend=FALSE,size=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("spatial_distribution_of_ctr_emb1_stage_7.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


# legend

i=7
ggplot(dataTemp[dataTemp$stage==i,],aes(x=X,y=Y,col=lineage))+
geom_point(cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("spatial_distribution_of_ctr_emb1_stage_legend.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))













