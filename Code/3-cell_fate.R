
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


##

get_mo<-function(cell){
	temp0<-match(cell,lineage_info)
	if(is.na(temp0)){
		temp<-unlist(strsplit(cell,""))
		return(paste0(head(temp,length(temp)-1),collapse=""))		
	}else{
		cell_mo_rowCol<-as.numeric(rowColRec(temp0,dim(lineage_info)[1],dim(lineage_info)[2]))+c(-1,0)
		if(cell_mo_rowCol[1]==0){
			return("")
		}else{
			return(lineage_info[cell_mo_rowCol[1],cell_mo_rowCol[2]])
		}
	}
}


#

mo_fate_update<-function(cells,groups){

	cells_extend<-do.call(rbind,lapply(cells,function(x){
		print(which(cells==x))
		focalInd<-which(lineage_info==x)
		if(length(focalInd)>0){
			rowColInd<-rowColRec(focalInd[1],dim(lineage_info)[1],dim(lineage_info)[2])
			focal_cells=c(lineage_info[rowColInd[1]:dim(lineage_info)[1],which(lineage_info[rowColInd[1],]==x)])
			focal_groups=rep(groups[cells==x],length(focal_cells))
			dataTemp<-cbind(cell=focal_cells,group=focal_groups)
			dataTemp<-na.omit(dataTemp)
			dataTemp<-dataTemp[!duplicated(apply(dataTemp,1,paste0,collapse="-")),]
			return(dataTemp)
		}else{
			dataTemp<-cbind(cell=x,group=groups[cells==x])
			dataTemp<-na.omit(dataTemp)	
			return(dataTemp)			
		}
		
	}))
	
	
	cells_extend<-cells_extend[!duplicated(apply(cells_extend,1,paste0,collapse="-")),]
	
	cells<-cells_extend[,1]
	groups<-cells_extend[,2]
	
	names(cells)<-cells
	names(groups)<-cells
	
	cells_ori<-cells[-1]
	groups_ori<-groups[-1]
	cells_update<-cells
	groups_update<-groups

	
	while(length(cells_ori)<length(cells_update)){
		cells_ori<-cells_update
		groups_ori<-groups_update
		
		temp_mo<-sapply(cells_ori,get_mo)
		names(temp_mo)<-names(cells_ori)
		temp_group<-groups_ori
		names(temp_group)<-names(cells_ori)
		
		temp<-cbind(c(cells_ori,temp_mo[groups_ori!="Death"]),c(groups_ori,temp_group[groups_ori!="Death"]))
		
		rownames(temp)<-names(c(cells_ori,temp_mo[groups_ori!="Death"]))
		temp<-temp[!duplicated(apply(temp,1,paste,collapse="_"))&temp[,1]!="",]
		
		cells_update<-temp[,1]
		groups_update<-temp[,2]
		names(cells_update)<-rownames(temp)
		names(groups_update)<-rownames(temp)
		
	}
	

	coins<-sapply(cells_update,function(x){
		temp<-which(cells_update==x)
		
		if(length(temp)==1){
			return(TRUE)
		}else{
			if(length(table(groups_update[temp]))==1){
				return(TRUE)
			}else{
				return(FALSE)
			}
		}		
	})
	
	temp<-cbind(cells_update[coins],groups_update[coins])
	
	rownames(temp)<-names(cells_update)[coins]
	
	temp<-temp[!duplicated(temp[,1]),]
	
	return(temp)

}


##

check_sis<-function(cells,groups){
	
	dataTemp<-data.frame(cell_mo_da1_da2[na.omit(unique(c(match(cells,cell_mo_da1_da2[,1]),match(cells,cell_mo_da1_da2[,2]),match(cells,cell_mo_da1_da2[,3])))),])

	dataTemp$mo_tissue<-sapply(dataTemp[,1],function(x){ifelse(length(which(cells==x))==0,"plu",groups[which(cells==x)])})
	dataTemp$da1_tissue<-sapply(dataTemp[,2],function(x){ifelse(length(which(cells==x))==0,"plu",groups[which(cells==x)])})
	dataTemp$da2_tissue<-sapply(dataTemp[,3],function(x){ifelse(length(which(cells==x))==0,"plu",groups[which(cells==x)])})
	
	
	errorInd<-na.omit(match(unique(dataTemp[setdiff(1:dim(dataTemp)[1],which((dataTemp$mo_tissue==dataTemp$da1_tissue&dataTemp$mo_tissue==dataTemp$da2_tissue)|dataTemp$mo_tissue=="plu"|(dataTemp$mo_tissue==dataTemp$da1_tissue&dataTemp$da2_tissue=="Death")|(dataTemp$da1_tissue=="Death"&dataTemp$mo_tissue==dataTemp$da2_tissue))),1]),cells))
	if(length(errorInd)>0){
		cells_new<-cells[-errorInd]
		groups_new<-groups[-errorInd]	
		return(cbind(cells_new,groups_new))
		
	}else{
		return(cbind(cells,groups))
	
	}

}

mo_fate_update_sis<-function(cells,groups){
	
	dataTemp<-check_sis(cells,groups)
	cells_new<-dataTemp[,1]
	groups_new<-dataTemp[,2]
	while(length(cells_new)!=length(cells)){
		cells<-cells_new
		groups<-groups_new
		dataTemp<-check_sis(cells,groups)
		cells_new<-dataTemp[,1]
		groups_new<-dataTemp[,2]		
	}
	
	return(cbind(cells_new,groups_new))
	
}


#######################################################################

#######################################################################

focal_projectFolderPath<-r"(E:\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/data/subsetDiff"))

#

loadResources(c("E:/papers/canalizationENV/functions/rowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData"),projectFolderPath=focal_projectFolderPath)


######################################################################################################

###

Cell_list_TF<-read.csv(paste0(focal_projectFolderPath,"/data/annotation/Cell_list_TF.csv"))
TissueFC_group_temp<-Cell_list_TF[,c(2,4)]
names(TissueFC_group_temp)<-c("Cell","Tissue")
TissueFC_group_temp[TissueFC_group_temp$Cell=="MSpppaap"|TissueFC_group_temp$Cell=="MSappaap",2]<-"Germline"

TissueFC_group_temp<-TissueFC_group_temp[TissueFC_group_temp[,2]!="",]

dim(TissueFC_group_temp)

Tissue_1<-mo_fate_update(TissueFC_group_temp[,1],TissueFC_group_temp[,2])
dim(Tissue_1)


Tissue_2<-mo_fate_update_sis(Tissue_1[,1],Tissue_1[,2])
dim(Tissue_2)


dataTemp<-data.frame(cell_mo_da1_da2[na.omit(unique(c(match(Tissue_2[,1],cell_mo_da1_da2[,1]),match(Tissue_2[,1],cell_mo_da1_da2[,2]),match(Tissue_2[,1],cell_mo_da1_da2[,3])))),])

dataTemp$mo_tissue<-sapply(dataTemp[,1],function(x){ifelse(length(which(Tissue_2[,1]==x))==0,"plu",Tissue_2[which(Tissue_2[,1]==x),2])})
dataTemp$da1_tissue<-sapply(dataTemp[,2],function(x){ifelse(length(which(Tissue_2[,1]==x))==0,"plu",Tissue_2[which(Tissue_2[,1]==x),2])})
dataTemp$da2_tissue<-sapply(dataTemp[,3],function(x){ifelse(length(which(Tissue_2[,1]==x))==0,"plu",Tissue_2[which(Tissue_2[,1]==x),2])})

dataTemp2<-dataTemp[-which(dataTemp$mo_tissue=="plu"|(dataTemp$mo_tissue==dataTemp$da1_tissue&dataTemp$mo_tissue==dataTemp$da2_tissue)),]
dataTemp2[order(dataTemp2$mo),]


Tissue_2[order(Tissue_2[,1]),]
colnames(Tissue_2)<-c("Cell","Tissue")

TissueFC_group<-data.frame(Tissue_2)

# combine FC and Tissue 

write.csv(TissueFC_group,paste0(focal_projectFolderPath,"/data/annotation/TissueFC_group.csv"),row.names=F)



