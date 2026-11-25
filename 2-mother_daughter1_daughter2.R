


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


focal_projectFolderPath<-r"(F:\wangj\Desktop\temp\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/data/basic"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/data/basic/all_cells.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info.RData"),projectFolderPath=focal_projectFolderPath)


###

dataTemp<-data.frame(mo=NULL,da1=NULL)

for(i in 1:(dim(lineage_info)[1]-1)){
	focal_table<-table(lineage_info[i,])
	focal_mo<-names(focal_table)
	for(j in 1:length(focal_mo)){
		focal_table_2<-table(lineage_info[i+1,which(lineage_info[i,]==focal_mo[j])])
		focal_da<-names(focal_table_2)
		dataTemp<-rbind(dataTemp,
			data.frame(mo=rep(focal_mo[j],length(focal_da)),da1=focal_da))
	}
}

cell_mo_da1_da2<-data.frame(dataTemp,da2=sapply(1:dim(dataTemp)[1],function(i){
	tempInd<-setdiff(which(dataTemp[,1]==dataTemp[i,1]),i)
	res<-NA
	if(length(tempInd)==1){
		res<-dataTemp[tempInd,2]
	}
	return(res)
}))


save(cell_mo_da1_da2,file="cell_mo_da1_da2.RData")


