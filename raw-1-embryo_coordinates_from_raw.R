

#dirname()
#basename()
#library(tools)
#file_path_sans_ext()
#file_ext()



### raw to coo

#

raw2coo<-function(focalFile){
	
	focal_data<-read.table(focalFile,header=T)
	
	focal_cells<-unique(focal_data$cell_name)
	
	focal_coo<-lapply(focal_cells,function(x){
		temp<-focal_data[which(focal_data$cell_name==x),c(3:5,1)]
		temp[,3]<-temp[,3]*4.5
		return(temp)
	})
	
	names(focal_coo)<-focal_cells
	
	return(focal_coo)
	
}

#

raw2coo2<-function(focalFile){
	
	focal_data<-read.table(focalFile,header=T)
	
	focal_cells<-unique(focal_data$cell_name)
	
	focal_coo<-lapply(focal_cells,function(x){
		temp<-focal_data[which(focal_data$cell_name==x),c(3:5,1)]
		return(temp)
	})
	
	names(focal_coo)<-focal_cells
	
	return(focal_coo)
	
}

### get coo from dir

#

coo_from_dir<-function(focalPath){
	files<-list.files(path=focalPath,full.name=T)
	embryos<-sapply(basename(files),function(x){unlist(strsplit(x,"_raw_data.txt"))})
	
	raw2coo_list<-lapply(1:length(files),function(i){
		print(i)
		return(raw2coo(files[i]))
	})
	
	names(raw2coo_list)<-embryos
	
	return(raw2coo_list)

}

#

coo_from_dir2<-function(focalPath){
	files<-list.files(path=focalPath,full.name=T)
	embryos<-sapply(basename(files),function(x){unlist(strsplit(x,"_raw_data.txt"))})
	
	raw2coo_list<-lapply(1:length(files),function(i){
		print(i)
		return(raw2coo2(files[i]))
	})
	
	names(raw2coo_list)<-embryos
	
	return(raw2coo_list)

}

### coo to quantile 

coo2quantile<-function(focal_coo){
	
	focal_coo2quantile<-lapply(1:length(focal_coo),function(iii){
		
		print(iii)
		
		focal_embryo<-lapply(focal_coo[[iii]],function(y){
	
			if(dim(y)[1]>1){
				temp<-cbind(X=approx(y[,4],y[,1],n=11)$y,Y=approx(y[,4],y[,2],n=11)$y,Z=approx(y[,4],y[,3],n=11)$y)
			}
			if(dim(y)[1]==1){
				temp<-cbind(X=rep(y[1,1],11),Y=rep(y[1,2],11),Z=rep(y[1,3],11))
			}

			rownames(temp)<-paste0("s",1:11)
	
			return(temp)
	
		})
		
		return(focal_embryo)
		
	})
	
	names(focal_coo2quantile)<-names(focal_coo)
	
	return(focal_coo2quantile)
}


### quantile to refine 


coo2quantile_refine<-function(focal_coo2quantile,all_cells=NA){
	
	if(is.na(all_cells[1])){
		all_cells<-unique(unlist(lapply(focal_coo2quantile,function(x){names(x)})))
	}

	focal_coo2quantile_refine<-lapply(1:length(focal_coo2quantile),function(ii){
		
		print(ii)
		
		temp<-list()
	
		for(i in 1:length(all_cells)){
			if(is.na(match(all_cells[i],names(focal_coo2quantile[[ii]])))){
				temp[[i]]<-matrix(nrow=11,ncol=3)
				colnames(temp[[i]])<-c("X","Y","Z")
				rownames(temp[[i]])<-paste0("s",1:11)
			}
			if(!is.na(match(all_cells[i],names(focal_coo2quantile[[ii]])))){
				temp[[i]]<-focal_coo2quantile[[ii]][[match(all_cells[i],names(focal_coo2quantile[[ii]]))]]
			}
		
		}
	
		names(temp)<-all_cells
	
		return(temp)
	
	})
	
	names(focal_coo2quantile_refine)<-names(focal_coo2quantile)
	return(focal_coo2quantile_refine)
	
}


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


#loadResources(c("F:/wangj/Desktop/temp/canalization/functions/rowColRecovery.R",
#	"F:/wangj/Desktop/temp/canalization/data/basic/lineage_info.RData"),
#	r"(E:\papers\canalization)")


############################################################################
#					extract raw data and conduct linear interpolation
############################################################################

focal_projectFolderPath<-r"(E:\canalizationENV)"


##

setwd(paste0(focal_projectFolderPath,"/data/interpolation"))

## raw 2 coo 

coo_Control<-coo_from_dir(paste0(focal_projectFolderPath,"/raw_data/Raw_data_Control/spatial"))
coo_RNAi<-coo_from_dir(paste0(focal_projectFolderPath,"/raw_data/Raw_data_RNAi/spatial"))
coo_Env<-coo_from_dir2(paste0(focal_projectFolderPath,"/raw_data/Raw_data_Env/spatial"))

coo_Control<-coo_Control[match(paste0("ctr_emb",1:105),names(coo_Control))]

save(coo_Control,file="coo_Control.RData")
save(coo_RNAi,file="coo_RNAi.RData")
save(coo_Env,file="coo_Env.RData")


## coo to quantile to refine 
coo_quantile<-coo2quantile(c(coo_Control,coo_RNAi,coo_Env))
coo_quantile_refine<-coo2quantile_refine(coo_quantile)

save(coo_quantile,file="coo_quantile.RData")
save(coo_quantile_refine,file="coo_quantile_refine.RData")


## initial stage matrix 

coo_quantile_refine_initial_x<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_initial_x)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_initial_x)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_initial_y<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_initial_y)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_initial_y)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_initial_z<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_initial_z)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_initial_z)<-names(coo_quantile_refine[[1]])

for(i in 1:length(coo_quantile_refine)){
	print(i)
	for(j in 1:length(coo_quantile_refine[[1]])){
		coo_quantile_refine_initial_x[i,j]<-coo_quantile_refine[[i]][[j]][1,1]
		coo_quantile_refine_initial_y[i,j]<-coo_quantile_refine[[i]][[j]][1,2]
		coo_quantile_refine_initial_z[i,j]<-coo_quantile_refine[[i]][[j]][1,3]
		
	}
}

save(coo_quantile_refine_initial_x,file="coo_quantile_refine_initial_x.RData")
save(coo_quantile_refine_initial_y,file="coo_quantile_refine_initial_y.RData")
save(coo_quantile_refine_initial_z,file="coo_quantile_refine_initial_z.RData")


### Cell350

cell_select<-read.csv(paste0(focal_projectFolderPath,"/data/basic/coordinates_wt.csv"))[,1]

coo_quantile_refine_initial_x_cell350<-coo_quantile_refine_initial_x[,match(cell_select,colnames(coo_quantile_refine_initial_x))]
coo_quantile_refine_initial_y_cell350<-coo_quantile_refine_initial_y[,match(cell_select,colnames(coo_quantile_refine_initial_y))]
coo_quantile_refine_initial_z_cell350<-coo_quantile_refine_initial_z[,match(cell_select,colnames(coo_quantile_refine_initial_z))]

nonNAInd<-which(apply(cbind(coo_quantile_refine_initial_x_cell350,coo_quantile_refine_initial_y_cell350,coo_quantile_refine_initial_z_cell350),1,function(x){length(which(is.na(x)))})==0)

coo_quantile_refine_initial_x_cell350<-coo_quantile_refine_initial_x_cell350[nonNAInd,]
coo_quantile_refine_initial_y_cell350<-coo_quantile_refine_initial_y_cell350[nonNAInd,]
coo_quantile_refine_initial_z_cell350<-coo_quantile_refine_initial_z_cell350[nonNAInd,]


save(coo_quantile_refine_initial_x_cell350,file="coo_quantile_refine_initial_x_cell350.RData")
save(coo_quantile_refine_initial_y_cell350,file="coo_quantile_refine_initial_y_cell350.RData")
save(coo_quantile_refine_initial_z_cell350,file="coo_quantile_refine_initial_z_cell350.RData")


### coo list 

coo_list_q1_cell350<-lapply(1:dim(coo_quantile_refine_initial_x_cell350)[1],function(i){
	temp<-cbind(coo_quantile_refine_initial_x_cell350[i,],coo_quantile_refine_initial_y_cell350[i,],coo_quantile_refine_initial_z_cell350[i,])
	rownames(temp)<-colnames(coo_quantile_refine_initial_x_cell350)
	colnames(temp)<-c("X","Y","Z")
	
	return(temp)
	
})

names(coo_list_q1_cell350)<-rownames(coo_quantile_refine_initial_x_cell350)


save(coo_list_q1_cell350,file="coo_list_q1_cell350.RData")

################################################################################

### final stage 


coo_quantile_refine_final_x<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_final_x)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_final_x)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_final_y<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_final_y)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_final_y)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_final_z<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_final_z)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_final_z)<-names(coo_quantile_refine[[1]])

for(i in 1:length(coo_quantile_refine)){
	print(i)
	for(j in 1:length(coo_quantile_refine[[1]])){
		coo_quantile_refine_final_x[i,j]<-coo_quantile_refine[[i]][[j]][11,1]
		coo_quantile_refine_final_y[i,j]<-coo_quantile_refine[[i]][[j]][11,2]
		coo_quantile_refine_final_z[i,j]<-coo_quantile_refine[[i]][[j]][11,3]
		
	}
}

save(coo_quantile_refine_final_x,file="coo_quantile_refine_final_x.RData")
save(coo_quantile_refine_final_y,file="coo_quantile_refine_final_y.RData")
save(coo_quantile_refine_final_z,file="coo_quantile_refine_final_z.RData")


### Cell350

cell_select<-read.csv(paste0(focal_projectFolderPath,"/data/basic/coordinates_wt.csv"))[,1]

coo_quantile_refine_final_x_cell350<-coo_quantile_refine_final_x[,match(cell_select,colnames(coo_quantile_refine_final_x))]
coo_quantile_refine_final_y_cell350<-coo_quantile_refine_final_y[,match(cell_select,colnames(coo_quantile_refine_final_y))]
coo_quantile_refine_final_z_cell350<-coo_quantile_refine_final_z[,match(cell_select,colnames(coo_quantile_refine_final_z))]

nonNAInd<-which(apply(cbind(coo_quantile_refine_final_x_cell350,coo_quantile_refine_final_y_cell350,coo_quantile_refine_final_z_cell350),1,function(x){length(which(is.na(x)))})==0)

coo_quantile_refine_final_x_cell350<-coo_quantile_refine_final_x_cell350[nonNAInd,]
coo_quantile_refine_final_y_cell350<-coo_quantile_refine_final_y_cell350[nonNAInd,]
coo_quantile_refine_final_z_cell350<-coo_quantile_refine_final_z_cell350[nonNAInd,]


save(coo_quantile_refine_final_x_cell350,file="coo_quantile_refine_final_x_cell350.RData")
save(coo_quantile_refine_final_y_cell350,file="coo_quantile_refine_final_y_cell350.RData")
save(coo_quantile_refine_final_z_cell350,file="coo_quantile_refine_final_z_cell350.RData")


### coo list 

coo_list_q11_cell350<-lapply(1:dim(coo_quantile_refine_final_x_cell350)[1],function(i){
	temp<-cbind(coo_quantile_refine_final_x_cell350[i,],coo_quantile_refine_final_y_cell350[i,],coo_quantile_refine_final_z_cell350[i,])
	rownames(temp)<-colnames(coo_quantile_refine_final_x_cell350)
	colnames(temp)<-c("X","Y","Z")
	
	return(temp)
	
})

names(coo_list_q11_cell350)<-rownames(coo_quantile_refine_final_x_cell350)


save(coo_list_q11_cell350,file="coo_list_q11_cell350.RData")



################################################################################

### stage - 2 


coo_quantile_refine_2_x<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_2_x)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_2_x)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_2_y<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_2_y)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_2_y)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_2_z<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_2_z)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_2_z)<-names(coo_quantile_refine[[1]])

for(i in 1:length(coo_quantile_refine)){
	print(i)
	for(j in 1:length(coo_quantile_refine[[1]])){
		coo_quantile_refine_2_x[i,j]<-coo_quantile_refine[[i]][[j]][2,1]
		coo_quantile_refine_2_y[i,j]<-coo_quantile_refine[[i]][[j]][2,2]
		coo_quantile_refine_2_z[i,j]<-coo_quantile_refine[[i]][[j]][2,3]
		
	}
}

save(coo_quantile_refine_2_x,file="coo_quantile_refine_2_x.RData")
save(coo_quantile_refine_2_y,file="coo_quantile_refine_2_y.RData")
save(coo_quantile_refine_2_z,file="coo_quantile_refine_2_z.RData")


### Cell350

cell_select<-read.csv(paste0(focal_projectFolderPath,"/data/basic/coordinates_wt.csv"))[,1]

coo_quantile_refine_2_x_cell350<-coo_quantile_refine_2_x[,match(cell_select,colnames(coo_quantile_refine_2_x))]
coo_quantile_refine_2_y_cell350<-coo_quantile_refine_2_y[,match(cell_select,colnames(coo_quantile_refine_2_y))]
coo_quantile_refine_2_z_cell350<-coo_quantile_refine_2_z[,match(cell_select,colnames(coo_quantile_refine_2_z))]

nonNAInd<-which(apply(cbind(coo_quantile_refine_2_x_cell350,coo_quantile_refine_2_y_cell350,coo_quantile_refine_2_z_cell350),1,function(x){length(which(is.na(x)))})==0)

coo_quantile_refine_2_x_cell350<-coo_quantile_refine_2_x_cell350[nonNAInd,]
coo_quantile_refine_2_y_cell350<-coo_quantile_refine_2_y_cell350[nonNAInd,]
coo_quantile_refine_2_z_cell350<-coo_quantile_refine_2_z_cell350[nonNAInd,]


save(coo_quantile_refine_2_x_cell350,file="coo_quantile_refine_2_x_cell350.RData")
save(coo_quantile_refine_2_y_cell350,file="coo_quantile_refine_2_y_cell350.RData")
save(coo_quantile_refine_2_z_cell350,file="coo_quantile_refine_2_z_cell350.RData")


### coo list 

coo_list_q2_cell350<-lapply(1:dim(coo_quantile_refine_2_x_cell350)[1],function(i){
	temp<-cbind(coo_quantile_refine_2_x_cell350[i,],coo_quantile_refine_2_y_cell350[i,],coo_quantile_refine_2_z_cell350[i,])
	rownames(temp)<-colnames(coo_quantile_refine_2_x_cell350)
	colnames(temp)<-c("X","Y","Z")
	
	return(temp)
	
})

names(coo_list_q2_cell350)<-rownames(coo_quantile_refine_2_x_cell350)


save(coo_list_q2_cell350,file="coo_list_q2_cell350.RData")


################################################################################

### stage - 3 


coo_quantile_refine_3_x<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_3_x)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_3_x)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_3_y<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_3_y)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_3_y)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_3_z<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_3_z)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_3_z)<-names(coo_quantile_refine[[1]])

for(i in 1:length(coo_quantile_refine)){
	print(i)
	for(j in 1:length(coo_quantile_refine[[1]])){
		coo_quantile_refine_3_x[i,j]<-coo_quantile_refine[[i]][[j]][3,1]
		coo_quantile_refine_3_y[i,j]<-coo_quantile_refine[[i]][[j]][3,2]
		coo_quantile_refine_3_z[i,j]<-coo_quantile_refine[[i]][[j]][3,3]
		
	}
}

save(coo_quantile_refine_3_x,file="coo_quantile_refine_3_x.RData")
save(coo_quantile_refine_3_y,file="coo_quantile_refine_3_y.RData")
save(coo_quantile_refine_3_z,file="coo_quantile_refine_3_z.RData")


### Cell350

cell_select<-read.csv(paste0(focal_projectFolderPath,"/data/basic/coordinates_wt.csv"))[,1]

coo_quantile_refine_3_x_cell350<-coo_quantile_refine_3_x[,match(cell_select,colnames(coo_quantile_refine_3_x))]
coo_quantile_refine_3_y_cell350<-coo_quantile_refine_3_y[,match(cell_select,colnames(coo_quantile_refine_3_y))]
coo_quantile_refine_3_z_cell350<-coo_quantile_refine_3_z[,match(cell_select,colnames(coo_quantile_refine_3_z))]

nonNAInd<-which(apply(cbind(coo_quantile_refine_3_x_cell350,coo_quantile_refine_3_y_cell350,coo_quantile_refine_3_z_cell350),1,function(x){length(which(is.na(x)))})==0)

coo_quantile_refine_3_x_cell350<-coo_quantile_refine_3_x_cell350[nonNAInd,]
coo_quantile_refine_3_y_cell350<-coo_quantile_refine_3_y_cell350[nonNAInd,]
coo_quantile_refine_3_z_cell350<-coo_quantile_refine_3_z_cell350[nonNAInd,]


save(coo_quantile_refine_3_x_cell350,file="coo_quantile_refine_3_x_cell350.RData")
save(coo_quantile_refine_3_y_cell350,file="coo_quantile_refine_3_y_cell350.RData")
save(coo_quantile_refine_3_z_cell350,file="coo_quantile_refine_3_z_cell350.RData")


### coo list 

coo_list_q3_cell350<-lapply(1:dim(coo_quantile_refine_3_x_cell350)[1],function(i){
	temp<-cbind(coo_quantile_refine_3_x_cell350[i,],coo_quantile_refine_3_y_cell350[i,],coo_quantile_refine_3_z_cell350[i,])
	rownames(temp)<-colnames(coo_quantile_refine_3_x_cell350)
	colnames(temp)<-c("X","Y","Z")
	
	return(temp)
	
})

names(coo_list_q3_cell350)<-rownames(coo_quantile_refine_3_x_cell350)


save(coo_list_q3_cell350,file="coo_list_q3_cell350.RData")



################################################################################


### stage - 4 


coo_quantile_refine_4_x<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_4_x)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_4_x)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_4_y<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_4_y)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_4_y)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_4_z<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_4_z)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_4_z)<-names(coo_quantile_refine[[1]])

for(i in 1:length(coo_quantile_refine)){
	print(i)
	for(j in 1:length(coo_quantile_refine[[1]])){
		coo_quantile_refine_4_x[i,j]<-coo_quantile_refine[[i]][[j]][4,1]
		coo_quantile_refine_4_y[i,j]<-coo_quantile_refine[[i]][[j]][4,2]
		coo_quantile_refine_4_z[i,j]<-coo_quantile_refine[[i]][[j]][4,3]
		
	}
}

save(coo_quantile_refine_4_x,file="coo_quantile_refine_4_x.RData")
save(coo_quantile_refine_4_y,file="coo_quantile_refine_4_y.RData")
save(coo_quantile_refine_4_z,file="coo_quantile_refine_4_z.RData")


### Cell350

cell_select<-read.csv(paste0(focal_projectFolderPath,"/data/basic/coordinates_wt.csv"))[,1]

coo_quantile_refine_4_x_cell350<-coo_quantile_refine_4_x[,match(cell_select,colnames(coo_quantile_refine_4_x))]
coo_quantile_refine_4_y_cell350<-coo_quantile_refine_4_y[,match(cell_select,colnames(coo_quantile_refine_4_y))]
coo_quantile_refine_4_z_cell350<-coo_quantile_refine_4_z[,match(cell_select,colnames(coo_quantile_refine_4_z))]

nonNAInd<-which(apply(cbind(coo_quantile_refine_4_x_cell350,coo_quantile_refine_4_y_cell350,coo_quantile_refine_4_z_cell350),1,function(x){length(which(is.na(x)))})==0)

coo_quantile_refine_4_x_cell350<-coo_quantile_refine_4_x_cell350[nonNAInd,]
coo_quantile_refine_4_y_cell350<-coo_quantile_refine_4_y_cell350[nonNAInd,]
coo_quantile_refine_4_z_cell350<-coo_quantile_refine_4_z_cell350[nonNAInd,]


save(coo_quantile_refine_4_x_cell350,file="coo_quantile_refine_4_x_cell350.RData")
save(coo_quantile_refine_4_y_cell350,file="coo_quantile_refine_4_y_cell350.RData")
save(coo_quantile_refine_4_z_cell350,file="coo_quantile_refine_4_z_cell350.RData")


### coo list 

coo_list_q4_cell350<-lapply(1:dim(coo_quantile_refine_4_x_cell350)[1],function(i){
	temp<-cbind(coo_quantile_refine_4_x_cell350[i,],coo_quantile_refine_4_y_cell350[i,],coo_quantile_refine_4_z_cell350[i,])
	rownames(temp)<-colnames(coo_quantile_refine_4_x_cell350)
	colnames(temp)<-c("X","Y","Z")
	
	return(temp)
	
})

names(coo_list_q4_cell350)<-rownames(coo_quantile_refine_4_x_cell350)


save(coo_list_q4_cell350,file="coo_list_q4_cell350.RData")


################################################################################

### stage - 5 


coo_quantile_refine_5_x<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_5_x)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_5_x)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_5_y<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_5_y)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_5_y)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_5_z<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_5_z)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_5_z)<-names(coo_quantile_refine[[1]])

for(i in 1:length(coo_quantile_refine)){
	print(i)
	for(j in 1:length(coo_quantile_refine[[1]])){
		coo_quantile_refine_5_x[i,j]<-coo_quantile_refine[[i]][[j]][5,1]
		coo_quantile_refine_5_y[i,j]<-coo_quantile_refine[[i]][[j]][5,2]
		coo_quantile_refine_5_z[i,j]<-coo_quantile_refine[[i]][[j]][5,3]
		
	}
}

save(coo_quantile_refine_5_x,file="coo_quantile_refine_5_x.RData")
save(coo_quantile_refine_5_y,file="coo_quantile_refine_5_y.RData")
save(coo_quantile_refine_5_z,file="coo_quantile_refine_5_z.RData")


### Cell350

cell_select<-read.csv(paste0(focal_projectFolderPath,"/data/basic/coordinates_wt.csv"))[,1]

coo_quantile_refine_5_x_cell350<-coo_quantile_refine_5_x[,match(cell_select,colnames(coo_quantile_refine_5_x))]
coo_quantile_refine_5_y_cell350<-coo_quantile_refine_5_y[,match(cell_select,colnames(coo_quantile_refine_5_y))]
coo_quantile_refine_5_z_cell350<-coo_quantile_refine_5_z[,match(cell_select,colnames(coo_quantile_refine_5_z))]

nonNAInd<-which(apply(cbind(coo_quantile_refine_5_x_cell350,coo_quantile_refine_5_y_cell350,coo_quantile_refine_5_z_cell350),1,function(x){length(which(is.na(x)))})==0)

coo_quantile_refine_5_x_cell350<-coo_quantile_refine_5_x_cell350[nonNAInd,]
coo_quantile_refine_5_y_cell350<-coo_quantile_refine_5_y_cell350[nonNAInd,]
coo_quantile_refine_5_z_cell350<-coo_quantile_refine_5_z_cell350[nonNAInd,]


save(coo_quantile_refine_5_x_cell350,file="coo_quantile_refine_5_x_cell350.RData")
save(coo_quantile_refine_5_y_cell350,file="coo_quantile_refine_5_y_cell350.RData")
save(coo_quantile_refine_5_z_cell350,file="coo_quantile_refine_5_z_cell350.RData")


### coo list 

coo_list_q5_cell350<-lapply(1:dim(coo_quantile_refine_5_x_cell350)[1],function(i){
	temp<-cbind(coo_quantile_refine_5_x_cell350[i,],coo_quantile_refine_5_y_cell350[i,],coo_quantile_refine_5_z_cell350[i,])
	rownames(temp)<-colnames(coo_quantile_refine_5_x_cell350)
	colnames(temp)<-c("X","Y","Z")
	
	return(temp)
	
})

names(coo_list_q5_cell350)<-rownames(coo_quantile_refine_5_x_cell350)


save(coo_list_q5_cell350,file="coo_list_q5_cell350.RData")



################################################################################

### stage - 6


coo_quantile_refine_6_x<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_6_x)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_6_x)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_6_y<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_6_y)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_6_y)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_6_z<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_6_z)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_6_z)<-names(coo_quantile_refine[[1]])

for(i in 1:length(coo_quantile_refine)){
	print(i)
	for(j in 1:length(coo_quantile_refine[[1]])){
		coo_quantile_refine_6_x[i,j]<-coo_quantile_refine[[i]][[j]][6,1]
		coo_quantile_refine_6_y[i,j]<-coo_quantile_refine[[i]][[j]][6,2]
		coo_quantile_refine_6_z[i,j]<-coo_quantile_refine[[i]][[j]][6,3]
		
	}
}

save(coo_quantile_refine_6_x,file="coo_quantile_refine_6_x.RData")
save(coo_quantile_refine_6_y,file="coo_quantile_refine_6_y.RData")
save(coo_quantile_refine_6_z,file="coo_quantile_refine_6_z.RData")


### Cell350

cell_select<-read.csv(paste0(focal_projectFolderPath,"/data/basic/coordinates_wt.csv"))[,1]

coo_quantile_refine_6_x_cell350<-coo_quantile_refine_6_x[,match(cell_select,colnames(coo_quantile_refine_6_x))]
coo_quantile_refine_6_y_cell350<-coo_quantile_refine_6_y[,match(cell_select,colnames(coo_quantile_refine_6_y))]
coo_quantile_refine_6_z_cell350<-coo_quantile_refine_6_z[,match(cell_select,colnames(coo_quantile_refine_6_z))]

nonNAInd<-which(apply(cbind(coo_quantile_refine_6_x_cell350,coo_quantile_refine_6_y_cell350,coo_quantile_refine_6_z_cell350),1,function(x){length(which(is.na(x)))})==0)

coo_quantile_refine_6_x_cell350<-coo_quantile_refine_6_x_cell350[nonNAInd,]
coo_quantile_refine_6_y_cell350<-coo_quantile_refine_6_y_cell350[nonNAInd,]
coo_quantile_refine_6_z_cell350<-coo_quantile_refine_6_z_cell350[nonNAInd,]


save(coo_quantile_refine_6_x_cell350,file="coo_quantile_refine_6_x_cell350.RData")
save(coo_quantile_refine_6_y_cell350,file="coo_quantile_refine_6_y_cell350.RData")
save(coo_quantile_refine_6_z_cell350,file="coo_quantile_refine_6_z_cell350.RData")


### coo list 

coo_list_q6_cell350<-lapply(1:dim(coo_quantile_refine_6_x_cell350)[1],function(i){
	temp<-cbind(coo_quantile_refine_6_x_cell350[i,],coo_quantile_refine_6_y_cell350[i,],coo_quantile_refine_6_z_cell350[i,])
	rownames(temp)<-colnames(coo_quantile_refine_6_x_cell350)
	colnames(temp)<-c("X","Y","Z")
	
	return(temp)
	
})

names(coo_list_q6_cell350)<-rownames(coo_quantile_refine_6_x_cell350)


save(coo_list_q6_cell350,file="coo_list_q6_cell350.RData")



################################################################################

### stage - 7


coo_quantile_refine_7_x<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_7_x)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_7_x)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_7_y<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_7_y)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_7_y)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_7_z<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_7_z)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_7_z)<-names(coo_quantile_refine[[1]])

for(i in 1:length(coo_quantile_refine)){
	print(i)
	for(j in 1:length(coo_quantile_refine[[1]])){
		coo_quantile_refine_7_x[i,j]<-coo_quantile_refine[[i]][[j]][7,1]
		coo_quantile_refine_7_y[i,j]<-coo_quantile_refine[[i]][[j]][7,2]
		coo_quantile_refine_7_z[i,j]<-coo_quantile_refine[[i]][[j]][7,3]
		
	}
}

save(coo_quantile_refine_7_x,file="coo_quantile_refine_7_x.RData")
save(coo_quantile_refine_7_y,file="coo_quantile_refine_7_y.RData")
save(coo_quantile_refine_7_z,file="coo_quantile_refine_7_z.RData")


### Cell350

cell_select<-read.csv(paste0(focal_projectFolderPath,"/data/basic/coordinates_wt.csv"))[,1]

coo_quantile_refine_7_x_cell350<-coo_quantile_refine_7_x[,match(cell_select,colnames(coo_quantile_refine_7_x))]
coo_quantile_refine_7_y_cell350<-coo_quantile_refine_7_y[,match(cell_select,colnames(coo_quantile_refine_7_y))]
coo_quantile_refine_7_z_cell350<-coo_quantile_refine_7_z[,match(cell_select,colnames(coo_quantile_refine_7_z))]

nonNAInd<-which(apply(cbind(coo_quantile_refine_7_x_cell350,coo_quantile_refine_7_y_cell350,coo_quantile_refine_7_z_cell350),1,function(x){length(which(is.na(x)))})==0)

coo_quantile_refine_7_x_cell350<-coo_quantile_refine_7_x_cell350[nonNAInd,]
coo_quantile_refine_7_y_cell350<-coo_quantile_refine_7_y_cell350[nonNAInd,]
coo_quantile_refine_7_z_cell350<-coo_quantile_refine_7_z_cell350[nonNAInd,]


save(coo_quantile_refine_7_x_cell350,file="coo_quantile_refine_7_x_cell350.RData")
save(coo_quantile_refine_7_y_cell350,file="coo_quantile_refine_7_y_cell350.RData")
save(coo_quantile_refine_7_z_cell350,file="coo_quantile_refine_7_z_cell350.RData")


### coo list 

coo_list_q7_cell350<-lapply(1:dim(coo_quantile_refine_7_x_cell350)[1],function(i){
	temp<-cbind(coo_quantile_refine_7_x_cell350[i,],coo_quantile_refine_7_y_cell350[i,],coo_quantile_refine_7_z_cell350[i,])
	rownames(temp)<-colnames(coo_quantile_refine_7_x_cell350)
	colnames(temp)<-c("X","Y","Z")
	
	return(temp)
	
})

names(coo_list_q7_cell350)<-rownames(coo_quantile_refine_7_x_cell350)


save(coo_list_q7_cell350,file="coo_list_q7_cell350.RData")



################################################################################

### stage - 8


coo_quantile_refine_8_x<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_8_x)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_8_x)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_8_y<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_8_y)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_8_y)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_8_z<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_8_z)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_8_z)<-names(coo_quantile_refine[[1]])

for(i in 1:length(coo_quantile_refine)){
	print(i)
	for(j in 1:length(coo_quantile_refine[[1]])){
		coo_quantile_refine_8_x[i,j]<-coo_quantile_refine[[i]][[j]][8,1]
		coo_quantile_refine_8_y[i,j]<-coo_quantile_refine[[i]][[j]][8,2]
		coo_quantile_refine_8_z[i,j]<-coo_quantile_refine[[i]][[j]][8,3]
		
	}
}

save(coo_quantile_refine_8_x,file="coo_quantile_refine_8_x.RData")
save(coo_quantile_refine_8_y,file="coo_quantile_refine_8_y.RData")
save(coo_quantile_refine_8_z,file="coo_quantile_refine_8_z.RData")


### Cell350

cell_select<-read.csv(paste0(focal_projectFolderPath,"/data/basic/coordinates_wt.csv"))[,1]

coo_quantile_refine_8_x_cell350<-coo_quantile_refine_8_x[,match(cell_select,colnames(coo_quantile_refine_8_x))]
coo_quantile_refine_8_y_cell350<-coo_quantile_refine_8_y[,match(cell_select,colnames(coo_quantile_refine_8_y))]
coo_quantile_refine_8_z_cell350<-coo_quantile_refine_8_z[,match(cell_select,colnames(coo_quantile_refine_8_z))]

nonNAInd<-which(apply(cbind(coo_quantile_refine_8_x_cell350,coo_quantile_refine_8_y_cell350,coo_quantile_refine_8_z_cell350),1,function(x){length(which(is.na(x)))})==0)

coo_quantile_refine_8_x_cell350<-coo_quantile_refine_8_x_cell350[nonNAInd,]
coo_quantile_refine_8_y_cell350<-coo_quantile_refine_8_y_cell350[nonNAInd,]
coo_quantile_refine_8_z_cell350<-coo_quantile_refine_8_z_cell350[nonNAInd,]


save(coo_quantile_refine_8_x_cell350,file="coo_quantile_refine_8_x_cell350.RData")
save(coo_quantile_refine_8_y_cell350,file="coo_quantile_refine_8_y_cell350.RData")
save(coo_quantile_refine_8_z_cell350,file="coo_quantile_refine_8_z_cell350.RData")


### coo list 

coo_list_q8_cell350<-lapply(1:dim(coo_quantile_refine_8_x_cell350)[1],function(i){
	temp<-cbind(coo_quantile_refine_8_x_cell350[i,],coo_quantile_refine_8_y_cell350[i,],coo_quantile_refine_8_z_cell350[i,])
	rownames(temp)<-colnames(coo_quantile_refine_8_x_cell350)
	colnames(temp)<-c("X","Y","Z")
	
	return(temp)
	
})

names(coo_list_q8_cell350)<-rownames(coo_quantile_refine_8_x_cell350)


save(coo_list_q8_cell350,file="coo_list_q8_cell350.RData")



################################################################################

### stage - 9


coo_quantile_refine_9_x<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_9_x)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_9_x)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_9_y<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_9_y)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_9_y)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_9_z<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_9_z)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_9_z)<-names(coo_quantile_refine[[1]])

for(i in 1:length(coo_quantile_refine)){
	print(i)
	for(j in 1:length(coo_quantile_refine[[1]])){
		coo_quantile_refine_9_x[i,j]<-coo_quantile_refine[[i]][[j]][9,1]
		coo_quantile_refine_9_y[i,j]<-coo_quantile_refine[[i]][[j]][9,2]
		coo_quantile_refine_9_z[i,j]<-coo_quantile_refine[[i]][[j]][9,3]
		
	}
}

save(coo_quantile_refine_9_x,file="coo_quantile_refine_9_x.RData")
save(coo_quantile_refine_9_y,file="coo_quantile_refine_9_y.RData")
save(coo_quantile_refine_9_z,file="coo_quantile_refine_9_z.RData")


### Cell350

cell_select<-read.csv(paste0(focal_projectFolderPath,"/data/basic/coordinates_wt.csv"))[,1]

coo_quantile_refine_9_x_cell350<-coo_quantile_refine_9_x[,match(cell_select,colnames(coo_quantile_refine_9_x))]
coo_quantile_refine_9_y_cell350<-coo_quantile_refine_9_y[,match(cell_select,colnames(coo_quantile_refine_9_y))]
coo_quantile_refine_9_z_cell350<-coo_quantile_refine_9_z[,match(cell_select,colnames(coo_quantile_refine_9_z))]

nonNAInd<-which(apply(cbind(coo_quantile_refine_9_x_cell350,coo_quantile_refine_9_y_cell350,coo_quantile_refine_9_z_cell350),1,function(x){length(which(is.na(x)))})==0)

coo_quantile_refine_9_x_cell350<-coo_quantile_refine_9_x_cell350[nonNAInd,]
coo_quantile_refine_9_y_cell350<-coo_quantile_refine_9_y_cell350[nonNAInd,]
coo_quantile_refine_9_z_cell350<-coo_quantile_refine_9_z_cell350[nonNAInd,]


save(coo_quantile_refine_9_x_cell350,file="coo_quantile_refine_9_x_cell350.RData")
save(coo_quantile_refine_9_y_cell350,file="coo_quantile_refine_9_y_cell350.RData")
save(coo_quantile_refine_9_z_cell350,file="coo_quantile_refine_9_z_cell350.RData")


### coo list 

coo_list_q9_cell350<-lapply(1:dim(coo_quantile_refine_9_x_cell350)[1],function(i){
	temp<-cbind(coo_quantile_refine_9_x_cell350[i,],coo_quantile_refine_9_y_cell350[i,],coo_quantile_refine_9_z_cell350[i,])
	rownames(temp)<-colnames(coo_quantile_refine_9_x_cell350)
	colnames(temp)<-c("X","Y","Z")
	
	return(temp)
	
})

names(coo_list_q9_cell350)<-rownames(coo_quantile_refine_9_x_cell350)


save(coo_list_q9_cell350,file="coo_list_q9_cell350.RData")



################################################################################

### stage - 10


coo_quantile_refine_10_x<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_10_x)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_10_x)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_10_y<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_10_y)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_10_y)<-names(coo_quantile_refine[[1]])

coo_quantile_refine_10_z<-matrix(nrow=length(coo_quantile_refine),ncol=length(coo_quantile_refine[[1]]))
rownames(coo_quantile_refine_10_z)<-names(coo_quantile_refine)
colnames(coo_quantile_refine_10_z)<-names(coo_quantile_refine[[1]])

for(i in 1:length(coo_quantile_refine)){
	print(i)
	for(j in 1:length(coo_quantile_refine[[1]])){
		coo_quantile_refine_10_x[i,j]<-coo_quantile_refine[[i]][[j]][10,1]
		coo_quantile_refine_10_y[i,j]<-coo_quantile_refine[[i]][[j]][10,2]
		coo_quantile_refine_10_z[i,j]<-coo_quantile_refine[[i]][[j]][10,3]
		
	}
}

save(coo_quantile_refine_10_x,file="coo_quantile_refine_10_x.RData")
save(coo_quantile_refine_10_y,file="coo_quantile_refine_10_y.RData")
save(coo_quantile_refine_10_z,file="coo_quantile_refine_10_z.RData")


### Cell350

cell_select<-read.csv(paste0(focal_projectFolderPath,"/data/basic/coordinates_wt.csv"))[,1]

coo_quantile_refine_10_x_cell350<-coo_quantile_refine_10_x[,match(cell_select,colnames(coo_quantile_refine_10_x))]
coo_quantile_refine_10_y_cell350<-coo_quantile_refine_10_y[,match(cell_select,colnames(coo_quantile_refine_10_y))]
coo_quantile_refine_10_z_cell350<-coo_quantile_refine_10_z[,match(cell_select,colnames(coo_quantile_refine_10_z))]

nonNAInd<-which(apply(cbind(coo_quantile_refine_10_x_cell350,coo_quantile_refine_10_y_cell350,coo_quantile_refine_10_z_cell350),1,function(x){length(which(is.na(x)))})==0)

coo_quantile_refine_10_x_cell350<-coo_quantile_refine_10_x_cell350[nonNAInd,]
coo_quantile_refine_10_y_cell350<-coo_quantile_refine_10_y_cell350[nonNAInd,]
coo_quantile_refine_10_z_cell350<-coo_quantile_refine_10_z_cell350[nonNAInd,]


save(coo_quantile_refine_10_x_cell350,file="coo_quantile_refine_10_x_cell350.RData")
save(coo_quantile_refine_10_y_cell350,file="coo_quantile_refine_10_y_cell350.RData")
save(coo_quantile_refine_10_z_cell350,file="coo_quantile_refine_10_z_cell350.RData")


### coo list 

coo_list_q10_cell350<-lapply(1:dim(coo_quantile_refine_10_x_cell350)[1],function(i){
	temp<-cbind(coo_quantile_refine_10_x_cell350[i,],coo_quantile_refine_10_y_cell350[i,],coo_quantile_refine_10_z_cell350[i,])
	rownames(temp)<-colnames(coo_quantile_refine_10_x_cell350)
	colnames(temp)<-c("X","Y","Z")
	
	return(temp)
	
})

names(coo_list_q10_cell350)<-rownames(coo_quantile_refine_10_x_cell350)


save(coo_list_q10_cell350,file="coo_list_q10_cell350.RData")


################################################################################


