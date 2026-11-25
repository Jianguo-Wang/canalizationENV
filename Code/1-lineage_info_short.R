

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


### get lineage of a cell 

#

get_lineage_AB<-function(end_list){
	
	return(lapply(end_list,function(x){
		temp<-unlist(strsplit(x,""))
		lineage_temp<-c("P0","AB",sapply(3:length(temp),function(i){
			paste0(temp[1:i],collapse="")
		}))
	
		return(lineage_temp)

	}))

}

#

get_lineage_MS<-function(end_list){
	
	return(lapply(end_list,function(x){
		temp<-unlist(strsplit(x,""))
		lineage_temp<-c("P0","P1","EMS",sapply(2:length(temp),function(i){
			paste0(temp[1:i],collapse="")
		}))
	
		return(lineage_temp)

	}))

}

#

get_lineage_E<-function(end_list){
	
	return(lapply(end_list,function(x){
		temp<-unlist(strsplit(x,""))
		lineage_temp<-c("P0","P1","EMS","E",sapply(2:length(temp),function(i){
			paste0(temp[1:i],collapse="")
		}))
	
		return(lineage_temp)

	}))

}

#

get_lineage_C<-function(end_list){
	
	return(lapply(end_list,function(x){
		temp<-unlist(strsplit(x,""))
		lineage_temp<-c("P0","P1","P2","C",sapply(2:length(temp),function(i){
			paste0(temp[1:i],collapse="")
		}))
	
		return(lineage_temp)

	}))

}

#

get_lineage_D<-function(end_list){
	
	return(lapply(end_list,function(x){
		temp<-unlist(strsplit(x,""))
		lineage_temp<-c("P0","P1","P2","P3","D",sapply(2:length(temp),function(i){
			paste0(temp[1:i],collapse="")
		}))
	
		return(lineage_temp)

	}))

}

#

get_lineage_Z2<-function(end_list){
	
	return(lapply(end_list,function(x){
		temp<-unlist(strsplit(x,""))
		lineage_temp<-c("P0","P1","P2","P3","P4","Z2",sapply(3:length(temp),function(i){
			paste0(temp[1:i],collapse="")
		}))
	
		return(lineage_temp)

	}))

}

#

get_lineage_Z3<-function(end_list){
	
	return(lapply(end_list,function(x){
		temp<-unlist(strsplit(x,""))
		lineage_temp<-c("P0","P1","P2","P3","P4","Z3",sapply(3:length(temp),function(i){
			paste0(temp[1:i],collapse="")
		}))
	
		return(lineage_temp)

	}))

}




###################################################################################

###################################################################################


focal_projectFolderPath<-r"(E:\canalizationENV)"


##

setwd(paste0(focal_projectFolderPath,"/data/basic"))


##

loadResources(c("E:/canalizationENV/data/interpolation/coo_Control.RData",
				"E:/canalizationENV/data/interpolation/coo_RNAi.RData",
				"E:/canalizationENV/data/interpolation/coo_Env.RData"),projectFolderPath=focal_projectFolderPath)


##

cell714<-read.csv("E:/canalizationENV/data/basic/coordinates_wt.csv")[,1]


save(cell714,file="cell714.RData")

all_cells = cell714

## AB

AB_end_10<-intersect(all_cells[grep("AB",all_cells)],all_cells[nchar(all_cells)==10])
AB_end_10_list<-get_lineage_AB(AB_end_10)

intersect(all_cells[grep("AB",all_cells)],setdiff(all_cells[nchar(all_cells)==9],unlist(AB_end_10_list)))


## MS 

MS_end_8<-intersect(all_cells[grep("MS",all_cells)],all_cells[nchar(all_cells)==8])
MS_end_8_list<-get_lineage_MS(MS_end_8)

MS_end_7<-intersect(all_cells[grep("MS",all_cells)],setdiff(all_cells[nchar(all_cells)==7],unlist(MS_end_8_list)))
MS_end_7_list<-get_lineage_MS(MS_end_7)

intersect(all_cells[grep("MS",all_cells)],setdiff(setdiff(all_cells[nchar(all_cells)==6],unlist(MS_end_7_list)),unlist(MS_end_8_list)))


## E

E_end_5<-intersect(all_cells[grep("E",all_cells)],all_cells[nchar(all_cells)==5])
E_end_5_list<-get_lineage_E(E_end_5)

setdiff(all_cells[grep("E",all_cells)],unlist(E_end_5_list))

## C

C_end_6<-intersect(all_cells[grep("C",all_cells)],all_cells[nchar(all_cells)==6])
C_end_6_list<-get_lineage_C(C_end_6)

C_end_5<-intersect(all_cells[grep("C",all_cells)],setdiff(all_cells[nchar(all_cells)==5],unlist(C_end_6_list)))
C_end_5_list<-get_lineage_C(C_end_5)

setdiff(all_cells[grep("C",all_cells)],union(unlist(C_end_6_list),unlist(C_end_5_list)))

## D

D_end_4<-intersect(all_cells[grep("D",all_cells)],all_cells[nchar(all_cells)==4])
D_end_4_list<-get_lineage_D(D_end_4)

setdiff(all_cells[grep("D",all_cells)],unlist(D_end_4_list))

## Z2 

Z2_end_2_list<-list(c("P0","P1","P2","P3","P4","Z2"))

setdiff(all_cells[grep("Z2",all_cells)],unlist(Z2_end_2_list))

## Z3

Z3_end_2_list<-list(c("P0","P1","P2","P3","P4","Z3"))

setdiff(all_cells[grep("Z3",all_cells)],unlist(Z3_end_2_list))


### generate lineage_info 

lineage_info_list<-c(AB_end_10_list,MS_end_8_list,MS_end_7_list,E_end_5_list,C_end_6_list,C_end_5_list,D_end_4_list,Z2_end_2_list,Z3_end_2_list)

lineage_info<-matrix(nrow=max(sapply(lineage_info_list,length)),ncol=length(lineage_info_list))
for(i in 1:length(lineage_info_list)){
	lineage_info[1:length(lineage_info_list[[i]]),i]<-lineage_info_list[[i]]
}

rownames(lineage_info)<-paste0("g",1:10)
colnames(lineage_info)<-c(AB_end_10,MS_end_8,MS_end_7,E_end_5,C_end_6,C_end_5,D_end_4,"Z2","Z3")

lineage_info[1,]<-"Zygote"

save(lineage_info,file="lineage_info_short.RData")


###############################################################################
                                   end
###############################################################################








