

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


focal_projectFolderPath<-r"(F:\wangj\Desktop\temp\canalizationENV)"

#

setwd(paste0(focal_projectFolderPath,"/figs"))

#

loadResources(c("F:/wangj/Desktop/temp/canalizationENV/functions/set_theme.R",
	"F:/wangj/Desktop/temp/canalizationENV/functions/RowColRecovery.R",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/cell_mo_da1_da2.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/lineage_info_short.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz.RData",
	"F:/wangj/Desktop/temp/canalizationENV/data/basic/q1_cell350_alignment_xyz_noise.RData",	
	"F:/wangj/Desktop/temp/canalizationENV/data/alignment/embryo_size.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_x_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_y_model_mother_coef.RData",
	"E:/papers/canalizationENV/data/spatial/xyz_z_model_mother_coef.RData"
),projectFolderPath=focal_projectFolderPath)



#

q1_cell350_alignment_xyz_noise_x<-q1_cell350_alignment_xyz_noise[,1:714]
q1_cell350_alignment_xyz_noise_y<-q1_cell350_alignment_xyz_noise[,(1+714):(714*2)]
q1_cell350_alignment_xyz_noise_z<-q1_cell350_alignment_xyz_noise[,(1+714*2):(714*3)]


cell_name<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise[,1:714]),"_x"))
cell_gen<-sapply(match(cell_name,lineage_info),function(x){rowColRec(x,dim(lineage_info)[1],dim(lineage_info)[2])[1]})


########################################################################################################

# x

cell_var_x<-apply(q1_cell350_alignment_xyz_noise_x,2,var)
var_x_obs<-cell_var_x
var_x_mean_obs<-sapply(4:10,function(i){mean(cell_var_x[cell_gen==i])})

cell_var_new<-cell_var_x

for(i in 1:714){
	focal_da1<-cell_name[i]
	focal_mo<-cell_mo_da1_da2[which(cell_mo_da1_da2[,2]==focal_da1),1]
	focal_mo_ind<-match(focal_mo,cell_name)
	
	if(!is.na(focal_mo_ind)){
		cell_var_new[i]<-var(residuals(lm(q1_cell350_alignment_xyz_noise_x[,i]~q1_cell350_alignment_xyz_noise_x[,focal_mo_ind])))
	}
	
}

cell_var_expected<-sapply(cell_name,function(x){
	focal_lineage<-lineage_info[1:rowColRec(which(lineage_info==x)[1],dim(lineage_info)[1],dim(lineage_info)[2])[1],rowColRec(which(lineage_info==x)[1],dim(lineage_info)[1],dim(lineage_info)[2])[2]]
	focalInd<-na.omit(match(focal_lineage,cell_name))
	if(length(focalInd)!=0){
		focal_exp<-sum(cell_var_new[focalInd])
	}else{
		focal_exp<-cell_var_new[match(x,cell_name)]
	}
	return(focal_exp)
})

var_x_exp<-cell_var_expected
var_x_mean_exp<-sapply(4:10,function(i){mean(cell_var_expected[cell_gen==i])})


# y

cell_var_y<-apply(q1_cell350_alignment_xyz_noise_y,2,var)
var_y_obs<-cell_var_y
var_y_mean_obs<-sapply(4:10,function(i){mean(cell_var_y[cell_gen==i])})

cell_var_new<-cell_var_y

for(i in 1:714){
	focal_da1<-cell_name[i]
	focal_mo<-cell_mo_da1_da2[which(cell_mo_da1_da2[,2]==focal_da1),1]
	focal_mo_ind<-match(focal_mo,cell_name)
	
	if(!is.na(focal_mo_ind)){
		cell_var_new[i]<-var(residuals(lm(q1_cell350_alignment_xyz_noise_y[,i]~q1_cell350_alignment_xyz_noise_y[,focal_mo_ind])))
	}
	
}

cell_var_expected<-sapply(cell_name,function(x){
	focal_lineage<-lineage_info[1:rowColRec(which(lineage_info==x)[1],dim(lineage_info)[1],dim(lineage_info)[2])[1],rowColRec(which(lineage_info==x)[1],dim(lineage_info)[1],dim(lineage_info)[2])[2]]
	focalInd<-na.omit(match(focal_lineage,cell_name))
	if(length(focalInd)!=0){
		focal_exp<-sum(cell_var_new[focalInd])
	}else{
		focal_exp<-cell_var_new[match(x,cell_name)]
	}
	return(focal_exp)
})

var_y_exp<-cell_var_expected
var_y_mean_exp<-sapply(4:10,function(i){mean(cell_var_expected[cell_gen==i])})



# z

cell_var_z<-apply(q1_cell350_alignment_xyz_noise_z,2,var)
var_z_obs<-cell_var_z
var_z_mean_obs<-c(sapply(4:10,function(i){mean(cell_var_z[cell_gen==i])}))

cell_var_new<-cell_var_z

for(i in 1:714){
	focal_da1<-cell_name[i]
	focal_mo<-cell_mo_da1_da2[which(cell_mo_da1_da2[,2]==focal_da1),1]
	focal_mo_ind<-match(focal_mo,cell_name)
	
	if(!is.na(focal_mo_ind)){
		cell_var_new[i]<-var(residuals(lm(q1_cell350_alignment_xyz_noise_z[,i]~q1_cell350_alignment_xyz_noise_z[,focal_mo_ind])))
	}
	
}

cell_var_expected<-sapply(cell_name,function(x){
	focal_lineage<-lineage_info[1:rowColRec(which(lineage_info==x)[1],dim(lineage_info)[1],dim(lineage_info)[2])[1],rowColRec(which(lineage_info==x)[1],dim(lineage_info)[1],dim(lineage_info)[2])[2]]
	focalInd<-na.omit(match(focal_lineage,cell_name))
	if(length(focalInd)!=0){
		focal_exp<-sum(cell_var_new[focalInd])
	}else{
		focal_exp<-cell_var_new[match(x,cell_name)]
	}
	return(focal_exp)
})

var_z_exp<-cell_var_expected
var_z_mean_exp<-sapply(4:10,function(i){mean(cell_var_expected[cell_gen==i])})


#

cell_var_exp_obs_x<-data.frame(cell=cell_name,cell_var_exp=var_x_exp,cell_var_obs=var_x_obs,generation=cell_gen)
cell_var_exp_obs_y<-data.frame(cell=cell_name,cell_var_exp=var_y_exp,cell_var_obs=var_y_obs,generation=cell_gen)
cell_var_exp_obs_z<-data.frame(cell=cell_name,cell_var_exp=var_z_exp,cell_var_obs=var_z_obs,generation=cell_gen)

save(cell_var_exp_obs_x,file=paste0(focal_projectFolderPath,"/data/withoutNegativeFeedback/cell_var_exp_obs_x.RData"))
save(cell_var_exp_obs_y,file=paste0(focal_projectFolderPath,"/data/withoutNegativeFeedback/cell_var_exp_obs_y.RData"))
save(cell_var_exp_obs_z,file=paste0(focal_projectFolderPath,"/data/withoutNegativeFeedback/cell_var_exp_obs_z.RData"))


cell_var_mean_exp_obs_x<-data.frame(loc=rep(1,14),generation=rep(c(10:4),2),cell_var_mean=c(var_x_mean_exp,var_x_mean_obs),type=c(rep("exp",7),rep("obs",7)))
cell_var_mean_exp_obs_y<-data.frame(loc=rep(1,14),generation=rep(c(10:4),2),cell_var_mean=c(var_y_mean_exp,var_y_mean_obs),type=c(rep("exp",7),rep("obs",7)))
cell_var_mean_exp_obs_z<-data.frame(loc=rep(1,14),generation=rep(c(10:4),2),cell_var_mean=c(var_z_mean_exp,var_z_mean_obs),type=c(rep("exp",7),rep("obs",7)))


save(cell_var_mean_exp_obs_x,file=paste0(focal_projectFolderPath,"/data/withoutNegativeFeedback/cell_var_mean_exp_obs_x.RData"))
save(cell_var_mean_exp_obs_y,file=paste0(focal_projectFolderPath,"/data/withoutNegativeFeedback/cell_var_mean_exp_obs_y.RData"))
save(cell_var_mean_exp_obs_z,file=paste0(focal_projectFolderPath,"/data/withoutNegativeFeedback/cell_var_mean_exp_obs_z.RData"))


##################################################################################################################3

#x

ggplot(data=cell_var_mean_exp_obs_x,aes(x=loc,y=generation))+
geom_point(size=sqrt(cell_var_mean_exp_obs_x$cell_var_mean)*0.2,stroke=0.5,pch=c(rep(1,7),rep(16,7)),col=c(rep("#00BFC4",7),rep("#F8766D",7)))+
labs(title="")

var_x_mean_exp/var_x_mean_obs

ggsave("withoutNF_variance_accumulation_x_average.pdf",width=unit(2.3,"in"),height=unit(2.4,"in"))

#y 

ggplot(data=cell_var_mean_exp_obs_y,aes(x=loc,y=generation))+
geom_point(size=sqrt(cell_var_mean_exp_obs_y$cell_var_mean)*0.2,stroke=0.5,pch=c(rep(1,7),rep(16,7)),col=c(rep("#00BFC4",7),rep("#F8766D",7)))+
labs(title="")

var_y_mean_exp/var_y_mean_obs

ggsave("withoutNF_variance_accumulation_y_average.pdf",width=unit(2.3,"in"),height=unit(2.4,"in"))

#z

ggplot(data=cell_var_mean_exp_obs_z,aes(x=loc,y=generation))+
geom_point(size=sqrt(cell_var_mean_exp_obs_z$cell_var_mean)*0.2,stroke=0.5,pch=c(rep(1,7),rep(16,7)),col=c(rep("#00BFC4",7),rep("#F8766D",7)))+
labs(title="")

var_z_mean_exp/var_z_mean_obs

ggsave("withoutNF_variance_accumulation_z_average.pdf",width=unit(2.3,"in"),height=unit(2.4,"in"))


####################################################################33


triangle_tree_10<-data.frame(x=c(1:2^9),y=rep(1,2^9))
triangle_tree_9<-data.frame(x=(triangle_tree_10[c(1:2^8)*2-1,1]+triangle_tree_10[c(1:2^8)*2,1])/2,y=rep(2,2^8))
triangle_tree_8<-data.frame(x=(triangle_tree_9[c(1:2^7)*2-1,1]+triangle_tree_9[c(1:2^7)*2,1])/2,y=rep(3,2^7))
triangle_tree_7<-data.frame(x=(triangle_tree_8[c(1:2^6)*2-1,1]+triangle_tree_8[c(1:2^6)*2,1])/2,y=rep(4,2^6))
triangle_tree_6<-data.frame(x=(triangle_tree_7[c(1:2^5)*2-1,1]+triangle_tree_7[c(1:2^5)*2,1])/2,y=rep(5,2^5))
triangle_tree_5<-data.frame(x=(triangle_tree_6[c(1:2^4)*2-1,1]+triangle_tree_6[c(1:2^4)*2,1])/2,y=rep(6,2^4))
triangle_tree_4<-data.frame(x=(triangle_tree_5[c(1:2^3)*2-1,1]+triangle_tree_5[c(1:2^3)*2,1])/2,y=rep(7,2^3))
triangle_tree_3<-data.frame(x=(triangle_tree_4[c(1:2^2)*2-1,1]+triangle_tree_4[c(1:2^2)*2,1])/2,y=rep(8,2^2))
triangle_tree_2<-data.frame(x=(triangle_tree_3[c(1:2^1)*2-1,1]+triangle_tree_3[c(1:2^1)*2,1])/2,y=rep(9,2^1))
triangle_tree_1<-data.frame(x=(triangle_tree_2[c(1:2^0)*2-1,1]+triangle_tree_2[c(1:2^0)*2,1])/2,y=rep(10,2^0))


triangle_tree_1$barcode="0"
triangle_tree_2$barcode=c("00","01")
triangle_tree_3$barcode=c(sapply(triangle_tree_2$barcode,function(x){paste0(x,c("0","1"))}))
triangle_tree_4$barcode=c(sapply(triangle_tree_3$barcode,function(x){paste0(x,c("0","1"))}))
triangle_tree_5$barcode=c(sapply(triangle_tree_4$barcode,function(x){paste0(x,c("0","1"))}))
triangle_tree_6$barcode=c(sapply(triangle_tree_5$barcode,function(x){paste0(x,c("0","1"))}))
triangle_tree_7$barcode=c(sapply(triangle_tree_6$barcode,function(x){paste0(x,c("0","1"))}))
triangle_tree_8$barcode=c(sapply(triangle_tree_7$barcode,function(x){paste0(x,c("0","1"))}))
triangle_tree_9$barcode=c(sapply(triangle_tree_8$barcode,function(x){paste0(x,c("0","1"))}))
triangle_tree_10$barcode=c(sapply(triangle_tree_9$barcode,function(x){paste0(x,c("0","1"))}))



triangle_tree<-rbind(triangle_tree_1,
	triangle_tree_2,
	triangle_tree_3,
	triangle_tree_4,
	triangle_tree_5,
	triangle_tree_6,
	triangle_tree_7,
	triangle_tree_8,
	triangle_tree_9,
	triangle_tree_10)


mo_da_link<-rbind(data.frame(mo=c(sapply(triangle_tree_1$barcode,rep,2)),da=triangle_tree_2$barcode),
				  data.frame(mo=c(sapply(triangle_tree_2$barcode,rep,2)),da=triangle_tree_3$barcode),
				  data.frame(mo=c(sapply(triangle_tree_3$barcode,rep,2)),da=triangle_tree_4$barcode),
				  data.frame(mo=c(sapply(triangle_tree_4$barcode,rep,2)),da=triangle_tree_5$barcode),
				  data.frame(mo=c(sapply(triangle_tree_5$barcode,rep,2)),da=triangle_tree_6$barcode),
				  data.frame(mo=c(sapply(triangle_tree_6$barcode,rep,2)),da=triangle_tree_7$barcode),
				  data.frame(mo=c(sapply(triangle_tree_7$barcode,rep,2)),da=triangle_tree_8$barcode),
				  data.frame(mo=c(sapply(triangle_tree_8$barcode,rep,2)),da=triangle_tree_9$barcode),
				  data.frame(mo=c(sapply(triangle_tree_9$barcode,rep,2)),da=triangle_tree_10$barcode))
				
mo_da_link$x_mo<-triangle_tree[match(mo_da_link$mo,triangle_tree$barcode),1]
mo_da_link$y_mo<-triangle_tree[match(mo_da_link$mo,triangle_tree$barcode),2]
mo_da_link$x_da<-triangle_tree[match(mo_da_link$da,triangle_tree$barcode),1]
mo_da_link$y_da<-triangle_tree[match(mo_da_link$da,triangle_tree$barcode),2]

	
#

all_cell<-unique(na.omit(c(lineage_info)))

all_cell_barcode_1<-data.frame(cell=unique(na.omit(c(lineage_info[1,]))),barcode="0")
all_cell_barcode_2<-data.frame(cell=unique(na.omit(c(lineage_info[2,]))),barcode=c("00","01"))

all_cell_barcode_3<-data.frame(cell=c(sapply(all_cell_barcode_2$cell,function(x){unique(lineage_info[3,lineage_info[2,]==x])})),
	barcode=c(sapply(all_cell_barcode_2$cell,function(x){
		temp<-length(na.omit(unique(lineage_info[3,lineage_info[2,]==x])));
		if(temp==2){
			return(paste0(all_cell_barcode_2$barcode[match(x,all_cell_barcode_2$cell)],c("0","1")))
		}
		if(temp==1){
			return(paste0(all_cell_barcode_2$barcode[match(x,all_cell_barcode_2$cell)],"0"))
		}
		if(temp==0){
			return(NA)
		}			
	})))


all_cell_barcode_4<-data.frame(cell=c(sapply(all_cell_barcode_3$cell,function(x){unique(lineage_info[4,lineage_info[3,]==x])})),
	barcode=c(sapply(all_cell_barcode_3$cell,function(x){
		temp<-length(na.omit(unique(lineage_info[4,lineage_info[3,]==x])));
		if(temp==2){
			return(paste0(all_cell_barcode_3$barcode[match(x,all_cell_barcode_3$cell)],c("0","1")))
		}
		if(temp==1){
			return(paste0(all_cell_barcode_3$barcode[match(x,all_cell_barcode_3$cell)],"0"))
		}
		if(temp==0){
			return(NA)
		}			
	})))

	
all_cell_barcode_5<-data.frame(cell=c(sapply(all_cell_barcode_4$cell,function(x){unique(lineage_info[5,lineage_info[4,]==x])})),
	barcode=c(sapply(all_cell_barcode_4$cell,function(x){
		temp<-length(na.omit(unique(lineage_info[5,lineage_info[4,]==x])));
		if(temp==2){
			return(paste0(all_cell_barcode_4$barcode[match(x,all_cell_barcode_4$cell)],c("0","1")))
		}
		if(temp==1){
			return(paste0(all_cell_barcode_4$barcode[match(x,all_cell_barcode_4$cell)],"0"))
		}
		if(temp==0){
			return(NA)
		}			
	})))


	
all_cell_barcode_6<-data.frame(cell=unlist(lapply(all_cell_barcode_5$cell,function(x){na.omit(unique(lineage_info[6,lineage_info[5,]==x]))})),
	barcode=unlist(lapply(all_cell_barcode_5$cell,function(x){
		temp<-length(na.omit(unique(lineage_info[6,lineage_info[5,]==x])));
		if(temp==2){
			return(paste0(all_cell_barcode_5$barcode[match(x,all_cell_barcode_5$cell)],c("0","1")))
		}
		if(temp==1){
			return(paste0(all_cell_barcode_5$barcode[match(x,all_cell_barcode_5$cell)],"0"))
		}
		if(temp==0){
			return(NA)
		}			
	})))

	
all_cell_barcode_7<-data.frame(cell=unlist(lapply(all_cell_barcode_6$cell,function(x){na.omit(unique(lineage_info[7,lineage_info[6,]==x]))})),
	barcode=unlist(lapply(all_cell_barcode_6$cell,function(x){
		temp<-length(na.omit(unique(lineage_info[7,lineage_info[6,]==x])));
		if(temp==2){
			return(paste0(all_cell_barcode_6$barcode[match(x,all_cell_barcode_6$cell)],c("0","1")))
		}
		if(temp==1){
			return(paste0(all_cell_barcode_6$barcode[match(x,all_cell_barcode_6$cell)],"0"))
		}
		if(temp==0){
			return(NULL)
		}			
	})))



all_cell_barcode_8<-data.frame(cell=unlist(lapply(all_cell_barcode_7$cell,function(x){na.omit(unique(lineage_info[8,lineage_info[7,]==x]))})),
	barcode=unlist(lapply(all_cell_barcode_7$cell,function(x){
		temp<-length(na.omit(unique(lineage_info[8,lineage_info[7,]==x])));
		if(temp==2){
			return(paste0(all_cell_barcode_7$barcode[match(x,all_cell_barcode_7$cell)],c("0","1")))
		}
		if(temp==1){
			return(paste0(all_cell_barcode_7$barcode[match(x,all_cell_barcode_7$cell)],"0"))
		}
		if(temp==0){
			return(NULL)
		}			
	})))


all_cell_barcode_9<-data.frame(cell=unlist(lapply(all_cell_barcode_8$cell,function(x){na.omit(unique(lineage_info[9,lineage_info[8,]==x]))})),
	barcode=unlist(lapply(all_cell_barcode_8$cell,function(x){
		temp<-length(na.omit(unique(lineage_info[9,lineage_info[8,]==x])));
		if(temp==2){
			return(paste0(all_cell_barcode_8$barcode[match(x,all_cell_barcode_8$cell)],c("0","1")))
		}
		if(temp==1){
			return(paste0(all_cell_barcode_8$barcode[match(x,all_cell_barcode_8$cell)],"0"))
		}
		if(temp==0){
			return(NULL)
		}			
	})))


all_cell_barcode_10<-data.frame(cell=unlist(lapply(all_cell_barcode_9$cell,function(x){na.omit(unique(lineage_info[10,lineage_info[9,]==x]))})),
	barcode=unlist(lapply(all_cell_barcode_9$cell,function(x){
		temp<-length(na.omit(unique(lineage_info[10,lineage_info[9,]==x])));
		if(temp==2){
			return(paste0(all_cell_barcode_9$barcode[match(x,all_cell_barcode_9$cell)],c("0","1")))
		}
		if(temp==1){
			return(paste0(all_cell_barcode_9$barcode[match(x,all_cell_barcode_9$cell)],"0"))
		}
		if(temp==0){
			return(NULL)
		}			
	})))


all_cell_barcode<-rbind(all_cell_barcode_1,
						all_cell_barcode_2,
						all_cell_barcode_3,
						all_cell_barcode_4,
						all_cell_barcode_5,
						all_cell_barcode_6,
						all_cell_barcode_7,
						all_cell_barcode_8,
						all_cell_barcode_9,
						all_cell_barcode_10)



					
triangle_tree_cut<-triangle_tree[match(all_cell_barcode$barcode,triangle_tree$barcode),]					
triangle_tree_cut$cell<-all_cell_barcode$cell
					
mo_da_link_cut<-mo_da_link[na.omit(sapply(1:dim(mo_da_link)[1],function(i){ifelse(length(intersect(mo_da_link[i,1:2],all_cell_barcode$barcode))==2,i,NA)})),]

save(mo_da_link_cut,file=paste0(focal_projectFolderPath,"/functions/mo_da_link_cut.RData"))					
save(triangle_tree_cut,file=paste0(focal_projectFolderPath,"/functions/triangle_tree_cut.RData"))					


#

mo_da_link_cut_script<-sapply(1:dim(mo_da_link_cut)[1],function(i){paste0("geom_line(data=data.frame(x=c(mo_da_link_cut[",i,",3],mo_da_link_cut[",i,",5]),y=c(mo_da_link_cut[",i,",4],mo_da_link_cut[",i,",6])),aes(x=x,y=y),col=\"grey\"",",lwd=0.1)+")})
write.table(mo_da_link_cut_script,paste0(focal_projectFolderPath,"/data/withoutNegativeFeedback/mo_da_link_cut_script.txt"),quote=FALSE,row.names=F,col.names=F)




plot_triangle_tree<-function(focal_triangle_tree_cut,focal_cells){

ggplot(data=focal_triangle_tree_cut,aes(x=x,y=y))+geom_point(size=0.02,stroke=0.1,col="grey")+
geom_line(data=data.frame(x=c(mo_da_link_cut[1,3],mo_da_link_cut[1,5]),y=c(mo_da_link_cut[1,4],mo_da_link_cut[1,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[2,3],mo_da_link_cut[2,5]),y=c(mo_da_link_cut[2,4],mo_da_link_cut[2,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[3,3],mo_da_link_cut[3,5]),y=c(mo_da_link_cut[3,4],mo_da_link_cut[3,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[4,3],mo_da_link_cut[4,5]),y=c(mo_da_link_cut[4,4],mo_da_link_cut[4,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[5,3],mo_da_link_cut[5,5]),y=c(mo_da_link_cut[5,4],mo_da_link_cut[5,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[6,3],mo_da_link_cut[6,5]),y=c(mo_da_link_cut[6,4],mo_da_link_cut[6,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[7,3],mo_da_link_cut[7,5]),y=c(mo_da_link_cut[7,4],mo_da_link_cut[7,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[8,3],mo_da_link_cut[8,5]),y=c(mo_da_link_cut[8,4],mo_da_link_cut[8,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[9,3],mo_da_link_cut[9,5]),y=c(mo_da_link_cut[9,4],mo_da_link_cut[9,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[10,3],mo_da_link_cut[10,5]),y=c(mo_da_link_cut[10,4],mo_da_link_cut[10,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[11,3],mo_da_link_cut[11,5]),y=c(mo_da_link_cut[11,4],mo_da_link_cut[11,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[12,3],mo_da_link_cut[12,5]),y=c(mo_da_link_cut[12,4],mo_da_link_cut[12,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[13,3],mo_da_link_cut[13,5]),y=c(mo_da_link_cut[13,4],mo_da_link_cut[13,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[14,3],mo_da_link_cut[14,5]),y=c(mo_da_link_cut[14,4],mo_da_link_cut[14,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[15,3],mo_da_link_cut[15,5]),y=c(mo_da_link_cut[15,4],mo_da_link_cut[15,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[16,3],mo_da_link_cut[16,5]),y=c(mo_da_link_cut[16,4],mo_da_link_cut[16,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[17,3],mo_da_link_cut[17,5]),y=c(mo_da_link_cut[17,4],mo_da_link_cut[17,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[18,3],mo_da_link_cut[18,5]),y=c(mo_da_link_cut[18,4],mo_da_link_cut[18,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[19,3],mo_da_link_cut[19,5]),y=c(mo_da_link_cut[19,4],mo_da_link_cut[19,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[20,3],mo_da_link_cut[20,5]),y=c(mo_da_link_cut[20,4],mo_da_link_cut[20,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[21,3],mo_da_link_cut[21,5]),y=c(mo_da_link_cut[21,4],mo_da_link_cut[21,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[22,3],mo_da_link_cut[22,5]),y=c(mo_da_link_cut[22,4],mo_da_link_cut[22,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[23,3],mo_da_link_cut[23,5]),y=c(mo_da_link_cut[23,4],mo_da_link_cut[23,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[24,3],mo_da_link_cut[24,5]),y=c(mo_da_link_cut[24,4],mo_da_link_cut[24,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[25,3],mo_da_link_cut[25,5]),y=c(mo_da_link_cut[25,4],mo_da_link_cut[25,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[26,3],mo_da_link_cut[26,5]),y=c(mo_da_link_cut[26,4],mo_da_link_cut[26,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[27,3],mo_da_link_cut[27,5]),y=c(mo_da_link_cut[27,4],mo_da_link_cut[27,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[28,3],mo_da_link_cut[28,5]),y=c(mo_da_link_cut[28,4],mo_da_link_cut[28,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[29,3],mo_da_link_cut[29,5]),y=c(mo_da_link_cut[29,4],mo_da_link_cut[29,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[30,3],mo_da_link_cut[30,5]),y=c(mo_da_link_cut[30,4],mo_da_link_cut[30,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[31,3],mo_da_link_cut[31,5]),y=c(mo_da_link_cut[31,4],mo_da_link_cut[31,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[32,3],mo_da_link_cut[32,5]),y=c(mo_da_link_cut[32,4],mo_da_link_cut[32,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[33,3],mo_da_link_cut[33,5]),y=c(mo_da_link_cut[33,4],mo_da_link_cut[33,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[34,3],mo_da_link_cut[34,5]),y=c(mo_da_link_cut[34,4],mo_da_link_cut[34,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[35,3],mo_da_link_cut[35,5]),y=c(mo_da_link_cut[35,4],mo_da_link_cut[35,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[36,3],mo_da_link_cut[36,5]),y=c(mo_da_link_cut[36,4],mo_da_link_cut[36,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[37,3],mo_da_link_cut[37,5]),y=c(mo_da_link_cut[37,4],mo_da_link_cut[37,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[38,3],mo_da_link_cut[38,5]),y=c(mo_da_link_cut[38,4],mo_da_link_cut[38,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[39,3],mo_da_link_cut[39,5]),y=c(mo_da_link_cut[39,4],mo_da_link_cut[39,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[40,3],mo_da_link_cut[40,5]),y=c(mo_da_link_cut[40,4],mo_da_link_cut[40,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[41,3],mo_da_link_cut[41,5]),y=c(mo_da_link_cut[41,4],mo_da_link_cut[41,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[42,3],mo_da_link_cut[42,5]),y=c(mo_da_link_cut[42,4],mo_da_link_cut[42,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[43,3],mo_da_link_cut[43,5]),y=c(mo_da_link_cut[43,4],mo_da_link_cut[43,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[44,3],mo_da_link_cut[44,5]),y=c(mo_da_link_cut[44,4],mo_da_link_cut[44,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[45,3],mo_da_link_cut[45,5]),y=c(mo_da_link_cut[45,4],mo_da_link_cut[45,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[46,3],mo_da_link_cut[46,5]),y=c(mo_da_link_cut[46,4],mo_da_link_cut[46,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[47,3],mo_da_link_cut[47,5]),y=c(mo_da_link_cut[47,4],mo_da_link_cut[47,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[48,3],mo_da_link_cut[48,5]),y=c(mo_da_link_cut[48,4],mo_da_link_cut[48,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[49,3],mo_da_link_cut[49,5]),y=c(mo_da_link_cut[49,4],mo_da_link_cut[49,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[50,3],mo_da_link_cut[50,5]),y=c(mo_da_link_cut[50,4],mo_da_link_cut[50,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[51,3],mo_da_link_cut[51,5]),y=c(mo_da_link_cut[51,4],mo_da_link_cut[51,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[52,3],mo_da_link_cut[52,5]),y=c(mo_da_link_cut[52,4],mo_da_link_cut[52,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[53,3],mo_da_link_cut[53,5]),y=c(mo_da_link_cut[53,4],mo_da_link_cut[53,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[54,3],mo_da_link_cut[54,5]),y=c(mo_da_link_cut[54,4],mo_da_link_cut[54,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[55,3],mo_da_link_cut[55,5]),y=c(mo_da_link_cut[55,4],mo_da_link_cut[55,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[56,3],mo_da_link_cut[56,5]),y=c(mo_da_link_cut[56,4],mo_da_link_cut[56,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[57,3],mo_da_link_cut[57,5]),y=c(mo_da_link_cut[57,4],mo_da_link_cut[57,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[58,3],mo_da_link_cut[58,5]),y=c(mo_da_link_cut[58,4],mo_da_link_cut[58,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[59,3],mo_da_link_cut[59,5]),y=c(mo_da_link_cut[59,4],mo_da_link_cut[59,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[60,3],mo_da_link_cut[60,5]),y=c(mo_da_link_cut[60,4],mo_da_link_cut[60,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[61,3],mo_da_link_cut[61,5]),y=c(mo_da_link_cut[61,4],mo_da_link_cut[61,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[62,3],mo_da_link_cut[62,5]),y=c(mo_da_link_cut[62,4],mo_da_link_cut[62,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[63,3],mo_da_link_cut[63,5]),y=c(mo_da_link_cut[63,4],mo_da_link_cut[63,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[64,3],mo_da_link_cut[64,5]),y=c(mo_da_link_cut[64,4],mo_da_link_cut[64,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[65,3],mo_da_link_cut[65,5]),y=c(mo_da_link_cut[65,4],mo_da_link_cut[65,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[66,3],mo_da_link_cut[66,5]),y=c(mo_da_link_cut[66,4],mo_da_link_cut[66,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[67,3],mo_da_link_cut[67,5]),y=c(mo_da_link_cut[67,4],mo_da_link_cut[67,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[68,3],mo_da_link_cut[68,5]),y=c(mo_da_link_cut[68,4],mo_da_link_cut[68,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[69,3],mo_da_link_cut[69,5]),y=c(mo_da_link_cut[69,4],mo_da_link_cut[69,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[70,3],mo_da_link_cut[70,5]),y=c(mo_da_link_cut[70,4],mo_da_link_cut[70,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[71,3],mo_da_link_cut[71,5]),y=c(mo_da_link_cut[71,4],mo_da_link_cut[71,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[72,3],mo_da_link_cut[72,5]),y=c(mo_da_link_cut[72,4],mo_da_link_cut[72,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[73,3],mo_da_link_cut[73,5]),y=c(mo_da_link_cut[73,4],mo_da_link_cut[73,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[74,3],mo_da_link_cut[74,5]),y=c(mo_da_link_cut[74,4],mo_da_link_cut[74,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[75,3],mo_da_link_cut[75,5]),y=c(mo_da_link_cut[75,4],mo_da_link_cut[75,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[76,3],mo_da_link_cut[76,5]),y=c(mo_da_link_cut[76,4],mo_da_link_cut[76,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[77,3],mo_da_link_cut[77,5]),y=c(mo_da_link_cut[77,4],mo_da_link_cut[77,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[78,3],mo_da_link_cut[78,5]),y=c(mo_da_link_cut[78,4],mo_da_link_cut[78,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[79,3],mo_da_link_cut[79,5]),y=c(mo_da_link_cut[79,4],mo_da_link_cut[79,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[80,3],mo_da_link_cut[80,5]),y=c(mo_da_link_cut[80,4],mo_da_link_cut[80,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[81,3],mo_da_link_cut[81,5]),y=c(mo_da_link_cut[81,4],mo_da_link_cut[81,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[82,3],mo_da_link_cut[82,5]),y=c(mo_da_link_cut[82,4],mo_da_link_cut[82,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[83,3],mo_da_link_cut[83,5]),y=c(mo_da_link_cut[83,4],mo_da_link_cut[83,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[84,3],mo_da_link_cut[84,5]),y=c(mo_da_link_cut[84,4],mo_da_link_cut[84,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[85,3],mo_da_link_cut[85,5]),y=c(mo_da_link_cut[85,4],mo_da_link_cut[85,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[86,3],mo_da_link_cut[86,5]),y=c(mo_da_link_cut[86,4],mo_da_link_cut[86,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[87,3],mo_da_link_cut[87,5]),y=c(mo_da_link_cut[87,4],mo_da_link_cut[87,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[88,3],mo_da_link_cut[88,5]),y=c(mo_da_link_cut[88,4],mo_da_link_cut[88,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[89,3],mo_da_link_cut[89,5]),y=c(mo_da_link_cut[89,4],mo_da_link_cut[89,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[90,3],mo_da_link_cut[90,5]),y=c(mo_da_link_cut[90,4],mo_da_link_cut[90,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[91,3],mo_da_link_cut[91,5]),y=c(mo_da_link_cut[91,4],mo_da_link_cut[91,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[92,3],mo_da_link_cut[92,5]),y=c(mo_da_link_cut[92,4],mo_da_link_cut[92,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[93,3],mo_da_link_cut[93,5]),y=c(mo_da_link_cut[93,4],mo_da_link_cut[93,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[94,3],mo_da_link_cut[94,5]),y=c(mo_da_link_cut[94,4],mo_da_link_cut[94,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[95,3],mo_da_link_cut[95,5]),y=c(mo_da_link_cut[95,4],mo_da_link_cut[95,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[96,3],mo_da_link_cut[96,5]),y=c(mo_da_link_cut[96,4],mo_da_link_cut[96,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[97,3],mo_da_link_cut[97,5]),y=c(mo_da_link_cut[97,4],mo_da_link_cut[97,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[98,3],mo_da_link_cut[98,5]),y=c(mo_da_link_cut[98,4],mo_da_link_cut[98,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[99,3],mo_da_link_cut[99,5]),y=c(mo_da_link_cut[99,4],mo_da_link_cut[99,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[100,3],mo_da_link_cut[100,5]),y=c(mo_da_link_cut[100,4],mo_da_link_cut[100,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[101,3],mo_da_link_cut[101,5]),y=c(mo_da_link_cut[101,4],mo_da_link_cut[101,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[102,3],mo_da_link_cut[102,5]),y=c(mo_da_link_cut[102,4],mo_da_link_cut[102,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[103,3],mo_da_link_cut[103,5]),y=c(mo_da_link_cut[103,4],mo_da_link_cut[103,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[104,3],mo_da_link_cut[104,5]),y=c(mo_da_link_cut[104,4],mo_da_link_cut[104,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[105,3],mo_da_link_cut[105,5]),y=c(mo_da_link_cut[105,4],mo_da_link_cut[105,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[106,3],mo_da_link_cut[106,5]),y=c(mo_da_link_cut[106,4],mo_da_link_cut[106,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[107,3],mo_da_link_cut[107,5]),y=c(mo_da_link_cut[107,4],mo_da_link_cut[107,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[108,3],mo_da_link_cut[108,5]),y=c(mo_da_link_cut[108,4],mo_da_link_cut[108,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[109,3],mo_da_link_cut[109,5]),y=c(mo_da_link_cut[109,4],mo_da_link_cut[109,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[110,3],mo_da_link_cut[110,5]),y=c(mo_da_link_cut[110,4],mo_da_link_cut[110,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[111,3],mo_da_link_cut[111,5]),y=c(mo_da_link_cut[111,4],mo_da_link_cut[111,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[112,3],mo_da_link_cut[112,5]),y=c(mo_da_link_cut[112,4],mo_da_link_cut[112,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[113,3],mo_da_link_cut[113,5]),y=c(mo_da_link_cut[113,4],mo_da_link_cut[113,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[114,3],mo_da_link_cut[114,5]),y=c(mo_da_link_cut[114,4],mo_da_link_cut[114,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[115,3],mo_da_link_cut[115,5]),y=c(mo_da_link_cut[115,4],mo_da_link_cut[115,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[116,3],mo_da_link_cut[116,5]),y=c(mo_da_link_cut[116,4],mo_da_link_cut[116,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[117,3],mo_da_link_cut[117,5]),y=c(mo_da_link_cut[117,4],mo_da_link_cut[117,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[118,3],mo_da_link_cut[118,5]),y=c(mo_da_link_cut[118,4],mo_da_link_cut[118,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[119,3],mo_da_link_cut[119,5]),y=c(mo_da_link_cut[119,4],mo_da_link_cut[119,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[120,3],mo_da_link_cut[120,5]),y=c(mo_da_link_cut[120,4],mo_da_link_cut[120,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[121,3],mo_da_link_cut[121,5]),y=c(mo_da_link_cut[121,4],mo_da_link_cut[121,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[122,3],mo_da_link_cut[122,5]),y=c(mo_da_link_cut[122,4],mo_da_link_cut[122,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[123,3],mo_da_link_cut[123,5]),y=c(mo_da_link_cut[123,4],mo_da_link_cut[123,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[124,3],mo_da_link_cut[124,5]),y=c(mo_da_link_cut[124,4],mo_da_link_cut[124,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[125,3],mo_da_link_cut[125,5]),y=c(mo_da_link_cut[125,4],mo_da_link_cut[125,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[126,3],mo_da_link_cut[126,5]),y=c(mo_da_link_cut[126,4],mo_da_link_cut[126,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[127,3],mo_da_link_cut[127,5]),y=c(mo_da_link_cut[127,4],mo_da_link_cut[127,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[128,3],mo_da_link_cut[128,5]),y=c(mo_da_link_cut[128,4],mo_da_link_cut[128,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[129,3],mo_da_link_cut[129,5]),y=c(mo_da_link_cut[129,4],mo_da_link_cut[129,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[130,3],mo_da_link_cut[130,5]),y=c(mo_da_link_cut[130,4],mo_da_link_cut[130,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[131,3],mo_da_link_cut[131,5]),y=c(mo_da_link_cut[131,4],mo_da_link_cut[131,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[132,3],mo_da_link_cut[132,5]),y=c(mo_da_link_cut[132,4],mo_da_link_cut[132,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[133,3],mo_da_link_cut[133,5]),y=c(mo_da_link_cut[133,4],mo_da_link_cut[133,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[134,3],mo_da_link_cut[134,5]),y=c(mo_da_link_cut[134,4],mo_da_link_cut[134,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[135,3],mo_da_link_cut[135,5]),y=c(mo_da_link_cut[135,4],mo_da_link_cut[135,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[136,3],mo_da_link_cut[136,5]),y=c(mo_da_link_cut[136,4],mo_da_link_cut[136,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[137,3],mo_da_link_cut[137,5]),y=c(mo_da_link_cut[137,4],mo_da_link_cut[137,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[138,3],mo_da_link_cut[138,5]),y=c(mo_da_link_cut[138,4],mo_da_link_cut[138,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[139,3],mo_da_link_cut[139,5]),y=c(mo_da_link_cut[139,4],mo_da_link_cut[139,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[140,3],mo_da_link_cut[140,5]),y=c(mo_da_link_cut[140,4],mo_da_link_cut[140,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[141,3],mo_da_link_cut[141,5]),y=c(mo_da_link_cut[141,4],mo_da_link_cut[141,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[142,3],mo_da_link_cut[142,5]),y=c(mo_da_link_cut[142,4],mo_da_link_cut[142,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[143,3],mo_da_link_cut[143,5]),y=c(mo_da_link_cut[143,4],mo_da_link_cut[143,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[144,3],mo_da_link_cut[144,5]),y=c(mo_da_link_cut[144,4],mo_da_link_cut[144,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[145,3],mo_da_link_cut[145,5]),y=c(mo_da_link_cut[145,4],mo_da_link_cut[145,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[146,3],mo_da_link_cut[146,5]),y=c(mo_da_link_cut[146,4],mo_da_link_cut[146,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[147,3],mo_da_link_cut[147,5]),y=c(mo_da_link_cut[147,4],mo_da_link_cut[147,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[148,3],mo_da_link_cut[148,5]),y=c(mo_da_link_cut[148,4],mo_da_link_cut[148,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[149,3],mo_da_link_cut[149,5]),y=c(mo_da_link_cut[149,4],mo_da_link_cut[149,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[150,3],mo_da_link_cut[150,5]),y=c(mo_da_link_cut[150,4],mo_da_link_cut[150,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[151,3],mo_da_link_cut[151,5]),y=c(mo_da_link_cut[151,4],mo_da_link_cut[151,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[152,3],mo_da_link_cut[152,5]),y=c(mo_da_link_cut[152,4],mo_da_link_cut[152,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[153,3],mo_da_link_cut[153,5]),y=c(mo_da_link_cut[153,4],mo_da_link_cut[153,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[154,3],mo_da_link_cut[154,5]),y=c(mo_da_link_cut[154,4],mo_da_link_cut[154,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[155,3],mo_da_link_cut[155,5]),y=c(mo_da_link_cut[155,4],mo_da_link_cut[155,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[156,3],mo_da_link_cut[156,5]),y=c(mo_da_link_cut[156,4],mo_da_link_cut[156,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[157,3],mo_da_link_cut[157,5]),y=c(mo_da_link_cut[157,4],mo_da_link_cut[157,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[158,3],mo_da_link_cut[158,5]),y=c(mo_da_link_cut[158,4],mo_da_link_cut[158,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[159,3],mo_da_link_cut[159,5]),y=c(mo_da_link_cut[159,4],mo_da_link_cut[159,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[160,3],mo_da_link_cut[160,5]),y=c(mo_da_link_cut[160,4],mo_da_link_cut[160,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[161,3],mo_da_link_cut[161,5]),y=c(mo_da_link_cut[161,4],mo_da_link_cut[161,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[162,3],mo_da_link_cut[162,5]),y=c(mo_da_link_cut[162,4],mo_da_link_cut[162,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[163,3],mo_da_link_cut[163,5]),y=c(mo_da_link_cut[163,4],mo_da_link_cut[163,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[164,3],mo_da_link_cut[164,5]),y=c(mo_da_link_cut[164,4],mo_da_link_cut[164,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[165,3],mo_da_link_cut[165,5]),y=c(mo_da_link_cut[165,4],mo_da_link_cut[165,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[166,3],mo_da_link_cut[166,5]),y=c(mo_da_link_cut[166,4],mo_da_link_cut[166,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[167,3],mo_da_link_cut[167,5]),y=c(mo_da_link_cut[167,4],mo_da_link_cut[167,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[168,3],mo_da_link_cut[168,5]),y=c(mo_da_link_cut[168,4],mo_da_link_cut[168,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[169,3],mo_da_link_cut[169,5]),y=c(mo_da_link_cut[169,4],mo_da_link_cut[169,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[170,3],mo_da_link_cut[170,5]),y=c(mo_da_link_cut[170,4],mo_da_link_cut[170,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[171,3],mo_da_link_cut[171,5]),y=c(mo_da_link_cut[171,4],mo_da_link_cut[171,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[172,3],mo_da_link_cut[172,5]),y=c(mo_da_link_cut[172,4],mo_da_link_cut[172,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[173,3],mo_da_link_cut[173,5]),y=c(mo_da_link_cut[173,4],mo_da_link_cut[173,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[174,3],mo_da_link_cut[174,5]),y=c(mo_da_link_cut[174,4],mo_da_link_cut[174,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[175,3],mo_da_link_cut[175,5]),y=c(mo_da_link_cut[175,4],mo_da_link_cut[175,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[176,3],mo_da_link_cut[176,5]),y=c(mo_da_link_cut[176,4],mo_da_link_cut[176,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[177,3],mo_da_link_cut[177,5]),y=c(mo_da_link_cut[177,4],mo_da_link_cut[177,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[178,3],mo_da_link_cut[178,5]),y=c(mo_da_link_cut[178,4],mo_da_link_cut[178,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[179,3],mo_da_link_cut[179,5]),y=c(mo_da_link_cut[179,4],mo_da_link_cut[179,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[180,3],mo_da_link_cut[180,5]),y=c(mo_da_link_cut[180,4],mo_da_link_cut[180,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[181,3],mo_da_link_cut[181,5]),y=c(mo_da_link_cut[181,4],mo_da_link_cut[181,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[182,3],mo_da_link_cut[182,5]),y=c(mo_da_link_cut[182,4],mo_da_link_cut[182,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[183,3],mo_da_link_cut[183,5]),y=c(mo_da_link_cut[183,4],mo_da_link_cut[183,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[184,3],mo_da_link_cut[184,5]),y=c(mo_da_link_cut[184,4],mo_da_link_cut[184,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[185,3],mo_da_link_cut[185,5]),y=c(mo_da_link_cut[185,4],mo_da_link_cut[185,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[186,3],mo_da_link_cut[186,5]),y=c(mo_da_link_cut[186,4],mo_da_link_cut[186,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[187,3],mo_da_link_cut[187,5]),y=c(mo_da_link_cut[187,4],mo_da_link_cut[187,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[188,3],mo_da_link_cut[188,5]),y=c(mo_da_link_cut[188,4],mo_da_link_cut[188,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[189,3],mo_da_link_cut[189,5]),y=c(mo_da_link_cut[189,4],mo_da_link_cut[189,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[190,3],mo_da_link_cut[190,5]),y=c(mo_da_link_cut[190,4],mo_da_link_cut[190,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[191,3],mo_da_link_cut[191,5]),y=c(mo_da_link_cut[191,4],mo_da_link_cut[191,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[192,3],mo_da_link_cut[192,5]),y=c(mo_da_link_cut[192,4],mo_da_link_cut[192,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[193,3],mo_da_link_cut[193,5]),y=c(mo_da_link_cut[193,4],mo_da_link_cut[193,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[194,3],mo_da_link_cut[194,5]),y=c(mo_da_link_cut[194,4],mo_da_link_cut[194,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[195,3],mo_da_link_cut[195,5]),y=c(mo_da_link_cut[195,4],mo_da_link_cut[195,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[196,3],mo_da_link_cut[196,5]),y=c(mo_da_link_cut[196,4],mo_da_link_cut[196,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[197,3],mo_da_link_cut[197,5]),y=c(mo_da_link_cut[197,4],mo_da_link_cut[197,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[198,3],mo_da_link_cut[198,5]),y=c(mo_da_link_cut[198,4],mo_da_link_cut[198,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[199,3],mo_da_link_cut[199,5]),y=c(mo_da_link_cut[199,4],mo_da_link_cut[199,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[200,3],mo_da_link_cut[200,5]),y=c(mo_da_link_cut[200,4],mo_da_link_cut[200,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[201,3],mo_da_link_cut[201,5]),y=c(mo_da_link_cut[201,4],mo_da_link_cut[201,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[202,3],mo_da_link_cut[202,5]),y=c(mo_da_link_cut[202,4],mo_da_link_cut[202,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[203,3],mo_da_link_cut[203,5]),y=c(mo_da_link_cut[203,4],mo_da_link_cut[203,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[204,3],mo_da_link_cut[204,5]),y=c(mo_da_link_cut[204,4],mo_da_link_cut[204,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[205,3],mo_da_link_cut[205,5]),y=c(mo_da_link_cut[205,4],mo_da_link_cut[205,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[206,3],mo_da_link_cut[206,5]),y=c(mo_da_link_cut[206,4],mo_da_link_cut[206,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[207,3],mo_da_link_cut[207,5]),y=c(mo_da_link_cut[207,4],mo_da_link_cut[207,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[208,3],mo_da_link_cut[208,5]),y=c(mo_da_link_cut[208,4],mo_da_link_cut[208,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[209,3],mo_da_link_cut[209,5]),y=c(mo_da_link_cut[209,4],mo_da_link_cut[209,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[210,3],mo_da_link_cut[210,5]),y=c(mo_da_link_cut[210,4],mo_da_link_cut[210,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[211,3],mo_da_link_cut[211,5]),y=c(mo_da_link_cut[211,4],mo_da_link_cut[211,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[212,3],mo_da_link_cut[212,5]),y=c(mo_da_link_cut[212,4],mo_da_link_cut[212,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[213,3],mo_da_link_cut[213,5]),y=c(mo_da_link_cut[213,4],mo_da_link_cut[213,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[214,3],mo_da_link_cut[214,5]),y=c(mo_da_link_cut[214,4],mo_da_link_cut[214,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[215,3],mo_da_link_cut[215,5]),y=c(mo_da_link_cut[215,4],mo_da_link_cut[215,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[216,3],mo_da_link_cut[216,5]),y=c(mo_da_link_cut[216,4],mo_da_link_cut[216,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[217,3],mo_da_link_cut[217,5]),y=c(mo_da_link_cut[217,4],mo_da_link_cut[217,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[218,3],mo_da_link_cut[218,5]),y=c(mo_da_link_cut[218,4],mo_da_link_cut[218,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[219,3],mo_da_link_cut[219,5]),y=c(mo_da_link_cut[219,4],mo_da_link_cut[219,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[220,3],mo_da_link_cut[220,5]),y=c(mo_da_link_cut[220,4],mo_da_link_cut[220,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[221,3],mo_da_link_cut[221,5]),y=c(mo_da_link_cut[221,4],mo_da_link_cut[221,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[222,3],mo_da_link_cut[222,5]),y=c(mo_da_link_cut[222,4],mo_da_link_cut[222,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[223,3],mo_da_link_cut[223,5]),y=c(mo_da_link_cut[223,4],mo_da_link_cut[223,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[224,3],mo_da_link_cut[224,5]),y=c(mo_da_link_cut[224,4],mo_da_link_cut[224,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[225,3],mo_da_link_cut[225,5]),y=c(mo_da_link_cut[225,4],mo_da_link_cut[225,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[226,3],mo_da_link_cut[226,5]),y=c(mo_da_link_cut[226,4],mo_da_link_cut[226,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[227,3],mo_da_link_cut[227,5]),y=c(mo_da_link_cut[227,4],mo_da_link_cut[227,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[228,3],mo_da_link_cut[228,5]),y=c(mo_da_link_cut[228,4],mo_da_link_cut[228,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[229,3],mo_da_link_cut[229,5]),y=c(mo_da_link_cut[229,4],mo_da_link_cut[229,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[230,3],mo_da_link_cut[230,5]),y=c(mo_da_link_cut[230,4],mo_da_link_cut[230,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[231,3],mo_da_link_cut[231,5]),y=c(mo_da_link_cut[231,4],mo_da_link_cut[231,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[232,3],mo_da_link_cut[232,5]),y=c(mo_da_link_cut[232,4],mo_da_link_cut[232,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[233,3],mo_da_link_cut[233,5]),y=c(mo_da_link_cut[233,4],mo_da_link_cut[233,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[234,3],mo_da_link_cut[234,5]),y=c(mo_da_link_cut[234,4],mo_da_link_cut[234,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[235,3],mo_da_link_cut[235,5]),y=c(mo_da_link_cut[235,4],mo_da_link_cut[235,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[236,3],mo_da_link_cut[236,5]),y=c(mo_da_link_cut[236,4],mo_da_link_cut[236,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[237,3],mo_da_link_cut[237,5]),y=c(mo_da_link_cut[237,4],mo_da_link_cut[237,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[238,3],mo_da_link_cut[238,5]),y=c(mo_da_link_cut[238,4],mo_da_link_cut[238,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[239,3],mo_da_link_cut[239,5]),y=c(mo_da_link_cut[239,4],mo_da_link_cut[239,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[240,3],mo_da_link_cut[240,5]),y=c(mo_da_link_cut[240,4],mo_da_link_cut[240,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[241,3],mo_da_link_cut[241,5]),y=c(mo_da_link_cut[241,4],mo_da_link_cut[241,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[242,3],mo_da_link_cut[242,5]),y=c(mo_da_link_cut[242,4],mo_da_link_cut[242,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[243,3],mo_da_link_cut[243,5]),y=c(mo_da_link_cut[243,4],mo_da_link_cut[243,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[244,3],mo_da_link_cut[244,5]),y=c(mo_da_link_cut[244,4],mo_da_link_cut[244,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[245,3],mo_da_link_cut[245,5]),y=c(mo_da_link_cut[245,4],mo_da_link_cut[245,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[246,3],mo_da_link_cut[246,5]),y=c(mo_da_link_cut[246,4],mo_da_link_cut[246,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[247,3],mo_da_link_cut[247,5]),y=c(mo_da_link_cut[247,4],mo_da_link_cut[247,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[248,3],mo_da_link_cut[248,5]),y=c(mo_da_link_cut[248,4],mo_da_link_cut[248,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[249,3],mo_da_link_cut[249,5]),y=c(mo_da_link_cut[249,4],mo_da_link_cut[249,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[250,3],mo_da_link_cut[250,5]),y=c(mo_da_link_cut[250,4],mo_da_link_cut[250,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[251,3],mo_da_link_cut[251,5]),y=c(mo_da_link_cut[251,4],mo_da_link_cut[251,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[252,3],mo_da_link_cut[252,5]),y=c(mo_da_link_cut[252,4],mo_da_link_cut[252,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[253,3],mo_da_link_cut[253,5]),y=c(mo_da_link_cut[253,4],mo_da_link_cut[253,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[254,3],mo_da_link_cut[254,5]),y=c(mo_da_link_cut[254,4],mo_da_link_cut[254,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[255,3],mo_da_link_cut[255,5]),y=c(mo_da_link_cut[255,4],mo_da_link_cut[255,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[256,3],mo_da_link_cut[256,5]),y=c(mo_da_link_cut[256,4],mo_da_link_cut[256,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[257,3],mo_da_link_cut[257,5]),y=c(mo_da_link_cut[257,4],mo_da_link_cut[257,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[258,3],mo_da_link_cut[258,5]),y=c(mo_da_link_cut[258,4],mo_da_link_cut[258,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[259,3],mo_da_link_cut[259,5]),y=c(mo_da_link_cut[259,4],mo_da_link_cut[259,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[260,3],mo_da_link_cut[260,5]),y=c(mo_da_link_cut[260,4],mo_da_link_cut[260,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[261,3],mo_da_link_cut[261,5]),y=c(mo_da_link_cut[261,4],mo_da_link_cut[261,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[262,3],mo_da_link_cut[262,5]),y=c(mo_da_link_cut[262,4],mo_da_link_cut[262,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[263,3],mo_da_link_cut[263,5]),y=c(mo_da_link_cut[263,4],mo_da_link_cut[263,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[264,3],mo_da_link_cut[264,5]),y=c(mo_da_link_cut[264,4],mo_da_link_cut[264,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[265,3],mo_da_link_cut[265,5]),y=c(mo_da_link_cut[265,4],mo_da_link_cut[265,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[266,3],mo_da_link_cut[266,5]),y=c(mo_da_link_cut[266,4],mo_da_link_cut[266,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[267,3],mo_da_link_cut[267,5]),y=c(mo_da_link_cut[267,4],mo_da_link_cut[267,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[268,3],mo_da_link_cut[268,5]),y=c(mo_da_link_cut[268,4],mo_da_link_cut[268,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[269,3],mo_da_link_cut[269,5]),y=c(mo_da_link_cut[269,4],mo_da_link_cut[269,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[270,3],mo_da_link_cut[270,5]),y=c(mo_da_link_cut[270,4],mo_da_link_cut[270,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[271,3],mo_da_link_cut[271,5]),y=c(mo_da_link_cut[271,4],mo_da_link_cut[271,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[272,3],mo_da_link_cut[272,5]),y=c(mo_da_link_cut[272,4],mo_da_link_cut[272,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[273,3],mo_da_link_cut[273,5]),y=c(mo_da_link_cut[273,4],mo_da_link_cut[273,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[274,3],mo_da_link_cut[274,5]),y=c(mo_da_link_cut[274,4],mo_da_link_cut[274,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[275,3],mo_da_link_cut[275,5]),y=c(mo_da_link_cut[275,4],mo_da_link_cut[275,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[276,3],mo_da_link_cut[276,5]),y=c(mo_da_link_cut[276,4],mo_da_link_cut[276,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[277,3],mo_da_link_cut[277,5]),y=c(mo_da_link_cut[277,4],mo_da_link_cut[277,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[278,3],mo_da_link_cut[278,5]),y=c(mo_da_link_cut[278,4],mo_da_link_cut[278,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[279,3],mo_da_link_cut[279,5]),y=c(mo_da_link_cut[279,4],mo_da_link_cut[279,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[280,3],mo_da_link_cut[280,5]),y=c(mo_da_link_cut[280,4],mo_da_link_cut[280,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[281,3],mo_da_link_cut[281,5]),y=c(mo_da_link_cut[281,4],mo_da_link_cut[281,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[282,3],mo_da_link_cut[282,5]),y=c(mo_da_link_cut[282,4],mo_da_link_cut[282,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[283,3],mo_da_link_cut[283,5]),y=c(mo_da_link_cut[283,4],mo_da_link_cut[283,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[284,3],mo_da_link_cut[284,5]),y=c(mo_da_link_cut[284,4],mo_da_link_cut[284,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[285,3],mo_da_link_cut[285,5]),y=c(mo_da_link_cut[285,4],mo_da_link_cut[285,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[286,3],mo_da_link_cut[286,5]),y=c(mo_da_link_cut[286,4],mo_da_link_cut[286,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[287,3],mo_da_link_cut[287,5]),y=c(mo_da_link_cut[287,4],mo_da_link_cut[287,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[288,3],mo_da_link_cut[288,5]),y=c(mo_da_link_cut[288,4],mo_da_link_cut[288,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[289,3],mo_da_link_cut[289,5]),y=c(mo_da_link_cut[289,4],mo_da_link_cut[289,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[290,3],mo_da_link_cut[290,5]),y=c(mo_da_link_cut[290,4],mo_da_link_cut[290,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[291,3],mo_da_link_cut[291,5]),y=c(mo_da_link_cut[291,4],mo_da_link_cut[291,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[292,3],mo_da_link_cut[292,5]),y=c(mo_da_link_cut[292,4],mo_da_link_cut[292,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[293,3],mo_da_link_cut[293,5]),y=c(mo_da_link_cut[293,4],mo_da_link_cut[293,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[294,3],mo_da_link_cut[294,5]),y=c(mo_da_link_cut[294,4],mo_da_link_cut[294,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[295,3],mo_da_link_cut[295,5]),y=c(mo_da_link_cut[295,4],mo_da_link_cut[295,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[296,3],mo_da_link_cut[296,5]),y=c(mo_da_link_cut[296,4],mo_da_link_cut[296,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[297,3],mo_da_link_cut[297,5]),y=c(mo_da_link_cut[297,4],mo_da_link_cut[297,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[298,3],mo_da_link_cut[298,5]),y=c(mo_da_link_cut[298,4],mo_da_link_cut[298,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[299,3],mo_da_link_cut[299,5]),y=c(mo_da_link_cut[299,4],mo_da_link_cut[299,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[300,3],mo_da_link_cut[300,5]),y=c(mo_da_link_cut[300,4],mo_da_link_cut[300,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[301,3],mo_da_link_cut[301,5]),y=c(mo_da_link_cut[301,4],mo_da_link_cut[301,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[302,3],mo_da_link_cut[302,5]),y=c(mo_da_link_cut[302,4],mo_da_link_cut[302,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[303,3],mo_da_link_cut[303,5]),y=c(mo_da_link_cut[303,4],mo_da_link_cut[303,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[304,3],mo_da_link_cut[304,5]),y=c(mo_da_link_cut[304,4],mo_da_link_cut[304,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[305,3],mo_da_link_cut[305,5]),y=c(mo_da_link_cut[305,4],mo_da_link_cut[305,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[306,3],mo_da_link_cut[306,5]),y=c(mo_da_link_cut[306,4],mo_da_link_cut[306,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[307,3],mo_da_link_cut[307,5]),y=c(mo_da_link_cut[307,4],mo_da_link_cut[307,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[308,3],mo_da_link_cut[308,5]),y=c(mo_da_link_cut[308,4],mo_da_link_cut[308,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[309,3],mo_da_link_cut[309,5]),y=c(mo_da_link_cut[309,4],mo_da_link_cut[309,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[310,3],mo_da_link_cut[310,5]),y=c(mo_da_link_cut[310,4],mo_da_link_cut[310,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[311,3],mo_da_link_cut[311,5]),y=c(mo_da_link_cut[311,4],mo_da_link_cut[311,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[312,3],mo_da_link_cut[312,5]),y=c(mo_da_link_cut[312,4],mo_da_link_cut[312,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[313,3],mo_da_link_cut[313,5]),y=c(mo_da_link_cut[313,4],mo_da_link_cut[313,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[314,3],mo_da_link_cut[314,5]),y=c(mo_da_link_cut[314,4],mo_da_link_cut[314,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[315,3],mo_da_link_cut[315,5]),y=c(mo_da_link_cut[315,4],mo_da_link_cut[315,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[316,3],mo_da_link_cut[316,5]),y=c(mo_da_link_cut[316,4],mo_da_link_cut[316,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[317,3],mo_da_link_cut[317,5]),y=c(mo_da_link_cut[317,4],mo_da_link_cut[317,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[318,3],mo_da_link_cut[318,5]),y=c(mo_da_link_cut[318,4],mo_da_link_cut[318,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[319,3],mo_da_link_cut[319,5]),y=c(mo_da_link_cut[319,4],mo_da_link_cut[319,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[320,3],mo_da_link_cut[320,5]),y=c(mo_da_link_cut[320,4],mo_da_link_cut[320,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[321,3],mo_da_link_cut[321,5]),y=c(mo_da_link_cut[321,4],mo_da_link_cut[321,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[322,3],mo_da_link_cut[322,5]),y=c(mo_da_link_cut[322,4],mo_da_link_cut[322,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[323,3],mo_da_link_cut[323,5]),y=c(mo_da_link_cut[323,4],mo_da_link_cut[323,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[324,3],mo_da_link_cut[324,5]),y=c(mo_da_link_cut[324,4],mo_da_link_cut[324,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[325,3],mo_da_link_cut[325,5]),y=c(mo_da_link_cut[325,4],mo_da_link_cut[325,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[326,3],mo_da_link_cut[326,5]),y=c(mo_da_link_cut[326,4],mo_da_link_cut[326,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[327,3],mo_da_link_cut[327,5]),y=c(mo_da_link_cut[327,4],mo_da_link_cut[327,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[328,3],mo_da_link_cut[328,5]),y=c(mo_da_link_cut[328,4],mo_da_link_cut[328,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[329,3],mo_da_link_cut[329,5]),y=c(mo_da_link_cut[329,4],mo_da_link_cut[329,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[330,3],mo_da_link_cut[330,5]),y=c(mo_da_link_cut[330,4],mo_da_link_cut[330,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[331,3],mo_da_link_cut[331,5]),y=c(mo_da_link_cut[331,4],mo_da_link_cut[331,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[332,3],mo_da_link_cut[332,5]),y=c(mo_da_link_cut[332,4],mo_da_link_cut[332,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[333,3],mo_da_link_cut[333,5]),y=c(mo_da_link_cut[333,4],mo_da_link_cut[333,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[334,3],mo_da_link_cut[334,5]),y=c(mo_da_link_cut[334,4],mo_da_link_cut[334,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[335,3],mo_da_link_cut[335,5]),y=c(mo_da_link_cut[335,4],mo_da_link_cut[335,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[336,3],mo_da_link_cut[336,5]),y=c(mo_da_link_cut[336,4],mo_da_link_cut[336,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[337,3],mo_da_link_cut[337,5]),y=c(mo_da_link_cut[337,4],mo_da_link_cut[337,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[338,3],mo_da_link_cut[338,5]),y=c(mo_da_link_cut[338,4],mo_da_link_cut[338,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[339,3],mo_da_link_cut[339,5]),y=c(mo_da_link_cut[339,4],mo_da_link_cut[339,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[340,3],mo_da_link_cut[340,5]),y=c(mo_da_link_cut[340,4],mo_da_link_cut[340,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[341,3],mo_da_link_cut[341,5]),y=c(mo_da_link_cut[341,4],mo_da_link_cut[341,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[342,3],mo_da_link_cut[342,5]),y=c(mo_da_link_cut[342,4],mo_da_link_cut[342,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[343,3],mo_da_link_cut[343,5]),y=c(mo_da_link_cut[343,4],mo_da_link_cut[343,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[344,3],mo_da_link_cut[344,5]),y=c(mo_da_link_cut[344,4],mo_da_link_cut[344,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[345,3],mo_da_link_cut[345,5]),y=c(mo_da_link_cut[345,4],mo_da_link_cut[345,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[346,3],mo_da_link_cut[346,5]),y=c(mo_da_link_cut[346,4],mo_da_link_cut[346,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[347,3],mo_da_link_cut[347,5]),y=c(mo_da_link_cut[347,4],mo_da_link_cut[347,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[348,3],mo_da_link_cut[348,5]),y=c(mo_da_link_cut[348,4],mo_da_link_cut[348,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[349,3],mo_da_link_cut[349,5]),y=c(mo_da_link_cut[349,4],mo_da_link_cut[349,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[350,3],mo_da_link_cut[350,5]),y=c(mo_da_link_cut[350,4],mo_da_link_cut[350,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[351,3],mo_da_link_cut[351,5]),y=c(mo_da_link_cut[351,4],mo_da_link_cut[351,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[352,3],mo_da_link_cut[352,5]),y=c(mo_da_link_cut[352,4],mo_da_link_cut[352,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[353,3],mo_da_link_cut[353,5]),y=c(mo_da_link_cut[353,4],mo_da_link_cut[353,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[354,3],mo_da_link_cut[354,5]),y=c(mo_da_link_cut[354,4],mo_da_link_cut[354,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[355,3],mo_da_link_cut[355,5]),y=c(mo_da_link_cut[355,4],mo_da_link_cut[355,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[356,3],mo_da_link_cut[356,5]),y=c(mo_da_link_cut[356,4],mo_da_link_cut[356,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[357,3],mo_da_link_cut[357,5]),y=c(mo_da_link_cut[357,4],mo_da_link_cut[357,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[358,3],mo_da_link_cut[358,5]),y=c(mo_da_link_cut[358,4],mo_da_link_cut[358,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[359,3],mo_da_link_cut[359,5]),y=c(mo_da_link_cut[359,4],mo_da_link_cut[359,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[360,3],mo_da_link_cut[360,5]),y=c(mo_da_link_cut[360,4],mo_da_link_cut[360,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[361,3],mo_da_link_cut[361,5]),y=c(mo_da_link_cut[361,4],mo_da_link_cut[361,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[362,3],mo_da_link_cut[362,5]),y=c(mo_da_link_cut[362,4],mo_da_link_cut[362,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[363,3],mo_da_link_cut[363,5]),y=c(mo_da_link_cut[363,4],mo_da_link_cut[363,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[364,3],mo_da_link_cut[364,5]),y=c(mo_da_link_cut[364,4],mo_da_link_cut[364,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[365,3],mo_da_link_cut[365,5]),y=c(mo_da_link_cut[365,4],mo_da_link_cut[365,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[366,3],mo_da_link_cut[366,5]),y=c(mo_da_link_cut[366,4],mo_da_link_cut[366,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[367,3],mo_da_link_cut[367,5]),y=c(mo_da_link_cut[367,4],mo_da_link_cut[367,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[368,3],mo_da_link_cut[368,5]),y=c(mo_da_link_cut[368,4],mo_da_link_cut[368,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[369,3],mo_da_link_cut[369,5]),y=c(mo_da_link_cut[369,4],mo_da_link_cut[369,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[370,3],mo_da_link_cut[370,5]),y=c(mo_da_link_cut[370,4],mo_da_link_cut[370,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[371,3],mo_da_link_cut[371,5]),y=c(mo_da_link_cut[371,4],mo_da_link_cut[371,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[372,3],mo_da_link_cut[372,5]),y=c(mo_da_link_cut[372,4],mo_da_link_cut[372,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[373,3],mo_da_link_cut[373,5]),y=c(mo_da_link_cut[373,4],mo_da_link_cut[373,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[374,3],mo_da_link_cut[374,5]),y=c(mo_da_link_cut[374,4],mo_da_link_cut[374,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[375,3],mo_da_link_cut[375,5]),y=c(mo_da_link_cut[375,4],mo_da_link_cut[375,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[376,3],mo_da_link_cut[376,5]),y=c(mo_da_link_cut[376,4],mo_da_link_cut[376,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[377,3],mo_da_link_cut[377,5]),y=c(mo_da_link_cut[377,4],mo_da_link_cut[377,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[378,3],mo_da_link_cut[378,5]),y=c(mo_da_link_cut[378,4],mo_da_link_cut[378,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[379,3],mo_da_link_cut[379,5]),y=c(mo_da_link_cut[379,4],mo_da_link_cut[379,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[380,3],mo_da_link_cut[380,5]),y=c(mo_da_link_cut[380,4],mo_da_link_cut[380,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[381,3],mo_da_link_cut[381,5]),y=c(mo_da_link_cut[381,4],mo_da_link_cut[381,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[382,3],mo_da_link_cut[382,5]),y=c(mo_da_link_cut[382,4],mo_da_link_cut[382,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[383,3],mo_da_link_cut[383,5]),y=c(mo_da_link_cut[383,4],mo_da_link_cut[383,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[384,3],mo_da_link_cut[384,5]),y=c(mo_da_link_cut[384,4],mo_da_link_cut[384,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[385,3],mo_da_link_cut[385,5]),y=c(mo_da_link_cut[385,4],mo_da_link_cut[385,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[386,3],mo_da_link_cut[386,5]),y=c(mo_da_link_cut[386,4],mo_da_link_cut[386,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[387,3],mo_da_link_cut[387,5]),y=c(mo_da_link_cut[387,4],mo_da_link_cut[387,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[388,3],mo_da_link_cut[388,5]),y=c(mo_da_link_cut[388,4],mo_da_link_cut[388,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[389,3],mo_da_link_cut[389,5]),y=c(mo_da_link_cut[389,4],mo_da_link_cut[389,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[390,3],mo_da_link_cut[390,5]),y=c(mo_da_link_cut[390,4],mo_da_link_cut[390,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[391,3],mo_da_link_cut[391,5]),y=c(mo_da_link_cut[391,4],mo_da_link_cut[391,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[392,3],mo_da_link_cut[392,5]),y=c(mo_da_link_cut[392,4],mo_da_link_cut[392,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[393,3],mo_da_link_cut[393,5]),y=c(mo_da_link_cut[393,4],mo_da_link_cut[393,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[394,3],mo_da_link_cut[394,5]),y=c(mo_da_link_cut[394,4],mo_da_link_cut[394,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[395,3],mo_da_link_cut[395,5]),y=c(mo_da_link_cut[395,4],mo_da_link_cut[395,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[396,3],mo_da_link_cut[396,5]),y=c(mo_da_link_cut[396,4],mo_da_link_cut[396,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[397,3],mo_da_link_cut[397,5]),y=c(mo_da_link_cut[397,4],mo_da_link_cut[397,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[398,3],mo_da_link_cut[398,5]),y=c(mo_da_link_cut[398,4],mo_da_link_cut[398,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[399,3],mo_da_link_cut[399,5]),y=c(mo_da_link_cut[399,4],mo_da_link_cut[399,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[400,3],mo_da_link_cut[400,5]),y=c(mo_da_link_cut[400,4],mo_da_link_cut[400,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[401,3],mo_da_link_cut[401,5]),y=c(mo_da_link_cut[401,4],mo_da_link_cut[401,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[402,3],mo_da_link_cut[402,5]),y=c(mo_da_link_cut[402,4],mo_da_link_cut[402,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[403,3],mo_da_link_cut[403,5]),y=c(mo_da_link_cut[403,4],mo_da_link_cut[403,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[404,3],mo_da_link_cut[404,5]),y=c(mo_da_link_cut[404,4],mo_da_link_cut[404,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[405,3],mo_da_link_cut[405,5]),y=c(mo_da_link_cut[405,4],mo_da_link_cut[405,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[406,3],mo_da_link_cut[406,5]),y=c(mo_da_link_cut[406,4],mo_da_link_cut[406,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[407,3],mo_da_link_cut[407,5]),y=c(mo_da_link_cut[407,4],mo_da_link_cut[407,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[408,3],mo_da_link_cut[408,5]),y=c(mo_da_link_cut[408,4],mo_da_link_cut[408,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[409,3],mo_da_link_cut[409,5]),y=c(mo_da_link_cut[409,4],mo_da_link_cut[409,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[410,3],mo_da_link_cut[410,5]),y=c(mo_da_link_cut[410,4],mo_da_link_cut[410,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[411,3],mo_da_link_cut[411,5]),y=c(mo_da_link_cut[411,4],mo_da_link_cut[411,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[412,3],mo_da_link_cut[412,5]),y=c(mo_da_link_cut[412,4],mo_da_link_cut[412,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[413,3],mo_da_link_cut[413,5]),y=c(mo_da_link_cut[413,4],mo_da_link_cut[413,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[414,3],mo_da_link_cut[414,5]),y=c(mo_da_link_cut[414,4],mo_da_link_cut[414,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[415,3],mo_da_link_cut[415,5]),y=c(mo_da_link_cut[415,4],mo_da_link_cut[415,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[416,3],mo_da_link_cut[416,5]),y=c(mo_da_link_cut[416,4],mo_da_link_cut[416,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[417,3],mo_da_link_cut[417,5]),y=c(mo_da_link_cut[417,4],mo_da_link_cut[417,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[418,3],mo_da_link_cut[418,5]),y=c(mo_da_link_cut[418,4],mo_da_link_cut[418,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[419,3],mo_da_link_cut[419,5]),y=c(mo_da_link_cut[419,4],mo_da_link_cut[419,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[420,3],mo_da_link_cut[420,5]),y=c(mo_da_link_cut[420,4],mo_da_link_cut[420,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[421,3],mo_da_link_cut[421,5]),y=c(mo_da_link_cut[421,4],mo_da_link_cut[421,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[422,3],mo_da_link_cut[422,5]),y=c(mo_da_link_cut[422,4],mo_da_link_cut[422,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[423,3],mo_da_link_cut[423,5]),y=c(mo_da_link_cut[423,4],mo_da_link_cut[423,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[424,3],mo_da_link_cut[424,5]),y=c(mo_da_link_cut[424,4],mo_da_link_cut[424,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[425,3],mo_da_link_cut[425,5]),y=c(mo_da_link_cut[425,4],mo_da_link_cut[425,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[426,3],mo_da_link_cut[426,5]),y=c(mo_da_link_cut[426,4],mo_da_link_cut[426,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[427,3],mo_da_link_cut[427,5]),y=c(mo_da_link_cut[427,4],mo_da_link_cut[427,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[428,3],mo_da_link_cut[428,5]),y=c(mo_da_link_cut[428,4],mo_da_link_cut[428,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[429,3],mo_da_link_cut[429,5]),y=c(mo_da_link_cut[429,4],mo_da_link_cut[429,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[430,3],mo_da_link_cut[430,5]),y=c(mo_da_link_cut[430,4],mo_da_link_cut[430,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[431,3],mo_da_link_cut[431,5]),y=c(mo_da_link_cut[431,4],mo_da_link_cut[431,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[432,3],mo_da_link_cut[432,5]),y=c(mo_da_link_cut[432,4],mo_da_link_cut[432,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[433,3],mo_da_link_cut[433,5]),y=c(mo_da_link_cut[433,4],mo_da_link_cut[433,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[434,3],mo_da_link_cut[434,5]),y=c(mo_da_link_cut[434,4],mo_da_link_cut[434,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[435,3],mo_da_link_cut[435,5]),y=c(mo_da_link_cut[435,4],mo_da_link_cut[435,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[436,3],mo_da_link_cut[436,5]),y=c(mo_da_link_cut[436,4],mo_da_link_cut[436,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[437,3],mo_da_link_cut[437,5]),y=c(mo_da_link_cut[437,4],mo_da_link_cut[437,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[438,3],mo_da_link_cut[438,5]),y=c(mo_da_link_cut[438,4],mo_da_link_cut[438,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[439,3],mo_da_link_cut[439,5]),y=c(mo_da_link_cut[439,4],mo_da_link_cut[439,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[440,3],mo_da_link_cut[440,5]),y=c(mo_da_link_cut[440,4],mo_da_link_cut[440,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[441,3],mo_da_link_cut[441,5]),y=c(mo_da_link_cut[441,4],mo_da_link_cut[441,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[442,3],mo_da_link_cut[442,5]),y=c(mo_da_link_cut[442,4],mo_da_link_cut[442,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[443,3],mo_da_link_cut[443,5]),y=c(mo_da_link_cut[443,4],mo_da_link_cut[443,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[444,3],mo_da_link_cut[444,5]),y=c(mo_da_link_cut[444,4],mo_da_link_cut[444,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[445,3],mo_da_link_cut[445,5]),y=c(mo_da_link_cut[445,4],mo_da_link_cut[445,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[446,3],mo_da_link_cut[446,5]),y=c(mo_da_link_cut[446,4],mo_da_link_cut[446,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[447,3],mo_da_link_cut[447,5]),y=c(mo_da_link_cut[447,4],mo_da_link_cut[447,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[448,3],mo_da_link_cut[448,5]),y=c(mo_da_link_cut[448,4],mo_da_link_cut[448,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[449,3],mo_da_link_cut[449,5]),y=c(mo_da_link_cut[449,4],mo_da_link_cut[449,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[450,3],mo_da_link_cut[450,5]),y=c(mo_da_link_cut[450,4],mo_da_link_cut[450,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[451,3],mo_da_link_cut[451,5]),y=c(mo_da_link_cut[451,4],mo_da_link_cut[451,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[452,3],mo_da_link_cut[452,5]),y=c(mo_da_link_cut[452,4],mo_da_link_cut[452,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[453,3],mo_da_link_cut[453,5]),y=c(mo_da_link_cut[453,4],mo_da_link_cut[453,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[454,3],mo_da_link_cut[454,5]),y=c(mo_da_link_cut[454,4],mo_da_link_cut[454,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[455,3],mo_da_link_cut[455,5]),y=c(mo_da_link_cut[455,4],mo_da_link_cut[455,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[456,3],mo_da_link_cut[456,5]),y=c(mo_da_link_cut[456,4],mo_da_link_cut[456,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[457,3],mo_da_link_cut[457,5]),y=c(mo_da_link_cut[457,4],mo_da_link_cut[457,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[458,3],mo_da_link_cut[458,5]),y=c(mo_da_link_cut[458,4],mo_da_link_cut[458,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[459,3],mo_da_link_cut[459,5]),y=c(mo_da_link_cut[459,4],mo_da_link_cut[459,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[460,3],mo_da_link_cut[460,5]),y=c(mo_da_link_cut[460,4],mo_da_link_cut[460,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[461,3],mo_da_link_cut[461,5]),y=c(mo_da_link_cut[461,4],mo_da_link_cut[461,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[462,3],mo_da_link_cut[462,5]),y=c(mo_da_link_cut[462,4],mo_da_link_cut[462,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[463,3],mo_da_link_cut[463,5]),y=c(mo_da_link_cut[463,4],mo_da_link_cut[463,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[464,3],mo_da_link_cut[464,5]),y=c(mo_da_link_cut[464,4],mo_da_link_cut[464,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[465,3],mo_da_link_cut[465,5]),y=c(mo_da_link_cut[465,4],mo_da_link_cut[465,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[466,3],mo_da_link_cut[466,5]),y=c(mo_da_link_cut[466,4],mo_da_link_cut[466,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[467,3],mo_da_link_cut[467,5]),y=c(mo_da_link_cut[467,4],mo_da_link_cut[467,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[468,3],mo_da_link_cut[468,5]),y=c(mo_da_link_cut[468,4],mo_da_link_cut[468,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[469,3],mo_da_link_cut[469,5]),y=c(mo_da_link_cut[469,4],mo_da_link_cut[469,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[470,3],mo_da_link_cut[470,5]),y=c(mo_da_link_cut[470,4],mo_da_link_cut[470,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[471,3],mo_da_link_cut[471,5]),y=c(mo_da_link_cut[471,4],mo_da_link_cut[471,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[472,3],mo_da_link_cut[472,5]),y=c(mo_da_link_cut[472,4],mo_da_link_cut[472,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[473,3],mo_da_link_cut[473,5]),y=c(mo_da_link_cut[473,4],mo_da_link_cut[473,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[474,3],mo_da_link_cut[474,5]),y=c(mo_da_link_cut[474,4],mo_da_link_cut[474,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[475,3],mo_da_link_cut[475,5]),y=c(mo_da_link_cut[475,4],mo_da_link_cut[475,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[476,3],mo_da_link_cut[476,5]),y=c(mo_da_link_cut[476,4],mo_da_link_cut[476,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[477,3],mo_da_link_cut[477,5]),y=c(mo_da_link_cut[477,4],mo_da_link_cut[477,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[478,3],mo_da_link_cut[478,5]),y=c(mo_da_link_cut[478,4],mo_da_link_cut[478,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[479,3],mo_da_link_cut[479,5]),y=c(mo_da_link_cut[479,4],mo_da_link_cut[479,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[480,3],mo_da_link_cut[480,5]),y=c(mo_da_link_cut[480,4],mo_da_link_cut[480,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[481,3],mo_da_link_cut[481,5]),y=c(mo_da_link_cut[481,4],mo_da_link_cut[481,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[482,3],mo_da_link_cut[482,5]),y=c(mo_da_link_cut[482,4],mo_da_link_cut[482,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[483,3],mo_da_link_cut[483,5]),y=c(mo_da_link_cut[483,4],mo_da_link_cut[483,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[484,3],mo_da_link_cut[484,5]),y=c(mo_da_link_cut[484,4],mo_da_link_cut[484,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[485,3],mo_da_link_cut[485,5]),y=c(mo_da_link_cut[485,4],mo_da_link_cut[485,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[486,3],mo_da_link_cut[486,5]),y=c(mo_da_link_cut[486,4],mo_da_link_cut[486,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[487,3],mo_da_link_cut[487,5]),y=c(mo_da_link_cut[487,4],mo_da_link_cut[487,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[488,3],mo_da_link_cut[488,5]),y=c(mo_da_link_cut[488,4],mo_da_link_cut[488,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[489,3],mo_da_link_cut[489,5]),y=c(mo_da_link_cut[489,4],mo_da_link_cut[489,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[490,3],mo_da_link_cut[490,5]),y=c(mo_da_link_cut[490,4],mo_da_link_cut[490,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[491,3],mo_da_link_cut[491,5]),y=c(mo_da_link_cut[491,4],mo_da_link_cut[491,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[492,3],mo_da_link_cut[492,5]),y=c(mo_da_link_cut[492,4],mo_da_link_cut[492,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[493,3],mo_da_link_cut[493,5]),y=c(mo_da_link_cut[493,4],mo_da_link_cut[493,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[494,3],mo_da_link_cut[494,5]),y=c(mo_da_link_cut[494,4],mo_da_link_cut[494,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[495,3],mo_da_link_cut[495,5]),y=c(mo_da_link_cut[495,4],mo_da_link_cut[495,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[496,3],mo_da_link_cut[496,5]),y=c(mo_da_link_cut[496,4],mo_da_link_cut[496,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[497,3],mo_da_link_cut[497,5]),y=c(mo_da_link_cut[497,4],mo_da_link_cut[497,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[498,3],mo_da_link_cut[498,5]),y=c(mo_da_link_cut[498,4],mo_da_link_cut[498,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[499,3],mo_da_link_cut[499,5]),y=c(mo_da_link_cut[499,4],mo_da_link_cut[499,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[500,3],mo_da_link_cut[500,5]),y=c(mo_da_link_cut[500,4],mo_da_link_cut[500,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[501,3],mo_da_link_cut[501,5]),y=c(mo_da_link_cut[501,4],mo_da_link_cut[501,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[502,3],mo_da_link_cut[502,5]),y=c(mo_da_link_cut[502,4],mo_da_link_cut[502,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[503,3],mo_da_link_cut[503,5]),y=c(mo_da_link_cut[503,4],mo_da_link_cut[503,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[504,3],mo_da_link_cut[504,5]),y=c(mo_da_link_cut[504,4],mo_da_link_cut[504,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[505,3],mo_da_link_cut[505,5]),y=c(mo_da_link_cut[505,4],mo_da_link_cut[505,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[506,3],mo_da_link_cut[506,5]),y=c(mo_da_link_cut[506,4],mo_da_link_cut[506,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[507,3],mo_da_link_cut[507,5]),y=c(mo_da_link_cut[507,4],mo_da_link_cut[507,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[508,3],mo_da_link_cut[508,5]),y=c(mo_da_link_cut[508,4],mo_da_link_cut[508,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[509,3],mo_da_link_cut[509,5]),y=c(mo_da_link_cut[509,4],mo_da_link_cut[509,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[510,3],mo_da_link_cut[510,5]),y=c(mo_da_link_cut[510,4],mo_da_link_cut[510,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[511,3],mo_da_link_cut[511,5]),y=c(mo_da_link_cut[511,4],mo_da_link_cut[511,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[512,3],mo_da_link_cut[512,5]),y=c(mo_da_link_cut[512,4],mo_da_link_cut[512,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[513,3],mo_da_link_cut[513,5]),y=c(mo_da_link_cut[513,4],mo_da_link_cut[513,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[514,3],mo_da_link_cut[514,5]),y=c(mo_da_link_cut[514,4],mo_da_link_cut[514,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[515,3],mo_da_link_cut[515,5]),y=c(mo_da_link_cut[515,4],mo_da_link_cut[515,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[516,3],mo_da_link_cut[516,5]),y=c(mo_da_link_cut[516,4],mo_da_link_cut[516,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[517,3],mo_da_link_cut[517,5]),y=c(mo_da_link_cut[517,4],mo_da_link_cut[517,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[518,3],mo_da_link_cut[518,5]),y=c(mo_da_link_cut[518,4],mo_da_link_cut[518,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[519,3],mo_da_link_cut[519,5]),y=c(mo_da_link_cut[519,4],mo_da_link_cut[519,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[520,3],mo_da_link_cut[520,5]),y=c(mo_da_link_cut[520,4],mo_da_link_cut[520,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[521,3],mo_da_link_cut[521,5]),y=c(mo_da_link_cut[521,4],mo_da_link_cut[521,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[522,3],mo_da_link_cut[522,5]),y=c(mo_da_link_cut[522,4],mo_da_link_cut[522,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[523,3],mo_da_link_cut[523,5]),y=c(mo_da_link_cut[523,4],mo_da_link_cut[523,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[524,3],mo_da_link_cut[524,5]),y=c(mo_da_link_cut[524,4],mo_da_link_cut[524,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[525,3],mo_da_link_cut[525,5]),y=c(mo_da_link_cut[525,4],mo_da_link_cut[525,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[526,3],mo_da_link_cut[526,5]),y=c(mo_da_link_cut[526,4],mo_da_link_cut[526,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[527,3],mo_da_link_cut[527,5]),y=c(mo_da_link_cut[527,4],mo_da_link_cut[527,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[528,3],mo_da_link_cut[528,5]),y=c(mo_da_link_cut[528,4],mo_da_link_cut[528,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[529,3],mo_da_link_cut[529,5]),y=c(mo_da_link_cut[529,4],mo_da_link_cut[529,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[530,3],mo_da_link_cut[530,5]),y=c(mo_da_link_cut[530,4],mo_da_link_cut[530,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[531,3],mo_da_link_cut[531,5]),y=c(mo_da_link_cut[531,4],mo_da_link_cut[531,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[532,3],mo_da_link_cut[532,5]),y=c(mo_da_link_cut[532,4],mo_da_link_cut[532,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[533,3],mo_da_link_cut[533,5]),y=c(mo_da_link_cut[533,4],mo_da_link_cut[533,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[534,3],mo_da_link_cut[534,5]),y=c(mo_da_link_cut[534,4],mo_da_link_cut[534,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[535,3],mo_da_link_cut[535,5]),y=c(mo_da_link_cut[535,4],mo_da_link_cut[535,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[536,3],mo_da_link_cut[536,5]),y=c(mo_da_link_cut[536,4],mo_da_link_cut[536,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[537,3],mo_da_link_cut[537,5]),y=c(mo_da_link_cut[537,4],mo_da_link_cut[537,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[538,3],mo_da_link_cut[538,5]),y=c(mo_da_link_cut[538,4],mo_da_link_cut[538,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[539,3],mo_da_link_cut[539,5]),y=c(mo_da_link_cut[539,4],mo_da_link_cut[539,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[540,3],mo_da_link_cut[540,5]),y=c(mo_da_link_cut[540,4],mo_da_link_cut[540,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[541,3],mo_da_link_cut[541,5]),y=c(mo_da_link_cut[541,4],mo_da_link_cut[541,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[542,3],mo_da_link_cut[542,5]),y=c(mo_da_link_cut[542,4],mo_da_link_cut[542,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[543,3],mo_da_link_cut[543,5]),y=c(mo_da_link_cut[543,4],mo_da_link_cut[543,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[544,3],mo_da_link_cut[544,5]),y=c(mo_da_link_cut[544,4],mo_da_link_cut[544,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[545,3],mo_da_link_cut[545,5]),y=c(mo_da_link_cut[545,4],mo_da_link_cut[545,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[546,3],mo_da_link_cut[546,5]),y=c(mo_da_link_cut[546,4],mo_da_link_cut[546,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[547,3],mo_da_link_cut[547,5]),y=c(mo_da_link_cut[547,4],mo_da_link_cut[547,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[548,3],mo_da_link_cut[548,5]),y=c(mo_da_link_cut[548,4],mo_da_link_cut[548,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[549,3],mo_da_link_cut[549,5]),y=c(mo_da_link_cut[549,4],mo_da_link_cut[549,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[550,3],mo_da_link_cut[550,5]),y=c(mo_da_link_cut[550,4],mo_da_link_cut[550,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[551,3],mo_da_link_cut[551,5]),y=c(mo_da_link_cut[551,4],mo_da_link_cut[551,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[552,3],mo_da_link_cut[552,5]),y=c(mo_da_link_cut[552,4],mo_da_link_cut[552,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[553,3],mo_da_link_cut[553,5]),y=c(mo_da_link_cut[553,4],mo_da_link_cut[553,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[554,3],mo_da_link_cut[554,5]),y=c(mo_da_link_cut[554,4],mo_da_link_cut[554,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[555,3],mo_da_link_cut[555,5]),y=c(mo_da_link_cut[555,4],mo_da_link_cut[555,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[556,3],mo_da_link_cut[556,5]),y=c(mo_da_link_cut[556,4],mo_da_link_cut[556,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[557,3],mo_da_link_cut[557,5]),y=c(mo_da_link_cut[557,4],mo_da_link_cut[557,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[558,3],mo_da_link_cut[558,5]),y=c(mo_da_link_cut[558,4],mo_da_link_cut[558,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[559,3],mo_da_link_cut[559,5]),y=c(mo_da_link_cut[559,4],mo_da_link_cut[559,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[560,3],mo_da_link_cut[560,5]),y=c(mo_da_link_cut[560,4],mo_da_link_cut[560,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[561,3],mo_da_link_cut[561,5]),y=c(mo_da_link_cut[561,4],mo_da_link_cut[561,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[562,3],mo_da_link_cut[562,5]),y=c(mo_da_link_cut[562,4],mo_da_link_cut[562,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[563,3],mo_da_link_cut[563,5]),y=c(mo_da_link_cut[563,4],mo_da_link_cut[563,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[564,3],mo_da_link_cut[564,5]),y=c(mo_da_link_cut[564,4],mo_da_link_cut[564,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[565,3],mo_da_link_cut[565,5]),y=c(mo_da_link_cut[565,4],mo_da_link_cut[565,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[566,3],mo_da_link_cut[566,5]),y=c(mo_da_link_cut[566,4],mo_da_link_cut[566,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[567,3],mo_da_link_cut[567,5]),y=c(mo_da_link_cut[567,4],mo_da_link_cut[567,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[568,3],mo_da_link_cut[568,5]),y=c(mo_da_link_cut[568,4],mo_da_link_cut[568,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[569,3],mo_da_link_cut[569,5]),y=c(mo_da_link_cut[569,4],mo_da_link_cut[569,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[570,3],mo_da_link_cut[570,5]),y=c(mo_da_link_cut[570,4],mo_da_link_cut[570,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[571,3],mo_da_link_cut[571,5]),y=c(mo_da_link_cut[571,4],mo_da_link_cut[571,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[572,3],mo_da_link_cut[572,5]),y=c(mo_da_link_cut[572,4],mo_da_link_cut[572,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[573,3],mo_da_link_cut[573,5]),y=c(mo_da_link_cut[573,4],mo_da_link_cut[573,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[574,3],mo_da_link_cut[574,5]),y=c(mo_da_link_cut[574,4],mo_da_link_cut[574,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[575,3],mo_da_link_cut[575,5]),y=c(mo_da_link_cut[575,4],mo_da_link_cut[575,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[576,3],mo_da_link_cut[576,5]),y=c(mo_da_link_cut[576,4],mo_da_link_cut[576,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[577,3],mo_da_link_cut[577,5]),y=c(mo_da_link_cut[577,4],mo_da_link_cut[577,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[578,3],mo_da_link_cut[578,5]),y=c(mo_da_link_cut[578,4],mo_da_link_cut[578,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[579,3],mo_da_link_cut[579,5]),y=c(mo_da_link_cut[579,4],mo_da_link_cut[579,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[580,3],mo_da_link_cut[580,5]),y=c(mo_da_link_cut[580,4],mo_da_link_cut[580,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[581,3],mo_da_link_cut[581,5]),y=c(mo_da_link_cut[581,4],mo_da_link_cut[581,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[582,3],mo_da_link_cut[582,5]),y=c(mo_da_link_cut[582,4],mo_da_link_cut[582,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[583,3],mo_da_link_cut[583,5]),y=c(mo_da_link_cut[583,4],mo_da_link_cut[583,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[584,3],mo_da_link_cut[584,5]),y=c(mo_da_link_cut[584,4],mo_da_link_cut[584,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[585,3],mo_da_link_cut[585,5]),y=c(mo_da_link_cut[585,4],mo_da_link_cut[585,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[586,3],mo_da_link_cut[586,5]),y=c(mo_da_link_cut[586,4],mo_da_link_cut[586,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[587,3],mo_da_link_cut[587,5]),y=c(mo_da_link_cut[587,4],mo_da_link_cut[587,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[588,3],mo_da_link_cut[588,5]),y=c(mo_da_link_cut[588,4],mo_da_link_cut[588,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[589,3],mo_da_link_cut[589,5]),y=c(mo_da_link_cut[589,4],mo_da_link_cut[589,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[590,3],mo_da_link_cut[590,5]),y=c(mo_da_link_cut[590,4],mo_da_link_cut[590,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[591,3],mo_da_link_cut[591,5]),y=c(mo_da_link_cut[591,4],mo_da_link_cut[591,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[592,3],mo_da_link_cut[592,5]),y=c(mo_da_link_cut[592,4],mo_da_link_cut[592,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[593,3],mo_da_link_cut[593,5]),y=c(mo_da_link_cut[593,4],mo_da_link_cut[593,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[594,3],mo_da_link_cut[594,5]),y=c(mo_da_link_cut[594,4],mo_da_link_cut[594,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[595,3],mo_da_link_cut[595,5]),y=c(mo_da_link_cut[595,4],mo_da_link_cut[595,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[596,3],mo_da_link_cut[596,5]),y=c(mo_da_link_cut[596,4],mo_da_link_cut[596,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[597,3],mo_da_link_cut[597,5]),y=c(mo_da_link_cut[597,4],mo_da_link_cut[597,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[598,3],mo_da_link_cut[598,5]),y=c(mo_da_link_cut[598,4],mo_da_link_cut[598,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[599,3],mo_da_link_cut[599,5]),y=c(mo_da_link_cut[599,4],mo_da_link_cut[599,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[600,3],mo_da_link_cut[600,5]),y=c(mo_da_link_cut[600,4],mo_da_link_cut[600,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[601,3],mo_da_link_cut[601,5]),y=c(mo_da_link_cut[601,4],mo_da_link_cut[601,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[602,3],mo_da_link_cut[602,5]),y=c(mo_da_link_cut[602,4],mo_da_link_cut[602,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[603,3],mo_da_link_cut[603,5]),y=c(mo_da_link_cut[603,4],mo_da_link_cut[603,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[604,3],mo_da_link_cut[604,5]),y=c(mo_da_link_cut[604,4],mo_da_link_cut[604,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[605,3],mo_da_link_cut[605,5]),y=c(mo_da_link_cut[605,4],mo_da_link_cut[605,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[606,3],mo_da_link_cut[606,5]),y=c(mo_da_link_cut[606,4],mo_da_link_cut[606,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[607,3],mo_da_link_cut[607,5]),y=c(mo_da_link_cut[607,4],mo_da_link_cut[607,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[608,3],mo_da_link_cut[608,5]),y=c(mo_da_link_cut[608,4],mo_da_link_cut[608,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[609,3],mo_da_link_cut[609,5]),y=c(mo_da_link_cut[609,4],mo_da_link_cut[609,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[610,3],mo_da_link_cut[610,5]),y=c(mo_da_link_cut[610,4],mo_da_link_cut[610,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[611,3],mo_da_link_cut[611,5]),y=c(mo_da_link_cut[611,4],mo_da_link_cut[611,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[612,3],mo_da_link_cut[612,5]),y=c(mo_da_link_cut[612,4],mo_da_link_cut[612,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[613,3],mo_da_link_cut[613,5]),y=c(mo_da_link_cut[613,4],mo_da_link_cut[613,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[614,3],mo_da_link_cut[614,5]),y=c(mo_da_link_cut[614,4],mo_da_link_cut[614,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[615,3],mo_da_link_cut[615,5]),y=c(mo_da_link_cut[615,4],mo_da_link_cut[615,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[616,3],mo_da_link_cut[616,5]),y=c(mo_da_link_cut[616,4],mo_da_link_cut[616,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[617,3],mo_da_link_cut[617,5]),y=c(mo_da_link_cut[617,4],mo_da_link_cut[617,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[618,3],mo_da_link_cut[618,5]),y=c(mo_da_link_cut[618,4],mo_da_link_cut[618,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[619,3],mo_da_link_cut[619,5]),y=c(mo_da_link_cut[619,4],mo_da_link_cut[619,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[620,3],mo_da_link_cut[620,5]),y=c(mo_da_link_cut[620,4],mo_da_link_cut[620,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[621,3],mo_da_link_cut[621,5]),y=c(mo_da_link_cut[621,4],mo_da_link_cut[621,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[622,3],mo_da_link_cut[622,5]),y=c(mo_da_link_cut[622,4],mo_da_link_cut[622,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[623,3],mo_da_link_cut[623,5]),y=c(mo_da_link_cut[623,4],mo_da_link_cut[623,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[624,3],mo_da_link_cut[624,5]),y=c(mo_da_link_cut[624,4],mo_da_link_cut[624,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[625,3],mo_da_link_cut[625,5]),y=c(mo_da_link_cut[625,4],mo_da_link_cut[625,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[626,3],mo_da_link_cut[626,5]),y=c(mo_da_link_cut[626,4],mo_da_link_cut[626,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[627,3],mo_da_link_cut[627,5]),y=c(mo_da_link_cut[627,4],mo_da_link_cut[627,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[628,3],mo_da_link_cut[628,5]),y=c(mo_da_link_cut[628,4],mo_da_link_cut[628,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[629,3],mo_da_link_cut[629,5]),y=c(mo_da_link_cut[629,4],mo_da_link_cut[629,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[630,3],mo_da_link_cut[630,5]),y=c(mo_da_link_cut[630,4],mo_da_link_cut[630,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[631,3],mo_da_link_cut[631,5]),y=c(mo_da_link_cut[631,4],mo_da_link_cut[631,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[632,3],mo_da_link_cut[632,5]),y=c(mo_da_link_cut[632,4],mo_da_link_cut[632,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[633,3],mo_da_link_cut[633,5]),y=c(mo_da_link_cut[633,4],mo_da_link_cut[633,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[634,3],mo_da_link_cut[634,5]),y=c(mo_da_link_cut[634,4],mo_da_link_cut[634,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[635,3],mo_da_link_cut[635,5]),y=c(mo_da_link_cut[635,4],mo_da_link_cut[635,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[636,3],mo_da_link_cut[636,5]),y=c(mo_da_link_cut[636,4],mo_da_link_cut[636,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[637,3],mo_da_link_cut[637,5]),y=c(mo_da_link_cut[637,4],mo_da_link_cut[637,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[638,3],mo_da_link_cut[638,5]),y=c(mo_da_link_cut[638,4],mo_da_link_cut[638,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[639,3],mo_da_link_cut[639,5]),y=c(mo_da_link_cut[639,4],mo_da_link_cut[639,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[640,3],mo_da_link_cut[640,5]),y=c(mo_da_link_cut[640,4],mo_da_link_cut[640,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[641,3],mo_da_link_cut[641,5]),y=c(mo_da_link_cut[641,4],mo_da_link_cut[641,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[642,3],mo_da_link_cut[642,5]),y=c(mo_da_link_cut[642,4],mo_da_link_cut[642,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[643,3],mo_da_link_cut[643,5]),y=c(mo_da_link_cut[643,4],mo_da_link_cut[643,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[644,3],mo_da_link_cut[644,5]),y=c(mo_da_link_cut[644,4],mo_da_link_cut[644,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[645,3],mo_da_link_cut[645,5]),y=c(mo_da_link_cut[645,4],mo_da_link_cut[645,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[646,3],mo_da_link_cut[646,5]),y=c(mo_da_link_cut[646,4],mo_da_link_cut[646,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[647,3],mo_da_link_cut[647,5]),y=c(mo_da_link_cut[647,4],mo_da_link_cut[647,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[648,3],mo_da_link_cut[648,5]),y=c(mo_da_link_cut[648,4],mo_da_link_cut[648,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[649,3],mo_da_link_cut[649,5]),y=c(mo_da_link_cut[649,4],mo_da_link_cut[649,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[650,3],mo_da_link_cut[650,5]),y=c(mo_da_link_cut[650,4],mo_da_link_cut[650,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[651,3],mo_da_link_cut[651,5]),y=c(mo_da_link_cut[651,4],mo_da_link_cut[651,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[652,3],mo_da_link_cut[652,5]),y=c(mo_da_link_cut[652,4],mo_da_link_cut[652,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[653,3],mo_da_link_cut[653,5]),y=c(mo_da_link_cut[653,4],mo_da_link_cut[653,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[654,3],mo_da_link_cut[654,5]),y=c(mo_da_link_cut[654,4],mo_da_link_cut[654,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[655,3],mo_da_link_cut[655,5]),y=c(mo_da_link_cut[655,4],mo_da_link_cut[655,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[656,3],mo_da_link_cut[656,5]),y=c(mo_da_link_cut[656,4],mo_da_link_cut[656,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[657,3],mo_da_link_cut[657,5]),y=c(mo_da_link_cut[657,4],mo_da_link_cut[657,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[658,3],mo_da_link_cut[658,5]),y=c(mo_da_link_cut[658,4],mo_da_link_cut[658,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[659,3],mo_da_link_cut[659,5]),y=c(mo_da_link_cut[659,4],mo_da_link_cut[659,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[660,3],mo_da_link_cut[660,5]),y=c(mo_da_link_cut[660,4],mo_da_link_cut[660,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[661,3],mo_da_link_cut[661,5]),y=c(mo_da_link_cut[661,4],mo_da_link_cut[661,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[662,3],mo_da_link_cut[662,5]),y=c(mo_da_link_cut[662,4],mo_da_link_cut[662,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[663,3],mo_da_link_cut[663,5]),y=c(mo_da_link_cut[663,4],mo_da_link_cut[663,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[664,3],mo_da_link_cut[664,5]),y=c(mo_da_link_cut[664,4],mo_da_link_cut[664,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[665,3],mo_da_link_cut[665,5]),y=c(mo_da_link_cut[665,4],mo_da_link_cut[665,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[666,3],mo_da_link_cut[666,5]),y=c(mo_da_link_cut[666,4],mo_da_link_cut[666,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[667,3],mo_da_link_cut[667,5]),y=c(mo_da_link_cut[667,4],mo_da_link_cut[667,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[668,3],mo_da_link_cut[668,5]),y=c(mo_da_link_cut[668,4],mo_da_link_cut[668,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[669,3],mo_da_link_cut[669,5]),y=c(mo_da_link_cut[669,4],mo_da_link_cut[669,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[670,3],mo_da_link_cut[670,5]),y=c(mo_da_link_cut[670,4],mo_da_link_cut[670,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[671,3],mo_da_link_cut[671,5]),y=c(mo_da_link_cut[671,4],mo_da_link_cut[671,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[672,3],mo_da_link_cut[672,5]),y=c(mo_da_link_cut[672,4],mo_da_link_cut[672,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[673,3],mo_da_link_cut[673,5]),y=c(mo_da_link_cut[673,4],mo_da_link_cut[673,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[674,3],mo_da_link_cut[674,5]),y=c(mo_da_link_cut[674,4],mo_da_link_cut[674,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[675,3],mo_da_link_cut[675,5]),y=c(mo_da_link_cut[675,4],mo_da_link_cut[675,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[676,3],mo_da_link_cut[676,5]),y=c(mo_da_link_cut[676,4],mo_da_link_cut[676,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[677,3],mo_da_link_cut[677,5]),y=c(mo_da_link_cut[677,4],mo_da_link_cut[677,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[678,3],mo_da_link_cut[678,5]),y=c(mo_da_link_cut[678,4],mo_da_link_cut[678,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[679,3],mo_da_link_cut[679,5]),y=c(mo_da_link_cut[679,4],mo_da_link_cut[679,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[680,3],mo_da_link_cut[680,5]),y=c(mo_da_link_cut[680,4],mo_da_link_cut[680,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[681,3],mo_da_link_cut[681,5]),y=c(mo_da_link_cut[681,4],mo_da_link_cut[681,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[682,3],mo_da_link_cut[682,5]),y=c(mo_da_link_cut[682,4],mo_da_link_cut[682,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[683,3],mo_da_link_cut[683,5]),y=c(mo_da_link_cut[683,4],mo_da_link_cut[683,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[684,3],mo_da_link_cut[684,5]),y=c(mo_da_link_cut[684,4],mo_da_link_cut[684,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[685,3],mo_da_link_cut[685,5]),y=c(mo_da_link_cut[685,4],mo_da_link_cut[685,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[686,3],mo_da_link_cut[686,5]),y=c(mo_da_link_cut[686,4],mo_da_link_cut[686,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[687,3],mo_da_link_cut[687,5]),y=c(mo_da_link_cut[687,4],mo_da_link_cut[687,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[688,3],mo_da_link_cut[688,5]),y=c(mo_da_link_cut[688,4],mo_da_link_cut[688,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[689,3],mo_da_link_cut[689,5]),y=c(mo_da_link_cut[689,4],mo_da_link_cut[689,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[690,3],mo_da_link_cut[690,5]),y=c(mo_da_link_cut[690,4],mo_da_link_cut[690,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[691,3],mo_da_link_cut[691,5]),y=c(mo_da_link_cut[691,4],mo_da_link_cut[691,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[692,3],mo_da_link_cut[692,5]),y=c(mo_da_link_cut[692,4],mo_da_link_cut[692,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[693,3],mo_da_link_cut[693,5]),y=c(mo_da_link_cut[693,4],mo_da_link_cut[693,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[694,3],mo_da_link_cut[694,5]),y=c(mo_da_link_cut[694,4],mo_da_link_cut[694,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[695,3],mo_da_link_cut[695,5]),y=c(mo_da_link_cut[695,4],mo_da_link_cut[695,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[696,3],mo_da_link_cut[696,5]),y=c(mo_da_link_cut[696,4],mo_da_link_cut[696,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[697,3],mo_da_link_cut[697,5]),y=c(mo_da_link_cut[697,4],mo_da_link_cut[697,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[698,3],mo_da_link_cut[698,5]),y=c(mo_da_link_cut[698,4],mo_da_link_cut[698,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[699,3],mo_da_link_cut[699,5]),y=c(mo_da_link_cut[699,4],mo_da_link_cut[699,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[700,3],mo_da_link_cut[700,5]),y=c(mo_da_link_cut[700,4],mo_da_link_cut[700,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[701,3],mo_da_link_cut[701,5]),y=c(mo_da_link_cut[701,4],mo_da_link_cut[701,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[702,3],mo_da_link_cut[702,5]),y=c(mo_da_link_cut[702,4],mo_da_link_cut[702,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[703,3],mo_da_link_cut[703,5]),y=c(mo_da_link_cut[703,4],mo_da_link_cut[703,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[704,3],mo_da_link_cut[704,5]),y=c(mo_da_link_cut[704,4],mo_da_link_cut[704,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[705,3],mo_da_link_cut[705,5]),y=c(mo_da_link_cut[705,4],mo_da_link_cut[705,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[706,3],mo_da_link_cut[706,5]),y=c(mo_da_link_cut[706,4],mo_da_link_cut[706,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[707,3],mo_da_link_cut[707,5]),y=c(mo_da_link_cut[707,4],mo_da_link_cut[707,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[708,3],mo_da_link_cut[708,5]),y=c(mo_da_link_cut[708,4],mo_da_link_cut[708,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[709,3],mo_da_link_cut[709,5]),y=c(mo_da_link_cut[709,4],mo_da_link_cut[709,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[710,3],mo_da_link_cut[710,5]),y=c(mo_da_link_cut[710,4],mo_da_link_cut[710,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[711,3],mo_da_link_cut[711,5]),y=c(mo_da_link_cut[711,4],mo_da_link_cut[711,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[712,3],mo_da_link_cut[712,5]),y=c(mo_da_link_cut[712,4],mo_da_link_cut[712,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[713,3],mo_da_link_cut[713,5]),y=c(mo_da_link_cut[713,4],mo_da_link_cut[713,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[714,3],mo_da_link_cut[714,5]),y=c(mo_da_link_cut[714,4],mo_da_link_cut[714,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[715,3],mo_da_link_cut[715,5]),y=c(mo_da_link_cut[715,4],mo_da_link_cut[715,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[716,3],mo_da_link_cut[716,5]),y=c(mo_da_link_cut[716,4],mo_da_link_cut[716,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[717,3],mo_da_link_cut[717,5]),y=c(mo_da_link_cut[717,4],mo_da_link_cut[717,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[718,3],mo_da_link_cut[718,5]),y=c(mo_da_link_cut[718,4],mo_da_link_cut[718,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[719,3],mo_da_link_cut[719,5]),y=c(mo_da_link_cut[719,4],mo_da_link_cut[719,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[720,3],mo_da_link_cut[720,5]),y=c(mo_da_link_cut[720,4],mo_da_link_cut[720,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_point(data=focal_triangle_tree_cut[match(focal_cells,focal_triangle_tree_cut$cell),],aes(x=x,y=y),size=sqrt(focal_triangle_tree_cut$var_exp[match(focal_cells,focal_triangle_tree_cut$cell)])*0.2,stroke=0.5,pch=1,col="#00BFC4")+
geom_point(data=focal_triangle_tree_cut[match(focal_cells,focal_triangle_tree_cut$cell),],aes(x=x,y=y),size=sqrt(focal_triangle_tree_cut$var_obs[match(focal_cells,focal_triangle_tree_cut$cell)])*0.2,stroke=0.1,pch=16,col="#F8766D")+
labs(title="")

}



# x  case

triangle_tree_cut_x<-triangle_tree_cut
triangle_tree_cut_x$var_exp<-rep(NA,dim(triangle_tree_cut)[1])
triangle_tree_cut_x$var_obs<-rep(NA,dim(triangle_tree_cut)[1])

triangle_tree_cut_x$var_exp[match(cell_var_exp_obs_x$cell,triangle_tree_cut$cell)]<-cell_var_exp_obs_x$cell_var_exp
triangle_tree_cut_x$var_obs[match(cell_var_exp_obs_x$cell,triangle_tree_cut$cell)]<-cell_var_exp_obs_x$cell_var_obs

# 

focal_cells<-na.omit(c(lineage_info[4:10,match(c("ABalaaaala","ABaraaaaaa","ABplaaaaaa","ABpraaaaaa","MSaaaaaa","Ealaa","Caaaaa","Daaa"),colnames(lineage_info))]))

plot_triangle_tree(triangle_tree_cut_x,focal_cells)


ggsave("withoutNF_variance_accumulation_x.pdf",width=unit(6.2,"in"),height=unit(2.4,"in"))



# y  case

triangle_tree_cut_y<-triangle_tree_cut

triangle_tree_cut_y$var_exp<-rep(NA,dim(triangle_tree_cut_y)[1])
triangle_tree_cut_y$var_obs<-rep(NA,dim(triangle_tree_cut_y)[1])

triangle_tree_cut_y$var_exp[match(cell_var_exp_obs_y$cell,triangle_tree_cut_y$cell)]<-cell_var_exp_obs_y$cell_var_exp
triangle_tree_cut_y$var_obs[match(cell_var_exp_obs_y$cell,triangle_tree_cut_y$cell)]<-cell_var_exp_obs_y$cell_var_obs

# 

focal_cells<-na.omit(c(lineage_info[4:10,match(c("ABalaaaala","ABaraaaaaa","ABplaaaaaa","ABpraaaaaa","MSaaaaaa","Ealaa","Caaaaa","Daaa"),colnames(lineage_info))]))


plot_triangle_tree(triangle_tree_cut_y,focal_cells)


ggsave("withoutNF_variance_accumulation_y.pdf",width=unit(6.2,"in"),height=unit(2.4,"in"))



# z  case

triangle_tree_cut_z<-triangle_tree_cut


triangle_tree_cut_z$var_exp<-rep(NA,dim(triangle_tree_cut_z)[1])
triangle_tree_cut_z$var_obs<-rep(NA,dim(triangle_tree_cut_z)[1])

triangle_tree_cut_z$var_exp[match(cell_var_exp_obs_z$cell,triangle_tree_cut_z$cell)]<-cell_var_exp_obs_z$cell_var_exp
triangle_tree_cut_z$var_obs[match(cell_var_exp_obs_z$cell,triangle_tree_cut_z$cell)]<-cell_var_exp_obs_z$cell_var_obs

# 

focal_cells<-na.omit(c(lineage_info[4:10,match(c("ABalaaaala","ABaraaaaaa","ABplaaaaaa","ABpraaaaaa","MSaaaaaa","Ealaa","Caaaaa","Daaa"),colnames(lineage_info))]))


plot_triangle_tree(triangle_tree_cut_z,focal_cells)


ggsave("withoutNF_variance_accumulation_z.pdf",width=unit(6.2,"in"),height=unit(2.4,"in"))














#total 

ggplot(data=triangle_tree_cut,aes(x=x,y=y))+geom_point(size=0.02,stroke=0.1,col="white")+
geom_line(data=data.frame(x=c(mo_da_link_cut[1,3],mo_da_link_cut[1,5]),y=c(mo_da_link_cut[1,4],mo_da_link_cut[1,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[2,3],mo_da_link_cut[2,5]),y=c(mo_da_link_cut[2,4],mo_da_link_cut[2,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[3,3],mo_da_link_cut[3,5]),y=c(mo_da_link_cut[3,4],mo_da_link_cut[3,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[4,3],mo_da_link_cut[4,5]),y=c(mo_da_link_cut[4,4],mo_da_link_cut[4,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[5,3],mo_da_link_cut[5,5]),y=c(mo_da_link_cut[5,4],mo_da_link_cut[5,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[6,3],mo_da_link_cut[6,5]),y=c(mo_da_link_cut[6,4],mo_da_link_cut[6,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[7,3],mo_da_link_cut[7,5]),y=c(mo_da_link_cut[7,4],mo_da_link_cut[7,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[8,3],mo_da_link_cut[8,5]),y=c(mo_da_link_cut[8,4],mo_da_link_cut[8,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[9,3],mo_da_link_cut[9,5]),y=c(mo_da_link_cut[9,4],mo_da_link_cut[9,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[10,3],mo_da_link_cut[10,5]),y=c(mo_da_link_cut[10,4],mo_da_link_cut[10,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[11,3],mo_da_link_cut[11,5]),y=c(mo_da_link_cut[11,4],mo_da_link_cut[11,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[12,3],mo_da_link_cut[12,5]),y=c(mo_da_link_cut[12,4],mo_da_link_cut[12,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[13,3],mo_da_link_cut[13,5]),y=c(mo_da_link_cut[13,4],mo_da_link_cut[13,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[14,3],mo_da_link_cut[14,5]),y=c(mo_da_link_cut[14,4],mo_da_link_cut[14,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[15,3],mo_da_link_cut[15,5]),y=c(mo_da_link_cut[15,4],mo_da_link_cut[15,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[16,3],mo_da_link_cut[16,5]),y=c(mo_da_link_cut[16,4],mo_da_link_cut[16,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[17,3],mo_da_link_cut[17,5]),y=c(mo_da_link_cut[17,4],mo_da_link_cut[17,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[18,3],mo_da_link_cut[18,5]),y=c(mo_da_link_cut[18,4],mo_da_link_cut[18,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[19,3],mo_da_link_cut[19,5]),y=c(mo_da_link_cut[19,4],mo_da_link_cut[19,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[20,3],mo_da_link_cut[20,5]),y=c(mo_da_link_cut[20,4],mo_da_link_cut[20,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[21,3],mo_da_link_cut[21,5]),y=c(mo_da_link_cut[21,4],mo_da_link_cut[21,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[22,3],mo_da_link_cut[22,5]),y=c(mo_da_link_cut[22,4],mo_da_link_cut[22,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[23,3],mo_da_link_cut[23,5]),y=c(mo_da_link_cut[23,4],mo_da_link_cut[23,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[24,3],mo_da_link_cut[24,5]),y=c(mo_da_link_cut[24,4],mo_da_link_cut[24,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[25,3],mo_da_link_cut[25,5]),y=c(mo_da_link_cut[25,4],mo_da_link_cut[25,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[26,3],mo_da_link_cut[26,5]),y=c(mo_da_link_cut[26,4],mo_da_link_cut[26,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[27,3],mo_da_link_cut[27,5]),y=c(mo_da_link_cut[27,4],mo_da_link_cut[27,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[28,3],mo_da_link_cut[28,5]),y=c(mo_da_link_cut[28,4],mo_da_link_cut[28,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[29,3],mo_da_link_cut[29,5]),y=c(mo_da_link_cut[29,4],mo_da_link_cut[29,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[30,3],mo_da_link_cut[30,5]),y=c(mo_da_link_cut[30,4],mo_da_link_cut[30,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[31,3],mo_da_link_cut[31,5]),y=c(mo_da_link_cut[31,4],mo_da_link_cut[31,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[32,3],mo_da_link_cut[32,5]),y=c(mo_da_link_cut[32,4],mo_da_link_cut[32,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[33,3],mo_da_link_cut[33,5]),y=c(mo_da_link_cut[33,4],mo_da_link_cut[33,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[34,3],mo_da_link_cut[34,5]),y=c(mo_da_link_cut[34,4],mo_da_link_cut[34,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[35,3],mo_da_link_cut[35,5]),y=c(mo_da_link_cut[35,4],mo_da_link_cut[35,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[36,3],mo_da_link_cut[36,5]),y=c(mo_da_link_cut[36,4],mo_da_link_cut[36,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[37,3],mo_da_link_cut[37,5]),y=c(mo_da_link_cut[37,4],mo_da_link_cut[37,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[38,3],mo_da_link_cut[38,5]),y=c(mo_da_link_cut[38,4],mo_da_link_cut[38,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[39,3],mo_da_link_cut[39,5]),y=c(mo_da_link_cut[39,4],mo_da_link_cut[39,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[40,3],mo_da_link_cut[40,5]),y=c(mo_da_link_cut[40,4],mo_da_link_cut[40,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[41,3],mo_da_link_cut[41,5]),y=c(mo_da_link_cut[41,4],mo_da_link_cut[41,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[42,3],mo_da_link_cut[42,5]),y=c(mo_da_link_cut[42,4],mo_da_link_cut[42,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[43,3],mo_da_link_cut[43,5]),y=c(mo_da_link_cut[43,4],mo_da_link_cut[43,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[44,3],mo_da_link_cut[44,5]),y=c(mo_da_link_cut[44,4],mo_da_link_cut[44,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[45,3],mo_da_link_cut[45,5]),y=c(mo_da_link_cut[45,4],mo_da_link_cut[45,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[46,3],mo_da_link_cut[46,5]),y=c(mo_da_link_cut[46,4],mo_da_link_cut[46,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[47,3],mo_da_link_cut[47,5]),y=c(mo_da_link_cut[47,4],mo_da_link_cut[47,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[48,3],mo_da_link_cut[48,5]),y=c(mo_da_link_cut[48,4],mo_da_link_cut[48,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[49,3],mo_da_link_cut[49,5]),y=c(mo_da_link_cut[49,4],mo_da_link_cut[49,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[50,3],mo_da_link_cut[50,5]),y=c(mo_da_link_cut[50,4],mo_da_link_cut[50,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[51,3],mo_da_link_cut[51,5]),y=c(mo_da_link_cut[51,4],mo_da_link_cut[51,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[52,3],mo_da_link_cut[52,5]),y=c(mo_da_link_cut[52,4],mo_da_link_cut[52,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[53,3],mo_da_link_cut[53,5]),y=c(mo_da_link_cut[53,4],mo_da_link_cut[53,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[54,3],mo_da_link_cut[54,5]),y=c(mo_da_link_cut[54,4],mo_da_link_cut[54,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[55,3],mo_da_link_cut[55,5]),y=c(mo_da_link_cut[55,4],mo_da_link_cut[55,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[56,3],mo_da_link_cut[56,5]),y=c(mo_da_link_cut[56,4],mo_da_link_cut[56,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[57,3],mo_da_link_cut[57,5]),y=c(mo_da_link_cut[57,4],mo_da_link_cut[57,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[58,3],mo_da_link_cut[58,5]),y=c(mo_da_link_cut[58,4],mo_da_link_cut[58,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[59,3],mo_da_link_cut[59,5]),y=c(mo_da_link_cut[59,4],mo_da_link_cut[59,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[60,3],mo_da_link_cut[60,5]),y=c(mo_da_link_cut[60,4],mo_da_link_cut[60,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[61,3],mo_da_link_cut[61,5]),y=c(mo_da_link_cut[61,4],mo_da_link_cut[61,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[62,3],mo_da_link_cut[62,5]),y=c(mo_da_link_cut[62,4],mo_da_link_cut[62,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[63,3],mo_da_link_cut[63,5]),y=c(mo_da_link_cut[63,4],mo_da_link_cut[63,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[64,3],mo_da_link_cut[64,5]),y=c(mo_da_link_cut[64,4],mo_da_link_cut[64,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[65,3],mo_da_link_cut[65,5]),y=c(mo_da_link_cut[65,4],mo_da_link_cut[65,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[66,3],mo_da_link_cut[66,5]),y=c(mo_da_link_cut[66,4],mo_da_link_cut[66,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[67,3],mo_da_link_cut[67,5]),y=c(mo_da_link_cut[67,4],mo_da_link_cut[67,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[68,3],mo_da_link_cut[68,5]),y=c(mo_da_link_cut[68,4],mo_da_link_cut[68,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[69,3],mo_da_link_cut[69,5]),y=c(mo_da_link_cut[69,4],mo_da_link_cut[69,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[70,3],mo_da_link_cut[70,5]),y=c(mo_da_link_cut[70,4],mo_da_link_cut[70,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[71,3],mo_da_link_cut[71,5]),y=c(mo_da_link_cut[71,4],mo_da_link_cut[71,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[72,3],mo_da_link_cut[72,5]),y=c(mo_da_link_cut[72,4],mo_da_link_cut[72,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[73,3],mo_da_link_cut[73,5]),y=c(mo_da_link_cut[73,4],mo_da_link_cut[73,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[74,3],mo_da_link_cut[74,5]),y=c(mo_da_link_cut[74,4],mo_da_link_cut[74,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[75,3],mo_da_link_cut[75,5]),y=c(mo_da_link_cut[75,4],mo_da_link_cut[75,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[76,3],mo_da_link_cut[76,5]),y=c(mo_da_link_cut[76,4],mo_da_link_cut[76,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[77,3],mo_da_link_cut[77,5]),y=c(mo_da_link_cut[77,4],mo_da_link_cut[77,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[78,3],mo_da_link_cut[78,5]),y=c(mo_da_link_cut[78,4],mo_da_link_cut[78,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[79,3],mo_da_link_cut[79,5]),y=c(mo_da_link_cut[79,4],mo_da_link_cut[79,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[80,3],mo_da_link_cut[80,5]),y=c(mo_da_link_cut[80,4],mo_da_link_cut[80,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[81,3],mo_da_link_cut[81,5]),y=c(mo_da_link_cut[81,4],mo_da_link_cut[81,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[82,3],mo_da_link_cut[82,5]),y=c(mo_da_link_cut[82,4],mo_da_link_cut[82,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[83,3],mo_da_link_cut[83,5]),y=c(mo_da_link_cut[83,4],mo_da_link_cut[83,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[84,3],mo_da_link_cut[84,5]),y=c(mo_da_link_cut[84,4],mo_da_link_cut[84,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[85,3],mo_da_link_cut[85,5]),y=c(mo_da_link_cut[85,4],mo_da_link_cut[85,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[86,3],mo_da_link_cut[86,5]),y=c(mo_da_link_cut[86,4],mo_da_link_cut[86,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[87,3],mo_da_link_cut[87,5]),y=c(mo_da_link_cut[87,4],mo_da_link_cut[87,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[88,3],mo_da_link_cut[88,5]),y=c(mo_da_link_cut[88,4],mo_da_link_cut[88,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[89,3],mo_da_link_cut[89,5]),y=c(mo_da_link_cut[89,4],mo_da_link_cut[89,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[90,3],mo_da_link_cut[90,5]),y=c(mo_da_link_cut[90,4],mo_da_link_cut[90,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[91,3],mo_da_link_cut[91,5]),y=c(mo_da_link_cut[91,4],mo_da_link_cut[91,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[92,3],mo_da_link_cut[92,5]),y=c(mo_da_link_cut[92,4],mo_da_link_cut[92,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[93,3],mo_da_link_cut[93,5]),y=c(mo_da_link_cut[93,4],mo_da_link_cut[93,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[94,3],mo_da_link_cut[94,5]),y=c(mo_da_link_cut[94,4],mo_da_link_cut[94,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[95,3],mo_da_link_cut[95,5]),y=c(mo_da_link_cut[95,4],mo_da_link_cut[95,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[96,3],mo_da_link_cut[96,5]),y=c(mo_da_link_cut[96,4],mo_da_link_cut[96,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[97,3],mo_da_link_cut[97,5]),y=c(mo_da_link_cut[97,4],mo_da_link_cut[97,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[98,3],mo_da_link_cut[98,5]),y=c(mo_da_link_cut[98,4],mo_da_link_cut[98,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[99,3],mo_da_link_cut[99,5]),y=c(mo_da_link_cut[99,4],mo_da_link_cut[99,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[100,3],mo_da_link_cut[100,5]),y=c(mo_da_link_cut[100,4],mo_da_link_cut[100,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[101,3],mo_da_link_cut[101,5]),y=c(mo_da_link_cut[101,4],mo_da_link_cut[101,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[102,3],mo_da_link_cut[102,5]),y=c(mo_da_link_cut[102,4],mo_da_link_cut[102,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[103,3],mo_da_link_cut[103,5]),y=c(mo_da_link_cut[103,4],mo_da_link_cut[103,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[104,3],mo_da_link_cut[104,5]),y=c(mo_da_link_cut[104,4],mo_da_link_cut[104,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[105,3],mo_da_link_cut[105,5]),y=c(mo_da_link_cut[105,4],mo_da_link_cut[105,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[106,3],mo_da_link_cut[106,5]),y=c(mo_da_link_cut[106,4],mo_da_link_cut[106,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[107,3],mo_da_link_cut[107,5]),y=c(mo_da_link_cut[107,4],mo_da_link_cut[107,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[108,3],mo_da_link_cut[108,5]),y=c(mo_da_link_cut[108,4],mo_da_link_cut[108,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[109,3],mo_da_link_cut[109,5]),y=c(mo_da_link_cut[109,4],mo_da_link_cut[109,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[110,3],mo_da_link_cut[110,5]),y=c(mo_da_link_cut[110,4],mo_da_link_cut[110,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[111,3],mo_da_link_cut[111,5]),y=c(mo_da_link_cut[111,4],mo_da_link_cut[111,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[112,3],mo_da_link_cut[112,5]),y=c(mo_da_link_cut[112,4],mo_da_link_cut[112,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[113,3],mo_da_link_cut[113,5]),y=c(mo_da_link_cut[113,4],mo_da_link_cut[113,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[114,3],mo_da_link_cut[114,5]),y=c(mo_da_link_cut[114,4],mo_da_link_cut[114,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[115,3],mo_da_link_cut[115,5]),y=c(mo_da_link_cut[115,4],mo_da_link_cut[115,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[116,3],mo_da_link_cut[116,5]),y=c(mo_da_link_cut[116,4],mo_da_link_cut[116,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[117,3],mo_da_link_cut[117,5]),y=c(mo_da_link_cut[117,4],mo_da_link_cut[117,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[118,3],mo_da_link_cut[118,5]),y=c(mo_da_link_cut[118,4],mo_da_link_cut[118,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[119,3],mo_da_link_cut[119,5]),y=c(mo_da_link_cut[119,4],mo_da_link_cut[119,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[120,3],mo_da_link_cut[120,5]),y=c(mo_da_link_cut[120,4],mo_da_link_cut[120,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[121,3],mo_da_link_cut[121,5]),y=c(mo_da_link_cut[121,4],mo_da_link_cut[121,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[122,3],mo_da_link_cut[122,5]),y=c(mo_da_link_cut[122,4],mo_da_link_cut[122,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[123,3],mo_da_link_cut[123,5]),y=c(mo_da_link_cut[123,4],mo_da_link_cut[123,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[124,3],mo_da_link_cut[124,5]),y=c(mo_da_link_cut[124,4],mo_da_link_cut[124,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[125,3],mo_da_link_cut[125,5]),y=c(mo_da_link_cut[125,4],mo_da_link_cut[125,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[126,3],mo_da_link_cut[126,5]),y=c(mo_da_link_cut[126,4],mo_da_link_cut[126,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[127,3],mo_da_link_cut[127,5]),y=c(mo_da_link_cut[127,4],mo_da_link_cut[127,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[128,3],mo_da_link_cut[128,5]),y=c(mo_da_link_cut[128,4],mo_da_link_cut[128,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[129,3],mo_da_link_cut[129,5]),y=c(mo_da_link_cut[129,4],mo_da_link_cut[129,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[130,3],mo_da_link_cut[130,5]),y=c(mo_da_link_cut[130,4],mo_da_link_cut[130,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[131,3],mo_da_link_cut[131,5]),y=c(mo_da_link_cut[131,4],mo_da_link_cut[131,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[132,3],mo_da_link_cut[132,5]),y=c(mo_da_link_cut[132,4],mo_da_link_cut[132,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[133,3],mo_da_link_cut[133,5]),y=c(mo_da_link_cut[133,4],mo_da_link_cut[133,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[134,3],mo_da_link_cut[134,5]),y=c(mo_da_link_cut[134,4],mo_da_link_cut[134,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[135,3],mo_da_link_cut[135,5]),y=c(mo_da_link_cut[135,4],mo_da_link_cut[135,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[136,3],mo_da_link_cut[136,5]),y=c(mo_da_link_cut[136,4],mo_da_link_cut[136,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[137,3],mo_da_link_cut[137,5]),y=c(mo_da_link_cut[137,4],mo_da_link_cut[137,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[138,3],mo_da_link_cut[138,5]),y=c(mo_da_link_cut[138,4],mo_da_link_cut[138,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[139,3],mo_da_link_cut[139,5]),y=c(mo_da_link_cut[139,4],mo_da_link_cut[139,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[140,3],mo_da_link_cut[140,5]),y=c(mo_da_link_cut[140,4],mo_da_link_cut[140,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[141,3],mo_da_link_cut[141,5]),y=c(mo_da_link_cut[141,4],mo_da_link_cut[141,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[142,3],mo_da_link_cut[142,5]),y=c(mo_da_link_cut[142,4],mo_da_link_cut[142,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[143,3],mo_da_link_cut[143,5]),y=c(mo_da_link_cut[143,4],mo_da_link_cut[143,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[144,3],mo_da_link_cut[144,5]),y=c(mo_da_link_cut[144,4],mo_da_link_cut[144,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[145,3],mo_da_link_cut[145,5]),y=c(mo_da_link_cut[145,4],mo_da_link_cut[145,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[146,3],mo_da_link_cut[146,5]),y=c(mo_da_link_cut[146,4],mo_da_link_cut[146,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[147,3],mo_da_link_cut[147,5]),y=c(mo_da_link_cut[147,4],mo_da_link_cut[147,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[148,3],mo_da_link_cut[148,5]),y=c(mo_da_link_cut[148,4],mo_da_link_cut[148,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[149,3],mo_da_link_cut[149,5]),y=c(mo_da_link_cut[149,4],mo_da_link_cut[149,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[150,3],mo_da_link_cut[150,5]),y=c(mo_da_link_cut[150,4],mo_da_link_cut[150,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[151,3],mo_da_link_cut[151,5]),y=c(mo_da_link_cut[151,4],mo_da_link_cut[151,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[152,3],mo_da_link_cut[152,5]),y=c(mo_da_link_cut[152,4],mo_da_link_cut[152,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[153,3],mo_da_link_cut[153,5]),y=c(mo_da_link_cut[153,4],mo_da_link_cut[153,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[154,3],mo_da_link_cut[154,5]),y=c(mo_da_link_cut[154,4],mo_da_link_cut[154,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[155,3],mo_da_link_cut[155,5]),y=c(mo_da_link_cut[155,4],mo_da_link_cut[155,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[156,3],mo_da_link_cut[156,5]),y=c(mo_da_link_cut[156,4],mo_da_link_cut[156,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[157,3],mo_da_link_cut[157,5]),y=c(mo_da_link_cut[157,4],mo_da_link_cut[157,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[158,3],mo_da_link_cut[158,5]),y=c(mo_da_link_cut[158,4],mo_da_link_cut[158,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[159,3],mo_da_link_cut[159,5]),y=c(mo_da_link_cut[159,4],mo_da_link_cut[159,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[160,3],mo_da_link_cut[160,5]),y=c(mo_da_link_cut[160,4],mo_da_link_cut[160,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[161,3],mo_da_link_cut[161,5]),y=c(mo_da_link_cut[161,4],mo_da_link_cut[161,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[162,3],mo_da_link_cut[162,5]),y=c(mo_da_link_cut[162,4],mo_da_link_cut[162,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[163,3],mo_da_link_cut[163,5]),y=c(mo_da_link_cut[163,4],mo_da_link_cut[163,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[164,3],mo_da_link_cut[164,5]),y=c(mo_da_link_cut[164,4],mo_da_link_cut[164,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[165,3],mo_da_link_cut[165,5]),y=c(mo_da_link_cut[165,4],mo_da_link_cut[165,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[166,3],mo_da_link_cut[166,5]),y=c(mo_da_link_cut[166,4],mo_da_link_cut[166,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[167,3],mo_da_link_cut[167,5]),y=c(mo_da_link_cut[167,4],mo_da_link_cut[167,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[168,3],mo_da_link_cut[168,5]),y=c(mo_da_link_cut[168,4],mo_da_link_cut[168,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[169,3],mo_da_link_cut[169,5]),y=c(mo_da_link_cut[169,4],mo_da_link_cut[169,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[170,3],mo_da_link_cut[170,5]),y=c(mo_da_link_cut[170,4],mo_da_link_cut[170,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[171,3],mo_da_link_cut[171,5]),y=c(mo_da_link_cut[171,4],mo_da_link_cut[171,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[172,3],mo_da_link_cut[172,5]),y=c(mo_da_link_cut[172,4],mo_da_link_cut[172,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[173,3],mo_da_link_cut[173,5]),y=c(mo_da_link_cut[173,4],mo_da_link_cut[173,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[174,3],mo_da_link_cut[174,5]),y=c(mo_da_link_cut[174,4],mo_da_link_cut[174,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[175,3],mo_da_link_cut[175,5]),y=c(mo_da_link_cut[175,4],mo_da_link_cut[175,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[176,3],mo_da_link_cut[176,5]),y=c(mo_da_link_cut[176,4],mo_da_link_cut[176,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[177,3],mo_da_link_cut[177,5]),y=c(mo_da_link_cut[177,4],mo_da_link_cut[177,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[178,3],mo_da_link_cut[178,5]),y=c(mo_da_link_cut[178,4],mo_da_link_cut[178,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[179,3],mo_da_link_cut[179,5]),y=c(mo_da_link_cut[179,4],mo_da_link_cut[179,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[180,3],mo_da_link_cut[180,5]),y=c(mo_da_link_cut[180,4],mo_da_link_cut[180,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[181,3],mo_da_link_cut[181,5]),y=c(mo_da_link_cut[181,4],mo_da_link_cut[181,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[182,3],mo_da_link_cut[182,5]),y=c(mo_da_link_cut[182,4],mo_da_link_cut[182,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[183,3],mo_da_link_cut[183,5]),y=c(mo_da_link_cut[183,4],mo_da_link_cut[183,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[184,3],mo_da_link_cut[184,5]),y=c(mo_da_link_cut[184,4],mo_da_link_cut[184,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[185,3],mo_da_link_cut[185,5]),y=c(mo_da_link_cut[185,4],mo_da_link_cut[185,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[186,3],mo_da_link_cut[186,5]),y=c(mo_da_link_cut[186,4],mo_da_link_cut[186,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[187,3],mo_da_link_cut[187,5]),y=c(mo_da_link_cut[187,4],mo_da_link_cut[187,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[188,3],mo_da_link_cut[188,5]),y=c(mo_da_link_cut[188,4],mo_da_link_cut[188,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[189,3],mo_da_link_cut[189,5]),y=c(mo_da_link_cut[189,4],mo_da_link_cut[189,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[190,3],mo_da_link_cut[190,5]),y=c(mo_da_link_cut[190,4],mo_da_link_cut[190,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[191,3],mo_da_link_cut[191,5]),y=c(mo_da_link_cut[191,4],mo_da_link_cut[191,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[192,3],mo_da_link_cut[192,5]),y=c(mo_da_link_cut[192,4],mo_da_link_cut[192,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[193,3],mo_da_link_cut[193,5]),y=c(mo_da_link_cut[193,4],mo_da_link_cut[193,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[194,3],mo_da_link_cut[194,5]),y=c(mo_da_link_cut[194,4],mo_da_link_cut[194,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[195,3],mo_da_link_cut[195,5]),y=c(mo_da_link_cut[195,4],mo_da_link_cut[195,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[196,3],mo_da_link_cut[196,5]),y=c(mo_da_link_cut[196,4],mo_da_link_cut[196,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[197,3],mo_da_link_cut[197,5]),y=c(mo_da_link_cut[197,4],mo_da_link_cut[197,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[198,3],mo_da_link_cut[198,5]),y=c(mo_da_link_cut[198,4],mo_da_link_cut[198,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[199,3],mo_da_link_cut[199,5]),y=c(mo_da_link_cut[199,4],mo_da_link_cut[199,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[200,3],mo_da_link_cut[200,5]),y=c(mo_da_link_cut[200,4],mo_da_link_cut[200,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[201,3],mo_da_link_cut[201,5]),y=c(mo_da_link_cut[201,4],mo_da_link_cut[201,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[202,3],mo_da_link_cut[202,5]),y=c(mo_da_link_cut[202,4],mo_da_link_cut[202,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[203,3],mo_da_link_cut[203,5]),y=c(mo_da_link_cut[203,4],mo_da_link_cut[203,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[204,3],mo_da_link_cut[204,5]),y=c(mo_da_link_cut[204,4],mo_da_link_cut[204,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[205,3],mo_da_link_cut[205,5]),y=c(mo_da_link_cut[205,4],mo_da_link_cut[205,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[206,3],mo_da_link_cut[206,5]),y=c(mo_da_link_cut[206,4],mo_da_link_cut[206,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[207,3],mo_da_link_cut[207,5]),y=c(mo_da_link_cut[207,4],mo_da_link_cut[207,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[208,3],mo_da_link_cut[208,5]),y=c(mo_da_link_cut[208,4],mo_da_link_cut[208,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[209,3],mo_da_link_cut[209,5]),y=c(mo_da_link_cut[209,4],mo_da_link_cut[209,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[210,3],mo_da_link_cut[210,5]),y=c(mo_da_link_cut[210,4],mo_da_link_cut[210,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[211,3],mo_da_link_cut[211,5]),y=c(mo_da_link_cut[211,4],mo_da_link_cut[211,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[212,3],mo_da_link_cut[212,5]),y=c(mo_da_link_cut[212,4],mo_da_link_cut[212,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[213,3],mo_da_link_cut[213,5]),y=c(mo_da_link_cut[213,4],mo_da_link_cut[213,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[214,3],mo_da_link_cut[214,5]),y=c(mo_da_link_cut[214,4],mo_da_link_cut[214,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[215,3],mo_da_link_cut[215,5]),y=c(mo_da_link_cut[215,4],mo_da_link_cut[215,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[216,3],mo_da_link_cut[216,5]),y=c(mo_da_link_cut[216,4],mo_da_link_cut[216,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[217,3],mo_da_link_cut[217,5]),y=c(mo_da_link_cut[217,4],mo_da_link_cut[217,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[218,3],mo_da_link_cut[218,5]),y=c(mo_da_link_cut[218,4],mo_da_link_cut[218,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[219,3],mo_da_link_cut[219,5]),y=c(mo_da_link_cut[219,4],mo_da_link_cut[219,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[220,3],mo_da_link_cut[220,5]),y=c(mo_da_link_cut[220,4],mo_da_link_cut[220,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[221,3],mo_da_link_cut[221,5]),y=c(mo_da_link_cut[221,4],mo_da_link_cut[221,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[222,3],mo_da_link_cut[222,5]),y=c(mo_da_link_cut[222,4],mo_da_link_cut[222,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[223,3],mo_da_link_cut[223,5]),y=c(mo_da_link_cut[223,4],mo_da_link_cut[223,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[224,3],mo_da_link_cut[224,5]),y=c(mo_da_link_cut[224,4],mo_da_link_cut[224,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[225,3],mo_da_link_cut[225,5]),y=c(mo_da_link_cut[225,4],mo_da_link_cut[225,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[226,3],mo_da_link_cut[226,5]),y=c(mo_da_link_cut[226,4],mo_da_link_cut[226,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[227,3],mo_da_link_cut[227,5]),y=c(mo_da_link_cut[227,4],mo_da_link_cut[227,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[228,3],mo_da_link_cut[228,5]),y=c(mo_da_link_cut[228,4],mo_da_link_cut[228,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[229,3],mo_da_link_cut[229,5]),y=c(mo_da_link_cut[229,4],mo_da_link_cut[229,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[230,3],mo_da_link_cut[230,5]),y=c(mo_da_link_cut[230,4],mo_da_link_cut[230,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[231,3],mo_da_link_cut[231,5]),y=c(mo_da_link_cut[231,4],mo_da_link_cut[231,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[232,3],mo_da_link_cut[232,5]),y=c(mo_da_link_cut[232,4],mo_da_link_cut[232,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[233,3],mo_da_link_cut[233,5]),y=c(mo_da_link_cut[233,4],mo_da_link_cut[233,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[234,3],mo_da_link_cut[234,5]),y=c(mo_da_link_cut[234,4],mo_da_link_cut[234,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[235,3],mo_da_link_cut[235,5]),y=c(mo_da_link_cut[235,4],mo_da_link_cut[235,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[236,3],mo_da_link_cut[236,5]),y=c(mo_da_link_cut[236,4],mo_da_link_cut[236,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[237,3],mo_da_link_cut[237,5]),y=c(mo_da_link_cut[237,4],mo_da_link_cut[237,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[238,3],mo_da_link_cut[238,5]),y=c(mo_da_link_cut[238,4],mo_da_link_cut[238,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[239,3],mo_da_link_cut[239,5]),y=c(mo_da_link_cut[239,4],mo_da_link_cut[239,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[240,3],mo_da_link_cut[240,5]),y=c(mo_da_link_cut[240,4],mo_da_link_cut[240,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[241,3],mo_da_link_cut[241,5]),y=c(mo_da_link_cut[241,4],mo_da_link_cut[241,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[242,3],mo_da_link_cut[242,5]),y=c(mo_da_link_cut[242,4],mo_da_link_cut[242,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[243,3],mo_da_link_cut[243,5]),y=c(mo_da_link_cut[243,4],mo_da_link_cut[243,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[244,3],mo_da_link_cut[244,5]),y=c(mo_da_link_cut[244,4],mo_da_link_cut[244,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[245,3],mo_da_link_cut[245,5]),y=c(mo_da_link_cut[245,4],mo_da_link_cut[245,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[246,3],mo_da_link_cut[246,5]),y=c(mo_da_link_cut[246,4],mo_da_link_cut[246,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[247,3],mo_da_link_cut[247,5]),y=c(mo_da_link_cut[247,4],mo_da_link_cut[247,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[248,3],mo_da_link_cut[248,5]),y=c(mo_da_link_cut[248,4],mo_da_link_cut[248,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[249,3],mo_da_link_cut[249,5]),y=c(mo_da_link_cut[249,4],mo_da_link_cut[249,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[250,3],mo_da_link_cut[250,5]),y=c(mo_da_link_cut[250,4],mo_da_link_cut[250,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[251,3],mo_da_link_cut[251,5]),y=c(mo_da_link_cut[251,4],mo_da_link_cut[251,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[252,3],mo_da_link_cut[252,5]),y=c(mo_da_link_cut[252,4],mo_da_link_cut[252,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[253,3],mo_da_link_cut[253,5]),y=c(mo_da_link_cut[253,4],mo_da_link_cut[253,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[254,3],mo_da_link_cut[254,5]),y=c(mo_da_link_cut[254,4],mo_da_link_cut[254,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[255,3],mo_da_link_cut[255,5]),y=c(mo_da_link_cut[255,4],mo_da_link_cut[255,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[256,3],mo_da_link_cut[256,5]),y=c(mo_da_link_cut[256,4],mo_da_link_cut[256,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[257,3],mo_da_link_cut[257,5]),y=c(mo_da_link_cut[257,4],mo_da_link_cut[257,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[258,3],mo_da_link_cut[258,5]),y=c(mo_da_link_cut[258,4],mo_da_link_cut[258,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[259,3],mo_da_link_cut[259,5]),y=c(mo_da_link_cut[259,4],mo_da_link_cut[259,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[260,3],mo_da_link_cut[260,5]),y=c(mo_da_link_cut[260,4],mo_da_link_cut[260,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[261,3],mo_da_link_cut[261,5]),y=c(mo_da_link_cut[261,4],mo_da_link_cut[261,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[262,3],mo_da_link_cut[262,5]),y=c(mo_da_link_cut[262,4],mo_da_link_cut[262,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[263,3],mo_da_link_cut[263,5]),y=c(mo_da_link_cut[263,4],mo_da_link_cut[263,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[264,3],mo_da_link_cut[264,5]),y=c(mo_da_link_cut[264,4],mo_da_link_cut[264,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[265,3],mo_da_link_cut[265,5]),y=c(mo_da_link_cut[265,4],mo_da_link_cut[265,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[266,3],mo_da_link_cut[266,5]),y=c(mo_da_link_cut[266,4],mo_da_link_cut[266,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[267,3],mo_da_link_cut[267,5]),y=c(mo_da_link_cut[267,4],mo_da_link_cut[267,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[268,3],mo_da_link_cut[268,5]),y=c(mo_da_link_cut[268,4],mo_da_link_cut[268,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[269,3],mo_da_link_cut[269,5]),y=c(mo_da_link_cut[269,4],mo_da_link_cut[269,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[270,3],mo_da_link_cut[270,5]),y=c(mo_da_link_cut[270,4],mo_da_link_cut[270,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[271,3],mo_da_link_cut[271,5]),y=c(mo_da_link_cut[271,4],mo_da_link_cut[271,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[272,3],mo_da_link_cut[272,5]),y=c(mo_da_link_cut[272,4],mo_da_link_cut[272,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[273,3],mo_da_link_cut[273,5]),y=c(mo_da_link_cut[273,4],mo_da_link_cut[273,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[274,3],mo_da_link_cut[274,5]),y=c(mo_da_link_cut[274,4],mo_da_link_cut[274,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[275,3],mo_da_link_cut[275,5]),y=c(mo_da_link_cut[275,4],mo_da_link_cut[275,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[276,3],mo_da_link_cut[276,5]),y=c(mo_da_link_cut[276,4],mo_da_link_cut[276,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[277,3],mo_da_link_cut[277,5]),y=c(mo_da_link_cut[277,4],mo_da_link_cut[277,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[278,3],mo_da_link_cut[278,5]),y=c(mo_da_link_cut[278,4],mo_da_link_cut[278,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[279,3],mo_da_link_cut[279,5]),y=c(mo_da_link_cut[279,4],mo_da_link_cut[279,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[280,3],mo_da_link_cut[280,5]),y=c(mo_da_link_cut[280,4],mo_da_link_cut[280,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[281,3],mo_da_link_cut[281,5]),y=c(mo_da_link_cut[281,4],mo_da_link_cut[281,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[282,3],mo_da_link_cut[282,5]),y=c(mo_da_link_cut[282,4],mo_da_link_cut[282,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[283,3],mo_da_link_cut[283,5]),y=c(mo_da_link_cut[283,4],mo_da_link_cut[283,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[284,3],mo_da_link_cut[284,5]),y=c(mo_da_link_cut[284,4],mo_da_link_cut[284,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[285,3],mo_da_link_cut[285,5]),y=c(mo_da_link_cut[285,4],mo_da_link_cut[285,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[286,3],mo_da_link_cut[286,5]),y=c(mo_da_link_cut[286,4],mo_da_link_cut[286,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[287,3],mo_da_link_cut[287,5]),y=c(mo_da_link_cut[287,4],mo_da_link_cut[287,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[288,3],mo_da_link_cut[288,5]),y=c(mo_da_link_cut[288,4],mo_da_link_cut[288,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[289,3],mo_da_link_cut[289,5]),y=c(mo_da_link_cut[289,4],mo_da_link_cut[289,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[290,3],mo_da_link_cut[290,5]),y=c(mo_da_link_cut[290,4],mo_da_link_cut[290,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[291,3],mo_da_link_cut[291,5]),y=c(mo_da_link_cut[291,4],mo_da_link_cut[291,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[292,3],mo_da_link_cut[292,5]),y=c(mo_da_link_cut[292,4],mo_da_link_cut[292,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[293,3],mo_da_link_cut[293,5]),y=c(mo_da_link_cut[293,4],mo_da_link_cut[293,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[294,3],mo_da_link_cut[294,5]),y=c(mo_da_link_cut[294,4],mo_da_link_cut[294,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[295,3],mo_da_link_cut[295,5]),y=c(mo_da_link_cut[295,4],mo_da_link_cut[295,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[296,3],mo_da_link_cut[296,5]),y=c(mo_da_link_cut[296,4],mo_da_link_cut[296,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[297,3],mo_da_link_cut[297,5]),y=c(mo_da_link_cut[297,4],mo_da_link_cut[297,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[298,3],mo_da_link_cut[298,5]),y=c(mo_da_link_cut[298,4],mo_da_link_cut[298,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[299,3],mo_da_link_cut[299,5]),y=c(mo_da_link_cut[299,4],mo_da_link_cut[299,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[300,3],mo_da_link_cut[300,5]),y=c(mo_da_link_cut[300,4],mo_da_link_cut[300,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[301,3],mo_da_link_cut[301,5]),y=c(mo_da_link_cut[301,4],mo_da_link_cut[301,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[302,3],mo_da_link_cut[302,5]),y=c(mo_da_link_cut[302,4],mo_da_link_cut[302,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[303,3],mo_da_link_cut[303,5]),y=c(mo_da_link_cut[303,4],mo_da_link_cut[303,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[304,3],mo_da_link_cut[304,5]),y=c(mo_da_link_cut[304,4],mo_da_link_cut[304,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[305,3],mo_da_link_cut[305,5]),y=c(mo_da_link_cut[305,4],mo_da_link_cut[305,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[306,3],mo_da_link_cut[306,5]),y=c(mo_da_link_cut[306,4],mo_da_link_cut[306,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[307,3],mo_da_link_cut[307,5]),y=c(mo_da_link_cut[307,4],mo_da_link_cut[307,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[308,3],mo_da_link_cut[308,5]),y=c(mo_da_link_cut[308,4],mo_da_link_cut[308,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[309,3],mo_da_link_cut[309,5]),y=c(mo_da_link_cut[309,4],mo_da_link_cut[309,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[310,3],mo_da_link_cut[310,5]),y=c(mo_da_link_cut[310,4],mo_da_link_cut[310,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[311,3],mo_da_link_cut[311,5]),y=c(mo_da_link_cut[311,4],mo_da_link_cut[311,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[312,3],mo_da_link_cut[312,5]),y=c(mo_da_link_cut[312,4],mo_da_link_cut[312,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[313,3],mo_da_link_cut[313,5]),y=c(mo_da_link_cut[313,4],mo_da_link_cut[313,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[314,3],mo_da_link_cut[314,5]),y=c(mo_da_link_cut[314,4],mo_da_link_cut[314,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[315,3],mo_da_link_cut[315,5]),y=c(mo_da_link_cut[315,4],mo_da_link_cut[315,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[316,3],mo_da_link_cut[316,5]),y=c(mo_da_link_cut[316,4],mo_da_link_cut[316,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[317,3],mo_da_link_cut[317,5]),y=c(mo_da_link_cut[317,4],mo_da_link_cut[317,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[318,3],mo_da_link_cut[318,5]),y=c(mo_da_link_cut[318,4],mo_da_link_cut[318,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[319,3],mo_da_link_cut[319,5]),y=c(mo_da_link_cut[319,4],mo_da_link_cut[319,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[320,3],mo_da_link_cut[320,5]),y=c(mo_da_link_cut[320,4],mo_da_link_cut[320,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[321,3],mo_da_link_cut[321,5]),y=c(mo_da_link_cut[321,4],mo_da_link_cut[321,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[322,3],mo_da_link_cut[322,5]),y=c(mo_da_link_cut[322,4],mo_da_link_cut[322,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[323,3],mo_da_link_cut[323,5]),y=c(mo_da_link_cut[323,4],mo_da_link_cut[323,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[324,3],mo_da_link_cut[324,5]),y=c(mo_da_link_cut[324,4],mo_da_link_cut[324,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[325,3],mo_da_link_cut[325,5]),y=c(mo_da_link_cut[325,4],mo_da_link_cut[325,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[326,3],mo_da_link_cut[326,5]),y=c(mo_da_link_cut[326,4],mo_da_link_cut[326,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[327,3],mo_da_link_cut[327,5]),y=c(mo_da_link_cut[327,4],mo_da_link_cut[327,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[328,3],mo_da_link_cut[328,5]),y=c(mo_da_link_cut[328,4],mo_da_link_cut[328,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[329,3],mo_da_link_cut[329,5]),y=c(mo_da_link_cut[329,4],mo_da_link_cut[329,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[330,3],mo_da_link_cut[330,5]),y=c(mo_da_link_cut[330,4],mo_da_link_cut[330,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[331,3],mo_da_link_cut[331,5]),y=c(mo_da_link_cut[331,4],mo_da_link_cut[331,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[332,3],mo_da_link_cut[332,5]),y=c(mo_da_link_cut[332,4],mo_da_link_cut[332,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[333,3],mo_da_link_cut[333,5]),y=c(mo_da_link_cut[333,4],mo_da_link_cut[333,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[334,3],mo_da_link_cut[334,5]),y=c(mo_da_link_cut[334,4],mo_da_link_cut[334,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[335,3],mo_da_link_cut[335,5]),y=c(mo_da_link_cut[335,4],mo_da_link_cut[335,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[336,3],mo_da_link_cut[336,5]),y=c(mo_da_link_cut[336,4],mo_da_link_cut[336,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[337,3],mo_da_link_cut[337,5]),y=c(mo_da_link_cut[337,4],mo_da_link_cut[337,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[338,3],mo_da_link_cut[338,5]),y=c(mo_da_link_cut[338,4],mo_da_link_cut[338,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[339,3],mo_da_link_cut[339,5]),y=c(mo_da_link_cut[339,4],mo_da_link_cut[339,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[340,3],mo_da_link_cut[340,5]),y=c(mo_da_link_cut[340,4],mo_da_link_cut[340,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[341,3],mo_da_link_cut[341,5]),y=c(mo_da_link_cut[341,4],mo_da_link_cut[341,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[342,3],mo_da_link_cut[342,5]),y=c(mo_da_link_cut[342,4],mo_da_link_cut[342,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[343,3],mo_da_link_cut[343,5]),y=c(mo_da_link_cut[343,4],mo_da_link_cut[343,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[344,3],mo_da_link_cut[344,5]),y=c(mo_da_link_cut[344,4],mo_da_link_cut[344,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[345,3],mo_da_link_cut[345,5]),y=c(mo_da_link_cut[345,4],mo_da_link_cut[345,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[346,3],mo_da_link_cut[346,5]),y=c(mo_da_link_cut[346,4],mo_da_link_cut[346,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[347,3],mo_da_link_cut[347,5]),y=c(mo_da_link_cut[347,4],mo_da_link_cut[347,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[348,3],mo_da_link_cut[348,5]),y=c(mo_da_link_cut[348,4],mo_da_link_cut[348,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[349,3],mo_da_link_cut[349,5]),y=c(mo_da_link_cut[349,4],mo_da_link_cut[349,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[350,3],mo_da_link_cut[350,5]),y=c(mo_da_link_cut[350,4],mo_da_link_cut[350,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[351,3],mo_da_link_cut[351,5]),y=c(mo_da_link_cut[351,4],mo_da_link_cut[351,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[352,3],mo_da_link_cut[352,5]),y=c(mo_da_link_cut[352,4],mo_da_link_cut[352,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[353,3],mo_da_link_cut[353,5]),y=c(mo_da_link_cut[353,4],mo_da_link_cut[353,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[354,3],mo_da_link_cut[354,5]),y=c(mo_da_link_cut[354,4],mo_da_link_cut[354,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[355,3],mo_da_link_cut[355,5]),y=c(mo_da_link_cut[355,4],mo_da_link_cut[355,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[356,3],mo_da_link_cut[356,5]),y=c(mo_da_link_cut[356,4],mo_da_link_cut[356,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[357,3],mo_da_link_cut[357,5]),y=c(mo_da_link_cut[357,4],mo_da_link_cut[357,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[358,3],mo_da_link_cut[358,5]),y=c(mo_da_link_cut[358,4],mo_da_link_cut[358,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[359,3],mo_da_link_cut[359,5]),y=c(mo_da_link_cut[359,4],mo_da_link_cut[359,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[360,3],mo_da_link_cut[360,5]),y=c(mo_da_link_cut[360,4],mo_da_link_cut[360,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[361,3],mo_da_link_cut[361,5]),y=c(mo_da_link_cut[361,4],mo_da_link_cut[361,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[362,3],mo_da_link_cut[362,5]),y=c(mo_da_link_cut[362,4],mo_da_link_cut[362,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[363,3],mo_da_link_cut[363,5]),y=c(mo_da_link_cut[363,4],mo_da_link_cut[363,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[364,3],mo_da_link_cut[364,5]),y=c(mo_da_link_cut[364,4],mo_da_link_cut[364,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[365,3],mo_da_link_cut[365,5]),y=c(mo_da_link_cut[365,4],mo_da_link_cut[365,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[366,3],mo_da_link_cut[366,5]),y=c(mo_da_link_cut[366,4],mo_da_link_cut[366,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[367,3],mo_da_link_cut[367,5]),y=c(mo_da_link_cut[367,4],mo_da_link_cut[367,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[368,3],mo_da_link_cut[368,5]),y=c(mo_da_link_cut[368,4],mo_da_link_cut[368,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[369,3],mo_da_link_cut[369,5]),y=c(mo_da_link_cut[369,4],mo_da_link_cut[369,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[370,3],mo_da_link_cut[370,5]),y=c(mo_da_link_cut[370,4],mo_da_link_cut[370,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[371,3],mo_da_link_cut[371,5]),y=c(mo_da_link_cut[371,4],mo_da_link_cut[371,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[372,3],mo_da_link_cut[372,5]),y=c(mo_da_link_cut[372,4],mo_da_link_cut[372,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[373,3],mo_da_link_cut[373,5]),y=c(mo_da_link_cut[373,4],mo_da_link_cut[373,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[374,3],mo_da_link_cut[374,5]),y=c(mo_da_link_cut[374,4],mo_da_link_cut[374,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[375,3],mo_da_link_cut[375,5]),y=c(mo_da_link_cut[375,4],mo_da_link_cut[375,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[376,3],mo_da_link_cut[376,5]),y=c(mo_da_link_cut[376,4],mo_da_link_cut[376,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[377,3],mo_da_link_cut[377,5]),y=c(mo_da_link_cut[377,4],mo_da_link_cut[377,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[378,3],mo_da_link_cut[378,5]),y=c(mo_da_link_cut[378,4],mo_da_link_cut[378,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[379,3],mo_da_link_cut[379,5]),y=c(mo_da_link_cut[379,4],mo_da_link_cut[379,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[380,3],mo_da_link_cut[380,5]),y=c(mo_da_link_cut[380,4],mo_da_link_cut[380,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[381,3],mo_da_link_cut[381,5]),y=c(mo_da_link_cut[381,4],mo_da_link_cut[381,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[382,3],mo_da_link_cut[382,5]),y=c(mo_da_link_cut[382,4],mo_da_link_cut[382,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[383,3],mo_da_link_cut[383,5]),y=c(mo_da_link_cut[383,4],mo_da_link_cut[383,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[384,3],mo_da_link_cut[384,5]),y=c(mo_da_link_cut[384,4],mo_da_link_cut[384,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[385,3],mo_da_link_cut[385,5]),y=c(mo_da_link_cut[385,4],mo_da_link_cut[385,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[386,3],mo_da_link_cut[386,5]),y=c(mo_da_link_cut[386,4],mo_da_link_cut[386,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[387,3],mo_da_link_cut[387,5]),y=c(mo_da_link_cut[387,4],mo_da_link_cut[387,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[388,3],mo_da_link_cut[388,5]),y=c(mo_da_link_cut[388,4],mo_da_link_cut[388,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[389,3],mo_da_link_cut[389,5]),y=c(mo_da_link_cut[389,4],mo_da_link_cut[389,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[390,3],mo_da_link_cut[390,5]),y=c(mo_da_link_cut[390,4],mo_da_link_cut[390,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[391,3],mo_da_link_cut[391,5]),y=c(mo_da_link_cut[391,4],mo_da_link_cut[391,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[392,3],mo_da_link_cut[392,5]),y=c(mo_da_link_cut[392,4],mo_da_link_cut[392,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[393,3],mo_da_link_cut[393,5]),y=c(mo_da_link_cut[393,4],mo_da_link_cut[393,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[394,3],mo_da_link_cut[394,5]),y=c(mo_da_link_cut[394,4],mo_da_link_cut[394,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[395,3],mo_da_link_cut[395,5]),y=c(mo_da_link_cut[395,4],mo_da_link_cut[395,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[396,3],mo_da_link_cut[396,5]),y=c(mo_da_link_cut[396,4],mo_da_link_cut[396,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[397,3],mo_da_link_cut[397,5]),y=c(mo_da_link_cut[397,4],mo_da_link_cut[397,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[398,3],mo_da_link_cut[398,5]),y=c(mo_da_link_cut[398,4],mo_da_link_cut[398,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[399,3],mo_da_link_cut[399,5]),y=c(mo_da_link_cut[399,4],mo_da_link_cut[399,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[400,3],mo_da_link_cut[400,5]),y=c(mo_da_link_cut[400,4],mo_da_link_cut[400,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[401,3],mo_da_link_cut[401,5]),y=c(mo_da_link_cut[401,4],mo_da_link_cut[401,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[402,3],mo_da_link_cut[402,5]),y=c(mo_da_link_cut[402,4],mo_da_link_cut[402,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[403,3],mo_da_link_cut[403,5]),y=c(mo_da_link_cut[403,4],mo_da_link_cut[403,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[404,3],mo_da_link_cut[404,5]),y=c(mo_da_link_cut[404,4],mo_da_link_cut[404,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[405,3],mo_da_link_cut[405,5]),y=c(mo_da_link_cut[405,4],mo_da_link_cut[405,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[406,3],mo_da_link_cut[406,5]),y=c(mo_da_link_cut[406,4],mo_da_link_cut[406,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[407,3],mo_da_link_cut[407,5]),y=c(mo_da_link_cut[407,4],mo_da_link_cut[407,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[408,3],mo_da_link_cut[408,5]),y=c(mo_da_link_cut[408,4],mo_da_link_cut[408,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[409,3],mo_da_link_cut[409,5]),y=c(mo_da_link_cut[409,4],mo_da_link_cut[409,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[410,3],mo_da_link_cut[410,5]),y=c(mo_da_link_cut[410,4],mo_da_link_cut[410,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[411,3],mo_da_link_cut[411,5]),y=c(mo_da_link_cut[411,4],mo_da_link_cut[411,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[412,3],mo_da_link_cut[412,5]),y=c(mo_da_link_cut[412,4],mo_da_link_cut[412,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[413,3],mo_da_link_cut[413,5]),y=c(mo_da_link_cut[413,4],mo_da_link_cut[413,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[414,3],mo_da_link_cut[414,5]),y=c(mo_da_link_cut[414,4],mo_da_link_cut[414,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[415,3],mo_da_link_cut[415,5]),y=c(mo_da_link_cut[415,4],mo_da_link_cut[415,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[416,3],mo_da_link_cut[416,5]),y=c(mo_da_link_cut[416,4],mo_da_link_cut[416,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[417,3],mo_da_link_cut[417,5]),y=c(mo_da_link_cut[417,4],mo_da_link_cut[417,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[418,3],mo_da_link_cut[418,5]),y=c(mo_da_link_cut[418,4],mo_da_link_cut[418,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[419,3],mo_da_link_cut[419,5]),y=c(mo_da_link_cut[419,4],mo_da_link_cut[419,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[420,3],mo_da_link_cut[420,5]),y=c(mo_da_link_cut[420,4],mo_da_link_cut[420,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[421,3],mo_da_link_cut[421,5]),y=c(mo_da_link_cut[421,4],mo_da_link_cut[421,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[422,3],mo_da_link_cut[422,5]),y=c(mo_da_link_cut[422,4],mo_da_link_cut[422,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[423,3],mo_da_link_cut[423,5]),y=c(mo_da_link_cut[423,4],mo_da_link_cut[423,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[424,3],mo_da_link_cut[424,5]),y=c(mo_da_link_cut[424,4],mo_da_link_cut[424,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[425,3],mo_da_link_cut[425,5]),y=c(mo_da_link_cut[425,4],mo_da_link_cut[425,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[426,3],mo_da_link_cut[426,5]),y=c(mo_da_link_cut[426,4],mo_da_link_cut[426,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[427,3],mo_da_link_cut[427,5]),y=c(mo_da_link_cut[427,4],mo_da_link_cut[427,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[428,3],mo_da_link_cut[428,5]),y=c(mo_da_link_cut[428,4],mo_da_link_cut[428,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[429,3],mo_da_link_cut[429,5]),y=c(mo_da_link_cut[429,4],mo_da_link_cut[429,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[430,3],mo_da_link_cut[430,5]),y=c(mo_da_link_cut[430,4],mo_da_link_cut[430,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[431,3],mo_da_link_cut[431,5]),y=c(mo_da_link_cut[431,4],mo_da_link_cut[431,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[432,3],mo_da_link_cut[432,5]),y=c(mo_da_link_cut[432,4],mo_da_link_cut[432,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[433,3],mo_da_link_cut[433,5]),y=c(mo_da_link_cut[433,4],mo_da_link_cut[433,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[434,3],mo_da_link_cut[434,5]),y=c(mo_da_link_cut[434,4],mo_da_link_cut[434,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[435,3],mo_da_link_cut[435,5]),y=c(mo_da_link_cut[435,4],mo_da_link_cut[435,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[436,3],mo_da_link_cut[436,5]),y=c(mo_da_link_cut[436,4],mo_da_link_cut[436,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[437,3],mo_da_link_cut[437,5]),y=c(mo_da_link_cut[437,4],mo_da_link_cut[437,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[438,3],mo_da_link_cut[438,5]),y=c(mo_da_link_cut[438,4],mo_da_link_cut[438,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[439,3],mo_da_link_cut[439,5]),y=c(mo_da_link_cut[439,4],mo_da_link_cut[439,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[440,3],mo_da_link_cut[440,5]),y=c(mo_da_link_cut[440,4],mo_da_link_cut[440,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[441,3],mo_da_link_cut[441,5]),y=c(mo_da_link_cut[441,4],mo_da_link_cut[441,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[442,3],mo_da_link_cut[442,5]),y=c(mo_da_link_cut[442,4],mo_da_link_cut[442,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[443,3],mo_da_link_cut[443,5]),y=c(mo_da_link_cut[443,4],mo_da_link_cut[443,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[444,3],mo_da_link_cut[444,5]),y=c(mo_da_link_cut[444,4],mo_da_link_cut[444,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[445,3],mo_da_link_cut[445,5]),y=c(mo_da_link_cut[445,4],mo_da_link_cut[445,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[446,3],mo_da_link_cut[446,5]),y=c(mo_da_link_cut[446,4],mo_da_link_cut[446,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[447,3],mo_da_link_cut[447,5]),y=c(mo_da_link_cut[447,4],mo_da_link_cut[447,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[448,3],mo_da_link_cut[448,5]),y=c(mo_da_link_cut[448,4],mo_da_link_cut[448,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[449,3],mo_da_link_cut[449,5]),y=c(mo_da_link_cut[449,4],mo_da_link_cut[449,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[450,3],mo_da_link_cut[450,5]),y=c(mo_da_link_cut[450,4],mo_da_link_cut[450,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[451,3],mo_da_link_cut[451,5]),y=c(mo_da_link_cut[451,4],mo_da_link_cut[451,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[452,3],mo_da_link_cut[452,5]),y=c(mo_da_link_cut[452,4],mo_da_link_cut[452,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[453,3],mo_da_link_cut[453,5]),y=c(mo_da_link_cut[453,4],mo_da_link_cut[453,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[454,3],mo_da_link_cut[454,5]),y=c(mo_da_link_cut[454,4],mo_da_link_cut[454,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[455,3],mo_da_link_cut[455,5]),y=c(mo_da_link_cut[455,4],mo_da_link_cut[455,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[456,3],mo_da_link_cut[456,5]),y=c(mo_da_link_cut[456,4],mo_da_link_cut[456,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[457,3],mo_da_link_cut[457,5]),y=c(mo_da_link_cut[457,4],mo_da_link_cut[457,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[458,3],mo_da_link_cut[458,5]),y=c(mo_da_link_cut[458,4],mo_da_link_cut[458,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[459,3],mo_da_link_cut[459,5]),y=c(mo_da_link_cut[459,4],mo_da_link_cut[459,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[460,3],mo_da_link_cut[460,5]),y=c(mo_da_link_cut[460,4],mo_da_link_cut[460,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[461,3],mo_da_link_cut[461,5]),y=c(mo_da_link_cut[461,4],mo_da_link_cut[461,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[462,3],mo_da_link_cut[462,5]),y=c(mo_da_link_cut[462,4],mo_da_link_cut[462,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[463,3],mo_da_link_cut[463,5]),y=c(mo_da_link_cut[463,4],mo_da_link_cut[463,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[464,3],mo_da_link_cut[464,5]),y=c(mo_da_link_cut[464,4],mo_da_link_cut[464,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[465,3],mo_da_link_cut[465,5]),y=c(mo_da_link_cut[465,4],mo_da_link_cut[465,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[466,3],mo_da_link_cut[466,5]),y=c(mo_da_link_cut[466,4],mo_da_link_cut[466,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[467,3],mo_da_link_cut[467,5]),y=c(mo_da_link_cut[467,4],mo_da_link_cut[467,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[468,3],mo_da_link_cut[468,5]),y=c(mo_da_link_cut[468,4],mo_da_link_cut[468,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[469,3],mo_da_link_cut[469,5]),y=c(mo_da_link_cut[469,4],mo_da_link_cut[469,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[470,3],mo_da_link_cut[470,5]),y=c(mo_da_link_cut[470,4],mo_da_link_cut[470,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[471,3],mo_da_link_cut[471,5]),y=c(mo_da_link_cut[471,4],mo_da_link_cut[471,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[472,3],mo_da_link_cut[472,5]),y=c(mo_da_link_cut[472,4],mo_da_link_cut[472,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[473,3],mo_da_link_cut[473,5]),y=c(mo_da_link_cut[473,4],mo_da_link_cut[473,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[474,3],mo_da_link_cut[474,5]),y=c(mo_da_link_cut[474,4],mo_da_link_cut[474,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[475,3],mo_da_link_cut[475,5]),y=c(mo_da_link_cut[475,4],mo_da_link_cut[475,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[476,3],mo_da_link_cut[476,5]),y=c(mo_da_link_cut[476,4],mo_da_link_cut[476,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[477,3],mo_da_link_cut[477,5]),y=c(mo_da_link_cut[477,4],mo_da_link_cut[477,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[478,3],mo_da_link_cut[478,5]),y=c(mo_da_link_cut[478,4],mo_da_link_cut[478,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[479,3],mo_da_link_cut[479,5]),y=c(mo_da_link_cut[479,4],mo_da_link_cut[479,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[480,3],mo_da_link_cut[480,5]),y=c(mo_da_link_cut[480,4],mo_da_link_cut[480,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[481,3],mo_da_link_cut[481,5]),y=c(mo_da_link_cut[481,4],mo_da_link_cut[481,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[482,3],mo_da_link_cut[482,5]),y=c(mo_da_link_cut[482,4],mo_da_link_cut[482,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[483,3],mo_da_link_cut[483,5]),y=c(mo_da_link_cut[483,4],mo_da_link_cut[483,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[484,3],mo_da_link_cut[484,5]),y=c(mo_da_link_cut[484,4],mo_da_link_cut[484,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[485,3],mo_da_link_cut[485,5]),y=c(mo_da_link_cut[485,4],mo_da_link_cut[485,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[486,3],mo_da_link_cut[486,5]),y=c(mo_da_link_cut[486,4],mo_da_link_cut[486,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[487,3],mo_da_link_cut[487,5]),y=c(mo_da_link_cut[487,4],mo_da_link_cut[487,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[488,3],mo_da_link_cut[488,5]),y=c(mo_da_link_cut[488,4],mo_da_link_cut[488,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[489,3],mo_da_link_cut[489,5]),y=c(mo_da_link_cut[489,4],mo_da_link_cut[489,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[490,3],mo_da_link_cut[490,5]),y=c(mo_da_link_cut[490,4],mo_da_link_cut[490,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[491,3],mo_da_link_cut[491,5]),y=c(mo_da_link_cut[491,4],mo_da_link_cut[491,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[492,3],mo_da_link_cut[492,5]),y=c(mo_da_link_cut[492,4],mo_da_link_cut[492,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[493,3],mo_da_link_cut[493,5]),y=c(mo_da_link_cut[493,4],mo_da_link_cut[493,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[494,3],mo_da_link_cut[494,5]),y=c(mo_da_link_cut[494,4],mo_da_link_cut[494,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[495,3],mo_da_link_cut[495,5]),y=c(mo_da_link_cut[495,4],mo_da_link_cut[495,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[496,3],mo_da_link_cut[496,5]),y=c(mo_da_link_cut[496,4],mo_da_link_cut[496,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[497,3],mo_da_link_cut[497,5]),y=c(mo_da_link_cut[497,4],mo_da_link_cut[497,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[498,3],mo_da_link_cut[498,5]),y=c(mo_da_link_cut[498,4],mo_da_link_cut[498,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[499,3],mo_da_link_cut[499,5]),y=c(mo_da_link_cut[499,4],mo_da_link_cut[499,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[500,3],mo_da_link_cut[500,5]),y=c(mo_da_link_cut[500,4],mo_da_link_cut[500,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[501,3],mo_da_link_cut[501,5]),y=c(mo_da_link_cut[501,4],mo_da_link_cut[501,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[502,3],mo_da_link_cut[502,5]),y=c(mo_da_link_cut[502,4],mo_da_link_cut[502,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[503,3],mo_da_link_cut[503,5]),y=c(mo_da_link_cut[503,4],mo_da_link_cut[503,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[504,3],mo_da_link_cut[504,5]),y=c(mo_da_link_cut[504,4],mo_da_link_cut[504,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[505,3],mo_da_link_cut[505,5]),y=c(mo_da_link_cut[505,4],mo_da_link_cut[505,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[506,3],mo_da_link_cut[506,5]),y=c(mo_da_link_cut[506,4],mo_da_link_cut[506,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[507,3],mo_da_link_cut[507,5]),y=c(mo_da_link_cut[507,4],mo_da_link_cut[507,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[508,3],mo_da_link_cut[508,5]),y=c(mo_da_link_cut[508,4],mo_da_link_cut[508,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[509,3],mo_da_link_cut[509,5]),y=c(mo_da_link_cut[509,4],mo_da_link_cut[509,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[510,3],mo_da_link_cut[510,5]),y=c(mo_da_link_cut[510,4],mo_da_link_cut[510,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[511,3],mo_da_link_cut[511,5]),y=c(mo_da_link_cut[511,4],mo_da_link_cut[511,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[512,3],mo_da_link_cut[512,5]),y=c(mo_da_link_cut[512,4],mo_da_link_cut[512,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[513,3],mo_da_link_cut[513,5]),y=c(mo_da_link_cut[513,4],mo_da_link_cut[513,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[514,3],mo_da_link_cut[514,5]),y=c(mo_da_link_cut[514,4],mo_da_link_cut[514,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[515,3],mo_da_link_cut[515,5]),y=c(mo_da_link_cut[515,4],mo_da_link_cut[515,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[516,3],mo_da_link_cut[516,5]),y=c(mo_da_link_cut[516,4],mo_da_link_cut[516,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[517,3],mo_da_link_cut[517,5]),y=c(mo_da_link_cut[517,4],mo_da_link_cut[517,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[518,3],mo_da_link_cut[518,5]),y=c(mo_da_link_cut[518,4],mo_da_link_cut[518,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[519,3],mo_da_link_cut[519,5]),y=c(mo_da_link_cut[519,4],mo_da_link_cut[519,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[520,3],mo_da_link_cut[520,5]),y=c(mo_da_link_cut[520,4],mo_da_link_cut[520,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[521,3],mo_da_link_cut[521,5]),y=c(mo_da_link_cut[521,4],mo_da_link_cut[521,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[522,3],mo_da_link_cut[522,5]),y=c(mo_da_link_cut[522,4],mo_da_link_cut[522,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[523,3],mo_da_link_cut[523,5]),y=c(mo_da_link_cut[523,4],mo_da_link_cut[523,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[524,3],mo_da_link_cut[524,5]),y=c(mo_da_link_cut[524,4],mo_da_link_cut[524,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[525,3],mo_da_link_cut[525,5]),y=c(mo_da_link_cut[525,4],mo_da_link_cut[525,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[526,3],mo_da_link_cut[526,5]),y=c(mo_da_link_cut[526,4],mo_da_link_cut[526,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[527,3],mo_da_link_cut[527,5]),y=c(mo_da_link_cut[527,4],mo_da_link_cut[527,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[528,3],mo_da_link_cut[528,5]),y=c(mo_da_link_cut[528,4],mo_da_link_cut[528,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[529,3],mo_da_link_cut[529,5]),y=c(mo_da_link_cut[529,4],mo_da_link_cut[529,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[530,3],mo_da_link_cut[530,5]),y=c(mo_da_link_cut[530,4],mo_da_link_cut[530,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[531,3],mo_da_link_cut[531,5]),y=c(mo_da_link_cut[531,4],mo_da_link_cut[531,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[532,3],mo_da_link_cut[532,5]),y=c(mo_da_link_cut[532,4],mo_da_link_cut[532,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[533,3],mo_da_link_cut[533,5]),y=c(mo_da_link_cut[533,4],mo_da_link_cut[533,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[534,3],mo_da_link_cut[534,5]),y=c(mo_da_link_cut[534,4],mo_da_link_cut[534,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[535,3],mo_da_link_cut[535,5]),y=c(mo_da_link_cut[535,4],mo_da_link_cut[535,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[536,3],mo_da_link_cut[536,5]),y=c(mo_da_link_cut[536,4],mo_da_link_cut[536,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[537,3],mo_da_link_cut[537,5]),y=c(mo_da_link_cut[537,4],mo_da_link_cut[537,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[538,3],mo_da_link_cut[538,5]),y=c(mo_da_link_cut[538,4],mo_da_link_cut[538,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[539,3],mo_da_link_cut[539,5]),y=c(mo_da_link_cut[539,4],mo_da_link_cut[539,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[540,3],mo_da_link_cut[540,5]),y=c(mo_da_link_cut[540,4],mo_da_link_cut[540,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[541,3],mo_da_link_cut[541,5]),y=c(mo_da_link_cut[541,4],mo_da_link_cut[541,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[542,3],mo_da_link_cut[542,5]),y=c(mo_da_link_cut[542,4],mo_da_link_cut[542,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[543,3],mo_da_link_cut[543,5]),y=c(mo_da_link_cut[543,4],mo_da_link_cut[543,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[544,3],mo_da_link_cut[544,5]),y=c(mo_da_link_cut[544,4],mo_da_link_cut[544,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[545,3],mo_da_link_cut[545,5]),y=c(mo_da_link_cut[545,4],mo_da_link_cut[545,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[546,3],mo_da_link_cut[546,5]),y=c(mo_da_link_cut[546,4],mo_da_link_cut[546,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[547,3],mo_da_link_cut[547,5]),y=c(mo_da_link_cut[547,4],mo_da_link_cut[547,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[548,3],mo_da_link_cut[548,5]),y=c(mo_da_link_cut[548,4],mo_da_link_cut[548,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[549,3],mo_da_link_cut[549,5]),y=c(mo_da_link_cut[549,4],mo_da_link_cut[549,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[550,3],mo_da_link_cut[550,5]),y=c(mo_da_link_cut[550,4],mo_da_link_cut[550,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[551,3],mo_da_link_cut[551,5]),y=c(mo_da_link_cut[551,4],mo_da_link_cut[551,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[552,3],mo_da_link_cut[552,5]),y=c(mo_da_link_cut[552,4],mo_da_link_cut[552,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[553,3],mo_da_link_cut[553,5]),y=c(mo_da_link_cut[553,4],mo_da_link_cut[553,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[554,3],mo_da_link_cut[554,5]),y=c(mo_da_link_cut[554,4],mo_da_link_cut[554,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[555,3],mo_da_link_cut[555,5]),y=c(mo_da_link_cut[555,4],mo_da_link_cut[555,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[556,3],mo_da_link_cut[556,5]),y=c(mo_da_link_cut[556,4],mo_da_link_cut[556,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[557,3],mo_da_link_cut[557,5]),y=c(mo_da_link_cut[557,4],mo_da_link_cut[557,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[558,3],mo_da_link_cut[558,5]),y=c(mo_da_link_cut[558,4],mo_da_link_cut[558,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[559,3],mo_da_link_cut[559,5]),y=c(mo_da_link_cut[559,4],mo_da_link_cut[559,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[560,3],mo_da_link_cut[560,5]),y=c(mo_da_link_cut[560,4],mo_da_link_cut[560,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[561,3],mo_da_link_cut[561,5]),y=c(mo_da_link_cut[561,4],mo_da_link_cut[561,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[562,3],mo_da_link_cut[562,5]),y=c(mo_da_link_cut[562,4],mo_da_link_cut[562,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[563,3],mo_da_link_cut[563,5]),y=c(mo_da_link_cut[563,4],mo_da_link_cut[563,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[564,3],mo_da_link_cut[564,5]),y=c(mo_da_link_cut[564,4],mo_da_link_cut[564,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[565,3],mo_da_link_cut[565,5]),y=c(mo_da_link_cut[565,4],mo_da_link_cut[565,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[566,3],mo_da_link_cut[566,5]),y=c(mo_da_link_cut[566,4],mo_da_link_cut[566,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[567,3],mo_da_link_cut[567,5]),y=c(mo_da_link_cut[567,4],mo_da_link_cut[567,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[568,3],mo_da_link_cut[568,5]),y=c(mo_da_link_cut[568,4],mo_da_link_cut[568,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[569,3],mo_da_link_cut[569,5]),y=c(mo_da_link_cut[569,4],mo_da_link_cut[569,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[570,3],mo_da_link_cut[570,5]),y=c(mo_da_link_cut[570,4],mo_da_link_cut[570,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[571,3],mo_da_link_cut[571,5]),y=c(mo_da_link_cut[571,4],mo_da_link_cut[571,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[572,3],mo_da_link_cut[572,5]),y=c(mo_da_link_cut[572,4],mo_da_link_cut[572,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[573,3],mo_da_link_cut[573,5]),y=c(mo_da_link_cut[573,4],mo_da_link_cut[573,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[574,3],mo_da_link_cut[574,5]),y=c(mo_da_link_cut[574,4],mo_da_link_cut[574,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[575,3],mo_da_link_cut[575,5]),y=c(mo_da_link_cut[575,4],mo_da_link_cut[575,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[576,3],mo_da_link_cut[576,5]),y=c(mo_da_link_cut[576,4],mo_da_link_cut[576,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[577,3],mo_da_link_cut[577,5]),y=c(mo_da_link_cut[577,4],mo_da_link_cut[577,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[578,3],mo_da_link_cut[578,5]),y=c(mo_da_link_cut[578,4],mo_da_link_cut[578,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[579,3],mo_da_link_cut[579,5]),y=c(mo_da_link_cut[579,4],mo_da_link_cut[579,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[580,3],mo_da_link_cut[580,5]),y=c(mo_da_link_cut[580,4],mo_da_link_cut[580,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[581,3],mo_da_link_cut[581,5]),y=c(mo_da_link_cut[581,4],mo_da_link_cut[581,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[582,3],mo_da_link_cut[582,5]),y=c(mo_da_link_cut[582,4],mo_da_link_cut[582,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[583,3],mo_da_link_cut[583,5]),y=c(mo_da_link_cut[583,4],mo_da_link_cut[583,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[584,3],mo_da_link_cut[584,5]),y=c(mo_da_link_cut[584,4],mo_da_link_cut[584,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[585,3],mo_da_link_cut[585,5]),y=c(mo_da_link_cut[585,4],mo_da_link_cut[585,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[586,3],mo_da_link_cut[586,5]),y=c(mo_da_link_cut[586,4],mo_da_link_cut[586,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[587,3],mo_da_link_cut[587,5]),y=c(mo_da_link_cut[587,4],mo_da_link_cut[587,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[588,3],mo_da_link_cut[588,5]),y=c(mo_da_link_cut[588,4],mo_da_link_cut[588,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[589,3],mo_da_link_cut[589,5]),y=c(mo_da_link_cut[589,4],mo_da_link_cut[589,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[590,3],mo_da_link_cut[590,5]),y=c(mo_da_link_cut[590,4],mo_da_link_cut[590,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[591,3],mo_da_link_cut[591,5]),y=c(mo_da_link_cut[591,4],mo_da_link_cut[591,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[592,3],mo_da_link_cut[592,5]),y=c(mo_da_link_cut[592,4],mo_da_link_cut[592,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[593,3],mo_da_link_cut[593,5]),y=c(mo_da_link_cut[593,4],mo_da_link_cut[593,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[594,3],mo_da_link_cut[594,5]),y=c(mo_da_link_cut[594,4],mo_da_link_cut[594,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[595,3],mo_da_link_cut[595,5]),y=c(mo_da_link_cut[595,4],mo_da_link_cut[595,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[596,3],mo_da_link_cut[596,5]),y=c(mo_da_link_cut[596,4],mo_da_link_cut[596,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[597,3],mo_da_link_cut[597,5]),y=c(mo_da_link_cut[597,4],mo_da_link_cut[597,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[598,3],mo_da_link_cut[598,5]),y=c(mo_da_link_cut[598,4],mo_da_link_cut[598,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[599,3],mo_da_link_cut[599,5]),y=c(mo_da_link_cut[599,4],mo_da_link_cut[599,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[600,3],mo_da_link_cut[600,5]),y=c(mo_da_link_cut[600,4],mo_da_link_cut[600,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[601,3],mo_da_link_cut[601,5]),y=c(mo_da_link_cut[601,4],mo_da_link_cut[601,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[602,3],mo_da_link_cut[602,5]),y=c(mo_da_link_cut[602,4],mo_da_link_cut[602,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[603,3],mo_da_link_cut[603,5]),y=c(mo_da_link_cut[603,4],mo_da_link_cut[603,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[604,3],mo_da_link_cut[604,5]),y=c(mo_da_link_cut[604,4],mo_da_link_cut[604,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[605,3],mo_da_link_cut[605,5]),y=c(mo_da_link_cut[605,4],mo_da_link_cut[605,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[606,3],mo_da_link_cut[606,5]),y=c(mo_da_link_cut[606,4],mo_da_link_cut[606,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[607,3],mo_da_link_cut[607,5]),y=c(mo_da_link_cut[607,4],mo_da_link_cut[607,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[608,3],mo_da_link_cut[608,5]),y=c(mo_da_link_cut[608,4],mo_da_link_cut[608,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[609,3],mo_da_link_cut[609,5]),y=c(mo_da_link_cut[609,4],mo_da_link_cut[609,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[610,3],mo_da_link_cut[610,5]),y=c(mo_da_link_cut[610,4],mo_da_link_cut[610,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[611,3],mo_da_link_cut[611,5]),y=c(mo_da_link_cut[611,4],mo_da_link_cut[611,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[612,3],mo_da_link_cut[612,5]),y=c(mo_da_link_cut[612,4],mo_da_link_cut[612,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[613,3],mo_da_link_cut[613,5]),y=c(mo_da_link_cut[613,4],mo_da_link_cut[613,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[614,3],mo_da_link_cut[614,5]),y=c(mo_da_link_cut[614,4],mo_da_link_cut[614,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[615,3],mo_da_link_cut[615,5]),y=c(mo_da_link_cut[615,4],mo_da_link_cut[615,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[616,3],mo_da_link_cut[616,5]),y=c(mo_da_link_cut[616,4],mo_da_link_cut[616,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[617,3],mo_da_link_cut[617,5]),y=c(mo_da_link_cut[617,4],mo_da_link_cut[617,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[618,3],mo_da_link_cut[618,5]),y=c(mo_da_link_cut[618,4],mo_da_link_cut[618,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[619,3],mo_da_link_cut[619,5]),y=c(mo_da_link_cut[619,4],mo_da_link_cut[619,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[620,3],mo_da_link_cut[620,5]),y=c(mo_da_link_cut[620,4],mo_da_link_cut[620,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[621,3],mo_da_link_cut[621,5]),y=c(mo_da_link_cut[621,4],mo_da_link_cut[621,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[622,3],mo_da_link_cut[622,5]),y=c(mo_da_link_cut[622,4],mo_da_link_cut[622,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[623,3],mo_da_link_cut[623,5]),y=c(mo_da_link_cut[623,4],mo_da_link_cut[623,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[624,3],mo_da_link_cut[624,5]),y=c(mo_da_link_cut[624,4],mo_da_link_cut[624,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[625,3],mo_da_link_cut[625,5]),y=c(mo_da_link_cut[625,4],mo_da_link_cut[625,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[626,3],mo_da_link_cut[626,5]),y=c(mo_da_link_cut[626,4],mo_da_link_cut[626,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[627,3],mo_da_link_cut[627,5]),y=c(mo_da_link_cut[627,4],mo_da_link_cut[627,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[628,3],mo_da_link_cut[628,5]),y=c(mo_da_link_cut[628,4],mo_da_link_cut[628,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[629,3],mo_da_link_cut[629,5]),y=c(mo_da_link_cut[629,4],mo_da_link_cut[629,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[630,3],mo_da_link_cut[630,5]),y=c(mo_da_link_cut[630,4],mo_da_link_cut[630,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[631,3],mo_da_link_cut[631,5]),y=c(mo_da_link_cut[631,4],mo_da_link_cut[631,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[632,3],mo_da_link_cut[632,5]),y=c(mo_da_link_cut[632,4],mo_da_link_cut[632,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[633,3],mo_da_link_cut[633,5]),y=c(mo_da_link_cut[633,4],mo_da_link_cut[633,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[634,3],mo_da_link_cut[634,5]),y=c(mo_da_link_cut[634,4],mo_da_link_cut[634,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[635,3],mo_da_link_cut[635,5]),y=c(mo_da_link_cut[635,4],mo_da_link_cut[635,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[636,3],mo_da_link_cut[636,5]),y=c(mo_da_link_cut[636,4],mo_da_link_cut[636,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[637,3],mo_da_link_cut[637,5]),y=c(mo_da_link_cut[637,4],mo_da_link_cut[637,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[638,3],mo_da_link_cut[638,5]),y=c(mo_da_link_cut[638,4],mo_da_link_cut[638,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[639,3],mo_da_link_cut[639,5]),y=c(mo_da_link_cut[639,4],mo_da_link_cut[639,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[640,3],mo_da_link_cut[640,5]),y=c(mo_da_link_cut[640,4],mo_da_link_cut[640,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[641,3],mo_da_link_cut[641,5]),y=c(mo_da_link_cut[641,4],mo_da_link_cut[641,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[642,3],mo_da_link_cut[642,5]),y=c(mo_da_link_cut[642,4],mo_da_link_cut[642,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[643,3],mo_da_link_cut[643,5]),y=c(mo_da_link_cut[643,4],mo_da_link_cut[643,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[644,3],mo_da_link_cut[644,5]),y=c(mo_da_link_cut[644,4],mo_da_link_cut[644,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[645,3],mo_da_link_cut[645,5]),y=c(mo_da_link_cut[645,4],mo_da_link_cut[645,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[646,3],mo_da_link_cut[646,5]),y=c(mo_da_link_cut[646,4],mo_da_link_cut[646,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[647,3],mo_da_link_cut[647,5]),y=c(mo_da_link_cut[647,4],mo_da_link_cut[647,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[648,3],mo_da_link_cut[648,5]),y=c(mo_da_link_cut[648,4],mo_da_link_cut[648,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[649,3],mo_da_link_cut[649,5]),y=c(mo_da_link_cut[649,4],mo_da_link_cut[649,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[650,3],mo_da_link_cut[650,5]),y=c(mo_da_link_cut[650,4],mo_da_link_cut[650,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[651,3],mo_da_link_cut[651,5]),y=c(mo_da_link_cut[651,4],mo_da_link_cut[651,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[652,3],mo_da_link_cut[652,5]),y=c(mo_da_link_cut[652,4],mo_da_link_cut[652,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[653,3],mo_da_link_cut[653,5]),y=c(mo_da_link_cut[653,4],mo_da_link_cut[653,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[654,3],mo_da_link_cut[654,5]),y=c(mo_da_link_cut[654,4],mo_da_link_cut[654,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[655,3],mo_da_link_cut[655,5]),y=c(mo_da_link_cut[655,4],mo_da_link_cut[655,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[656,3],mo_da_link_cut[656,5]),y=c(mo_da_link_cut[656,4],mo_da_link_cut[656,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[657,3],mo_da_link_cut[657,5]),y=c(mo_da_link_cut[657,4],mo_da_link_cut[657,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[658,3],mo_da_link_cut[658,5]),y=c(mo_da_link_cut[658,4],mo_da_link_cut[658,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[659,3],mo_da_link_cut[659,5]),y=c(mo_da_link_cut[659,4],mo_da_link_cut[659,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[660,3],mo_da_link_cut[660,5]),y=c(mo_da_link_cut[660,4],mo_da_link_cut[660,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[661,3],mo_da_link_cut[661,5]),y=c(mo_da_link_cut[661,4],mo_da_link_cut[661,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[662,3],mo_da_link_cut[662,5]),y=c(mo_da_link_cut[662,4],mo_da_link_cut[662,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[663,3],mo_da_link_cut[663,5]),y=c(mo_da_link_cut[663,4],mo_da_link_cut[663,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[664,3],mo_da_link_cut[664,5]),y=c(mo_da_link_cut[664,4],mo_da_link_cut[664,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[665,3],mo_da_link_cut[665,5]),y=c(mo_da_link_cut[665,4],mo_da_link_cut[665,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[666,3],mo_da_link_cut[666,5]),y=c(mo_da_link_cut[666,4],mo_da_link_cut[666,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[667,3],mo_da_link_cut[667,5]),y=c(mo_da_link_cut[667,4],mo_da_link_cut[667,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[668,3],mo_da_link_cut[668,5]),y=c(mo_da_link_cut[668,4],mo_da_link_cut[668,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[669,3],mo_da_link_cut[669,5]),y=c(mo_da_link_cut[669,4],mo_da_link_cut[669,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[670,3],mo_da_link_cut[670,5]),y=c(mo_da_link_cut[670,4],mo_da_link_cut[670,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[671,3],mo_da_link_cut[671,5]),y=c(mo_da_link_cut[671,4],mo_da_link_cut[671,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[672,3],mo_da_link_cut[672,5]),y=c(mo_da_link_cut[672,4],mo_da_link_cut[672,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[673,3],mo_da_link_cut[673,5]),y=c(mo_da_link_cut[673,4],mo_da_link_cut[673,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[674,3],mo_da_link_cut[674,5]),y=c(mo_da_link_cut[674,4],mo_da_link_cut[674,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[675,3],mo_da_link_cut[675,5]),y=c(mo_da_link_cut[675,4],mo_da_link_cut[675,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[676,3],mo_da_link_cut[676,5]),y=c(mo_da_link_cut[676,4],mo_da_link_cut[676,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[677,3],mo_da_link_cut[677,5]),y=c(mo_da_link_cut[677,4],mo_da_link_cut[677,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[678,3],mo_da_link_cut[678,5]),y=c(mo_da_link_cut[678,4],mo_da_link_cut[678,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[679,3],mo_da_link_cut[679,5]),y=c(mo_da_link_cut[679,4],mo_da_link_cut[679,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[680,3],mo_da_link_cut[680,5]),y=c(mo_da_link_cut[680,4],mo_da_link_cut[680,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[681,3],mo_da_link_cut[681,5]),y=c(mo_da_link_cut[681,4],mo_da_link_cut[681,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[682,3],mo_da_link_cut[682,5]),y=c(mo_da_link_cut[682,4],mo_da_link_cut[682,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[683,3],mo_da_link_cut[683,5]),y=c(mo_da_link_cut[683,4],mo_da_link_cut[683,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[684,3],mo_da_link_cut[684,5]),y=c(mo_da_link_cut[684,4],mo_da_link_cut[684,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[685,3],mo_da_link_cut[685,5]),y=c(mo_da_link_cut[685,4],mo_da_link_cut[685,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[686,3],mo_da_link_cut[686,5]),y=c(mo_da_link_cut[686,4],mo_da_link_cut[686,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[687,3],mo_da_link_cut[687,5]),y=c(mo_da_link_cut[687,4],mo_da_link_cut[687,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[688,3],mo_da_link_cut[688,5]),y=c(mo_da_link_cut[688,4],mo_da_link_cut[688,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[689,3],mo_da_link_cut[689,5]),y=c(mo_da_link_cut[689,4],mo_da_link_cut[689,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[690,3],mo_da_link_cut[690,5]),y=c(mo_da_link_cut[690,4],mo_da_link_cut[690,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[691,3],mo_da_link_cut[691,5]),y=c(mo_da_link_cut[691,4],mo_da_link_cut[691,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[692,3],mo_da_link_cut[692,5]),y=c(mo_da_link_cut[692,4],mo_da_link_cut[692,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[693,3],mo_da_link_cut[693,5]),y=c(mo_da_link_cut[693,4],mo_da_link_cut[693,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[694,3],mo_da_link_cut[694,5]),y=c(mo_da_link_cut[694,4],mo_da_link_cut[694,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[695,3],mo_da_link_cut[695,5]),y=c(mo_da_link_cut[695,4],mo_da_link_cut[695,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[696,3],mo_da_link_cut[696,5]),y=c(mo_da_link_cut[696,4],mo_da_link_cut[696,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[697,3],mo_da_link_cut[697,5]),y=c(mo_da_link_cut[697,4],mo_da_link_cut[697,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[698,3],mo_da_link_cut[698,5]),y=c(mo_da_link_cut[698,4],mo_da_link_cut[698,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[699,3],mo_da_link_cut[699,5]),y=c(mo_da_link_cut[699,4],mo_da_link_cut[699,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[700,3],mo_da_link_cut[700,5]),y=c(mo_da_link_cut[700,4],mo_da_link_cut[700,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[701,3],mo_da_link_cut[701,5]),y=c(mo_da_link_cut[701,4],mo_da_link_cut[701,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[702,3],mo_da_link_cut[702,5]),y=c(mo_da_link_cut[702,4],mo_da_link_cut[702,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[703,3],mo_da_link_cut[703,5]),y=c(mo_da_link_cut[703,4],mo_da_link_cut[703,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[704,3],mo_da_link_cut[704,5]),y=c(mo_da_link_cut[704,4],mo_da_link_cut[704,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[705,3],mo_da_link_cut[705,5]),y=c(mo_da_link_cut[705,4],mo_da_link_cut[705,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[706,3],mo_da_link_cut[706,5]),y=c(mo_da_link_cut[706,4],mo_da_link_cut[706,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[707,3],mo_da_link_cut[707,5]),y=c(mo_da_link_cut[707,4],mo_da_link_cut[707,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[708,3],mo_da_link_cut[708,5]),y=c(mo_da_link_cut[708,4],mo_da_link_cut[708,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[709,3],mo_da_link_cut[709,5]),y=c(mo_da_link_cut[709,4],mo_da_link_cut[709,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[710,3],mo_da_link_cut[710,5]),y=c(mo_da_link_cut[710,4],mo_da_link_cut[710,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[711,3],mo_da_link_cut[711,5]),y=c(mo_da_link_cut[711,4],mo_da_link_cut[711,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[712,3],mo_da_link_cut[712,5]),y=c(mo_da_link_cut[712,4],mo_da_link_cut[712,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[713,3],mo_da_link_cut[713,5]),y=c(mo_da_link_cut[713,4],mo_da_link_cut[713,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[714,3],mo_da_link_cut[714,5]),y=c(mo_da_link_cut[714,4],mo_da_link_cut[714,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[715,3],mo_da_link_cut[715,5]),y=c(mo_da_link_cut[715,4],mo_da_link_cut[715,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[716,3],mo_da_link_cut[716,5]),y=c(mo_da_link_cut[716,4],mo_da_link_cut[716,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[717,3],mo_da_link_cut[717,5]),y=c(mo_da_link_cut[717,4],mo_da_link_cut[717,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[718,3],mo_da_link_cut[718,5]),y=c(mo_da_link_cut[718,4],mo_da_link_cut[718,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[719,3],mo_da_link_cut[719,5]),y=c(mo_da_link_cut[719,4],mo_da_link_cut[719,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_line(data=data.frame(x=c(mo_da_link_cut[720,3],mo_da_link_cut[720,5]),y=c(mo_da_link_cut[720,4],mo_da_link_cut[720,6])),aes(x=x,y=y),col="grey",lwd=0.1)+
geom_point(size=sqrt(triangle_tree_cut$var_exp)*0.01,stroke=0.1,pch=1,col="#00BFC4")+
geom_point(size=sqrt(triangle_tree_cut$var_obs)*0.01,stroke=0.1,pch=16,col="#F8766D")+
labs(title="")


ggsave("F:/wangj/Desktop/temp/canalization/figs/triangle_tree_noise_accumulate_suppress_total.pdf",width=unit(6.4,"in"),height=unit(1.9,"in"))



############################################## 
               without negative feedback
##############################################

setwd("F:/wangj/Desktop/temp/canalization/figs")

source("F:/wangj/Desktop/temp/canalization/functions/set_theme.R")
source("F:/wangj/Desktop/temp/canalization/functions/rowColRecovery.R")

load("F:/wangj/Desktop/temp/canalization/data/basic/xyz_aligned_sizeControl_t_center.RData")

load("F:/wangj/Desktop/temp/canalization/data/spatialTimePoints/q1_cell350_alignment_xyz_noise.RData")
load("F:/wangj/Desktop/temp/canalization/data/basic/cell_mo_da1_da2.RData")
load("F:/wangj/Desktop/temp/canalization/data/basic/lineage_info.RData")
load("F:/wangj/Desktop/temp/canalization/data/spatialTimePoints/coo_list_q6_centroid_rot_sizeControl.RData")


q6_xyz_aligned_sizeControl_t<-t(sapply(1:length(coo_list_q6_centroid_rot_sizeControl),function(i){c(coo_list_q6_centroid_rot_sizeControl[[i]][,1],coo_list_q6_centroid_rot_sizeControl[[i]][,2],coo_list_q6_centroid_rot_sizeControl[[i]][,3])}))

library(foreach)
library(glmnet)
library(parallel)
library(doParallel)


q6_x<-q1_cell350_alignment_xyz_noise[,1:714]
q6_y<-q1_cell350_alignment_xyz_noise[,c(1:714)+714]
q6_z<-q1_cell350_alignment_xyz_noise[,c(1:714)+714*2]


#

cell_name<-unlist(strsplit(colnames(q1_cell350_alignment_xyz_noise[,1:714]),"_x"))
cell_gen<-sapply(match(cell_name,lineage_info),function(x){rowColRec(x,dim(lineage_info)[1],dim(lineage_info)[2])[1]})

#


### x 

q6_x_without_nf<-q6_x

for(i in 5:10){
	focal_da_s<-cell_name[cell_gen==i]
	for(j in 1:length(focal_da_s)){
		
		focal_da<-focal_da_s[j]
		focal_mo<-cell_mo_da1_da2[which(cell_mo_da1_da2[,2]==focal_da),1]
		
		focal_mo_ind<-match(focal_mo,cell_name)
		focal_da_ind<-match(focal_da,cell_name)
		
		set.seed(1)
		q6_x_without_nf[,focal_da_ind]<-q6_x_without_nf[,focal_mo_ind]+residuals(lm(q6_x[,focal_da_ind]~q6_x[,focal_mo_ind]))
		
	}
	
}

q6_x_without_nf_pos<-q6_x_without_nf

for(i in 1:714){
	q6_x_without_nf_pos[,i]<-q6_x_without_nf[,i]+mean(q6_xyz_aligned_sizeControl_t[,i])
}


### y 

q6_y_without_nf<-q6_y

for(i in 5:10){
	focal_da_s<-cell_name[cell_gen==i]
	for(j in 1:length(focal_da_s)){
		
		focal_da<-focal_da_s[j]
		focal_mo<-cell_mo_da1_da2[which(cell_mo_da1_da2[,2]==focal_da),1]
		
		focal_mo_ind<-match(focal_mo,cell_name)
		focal_da_ind<-match(focal_da,cell_name)
		
		set.seed(1)
		q6_y_without_nf[,focal_da_ind]<-q6_y_without_nf[,focal_mo_ind]+residuals(lm(q6_y[,focal_da_ind]~q6_y[,focal_mo_ind]))
		
	}
	
}

q6_y_without_nf_pos<-q6_y_without_nf

for(i in 1:714){
	q6_y_without_nf_pos[,i]<-q6_y_without_nf[,i]+mean(q6_xyz_aligned_sizeControl_t[,i+714])
}


### z  

q6_z_without_nf<-q6_z

for(i in 5:10){
	focal_da_s<-cell_name[cell_gen==i]
	for(j in 1:length(focal_da_s)){
		
		focal_da<-focal_da_s[j]
		focal_mo<-cell_mo_da1_da2[which(cell_mo_da1_da2[,2]==focal_da),1]
		
		focal_mo_ind<-match(focal_mo,cell_name)
		focal_da_ind<-match(focal_da,cell_name)
		
		set.seed(1)
		q6_z_without_nf[,focal_da_ind]<-q6_z_without_nf[,focal_mo_ind]+residuals(lm(q6_z[,focal_da_ind]~q6_z[,focal_mo_ind]))
		
	}
	
}


q6_z_without_nf_pos<-q6_z_without_nf

for(i in 1:714){
	q6_z_without_nf_pos[,i]<-q6_z_without_nf[,i]+mean(q6_xyz_aligned_sizeControl_t[,i+714*2])
}


###

load("F:/wangj/desktop/temp/canalization/data/basic/cell_stage_df.RData")

### without nf 

cell_stage_df[,1]<-q6_x_without_nf_pos[1,]
cell_stage_df[,2]<-q6_y_without_nf_pos[1,]
cell_stage_df[,3]<-q6_z_without_nf_pos[1,]

dataTemp<-cell_stage_df

####

i=1
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("without_nf_lineage_xy_stage_1-1.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=2
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("without_nf_lineage_xy_stage_1-2.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=3
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("without_nf_lineage_xy_stage_1-3.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=4
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("without_nf_lineage_xy_stage_1-4.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=5
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("without_nf_lineage_xy_stage_1-5.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=6
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("without_nf_lineage_xy_stage_1-6.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=7
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("without_nf_lineage_xy_stage_1-7.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))



### with nf 

cell_stage_df[,1]<-q6_xyz_aligned_sizeControl_t[1,1:714]
cell_stage_df[,2]<-q6_xyz_aligned_sizeControl_t[1,c(1:714)+714]
cell_stage_df[,3]<-q6_xyz_aligned_sizeControl_t[1,c(1:714)+714*2]

dataTemp<-cell_stage_df

####

i=1
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("with_nf_lineage_xy_stage_1-1.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=2
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("with_nf_lineage_xy_stage_1-2.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=3
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("with_nf_lineage_xy_stage_1-3.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=4
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("with_nf_lineage_xy_stage_1-4.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=5
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("with_nf_lineage_xy_stage_1-5.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=6
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("with_nf_lineage_xy_stage_1-6.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


i=7
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("with_nf_lineage_xy_stage_1-7.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))


### MT1

### without nf 
j=106

cell_stage_df[,1]<-q6_x_without_nf_pos[j,]
cell_stage_df[,2]<-q6_y_without_nf_pos[j,]
cell_stage_df[,3]<-q6_z_without_nf_pos[j,]

dataTemp<-cell_stage_df

i=7
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("without_nf_lineage_xy_stage_2-7_ABHD-11.2_emb1.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))

### with nf 

j=106

cell_stage_df[,1]<-q6_xyz_aligned_sizeControl_t[j,1:714]
cell_stage_df[,2]<-q6_xyz_aligned_sizeControl_t[j,c(1:714)+714]
cell_stage_df[,3]<-q6_xyz_aligned_sizeControl_t[j,c(1:714)+714*2]

dataTemp<-cell_stage_df

i=7
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("with_nf_lineage_xy_stage_2-7_ABHD-11.2_emb1.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))



### MT2

### without nf 
j=121

cell_stage_df[,1]<-q6_x_without_nf_pos[j,]
cell_stage_df[,2]<-q6_y_without_nf_pos[j,]
cell_stage_df[,3]<-q6_z_without_nf_pos[j,]

dataTemp<-cell_stage_df

i=7
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("without_nf_lineage_xy_stage_3-7_ACL-11_emb1.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))

### with nf 

j=121

cell_stage_df[,1]<-q6_xyz_aligned_sizeControl_t[j,1:714]
cell_stage_df[,2]<-q6_xyz_aligned_sizeControl_t[j,c(1:714)+714]
cell_stage_df[,3]<-q6_xyz_aligned_sizeControl_t[j,c(1:714)+714*2]

dataTemp<-cell_stage_df

i=7
ggplot(dataTemp[dataTemp$stage==i,],aes(x=x,y=y,col=lineage))+
geom_point(show.legend=FALSE,cex=0.5,stroke=0.1)+
labs(x="",y="",title="")+
scale_x_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))+
scale_y_continuous(limits=c(min(dataTemp[,1:2]),max(dataTemp[,1:2])))

ggsave("with_nf_lineage_xy_stage_3-7_ACL-11_emb1.pdf",width=unit(1.5,"in"),height=unit(1.6,"in"))







