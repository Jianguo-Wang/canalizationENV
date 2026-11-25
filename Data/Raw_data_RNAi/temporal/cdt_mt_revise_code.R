



############################################################################
#					extract raw data
############################################################################

focal_projectFolderPath<-r"(E:\canalizationENV)"

##

setwd(paste0(focal_projectFolderPath,"/data/basic"))

load(paste0(focal_projectFolderPath,"/data/basic/q1_cell350_alignment_xyz.RData"))










#cdt_mt (artificial modification two RNAi embryos)

cdt_mt<-read.csv(paste0(focal_projectFolderPath,"/raw_data/Raw_data_RNAi/temporal/cell_division_time_mt.csv"),header=F)
cdt_mt_colnames<-as.character(cdt_mt[1,-1])
cdt_mt_rownames<-as.character(cdt_mt[-1,1])
cdt_mt<-read.csv(paste0(focal_projectFolderPath,"/raw_data/Raw_data_RNAi/temporal/cell_division_time_mt.csv"))
cdt_mt<-as.matrix(cdt_mt[,-1])
colnames(cdt_mt)<-cdt_mt_colnames
rownames(cdt_mt)<-cdt_mt_rownames

temp<-read.table(paste0(focal_projectFolderPath,"/Cell_division_timing_matrix_E03H4.8_and_NED-8.txt"),header=T,sep="\t")
temp<-temp[match(rownames(cdt_mt),temp[,1]),]
#write.csv(temp,paste0(focal_projectFolderPath,"/raw_data/Raw_data_RNAi/temporal/Cell_division_timing_matrix_E03H4.8_and_NED-8_refine.csv"),row.names=F)

cdt_mt[,match(c("E03H4.8_emb1","NED-8_emb1"),colnames(cdt_mt))]<-as.matrix(temp[,2:3])


write.csv(cdt_mt,paste0(focal_projectFolderPath,"/raw_data/Raw_data_RNAi/temporal/cell_division_time_mt_revise.csv"),row.names=TRUE)
