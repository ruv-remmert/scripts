heatmap<-read.delim("TES_nuc_regularity/values_Heatmap.txt", skip=1, header = F)
geneList<-heatmap[,1:6]
colnames(geneList)<-c("chr", "start", "end", "name", "nucPositioning","strand")
heatmapVal<-heatmap[,-(1:6)]

heatmapHEAD<-read.delim("TES_nuc_regularity/values_Heatmap.txt", nrows=1, header = F)

features<-heatmapHEAD %>% str_split("],")
#get positions
binSize<-grep("bin size",features[[1]], value = T) %>% str_split(",") %>% 
    unlist() %>%  first() %>% str_split_i("\\:\\[",2) %>% as.numeric() 
upstream<-grep("upstream",features[[1]], value = T) %>% str_split(",") %>% 
    unlist() %>%  first() %>% str_split_i("\\:\\[",2) %>% as.numeric() 
downstream<-grep("downstream",features[[1]], value = T) %>% str_split(",") %>% 
    unlist() %>%  first() %>% str_split_i("\\:\\[",2) %>% as.numeric() 

ncol_sample<-(upstream+downstream)/binSize
pos<-seq(c(-1)*(upstream-binSize/2),downstream-binSize/2,binSize)

samples<-grep("sample_labels",features[[1]], value = T) %>% 
     str_split_i("\\:\\[",2) %>% str_split(",") %>% unlist() %>% 
    str_split_i("Window_",2)

getHeatSamp<-function(name,ncol_sample=ncol_sample) {
    idx<-which(samples==name)
    idxHeat<-seq(((idx-1)*ncol_sample+1),idx*ncol_sample)
    heat<-heatmapVal[,idxHeat]
    colnames(heat)<-pos
    return(heat)
}


heat.tko.rap<-getHeatSamp(name="_TKO_Rpb1_Rapamycin_120_min",ncol_sample=ncol_sample)
heat.tko.ino80.rap.r1<-getHeatSamp(name="_TKO_Rpb1_Ino80_Rapamycin_120_min_Rep_1",ncol_sample=ncol_sample)
heat.wt.r1<-getHeatSamp(name="_Wild_type_Rep_1",ncol_sample=ncol_sample)
heat.tko.r1<-getHeatSamp(name="_TKO_Rep_1",ncol_sample=ncol_sample)

idx.selection<-which(pos<0 & pos>-300)

reg.tko.rap<-apply(heat.tko.rap[,idx.selection],1,sum)
reg.tko.ino80.rap.r1<-apply(heat.tko.ino80.rap.r1[,idx.selection],1,sum)
reg.wt.r1<-apply(heat.wt.r1[,idx.selection],1,sum)
reg.tko.r1<-apply(heat.tko.r1[,idx.selection],1,sum)

medianReg.tes<-apply(cbind(reg.tko.rap,reg.tko.ino80.rap.r1,reg.wt.r1,reg.tko.r1),1,median)

names(medianReg.tes)<-geneList$name

sortGenesReg<-sort(medianReg.tes, decreasing = T)

#import original gene list
genesStart<-read.delim("TES_nuc.bed", header=F)
colnames(genesStart)<-c("chr", "start", "end", "name", "nucPositioning","strand")
genesSorted<-genesStart[match(names(sortGenesReg),genesStart$name),]
genesSorted$nucPositioning<-sortGenesReg

write.table(genesSorted, file ="TES_nuc_regSorted.bed",
            col.names = F, row.names=F, quote=F, sep="\t")
