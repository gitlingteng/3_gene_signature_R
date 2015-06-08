load("/home/lteng/Desktop/Projects/9Raffyd_nano/nano.rdata")
load("/home/lteng/Desktop/Projects/9Raffyd_nano/orig_data/williamson/cellfile/williamgene.rdata")
save.image("/home/lteng/Desktop/Projects/9Raffyd_nano/nano.rdata")  
source("./project_function.r")  # function file

nano1=read.delim("Nanstr_all_with.csv",header=TRUE,row.names =1,sep="\t")
sampl.status =read.delim("Sample_status.csv",header=TRUE,sep="\t")

#nano1
colnames =c("X509", "X4463","X5700", "X6152", "X10012","X10440", "X16655", "X22707", "X2574AP", "X12251AP", "RD.FFPE..ERMS.","RH28.FFPE..ARMSp.", "X2369", "X2574" , "X6237","X6468", "X10191", "X12240", "X12251","X13413", "X15438", "X24130", "IRS.12", "IRS.19", "IRS.44","IRS.58", "IRS.66" , "IRS.73" , "IRS.78" ,"IRS.83", "IRS.93", "IRS.135", "IRS.136" )
#transpose the original data frame
nano1.t=as.data.frame(t(nano1[,-1]))

#change the rownames
rownames(nano1.t) =as.character(sampl.status[[1]])
nano1.t$Fusionstatus =sampl.status[[2]]

                                    ##PCA analysis##
pdf("PCA plot of nanostring data with all samples and all cell lines.pdf")
nano1pca =pcaplots(nano1.t,54,"Nanostring all sample with cell lines",TRUE)

dev.off()
#select genes 
nano1gene =nano1pca[order( nano1pca[,1],decreasing=TRUE), ]
nano1gene02 =nano1pca[order( nano1pca[,2],decreasing=TRUE), ]
genelist =rownames(nano1gene)
genelist02 =rownames(nano1gene02)
selgenes =genelist[seq(1,5,step=1)]
selgenes02 =genelist02[seq(1,5,step=1)]
#wantgenes =c(selgenes,selgenes02)

wilsongenes =as.character(na.omit(read.delim("wilsongenes.txt",header=FALSE,skip=0))[,1])
 
wantgenes =wilsongenes[c(seq(1,5))]
wantgenes02 =wantgenes[c(-1,-2)]

selnano1 = nano1.t[,c(which(colnames(nano1.t)[-54] %in% wantgenes02),54) ]
pdf("PCA plot of nanostring data with all samples and all cell lines of selected own genes03.pdf")
pca_plots(selnano1 ,length(selnano1),"Nanostring all sample with cell lines of selected own genes",TRUE)

dev.off()


seq(1,31,step=1)
##nano2
pdf("PCA plot of nanostring data with all samples  without cell lines.pdf")
pca_plots(nano2.t,54,"Nanostring all sample without cell lines",TRUE)
dev.off()
##nano3
pdf("PCA plot of nanostring data with annotated samples  with cell lines.pdf")
pca_plots(nano3.t,54,"Nanostring annotated sample with cell lines",TRUE)
dev.off()
##nano4
pdf("PCA plot of nanostring data with annotated samples  without cell lines.pdf")
pca_plots(nano4.t,54,"Nanostring annotated sample without cell lines",TRUE)
dev.off()


names427=c("PIPOX","FZD7","RBM13","FBN2","SPRED2","HMGA2"  ,"TFAP2B" )
totallist =list(v1=c("PIPOX","FZD7"),v2=c("PIPOX","RBM13"),v3=c("PIPOX","SPRED2"),v4=c("TFAP2B","FBN2"),v5=c("PIPOX","HMGA2"))
totallist[[1]]
                       ##########################
                       #      ##FUNCTIONS##     #
                       ##########################
                                        #pcadata =iris2; ptitle ="iris2";k=c(5,6)
pcadata =iris; ptitle ="iris";k=5;showlabel=FALSE
pca_plots(iris,5,"iris",FALSE)
#dev.off()

                                       # #================================================
                                          ##  final code,after change the status of "10191"
                                       # #================================================

nanodata.copy =nano1.t[order(nano1.t$Fusionstatus),]
nanodata.copy[which(rownames(nanodata.copy) %in% "10191"),]$Fusionstatus ="Fusion Negative"
selnano = nanodata.copy
selnanostatus = selnano[-which(selnano$Fusionstatus =="Unknown"),]

pdf("PCA plot of nanostring samples with status and all genes.pdf")
pca_plots(selnanostatus ,length(selnanostatus),"Nanostring sample with status and all genes",TRUE)

selnano2gene =selnanostatus[,c(17,40,54)]
str(selnano2gene)
pca_plots(selnano2gene ,length(selnano2gene),"Nanostring sample with status and all genes",TRUE)
tbltoplot=selnano2gene
s=3
groupplot(tbltoplot,s,"Nanostring sample with status and  gene FZD7 and PIPOX")
dev.off()

s=3
title="plot of FZD7 and PIPOX"
groupplot(tbltoplot,s,title)

#####
pdf("PCA plot of nanostring data with all samples and all cell lines of all genes.pdf")
pca_plots(selnano1 ,length(selnano1),"Nanostring all sample with cell lines of all genes",TRUE)
dev.off()

str(nanodata.copy[,c(6,10,17,27,31,40,46,54)])
selnano01 = nanodata.copy[,c(17,40,54)]  ##11111
pdf("PCA plot of nanostring data with all samples and own selected genes.pdf")
myheatmap(selnano01,3,"heatmap for FZD7 PIPOX ")
pca_plots(selnano01 ,length(selnano01),"Nanostring all sample with cell lines of own selected FZD7 PIPOX",TRUE)
dev.off()

selnano1 = nanodata.copy[,c(3,4,8,13,5,40,49,54)]  ##11111
metaHeatmap(res)
heatmap.2(selectedgenes, col=greenred)
pdf("PCA plot of nanostring data with all samples and NMF selected genes.pdf")
myheatmap(selnano1,8,"heatmap for  NMFselected genes  ")
pca_plots(selnano1 ,length(selnano1),"Nanostring all sample with cell lines of NMFselected genes",TRUE)
dev.off()

selnano2 = nanodata.copy[,c(14,21,54)]  ##11111
pdf("PCA plot of nanostring data with all samples and NMF selected genes.pdf")

myheatmap(selnano1,8,"heatmap for  NMFselected genes  ")
pca_plots(selnano1 ,length(selnano1),"Nanostring all sample with cell lines of NMFselected genes",TRUE)
dev.off()


myheatmap(th,5,"Fig 1 heatmap for 4 biological process")
myheatmap(selnano1,length(selnano1),"Fig 2 heatmap for  NMFselected genes(Bio process 4)  ")
pca_plots(selnano1 ,length(selnano1),"Fig 3 Nanostring all sample with NMFselected genes(Bio process 4)",TRUE)
myheatmap(selnano2,length(selnano2),"Fig 4 heatmap for  NMFselected genes(Bio process 2)  ")
pca_plots(selnano2 ,length(selnano2),"Fig 5 Nanostring all sample with NMFselected genes(Bio process 2)",TRUE)

dev.off()



######################
##run the NB Gene Classifier
######################
                                        #nbdata.copy02=nbdata.copy[,c(1,2,3,4, 5,6, 7,8,109)]  

  origdata=nanodata.copy
  k=length(origdata);
#  beginlist =which(colnames(origdata) %in% wantgenes02)
##
                                        # # torunlist01=list(); beginlist=c(17,40); torunlist01[[1]]=beginlist
  
  classify.list=list()
   totaltorunlist=list()
  runlist=(1:k)[-k]
  for(i in 1:length(runlist))
    {
                                        #i=1
      ids = runlist[[i]]
      idlist =(1:ncol(origdata))[-c(ids,k)]
      torunlist =lapply(idlist,function(x) c(ids,x))
      h=length(totaltorunlist)
      totaltorunlist[(h+1):(h+length(torunlist))]=torunlist                                      
    }
                                        #clear the list, filter out duplicate ones
  torunlist01 =lapply( totaltorunlist,function(x) sort(x)) 
  torunlist01=unique(torunlist01); str( torunlist01)
  combination2.result=runclassifier(origdata,torunlist01,"nb")
  runlist02=combination2.result[[1]][1:60]
                                        #  combination2.result=re.model
  m=length(classify.list)
  for(j in 1:60)
    {
      classify.list[[m+j]]=list(id=combination2.result[[1]][[j]],value =combination2.result[[2]][[j]] )
    }
  save.image("nano.rdata")

                                        #classify.list=list()
  finalresult =classifier.nb(origdata,runlist02[1],k,classify.list[1],5)
  runlist03=lapply(finalresult[22:31],function(x) x$id)
  finalresult03 =classifier.nb(origdata,runlist03,k,finalresult,29)
  runlist04=lapply(finalresult03[1:20],function(x) x$id)


valuelist04=lapply(finalresult03,function(x) x$value)
valueorder = order(unlist(valuelist04),decreasing =TRUE)
torunlist04.reorder=lapply(finalresult03[valueorder],function(x) x$id)#6,13,14

##plot of combinations based on 17,40
pdf("PCA plot of nanostring data with all samples and all cell lines of selected own genes03.pdf")
for(i in 1:length(torunlist04.reorder)) {
    selnano2 = nanodata.copy[, c(torunlist04.reorder[[i]],54) ] #
    subtitle =paste("Nanostring all sample with cell lines of selected own genes ",i,sep="--")

pca_plots(selnano2 ,length(selnano2),subtitle ,TRUE)
      
}
dev.off()

