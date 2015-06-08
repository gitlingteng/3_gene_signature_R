library(gplots)	
library(caret)	
library(klaR)
library(limma)
library("combinat")
 # library(car)
library(scatterplot3d)
Plot.Pca.Heatmap =function(plotdata,k,stitle)
  {
    pcresult =pcaplots(plotdata,k,stitle,FALSE,TRUE)
    myheatmap(plotdata,k,stitle)
    return(pcresult)
  }

##=========================================================================================
###pcadta:dataframe to be PCA analyzed,ptitle:name to add astitle, k,colnum of status
#pdf(file=fileName,paper="a4r",width =15,height=7, fonts="Times")  dev.off()  

###PCA Analysis(against ROWS)##################
###############################################
##pcadta:dataframe to be PCA analyzed,ptitle:name to add astitle, k,colnum of status
#pdf(file=fileName,paper="a4r",width =15,height=7, fonts="Times")  dev.off()  
                                      
                                        #========================
                                        #pcaplots()
                                        #=========================
pcaplots=function(pcadata,k,ptitle,showlabel,show3d){
  

    ## k=length(t.myExprs)
  ## pcadata =t.myExprs
  ## ptitle="RG2011"
  ## showlabel =1
  shape.choice =c(17,4,13,9,18,9,10,3,5)
                                        # pr = prcomp(pcadata[,-k],scale=TRUE)  
  pr = prcomp(pcadata[,-k],retex=TRUE)  
                                        #textplot(capture.output(pcadata),halign="center",  valign="top",cex=0.55, show.rownames =TRUE)			    
                                        # textplot(capture.output(summary(pr)),halign="left",  valign="top",cex=0.7, show.rownames =TRUE)	      
                                        # plot(pr) 
                                        # biplot(pr,var.axes = TRUE,cex =0.4,arrow.len =0.03)       		   					        

  ##against samples
  preddata =as.data.frame(pr$x)
  pc1 <- preddata[, 1]
  pc2 <- preddata[, 2]
  if(length(k)==1)  preddata$Status1 =pcadata[[k]]  else {preddata$Status1 =pcadata[[k[1]]]; preddata$Status2 =pcadata[[k[2]]]}
                                        #textplot(capture.output(preddata),halign="center",  valign="center", show.rownames =TRUE)

                                     
                                        #Define shapes and colors of the plots based on the Status
  if(length(k)==1) index=1 else index=sapply(preddata$Status2, pmatch, levels(preddata$Status2) ) 
  shapes=shape.choice[index]
  
  colors =sapply(preddata$Status1 ,pmatch, levels(factor(preddata$Status1)   ) )+1
  colors.fac=factor(colors)

                                        #make 2d plot
  if(length(k)==1) legendlist= levels(droplevels(factor(preddata$Status1)))  else legendlist= c(levels(droplevels(factor(preddata$Status1))),  levels(droplevels(factor(preddata$Status2))))
  plot(pc1,pc2,main=ptitle, xlab="pc1",ylab="pc2",asp = 1,pch =shapes,col=colors, cex = 0.8)
  legend(x="topleft",cex=1.0, legend=legendlist, col=c(levels(colors.fac),rep(1,length(unique(shapes)))),pch=c(rep(15,length(unique(colors))),unique(shapes)))

  if(showlabel) text(pc1,pc2,col="grey", rownames(pcadata), cex = 1   )
  
  ##make 3d plot  
  if(show3d) {
          #Define pc1,pc3,pc3
    
        pc3 <- preddata[, 3]

       s3d <- scatterplot3d(pc3,pc1,pc2, cex.symbols=0.6, pch=shapes,color = colors, main=paste("3D Scatterplot of",ptitle,sep=" "))
  legend(x="topleft",cex=1.0, legend=legendlist, col=c(levels(colors.fac),rep(1,length(unique(shapes)))),pch=c(rep(15,length(unique(colors))),unique(shapes)))
     #  if(interactive() && require(rgl) && require(mgcv)){
                                       # #scatter3d(PC1 ~ PC2 + PC3 | Status,revolutions=3, data=preddata,main =ptitle)
                                    #    print(c(preddata$PC1,preddata$PC2,preddata$PC3))
                                     #   scatter3d(PC1 ~ PC2 + PC3 | Status, data=preddata,main =ptitle)
                                      #  rgl.bringtotop()
                                      #  rgl.viewpoint(0,20)
                                       # rgl.snapshot(paste(ptitle,".png",sep=""))
                                      #  Sys.sleep(5) # wait 5 seconds
                                    #    }
      # scatter3d(PC1 ~ PC2 + PC3 | Status, data=preddata,main =ptitle) 
     #  rgl.snapshot(paste(ptitle,".png",sep=""))
     #  rgl.text(ptilte)
                                        
  }                                  
 

                                        #against column namess
  a1 <- pr$rotation[, 1]
  a2 <- pr$rotation[, 2]
                                        # plot(a1,a2,main="PCA analysis of variables' importance ",col = 14, pch = 16, asp = 1, xlab = "PC      1",ylab = "PC 2", cex = 0.8)
                                        # colnames(pcadata)
                                        # text(a1+0.02,a2, colnames(pcadata), cex = 0.60)
  loading =pr$rotation
                                        #sort the matrix by PC1 vallue
  newloadpc1=loading[order(loading[,1],decreasing=TRUE), ][,1]
  newloadpc2=loading[order(loading[,2],decreasing=TRUE), ][,2]
  dir.create("outputs")                                      #export rotation  to csv file
  ## filen=paste(ptitle,".csv",sep="")
  ## write.table(newloadpc1,file=paste("./outputs/pca_load_pc1_",filen,sep=""),quote=FALSE,row.names=TRUE,col.names=FALSE,sep=",")
  ## write.table(newloadpc2,file=paste("./outputs/pca_load_pc2_",filen,sep=""),quote=FALSE,row.names=TRUE,col.names=FALSE,sep=",")
  return(pr$rotation[,(1:2)]) 
}
                                        
                                        


		
#######end of function PCA######## 
##=====================================================================================

##========================Normalize function mapstd=====================================
mapstd <- function(X, kmean, ksd) {
  if(missing(kmean)||missing(ksd)){
    kmean <- 0
    ksd <- 1
  }
  return((X-mean(X))*(ksd/sd(X))+kmean)
}

##=====================================================================================
##hirarichal cluster

#hiercluster(data,k,title)
hiercluster=function(data,k,title) {
  d.t=dist(data[,-k])
  labelnames =paste(rownames(data),data[[k]],sep="_")
  plot(hclust(d.t), main=paste("Hierarchical clustering for ",title,sep=""), labels=labelnames)
  }


##=====================================================================================
##heatmap###

                                        #myheatmap(heatdata,k,ptitle)

myheatmap=function(heatdata,k,ptitle){
  heatdata=heatdata[order(heatdata[[k]]),]
  x =as.matrix(heatdata[,-k])
                                        #textplot("Heat map",cex =1)
  fac=factor(heatdata[[k]]);fac
  vec = rainbow(nlevels(fac),start=.1,end=.7);vec
  index =sapply(heatdata[[k]], pmatch, levels(fac))
  vec.color=vec[index]
  legendtitle="Status"
                                        #heatmap with out row dendrogram
  ##    #heatmap with no row dendrogram
  ## hv <- heatmap.2(x, col=greenred,Rowv=FALSE,
  ##                 RowSideColors= vec.color, 
  ##                 main=paste("heatmap",ptitle,sep=" "),
  ##                 tracecol=NULL
  ##                 )
  ## y <- legend(x="bottomleft",legend=levels(fac), box.col="grey",fill=vec,title=legendtitle,cex=1.0)
  ##heatmap with  row dendrogram
  hv.02 <- heatmap.2(x, col=greenred,
                     RowSideColors= vec.color, 
                     main=paste("heatmap",ptitle,sep=" "),
                     tracecol=NULL,dendrogram="column"
                     )
  y.02 <- legend(x="bottomleft",1,legend=levels(fac), box.col="grey",fill=vec,title=legendtitle,cex=0.45)

}

##=====================================================================================
##heatmap for gene differentiation analysis
#======================================================================================

myheatmap.diff=function(heatdata,k,ptitle){
                                        # heatdata =seldata;k=length(seldata);ptitle=""
  heatdata=heatdata[order(heatdata[[k]]),]
                                        #textplot("Heat map",cex =1)
  fac=factor(heatdata[[k]]);fac
  unique(heatdata[[k]])
  vec = rainbow(nlevels(fac),start=.1,end=.7);vec
  index =sapply(heatdata[[k]], pmatch, levels(fac))
  vec.color=vec[index]
  legendtitle="Status"
                                        #heatmap with out row dendrogram
  ##    #heatmap with no row dendrogram
  ## hv <- heatmap.2(x, col=greenred,Rowv=FALSE,
  ##                 RowSideColors= vec.color, 
  ##                 main=paste("heatmap",ptitle,sep=" "),
  ##                 tracecol=NULL
  ##                 )
  ## y <- legend(x="bottomleft",legend=levels(fac), box.col="grey",fill=vec,title=legendtitle,cex=1.0)
  ##heatmap with  row dendrogram
  ##Calculate mean,std, construct a new matrix for plot
  x =as.matrix(heatdata[,-k])
  V=t(as.matrix(heatdata[,-k]))
                                        #col.labels = 0
  n.rows <- length(V[,1])
  n.cols <- length(V[1,])
  row.mean <- apply(V, MARGIN=1, FUN=mean)
  row.sd <- apply(V, MARGIN=1, FUN=sd)
  row.n <- length(V[,1])
  for (i in 1:n.rows) {
    if (row.sd[i] == 0) {
      V[i,] <- 0
    } else {
      V[i,] <- (V[i,] - row.mean[i])/(0.5 * row.sd[i])
    }
    V[i,] <- ifelse(V[i,] < -6, -6, V[i,])
    V[i,] <- ifelse(V[i,] > 6, 6, V[i,])
  }


  mid.range.V <- mean(range(V)) - 0.1
  ## heatm <- matrix(0, nrow = n.rows, ncol = n.cols)
  ## heatm[1:n.rows,] <- V[seq(n.rows, 1, -1),]
  ## heatm[n.rows + 1,] <- ifelse(col.labels == 0, 7, -7)

  ##ready for heatmap
  hv.01 <- heatmap.2(V, col=bluered,
                     ColSideColors= vec.color,key=FALSE,
                     main=paste("Heatmap For",ptitle,sep=" "),
                     tracecol=NULL, dendrogram="none",cexRow=0.8
                     )  #Rowv=FALSE,Colv=FALSE,cexCol =0.5,key=FALSE,
  statuslabel =paste(unique(heatdata[[k]]),collapse="                    ")

  ## mtext(statuslabel  , 3, line=-5)
  ## axis(3, at=c(1,4), labels=col.classes, tick=FALSE, las = 1, cex.axis=1.25, font.axis=2, line=-1)
                                        # heatmap.2( (heatm),col=bluered,Rowv=FALSE, Colv=FALSE, dendrogram="none",tracecol=NULL) #
  y.02 <- legend(x="topright",-4,legend=levels(fac), box.col="grey",fill=vec,title=legendtitle,cex=0.45)
   return(V)
}

###===========Heatmap for simple differentiation==========================
myheatmap.simdiff =function(seldata) {    #colnames are sample names, rows are gene feartures
  V=as.matrix(seldata)
  n.rows <- length(V[,1])
  n.cols <- length(V[1,])
  row.mean <- apply(V, MARGIN=1, FUN=mean)
  row.sd <- apply(V, MARGIN=1, FUN=sd)
  row.n <- length(V[,1])
  for (i in 1:n.rows) {
    if (row.sd[i] == 0) {
      V[i,] <- 0
    } else {
      V[i,] <- (V[i,] - row.mean[i])/(0.5 * row.sd[i])
    }
    V[i,] <- ifelse(V[i,] < -6, -6, V[i,])
    V[i,] <- ifelse(V[i,] > 6, 6, V[i,])
  }

  heatmap.2(V, col=bluered,key=FALSE,Rowv=FALSE, main="Simple Heatmap",tracecol=NULL, dendrogram="none",cexRow =0.1 )
  }
#=========================================================================================
#function :plot.pca.heat.density(final.data) make pca plot, heatmap , and density plot, suppose status is last column

plot.pca.heat.density =function(final.data) {

  k=length(final.data)
  pcaplots(final.data,k,"select_ovarian",FALSE)
  myheatmap(final.data,k,"select_ovarian")
                                        #plot the density of sensi and resist respectively
  test=final.data
  k=length(test)
  test$status =as.factor(as.character(test$status))
  str(test)
  for(i in 1:(length(test)-1)) {
                                        #i=1
    x=test[[i]]
    name =names(test)[i]
    relist01=aggregate(x,by=list(test[[k]]),FUN=function(y) return(y))
    len =nrow(relist01)
    colors=c("green","red","blue","yellow")
    for(i in 1:len)
      {###i=1   as.numeric(relist[i,][[2]]) 
        if(i==1) plot(density( unlist(relist01[i,][[2]])),col=colors[i], xlab="Expression value",main= name ) else
        {
          r=density( unlist(relist01[i,][[2]]))
          points(r$x,r$y,col=colors[i],pch=".")
        }      
      }
    legend(x="topright",cex=1  , legend=as.character(relist01[,1]),col= colors[1:len] ,pch=15) 
  }
}




##====================================================================================
##Supervised Classifier scoring algorithm
##=============function: runclassifier================================================
##====================================================================================

runclassifier =function(nbdata,newlist,classifier,k){
#library(caret)	
#library(klaR)
    ## nbdata =origdata;newlist =torunlist01;classifier="nb"
 #   classifier="svmRadial";classifier="classifier="lvq"
   # k=length(nbdata)
  #  classifier="nb"
    acculistcv=list()     
    for(j in 1:length(newlist)) {
    #  print(j)                                        #j=1
      inputdata=nbdata[,unlist(c(k,newlist[[j]]))]#;str(inputdata)
      
      fitcv<- train(inputdata[,-1],factor(inputdata[[1]]),classifier ,trControl = trainControl(method = "cv", number=nrow(inputdata)))      
      acculistcv[[j]]=max(fitcv$results$Accuracy)     
    }
    fitcv$results$Accuracy
    acculistcv = unlist(acculistcv)
    order.acculist = rev(order(acculistcv))
    rev(order.acculist )
   # plot(density(order.acculist))  
    matchid01=newlist[order.acculist ]  
    modelresult= list(ids=matchid01, value=acculistcv[order.acculist])
    return(modelresult)
  }

  ##===========function classifier.nb
  ##===========================
## runlist0202=combin02.result[[1]][1:60]
## classify.list=list()                                        #  combination2
## m=length(classify.list)
## for(j in 1:60)
##   {
##     classify.list[[m+j]]=list(id=combin02.result[[1]][[j]],value =combin02.result[[2]][[j]] )
##   }
## finalresult =classifier.nb(origdata,runlist0202,k,classify.list,5) --5 means maximam length
##runlist03=lapply(finalresult[22:31],function(x) x$id)
##finalresult03 =classifier.nb(origdata,runlist03,k,finalresult,29)
##runlist04=lapply(finalresult03[1:20],function(x) x$id)


##valuelist04=lapply(finalresult03,function(x) x$value)
##valueorder = order(unlist(valuelist04),decreasing =TRUE)
##torunlist04.reorder=lapply(finalresult03[valueorder],function(x) x$id)#6,13,14
  
##  save.image("nano.rdata")
##plot of combinations based on 17,40

#pdf("PCA plot of nanostring data with all samples and all cell lines of selected own genes03.pdf")
#for(i in 1:length(torunlist04.reorder)) {
    #selnano2 = nanodata.copy[, c(torunlist04.reorder[[i]],54) ] #
    #subtitle =paste("Nanostring all sample with cell lines of selected own genes ",i,sep="--")

#pca_plots(selnano2 ,length(selnano2),subtitle ,TRUE)
      
#}
#dev.off()

##*************************************************************##
#testfinal01 = classifier.nb(origdata,list(),k,list(),2,100)  #might take 24 hours #6000sets/hour,
    ## runlist =list(); classify.list=list()
    ## origdata =  highsel.exprs;k=length(highsel.exprs);maxilen=2; selnum =100(normally, more than 100 will take longer)
##
#runlist02=lapply(testfinal01[380:479],function(x) x$id)  
#testfinal01 = classifier.nb(origdata,runlist02,k,testfinal01,4,100)  

classifier.nb=function(origdata,runlist,k,classify.list,maxilen,selnum)
  { #6000sets/hour,
    ## runlist =list(); classify.list=list()
    ## origdata =  highsel.exprs;k=length(highsel.exprs);maxilen=2; selnum =100(normally, more than 100 will take longer)
   if(length(runlist)==0) nn=0 else nn=length(runlist[[1]])
    print("........Combination of ")
    print(nn)
    print("Numer of input set:")
    print(length(runlist))
    totaltorunlist =list()
    torunlist02 =list()
    torunlist01=list()
    
    if(nn<maxilen)
      {
        if(length(runlist)==0) 
          {
            columlist =(1:length(origdata))[-k]
            for(i in 1:length(columlist)) torunlist01[[i]]=columlist[i]    
             for(i in 1:length(torunlist01))
            {
                                        #i=1
              ids = torunlist01[[i]]
              idlist =(1:ncol(origdata))[-c(ids,k)]
              #add new id to runlist, output new  totaltorunlist to runclassifier
              torunlist =lapply(idlist,function(x) c(ids,x))
              h=length(totaltorunlist)
              totaltorunlist[(h+1):(h+length(torunlist))]=torunlist                                      
                                        #  classify.compare[[h+1]]=re.model
            }
           torunlist02 =lapply( totaltorunlist,function(x) sort(x)) 
        #  print(length( torunlist02 ))
          torunlist02=unique(torunlist02)             
          } else 
        {
          for(i in 1:length(runlist))
            {
                                        #i=1
              ids = runlist[[i]]
              idlist =(1:ncol(origdata))[-c(ids,k)]
              #add new id to runlist, output new  totaltorunlist to runclassifier
              torunlist =lapply(idlist,function(x) c(ids,x))
              h=length(totaltorunlist)
              totaltorunlist[(h+1):(h+length(torunlist))]=torunlist                                      
                                        #  classify.compare[[h+1]]=re.model
            }
                                        #clear the list, filter out duplicate ones
          torunlist02 =lapply( totaltorunlist,function(x) sort(x)) 
        #  print(length( torunlist02 ))
          torunlist02=unique(torunlist02)   
           
          
        }
        print("Numer of running set:")
        print(length( torunlist02 ))
        print(paste("Estimated running time: ", length(( torunlist02 ))/6000,"  hours",sep=""))
                                        # summary02 = runclassifier(torunlist,"nb")
        print("Begin the classifier......")
        re.model=runclassifier(origdata,torunlist02,"nb",k)
        print("End the classifier......")
                                        # re.model$ids[which(re.model$value >= 0.9)]
        
                                        #runlist02 = re.model[[1]][which(re.model$value >= 0.9)]
                                        
        if(length(re.model[[1]]) <selnum ) highnum =length(re.model[[1]]) else highnum =selnum                               
        runlist02 = re.model[[1]][1:highnum]
        m=length(classify.list)
        if(length(re.model[[1]]) <500) runs =length(re.model[[1]]) else runs =500
        for(j in 1:runs)
          {
            classify.list[[m+j]]=list(id=re.model[[1]][[j]],value =re.model[[2]][[j]] )
          }
        
                                        #   outputlist=list(runlist = runlist02,classlist=classify.list)
        resultlist= classifier.nb(origdata,runlist02,k, classify.list,maxilen,selnum)
        
                                        #   classify.list[[i]]=re.model      
      } 
    else  {return(classify.list)  } 
    
    return(resultlist)
  }

 


##=====================================================

#FUNCTION---naiveBayes,output to a pdf file.
naiveBayes_fun =function(nbdata,k,stitle )
{
  library(gplots)	
  library(caret)	
  library(klaR)
  nbdata[[k]] =as.factor(as.character(nbdata[[k]]))
  num = ncol(nbdata)-1
                                        #list of possible combination of antibody names
  newlist=list()
  for(n in 1:num) {
    mat=combinations(num,n,((1:ncol(nbdata))[-k]));mat
    for(h in 1:nrow(mat)) {name=paste(n,h,sep="num");newlist[[name]]=mat[h,]}
  };newlist
                                        #combination of LGOCV(leave 4 out), randomly selected 800 of the total number
  rownum =nrow(nbdata)
  wholelist   =combinations(rownum,4)
  partmatrix =wholelist[resample(1:nrow(wholelist),800),]
  partlist =list()
  acculist=list()
  acculistcv=list()
  for(i in 1:nrow(partmatrix)) {partlist[[i]]=(1:rownum)[-partmatrix[i,]]};partlist
                                        #train NA classifier based on partlist,which has 800 combinations
                                        #			
  for(j in 1:length(newlist)) {
    print(j)
    inputdata=nbdata[,c(k,newlist[[j]])];str(inputdata)
    
    fitcv<- train(inputdata[,-1],inputdata[[1]], "nb",trControl = trainControl(method = "cv", number=10))
    acculistcv[[j]]=max(fitcv$results$Accuracy)			
  }
                                        #get index of the antiboies' names, and the best combination data that has the highest scores	
  alist=unlist(acculistcv)  
  maxid = which(alist == max(alist));print(maxid)
  
  bestcomb=nbdata[,c(k,newlist[[maxid]])];print(bestcomb)
  bestfit <- train(bestcomb[,-1],bestcomb[[1]], "nb",trControl = trainControl(method = "LGOCV",index=partlist));bestfit	
                                        #pdf(paste(stitle ,"PCA_NaiveBayes_Heatmap.pdf",sep=" "),paper="a4",width =15,height=7, fonts="Times") 
  textplot(paste(" Naive Bayes Classifier and pca plot ,heatmap",stitle,sep=" "))
                                        # textplot("Original data",cex =1)  
                                        # textplot(capture.output(head(nbdata)), halign="center",   valign="top", show.rownames = TRUE)  
  textplot(c("Best Combination: ",names(bestcomb)[-1]),cex =1)  
# textplot("Best Combination",cex =1)
  textplot(capture.output( head(bestcomb)), halign="center", cex=0.55,  valign="top", show.rownames = TRUE)				
  
  textplot(c("Naive Bayes Classifier Summary","-----------------------------",capture.output(bestfit)), halign="center",        valign="top", show.rownames = TRUE)
                                        # finalclassifier = fit$finalModel
                                        # finalclassifier
                                        # predicttb = predict(finalclassifier, nbdata[,-k])
                                        # predicttb
                                        # ntb=table(nbdata[,k],predicttb$class)
                                        # ntb
                                        # summary(ntb)
  pcaplots( nbdata,k,stitle ,FALSE    )
  myheatmap(  nbdata,k,stitle)
   revalue=names(bestcomb)[-1]
   return(revalue)
}
                   
 #End of function#####
 
##=======================================================================================
## ========================= PLOT FUNCTIONS   ===========================================
##=======================================================================================

plot.line=function(x)
  {
    ncolor=ncol(x)
    n= nrow(x)
    colors <- rainbow(8)
    linetype <- c(1:ncolor)
    plotchar <- seq(18,18+ncolor,1)
 plot(1:n, seq(min(x),max(x),length=n), type="n", xlab="Samples",   ylab="expression" )
    for(i in 1:ncol(x))
      {
        lines(1:length(x[[i]]), x[[i]], type="b", lwd=1.5,lty=linetype[i], col=colors[i], pch=plotchar[i])
      }
    }
 #End of function#####
##=======================================================================================


##scatterplot with colored group
##=======================================================================================
#s=3
#title="plot of FZD7 and PIPOX"
#name4text = unlist(lapply(rownames(tbltoplot), function(x) unlist(strsplit(x,split="_U",fix=TRUE))[1]   ))
#groupplot(tbltoplot,s,title)

#####
groupplot =function(tbltoplot,s,title,name4text)
  {
    statuscol =tbltoplot[[s]]
    pc1=tbltoplot[,-s][[1]]
    pc2=tbltoplot[,-s][[2]]
    labnames=names(tbltoplot[,-s])
    colors =sapply(statuscol , pmatch, levels(statuscol) )+1
    colors.fac=factor(colors)
    shapes=20

    legendlist= levels(droplevels(statuscol ))
    plot(pc1,pc2,main=title,xlab=labnames[1],ylab=labnames[2], asp = 1,pch =shapes,col=colors, cex = 1.3)
    legend(x="topright",cex=1.0, legend=legendlist, col=c(levels(colors.fac),rep(1,length(unique(shapes)))),pch=c(rep(15,length(unique(colors)))))
    
    text(pc1,pc2,col="grey", name4text, cex = 0.60)	
    }

##===============================================================================
#==
densityplot =function(final.data) {

  test=final.data
  k=length(test)
  test$status =as.factor(as.character(test$status))
  str(test)
  for(i in 1:(length(test)-1)) {
                                        #i=1
    x=test[[i]]
    name =names(test)[i]
    relist01=aggregate(x,by=list(test[[k]]),FUN=function(y) return(y))
    len =nrow(relist01)
    colors=c("green","red","blue","yellow")
    for(i in 1:len)
      {###i=1   as.numeric(relist[i,][[2]]) 
        if(i==1) plot(density( unlist(relist01[i,][[2]])),col=colors[i], xlab="Expression value",main= name ) else
        {
          r=density( unlist(relist01[i,][[2]]))
          points(r$x,r$y,col=colors[i],pch=".")
        }      
      }
    legend(x="topright",cex=1  , legend=as.character(relist01[,1]),col= colors[1:len] ,pch=15) 
  }
}

##=========Normalize functions=============================
NormalizeRows <- function(V) { 
#
# Stardardize rows of a gene expression matrix
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.

        row.mean <- apply(V, MARGIN=1, FUN=mean)
               row.sd <- apply(V, MARGIN=1, FUN=sd)
        row.n <- length(V[,1])
        for (i in 1:row.n) {
             if (row.sd[i] == 0) {
                  V[i,] <- 0
           } else {
              V[i,] <- (V[i,] - row.mean[i])/row.sd[i]
           }
        }
        return(V)
}

NormalizeCols <- function(V) { 
#
# Stardardize columns of a gene expression matrix
#


        col.mean <- apply(V, MARGIN=2, FUN=mean)
               col.sd <- apply(V, MARGIN=2, FUN=sd)
        col.n <- length(V[1,])
        for (i in 1:col.n) {
             if (col.sd[i] == 0) {
                  V[i,] <- 0
           } else {
              V[,i] <- (V[,i] - col.mean[i])/col.sd[i]
           }
        }
        return(V)
}

# end of auxiliary functions


##=========functions  avedupTblRows ================

avedupTblRows =function(mergtbl,h,d) { #h(column id which need to average) d(column which non numeric and can't be averaged)

#   mergtbl =merg.genemap; h=length(merg.genemap);d=c(1)
 dupnames =unique(mergtbl[[h]][which(duplicated(mergtbl[[h]]))] )
  
  avgtable =as.data.frame(sapply(dupnames,function(x) unlist(colMeans((mergtbl[grep(x,mergtbl[[h]]),-c(d,h)]))))     )
  names(avgtable ) =dupnames
  mergtbl =  mergtbl[-which(mergtbl[[h]] %in% dupnames),]
  ##mergtbl02 with tvalue and pvalue
  str(mergtbl)
  rownames(mergtbl) =mergtbl[[h]]
  mergtbl02 =as.data.frame(t(mergtbl[,-c(h,d)]))
  mergtbl02 = as.data.frame(t(cbind(mergtbl02,avgtable)))
  return(mergtbl02)
  }

##===========Affymatrix exon gene exression analysis =========================== 

getexprs.exon.gene =function(d.gene,gene.dabg,names.sub01,names.sub02,statusname,path,prjname,geneassignment) {
        # Get DABG values for each sample
        # At least 25% of probes in each sample have dabg p-value less than 0.05
        all(apply(gene.dabg, 2, quantile)[2,] < 0.05)

        # Get DABG values for each subtype
        gene.dabg.sub01 <- gene.dabg[,colnames(gene.dabg) %in% names.sub01]
        gene.dabg.sub02 <- gene.dabg[,colnames(gene.dabg) %in% names.sub02]

        gene.dabg.quantile.sub01 <- t(apply(gene.dabg.sub01, 1, quantile))
        gene.dabg.quantile.sub02 <- t(apply(gene.dabg.sub02, 1, quantile))
        # Retain probesets with p-value less than k.exon.dabgpval in at least 50% of samples
        # When k.gene.dabgpval = 0.05, 20038 probesets are retained.
        # When k.gene.dabgpval = 0.01, 18589 probesets are retained.
        k.gene.dabgpval <- 0.05
        gene.dabg.fltr.sub01 <- rownames(gene.dabg.sub01)[gene.dabg.quantile.sub01[,3] < k.gene.dabgpval]
        gene.dabg.fltr.sub02 <- rownames(gene.dabg.sub02)[gene.dabg.quantile.sub02[,3] < k.gene.dabgpval]
        gene.dabg.fltr <- unique(c(gene.dabg.fltr.sub01 , gene.dabg.fltr.sub02 )) #20234
        #  str(gene.dabg.fltr)
        # Filter expression matrix by DABG p-values
        ##20234*56
        gene.rma.dabg <- d.gene[rownames(d.gene) %in% gene.dabg.fltr,names(d.gene) %in% c(names.sub01,names.sub02)]
  
      ##Normalize
        gene.exprs=gene.rma.dabg
   #     gene.exprs.orig=t(gene.exprs)
        gene.exprs.norm =t(apply(gene.exprs,2,mapstd))#normalize again, proved to be necessary. "german_Ewing Gene Expression PCA:after DABG filtering.pdf"
        names.full=c(names.sub01,names.sub02)
        np.full=c(rep(statusname[1],length(names.sub01)),rep(rep(statusname[2],length(names.sub02))))
        gene.exprs.sel= gene.exprs.norm[match(names.full,rownames(gene.exprs.norm)),]
        gene.exprs.sel = cbind(np.full,as.data.frame(gene.exprs.sel))  #(11*21368)
        
        ##=====================================================================
        ##STEP 3 PCA PLOT
        ##=====================================================================
      pdf(paste(path,prjname," Gene Expression PCA:after DABG filtering.pdf"))
      pcaplots(gene.exprs.sel,1,"Gene expression renormalized:after DABG filtering",1,0)
      dev.off()
      ###
    #
    # HIERARCHICAL CLUSTERING
    #
    ###

    ##hirarichal cluster
    index =sapply(gene.exprs.sel[[1]], pmatch, levels(gene.exprs.sel[[1]]) )
    label.str=statusname[index]
    d.t=dist(gene.exprs.sel[,-1])
    pdf(paste(path,prjname,"  hiercluster.pdf"))
    plot(hclust(d.t), main="Hierarchical clustering for gene:after DABG filtering", labels=label.str)
    dev.off()
    ##=====================================================================
    ##STEP 4  Differentiation Analysis
    ##=====================================================================
    library("limma")
    #str(gene.exprs.sel)
    exprs.matrix =as.matrix(t(gene.exprs.sel[,-1]))
    str(exprs.matrix)
    design <- cbind(Grp1=1,Grp2vs1=c(rep(0,length(names.sub01)),rep(1,length(names.sub02)))) # "nonrelapse" control --0,"relapse"--1 event
    arrayw <- arrayWeights(exprs.matrix)

    pdf(paste(path,prjname,"  differentiation plot.pdf"))
    barplot(arrayw, xlab="Array", ylab="Weight", col="white", las=2)
    #classic fit
    fit <- lmFit(exprs.matrix ,design)
    fit <- eBayes(fit)
   # fit
    as.data.frame(fit[1:10,2])
    qqt(fit$t[,2],df=fit$df.residual+fit$df.prior)
    abline(0,1)
    volcanoplot(fit,coef=2,highlight=100)
    dev.off()

    options(digits=2)
    survival.top = topTable(fit ,coef=2,adjust.method="BH", n=2000)  #p <0.01 2303
#str(survival.top)
    survival.bigtop =survival.top[survival.top$P.Value<=0.05,]
    #map the probeset id to gene through annotation
    #str(geneassignment)
    survival.bigtop$genename =unlist(lapply(survival.bigtop[[1]],function(x) 
    {#   x=smalltop[[1]][21]
        index =grep(x,names(geneassignment))
       if(length(index)>0) { gname = geneassignment[[grep(x,names(geneassignment))]]}  else { gname =paste("probeset_",x,sep="")}
##    print(gname)
       return(gname)
}))

    survival.bigtop=survival.bigtop[order(survival.bigtop$t),]
    survival.bigtop=survival.bigtop[-grep("probe*",survival.bigtop[[7]]),]
    survival.bigsel =gene.exprs.sel[,c(1,match(survival.bigtop[[1]],colnames(gene.exprs.sel)))]#pvalue <0.05
    names(survival.bigsel) =c("status",survival.bigtop$genename)
    survival.smallsel =gene.exprs.sel[,c(1,match(survival.bigtop[survival.bigtop$P.Value <=0.01,][[1]],colnames(gene.exprs.sel)))]  #pvalue <0.01
    names(survival.smallsel) =c("status",survival.bigtop[survival.bigtop$P.Value <=0.01,]$genename)
    pdf(paste(path,prjname,"  PCA plot and heatmap after limma differentiation analysis.pdf"))
    pcaplots(survival.bigsel,1,"PCA  <0.05",1,0)
    myheatmap(survival.bigsel,1,"Heatmap <0.05")
    pcaplots(survival.smallsel,1,"PCA <0.01",1,0)
    myheatmap(survival.bigsel,1,"Heatmap <0.01")
    dev.off()
     print("1,normalized geneexprs;2,fit result,topTable(fit ,coef=2,adjust.method=BH, n=2000),3,top 2000 differentiated genes table(control=0||-1,REL=1,UP regulate if t>0 for REL);4,top <0.05 genes expression;5,top <0.01 genes expression")
    return(list(gene.exprs.sel,fit,survival.top,survival.bigsel,survival.smallsel) )  #gene expression table, limma differentiation result
}

