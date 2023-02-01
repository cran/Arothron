#' permutangle
#'
#' Create palettes from an image
#' @param mat array: rgb array 
#' @param var numeric: desidered resize factor
#' @param group1 logical: if TRUE each color is counted once
#' @param group2 logical: if TRUE (color) variables are scaled 
#' @param scale numeric: desidered number of clusters (i.e., number of palettes)
#' @param iter numeric: length of the color vector of each palette
#' @param cex1 numeric: desidered PLS axis
#' @param cex2 numeric: size of colored squares
#' @param cex3 numeric: size of color names
#' @param cex4 numeric: size of color names
#' @param labels numeric: size of color names
#' @param pch1 numeric: size of color names
#' @param pch2 numeric: size of color names
#' @param pch3 numeric: size of color names
#' @param col1 numeric: size of color names
#' @param col2 numeric: size of color names
#' @return angle list: color palettes arranged in a list
#' @return permangles list: color palettes arranged in a list
#' @return angle list: color palettes arranged in a list
#' @return iterangles list: color palettes arranged in a list
#' @return p-value list: color palettes arranged in a list
#' @return PCA_angle list: color palettes arranged in a list
#' @return PCA_interangles list: color palettes arranged in a list
#' @return PCA_p-value list: color palettes arranged in a list
#' @author Antonio Profico
#' @examples
#' \dontrun{
#' require(shapes)
#' require(Morpho)
#' data("gorf.dat")
#' data("gorm.dat")
#' Array<-bindArr(gorf.dat,gorm.dat,along=3)
#' CS<-apply(Array,3,cSize)
#' Sex<-c(rep("F",dim(gorf.dat)[3]),rep("M",dim(gorm.dat)[3]))
#' 
#' #Shape and size space
#' AllTrajFB<-permutangle(procSym(Array,scale=FALSE,CSinit = FALSE)$PCscores,
#' var=CS,group1=which(Sex=="F"),group2=which(Sex=="M"),scale=FALSE,iter=50)
#' hist(AllTrajFB$iterangles,breaks = 100,xlim=c(0,90))
#' abline(v=AllTrajFB$angle,lwd=2,col="red")
#' hist(AllTrajFB$PCA_iterangles,breaks = 100,xlim=c(0,90))
#' abline(v=AllTrajFB$PCA_angle,lwd=2,col="red")
#' 
#' #Shape space
#' AllTrajFB<-permutangle(procSym(Array)$PCscores,
#' var=CS,group1=which(Sex=="F"),group2=which(Sex=="M"),scale=FALSE,iter=50)
#' hist(AllTrajFB$iterangles,breaks = 100,xlim=c(0,90))
#' abline(v=AllTrajFB$angle,lwd=2,col="red")
#' hist(AllTrajFB$PCA_iterangles,breaks = 100,xlim=c(0,90))
#' abline(v=AllTrajFB$PCA_angle,lwd=2,col="red")
#' }
#' @export

permutangle<-function(mat,var,group1,group2,scale=FALSE,iter=100,cex1=range01(var[group1]+1),
cex2=range01(var[group2]+1),cex3=0.7,cex4=1.2,labels=c("stgr1","stgr2","endgr1","endgr2"),pch1=19, pch2=19,pch3=19,
col1="red",col2="blue"){
  
  range01 <- function(x) {
    (x - min(x))/(max(x) - min(x))
  }
  
  predict.var<-function(mat,var,value){ 
    mylm<-lm(mat~var)
    new <- data.frame(var = value)
    pred<-predict(mylm,new)
    return(pred)
  }
  set1<-as.matrix(mat[group1,])
  set2<-as.matrix(mat[group2,]) 
  sett<-rbind(set1,set2) 
  indep<-c(var[group1],var[group2])
  
  pca<-prcomp(sett,scale. = scale,center = TRUE)
  values <- 0
  eigv <- pca$sdev^2
  values <- eigv[which(eigv > 1e-16)]
  lv <- length(values)
  PCs <- pca$rotation[, 1:lv]
  PCscores <- as.matrix(pca$x[, 1:lv])
  Variance <- cbind(sqrt(eigv), eigv/sum(eigv), cumsum(eigv)/sum(eigv)) * 100
  Variance <- Variance[1:lv, ]
  
  y<-indep[c(1:dim(set1)[1])]
  x<-PCscores[c(1:dim(set1)[1]),]
  coeff1<-lm(x~y)[[1]][2,]
  
  y<-indep[c((dim(set1)[1]+1):dim(sett)[1])]
  x<-PCscores[c((dim(set1)[1]+1):dim(sett)[1]),]
  coeff2<-lm(x~y)[[1]][2,]
  
  Radangle<-angleTest(coeff1,coeff2) 
  Angles<-Radangle$angle*180/pi
  
  tot<-sett
  halfspecs<-round(dim(tot)[1]/2)
  angles<-NULL
  angles2<-NULL
  for(i in 1:iter){
    group1i<-sample(1:dim(tot)[1],halfspecs)
    group2i<-(1:dim(tot)[1])[-group1i]
    
    y<-indep[group1i]
    x<-PCscores[group1i,]
    coeffg1<-lm(x~y)[[1]][2,]
    
    y<-indep[group2i]
    x<-PCscores[group2i,]
    coeffg2<-lm(x~y)[[1]][2,]
    
    radangle<-angleTest(coeffg1,coeffg2)
    angles[i]<-radangle$angle*180/pi
    
    sFi<-predict.var(as.matrix(PCscores[group1i,]),indep[group1i],min(indep[group1i]))
    bFi<-predict.var(as.matrix(PCscores[group1i,]),indep[group1i],max(indep[group1i]))
    sMi<-predict.var(as.matrix(PCscores[group2i,]),indep[group2i],min(indep[group2i]))
    bMi<-predict.var(as.matrix(PCscores[group2i,]),indep[group2i],max(indep[group2i]))
    Rotsi<-rbind(sFi,sMi,bFi,bMi)
    Rots_pcai <- prcomp(Rotsi, center = TRUE,scale=scale)
    matangle<-Rots_pcai$x[c(1,3),1:2]
    matangle[,1]<-matangle[,1]-matangle[1,1]
    matangle[,2]<-matangle[,2]-matangle[1,2]
    matangle2<-Rots_pcai$x[c(2,4),1:2]
    matangle2[,1]<-matangle2[,1]-matangle2[1,1]
    matangle2[,2]<-matangle2[,2]-matangle2[1,2]
    Radangle<-angleTest(matangle,matangle2)
    angles2[i]<-Radangle$angle*180/pi
  }
  
  perpvalue<-((length(which(angles>=Angles))+1)/length(angles))
  
  sF<-predict.var(as.matrix(set1),indep[c(1:dim(set1)[1])],min(indep[c(1:dim(set1)[1])]))
  bF<-predict.var(as.matrix(set1),indep[c(1:dim(set1)[1])],max(indep[c(1:dim(set1)[1])]))
  sM<-predict.var(as.matrix(set2),indep[c((dim(set1)[1]+1):dim(sett)[1])],min(indep[c((dim(set1)[1]+1):dim(sett)[1])]))
  bM<-predict.var(as.matrix(set2),indep[c((dim(set1)[1]+1):dim(sett)[1])],max(indep[c((dim(set1)[1]+1):dim(sett)[1])]))
  
  Rots<-rbind(sF,sM,bF,bM)
  Rots_pca <- prcomp(Rots, center = TRUE,scale=scale)
  eigv <- Rots_pca$sdev^2
  
  Variance <- cbind(sqrt(eigv), eigv/sum(eigv), cumsum(eigv)/sum(eigv)) * 100
  matangle<-Rots_pca$x[c(1,3),1:2]
  matangle[,1]<-matangle[,1]-matangle[1,1]
  matangle[,2]<-matangle[,2]-matangle[1,2]
  matangle2<-Rots_pca$x[c(2,4),1:2]
  matangle2[,1]<-matangle2[,1]-matangle2[1,1]
  matangle2[,2]<-matangle2[,2]-matangle2[1,2]
  Radangle<-angleTest(matangle,matangle2)
  Angles_plot<-Radangle$angle*180/pi
  
  perpvalue2<-((length(which(angles2>=Angles_plot))+1)/length(angles2))
  
  pred <- predict(Rots_pca, newdata=sett)
  regr1<-indep[c(1:dim(set1)[1])]
  regr2<-indep[c((dim(set1)[1]+1):dim(sett)[1])]
  plot(Rots_pca$x,asp=1,xlim=range(pred[,1]),ylim=range(pred[,2]),pch=pch3,cex=cex3,
       xlab=paste("PC1 (",round(Variance[1,2],2),"%)",sep=""),
       ylab=paste("PC2 (",round(Variance[2,2],2),"%)",sep=""))
  text(Rots_pca$x,labels=labels,pos=1,cex=cex4)
  arrows(Rots_pca$x[1,1],Rots_pca$x[1,2],
         Rots_pca$x[3,1],Rots_pca$x[3,2],lwd=2)
  arrows(Rots_pca$x[2,1],Rots_pca$x[2,2],
         Rots_pca$x[4,1],Rots_pca$x[4,2],lwd=2)
  abline(v=0,h=0)
  points(pred[c(1:dim(set1)[1]),],col=col1,cex=cex1,pch=pch1)
  points(pred[c((dim(set1)[1]+1):dim(sett)[1]),],col=col2,pch=pch2,cex=cex2)
  
  return(list("angle"=Angles,"iterangles"=angles,"p-value"=perpvalue,
              "PCA_angle"=Angles_plot,"PCA_iterangles"=angles2,"PCA_p-value"=perpvalue2))
}

