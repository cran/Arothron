#' image2palettes
#'
#' Create palettes from an image
#' @param array array: rgb array 
#' @param resize numeric: desidered resize factor
#' @param unique logical: if TRUE each color is counted once
#' @param scale logical: if TRUE (color) variables are scaled 
#' @param k numeric: desidered number of clusters (i.e., number of palettes)
#' @param lcols numeric: length of the color vector of each palette
#' @param plsaxis numeric: desidered PLS axis
#' @param cex numeric: size of colored squares
#' @param cext numeric: size of color names
#' @return paletteslist list: color palettes arranged in a list
#' @author Antonio Profico
#' @examples
#' \dontrun{
#' require(jpeg)
#' require(Morpho)
#' data("Altapic")
#' image2palettes(Altapic,resize=1,unique=T,scale=T,k=3,lcols=5,plsaxis=1,cext=0.5)
#' }
#' @export

image2palettes<-function(array,resize=4,unique=FALSE,scale=F,k=3,lcols=7,
                         plsaxis=1,cex=5,cext=0.5)
{
nrow<-dim(array)[1]
ncol<-dim(array)[2]
imageR<-array[round(seq(1,nrow,length.out=nrow/resize)),
              round(seq(1,ncol,length.out=ncol/resize)),]
imageR[,,1]<-imageR[dim(imageR)[1]:1,dim(imageR)[2]:1,1]
imageR[,,2]<-imageR[dim(imageR)[1]:1,dim(imageR)[2]:1,2]
imageR[,,3]<-imageR[dim(imageR)[1]:1,dim(imageR)[2]:1,3]

cols<-(t(vecx(imageR)))
if(isTRUE(unique)) cols<-unique(cols)
rgbtoB<-function(rgbval){
out<-as.numeric(sum(rgbval*c(0.2126,0.7152,0.0722)))
}
RGB2XYZ<-function(mat){
  # array<-image
  # dim(array)
  selr1<-which(mat[,1]>0.04045,arr.ind = T)
  selg1<-which(mat[,2]>0.04045,arr.ind = T)
  selb1<-which(mat[,3]>0.04045,arr.ind = T)
  selr2<-which(mat[,1]<=0.04045,arr.ind = T)
  selg2<-which(mat[,2]<=0.04045,arr.ind = T)
  selb2<-which(mat[,3]<=0.04045,arr.ind = T)
  mat[,1][selr1]<-((mat[,1][selr1] + 0.055 ) / 1.055 ) ^ 2.4
  mat[,1][selr2]<-mat[,1][selr2]/ 12.92
  mat[,2][selg1]<-((mat[,2][selg1] + 0.055 ) / 1.055 ) ^ 2.4
  mat[,2][selg2]<-mat[,2][selg2]/ 12.92
  mat[,3][selb1]<-((mat[,3][selb1] + 0.055 ) / 1.055 ) ^ 2.4
  mat[,3][selb2]<-mat[,3][selb2]/ 12.92
  
  mat<-mat*100
  X = mat[,1] * 0.4124 + mat[,2] * 0.3576 + mat[,3] * 0.1805
  Y = mat[,1] * 0.2126 + mat[,2] * 0.7152 + mat[,3] * 0.0722
  Z = mat[,1] * 0.0193 + mat[,2] * 0.1192 + mat[,3] * 0.9505
  out<-cbind(c(X),c(Y),c(Z))
  return(out)
}
XYZ2Lab<-function(mat){
  # mat<-RGB2XYZ(image)
  # dim(mat)
  selr1<-which(mat[,1]>0.008856 ,arr.ind = T)
  selg1<-which(mat[,2]>0.008856 ,arr.ind = T)
  selb1<-which(mat[,3]>0.008856 ,arr.ind = T)
  selr2<-which(mat[,1]<=0.008856 ,arr.ind = T)
  selg2<-which(mat[,2]<=0.008856 ,arr.ind = T)
  selb2<-which(mat[,3]<=0.008856 ,arr.ind = T)
  mat[,1][selr1]<-mat[selr1,1] ^ c(1/3)
  mat[,1][selr2]<-(mat[selr2,1]*7.787) + ( 16 / 116 )
  mat[,2][selg1]<-mat[selg1,2] ^ c(1/3)
  mat[,2][selg2]<-(mat[selg2,2]*7.787) + ( 16 / 116 )
  mat[,3][selb1]<-mat[selb1,3] ^ c(1/3)
  mat[,3][selb2]<-(mat[selb2,3]*7.787) + ( 16 / 116 )
  
  CIEL<-( 116 * mat[,2] ) - 16
  CIEa<-500 * ( mat[,1] - mat[,2] )
  CIEb<-200 * ( mat[,2] - mat[,3] )
  out<-cbind(c(CIEL),c(CIEa),c(CIEb))
  return(out)
}
allcols3<-c(apply(cols,1,function(x) rgb(red=x[1],green=x[2],blue=x[3])))
HSV<-XYZ2Lab(RGB2XYZ(cols))
H<-HSV[,1];S<-HSV[,2];V<-HSV[,3];B<-unlist(apply(cols,1,rgbtoB));G<-rowMeans(cols)
allcols4<-cbind(H,S,V,B,G)
selmat2<-which(apply(as.matrix(allcols4),2,sd)>0)
pcai<-(prcomp(cbind(cols,allcols4[,selmat2]),scale. = scale)$x)
kmpls<-kmeans(pcai,centers = k,iter.max = 1000)
paletteslist<-list()
paletteslistg<-list()
for(i in 1:k){
selmat3<-which(apply(as.matrix(allcols4[which(kmpls$cluster==i),selmat2]),2,sd)>0)
plsii<-pls2B(prcomp(cols[which(kmpls$cluster==i),],scale. = scale)$x,prcomp(allcols4[which(kmpls$cluster==i),selmat2[selmat3]],scale. = scale)$x,rounds = 0)
selk<-which(kmpls$cluster==i)
paletteslist[[i]]<-allcols3[selk[order(plsii$Xscores[,plsaxis])[seq(1,length(plsii$Xscores[,plsaxis]),length.out=lcols)]]]
ciccio<-cbind(rowMeans(t(col2rgb(paletteslist[[i]])/255)),rowMeans(t(col2rgb(paletteslist[[i]])/255)),rowMeans(t(col2rgb(paletteslist[[i]])/255)))
ciccia<-c(apply(ciccio,1,function(x) rgb(red=x[1],green=x[2],blue=x[3])))
paletteslistg[[i]]<-ciccia
}

plot(NA,xlim = c(0, dim(imageR)[2]),
     ylim = c(0, dim(imageR)[1]), asp=1,axes = F,xaxt='n',yaxt='n',ann=FALSE)
rasterImage(imageR, dim(imageR)[2], dim(imageR)[1],0, 0)

plot(NA,xlim=c(0,lcols+2),ylim=c(0,k+2),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",main="")
for(i in 1:k){
  cols<-rep(1,length(paletteslist[[i]]))
  cols[which(colMeans(col2rgb(paletteslist[[i]])/255)<0.5)]<-"white"
  points(rep(i,lcols),col=paletteslist[[i]],pch=15,cex=cex)  
  text(rep(i,lcols),labels=paletteslist[[i]],cex=cext,col=cols)
}

return(paletteslist)
}
