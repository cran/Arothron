#' listtoarray
#' convert a list into an array 
#' @param mylist a list
#' @return a kx3xn array with landmark coordinates
#' @author Antonio Profico, Costantino Buzi, Marina Melchionna, Paolo Piras, Pasquale Raia, Alessio Veneziano
#' @export

listtoarray<-function(mylist){
  final<-NULL
  for(i in 1:length(mylist)){
    temp<-array(mylist[[i]],dim=c(nrow(mylist[[i]]),ncol(mylist[[i]]),1))
    final<-abind::abind(final,temp)
  }
  return(final)
}