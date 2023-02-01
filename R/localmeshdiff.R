#' localmeshdiff
#' Calculate and Visualize local differences between two meshes
#' @param mesh1 reference mesh: object of class "mesh3d"
#' @param mesh2 target mesh: object of class "mesh3d"
#' @param ploton numeric: define which mesh will be used to visualize local differences
#' @param diffarea formula: define how calculating differences in area. area_shape1 refers to mesh1, area_shape2 refers to mesh2
#' @param paltot character vector: specify the colors which are used to create a color palette
#' @param from numeric: minimum distance to be colorised
#' @param to numeric: maximum distance to be colorised
#' @param n.int numeric: determines break points for color palette
#' @param out.rem logical: if TRUE outliers will be removed
#' @param fact numeric: factor k of the interquartile range 
#' @param visual numeric: if equals to 1 the mesh is plotted without a wireframe, if set on 2 a wireframe is added  
#' @param scale01 logical: if TRUE the vector of distances is scaled from 0 to 1
#' @param colwire character: color of the wireframe
#' @return vect numeric vector containing local differeces in area between the reference and the target mesh
#' @author Antonio Profico, Costantino Buzi, Silvia Castiglione, Marina Melchionna, Paolo Piras, Pasquale Raia, Alessio Veneziano
#' @references Melchionna, M., Profico, A., Castiglione, S., Sansalone, G., Serio, C., Mondanaro, A., ... & Manzi, G. (2020). 
#' From smart apes to human brain boxes. A uniquely derived brain shape in late hominins clade. Frontiers in Earth Science, 8, 273.
#' @examples
#' \dontrun{
#' library(Arothron)
#' library(rgl)
#' data("primendoR")
#' neaset<-primendoR$sets[,,11]
#' sapset<-primendoR$sets[,,14]
#' #defining a mesh for the neanderthal right hemisphere
#' neasur<-list("vb"=t(cbind(neaset,1)),"it"=primendoR$sur$it)
#' class(neasur)<-"mesh3d"
#' #defining a mesh for the modern human right hemisphere
#' sapsur<-list("vb"=t(cbind(sapset,1)),"it"=primendoR$sur$it)
#' class(neasur)<-"mesh3d"
#' layout3d(t(c(1,2)),sharedMouse = TRUE)
#' localmeshdiff(sapsur,neasur,1,scale01 = TRUE,
#' paltot=c("darkred","red","orange","white","lightblue","blue","darkblue"))
#' next3d()
#' localmeshdiff(neasur,sapsur,1,scale01 = TRUE,
#' paltot=c("darkred","red","orange","white","lightblue","blue","darkblue"))
#' }
#' @export

localmeshdiff<-function (mesh1, mesh2,ploton=1, diffarea=((area_shape1 - area_shape2)/area_shape2)*100, paltot = rainbow(200), from = NULL, 
to = NULL, n.int = 200, out.rem = TRUE, fact = 1.5, visual = 1, 
scale01 = TRUE, colwire = "pink") 
{
  range01 <- function(x) {
    (x - min(x))/(max(x) - min(x))
  }
  area_shape1 <- vcgArea(mesh1, perface = T)$pertriangle
  area_shape2 <- vcgArea(mesh2, perface = T)$pertriangle
  
  diff_areas <- diffarea
  
  sel <- which(is.na(diff_areas))
  if (length(sel) > 0) {
    mesh1$it <- mesh1$it[, -sel]
    mesh2$it <- mesh2$it[, -sel]
    mesh1 <- rmUnrefVertex(mesh1)
    mesh2 <- rmUnrefVertex(mesh2)
    area_shape1 <- vcgArea(mesh1, perface = T)$pertriangle
    area_shape2 <- vcgArea(mesh2, perface = T)$pertriangle
    diff_areas <- (area_shape1 - area_shape2)/area_shape1
  }
  if (out.rem == TRUE) {
    x = diff_areas
    qq <- quantile(x, c(1, 3)/4, names = FALSE)
    r <- diff(qq) * fact
    tst <- x < qq[1] - r | x > qq[2] + r
    tstp <- qq[2] + r
    tstn <- qq[1] - r
    diff_areas[x > tstp] <- tstp
    diff_areas[x < tstn] <- tstn
  }
  else {
    diff_areas = diff_areas
  }
  if (scale01 == TRUE) {
    diff_areas <- range01(diff_areas)
  }
  cat("the range of diff_areas is ", range(diff_areas), sep = "\n")
  if (is.null(to) == TRUE) {
    to <- max(diff_areas) * 1.01
  }
  if (is.null(from) == TRUE) {
    from <- min(diff_areas) * 1.01
  }
  selfromto <- which(diff_areas < to & diff_areas >= from)
  diff_areas_fromto <- diff_areas[selfromto]
  if (ploton == 1) {
    meshfromto <- mesh1
    meshwhite <- mesh1
  }
  if (ploton == 2) {
    meshfromto <- mesh2
    meshwhite <- mesh2
  }
  meshfromto$it <- meshfromto$it[, selfromto]
  meshwhite$it <- meshwhite$it[, -selfromto]
  colmap_tot <- colorRampPalette(paltot)
  breaks_tot <- cut(c(from, diff_areas_fromto, to), n.int)
  cols_tot <- colmap_tot(n.int)[breaks_tot]
  cols_tot <- cols_tot[-c(1, length(cols_tot))]
  plot(density(c(from, diff_areas, to)), main = "", xlab = "", 
       ylab = "")
  abline(v = seq(from, to, length.out = n.int), col = colmap_tot(n.int), 
         lwd = 5)
  points(density(diff_areas), type = "l", lwd = 2)
  if (visual == 1) {
    triangles3d(t(meshfromto$vb[, meshfromto$it]), col = rep(cols_tot, 
                                                             each = 3), alpha = 1, lit = T, specular = "black")
    triangles3d(t(meshwhite$vb[, meshwhite$it]), col = "grey", 
                alpha = 1, lit = T, specular = "black")
  }
  if (visual == 2) {
    triangles3d(t(meshfromto$vb[, meshfromto$it]), col = rep(cols_tot, 
                                                             each = 3), alpha = 1, lit = F, specular = "black")
    triangles3d(t(meshwhite$vb[, meshwhite$it]), col = "grey", 
                alpha = 1, lit = F, specular = "black")
    wire3d(meshfromto, col = colwire, lit = F, lwd = 2)
  }
  out <- list(vect = diff_areas)
  return(out)
}