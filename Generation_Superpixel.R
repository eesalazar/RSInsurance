#GENERATION OF SUPERPIXELS
#This script supports the article entitled "Managing high-resolution, high-frequency satellite data for enhanced analyses: a machine learning superpixel-based approach"

#REQUIRED PACKAGES
library(raster)
library(maptools)
library(sf)
library(dplyr)
library(rgeoda)

#REQUIRED FUNCTIONS
contar<-function (x, ..., FUN = sum) 
{
  if (missing(...)) 
    x[] <- FUN(x)
  else {
    g <- interaction(...)
    split(x, g) <- lapply(split(x, g), FUN)
  }
  x
}

df<-readRDS("C:/PHD_REGIONES/VALIDACION_SUPERPIXEL/SUPERPIXEL1/Serie_Superpixel1_250.rds")

#INPUT INFORMATION
df<- #Input dataframe. First 2 columns--> X-Y coordinates and rest of the columns->Time series
PCA_components<- #Desired components in the PCA
kmeans_clusters<- #Desired clusters in the k-means
minimum_size<- #Minimum size of superpixels
Nsuperpixel<- #Number of superpixels if SKATER is used
CPU<- #Number of threads for parallel computation
  
#GENERATE PCA 
serie<-df[,3:ncol(df)]
coord<-df[,1:2]
respca<-prcomp(df,rank=PCA_components)

#NON-SPATIAL CLUSTERING
PCA<-as.data.frame(respca$x)
cluster_kmeans<-kmeans(PCA,kmeans_clusters)
cluster_pixel<-as.data.frame(cluster_kmeans$cluster)
df_cluster<-cbind(coord,PCA,cluster_pixel)
columns<-colnames(df_cluster)[1:(ncol(df_cluster)-1)]
columns<-append(columns,"cluster")
colnames(df_cluster)<-columns
coord_cluster<-cbind(coord,cluster_pixel)
coord_pca<-cbind(coord,PCA)

#IDENTIFYING ISOLATED PIXELS & JOINING PIXELS THAT BELONG TO THE SAME CLUSTER
df_group1<-df_cluster[df_cluster$cluster==1,] 
r1<-rasterFromXYZ(df_group1[,1:3], crs="+proj=longlat +datum=WGS84 +no_defs")  #Esto sirve para hacer el clump
r1_clump<-clump(r1,direction=8) #Direction=8->Queens case (Diagonal pixels considered) // Direction=4-> Rooks case 
df_r1_clump<-as.data.frame(r1_clump,xy=TRUE)
df_r1_clump<-df_r1_clump[rowSums(is.na(df_r1_clump))<1,]
df_r1_clump<-cbind(df_r1_clump,dplyr::select(df_group1,-(x:y)))
for (h in 1:(ncol(df_r1_clump)-3)){#Assign the mean value of the feature
  df_r1_clump[,paste0("V",h,"m")]<-ave(df_r1_clump[,3+h],df_r1_clump$clumps)
}
df_r1_clump$number<-1
df_r1_clump$px<-contar(df_r1_clump$number,df_r1_clump$clumps) #Count number of pixels of each island
df_r1_med<-dplyr::select(df_r1_clump,-(PC1:PC200))
df_r1_med<-dplyr::select(df_r1_med,-(clumps))
df_r1_med<-dplyr::select(df_r1_med,-(number))
raster_cluster1<-rasterFromXYZ(df_r1_med, crs="+proj=longlat +datum=WGS84 +no_defs")
pol_1<-rasterToPolygons(raster_cluster1)
columnas<-colnames(df_r1_med)[3:ncol(df_r1_med)]
pol<-pol_1
for (i in 2:max(df_cluster["cluster"])){#Repeat for each cluster generated in kmeans
  df_group1<-df_cluster[df_cluster$cluster==i,] 
  r1<-rasterFromXYZ(df_group1[,1:3], crs="+proj=longlat +datum=WGS84 +no_defs")  #Esto sirve para hacer el clump
  r1_clump<-clump(r1,direction=8) #Direction=8->Queens case (Diagonal pixels considered) // Direction=4-> Rooks case 
  df_r1_clump<-as.data.frame(r1_clump,xy=TRUE)
  df_r1_clump<-df_r1_clump[rowSums(is.na(df_r1_clump))<1,]
  df_r1_clump<-cbind(df_r1_clump,dplyr::select(df_group1,-(x:y)))
  for (h in 1:(ncol(df_r1_clump)-3)){#Assign the mean value of the feature
    df_r1_clump[,paste0("V",h,"m")]<-ave(df_r1_clump[,3+h],df_r1_clump$clumps)
  }
  df_r1_clump$number<-1
  df_r1_clump$px<-contar(df_r1_clump$number,df_r1_clump$clumps) #Count number of pixels of each island
  df_r1_med<-dplyr::select(df_r1_clump,-(PC1:PC200))
  df_r1_med<-dplyr::select(df_r1_med,-(clumps))
  df_r1_med<-dplyr::select(df_r1_med,-(number))
  raster_cluster1<-rasterFromXYZ(df_r1_med, crs="+proj=longlat +datum=WGS84 +no_defs")
  pol_1<-rasterToPolygons(raster_cluster1)
  columnas<-colnames(df_r1_med)[3:ncol(df_r1_med)]
  pol<-pol_1
  pol<-rbind(pol,pol_1)
}

#SPATIAL CONSTRAINED CLUSTERING
guerry2<-pol
queen_w2<-queen_weights(guerry2)
guerry_geometry<-guerry2["geometry"]
data<-dplyr::select(as.data.frame(guerry2),(PC1:PC200))
px<-as.numeric(unlist(as.data.frame(guerry2)["px"]))
bound_variable<-cbind(guerry_geometry,px)
      #MAX-P
result <- rgeoda::maxp_greedy(w=queen_w2,df=data, bound_variable, min_bound=minimum_size,cpu_threads=CPU)
      #SKATER 
result <- rgeoda::skater(k=Nsuperpixel, w=queen_w2,df=data, bound_variable, min_bound=minimum_size,cpu_threads=CPU)

group<-as.data.frame(result$Clusters)
group$cluster<-result$Clusters
group<-dplyr::select(group,(cluster))
df_cluster<-cbind(guerry_geometry,group)  #This shows each pixel which superpixel it belongs to
plot(df_cluster)

#CONVERT TO RASTER (Optional)

df<-readRDS("C:/PHD_REGIONES/VALIDACION_SUPERPIXEL/SUPERPIXEL1/Serie_Superpixel1_250.rds")
coord<-df[,1:2]
serie<-df[,3:946]
poligono<-as_Spatial(df_cluster)
raster_base<- coord
raster_base$num<-0
raster_base<-rasterFromXYZ(raster_base, crs="+proj=longlat +datum=WGS84 +no_defs")  #Esto sirve para hacer el clump
raster_final<-rasterize(x=poligono,y=raster_base, field="cluster") #raster_base solo aporta el molde
saveRDS(raster_final,paste0("C:/PHD_REGIONES/VALIDACION_SUPERPIXEL/RDS_Raster_SuperPixel_1000_SKATER.rds"))
writeRaster(raster_final,paste0("C:/PHD_REGIONES/VALIDACION_SUPERPIXEL/Raster_SuperPixel_1000_SKATER.tif"))

#ASSIGN TO EACH PIXEL THE AVERAGE OF THE SUPERPIXEL
poligon<-as_Spatial(df_cluster)
raster_base<- coord
raster_base$num<-0
raster_base<-rasterFromXYZ(raster_base, crs="+proj=longlat +datum=WGS84 +no_defs")  
raster_final<-rasterize(x=poligon,y=raster_base, field="cluster") 
clusters<-raster_final 
clusters<-as.data.frame(clusters,xy=TRUE)
clusters<-clusters[rowSums(is.na(clusters))<1,]
clusters_serie<-cbind(clusters,df)
clusters_serie[,4:5]<-NULL
clusters_mean<-clusters_serie[,1:3]
for (h in 1:(ncol(clusters_serie)-3)){
  clusters_mean[,paste0("X",h)]<-ave(clusters_serie[,3+h],clusters_serie$layer)
}
df_group1<-clusters_mean[clusters_mean$layer==1,] 
df_group1<-dplyr::select(df_group1,-(layer))
raster_cluster1<-rasterFromXYZ(df_group1, crs="+proj=longlat +datum=WGS84 +no_defs")
pol_1<-rasterToPolygons(raster_cluster1)
columnas<-colnames(df_group1)[3:ncol(df_group1)]
pol_1<-raster::aggregate(pol_1,by=columnas)
pol<-pol_1
for (i in 2:max(clusters_mean["layer"])){
  df_group1<-clusters_mean[clusters_mean$layer==i,] 
  df_group1<-dplyr::select(df_group1,-(layer))
  raster_cluster1<-rasterFromXYZ(df_group1,crs="+proj=longlat +datum=WGS84 +no_defs")
  pol_1<-rasterToPolygons(raster_cluster1)
  columnas<-colnames(df_group1)[3:ncol(df_group1)]
  pol_1<-raster::aggregate(pol_1,by=columnas)
  pol<-rbind(pol,pol_1)
}

FINAL_RESULT<-pol
