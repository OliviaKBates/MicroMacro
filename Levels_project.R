## R code for "Soil microclimates predict smaller niche shifts than macroclimates in introduced ant species"

####	####	####	####
# this file produces:  
#an antmaps_occ of just species with  10 occurence points in each range - saved as '~/Desktop/Antmaps_data/antmapsData.RData'
# Microclimate datasets cropped and stacked on eachother saved in the nas '"/Volumes/FBM/DEE/GROUPS/gr-Bertelsmeier/PRIVE/obates/Data/Microclimate_data.RData'
rm(list = ls())

# Load libraries 
library(ade4); library(sf)
library(sp); library(rgdal)
library(rgeos); library(dismo)
library(ecospat); library(ggplot2)
library(dplyr); library(tidyr)
library(reshape2); library(raster)
library(rworldmap); library(RStoolbox)
library(viridis); library(ggridges)
library(geodata); library(maptools)
library(ggnewscale); library(geosphere)
library(outliers); library(class)
library(scales); library(ggpubr)
library(grDevices); library(extrafont)
library(dunn.test); library(aqp)
# Additional configuration
sf_use_s2(FALSE)


#### Load data
lcl <- ""
setwd(lcl)
load('MicroMacro.RData')
#load climate datasets 
#save(species_data, Species, ordination_results, hypervolume_results,ordination_restrict, hypervolume_restrict, nesting, landcover, w_rworldmap, sample_size, file='MicroMacro.RData')

## Step 1 - extracts raw climate data for each occurence point. For this use 
      #species_data <--  Species for which there were > 20 occurence points in both native and introduced range, after thinning to 1km resolution 
      ## to have it all in one dataframe if you need use code: all_ants <-  bind_rows(species_data,  .id = "column_label")
      #Species  <-- list of species 
      #r  <--  Chelsa dataset at 1km 
      # SBio <-- SoilTemp dataset at 1km 
#pca_SoilTemp <-- pca of soiltemp variables, created via: pca_soiltemp<- rasterPCA(SBio ,spca=TRUE, nComp=4, maskCheck=F) # this takes a long time 
#pca_chelsa <- pca of the chelsa variables, created via:  pca_xhelsa<- rasterPCA(r ,spca=TRUE, nComp=4, maskCheck=F) # this takes a long time 


## Step 2 - niche shift analysis, produces 
      #ordination_results  <-- statistics of niche shifts using ordination approach 
          #for each species, D overlap, expansion, stability, unfilling results are presented (distinguished between datasets by Chelsa_ and soil_)
          # bio1...bio11 display the correlation to each bio variable for Chelsa Dataset during between pca 
          # soilbio1...soilbio11 display the correlation to each bio variable for the SoilTemp dataset during between pca 
      #hypervolume_results  <-- statistics of niche shifts using hypervolume approach 
              # statistics 
                  #'Determinant ratio (Determinant_ratio), Mahalanobis distance (Mahalanobis_distance ), Bhattacharyya distance (Bhattacharyya_distance), centroid distance (dis_centroid), Jaccard similarity (jaccard ), niche volume of native range (volume_nat), niche volume of invasive range  (volume_inv)
      #ordination_restrict #for restrited continent analysis (see SM)
      #hypervolume_restrict  # for restrited continent analysis (see SM)

### Step 3 - analysis stages
      # figure 4 
          # tab contains the correlation of occurence points for each bio variable for each ('Species') and ('Random')
      #nesting # nesting catefories for species 
          ## categorised into: (Categories)
                  # C == aboreal cavities and all structure present in living trees/plants/ dead trees/ 
                  # B == aboreal builds 
                  # F == in fungus 
                  # L == ground leaf litter 
                  # M == under moss/wood
                  # R == rocks/crevicies 
                  # S = in the soil/under rocks/objects on ground 
                  # W == dead wood, stumps on floor 
          # Cat_sim  cattegorieses these further into Aboreal (high), ground (low), or generalist (wide)
      #landcover ## landcover data ‘HYBMAP’, a Global hybrid land-cover map of categorised land-cover data from the year 2013 (Zhu et al., 2023) 
      #species_data #used for figures 4 and 5 (see step 1)
      #w_rworldmap # rworldmap
      #sample_size #sample size for each species, split by native sample size (nativeSS), introduced range sample size (invSS), and species name (Species)


##### Below code produces variable 'species_data' 
### Extracting out the values for each of the occurence points 
gridnumbers <- NULL
niche_expansion_1 <- NULL
sample_size_1 <- NULL 
## creating a list with datasets for each species 
species_data <- list()
SBio <- raster::stack(SBio)
r <- raster::stack(r)
####  Niche shift analysis 
###Micrcoclimate levels - climate data levels PART 2
#comparing between different levels of climate data to observe niche shifts. 
####	####	####	####
##### ##### ##### ##### ##### 
# functions 
##### ##### ##### ##### ##### 

##initialise the graph code for every run
"map_overlap" <- function(double=F, y, y2=NULL, status, species_name, graph_path, axis_label=NULL){
  # argument double: to plot two datasets in one go (i.e. Macro and Microclimate )
  if (double==F){
    mypath3 <- file.path(graph_path) #path to save file 
    png(mypath3 )
    y_E <- y[status == "E"] #specify non-native 
    y_N <- y[status == "N"]#specify native
    dens_E <- density(y_E)  #density of occurence
    dens_N <- density(y_N)
    xlim <- range(dens_E$x, dens_N$x) #x axis range
    ylim <- range(0, dens_E$y, dens_N$y)
    ylim[2] <- ylim[2] + 0.1
    hist(y_E, proba = T,  xlim = xlim, ylim = ylim, col= 'white', density = 10, angle = 135, main = species_name, xlab = "Between-Class PCA 1", ylab = "Density", cex.lab=1.5, cex.axis=1.5)
    hist(y_N, proba = T, add = T, col = 'white', density = 10, angle = 45)
    polygon(dens_E, density = -1, col = col_E)
    polygon(dens_N, density = -1, col = col_N)
    mat <- cbind(dens_E$x, dens_E$y)
    mat <- rbind(mat, mat[1,])
    pol_E <- st_polygon(list(mat))
    mat <- cbind(dens_N$x, dens_N$y)
    mat <- rbind(mat, mat[1,])
    pol_N <- st_polygon(list(mat))
    pol_inter <- st_intersection(st_buffer(pol_E,0),  st_buffer(pol_N,0))
    plot(pol_inter, add = T, col = "grey", ylim=c(0,1))
    rug(y_N, side = 1, line = -0.15, col = col_N, tck = 0.03, lwd=2)
    rug(y_E, side = 1, line = -0.15, col = col_E, tck = 0.03, lwd=2)
    dev.off()
  }else { #this is if double=T , i.e. present two datasets on one file 
    mypath3 <- file.path(graph_path) #path to save file 
    png(  mypath3 )
    #world
    col_E <- rgb(1,0,0,0.4) #non-native colours
    col_N <- rgb(0,0,1,0.35)
    y_E <- y[status == "E"]
    y_N <- y[status == "N"]
    dens_E <- density(y_E) 
    dens_N <- density(y_N)
    #Soil
    col2_E <- rgb(1,0.6, 0,0.5) #non-native colours
    col2_N <- rgb(0.2,0.6,0.6,0.4)
    y2_E <- y2[status == "E"]
    y2_N <- y2[status == "N"]
    dens2_E <- density(y2_E) 
    dens2_N <- density(y2_N)
    #limits
    xlim <- range(c(dens_E$x, dens_N$x,   dens2_E$x,     dens2_N$x) )#x axis range
    ylim <- range(0, dens_E$y, dens_N$y,  dens2_E$y, dens2_N$y)
    ylim[2] <- ylim[2] + 0.1
    hist(y_E, proba = T,  xlim = xlim, ylim = ylim, col= 'white', density = 10, angle = 135, main = species_name, xlab = axis_label, ylab = "Density", cex.lab=1.5, cex.axis=1.5)
    hist(y_N, proba = T, add = T, col = 'white', density = 10, angle = 45)
    hist(y2_E, proba = T, add = T, col = 'white', density = 10, angle = 45)
    hist(y2_N, proba = T, add = T, col = 'white', density = 10, angle = 45)
    polygon(dens_E, density = -1, col = col_E, border='darkblue')
    polygon(dens_N, density = -1, col = col_N, border= 'darkblue')
    polygon(dens2_E, density = -1, col = col2_E, border= 'darkgreen')
    polygon(dens2_N, density = -1, col = col2_N, border='darkgreen')
    #plot intersection and rugs world
    mat <- cbind(dens_E$x, dens_E$y)
    mat <- rbind(mat, mat[1,])
    pol_E <- st_polygon(list(mat))
    mat <- cbind(dens_N$x, dens_N$y)
    mat <- rbind(mat, mat[1,])
    pol_N <- st_polygon(list(mat))
    pol_inter <- st_intersection(st_buffer(pol_E,0),  st_buffer(pol_N,0))
    # plot(pol_inter, add = T, col = "grey", ylim=c(0,1))
    rug(y_N, side = 1, line = -0.15, col = col_N, tck = 0.02, lwd=0.5)
    rug(y_E, side = 1, line = -0.15, col = col_E, tck = 0.02, lwd=0.5)
    #plot intersection and rugs soil
    mat <- cbind(dens2_E$x, dens2_E$y)
    mat <- rbind(mat, mat[1,])
    pol_E <- st_polygon(list(mat))
    mat <- cbind(dens2_N$x, dens2_N$y)
    mat <- rbind(mat, mat[1,])
    pol_N <- st_polygon(list(mat))
    pol_inter <- st_intersection(st_buffer(pol_E,0),  st_buffer(pol_N,0))
    # plot(pol_inter, add = T, col = "grey", ylim=c(0,1))
    rug(y2_N, side = 1, line = -0.15, col = col2_N, tck = 0.03, lwd=0.5)
    rug(y2_E, side = 1, line = -0.15, col = col2_E, tck = 0.03, lwd=0.5) 
    legend('topright', pch=rep(19, 4), col=c(col_N, col_E, col2_N, col2_E ), 
           legend=c('Native Chelsa', 'Non-native Chelsa', 'Native Soiltemp', 'Non-native Soiltemp'), 
           pt.cex=2, bty="n")
    dev.off()
  }
}
#### Extracting occurences and calculating Overlap Niche expansion and Equivilency Test
"ordination_shift" <- function(species = species, data =species_data,  col_names= NULL, nf=5, graph_path = NULL  ) {
  par(mfrow=c(1,1))
  par(mar=c(5, 5, 5, 5))
  Native_exotic_df <- species_data[[species]]
  status <- as.factor(Native_exotic_df$status) #whether point is non-native or native 
  ### PCA and between class analysis
  pca.env <- dudi.pca(Native_exotic_df[,col_names], scannf = F, nf =nf) ##so can first make a pca of all, look at eigen values ect. 
  status <- as.factor(status)
  bet1 <- bca(pca.env, status, scan = FALSE, nf = nf)    ## perform between-class analysis : highlihgts the differences between groups 
  #bet1$ratio #between-class ratio
  ## explore correlaitons to each bio variable to axis
  correlation_to_pca.env <- as.data.frame (bet1$as)#see the correlation with the normal pca for the axis. 
  correlation_to_pca <- t(correlation_to_pca.env) #this is a data frame of the contribution of each origional pca axis to the between-pca axis
  rownames(correlation_to_pca) <-  c(species)
  correlation_to_bio_variables <- as.data.frame (bet1$co) #to see the contribution of each worldclim BIO variable
  correlation_to_bio_variables <- t(correlation_to_bio_variables) 
  rownames(correlation_to_bio_variables) <-  c(species)
  
  ## save interesting values  
  between_class <- as.data.frame(bet1$ls) #data frame of the first axis of between class analysis
  data <-  cbind(status, between_class) #binding the two together
  #the scores for the whole environment
  scores.globclim <- data$CS1
  #scores for the presences 
  scores.sp.ex <- data$CS1[which(data$status=='E')] #only non-native points
  scores.sp.nat <-data$CS1[which(data$status=='N')]  # only native points 
  #scores for the whole study region, same as presences as no considering whole study area 
  scores.clim.ex <- data$CS1[which(data$status=='E')] # only non-native points
  scores.clim.nat <-data$CS1[which(data$status=='N')] #only native points 
  
  ### Evaluate niche shifts
  ## Grid the whole environment 
  grid.clim.nat <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1= scores.sp.nat, sp=scores.sp.nat, R=100, th.sp=0)
  #glob= all background dataframe, glob1= of the native/invasive envinoment range of the species (for each one), envirnomental variable for occurances only, 
  grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.sp.ex , sp=scores.sp.ex, R=100, th.sp=0)
  ## D Overlap
  D.overlap <- ecospat.niche.overlap (grid.clim.nat, grid.clim.inv, cor=F)$D #niche overlap. 
  # Niche expansion : delimiting niche categories and quantifying niche dynamics in analogue climates
  niche.dyn.whole <- ecospat.niche.dyn.index (grid.clim.nat, grid.clim.inv, intersection = NA) #if intersection= NA means no intersection of analogous envirnoment
  niche.dyn.whole$dynamic.index.w 
  ## save results
  row <- cbind.data.frame(D.overlap,  niche.dyn.whole$dynamic.index.w[[1]],niche.dyn.whole$dynamic.index.w[[2]], niche.dyn.whole$dynamic.index.w[[3]])
  names(row) <- c('D.overlap',  'expansion', 'stability', 'Unfilling')
  rownames(row) <- c(species)
  niche_expansion_1 <- cbind ( row,  correlation_to_pca,  correlation_to_bio_variables) #add the correlations of bio and pca axis to each row 
  return( niche_expansion_1)
  ## plot results
  map_overlap(bet1$ls[,1], status, species)
  #dev.copy(png, file=mypath3,  width = 7.7, height = 6.5, units='in', res=500)
  dev.off() #no more added to the file
  
}

'coords2continent' <-  function(points){  
  # The single argument to this function, points, is a data.frame in which:
  #   - column 1 contains the longitude in degrees
  #   - column 2 contains the latitude in degrees
  countriesSP <- getMap(resolution='low')
  #countriesSP <- getMap(resolution='high') #you could use high res map from rworldxtra if you were concerned about detail
  # converting points to a SpatialPoints object
  pointsSP = SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))  
  # use 'over' to get indices of the Polygons object containing each point 
  indices = over(pointsSP, countriesSP)
  #indices$continent   # returns the continent (6 continent model)
  indices$REGION   # returns the continent (7 continent model)
}

### per axis function 
"perAxis_shift" <- function(species = species, data =species_data,  vectors=vectors, nf=5, graph_path = NULL , dataset_name=NULL  ) {
  Native_exotic_df <- species_data[[species]]
  Native_exotic_df <-   Native_exotic_df %>% drop_na (SBIO3_0_5cm_Isothermality) 
  status <-  as.factor(Native_exotic_df$status)
  niche_expansion_1  <- data.frame('colname'=c())
  niche_expansion_2  <- data.frame('colname'=c())
  for (v in 1:nrow(vectors)){
    i = vectors[v, 1]
    ### PCA and between class analysis
    Native_exotic_df[,i]
    scores.globclim <- Native_exotic_df[,i]
    #scores for the presences 
    scores.sp.ex <- Native_exotic_df[,i][which(Native_exotic_df$status=='E')]
    scores.sp.nat <-Native_exotic_df[,i][which(Native_exotic_df$status=='N')] 
    #scores for the whole study region, same as presences as no considering whole study area 
    scores.clim.ex <- Native_exotic_df[,i][which(Native_exotic_df$status=='E')]
    scores.clim.nat <-Native_exotic_df[,i][which(Native_exotic_df$status=='N')] 
    ### Evaluate niche shifts
    ## Grid the whole environment 
    grid.clim.nat <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1= scores.sp.nat, sp=scores.sp.nat, R=100, th.sp=0)
    #glob= all background dataframe, glob1= of the native/invasive envinoment range of the species (for each one), envirnomental variable for occurances only, 
    grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.sp.ex , sp=scores.sp.ex, R=100, th.sp=0)
    
    ## D Overlap
    D.overlap <- ecospat.niche.overlap (grid.clim.nat, grid.clim.inv, cor=F)$D #niche overlap. 
    # Niche expansion : delimiting niche categories and quantifying niche dynamics in analogue climates
    niche.dyn.whole <- ecospat.niche.dyn.index (grid.clim.nat, grid.clim.inv, intersection = NA) #if intersection= NA means no intersection of analogous envirnoment
    niche.dyn.whole$dynamic.index.w 
    ## save results
    row <- cbind.data.frame(D.overlap,  niche.dyn.whole$dynamic.index.w[[1]],niche.dyn.whole$dynamic.index.w[[2]], niche.dyn.whole$dynamic.index.w[[3]])
    names(row) <- c('D.overlap',  'expansion', 'stability', 'Unfilling')
    axis <- colnames(Native_exotic_df[i])
    colnames(row) <- paste0(colnames(row), "_", dataset_name[1])
    row$axis <- axis
    niche_expansion_1 <- rbind(niche_expansion_1, row)
    world_dat <-  Native_exotic_df[,i]
    j <-  vectors[v, 2]
    ### PCA and between class analysis
    Native_exotic_df[,j]
    scores.globclim <- Native_exotic_df[,j]
    #scores for the presences 
    scores.sp.ex <- Native_exotic_df[,j][which(Native_exotic_df$status=='E')]
    scores.sp.nat <-Native_exotic_df[,j][which(Native_exotic_df$status=='N')] 
    #scores for the whole study region, same as presences as no considering whole study area 
    scores.clim.ex <- Native_exotic_df[,j][which(Native_exotic_df$status=='E')]
    scores.clim.nat <-Native_exotic_df[,j][which(Native_exotic_df$status=='N')] 
    ### Evaluate niche shifts
    ## Grid the whole environment 
    grid.clim.nat <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1= scores.sp.nat, sp=scores.sp.nat, R=100, th.sp=0)
    #glob= all background dataframe, glob1= of the native/invasive envinoment range of the species (for each one), envirnomental variable for occurances only, 
    grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.sp.ex , sp=scores.sp.ex, R=100, th.sp=0)
    
    ## D Overlap
    D.overlap <- ecospat.niche.overlap (grid.clim.nat, grid.clim.inv, cor=F)$D #niche overlap. 
    # Niche expansion : delimiting niche categories and quantifying niche dynamics in analogue climates
    niche.dyn.whole <- ecospat.niche.dyn.index (grid.clim.nat, grid.clim.inv, intersection = NA) #if intersection= NA means no intersection of analogous envirnoment
    niche.dyn.whole$dynamic.index.w 
    
    ## save results
    row <- cbind.data.frame(D.overlap,  niche.dyn.whole$dynamic.index.w[[1]],niche.dyn.whole$dynamic.index.w[[2]], niche.dyn.whole$dynamic.index.w[[3]])
    names(row) <- c('D.overlap',  'expansion', 'stability', 'Unfilling')
    axis <- colnames(Native_exotic_df[j])
    colnames(row) <- paste0(colnames(row), "_", dataset_name[2])
    row$axis1 <- axis
    niche_expansion_2 <- rbind(niche_expansion_2, row)
    soil_dat <-  Native_exotic_df[,j]
    map_overlap(  double=T, y= world_dat,   y2=soil_dat , status, species_name=species, axis_label=axis, graph_path= paste0(graph_path, species, axis, ".jpg"))
  }
  ## plot results
  #add the correlations of bio and pca axis to each row 
  #niche_expansion_1  <-  niche_expansion_1   %>% pivot_wider(names_from = axis, values_from = c(D.overlap, expansion, stability, Unfilling)) %>% as.data.frame()
  niche_expansion_1$species <-  species 
  niche_expansion_2$species <-  species 
  niche_expansion <- cbind( niche_expansion_2,  niche_expansion_1 )
  return(     niche_expansion)
}

##### ##### ##### ##### ##### 
# Run the overlap for all species  per axis 
##### ##### ##### ##### #####

##### ##### ##### ##### ##### 
# Run the overlap for all sepcies 
##### ##### ##### ##### #####

## All continents included model 
niche_results_temp <- NULL
for (j in 1:length(Species)) {
  col_E <- rgb(1,0,0,0.2) #non-native colours
  col_N <- rgb(0,0,1,0.2)
  graph_path <- paste0("../Graphs/Density_plots/Chelsa_temp/World/", Species[j], ".jpg")
  niche_expansion_worldclim <- ordination_shift(Species[j], col_names= c(15:25), graph_path=graph_path)
  colnames(   niche_expansion_worldclim ) <- c("Chelsa_D.overlap", "Chelsa_expansion", "Chelsa_stability", "Chelsa_Unfilling", "Chelsa_Axis1", "Chelsa_Axis2" , "Chelsa_Axis3" ,"Chelsa_Axis4", "Chelsa_Axis5", 'bio1', 'bio10', 'bio11', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9' )
  col_E <- rgb(0.55,0.2,0.65,0.3) #non-native colours
  col_N <- rgb(0.2,0.6,0.6,0.3)
  graph_path <- paste0("../Graphs/Density_plots/Soiltemps/World/", "density_all", Species[j], ".jpg")
  niche_expansion_soil <- ordination_shift(Species[j], col_names= c(4:14), graph_path=graph_path)
  colnames(  niche_expansion_soil) <- c("soil_D.overlap", "soil_expansion", "soil_stability", "soil_Unfilling", "soil_Axis1",     "soil_Axis2" ,    "soil_Axis3" ,    "soil_Axis4"   ,  "soil_Axis5", 'soilbio1', 'soilbio10', 'soilbio11', 'soilbio2', 'soilbio3', 'soilbio4', 'soilbio5', 'soilbio6', 'soilbio7', 'soilbio8', 'soilbio9'  )
  #bind together and save data 
  niche_expansion_1 <- cbind( niche_expansion_worldclim ,  niche_expansion_soil)
  niche_results_temp <- rbind(niche_results_temp, niche_expansion_1)
}
## Datasets are saved in this for script part 3
#write.csv(niche_results_temp , file=paste0(lcl, "Niche_result_soil_world.csv") ) 

###### ###### ###### ###### ###### ###### ###### ###### ###### 
### Hypervolume approach 
### If you do not want to run this part, you can use directly the hypervolume_soiltemp hypervolume_chelsa objects 
###### ###### ###### ###### ###### ###### ###### ###### ###### 
pca_chelsa <- rast(pca_chelsa)
pca_soiltemp <- rast(pca_soiltemp )
## Chelsa data 
#### NON FIXED Bandwidth method - to calcluate optimum bandwidths
setwd(lcl)
bandwidth_details <- NULL 
volume <- NULL 
climate_data <- pca_chelsa
for (i in 1:length(Species)) {
  ##extract the points for each species 
  #species <- data[which(data$Species=='a_taeniatulus'), ]
  species <-species_data[[as.character(Species[i])]]
  # for restrictive continent model run   species <- species %>% filter(species$continent=='Europe'| species$continent=='North America')
  ##seperate the distirbutions by native and exotic! (and then maybe by continant) 
  lat <-   species$lat
  lon <- species$lon
  coords <- data.frame(x=lon,y=lat) #creating a data frame of the latitudes and longitudes 
  points <- SpatialPointsDataFrame(coords, data=species, proj4string = climate_data@crs) #this converts the lat and longs into spatial data
  values <- raster::extract(climate_data, points, cellnumbers=T) #Ω` extracts the bioclim data for the coodrinate points we have 
  local_df <- cbind.data.frame(coordinates(points),  values, points$status)  #binds the bioclim data 
  local_df <- local_df[complete.cases(local_df),]
  #make a hypervolume, as based on the same pca axis then it should be comparable between species 
  #estimate bandwidth for each species, then use the highest one. 
  hyper_world <- hypervolume_gaussian(  local_df [,4:7], chunk.size=500) #gaussian method is better as the box method is flat therefore has errors 
  #dataset with badwidth details 
  bandwidth_sp <- hyper_world@Parameters$kde.bandwidth # this gives us the bandwidth. 
  bandwidth_sp <- as.data.frame(bandwidth_sp)
  bandwidth_sp <- t(bandwidth_sp)
  rownames(bandwidth_sp) <- as.character(Species[i])
  bandwidth_details <- rbind(bandwidth_details, bandwidth_sp)
  #dataset with the volume details 
  volume_sp <- get_volume(hyper_world)
  volume_sp <- as.data.frame(volume_sp)
  rownames(volume_sp) <- as.character(Species[i])
  volume <- rbind(volume, volume_sp)
  rm(hyper_world)
}
######
#remove the outliers ! 
bandwidth_details <- as.data.frame(bandwidth_details)
colnames(bandwidth_details) <- c('PC1', 'PC2', 'PC3', 'PC4')
while (grubbs.test(bandwidth_details$PC1)$p.value <0.05) {
  bandwidth_details$PC1[which(bandwidth_details$PC1==outlier(bandwidth_details$PC1))] <- NA
}
pca1 <- max(na.omit(bandwidth_details$PC1)) # 2,28
while (grubbs.test(bandwidth_details$PC2)$p.value <0.05) {
  bandwidth_details$PC2[which(bandwidth_details$PC2==outlier(bandwidth_details$PC2))] <- NA
}
pca2 <- max(na.omit(bandwidth_details$PC2)) #  1.28
while (grubbs.test(bandwidth_details$PC3)$p.value <0.05) {
  bandwidth_details$PC3[which(bandwidth_details$PC3==outlier(bandwidth_details$PC3))] <- NA
}
pca3 <-max(na.omit(bandwidth_details$PC3)) # 1.32
while (grubbs.test(bandwidth_details$PC4)$p.value <0.05) {
  bandwidth_details$PC4[which(bandwidth_details$PC4==outlier(bandwidth_details$PC4))] <- NA
}
pca4 <- max(na.omit(bandwidth_details$PC4)) # 0.67
bwdths <- c(pca1, pca2, pca3, pca4)

## then re-run using fixed bandwidth of =previous section
hypervolume_chelsa <- list()
length(Species)
for (i in 1:length(Species)) {
  species <-species_data[[as.character(Species[i])]]
  # for restrictive continent model run   species <- species %>% filter(species$continent=='Europe'| species$continent=='North America')
  ##seperate the distirbutions by native and exotic! 
  lat <-   species$lat
  lon <- species$lon
  coords <- data.frame(x=lon,y=lat) #creating a data frame of the latitudes and longitudes 
  points <- SpatialPointsDataFrame(coords, data=species, proj4string = climate_data@crs) #this converts the lat and longs into spatial data
  values <- raster::extract(climate_data, points, cellnumbers=F) #extracts the bioclim data for the coodrinate points we have 
  local_df <- cbind.data.frame(coordinates(points),  values, points$status)  #binds the bioclim data 
  local_df <- local_df[complete.cases(local_df),]
  # make a hypervolume, as based on the same pca axis then it should be comparable between species 
  native <- local_df[which(local_df$`points$status`=='N'),]
  hyper_native <- hypervolume_gaussian(native[,3:6], chunk.size=500, weight=NULL, kde.bandwidth=estimate_bandwidth(local_df [,4:7]), name='native' )#gaussian method is better as the box method is flat therefore has errors
  invasive <-  local_df[which(local_df$`points$status`=='E'),]
  hyper_invasive <- hypervolume_gaussian( invasive[,3:6], chunk.size=500, kde.bandwidth=estimate_bandwidth(local_df [,4:7]), name='invasive' ) #gaussian method is better as the box method is flat therefore has errors
  # join all together 
  sp_hypervolume <- hypervolume_join(   hyper_native,  hyper_invasive, names=c('native', 'invasive') )
  hypervolume_chelsa[[Species[i]]] <-  sp_hypervolume
  rm(hyper_invasive, hyper_native)
  print(paste0(i, ' of 100'))
}
names(hypervolume_chelsa)

##  SoilTemp data hypervolume ##
bandwidth_details <- NULL 
volume <- NULL 
climate_data <- pca_soiltemp
for (i in 1:length(Species)) {
  ##extract the points for each species 
  species <-species_data[[as.character(Species[i])]]
  # for restrictive continent model run   species <- species %>% filter(species$continent=='Europe'| species$continent=='North America')
  ##seperate the distirbutions by native and exotic! 
  lat <-   species$lat
  lon <- species$lon
  coords <- data.frame(x=lon,y=lat) #creating a data frame of the latitudes and longitudes 
  points <- SpatialPointsDataFrame(coords, data=species, proj4string = climate_data@crs) #this converts the lat and longs into spatial data
  values <- raster::extract(climate_data, points, cellnumbers=T) #Ω` extracts the bioclim data for the coodrinate points we have 
  local_df <- cbind.data.frame(coordinates(points),  values, points$status)  #binds the bioclim data 
  local_df <- local_df[complete.cases(local_df),]
  #make a hypervolume, as based on the same pca axis then it should be comparable between species 
  #estimate bandwidth for each species, then use the highest one. 
  hyper_world <- hypervolume_gaussian(  local_df [,4:7], chunk.size=500) #gaussian method is better as the box method is flat therefore has errors 
  #dataset with badwidth details 
  bandwidth_sp <- hyper_world@Parameters$kde.bandwidth # this gives us the bandwidth. 
  bandwidth_sp <- as.data.frame(bandwidth_sp)
  bandwidth_sp <- t(bandwidth_sp)
  rownames(bandwidth_sp) <- as.character(Species[i])
  bandwidth_details <- rbind(bandwidth_details, bandwidth_sp)
  #dataset with the volume details 
  volume_sp <- get_volume(hyper_world)
  volume_sp <- as.data.frame(volume_sp)
  rownames(volume_sp) <- as.character(Species[i])
  volume <- rbind(volume, volume_sp)
  rm(hyper_world)
}
######
#remove the outliers ! 
bandwidth_details <- as.data.frame(bandwidth_details)
colnames(bandwidth_details) <- c('PC1', 'PC2', 'PC3', 'PC4')
while (grubbs.test(bandwidth_details$PC1)$p.value <0.05) {
  bandwidth_details$PC1[which(bandwidth_details$PC1==outlier(bandwidth_details$PC1))] <- NA
}
pca1 <- max(na.omit(bandwidth_details$PC1)) # 2,28
while (grubbs.test(bandwidth_details$PC2)$p.value <0.05) {
  bandwidth_details$PC2[which(bandwidth_details$PC2==outlier(bandwidth_details$PC2))] <- NA
}
pca2 <- max(na.omit(bandwidth_details$PC2)) #  1.28
while (grubbs.test(bandwidth_details$PC3)$p.value <0.05) {
  bandwidth_details$PC3[which(bandwidth_details$PC3==outlier(bandwidth_details$PC3))] <- NA
}
pca3 <-max(na.omit(bandwidth_details$PC3)) # 1.32
while (grubbs.test(bandwidth_details$PC4)$p.value <0.05) {
  bandwidth_details$PC4[which(bandwidth_details$PC4==outlier(bandwidth_details$PC4))] <- NA
}
pca4 <- max(na.omit(bandwidth_details$PC4)) # 0.67
bwdths <- c(pca1, pca2, pca3, pca4)

hypervolume_soiltemp<- list()
climate_data <- pca_soiltemp
for (i in 1:length(Species)) {
  species <-species_data[[as.character(Species[i])]]
  # for restrictive continent model run   species <- species %>% filter(species$continent=='Europe'| species$continent=='North America')
  ##seperate the distirbutions by native and exotic!
  lat <-   species$lat
  lon <- species$lon
  coords <- data.frame(x=lon,y=lat) #creating a data frame of the latitudes and longitudes 
  points <- SpatialPointsDataFrame(coords, data=species, proj4string = climate_data@crs) #this converts the lat and longs into spatial data
  values <- raster::extract(climate_data, points, cellnumbers=F) #extracts the bioclim data for the coodrinate points we have 
  local_df <- cbind.data.frame(coordinates(points),  values, points$status)  #binds the bioclim data 
  local_df <- local_df[complete.cases(local_df),]
  #make a hypervolume, as based on the same pca axis then it should be comparable between species 
  native <-  local_df[which(local_df$`points$status`=='N'),]
  hyper_native <- hypervolume_gaussian(native[,3:6], chunk.size=500, kde.bandwidth=estimate_bandwidth(local_df [,4:7]), name='native' )#gaussian method is better as the box method is flat therefore has errors 
  invasive <-  local_df[which(local_df$`points$status`=='E'),]
  hyper_invasive <- hypervolume_gaussian( invasive[,3:6], chunk.size=500, kde.bandwidth=estimate_bandwidth(local_df [,4:7]), name='invasive' ) #gaussian method is better as the box method is flat therefore has errors 
  # join it all together 
  sp_hypervolume <- hypervolume_join(   hyper_native,  hyper_invasive, names=c('native', 'invasive') )
  hypervolume_soiltemp[[Species[i]]] <-  sp_hypervolume
  print(paste0(i, ' of 100'))
}
names(hypervolume_soiltemp) # check its okay 

####### ###
###  calculate the overlap between the two datasets for hypervolume by assessing niche shifts for each seperatly 
####### ###
stats <- NULL
hypervolume <- hypervolume_chelsa
for (j in 1:length(Species)) {
  spe <-   hypervolume[[j]] 
  name <-  names( hypervolume[j] )
  ##first, have the sizes of the different hypervolumes 
  # get volume for each hypervolume 
  volume_sp <- get_volume(spe) %>% as.data.frame() %>% t()
  colnames(  volume_sp) <- c('volume_nat', 'volume_inv')
  ##then assess centroid distance and overlap  
  hv_set <- hypervolume_set(spe[[1]],spe[[2]], check.memory=FALSE)
  overlap <-   hypervolume_overlap_statistics(hv_set)%>% as.data.frame() %>% t()
  ## MD and DR 
  dis_centroid  <-  hypervolume_distance (spe[[1]],spe[[2]], type='centroid')
  names( dis_centroid ) <- 'dis_centroid'
  ### calculate the distances in hypervolume space between the different cliamtes 
  Bhattacharyya <- MVNH_dissimilarity(db1 = spe[[1]]@Data, db2=spe[[2]]@Data) #multivariate normal hypervolume framework - 
  Bhattacharyya_distance <-  Bhattacharyya$Bhattacharyya_distance[1]
  names(Bhattacharyya_distance) <- 'Bhattacharyya_distance'
  Mahalanobis_distance <-  Bhattacharyya$Mahalanobis_distance[1]
  names( Mahalanobis_distance) <- ' Mahalanobis_distance'
  Determinant_ratio <-  Bhattacharyya$Determinant_ratio[1]
  names(Determinant_ratio) <- 'Determinant_ratio'
  #make a dataframe of it 
  local <-  cbind.data.frame(  volume_sp,   overlap ,  dis_centroid ,   Bhattacharyya_distance ,  Mahalanobis_distance,   Determinant_ratio)
  local$species <- name
  stats <- rbind(stats, local)
}
chelsa_stats <- stats
chelsa_stats$Dataset <- 'Chelsa'

stats <- NULL
hypervolume <- hypervolume_soiltemp
for (j in 1:length(Species)) {
  spe <-   hypervolume[[j]] 
  name <-  names( hypervolume[j] )
  ##first, have the sizes of the different hypervolumes 
  # get volume for each hypervolume 
  volume_sp <- get_volume(spe) %>% as.data.frame() %>% t()
  colnames(  volume_sp) <- c('volume_nat', 'volume_inv')
  ##then assess centroid distance and overlap  
  hv_set <- hypervolume_set(spe[[1]],spe[[2]], check.memory=FALSE)
  overlap <-   hypervolume_overlap_statistics(hv_set)%>% as.data.frame() %>% t()
  ## MD and DR 
  dis_centroid  <-  hypervolume_distance (spe[[1]],spe[[2]], type='centroid')
  names( dis_centroid ) <- 'dis_centroid'
  bat_dissim <- kernel.beta(spe) 
  ### calculate the distances in hypervolume space between the different cliamtes 
  Bhattacharyya <- MVNH_dissimilarity(db1 = spe[[1]]@Data, db2=spe[[2]]@Data) #multivariate normal hypervolume framework - 
  # #can assess relarive contributions of consiruent components in driving niche variation- integrate size and dissimilarity
  Bhattacharyya_distance <-  Bhattacharyya$Bhattacharyya_distance[1]
  names(Bhattacharyya_distance) <- 'Bhattacharyya_distance'
  Mahalanobis_distance <-  Bhattacharyya$Mahalanobis_distance[1]
  names( Mahalanobis_distance) <- ' Mahalanobis_distance'
  Determinant_ratio <-  Bhattacharyya$Determinant_ratio[1]
  names(Determinant_ratio) <- 'Determinant_ratio'

  #make a dataframe of it 
  local <-  cbind.data.frame(  volume_sp,   overlap ,  dis_centroid ,    Bhattacharyya_distance ,  Mahalanobis_distance,   Determinant_ratio)
  local$species <- name
  stats <- rbind(stats, local)
}
soil_stats <- stats
soil_stats$Dataset <- 'Soiltemps'

hypervolume_results <- rbind(soil_stats, chelsa_stats)
###### ###### ###### ###### ###### ###### ###### ###### ###### 


#######################################

### Step 3 - figures and analysis 
###Micrcoclimate levels - climate data levels PART 3
## ordination data 
#ordination_results  <- ordination_restrict  # to use if you want the restricted model 
ordination_results  <-  ordination_results  %>%   dplyr::select(X, Chelsa_D.overlap,  Chelsa_expansion,  soil_D.overlap, soil_expansion)
# X = Species names 
# worldclim.D.overlap - D overlap calculated using WorldClim
# worldclim_expansion - Expansion metric calculated using WorldClim
#"soil_D.overlap"    - D overlap calculated using SoilTemp
#"soil_expansion"     - - Expansion metric calculated using SoilTemp

## hypervolume data 
#hypervolume_results  <- hypervolume_restrict # to use if you want the restricted model 
hypervolume_results <-  hypervolume_results %>% dplyr::select(species, Dataset, Determinant_ratio, Mahalanobis_distance, Bhattacharyya_distance, 
                                                              dis_centroid, jaccard, volume_nat, volume_inv)
hypervolume_results <- hypervolume_results  %>%
              pivot_wider(names_from = Dataset, 
              values_from = c(Determinant_ratio, Mahalanobis_distance, Bhattacharyya_distance, dis_centroid, jaccard, volume_nat, volume_inv ))

############## ######
### Figure 2 ###
############## ######
#D overlap   
shapiro.test(ordination_results$Chelsa_D.overlap) #significant - need to use non-paremetric tests
shapiro.test(ordination_results$soil_D.overlap) #significant - need to use non-paremetric tests
cor.test( log(ordination_results$Chelsa_D.overlap), log(ordination_results$soil_D.overlap), method=c("kendall")) # significant correlation
d <- ggscatter(ordination_results, x = 'Chelsa_D.overlap', y = 'soil_D.overlap',# label='Species', 
               add = "reg.line", conf.int = TRUE, 
                cor.coef = F, cor.method = "kendall", digits = 3, size=0.6, cor.coef.size=3.3, repel = TRUE, font.label=c(10, 'plain', 'black'), font.family='Helvetica') +
          expand_limits(x =c( 0, 1), y = c( 0, 1)) + 
          geom_abline(intercept = 0, slope = 1, color="red", linewidth=0.6, alpha=0.5) +
          theme(text = element_text(family = "Helvetica", size=10)) +
          stat_cor( p.accuracy = 0.001, r.accuracy = 0.01,  method = "kendall", digits=2, label.x= 0, label.y=1.1, size = 3) + 
          coord_fixed() + 
          scale_y_continuous(limits=c(0,1))+
          scale_x_continuous(limits=c(0,1))+
          xlab('Macroclimate D Overlap') +  ylab(expression("Microclimate D overlap "))  + 
          labs(color = "Percentage of \nnative niche")   +
          annotate(geom = "polygon", x = c(-Inf, Inf, -Inf), y = c(-Inf, Inf, Inf), fill = "#458769", alpha = 0.2) +
          annotate(geom = "polygon", x = c(-Inf, Inf, Inf), y = c(-Inf, Inf, -Inf), fill = "#A2B9C2", alpha = 0.2 )   

#Expansion 
shapiro.test(ordination_results$Chelsa_D.overlap) #significant - need to use non-paremetric tests
shapiro.test(ordination_results$soil_expansion) #significant - need to use non-paremetric tests
cor.test(ordination_results$Chelsa_expansion, ordination_results$soil_expansion, method=c("kendall")) # there is a significant correlation between expansions
exp <- ggscatter(ordination_results, x = 'Chelsa_expansion', y = 'soil_expansion', 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = F, cor.method = "kendall", size=0.6, cor.coef.size=3.3, font.label=c(10, 'plain', 'black'), font.family='Arial') +
            theme(text = element_text(family = "Helvetica", size=10)) +
            stat_cor( p.accuracy = 0.001, r.accuracy = 0.01, method = "kendall", digits=2, label.x= 0, label.y=1.1, size = 3) + 
            scale_y_continuous(limits=c(0,1))+
            scale_x_continuous(limits=c(0,1))+
            expand_limits(x = c(0,1), y =c( 0,1)) + 
            coord_fixed() + 
            xlab('Macroclimate Expansion') +  ylab(expression("Microclimate Expansion"))  + 
            labs(color = "Percentage of \nnative niche") + 
            annotate(geom = "polygon", x = c(-Inf, Inf, -Inf), y = c(-Inf, Inf, Inf), fill ="#A2B9C2", alpha = 0.2) +
            annotate(geom = "polygon", x = c(-Inf, Inf, Inf), y = c(-Inf, Inf, -Inf), fill = "#458769", alpha = 0.2 ) +
            annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = Inf, 
                     colour = "red", size=0.6, alpha=0.5)
### Hypervolume methods correaltion in niche metrics between datasets 
#Jaccard Similarity 
shapiro.test(hypervolume_results$jaccard_Chelsa) 
shapiro.test(hypervolume_results$jaccard_Soiltemps) 
cor.test(hypervolume_results$jaccard_Chelsa, hypervolume_results$jaccard_Soiltemps, method=c("kendall")) # significant correlation too 
Jac <- ggscatter(hypervolume_results, x = 'jaccard_Chelsa', y = 'jaccard_Soiltemps',# label='Species', 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = F, cor.method = "pearson", size=0.6, cor.coef.size=3.3, repel = TRUE, font.label=c(10, 'plain', 'black'), font.family='Helvetica') +
          stat_cor( p.accuracy = 0.001, r.accuracy = 0.01, method = "kendall", digits=2, label.x= 0, label.y.npc = "top", size = 3) + 
          expand_limits(x = c(0,1), y =c( 0,1)) + 
          geom_abline(intercept = 0, slope = 1, color="red", size=0.6,  alpha=0.5) +
          theme(text = element_text(family = "Helvetica", size=10)) +
          xlab('Macroclimate  Jaccard Similarity') +  ylab(expression("Microclimate Jaccard Similarity"))  + 
          labs(color = "Percentage of \nnative niche")   +
          coord_fixed() + 
          annotate(geom = "polygon", x = c(-Inf, Inf, -Inf), y = c(-Inf, Inf, Inf), fill ="#458769", alpha = 0.2) +
          annotate(geom = "polygon", x = c(-Inf, Inf, Inf), y = c(-Inf, Inf, -Inf), fill = "#A2B9C2", alpha = 0.2 ) +
          annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = Inf,  colour = "red", size=0.6, alpha=0.5) 
# boxplots 
positive  <- function(x){ # functions determines if number is positive, negative or zero 
                if(x > 0) {
                  base::print("Positive")
                } else {
                  if(x == 0) {
                    base::print("Zero")
                  } else {
                    base::print("Negative")
                  }
                }}

### D overlap - Figure 2D
ordination_results$Ddifference <- ordination_results$Chelsa_D.overlap - ordination_results$soil_D.overlap
new <- lapply(round(ordination_results$Ddifference, digits = 2), FUN=positive) %>% unlist() %>% as.vector()
ordination_results$is_increasing <-new
plot_data <- melt(ordination_results, id.vars=c('species', 'is_increasing'), measure.vars=c('Chelsa_D.overlap', 'soil_D.overlap'))
kruskal.test(plot_data$value, plot_data$variable)
compare_means(value ~ variable, data = plot_data, paired = TRUE, method='wilcox.test')
length(which(          ordination_results$is_increasing=='Negative'))/95

d_box <- ggplot(plot_data, aes(x=variable, y=value, fill=variable, colour=variable)) +
                  geom_boxplot(aes(x = rep(c(1, 2), each = 95), group = variable),  outlier.shape=NA, width=0.2, show.legend = F) + 
                  scale_fill_manual(values=c( "#A2B9C2", "#8CC2AA"), labels=c('Chelsa', "SoilTemp")) + 
                  scale_colour_manual(values=c( "#50808E", "#6B917F") )+ 
                  geom_point(aes(group=species, fill=variable, x = rep(c(1.25, 1.75), each = 95)), shape=19,  
                             position = position_dodge(0.01), size=0.7, colour='black', show.legend = F) + 
                  new_scale_color() +
                  scale_colour_manual(values = c("#3E541E", "#BB6B59", 'grey34')) + 
                  geom_line(aes(group=species, col = is_increasing, x = rep(c(1.25, 1.75), each = 95)), alpha=0.5, size=0.3, position = position_dodge(0.01))  + 
                  stat_compare_means(paired = TRUE, label='p.signif',   label.x.npc = "middle") + 
                  xlab('Dataset') +  ylab('D overlap ') + 
                  scale_x_discrete(labels=c('Macroclimate', "Microclimate")) + 
                  theme_classic()  + 
                  labs(colour='Direction of change') + 
                  theme(axis.text = element_text(family = "Helvetica", size=10)) + 
                  theme(text = element_text(family = "Helvetica", size=10)) + 
                  theme(axis.text.x=element_text(colour="black")) 
### Jaccard dissimilarity - Figure 2F
hypervolume_results$Jdifference <- hypervolume_results$jaccard_Chelsa- hypervolume_results$jaccard_Soiltemps
new <- lapply(round(hypervolume_results$Jdifference, digits = 2), FUN=positive)%>% unlist() %>% as.vector()
hypervolume_results$is_increasing <-new
length(which(           hypervolume_results$is_increasing=='Negative'))/95
#statistics
plot_data <-  reshape2::melt(hypervolume_results, id.vars=c('species', 'is_increasing'), 
                             measure.vars=c('jaccard_Chelsa', 'jaccard_Soiltemps'))
compare_means(value ~ variable, data = plot_data, paired = TRUE, method='wilcox.test')
#plot
jac_box <-  ggplot( plot_data, aes(x=variable, y=value, fill=variable, colour=variable)) +
                  geom_boxplot(aes(x = rep(c(1, 2), each = 95), group = variable),  outlier.shape=NA, width=0.2, show.legend = F) + 
                  scale_fill_manual(values=c( "#A2B9C2", "#8CC2AA"), labels=c('Chelsa', "SoilTemp")) + 
                  scale_colour_manual(values=c( "#50808E", "#6B917F") )+ 
                  geom_point(aes(group=species, fill=variable, x = rep(c(1.25, 1.75), each = 95)), shape=19,  
                             position = position_dodge(0.01), size=0.7, colour='black', show.legend = F) + 
                  new_scale_color() +
                  scale_colour_manual(values = c("#3E541E", "#BB6B59", 'grey34')) + 
                  geom_line(aes(group=species, col = is_increasing, x = rep(c(1.25, 1.75), each = 95)), alpha=0.5, size=0.3, position = position_dodge(0.01))  + 
                  stat_compare_means(paired = TRUE, label='p.signif',   label.x.npc = "middle") + 
                  xlab('Dataset') +  ylab('Jaccard Similarity ') + 
                  scale_x_discrete(labels=c('Macroclimate', "Microclimate")) + 
                  theme_classic()  + 
                  labs(colour='Direction of change') + 
                  theme(text = element_text(family = "Helvetica", size=10)) + 
                  theme(axis.text = element_text(family = "Helvetica", size=10)) + 
                  theme(axis.text.x=element_text(colour="black")) 
##  Expansion - significance/non significance  - Figure 2E
## raw differences 
ordination_results$Exp_difference <- ordination_results$Chelsa_expansion - ordination_results$soil_expansion
length(which( ordination_results$Exp_difference == 0))  # 12 species had no change 
# similar numbers of species have significant expansions, but the identity of those species change. 
ordination_results$exp_change <- ifelse(ordination_results$Chelsa_expansion > 0.1 & ordination_results$soil_expansion > 0.1, 'still sig', 'change' )
ordination_results$exp_change_2 <- ifelse(ordination_results$Chelsa_expansion < 0.1 & ordination_results$soil_expansion < 0.1, 'still nonsig', 'change' )
## see change in each (Figure 2E)
ordination_results$exp_cat_world <- ifelse( ordination_results$Chelsa_expansion > 0.1, 'sig', 'nonsig')
ordination_results$exp_cat_soil <- ifelse( ordination_results$soil_expansion> 0.1, 'sig', 'nonsig')
plot_data1  <-  rbind( as.data.frame(table(ordination_results$exp_cat_world)) ,  as.data.frame(table(ordination_results$exp_cat_soil)))
plot_data1$dataset <- c('Chelsa', 'Chelsa', 'SoilTemp', 'SoilTemp')
levels(plot_data1$dataset) <- c('Chelsa', 'SoilTemp')
colours <- c( "#44987D", '#D57578')
## statistical test 
test <-  plot_data1 %>%  dcast(dataset~Var1, value.var='Freq')
chisq.test(test[,c(2, 3)])  ###Significant
# plot figure 2E
fig2e <-  ggplot(plot_data1 , aes(y=Freq, fill=Var1, x=dataset)) +
                geom_bar(position="stack", stat='identity') + 
                scale_fill_manual(values= colours, labels=c("0-10% ", '10-100%') ) +
                labs(fill = "Expansion \nconclusion") + 
                theme_classic()+
                annotate("text", x = 1.5, y = 100,   label = "*")+ 
                ylab('Number of Species')+ xlab('Dataset') + 
                theme(axis.text = element_text(family = "Helvetica", size=10, colour='black'))  + 
                scale_x_discrete(limits = levels(plot_data1$dataset)) + 
                theme(legend.text = element_text(family = "Helvetica", size=10), legend.title=element_text(family = "Helvetica", size=10)) 
figure2_final <- ggarrange( d,    exp,   Jac,d_box, fig2e, jac_box, common.legend = T, labels=c('A', "B", "C", 'D', 'E', 'F'))
# path = file.path("Graphs/fig2_1.pdf")
# ggsave(filename = path, plot=figure2_final, width=7.4, height=6 , units = "in", device='pdf')


### #### #### #### #### 
## Figure 3 ####
######## #### #### #### 
# combine 
index <- match(hypervolume_results$species, ordination_results$species)
hypervolume_results$Ddifference <- ordination_results$Ddifference[index]
hypervolume_results$expcategories <- ordination_results$expcategories[index]
hypervolume_results  <- hypervolume_results  %>% arrange(species)
### expansion figure ###
#D overlap
hypervolume_results  <- hypervolume_results  %>% arrange(desc(Ddifference))
hypervolume_results$species <- factor(  hypervolume_results$species, levels=hypervolume_results$species)

background_area_name <- c("Soiltemps", "Chelsa")
background_data <- data.frame(xmin = -Inf, xmax = Inf,
                              ymin = c(-Inf, 0),
                              ymax = c(0, Inf),
                              fill = factor(background_area_name,
                                            background_area_name))
dplot <- ggplot(data = hypervolume_results, aes(x = species, y=Ddifference)) +
                geom_bar(stat='identity', fill='gray15') +
                geom_rect(data = background_data, 
                          aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax, fill = fill),
                          inherit.aes = FALSE, alpha = 0.5) +
                scale_fill_manual(values= c( '#458769', '#A2B9C2' )) + 
                geom_bar(stat='identity', fill='gray15') +
                theme_classic()+
                coord_flip() + 
                ylab("D overlap difference \n(Macro-Micro)") + xlab("Species") + 
                rremove("ylab") + 
                theme(axis.text.y = element_text(face = "italic")) + 
                theme(plot.background = element_rect(fill="transparent", color=NA),
                      panel.background = element_rect(fill="transparent", color=NA)) + 
                theme(text = element_text(family = "Helvetica", size=10)) +
                theme(text=element_text(color="gray15"), axis.text=element_text(color="grey15"), axis.line=element_line(colour='grey15')) + 
                geom_hline(yintercept=0.1, lty=2, colour='#ffb6c1') +
                geom_hline(yintercept=-0.1, lty=2, colour='#ffb6c1') 
dplot
Ddifference <- ordination_results %>%  dplyr::select(Ddifference, X)

# Jaccard differences 
hypervolume_results$species <- factor(  hypervolume_results$species, levels=hypervolume_results$species)
background_area_name <- c("Soiltemps", "Chelsa")
background_data <- data.frame(xmin = -Inf, xmax = Inf,
                              ymin = c(-Inf, 0),
                              ymax = c(0, Inf),
                              fill = factor(background_area_name,
                                            background_area_name))
Jplot <- ggplot(data = hypervolume_results, aes(x = species, y=Jdifference)) +
                geom_bar(stat='identity', fill='gray15') +
                geom_rect(data = background_data, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax, fill = fill),
                          inherit.aes = FALSE, alpha = 0.5) +
                scale_fill_manual(values= c(  '#458769',  '#A2B9C2')) + 
                geom_bar(stat='identity', fill='gray15') +
                theme_classic()+
                coord_flip() + 
                rremove("ylab") + 
                theme( axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
                ylab("Jaccard similarity\n(Macro-Micro)") + xlab("Species") + 
                theme(plot.background = element_rect(fill="transparent", color=NA),
                      panel.background = element_rect(fill="transparent", color=NA)) + 
                theme(text = element_text(family = "Helvetica", size=10)) +
                theme(text=element_text(color="gray15"), axis.text=element_text(color="grey15"), axis.line=element_line(colour='grey15'))
Jplot # for more species might higher similarity for jaccard in Soil 

## save the figure together as one 
plot3 <- ggarrange(  dplot, Jplot,  labels=c("A", "B" ), common.legend = T, widths=c(2,1.5))
plot3
#path = file.path("Graphs/fig3.pdf")
#ggsave(path, plot = plot3, width = 180, height = 260, units = "mm")

########
### for Figure 4 ###
########
#### Load data
## Loading in the species and world data
Species <- names(species_data)

## for each species, compare the correlation between chelsa and sbio, and the varience between chelsa and sbio
table1 <- NULL 
tableVar <- NULL
for (i in 1:length(species_data)) {
  local <- species_data[[i]]
  spe <-names(species_data[i])
  bio1 <- cor.test(local$SBIO1_0_5cm_Annual_Mean_Temperature, local$CHELSA_bio10_01, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio2 <- cor.test(local$SBIO2_0_5cm_mean_diurnal_range, local$CHELSA_bio10_02, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio3 <-cor.test(local$SBIO3_0_5cm_Isothermality, local$CHELSA_bio10_03, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio4 <-cor.test(local$SBIO4_0_5cm_Temperature_Seasonality, local$CHELSA_bio10_04, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio5 <-cor.test(local$SBIO5_0_5cm_MaxT_warmestMonth, local$CHELSA_bio10_05, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio6 <-cor.test(local$SBIO6_0_5cm_MinT_coldestMonth, local$CHELSA_bio10_06, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio7 <-cor.test(local$SBIO7_0_5cm_annual_range, local$CHELSA_bio10_07, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio8 <-cor.test(local$SBIO8_0_5cm_meanT_wettestQ, local$CHELSA_bio10_08, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio9 <-cor.test(local$SBIO9_0_5cm_meanT_driestQ, local$CHELSA_bio10_09, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio10 <-cor.test(local$SBIO10_0_5cm_meanT_warmestQ, local$CHELSA_bio10_10, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio11 <-cor.test(local$SBIO11_0_5cm_meanT_coldestQ, local$CHELSA_bio10_11, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  #join all the data together 
  es <- cbind(bio1$estimate, bio2$estimate, bio3$estimate, bio4$estimate, bio5$estimate, bio6$estimate, bio7$estimate, bio8$estimate, bio9$estimate, bio10$estimate, bio11$estimate)
  colnames(es) <- c('Bio1', 'Bio2', 'Bio3', 'Bio4', 'Bio5', 'Bio6', 'Bio7', 'Bio8', 'Bio9', 'Bio10', 'Bio11')
  es
  rownames(es) <- spe 
  table1 <- rbind.data.frame(table1, es)
  
  # varience 
  vbio1 <- var(local$SBIO1_0_5cm_Annual_Mean_Temperature, local$CHELSA_bio10_01) # sinnificant correlation too 
  vbio2 <- var(local$SBIO2_0_5cm_mean_diurnal_range, local$CHELSA_bio10_02) # sinnificant correlation too 
  vbio3 <-var(local$SBIO3_0_5cm_Isothermality, local$CHELSA_bio10_03) # sinnificant correlation too 
  vbio4 <-var(local$SBIO4_0_5cm_Temperature_Seasonality, local$CHELSA_bio10_04) # sinnificant correlation too 
  vbio5 <-var(local$SBIO5_0_5cm_MaxT_warmestMonth, local$CHELSA_bio10_05) # sinnificant correlation too 
  vbio6 <-var(local$SBIO6_0_5cm_MinT_coldestMonth, local$CHELSA_bio10_06) # sinnificant correlation too 
  vbio7 <-var(local$SBIO7_0_5cm_annual_range, local$CHELSA_bio10_07) # sinnificant correlation too 
  vbio8 <-var(local$SBIO8_0_5cm_meanT_wettestQ, local$CHELSA_bio10_08) # sinnificant correlation too 
  vbio9 <-var(local$SBIO9_0_5cm_meanT_driestQ, local$CHELSA_bio10_09) # sinnificant correlation too 
  vbio10 <-var(local$SBIO10_0_5cm_meanT_warmestQ, local$CHELSA_bio10_10) # sinnificant correlation too 
  vbio11 <-var(local$SBIO11_0_5cm_meanT_coldestQ, local$CHELSA_bio10_11) # sinnificant correlation too 
  #join all the data together 
  va <- cbind(vbio1, vbio2, vbio3, vbio4, vbio5, vbio6, vbio7, vbio8, vbio9, vbio10, vbio11)
  colnames(va) <- c('Bio1', 'Bio2', 'Bio3', 'Bio4', 'Bio5', 'Bio6', 'Bio7', 'Bio8', 'Bio9', 'Bio10', 'Bio11')
  rownames(va ) <- spe 
  tableVar <- rbind.data.frame(  tableVar,   va)
  
}
#write.csv(table1, '../Data/Results/correlationMacroMicroOcc.csv')
table2 <- gather(table1 , variable, correlation, Bio1:Bio11)
table2$variable <- factor(table2$variable, levels=c('Bio1','Bio2', 'Bio3', 'Bio4', 'Bio5', 'Bio6', 'Bio7', 'Bio8', 'Bio9', 'Bio10', 'Bio11'))

max(table1[,1:11]); min(table1[,1:11])
median(table1$Bio1); range(table1$Bio1)
median(table1$Bio2); range(table1$Bio2)
median(table1$Bio3);range(table1$Bio3)
median(table1$Bio4); range(table1$Bio4)
median(table1$Bio5); range(table1$Bio5)
median(table1$Bio6); range(table1$Bio6)
median(table1$Bio7); range(table1$Bio7)
median(table1$Bio8); range(table1$Bio8)
median(table1$Bio9); range(table1$Bio9)
median(table1$Bio10); range(table1$Bio10)
median(table1$Bio11); range(table1$Bio11)

hist(table1$Bio11)
table1$species <- rownames(table1)
#statistical differences between observed correlations  - correlation
test <-  melt( table1, id.vars=c("species"))
kruskal.test(test$value, test$variable) # significant 
pairwise.wilcox.test(test$value, test$variable, p.adjust.method='BH')
#statistical differences between observed correlations  - varitation
tableVar$species <- rownames( tableVar)
test <-  melt(  tableVar, id.vars=c("species"))
kruskal.test(test$value, test$variable) # significant 
pairwise.wilcox.test(test$value, test$variable, p.adjust.method='BH')


#### then compared observed correlaitons/ varience to randomly generated variables , this takes a long time to run can use the vairbale 'tab' instead
ok <- st_as_sf(w_rworldmap)
#SBio <- rast(SBio)
#r <- rast(r)
# load in chelsa data downloaded from https://chelsa-climate.org/exchelsa-extended-bioclim/
rastlist <- list.files(path = "~/Desktop/Chelsa/", pattern='.tif',  all.files=TRUE, full.names=T)
r <- stack(rastlist)
crs(r) <- crs(' +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
r <- r[[1:11]]
r <- rast(r)
# load in data downloaded from https://zenodo.org/records/7134169 (version 2)
rastlist <- list.files(path = "SBIO/Raw", pattern='.tif',  all.files=TRUE, full.names=T)
SBio <- stack(rastlist)
crs(SBio) <- crs(' +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
SBio <- rast(SBio)
SBio
# making sure the two rasters are exactly the same 
r <- crop(r, SBio)
r <- terra::project(r, SBio)
SBio <- terra::project( SBio, r)
r <- terra::mask(r, SBio)
all.equal(r, SBio)


tabler <- NULL 
tableVar2 <- NULL 
# do correlation of ALL pixels 
#example sepceis with the plot of the whole world, xy on and aggregated map and then plot one random dample where they are on the correlartion plot then plot one species, and two variables- bio 1 and 3 and see the effect to illustrate whats happning 
for (i in 1:1000) {
  set.seed(i)
  points <- st_sample(ok, size = 1000, type = "random", crs = st_crs(4326)) 
  coordinatespoints <-  st_coordinates(points) %>% 
    as.data.frame()
  df <- raster::extract(r, coordinatespoints)
  df2 <- raster::extract(SBio, coordinatespoints)
  df3 <- cbind(df, df2)
  local<- df3[complete.cases(df3), ] %>% as.data.frame()
  bio1 <- cor.test(local$SBIO1_0_5cm_Annual_Mean_Temperature, local$CHELSA_bio10_01, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio2 <- cor.test(local$SBIO2_0_5cm_mean_diurnal_range, local$CHELSA_bio10_02, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio3 <-cor.test(local$SBIO3_0_5cm_Isothermality, local$CHELSA_bio10_03, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio4 <-cor.test(local$SBIO4_0_5cm_Temperature_Seasonality, local$CHELSA_bio10_04, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio5 <-cor.test(local$SBIO5_0_5cm_MaxT_warmestMonth, local$CHELSA_bio10_05, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio6 <-cor.test(local$SBIO6_0_5cm_MinT_coldestMonth, local$CHELSA_bio10_06, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio7 <-cor.test(local$SBIO7_0_5cm_annual_range, local$CHELSA_bio10_07, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio8 <-cor.test(local$SBIO8_0_5cm_meanT_wettestQ, local$CHELSA_bio10_08, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio9 <-cor.test(local$SBIO9_0_5cm_meanT_driestQ, local$CHELSA_bio10_09, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio10 <-cor.test(local$SBIO10_0_5cm_meanT_warmestQ, local$CHELSA_bio10_10, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  bio11 <-cor.test(local$SBIO11_0_5cm_meanT_coldestQ, local$CHELSA_bio10_11, method=c("kendall"), alternative='two.sided') # sinnificant correlation too 
  #join all the data together 
  es <- cbind(bio1$estimate, bio2$estimate, bio3$estimate, bio4$estimate, bio5$estimate, bio6$estimate, bio7$estimate, bio8$estimate, bio9$estimate, bio10$estimate, bio11$estimate)
  colnames(es) <- c('Bio1', 'Bio2', 'Bio3', 'Bio4', 'Bio5', 'Bio6', 'Bio7', 'Bio8', 'Bio9', 'Bio10', 'Bio11')
  es
  tabler <- rbind.data.frame(tabler, es)
  
  
  # varience 
  vbio1 <- var(local$SBIO1_0_5cm_Annual_Mean_Temperature, local$CHELSA_bio10_01) # sinnificant correlation too 
  vbio2 <- var(local$SBIO2_0_5cm_mean_diurnal_range, local$CHELSA_bio10_02) # sinnificant correlation too 
  vbio3 <-var(local$SBIO3_0_5cm_Isothermality, local$CHELSA_bio10_03) # sinnificant correlation too 
  vbio4 <-var(local$SBIO4_0_5cm_Temperature_Seasonality, local$CHELSA_bio10_04) # sinnificant correlation too 
  vbio5 <-var(local$SBIO5_0_5cm_MaxT_warmestMonth, local$CHELSA_bio10_05) # sinnificant correlation too 
  vbio6 <-var(local$SBIO6_0_5cm_MinT_coldestMonth, local$CHELSA_bio10_06) # sinnificant correlation too 
  vbio7 <-var(local$SBIO7_0_5cm_annual_range, local$CHELSA_bio10_07) # sinnificant correlation too 
  vbio8 <-var(local$SBIO8_0_5cm_meanT_wettestQ, local$CHELSA_bio10_08) # sinnificant correlation too 
  vbio9 <-var(local$SBIO9_0_5cm_meanT_driestQ, local$CHELSA_bio10_09) # sinnificant correlation too 
  vbio10 <-var(local$SBIO10_0_5cm_meanT_warmestQ, local$CHELSA_bio10_10) # sinnificant correlation too 
  vbio11 <-var(local$SBIO11_0_5cm_meanT_coldestQ, local$CHELSA_bio10_11) # sinnificant correlation too 
  va <- cbind(vbio1, vbio2, vbio3, vbio4, vbio5, vbio6, vbio7, vbio8, vbio9, vbio10, vbio11)
  colnames(va) <- c('Bio1', 'Bio2', 'Bio3', 'Bio4', 'Bio5', 'Bio6', 'Bio7', 'Bio8', 'Bio9', 'Bio10', 'Bio11')
  tableVar2 <- rbind.data.frame(  tableVar2,   va)
  
  base::print(paste0('done', as.numeric(i)))
}

tabler2 <- gather(tabler , variable, correlation, Bio1:Bio11)
tabler2$variable <- factor(tabler2$variable, levels=c('Bio1','Bio2', 'Bio3', 'Bio4', 'Bio5', 'Bio6', 'Bio7', 'Bio8', 'Bio9', 'Bio10', 'Bio11'))
#wide to long format for plotting
tableVar2 <- gather( tableVar, variable, correlation, Bio1:Bio11)
tableVar2$variable <- factor(tableVar2$variable, levels=c('Bio1','Bio2', 'Bio3', 'Bio4', 'Bio5', 'Bio6', 'Bio7', 'Bio8', 'Bio9', 'Bio10', 'Bio11'))

### figure 4
table2$data <- 'Species'
tabler2$data <- 'Random'
#tabler2$species <- 'random'
tab <- rbind(table2, tabler2)

tab$merged <- factor(paste( tab$data, tab$variable,sep='_'),
                     levels=c("Random_Bio1" ,  "Random_Bio2" ,"Random_Bio3" ,"Random_Bio4","Random_Bio5",  
                              "Random_Bio6","Random_Bio7","Random_Bio8","Random_Bio9","Random_Bio10" ,"Random_Bio11",
                              "Species_Bio1" ,"Species_Bio2" ,"Species_Bio3" ,"Species_Bio4","Species_Bio5",  
                              "Species_Bio6", "Species_Bio7" ,"Species_Bio8" ,
                              "Species_Bio9","Species_Bio10","Species_Bio11" ))
tab1 <- tab %>%  filter(variable=='Bio1') 
kruskal.test(tab1$correlation, tab1$data) # significant 
tab1 <-tab %>%  filter(variable=='Bio2')  
kruskal.test(tab1$correlation, tab1$data) # significant 
tab1 <-tab %>%  filter(variable=='Bio3') 
kruskal.test(tab1$correlation, tab1$data) # significant 
tab1 <-tab %>%  filter(variable=='Bio4') 
kruskal.test(tab1$correlation, tab1$data) # not significant 
tab1 <-tab %>%  filter(variable=='Bio5') 
kruskal.test(tab1$correlation, tab1$data) # significant 
tab1 <-tab %>%  filter(variable=='Bio6') 
kruskal.test(tab1$correlation, tab1$data) # significant 
tab1 <-tab %>%  filter(variable=='Bio7')
kruskal.test(tab1$correlation, tab1$data) # significant 
tab1 <-tab %>%  filter(variable=='Bio8') 
kruskal.test(tab1$correlation, tab1$data) # significant 
tab1 <-tab %>%  filter(variable=='Bio9')
kruskal.test(tab1$correlation, tab1$data) # significant 
tab1 <-tab %>%  filter(variable=='Bio10')
kruskal.test(tab1$correlation, tab1$data) # significant 
tab1 <-tab %>%  filter(variable=='Bio11') 
kruskal.test(tab1$correlation, tab1$data) # significant 

colp <- c(rep('grey20',11),viridis(11,begin=0.2, end=0.85))

plot4 <- ggplot(tab, aes(x=variable, y=correlation, fill=merged, col=merged) )+ 
                  ggdist::stat_halfeye(normalize = "groups",adjust = .5,  width = .7, .width = 0, 
                                       justification = -.2, point_colour = NA, alpha=0.7, outline_bars=T) + 
                  geom_jitter(width = .05, alpha = .2) + 
                  geom_boxplot(width = .08, outlier.shape = NA, alpha=0.7, position=position_dodge(0.2), fill='white') + #, colour='black') + 
                  # geom_density_ridges(aes(height =..scaled..)) + 
                  theme_classic()  + 
                  #expand_limits(y =c( -0.5, 1)) + 
                  scale_x_discrete(limits = rev(levels(tab$variable))) + 
                  ylab('Per Species Correlation coefficient') +
                  xlab('Bioclimatic Variable') +
                  scale_fill_manual(values=colp) +
                  scale_colour_manual(values=colp) +
                  theme(legend.position='none', text = element_text(family = "Helvetica", size=10)) + 
                  theme(axis.text.x=element_text(colour="black")) + 
                  coord_flip()
path = file.path("/Graphs/fig4.pdf")
ggsave(filename = path, plot=plo4t, width=6, height=6 , units = "in", device='pdf')


### Figure 5 ###
# comparing Correlation between metrics 
test_data <- hypervolume_results %>% dplyr::select(species,BDdifference, disCdifference, Ddifference, expcategories, 
                                                   Jdifference  ) #jdifference too 
test_data$species <- gsub( " ", "_",   test_data$species )

### Ant nesting types 
index <- match(test_data$species, nesting$species)
test_data$NestingTypes <-  nesting$Cat_simp[index]
#plot it 
d <- ggplot( test_data[complete.cases(test_data$NestingTypes),], aes(x=NestingTypes, y=Ddifference, fill=NestingTypes)) +
  geom_boxplot(color='black', outlier.shape = NA) + 
  geom_jitter(shape=19, position=position_jitter(0.2), size=0.7) + 
  scale_fill_manual(values=c("#55853D", "#764108", "#DB5E48"), labels=c("Aboreal", "Ground", "Generalist") ) + 
  xlab('Nesting type') +  ylab('D Overlap difference ') + 
  scale_x_discrete(labels = c("Aboreal", "Ground", "Generalist")) + 
  theme_classic()  + 
  theme(legend.position = 'none') + 
  theme(text = element_text(family = "Helvetica", size=10)) + 
  theme(axis.text.x=element_text(colour="black")) 
kruskal.test(test_data$NestingTypes, test_data$Ddifference)

jac <- ggplot(test_data[complete.cases(test_data$NestingTypes),], aes(x=NestingTypes, y=Jdifference, fill=NestingTypes)) +
  geom_boxplot(color='black', outlier.shape = NA) + 
  scale_fill_manual(values=c("#55853D", "#764108",  "#DB5E48"), labels=c("Aboreal", "Ground", "Generalist") ) + 
  geom_jitter(shape=19, position=position_jitter(0.2), size=0.7) + 
  xlab('Nesting type') +  ylab('Jaccard Similarity difference') + 
  scale_x_discrete(labels = c("Aboreal", "Ground", "Generalist")) + 
  theme_classic()  + 
  theme(legend.position = 'none') + 
  theme(text = element_text(family = "Helvetica", size=10)) + 
  theme(axis.text.x=element_text(colour="black")) 
kruskal.test(test_data$Jdifference, test_data$NestingTypes)

exp <- ggplot(test_data[complete.cases(test_data$NestingTypes),] , aes(fill=NestingTypes, x=expcategories)) +
  geom_bar() + 
  scale_fill_manual(values=c("#55853D", "#764108",  "#DB5E48"), labels=c("Aboreal", "Ground", "Generalist") )+
  theme_classic()  + 
  theme(legend.position =  c(0.2, 0.9)) + 
  scale_x_discrete(labels=c('Chelsa exp, \nSoilTemp no exp', 'Chelsa no exp, \nSoilTemp exp', 'No change \n(no exp', 'No change \n(exp)')) + 
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) + 
  theme(text = element_text(family = "Helvetica", size=10)) + 
  xlab('Expansion categories') +  ylab('Species (n) ') + 
  theme(axis.text.x=element_text(colour="black")) 
test.conclusion <- test_data  %>% dplyr::group_by(NestingTypes, expcategories) %>% 
  summarise(Freq = n())  %>%  as.data.frame()
test.conclusion <- spread(test.conclusion, expcategories , Freq)
test.conclusion <- test.conclusion [-which(is.na(test.conclusion$NestingTypes )),]

test.conclusion [is.na(test.conclusion )] <- 0
chisq.test(t(test.conclusion[,2:5]) ) # not significant 

## Habitat types 
Species <- names(species_data)
SpeciesLandconver <- NULL 
for (i in 1:length(Species)) {
  species <-species_data[[as.character(Species[i])]]
  ##seperate the distirbutions by native and exotic! (and then maybe by continant) 
  lat <-   species$lat
  lon <- species$lon
  points <- vect(species, geom=c("lon", "lat"), crs=crs(landcover), keepgeom=FALSE)
  #coords <- data.frame(lon=lon,lat=lat |> vect(crs = landcover@crs) )#creating a data frame of the latitudes and longitudes 
  #points <- SpatialPointsDataFrame(coords, data=species, proj4string = landcover@crs) #this converts the lat and longs into spatial data
  values <- terra::extract(landcover, points, bind=T) #extracts the bioclim data for the coodrinate points we have 
  # 1-5 == forests 
  # 6-7 == shrubs 
  # 8-10 # grasslands
  #percentages of all 
  lndcov <- values %>%  as.data.frame() %>% group_by(HYBMAP_IGBP_2013_LC) %>% 
    summarise(n= n()) %>% 
    mutate( perc = n/sum(n) *100)
  # percentages of forests
  forests <- lndcov[lndcov$HYBMAP_IGBP_2013_LC %in% c('1','2','3', '4', '5'),]  # just forests 
  forestcov <- sum(forests$perc)
  anyveg <- lndcov[lndcov$HYBMAP_IGBP_2013_LC %in% c('1','2','3', '4', '5', '6', '7'),]  # forests and shurblands
  anyvegcov <- sum(  anyveg$perc)
  urban <- lndcov[lndcov$HYBMAP_IGBP_2013_LC %in% c('13'),]  # urban areas 
  urbancov <- sum( urban$perc)
  #make df
  loc <- cbind.data.frame( forestcov,  anyvegcov ,   urbancov )
  
  ## then the same for just native 
  lndcov <- values[which(values$status=='N')] %>%  as.data.frame() %>% group_by(HYBMAP_IGBP_2013_LC) %>% 
    summarise(n= n()) %>% 
    mutate( perc = n/sum(n) *100)
  forests <- lndcov[lndcov$HYBMAP_IGBP_2013_LC %in% c('1','2','3', '4', '5'),] 
  forestcov <- sum(forests$perc)
  anyveg <- lndcov[lndcov$HYBMAP_IGBP_2013_LC %in% c('1','2','3', '4', '5', '6', '7'),] 
  anyvegcov <- sum(  anyveg$perc)
  urban <- lndcov[lndcov$HYBMAP_IGBP_2013_LC %in% c('13'),] 
  urbancov <- sum( urban$perc)
  #make df
  loc2 <- cbind.data.frame( forestcov,  anyvegcov ,   urbancov )
  colnames(loc2) <- c( 'forestcovNAT',  'anyvegcovNAT' ,   'urbancovNAT' )
  
  ## then the same for just Introduced
  lndcov <- values[which(values$status=='E')] %>%  as.data.frame() %>% group_by(HYBMAP_IGBP_2013_LC) %>% 
    summarise(n= n()) %>% 
    mutate( perc = n/sum(n) *100)
  forests <- lndcov[lndcov$HYBMAP_IGBP_2013_LC %in% c('1','2','3', '4', '5'),] 
  forestcov <- sum(forests$perc)
  anyveg <- lndcov[lndcov$HYBMAP_IGBP_2013_LC %in% c('1','2','3', '4', '5', '6', '7'),] 
  anyvegcov <- sum(  anyveg$perc)
  urban <- lndcov[lndcov$HYBMAP_IGBP_2013_LC %in% c('13'),] 
  urbancov <- sum( urban$perc)
  #make df
  loc3 <- cbind.data.frame( forestcov,  anyvegcov ,   urbancov )
  colnames(loc3) <- c( 'forestcovIN',  'anyvegcovIN' ,   'urbancovIN' )
  
  #join everything together 
  local <- cbind(loc, loc2, loc3)
  local$Species <- Species[i]
  SpeciesLandconver <- rbind(  SpeciesLandconver,   local)
}

# list of all species, the amount of landcover for all range,  and the native and invasive ranges seperatly 
SpeciesLandconver$Species <- gsub(" ", "_", SpeciesLandconver$Species)
index <- match(test_data$species,  as.factor(SpeciesLandconver$Species))
test_data$forestcov <-   SpeciesLandconver $forestcov[index]
test_data$forestcovDIFF <-  SpeciesLandconver $forestcovNAT[index] - SpeciesLandconver$forestcovIN[index]

# d overlap 
f_d <-ggscatter( test_data, x = 'forestcov', y = 'Ddifference',# label='Species', 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = F, cor.method = "pearson", size=0.4, cor.coef.size=1, ## correlation coefficient is Kendalls correlation Tau because data is non-normal
                 #  cor.coef.coord = c(0, 1.1),
                 repel = TRUE, font.label=c(10, 'plain', 'black'), 
                 font.family='Helvetica') +
                  stat_cor(method='pearson',
                  aes(label = paste(..r.label.., ..rr.label.., ..p.label.., sep = "~`,`~")), 
                  label.x = 3, size=3) + 
                  theme_minimal()+
                  theme(text = element_text(family = "Helvetica", size=10)) +
                  xlab('% under forest cover') +  ylab("D  Overlap difference")
cor(test_data$forestcovDIFF, test_data$Ddifference, method='pearson')

# expansion
f_e  <- ggplot(test_data, aes(x=expcategories, y=forestcov, fill=expcategories)) +
                geom_boxplot(color='black', outlier.shape = NA) + 
                geom_jitter(shape=19, position=position_jitter(0.2), size=0.7) + 
                xlab('Expansion category') +  ylab('% under forest cover') + 
                theme_classic()  + 
                scale_fill_brewer(palette=4) +
                theme(legend.position = 'none') + 
                theme(axis.text.x = element_text(angle = 40, hjust = 1)) + 
                scale_x_discrete(labels=c('Chelsa exp, \nSoilTemp no exp', 'Chelsa no exp, \nSoilTemp exp', 'No change \n(no exp)', 'No change \n(exp')) + 
                theme(text = element_text(family = "Helvetica", size=10)) + 
                theme(axis.text.x=element_text(colour="black")) 
kruskal.test(test_data$forestcovDIFF, test_data$expcategories)

# jaccard 
mod <- lm(Jdifference~forestcov, data=test_data)
hist(residuals.lm(mod)) # should i not be taking the residuals
simulationOutput <- simulateResiduals(fittedModel = mod, plot = F)
plot(simulationOutput)
summary(mod) #significant very weak correlation x
residuals.lm(mod) %>% as.vector
f_j <-  ggscatter( test_data, x = 'forestcov', y = 'Jdifference',# label='Species', 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = F, cor.method = "pearson", size=0.3, cor.coef.size=1, ## correlation coefficient is Kendalls correlation Tau because data is non-normal
                   #  cor.coef.coord = c(0, 1.1),
                   repel = TRUE, font.label=c(10, 'plain', 'black'), 
                   font.family='Helvetica') +
                    stat_cor(method='pearson', 
                             aes(label = paste(..r.label.., ..rr.label.., ..p.label.., sep = "~`,`~")), 
                             label.x = 3, size=3) + 
                    theme_minimal()+
                    theme(text = element_text(family = "Helvetica", size=10)) +
                    xlab('% under forest cover') +  ylab("Jaccard similarity difference") 
cor(test_data$forestcovDIFF, test_data$Jdifference, method='pearson')

fig5<- ggarrange(d,  exp, jac, f_d,   f_e,   f_j, ncol=3, nrow=2,  labels=c("A", "B","C", "D","E", "F"), common.legend = T)
ggsave(fig5, filename = "Graphs/fig5.pdf", width=7, height=6, units = "in", device='pdf', bg='white')

## check sample size 
test_data$nativeSS <-  sample_size$V1
test_data$invSS <-  sample_size$V2
test_data$Samplesize <- test_data$nativeSS+ test_data$invSS 

cor.test(log(test_data$Samplesize),   test_data$Ddifference, method='kendall')
cor.test(log(test_data$Samplesize),   test_data$BDdifference, method='kendall')
test_data$Exp_difference
cor.test(log(test_data$Samplesize),   test_data$Jdifference, method='kendall')
cor.test(log(test_data$Samplesize),   test_data$disCdifference, method='kendall')


