##%######################################################%##
#                                                          #
####        Calculating species' traits                 ####
#                                                          #
##%######################################################%##

{
  require(dplyr)
  require(terra)
  require(ggplot2)
  require(patchwork)
  require(here)
  require(sf)
  require(readxl)
  require(dplyr)
  require(janitor)
  require(readr)
  require(stringr)
  require(textclean)
  require(raster)
  require(tidyverse)
  require(fuzzySim)
  library(adehabitatHR)
  library(ecospat)
  library(BAMMtools)
  library(GGally)
  require(flexsdm)
  require(landscapemetrics)
}

##%######################################################%##
#                                                          #
####        Species and environmental data              ####
#                                                          #
##%######################################################%##

# California floristic provinces
cfp <- terra::vect("./1_Inputs/4_Shapefiles/JepsonRegions/JepsonRegions_CFP.shp")

# Read occurrences databases
occ <- data.table::fread('1_Inputs/1_Occurrences/3_OccPartition/occ_part_block.gz') %>% tibble()
occ2 <-  data.table::fread('1_Inputs/1_Occurrences/3_OccPartition/occ_part_band.gz') %>% tibble()

occ_df <- bind_rows(occ, occ2) %>%
  filter(species != 'Mimulus cardinalis' &
           species != 'Quercus dumosa' &
           species != 'Arctostaphylos mewukka')

# presence only occurrences
occ3 <- data.table::fread('1_Inputs/1_Occurrences/4_PresenceOnly/2_occ_presabs_part_allsp.gz') %>% tibble()
occ3["species"][occ3["species"] == "Calochortus obispnsis"] <- "Calochortus obispoensis"

# all occurrence records
occ_final <- bind_rows(occ_df, occ3)

# Environmental variables - Current conditions
env <- file.path('1_Inputs/2_Predictors/1_Current') %>%
  list.files(., full.names = T, pattern = '.tif$') %>%
  terra::rast()
env <- flexsdm::homogenize_na(env)
env_names <- names(env)

# List of data frames with each species p/a
sp_list <- split(occ_final, occ_final$species)

p_points <- list()

for(i in 1:length(sp_list)){
  p_points[[i]] <- sp_list[[i]] %>%
    dplyr::filter(pr_ab == 1) %>%
    terra::vect(crs = crs(env),
                geom = c('x', 'y'))
}

names(p_points) <- names(sp_list)


##%######################################################%##
#                                                          #
####                     Abundance                      ####
#                                                          #
##%######################################################%##

# Dmax =The largest distance between any pair of occurrence points of each species
# Dmin = Average of minimum distances between occurrence points for each point
# Local abundance = Dmax/Dmin

abundance <- list()

for (i in 1:length(p_points)) {
  dist <- terra::distance(p_points[[i]])
  distm <- as.matrix(dist)
  diag(distm) <- NA
  dmax <- distm[which.max(distm)]
  dmin <- apply(distm, 1, min, na.rm = TRUE)
  mdmin <- mean(dmin)
  abundance[[i]] <-
    data.frame(dmax,
               mdmin,
               abundance = dmax / mdmin,
               species = names(p_points[i]))
}

abundance <- bind_rows(abundance)

##%######################################################%##
#                                                          #
####                    Range Size                      ####
#                                                          #
##%######################################################%##

range_size100 <- list()
range_size95 <- list()

for (i in 1:length(p_points)) {
  print(names(p_points[i]))
  range_size100[[i]] <-
    mcp(
      xy = as(p_points[[i]], "Spatial"),
      percent = 100,
      unin = c("m"),
      unout = c("km2")
    )
  range_size95[[i]] <-
    mcp(
      xy = as(p_points[[i]], "Spatial"),
      percent = 95,
      unin = c("m"),
      unout = c("km2")
    )
  range_size100[[i]] <-
    as_tibble(range_size100[[i]]) %>%
    dplyr::mutate(species = names(p_points[i])) %>%
    dplyr::rename(area100 = area)
  range_size95[[i]] <-
    as_tibble(range_size95[[i]]) %>%
    dplyr::mutate(species = names(p_points[i])) %>%
    dplyr::rename(area95 = area)
}

range_size100 <- bind_rows(range_size100) %>% dplyr::select(-id)
range_size95 <- bind_rows(range_size95) %>% dplyr::select(-id)
range_size <- left_join(range_size100, range_size95, by = 'species')

df <- left_join(abundance, range_size, by = 'species')

df <- df[order(df$area100),]

hist(df$area100, xlab = "Range size", main = "Freabundquency distribution of range size")


##%######################################################%##
#                                                          #
####                  Habitat Specificity               ####
#                                                          #
##%######################################################%##


# niche breadth calculation (Diaz et al. 2020)
# 1) standardize environmental variables to z scores (mean = 0, SD= 1) to control for differences in variation
# and different measurement units and extract those values at spp. locations

p_points <- lapply(p_points, as.data.frame)

env_points <- sdm_extract(
  data = bind_rows(p_points),
  x = 'x',
  y = 'y',
  env_layer = env[[-9]],
  filter_na = TRUE
)

env_points <- split(env_points, env_points$species)


# standardized environmental variables
z_rast <- terra::scale(env[[-9]])

# Principal components analysis of environmental variables at occurrence locations
# climate variables principal components that sum up to 95% of the variation of the climatic variables (3 selected)
clim_z <- subset(z_rast, c(1,3,7,8,9))
clim_pca <- correct_colinvar(clim_z, method = 'pca')

# ordinal PC plot
db2 <- clim_pca$coefficients %>%
  dplyr::select(variable:PC3) %>%
  mutate(x0 = 0, y0 = 0)


clim_pca_plot <- ggplot(db2, aes(PC1, PC2, label = variable)) +
  geom_hline(yintercept = 0, col = 'gray50') +
  geom_vline(xintercept = 0, col = 'gray50') +
  geom_segment(
    aes(
      x = x0,
      y = y0,
      xend = PC1,
      yend = PC2
    ),
    arrow = arrow(type = "closed", angle = 15, unit(0.30, "cm")),
    col = 'blue'
  ) +
  ggrepel::geom_text_repel() +
  theme_minimal() +
  theme(text = element_text(size = 25))

ggsave(
  clim_pca_plot,
  filename = "G:/My Drive/Dissertation/002_Biogeography_exposure/figures/Sup_climate_PCA.png",
  dpi = 500,
  height = 25,
  width = 25,
  unit = 'cm'
)


clim_pca_r <- clim_pca$env_layer
names(clim_pca_r) <- c('PC1_clim', 'PC2_clim', 'PC3_clim')

# PCA of soil variables (4 selected)
soil_z <- subset(z_rast, c(2,4,5,6))
soil_pca <- correct_colinvar(soil_z, method = 'pca')

# ordinal PC plot
db2 <- soil_pca$coefficients %>%
  dplyr::select(variable:PC4) %>%
  mutate(x0 = 0, y0 = 0)

soil_pca_plot <- ggplot(db2, aes(PC1, PC2, label = variable)) +
  geom_hline(yintercept = 0, col = 'gray50') +
  geom_vline(xintercept = 0, col = 'gray50') +
  geom_segment(
    aes(
      x = x0,
      y = y0,
      xend = PC1,
      yend = PC2
    ),
    arrow = arrow(type = "closed", angle = 15, unit(0.30, "cm")),
    col = 'blue'
  ) +
  ggrepel::geom_text_repel() +
  theme_minimal() +
  theme(text = element_text(size = 25))

ggsave(
  soil_pca_plot,
  filename = "G:/My Drive/Dissertation/002_Biogeography_exposure/figures/Sup_soil_PCA.png",
  dpi = 500,
  height = 25,
  width = 25,
  unit = 'cm'
)


soil_pca_r <- soil_pca$env_layer
names(soil_pca_r) <- c('PC1_soil', 'PC2_soil', 'PC3_soil', 'PC4_soil')


# combining pca components into raster
env_pca <- c(clim_pca_r, soil_pca_r)

z_points <- sdm_extract(
  data = bind_rows(p_points),
  x = 'x',
  y = 'y',
  env_layer = env_pca,
  filter_na = TRUE
)

z_points <- split(z_points, z_points$species)

# PCA plot of niche breadth components


#2) mean value of each environmental variable across all plots occupied by a
# species as the sum of values for each individual plant divided by the total
# number of individuals of that species

env_var <- names(env_pca)

xser <- list()
ssqd <- list()
nb <- list()

for (i in 1:length(z_points)) {
  for (k in 1:length(env_var)) {
    xser[k] <-
      sum(z_points[[i]][paste(env_var[k])]) / nrow(z_points[[i]])
    
    x <- c(1:nrow(z_points[[i]]))
    
    ssqd[[k]] <-
      (z_points[[i]][paste(env_var[k])][c(x), ] - xser[k]) ^ 2
  }
  
  ssqd <- sum(ssqd %>% unlist())
  
  nb[[i]] <-
    data.frame(
      niche_breadth = ssqd/ nrow(z_points[[i]]),
      species = names(z_points[i]),
      pr = nrow(env_points[[i]]),
      tmn_range = diff(range(env_points[[i]]$tmn)),
      ppt_djf_range = diff(range(env_points[[i]]$ppt_djf)),
      ppt_jja_range = diff(range(env_points[[i]]$ppt_jja)),
      clay_range = diff(range(env_points[[i]]$pct_clay)),
      cwd_range = diff(range(env_points[[i]]$cwd)),
      aet_range = diff(range(env_points[[i]]$aet)),
      awc_range = diff(range(env_points[[i]]$awc)),
      depth_range = diff(range(env_points[[i]]$depth)),
      ph_range = diff(range(env_points[[i]]$ph))
    )
  
}

niche_breadth <- bind_rows(nb)
niche_breadth <- niche_breadth[order(niche_breadth$niche_breadth),]

hist(niche_breadth$niche_breadth, xlab = "Niche breadth", main = "Frequency distribution of niche breadth")


##%######################################################%##
#                                                          #
####      Calculating fragmentaiton traits              ####
#                                                          #
##%######################################################%##

# number of patches and patch isolation
# function with 1) m1 = habitat map, 2) threshold
# 3) species name and 4) model type

frag_traits <- function(m1, thr, species, model) {
  
  # binary map
  b1 <- m1
  b1[b1 >= thr] <- 1
  b1[b1 <= thr] <- 0

  # Current range traits
  
  c.landscape <-
    rbind(lsm_c_np(b1),
          # Number of patches
          lsm_c_ca(b1),
          # Total area
          lsm_c_enn_mn(b1)) %>%
          # Mean of euclidean nearest-neighbor distance (GISFrag)) %>%
          filter(class == 1) %>% # selecting only part of the landscape that is suitable
            dplyr::select(metric, value, class) %>%
            pivot_wider(names_from = metric,
                        values_from = value,
                        names_glue = "c_{metric}")
          
          range_traits <- tibble(species,
                                 model,
                                 c.landscape,
                                 occupied_area = c_ca/100)
          
          return(range_traits)
}


# thresholds for SDMs
thr <-
  readr::read_tsv("./2_Outputs/0_Model_performance/00_model_performance.txt") %>%
  filter(species != 'Mimulus cardinalis')
thr2 <- readr::read_tsv("./2_Outputs/0_Model_performance/00_model_performance_pres_only.txt")
thr <- bind_rows(thr, thr2) # combining presence-absence and presence-only species

# LU masked maps for only measuring impact of climate change
lu_dir <- "./2_Outputs/9_Final_SDM/4_SDM_LULC_ONLY" %>% list.dirs(., recursive = FALSE)
names(lu_dir) <- basename(lu_dir)
sp <- names(lu_dir)

c_range_traits <- list()

for (i in 1:length(lu_dir)) {
  
  print(sp[i])
  
  dir_save2 <- file.path("./2_Outputs/9_Final_SDM/2_SDM_LU", sp[i])
  dir_save2 <- dir_save2 %>% list.dirs(recursive = FALSE)
  files <- lapply(dir_save2, function(x) list.files(x, pattern = ".asc$", full.names = TRUE))
  names(files) <- basename(dir_save2)
  
  # threshold for species meanthr model
  thr2 <- thr %>%
    dplyr::filter(model == 'meanthr' ,
                  species == sp[i]) %>%
    dplyr::select(thr_value)
  
  
  # Current condition
  lf <- files$`01_current`
  r <- terra::rast(lf)
  
  c_range_traits[[i]] <- frag_traits(m1 = r, thr = thr2$thr_value, sp = sp[i], model = 'meanthr') %>% dplyr::select(-c_ca)

  
}

c_range_traits2 <- bind_rows(c_range_traits)
data.table::fwrite(c_range_traits2, '2_Outputs/3_Exposure/species_frag_data_raw.csv')



##%######################################################%##
#                                                          #
####                  Rarity Framework                  ####
#                                                          #
##%######################################################%##

# We classified each species into  the rarity framework using the average value of  all
# species for each index: Dmax, Dmax/Dmin, and habitat specificity.

df2 <- left_join(df, niche_breadth, by = 'species')
df3 <- left_join(df2, niche_pos, by = 'species')
df3$species<-gsub(" ", "_", df3$species)

# Range size: Narrow vs. wide
# Habitat specificity: Broad vs. restricted
# Abundance: Large vs. small

# explore natural breaks
getJenksBreaks(df3$area100, k = 5)
getJenksBreaks(df3$niche_breadth, k = 5)
getJenksBreaks(df3$abundance, k = 5)

# Range size: Narrow vs. wide
# Habitat specificity: Broad vs. restricted
# Abundance: Large vs. small

rarity_df <-
  df3 %>% dplyr::mutate(
    range_size = ifelse(area100 > 20000, 'w', 'n'),
    habitat = ifelse(niche_breadth >= 4, 'b', 'r'),
    abun_meas = ifelse(abundance >= 200, 'l', 's')
  )

rarity_df$rabin_class <- apply( rarity_df[19:21] , 1 , paste , collapse = "-" )

data.table::fwrite(rarity_df, '2_Outputs/3_Exposure/species_rarity_data_raw.csv')

rarity_df <- data.table::fread('2_Outputs/3_Exposure/species_rarity_data_raw.csv') %>% as_tibble()


##%######################################################%##
#                                                          #
####              Geographic traits                     ####
#                                                          #
##%######################################################%##

# Topographic heterogeneity of each species range
topo <- terra::rast("G:/My Drive/Franklin_grant/project/data/Topo_hetero/topo_hetero_cfp.tif")

# Average distance to coast of each species range
d_coast <- terra::rast("G:/My Drive/Franklin_grant/project/data/CoastDistance/coast_distance_cfp.tif")

# elevation
elev <- terra::rast("G:/My Drive/Franklin_grant/project/data/Elevation/dem_cfp.tif")

# stack biogeographic factors
bio_traits <- c(topo, d_coast, elev)
names(bio_traits) <- c('topo_het', 'd_coast', 'elevation')
bio_traits <- c(bio_traits, env)

# species occurrence points
p_points <- lapply(p_points, as.data.frame)

bio_points <- sdm_extract(
  data = bind_rows(p_points),
  x = 'x',
  y = 'y',
  env_layer = bio_traits,
  filter_na = TRUE
)

bio_df <- bio_points %>%
  group_by(species) %>%
  dplyr::summarise(mean_topo = mean(topo_het),
                   mean_d_coast = mean(d_coast),
                   mean_elev = mean(elevation))

data.table::fwrite(bio_df, '2_Outputs/3_Exposure/occ_based_biogeographic_traits_data.csv')
bio_df <- data.table::fread('2_Outputs/3_Exposure/occ_based_biogeographic_traits_data.csv') %>% as_tibble() %>%
  dplyr::select(species, mean_topo, mean_d_coast, mean_elev)

# merging data
rarity_df$species <-gsub("_"," ", rarity_df$species)
full_traits_df <- left_join(rarity_df, bio_df, by = c('species'))
full_traits_df2 <- left_join(full_traits_df, c_range_traits2, by = c('species'))


data.table::fwrite(full_traits_df2, '2_Outputs/3_Exposure/sp_traits_data.csv')
full_traits_df2 <- data.table::fread('2_Outputs/3_Exposure/sp_traits_data.csv') %>% as_tibble()
