##%######################################################%##
#                                                          #
####                Data Analysis                       ####
#                                                          #
##%######################################################%##

{
  require(jtools)
  require(terra)
  require(dplyr)
  require(flexsdm)
  require(here)
  require(ggplot2)
  require(viridis)
  require(tidyr)
  require(landscapemetrics)
  require(raster)
  require(landscapetools)
  require(data.table)
  library(BAMMtools)
  library(GGally)
  library(ggcorrplot)
  require(gridExtra)
  require(gam)
  library(sjPlot)
  library(sjmisc)
  library(sjlabelled)
  library(AER)
  library(plyr)
  require(patchwork)
  library(glmm)
  library(ggeffects)
  library(lme4)
  library(nlme)
  library(arm)
  library(gamlss.mx)
  library(tidyverse)
  library(stats)
  library(vegan)
  library(adespatial)
  library(rpart)
  library(rpart.plot)
}

##%######################################################%##
#                                                          #
####                      Data                          ####
#                                                          #
##%######################################################%##

# merging all data

# traits data
sp_df <-
  data.table::fread('sp_traits_data.csv') %>% as_tibble()

# exposure data
exp_df <-
  data.table::fread('exposure_data.csv') %>% as_tibble()

# Data set for modeling with values averaged by GCM
m_df <-
  exp_df %>% dplyr::group_by(sp, rcp, change_driver, year) %>% dplyr::summarise_all(funs(mean)) %>%
  dplyr::select(-model)

m_data <- left_join(sp_df, m_df, by = c('species' = 'sp')) %>% 
  dplyr::select(-scen)

# converting patch isolation and number of patches to log scale
m_data$log10_patch_isolation <- log10(m_data$patch_isolation)
m_data$log10_num_patches <- log10(m_data$num_patches)

##%######################################################%##
#                                                          #
####                       MODELS                       ####
#                                                          #
##%######################################################%##

# predictors
Vars <-
  c(
    "extent",
    "occupied_area",
    "abundance",
    "niche_breadth",
    "log10_num_patches",
    "log10_patch_isolation",
    "mean_topo",
    "mean_elev",
    "mean_d_coast"
  )
var_names <-
  c(
    expression('Extent' ~ (km ^ 2)),
    expression('Occupied area' ~ (km ^ 2)),
    "Abundance",
    "Niche Breadth",
    "log10(number of patches)",
    "log10(patch isolation)",
    "Topographic heterogeneity",
    "Elevation (m)",
    "Distance to coast (km)"
  )

# models arguments
Vars_m <- as.list(
  c(
    "rcp*poly(extent,2)",
    "rcp*occupied_area",
    "rcp*poly(abundance,2)",
    "rcp*niche_breadth",
    "rcp*log10_num_patches",
    "rcp*log10_patch_isolation",
    "rcp*mean_topo",
    "rcp*mean_elev",
    "rcp*mean_d_coast"
  )
)

# Climate change exposure
cc_data <- m_data %>% filter(change_driver == 'climate_change')

# Land use change exposure
luc_data <- m_data %>% filter(change_driver == 'lulc')
luc_data$exposure[luc_data$exposure <= 0] <- 0

# Climate and land use change exposure
comb_data <-
  m_data %>% filter(change_driver == 'climate_change + lulc')

# Full model formula (for selecting a distribution family)
# We selected extent as the measure of range size for our final models because it explains
# more variation in exposure than occupied area in exploratory analyses
fmula <-
  formula(
    exposure ~
      rcp +
      extent +
      rcp * extent +
      poly(abundance, 2) +
      niche_breadth +
      rcp * niche_breadth +
      log10_patch_isolation +
      log10_num_patches +
      rcp * log10_num_patches +
      mean_topo +
      rcp * mean_topo +
      mean_elev +
      re(random =  ~ 1 | species)
  )

# Selecting distribution with GAMLSS
library(gamlss)

m0 <- fitDist(cc_data$exposure, type = "realline") 

# histogram of response variable and chosen distribution family (SEP2)
fh <-
  histDist(
    cc_data$exposure,
    family = NO,
    nbins = 30,
    line.col = "black"
  )
fh0 <-
  histDist(
    cc_data$exposure,
    family = SN2,
    nbins = 30,
    line.col = "black"
  )
fh <-
  histDist(
    cc_data$exposure,
    family = SN1,
    nbins = 30,
    line.col = "black"
  )
fh2 <-
  histDist(
    cc_data$exposure,
    family = SHASHo2,
    nbins = 30,
    line.col = "black"
  )

# SN1 distribution family
m1 <-
  gamlss(fmula,
         data = cc_data,
         family = SN1)

plot(m1)
wp(m1) # worm plot
term.plot(m1, pages = 1)

# SN1 distribution family
m2 <-
  gamlss(
    fmula,
    data = cc_data,
    family = SN1
  )
plot(m2)
wp(m2) # worm plot
term.plot(m2, pages = 1)

# SHASHo2 family
m3 <-
  gamlss(
    fmula,
    data = cc_data,
    family = SHASHo2
  )
plot(m3)
wp(m3) # worm plot
term.plot(m3, pages = 1)

# NO family
m4 <-
  gamlss(
    fmula,
    data = cc_data,
    family = NO
  )
plot(m4)
wp(m4) # worm plot
term.plot(m4, pages = 1)

AIC(m1, m2, m3, m4)
# select m1 (family SN1)

# We selected the SN1 distribution family
# because of a relatively low AIC when compared to other distribution
# families. We also checked the model residuals which were
# better than those produced by models built with the other
# distributions. Although the SHASHO2 model had a lower AIC, the residuals were
# worse than the model produced with SN1 distribution family.
# Adding nu and sigma parameters improved model residuals when modeling
# ~exposure as a function of each species' trait

allModelsList <-
  lapply(paste("exposure ~", Vars_m, '+ re(random = ~ 1 |
                                                     species)'),
         as.formula)

nu <-
  lapply(paste("~rcp +", Vars),
         as.formula)
mu <-
  lapply(paste("~rcp +", Vars),
         as.formula)
sigma <-
  lapply(paste("~", Vars),
         as.formula)

allModelsResults <- list()

for (i in 1:length(allModelsList)) {
  allModelsResults[[i]] <-
    gamlss(
      allModelsList[[i]],
      data = cc_data,
      family = SN1,
      nu.fo = nu[[i]],
      sigma.fo = sigma[[i]]
    )
}

allModelTibble <- lapply(allModelsResults, function(x)
  as_tibble(unlist(summary(x)), rownames = 'term'))

# Calculating R-squared
allModelsList11 <-
  lapply(paste("exposure ~", Vars_m),
         as.formula)

nu <-
  lapply(paste("~rcp +", Vars),
         as.formula)
mu <-
  lapply(paste("~rcp +", Vars),
         as.formula)
sigma <-
  lapply(paste("~", Vars),
         as.formula)

allModelsResults11 <- list()

for (i in 1:length(allModelsList11)) {
  allModelsResults11[[i]] <-
    gamlss(
      allModelsList11[[i]],
      data = cc_data,
      family = SN1,
      nu.fo = nu[[i]],
      # mu.fo = mu[[i]],
      sigma.fo = sigma[[i]]
    )
}

lapply(allModelsResults11, Rsq)


allModelsList12 <-
  lapply(paste("exposure ~", Vars),
         as.formula)

nu <-
  lapply(paste("~", Vars),
         as.formula)
sigma <-
  lapply(paste("~", Vars),
         as.formula)

allModelsResults12 <- list()

for (i in 1:length(allModelsList12)) {
  allModelsResults12[[i]] <-
    gamlss(
      allModelsList12[[i]],
      data = cc_data,
      family = SN1,
      nu.fo = nu[[i]],
      # mu.fo = mu[[i]],
      sigma.fo = sigma[[i]]
    )
}

lapply(allModelsResults12, Rsq)


##%######################################################%##
#                                                          #
####               Land use change model                ####
#                                                          #
##%######################################################%##

luc_data <- m_data %>% filter(change_driver == 'lulc')
luc_data$exposure[luc_data$exposure <= 0.00] <- 0

# ZAGA
m1 <-
  gamlss(fmula,
         data = luc_data,
         family = ZAGA)
plot(m1)
wp(m1) # worm plot
term.plot(m1, pages = 1)

# BEINF0
m2 <-
  gamlss(
    fmula,
    data = luc_data,
    family = BEINF0
  )
plot(m2)
wp(m2) # worm plot
term.plot(m2, pages = 1)

AIC(m1, m2)
# Selected BEINF0 because of lower AIC and better residuals
# Adding nu, mu, and sigma parameters improved model residuals when modeling
# ~exposure as a function of each species' trait

# Loop for creating models for ~ exposure and each species trait
allModelsList2 <-
  lapply(paste("exposure ~", Vars_m, '+ re(random = ~ 1 |
                                                     species)'),
         as.formula)

allModelsResults2 <- list()

for (i in 1:length(allModelsList2)) {
  allModelsResults2[[i]] <-
    gamlss(
      allModelsList2[[i]],
      data = luc_data,
      family = BEINF0,
      nu.fo = nu[[i]],
      mu.fo = mu[[i]],
      sigma.fo = sigma[[i]]
    )
}

allModelTibble2 <- lapply(allModelsResults2, function(x)
  as_tibble(unlist(summary(x)), rownames = 'term'))


# Calculating R-squared
allModelsList21 <-
  lapply(paste("exposure ~", Vars_m),
         as.formula)

nu <-
  lapply(paste("~rcp +", Vars),
         as.formula)
mu <-
  lapply(paste("~rcp +", Vars),
         as.formula)
sigma <-
  lapply(paste("~", Vars),
         as.formula)

allModelsResults21 <- list()

for (i in 1:length(allModelsList21)) {
  allModelsResults21[[i]] <-
    gamlss(
      allModelsList21[[i]],
      data = luc_data,
      family = BEINF0,
      nu.fo = nu[[i]],
      # mu.fo = mu[[i]],
      sigma.fo = sigma[[i]]
    )
}

lapply(allModelsResults21, Rsq)

##%######################################################%##
#                                                          #
####         Land use change and climate change         ####
#                                                          #
##%######################################################%##

dist3 <- fitDist(comb_data$exposure, type = "realAll")
dist3$fits

fh2 <-
  histDist(
    comb_data$exposure,
    family = SN2,
    nbins = 30,
    line.col = "black"
  )
fh2 <-
  histDist(
    comb_data$exposure,
    family = SN1,
    nbins = 30,
    line.col = "black"
  )
fh2 <-
  histDist(
    comb_data$exposure,
    family = SEP2,
    nbins = 30,
    line.col = "black"
  )
fh2 <-
  histDist(
    comb_data$exposure,
    family = NO,
    nbins = 30,
    line.col = "black"
  )

# SN2
m1 <-
  gamlss(
    fmula,
    data = comb_data,
    family = SN2
  ) 

plot(m1)
wp(m1)
term.plot(m1, pages = 1)

# SN1
m2 <-
  gamlss(fmula,
         data = comb_data,
         family = SN1) 

plot(m2)
wp(m2)
term.plot(m2, pages = 1)

m3 <-
  gamlss(fmula,
         data = comb_data,
         family = NO) 

plot(m3)
wp(m3)
term.plot(m3, pages = 1)

AIC(m1, m2, m3) # select SN1

# overall, adding shape and location parameters improved model performance (AIC and residuals)

allModelsList3 <-
  lapply(paste("exposure ~", Vars_m, '+ re(random = ~ 1 |
                                                     species)'),
         as.formula)

allModelsResults3 <- list()

for (i in 1:length(allModelsList3)) {
  allModelsResults3[[i]] <-
    gamlss(
      allModelsList3[[i]],
      data = comb_data,
      family = SN1,
      sigma.fo = sigma[[i]]
    )
}

allModelTibble3 <- lapply(allModelsResults3, function(x)
  as_tibble(unlist(summary(x)), rownames = 'term'))

# all model estimates
allModels3 <- bind_rows(allModelTibble3)

# Calculating R-squared
allModelsList31 <-
  lapply(paste("exposure ~", Vars_m),
         as.formula)

nu <-
  lapply(paste("~rcp +", Vars),
         as.formula)
mu <-
  lapply(paste("~rcp +", Vars),
         as.formula)
sigma <-
  lapply(paste("~", Vars),
         as.formula)

allModelsResults31 <- list()

for (i in 1:length(allModelsList31)) {
  allModelsResults31[[i]] <-
    gamlss(
      allModelsList31[[i]],
      data = comb_data,
      family = SN1,
      nu.fo = nu[[i]],
      sigma.fo = sigma[[i]]
    )
}

lapply(allModelsResults31, Rsq)


##%######################################################%##
#                                                          #
####    Phylogenetic dependency tests                   ####
#                                                          #
##%######################################################%##

# Following methods performed by Estrada et al. 2015, we compared
# Moran's I phylogenetic correlograms for the response variable (exposure)
# and the residuals of each gamlss model. "This approach determines whether
# phylogenetic autocorrelation in a response variable has been captured by the
# model predictions. We also compared the significance of model coefficients
# of the gamlss models and from phylogenetic generalized least squares (PGLS)

# phylogenetic tree for California plants:
# Downloaded from https://datadryad.org/stash/dataset/doi:10.6078/D1VD4P
# Thornhill, A. H., B. G. Baldwin, W. A. Freyman, and S. Nosratinia. 2017a. Sequence matrix and tree files for Spatial phylogenetics of the native California flora.

# OTU/tip label data downloaded from: https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-017-0435-x
# Thornhill, A. H., B. G. Baldwin, W. A. Freyman, S. Nosratinia, M. M. Kling, N. Morueta-Holme, T. P. Madsen, D. D. Ackerly, and B. D. Mishler. 2017b. Spatial phylogenetics of the native California flora. BMC biology 15:96.

library(castor)
library(phylosignal)
library(phylobase)

# character strings for loops
var_names2 <-
  c(
    'Extent',
    'Abundance',
    'Niche breadth',
    'Number of patches (log10)',
    'Patch isolation (log10)',
    'Topographic heterogeneity',
    'Elevation',
    'Distance to coast'
  )

var_res <-
  c(
    'ext_res',
    'abun_res',
    'niche_res',
    'nump_res',
    'iso_res',
    'topo_res',
    'elev_res',
    'dc_res'
  )
var_names3 <-
  c('extent',
    'abun',
    'niche',
    'np',
    'patch_iso',
    'topo',
    'elev',
    'dis_c')

# Model data with gamlss residuals and operational taxon unit tip labels (from Thornhill phylgeny)
phy_mod_data <- data.table::fread('model_data_phylogenetic_residuals.csv') %>% tibble()

# subset the tree based on species we need
tree <- ape::read.nexus("Californian_clades_tree_final.nex") # Can be found here: https://datadryad.org/stash/dataset/doi:10.6078/D1VD4P
subtree <-
  get_subtree_with_tips(tree, only_tips = m_data_phy$OTU)$subtree

# list of change drivers & rcps for loops
drivers <- c('climate_change', 'lulc','climate_change + lulc')
rcps <- c('rcp45', 'rcp85')

pc <- list()
plot_param <- list()


# quick method for producing phylogenetic correlograms
for (k in 1:length(drivers)) {
  for (j in 1:length(rcps)) {
    df <-
      phylo4d(subtree,
              tip.data = phy_mod_data %>% filter(change_driver == drivers[k] &
                                                   rcp == rcps[j]))
    pc <- list() # list for correlogram plots
    
    for (i in 1:length(var_res)) {
      pc[[i]] <- phyloCorrelogram(df, trait = var_res[i])
      # insert saving information here
      plot(pc[[i]])
    }
  }
}

#loop for saving phylogenetic correlograms (can be inserted in loop above)
# for (i in 1:length(var_res)) {
#   print(i)
#   
#   pc[[i]] <- phyloCorrelogram(df, trait = var_res[i])
#   
#   png(
#     paste0(
#       "filename",
#       var_names3[i],
#       ".png"
#     ),
#     height = 12,
#     width = 12,
#     units = 'cm',
#     res = 500
#   )
#   
#   plot(
#     pc[[i]],
#     main = "",
#     xlab = "",
#     ylab = "",
#     axes = FALSE
#   )
#   axis(
#     side = 2,
#     at = seq(-0.4, 0.8, 0.2),
#     pos = 0,
#     tck = 0.02
#   )
#   axis(side = 1, tck = 0.02)
#   mtext("Phylogenetic Distance", side = 1, line = 2.5)
#   mtext("Correlation", side = 2, line = 2.5)
#   mtext(paste0(var_names2[i]), side = 3, line = 1)
#   
#   dev.off()
#   
# }


##%######################################################%##
#                                                          #
####                   Decision trees                   ####
#                                                          #
##%######################################################%##

dt_fmula <-
  formula(
    exposure ~
      extent +
      abundance +
      niche_breadth +
      mean_topo +
      mean_elev +
      mean_d_coast +
      num_patches +
      patch_isolation
  )

# ##%######################################################%##
#                                                          #
####                   Climate change                   ####
#                                                          #
##%######################################################%##

# RCP 4.5
cc_data45 <- cc_data %>% filter(rcp == 'rcp45')
tree_m1 <-
  rpart(dt_fmula,
        data = cc_data45,
        method = 'anova')

rpart.plot(
  tree_m1,
  box.palette = "BuRd",
  shadow.col = "gray",
  nn = T,
  extra = 1
)

# evaluate the tree
printcp(tree_m1)	# display cp table
plotcp(tree_m1) # cross-validation results
rsq.rpart(tree_m1)

# variable importance plot for full model
df <- data.frame(imp = tree_m1$variable.importance)

# RCP 8.5
cc_data85 <- cc_data %>% filter(rcp == 'rcp85')
tree_m2 <- rpart(dt_fmula, data = cc_data85, method = 'anova')

rpart.plot(
  tree_m2,
  box.palette = "BuRd",
  shadow.col = "gray",
  nn = T,
  extra = 1
)

# variable importance plot for full model
df <- data.frame(imp = tree_m2$variable.importance)

# evaluate the tree
printcp(tree_m2)	# display cp table
plotcp(tree_m2) # cross-validation results
rsq.rpart(tree_m2)


##%######################################################%##
#                                                          #
####                  Land use change                   ####
#                                                          #
##%######################################################%##

# RCP 4.5
luc_data45 <- luc_data %>% filter(rcp == 'rcp45')

tree_m3 <- rpart(dt_fmula, data = luc_data45, method = 'anova')

rpart.plot(
  tree_m3,
  box.palette = "BuRd",
  shadow.col = "gray",
  nn = T,
  extra = 1
)

# variable importance plot for full model
df <- data.frame(imp = tree_m3$variable.importance)

# evaluate the tree
printcp(tree_m3)	# display cp table
plotcp(tree_m3) # cross-validation results
rsq.rpart(tree_m3)

# RCP 8.5
luc_data85 <- luc_data %>% filter(rcp == 'rcp85')
tree_m4 <- rpart(dt_fmula, data = luc_data85, method = 'anova')

rpart.plot(
  tree_m4,
  box.palette = "BuRd",
  shadow.col = "gray",
  nn = T,
  extra = 1
)

# variable importance plot for full model
df <- data.frame(imp = tree_m4$variable.importance)

# evaluate the tree
printcp(tree_m4)	# display cp table
plotcp(tree_m4) # cross-validation results
rsq.rpart(tree_m4)


##%######################################################%##
#                                                          #
####      Climate and Land use change                   ####
#                                                          #
##%######################################################%##

# not included in manuscript (very similar to tree for climate change only)

# RCP 4.5
comb_data45 <- comb_data %>% filter(rcp == 'rcp45')

tree_m5 <- rpart(dt_fmula, data = comb_data45, method = 'anova')

rpart.plot(
  tree_m5,
  box.palette = "BuRd",
  shadow.col = "gray",
  nn = T,
  extra = 1
)

# variable importance plot for full model
df <- data.frame(imp = tree_m5$variable.importance)

# evaluate the tree
printcp(tree_m5)	# display cp table
plotcp(tree_m5) # cross-validation results
rsq.rpart(tree_m5)

# RCP 8.5
comb_data85 <- comb_data %>% filter(rcp == 'rcp85')
tree_m6 <- rpart(dt_fmula, data = comb_data85, method = 'anova')

rpart.plot(
  tree_m6,
  box.palette = "BuRd",
  shadow.col = "gray",
  nn = T,
  extra = 1
)

# variable importance plot for full model
df <- data.frame(imp = tree_m6$variable.importance)

# evaluate the tree
printcp(tree_m6)	# display cp table
plotcp(tree_m6) # cross-validation results
rsq.rpart(tree_m6)


##%######################################################%##
#                                                          #
####               Variance partitioning                ####
#                                                          #
##%######################################################%##

# adding polynomial of extent and abundance for variance paritioning
m_data$extent_poly <- poly(m_data$extent, 2)
m_data$abundance_poly <- poly(m_data$abundance, 2)


cc_data <- m_data %>% filter(change_driver == 'climate_change')
luc_data <- m_data %>% filter(change_driver == 'lulc')
comb_data <- m_data %>% filter(change_driver == 'climate_change + lulc')

# Rarity traits
rarity <-
  c('extent_poly',
    'niche_breadth',
    'abundance_poly',
    'patch_isolation',
    'num_patches')

# Geographic traits
geo <- c('mean_topo', 'mean_elev', 'mean_d_coast')

# RCP
rcp <- c('rcp')

# Climate change
mod <-
  varpart(Y = cc_data$exposure, cc_data[rarity], cc_data[geo], cc_data[rcp])

plot(
  mod,
  Xnames = c("Rarity traits", "Geographic traits", "RCP"),
  # name the partitions
  bg = c("#440154", "#30678D", "#B3DC2B"),
  alpha = 80,
  # colour the circles
  digits = 2,
  # only show 2 digits
  cex = 1.5
)


# Land use change
mod2 <-
  varpart(Y = luc_data$exposure, luc_data[rarity], luc_data[geo], luc_data[rcp])

plot(
  mod2,
  Xnames = c("Rarity traits", "Geographic traits", "RCP"),
  # name the partitions
  bg = c("#440154", "#30678D", "#B3DC2B"),
  alpha = 80,
  # color the circles
  digits = 2,
  # only show 2 digits
  cex = 1.5
)


