reticulate::use_condaenv("base")

library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)
# library(ggnewscale)
# library(ggtern)
# library(tidyverse)
#library(rgdal)
# library(lightgbm)
library(keras)
library(tensorflow)
library(tfprobability)
library(scoringRules)
library(scoringutils)
library(ranger)
library(gtools)
library(ggfortify)
library(compositions)
library(terra)
# library(tidyterra)
library(parallel)

options(scipen = 6000)

# Set tensorflow mixed precision:
mixed_precision <- keras$mixed_precision$Policy('mixed_float16')

keras$mixed_precision$set_global_policy(mixed_precision)

tfd_multivariate_normal_diag <- function(loc = NULL,
                                         scale_diag = NULL,
                                         validate_args = FALSE,
                                         allow_nan_stats = TRUE,
                                         name = "MultivariateNormalDiag") {
  args <- list(
    loc = loc,
    scale_diag = scale_diag,
    validate_args = validate_args,
    allow_nan_stats = allow_nan_stats,
    name = name
  )
  
  do.call(tfp$distributions$MultivariateNormalDiag, args)
}

# Load elevation data:
time <- Sys.time()
elev <- rast("data/British Isles 250m Copernicus DEM.tif")
# bielev <- aggregate(elev, fact = 2)
bielev <- elev
names(bielev) <- "alt"
elevdat <- as.data.table(terra::as.data.frame(bielev, xy = TRUE, na.rm = TRUE))
elevdat
Sys.time() - time

# Plot of UK elevation data
ggplot(elevdat) + geom_raster(aes(x = x, y = y, fill = alt)) + theme_bw() +
  scale_fill_viridis_c(option = "cividis", begin = 0.2, end = 1) + coord_equal() +
  theme(legend.position = "none") + labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)") +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks = seq(0,1000000, length.out = 6), expand = c(0,0), limits = c(0, 1225000)) +
  scale_x_continuous(breaks = seq(0,600000, length.out = 4), expand = c(0,0), limits = c(-200000, 664000))
ggsave(paste0("paperplots/UKelevation.png"), width = 144, height = 200, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

# Load geochemistry data ####

data <- fread("data/streamsedimentgeochemistryUKandIE.csv")

# Extract sample elevation from elevation data

time <- Sys.time()
data[, elevation := extract(bielev, methods = "bilinear", data[, c("Easting_BNG", "Northing_BNG"), with = FALSE])$alt, ]
data <- data[!is.na(data$elevation)]
Sys.time() - time

#data <- data[sample(1:nrow(data), 4000), ]

# Choose elements to map as clr compositon rgb ####

compelem <- c("K2O", "Fe2O3", "CaO")

R <- 1000000 - rowSums(data[, ..compelem])
# comp <- data.table(clr(acomp(cbind(data[, ..compelem], R))))
comp <- data.table(clr(acomp(data[, ..compelem])))

mapdat <- cbind(data[, c("Easting_BNG", "Northing_BNG", "elevation")], comp)
mapdat <- mapdat[-which(K2O == 0 | Fe2O3 == 0 | CaO == 0 | R == 0), ]

mapdat[mapdat == 0] <- NA
mapdat <- na.omit(mapdat)

scurve <- function(x, contrast){
  # return(1/(1+(x/(1-x))^-contrast)) this was a bug - the curve went the wrong way, reducing contrast!
  return(1/(1+exp(-((x*2)-1)/(1-contrast))))
}

# mapdat_rgb <- t(apply(mapdat[, ..compelem], 1, function(x) (x-predmapmean)/(pmax(predmapmax-predmapmean, predmapmean-predmapmin)*2)+0.5))

mapdat_rgb <- apply(mapdat[, ..compelem], 2, function(x) (x-mean(x))/(pmax(max(x)-mean(x), mean(x)-min(x))*2)+0.5)
mapdat_rgb <- apply(mapdat_rgb, 2, function(x) scurve(x, 0.75))
mapdat_lab <- convertColor(mapdat_rgb, from = "CIE RGB", to = "Lab")
mapdat_lab <- sweep(mapdat_lab, 2, c(0,12,4), FUN = "+")
mapdat_lum <- (mapdat_lab[,1]/100)
mapdat_lab[,1] <- (scurve(mapdat_lum, 0.8))*100
mapdat_back <- convertColor(mapdat_lab, from = "Lab", to = "CIE RGB", clip = TRUE)
mapdat_col <- rgb(mapdat_back)

mapdat$col <- mapdat_col

# Load bedrock shapefile ####

bedrock <- vect("data/digmap625_bedrock_arc/625k_V5_BEDROCK_Geology_Polygons.shp")
bedrock

bedrock_sf <- sf::st_as_sf(bedrock)
bedrock_sf <- sf::st_transform(bedrock_sf, sf::st_crs(elev))
crs(bedrock_sf)

bedrock_sf$N_vertices <- as.numeric(lapply(bedrock_sf$geometry, FUN = function(x) length(unlist(x))))
bedrock_sf$polygon_ID <- 1:nrow(bedrock_sf)

# Plot bedrock shapefile, coloured by number of vertices per polygon:
ggplot(data = bedrock_sf) + geom_sf(data = bedrock_sf, aes(fill = 10 - log(N_vertices)), lwd = 0) + theme_bw() +
  coord_sf(datum = sf::st_crs(27700)) + theme(legend.position="none") +
  labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)") +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks = seq(0,1000000, length.out = 6), expand = c(0,0), limits = c(0, 1225000)) +
  scale_x_continuous(breaks = seq(0,600000, length.out = 4), expand = c(0,0), limits = c(-200000, 664000))
ggsave(paste0("paperplots/625K_bedrock_polygons_by_n_vertices.png"), width = 144, height = 200, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

# Plot bedrock shapefile boundaries:
ggplot(data = bedrock_sf) + geom_sf(lwd = 0.2, fill = "gray96") + theme_bw() +
  coord_sf(datum = sf::st_crs(27700)) + theme(legend.position="none") +
  labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)") +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks = seq(0,1000000, length.out = 6), expand = c(0,0), limits = c(0, 1225000)) +
  scale_x_continuous(breaks = seq(0,600000, length.out = 4), expand = c(0,0), limits = c(-200000, 664000))
ggsave(paste0("paperplots/625K_bedrock_polygons_boundaries.png"), width = 144, height = 200, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

hist(1 - 1 / log(bedrock_sf$N_vertices) - 0.25, breaks = 100)
bedrock_sf$boundary_colour <- alpha("black", 1 - 1 / log(bedrock_sf$N_vertices) - 0.25)
bedrock_sf$boundary_colour[1:10]

ggplot(data = bedrock_sf) + geom_sf(lwd = 0.2, col = bedrock_sf$boundary_colour, fill = "gray96") + theme_bw() +
  scale_color_identity() +
  coord_sf(datum = sf::st_crs(27700)) + theme(legend.position="none") +
  labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)") +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks = seq(0,1000000, length.out = 6), expand = c(0,0), limits = c(0, 1225000)) +
  scale_x_continuous(breaks = seq(0,600000, length.out = 4), expand = c(0,0), limits = c(-200000, 664000))
ggsave(paste0("paperplots/625K_bedrock_polygons_boundaries_alpha_nvert.png"), width = 144, height = 200, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

# Plot point samples by composition, with bedrock linework overlain:
ggplot(data = mapdat[sample(1:nrow(mapdat)),]) + geom_point(aes(x = Easting_BNG, y = Northing_BNG, col = col), size = 0.5, shape = 16, alpha = 1) +
  geom_sf(data = bedrock_sf, fill = NA, col = alpha("black", 0.25)) + coord_sf(datum = sf::st_crs(27700)) +
  theme_bw() + scale_colour_identity() + labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)") +
  theme(legend.justification = c(1, 1), legend.position = c(0.98, 0.90), legend.background=element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks = seq(0,1000000, length.out = 6), expand = c(0,0), limits = c(0, 1225000)) +
  scale_x_continuous(breaks = seq(0,600000, length.out = 4), expand = c(0,0), limits = c(-200000, 664000))
ggsave(paste0("paperplots/rgb_", paste0(compelem, collapse = "_"), "_points_clr_curve_bedrock_lines_black.png"), width = 144, height = 200, units = "mm", type = "cairo", dpi = 300, scale = 1.375)


# Prepare data for BNN training: Extract sample-centred terrain images, of specified dimension and resolution ####
imagedim <- 27
imageres <- 250

rm(imgs)
imgs <- array(dim = c(nrow(mapdat), imagedim, imagedim))

seq.int(from = imageres/2, by = imageres, length.out = imagedim)-(imageres*imagedim)/2

cells <- as.data.table(expand.grid(x = seq.int(from = imageres/2, by = imageres, length.out = imagedim)-(imageres*imagedim)/2, 
                                   y = seq.int(from = imageres/2, by = imageres, length.out = imagedim)-(imageres*imagedim)/2))
cells[ ,coordx := rep(1:imagedim, imagedim),]
cells[ ,coordy := rep(1:imagedim, each = imagedim),]

# The sacred line (parallel processing for image extraction):
time <- Sys.time()
imgs <- array(as.numeric(unlist(mclapply(1:nrow(cells), FUN = function(i) extract(bielev, method = "bilinear", cbind(mapdat$Easting_BNG + cells[i,]$x, mapdat$Northing_BNG + cells[i,]$y))$alt,
                                         mc.cores = 8L))), dim = c(nrow(mapdat), imagedim, imagedim))
imgs[is.na(imgs)] = 0
Sys.time() - time

# Plot an image at random:
ggplot(reshape2::melt(imgs[sample(1:nrow(mapdat), 1),,]), aes(Var1,Var2, fill=value)) + 
  geom_raster() + theme_bw() + scale_fill_viridis_c(option = "B", name = "elevation") +
  coord_fixed()

# Scale sample-centered terrain images for use in neural network ####

mean = mean(as.data.frame(bielev, na.rm = TRUE)[,1])
sd =  sd(as.data.frame(bielev, na.rm = TRUE)[,1])

imgs_ann <- imgs - mapdat$elevation
imgs_ann <- (imgs_ann)/sd
imgs_ann[is.na(imgs_ann)] = 0
#imgs_ann <- imgs_ann - rowMeans(imgs_ann, na.rm = TRUE)

# Plot a normalised image at random
ggplot(reshape2::melt(imgs_ann[sample(1:nrow(mapdat), 1),,]), aes(Var1,Var2, fill=value)) + 
  geom_raster() + theme_bw() + scale_fill_viridis_c(option = "B", name = "normalised \n elevation") +
  coord_fixed()

# Extract easting, northing and elevation as location variables

loc <- mapdat[, c("Easting_BNG", "Northing_BNG", "elevation"), with = FALSE]

locmean <- apply(loc, 2, mean)
locsd <- apply(loc, 2, sd)

loc_ann <- t(apply(loc, 1, function(x) (x - locmean)/locsd))

# Plot and save three random terrain images to use as figure in paper
set.seed(678)
img1 <- reshape2::melt(imgs_ann[sample(1:nrow(mapdat), 1),,])
img2 <- reshape2::melt(imgs_ann[sample(1:nrow(mapdat), 1),,])
img3 <- reshape2::melt(imgs_ann[sample(1:nrow(mapdat), 1),,])

ggplot(img1, aes(Var1,Var2, fill=value)) + 
  geom_raster() + theme_bw() + 
  scale_fill_viridis_c(option = "cividis", limits = c(min(rbind(img1,img2,img3)$value),max(rbind(img1,img2,img3)$value))) +
  coord_fixed() + labs(x = "", y = "") + theme(legend.position = "none")
ggsave(paste0("plots/img1.png"), width = 45, height = 45, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

ggplot(img2, aes(Var1,Var2, fill=value)) + 
  geom_raster() + theme_bw() +
  scale_fill_viridis_c(option = "cividis", limits = c(min(rbind(img1,img2,img3)$value),max(rbind(img1,img2,img3)$value))) +
  coord_fixed() + labs(x = "", y = "") + theme(legend.position = "none")
ggsave(paste0("plots/img2.png"), width = 45, height = 45, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

ggplot(img3, aes(Var1,Var2, fill=value)) + 
  geom_raster() + theme_bw() + 
  scale_fill_viridis_c(option = "cividis", limits = c(min(rbind(img1,img2,img3)$value),max(rbind(img1,img2,img3)$value))) +
  coord_fixed() + labs(x = "", y = "") + theme(legend.position = "none")
ggsave(paste0("plots/img3.png"), width = 45, height = 45, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

# Use neural network to learn relationship between terrain features and geochemistry ####

set.seed(321)
fold_size = nrow(mapdat)/20
test <- sample(1:nrow(mapdat), fold_size)
val <- sample(which(!1:nrow(mapdat) %in% test), fold_size)
#train <- sample(which(!1:nrow(mapdat) %in% c(test, val)), length(which(!1:nrow(mapdat) %in% c(test, val)))*0.1)
train <- which(!1:nrow(mapdat) %in% c(test, val))

# Data Preparation --------------------------------------------------------

x_train <- list(imgs_ann[train,,], loc_ann[train, ])
x_val <- list(imgs_ann[val,,], loc_ann[val, ])
x_test <- list(imgs_ann[test,,], loc_ann[test, ])

dim(x_train[[1]]) <- c(nrow(x_train[[1]]), imagedim, imagedim, 1)
dim(x_val[[1]]) <- c(nrow(x_val[[1]]), imagedim, imagedim, 1)
dim(x_test[[1]]) <- c(nrow(x_test[[1]]), imagedim, imagedim, 1)

y_train <- as.matrix(mapdat[, c(4,5,6)][train])
y_val <- as.matrix(mapdat[, c(4,5,6)][val])
y_test <- as.matrix(mapdat[, c(4,5,6)][test])

# Defining Model ----------------------------------------------------------
time <- Sys.time()

dropratespat <- 0.8 #match what model was trained with
dropratedense <- 0.4 #match what model was trained with

# Convolutional stack:
conv_input <- layer_input(shape = c(imagedim, imagedim, 1))

conv_output <- conv_input %>%
  layer_conv_2d(filter = 256, kernel_size = c(3,3), dilation_rate = 1, strides = 3, kernel_initializer = initializer_he_normal()) %>%
  activation_relu() %>%
  layer_batch_normalization(momentum = 0.75) %>%
  layer_spatial_dropout_2d(rate = dropratespat) %>%
  layer_dropout(rate = dropratedense) %>%
  
  layer_conv_2d(filter = 256, kernel_size = c(3,3), dilation_rate = 1, strides = 3, kernel_initializer = initializer_he_normal()) %>%
  activation_relu() %>%
  layer_batch_normalization(momentum = 0.75) %>%
  layer_spatial_dropout_2d(rate = dropratespat) %>%
  layer_dropout(rate = dropratedense) %>%
  
  layer_conv_2d(filter = 256, kernel_size = c(3,3), dilation_rate = 1, strides = 1, kernel_initializer = initializer_he_normal()) %>%
  activation_relu() %>%
  layer_batch_normalization(momentum = 0.75) %>%
  layer_spatial_dropout_2d(rate = dropratespat) %>%
  layer_dropout(rate = dropratedense) %>%
  
  layer_flatten()

# Auxiliary input:
aux_input <- layer_input(shape = c(3))

aux_output <- aux_input %>%
  layer_dense(3840, kernel_initializer = initializer_he_normal()) %>%
  activation_relu() %>%
  layer_batch_normalization(momentum = 0.75) %>%
  layer_dropout(rate = dropratedense)

# Main process:
main_process <- layer_concatenate(c(conv_output, aux_output)) %>%
  layer_dense(2048, kernel_initializer = initializer_he_normal()) %>%
  activation_relu() %>%
  layer_batch_normalization(momentum = 0.75) %>%
  layer_dropout(rate = dropratedense) %>%
  layer_dense(1024, kernel_initializer = initializer_he_normal()) %>%
  activation_relu() %>%
  layer_batch_normalization(momentum = 0.75) %>%
  layer_dropout(rate = dropratedense)

# main_output <- layer_concatenate(c(main_process, conv_output)) %>%
main_output <- main_process %>%
  layer_dense(units = 6, activation = "linear", dtype = "float32", name = "dist_param") %>%
  layer_distribution_lambda(function(x)
    tfd_multivariate_normal_diag(loc = x[, 1:3, drop = FALSE],
                                 scale_diag = tf$math$softplus(0.1*x[, 4:6, drop = FALSE]))
  )

rm(model)
model <- keras_model(
  inputs = c(conv_input, aux_input),
  outputs = main_output
)

summary(model)

negloglik <- function(y, model) - (model %>% tfd_log_prob(y))

opt <- optimizer_nadam(learning_rate = 0.004, beta_1 = 0.95, use_ema = TRUE, ema_momentum = 0.95)

model %>% compile(
  loss = negloglik,
  optimizer = opt
)

# Training ----------------------------------------------------------------
batch_size <- 2^12
epochs <- 5000

# set.seed(999)
# 
# model %>% fit(
#   x_train, y_train,
#   batch_size = batch_size,
#   epochs = 50,
#   validation_data = list(x_val, y_val),
#   shuffle = TRUE,
#   verbose = 2)
# 
# history <- model %>% fit(
#   x_train, y_train,
#   batch_size = batch_size,
#   epochs = epochs,
#   validation_data = list(x_val, y_val),
#   shuffle = TRUE,
#   verbose = 2,
#   callbacks = list(callback_early_stopping(monitor = "val_loss", patience = 1500),
#                    callback_model_checkpoint(monitor = "val_loss", save_best_only = TRUE, save_weights_only = TRUE,
#                                              filepath = paste0(getwd(), "/models/", paste0(compelem, collapse = "_"), "_4090_250m_bn_tf215_modelweights_digital_geoscience_best.hdf5")))
# )
# min(na.omit(history$metrics$val_loss))
# Sys.time() - time
# 
# save_model_weights_hdf5(object = model, filepath = paste0(getwd(), "/models/", paste0(compelem, collapse = "_"), "_4090_250m_bn_tf215_modelweights_digital_geoscience_end.hdf5"))
# 
# min(na.omit(history$metrics$val_loss))
# which(na.omit(history$metrics$val_loss) == min(na.omit(history$metrics$val_loss)))
# 
# training_plot <- data.table(epoch = 1:length(history$metrics[[1]]), training = history$metrics$loss, validation = history$metrics$val_loss)
# training_plot <- as.data.table(reshape2::melt(training_plot, id.vars = "epoch", variable.name = "split", value.name = "negloglik"))
# 
# ggplot(training_plot[epoch <= 4000, ]) + theme_bw() +
#   geom_point(aes(x = epoch, y = negloglik, group = split, col = split), alpha = 0.4, shape = 16, size = 0.75) +
#   theme(legend.position = c(0.98, 0.98), legend.justification = c(1,1)) +
#   labs(x = "Training epoch", y = "Negative log-likelihood")
# ggsave(paste0("plots/BNN_training_progress_2.png"), width = 152, height = 90, units = "mm", type = "cairo", dpi = 300, scale = 1.375)


load_model_weights_hdf5(object = model, filepath = paste0(getwd(), "/models/", paste0(compelem, collapse = "_"), "_4090_250m_bn_tf215_modelweights_digital_geoscience_best.hdf5"))

meanmodel <- keras_model(
  inputs = model$input, 
  outputs = get_layer(model, "dist_param")$output
)

# Create maps of target variable using the trained bnn

# National scale maps ####
rm(predimgs, predimgs_ann, predloc, predloc_ann)

bielevcoarse <- aggregate(rast("data/British Isles 250m Copernicus DEM.tif"), fact=1) #500m grid cell map takes ~3.5mins to generate terrain patches for ANN (which seems a good trade-off)
predgrid <- as.data.table(as.data.frame(bielevcoarse, xy = TRUE, na.rm = TRUE))[y > 0]
names(predgrid)[3] <- "elevation"
nrow(predgrid)

# ggplot(predgrid) + geom_raster(aes(x = x, y = y, fill = elevation)) + theme_bw() + scale_fill_viridis_c(option = "B") + coord_equal()

# # The sacred line (parallel processing of images, or if already complete just load the RDS files):
# time <- Sys.time()
# predimgs <- array(as.numeric(unlist(mclapply(1:nrow(cells), FUN = function(i) extract(bielev, method = "bilinear", cbind(predgrid$x + cells[i,]$x, predgrid$y + cells[i,]$y))$alt,
#                                              mc.cores = 4L))), dim = c(nrow(predgrid), imagedim, imagedim))
# predimgs[is.na(predimgs)] = 0
# Sys.time() - time
# 
# predimgs_ann <- predimgs - predgrid$elevation
# predimgs_ann <- (predimgs_ann)/sd
# predimgs_ann[is.na(predimgs_ann)] = 0
# #predimgs_ann <- predimgs_ann - rowMeans(predimgs_ann, na.rm = TRUE)
# dim(predimgs_ann) <- c(nrow(predgrid), imagedim, imagedim, 1)
# 
# saveRDS(predimgs_ann, "data/predimgs_ann_27x27x250_at_250.RDS", compress = FALSE)
# 
# # Location information:
# predloc <- predgrid[, c("x", "y", "elevation"), with = FALSE]
# predlocmean <- apply(loc, 2, mean)
# predlocsd <- apply(loc, 2, sd)
# predloc_ann <- t(apply(predloc, 1, function(x) (x - predlocmean)/predlocsd))
# 
# saveRDS(predloc_ann, "data/predloc_ann_250.RDS", compress = FALSE)

predimgs_ann <- readRDS("data/predimgs_ann_27x27x250_at_250.RDS")
predloc_ann <- readRDS("data/predloc_ann_250.RDS")

# with(tf$device('/CPU:0'),
# preds <- cbind(predgrid, predict(meanmodel, list(predimgs_ann, predloc_ann), batch_size = 512)[,1:3])
# )

preds <- cbind(predgrid, predict(meanmodel, list(predimgs_ann, predloc_ann), batch_size = 512)[,1:3])
names(preds)[c(4,5,6)] <- compelem

predmapmean <- apply(preds[,4:6], 2, mean)
predmapsd <- apply(preds[,4:6], 2, sd)
predmapmax <- apply(preds[,4:6], 2, max)
predmapmin <- apply(preds[,4:6], 2, min)

# rgb_adj <- apply(preds[,4:6], 2, function(x) (x-mean(x))/(pmax(max(x)-mean(x), mean(x)-min(x))*2)+0.5)*0.95

rgb_adj <- t(apply(preds[,4:6], 1, function(x) (x-predmapmean)/(pmax(predmapmax-predmapmean, predmapmean-predmapmin)*2)+0.5))*0.95

col_lab <- convertColor(rgb_adj, from = "CIE RGB", to = "Lab")
col_lab <- sweep(col_lab, 2, c(0,4,2), FUN = "+")
luminance <- (col_lab[,1]/100)
col_lab[,1] <- scurve(0*ecdf(luminance)(luminance)+1*luminance, 0.82)*100
# col_lab[,1] <- (1*scurve(luminance, 0.85)-0.5*scurve(luminance, 0.95))*100
# hist(col_lab[,1], breaks = 100)
back_col <- convertColor(col_lab, from = "Lab", to = "CIE RGB", clip = TRUE)
simcol <- rgb(back_col)

ggplot(preds) + geom_raster(aes(x = x, y = y, fill = simcol), interpolate = TRUE) + theme_bw() + coord_equal() +
  scale_fill_identity() + labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)") +
  theme(legend.justification = c(1, 1), legend.position = c(0.98, 0.90), legend.background=element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks = seq(0,1000000, length.out = 6), expand = c(0,0), limits = c(0, 1225000)) +
  scale_x_continuous(breaks = seq(0,600000, length.out = 4), expand = c(0,0), limits = c(-200000, 664000))
ggsave(paste0("paperplots/rgb_", paste0(compelem, collapse = "_"), "_bnn_map_clr_curve.png"), width = 144, height = 200, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

# Colour polygons by averaging their geochemical composition of samples they contain ####
# Intersect points with polygons (which polygon is each geochemical sample in?)
mapdat_sf <- sf::st_as_sf(mapdat, coords = c("Easting_BNG", "Northing_BNG", "elevation"), remove = FALSE)
sf::st_crs(mapdat_sf) <- sf::st_crs(elev)

intersection <- sf::st_intersects(bedrock_sf, mapdat_sf, sparse = FALSE)
samples_polygons <- as.numeric(apply(intersection, 2, FUN = function(x) which(x == TRUE)))

mapdat$polygon <- samples_polygons

# Average geochemical composition within each polygon #####
bedrock_dt <- as.data.table(bedrock_sf)

r <- 1
for(r in 1:nrow(bedrock_dt)){
  bedrock_dt[r, c("K2O_sample_avg", "Fe2O3_sample_avg", "CaO_sample_avg") := as.list(colMeans(mapdat[train,][polygon == r, c("K2O", "Fe2O3", "CaO")])), ]
}

# Record which columns don't contain NAs (we can only pass non-NA rows through the RGB conversion, then have to map these back)
bedrock_dt_col_non_na_rows <- which(rowSums(is.na(bedrock_dt[, c("K2O_sample_avg", "Fe2O3_sample_avg", "CaO_sample_avg")])) == 0)

bedrock_rgb <- t(apply(na.omit(bedrock_dt[, c("K2O_sample_avg", "Fe2O3_sample_avg", "CaO_sample_avg")]), 1, function(x) (x-predmapmean)/(pmax(predmapmax-predmapmean, predmapmean-predmapmin)*2)+0.5))*0.95

# bedrock_rgb <- apply(bedrock_rgb, 2, function(x) scurve(x, 0.82))
bedrock_lab <- convertColor(bedrock_rgb, from = "CIE RGB", to = "Lab")
bedrock_lab <- sweep(bedrock_lab, 2, c(0,4,2), FUN = "+")
bedrock_lum <- (bedrock_lab[,1]/100)
bedrock_lab[,1] <- (scurve(bedrock_lum, 0.82))*100
bedrock_back <- convertColor(bedrock_lab, from = "Lab", to = "CIE RGB", clip = TRUE)
bedrock_sample_col <- rgb(bedrock_back)

bedrock_dt[bedrock_dt_col_non_na_rows, sample_avg_col := bedrock_sample_col, ]
bedrock_dt[!bedrock_dt_col_non_na_rows, sample_avg_col := "#000000", ]

bedrock_sf <- sf::st_as_sf(bedrock_dt)
# bedrock_sf <- sf::st_transform(bedrock_sf, sf::st_crs(elev))

ggplot(data = bedrock_sf) + geom_sf(data = bedrock_sf, aes(fill = sample_avg_col), color = NA, lwd = 1) + theme_bw() +
  scale_fill_identity() + coord_sf(datum = sf::st_crs(27700)) + theme(legend.position="none") +
  labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)") +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks = seq(0,1000000, length.out = 6), expand = c(0,0), limits = c(0, 1225000)) +
  scale_x_continuous(breaks = seq(0,600000, length.out = 4), expand = c(0,0), limits = c(-200000, 664000))
ggsave(paste0("paperplots/rgb_", paste0(compelem, collapse = "_"), "_polygon_sample_avgs_bnn_based_clr_curve.png"), width = 144, height = 200, units = "mm", type = "cairo", dpi = 300, scale = 1.375)


# Colour polygons by averaging their geochemical composition of BNN PREDICTIONS they contain ####
# Intersect points with polygons (which polygon is each geochemical sample in?)

sampleofpreds <- sample(1:nrow(preds), nrow(preds)/10)

preds_sf <- sf::st_as_sf(preds[sampleofpreds,], coords = c("x", "y", "elevation"), remove = FALSE)
sf::st_crs(preds_sf) <- sf::st_crs(elev)

preds_intersection <- sf::st_intersects(bedrock_sf, preds_sf, sparse = FALSE)
preds_samples_polygons <- as.numeric(apply(preds_intersection, 2, FUN = function(x) which(x == TRUE)))

preds[sampleofpreds, polygon := preds_samples_polygons, ]
preds

# Average geochemical composition within each polygon #####

r <- 1
for(r in 1:nrow(bedrock_dt)){
  bedrock_dt[r, c("K2O_bnn_avg", "Fe2O3_bnn_avg", "CaO_bnn_avg") := as.list(colMeans(preds[polygon == r, c("K2O", "Fe2O3", "CaO")])), ]
}

# Record which columns don't contain NAs (we can only pass non-NA rows through the RGB conversion, then have to map these back)
bedrock_dt_col_non_na_rows <- which(rowSums(is.na(bedrock_dt[, c("K2O_bnn_avg", "Fe2O3_bnn_avg", "CaO_bnn_avg")])) == 0)

bedrock_rgb <- t(apply(na.omit(bedrock_dt[, c("K2O_bnn_avg", "Fe2O3_bnn_avg", "CaO_bnn_avg")]), 1, function(x) (x-predmapmean)/(pmax(predmapmax-predmapmean, predmapmean-predmapmin)*2)+0.5))*0.95

# bedrock_rgb <- apply(bedrock_rgb, 2, function(x) scurve(x, 0.82))
bedrock_lab <- convertColor(bedrock_rgb, from = "CIE RGB", to = "Lab")
bedrock_lab <- sweep(bedrock_lab, 2, c(0,4,2), FUN = "+")
bedrock_lum <- (bedrock_lab[,1]/100)
bedrock_lab[,1] <- (scurve(bedrock_lum, 0.82))*100
bedrock_back <- convertColor(bedrock_lab, from = "Lab", to = "CIE RGB", clip = TRUE)
bedrock_bnn_col <- rgb(bedrock_back)

bedrock_dt[bedrock_dt_col_non_na_rows, bnn_avg_col := bedrock_bnn_col, ]
bedrock_dt[!bedrock_dt_col_non_na_rows, bnn_avg_col := "#000000", ]

bedrock_sf <- sf::st_as_sf(bedrock_dt)
# bedrock_sf <- sf::st_transform(bedrock_sf, sf::st_crs(elev))

ggplot(data = bedrock_sf) + geom_sf(data = bedrock_sf, aes(fill = bnn_avg_col), color = NA, lwd = 1) + theme_bw() +
  scale_fill_identity() + coord_sf(datum = sf::st_crs(27700)) + theme(legend.position="none") +
  labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)") +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  scale_y_continuous(breaks = seq(0,1000000, length.out = 6), expand = c(0,0), limits = c(0, 1225000)) +
  scale_x_continuous(breaks = seq(0,600000, length.out = 4), expand = c(0,0), limits = c(-200000, 664000))
ggsave(paste0("paperplots/rgb_", paste0(compelem, collapse = "_"), "_polygon_preds_avgs_bnn_based_clr_curve.png"), width = 144, height = 200, units = "mm", type = "cairo", dpi = 300, scale = 1.375)

# Assess performance of different 'models' on held out test set:

  holdout <- as.data.table(data.frame(Element = factor(rep(gsub("2.*|O.*", "", compelem), each = nrow(y_test)), levels = gsub("2.*|O.*", "", compelem)), obs = as.numeric(y_test), preds = as.numeric(predict(meanmodel, x_test)[,1:3])))
  
  ggplot(holdout[sample(1:nrow(holdout), nrow(holdout))]) + theme_bw() + geom_abline(slope = 1, intercept = 0) + theme(aspect.ratio = 1) + 
    coord_fixed(ratio = 1, xlim = c(-4,3), ylim = c(-4,3)) +
    geom_point(aes(x = obs, y = preds, col = Element), shape = 16, size = 0.45, alpha = 0.25) +
    labs(x = "Observed (centred log-ratio)", y = "Predicted (centred log-ratio)",
         subtitle = paste0("R\U00B2 for ", gsub("2.*|O.*", "", compelem)[1], " = ", round(cor(holdout[Element == gsub("2.*|O.*", "", compelem)[1], c(2:3)])[2]^2, 2), "   ", 
                           "R\U00B2 for ", gsub("2.*|O.*", "", compelem)[2], " = ", round(cor(holdout[Element == gsub("2.*|O.*", "", compelem)[2], c(2:3)])[2]^2, 2), "   ", 
                           "R\U00B2 for ", gsub("2.*|O.*", "", compelem)[3], " = ", round(cor(holdout[Element == gsub("2.*|O.*", "", compelem)[3], c(2:3)])[2]^2, 2), "\n",
                           "RMSE K = ", round(sqrt(mean((holdout[Element == gsub("2.*|O.*", "", compelem)[1]]$preds-holdout[Element == gsub("2.*|O.*", "", compelem)[1]]$obs)^2)), 2), "   ",
                           "RMSE Fe = ", round(sqrt(mean((holdout[Element == gsub("2.*|O.*", "", compelem)[2]]$preds-holdout[Element == gsub("2.*|O.*", "", compelem)[2]]$obs)^2)), 2), "   ",
                           "RMSE Ca = ", round(sqrt(mean((holdout[Element == gsub("2.*|O.*", "", compelem)[3]]$preds-holdout[Element == gsub("2.*|O.*", "", compelem)[3]]$obs)^2)), 2))) +
    theme(legend.position = c(0.05, 0.95), legend.justification = c(0,1), plot.subtitle = element_text(hjust = 0.5), plot.margin = margin(4,10,4,4))
  ggsave(paste0("paperplots/", paste0(compelem, collapse = "_"), "_BNN_mean_holdout_zoomed_RMSE.png"), width = 89, height = 89, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
  
  # Now for each element individually (as requested by reviewer):
    ggplot(holdout[sample(1:nrow(holdout), nrow(holdout))][Element == "K"]) + theme_bw() + geom_abline(slope = 1, intercept = 0) + theme(aspect.ratio = 1) + 
      coord_fixed(ratio = 1, xlim = c(-3,2), ylim = c(-3,2), expand = c(0,0)) +
      geom_point(aes(x = obs, y = preds), col = "#F8766D", shape = 16, size = 0.65, alpha = 0.5) +
      labs(x = "Observed (centred log-ratio)", y = "Predicted (centred log-ratio)",
           subtitle = paste0("R\U00B2 for ", gsub("2.*|O.*", "", compelem)[1], " = ", round(cor(holdout[Element == gsub("2.*|O.*", "", compelem)[1], c(2:3)])[2]^2, 2), "   ", 
                             "RMSE K = ", round(sqrt(mean((holdout[Element == gsub("2.*|O.*", "", compelem)[1]]$preds-holdout[Element == gsub("2.*|O.*", "", compelem)[1]]$obs)^2)), 2))) +
      theme(legend.position = c(0.05, 0.95), legend.justification = c(0,1), plot.subtitle = element_text(hjust = 0.5), plot.margin = margin(4,10,4,4))
    ggsave(paste0("paperplots/", paste0(compelem[1], collapse = "_"), "_BNN_mean_holdout_zoomed_RMSE.png"), width = 89, height = 89, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
    
    ggplot(holdout[sample(1:nrow(holdout), nrow(holdout))][Element == "Fe"]) + theme_bw() + geom_abline(slope = 1, intercept = 0) + theme(aspect.ratio = 1) + 
      coord_fixed(ratio = 1, xlim = c(-2,3), ylim = c(-2,3), expand = c(0,0)) +
      geom_point(aes(x = obs, y = preds), col = "#00BA38", shape = 16, size = 0.65, alpha = 0.5) +
      labs(x = "Observed (centred log-ratio)", y = "Predicted (centred log-ratio)",
           subtitle = paste0("R\U00B2 for ", gsub("2.*|O.*", "", compelem)[2], " = ", round(cor(holdout[Element == gsub("2.*|O.*", "", compelem)[2], c(2:3)])[2]^2, 2), "   ", 
                             "RMSE Fe = ", round(sqrt(mean((holdout[Element == gsub("2.*|O.*", "", compelem)[2]]$preds-holdout[Element == gsub("2.*|O.*", "", compelem)[2]]$obs)^2)), 2))) +
      theme(legend.position = c(0.05, 0.95), legend.justification = c(0,1), plot.subtitle = element_text(hjust = 0.5), plot.margin = margin(4,10,4,4))
    ggsave(paste0("paperplots/", paste0(compelem[2], collapse = "_"), "_BNN_mean_holdout_zoomed_RMSE.png"), width = 89, height = 89, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
    
    ggplot(holdout[sample(1:nrow(holdout), nrow(holdout))][Element == "Ca"]) + theme_bw() + geom_abline(slope = 1, intercept = 0) + theme(aspect.ratio = 1) + 
      coord_fixed(ratio = 1, xlim = c(-4,3), ylim = c(-4,3), expand = c(0,0)) +
      geom_point(aes(x = obs, y = preds), col = "#619CFF", shape = 16, size = 0.65, alpha = 0.5) +
      labs(x = "Observed (centred log-ratio)", y = "Predicted (centred log-ratio)",
           subtitle = paste0("R\U00B2 for ", gsub("2.*|O.*", "", compelem)[3], " = ", round(cor(holdout[Element == gsub("2.*|O.*", "", compelem)[3], c(2:3)])[2]^2, 2), "   ", 
                             "RMSE Ca = ", round(sqrt(mean((holdout[Element == gsub("2.*|O.*", "", compelem)[3]]$preds-holdout[Element == gsub("2.*|O.*", "", compelem)[3]]$obs)^2)), 2))) +
      theme(legend.position = c(0.05, 0.95), legend.justification = c(0,1), plot.subtitle = element_text(hjust = 0.5), plot.margin = margin(4,10,4,4))
    ggsave(paste0("paperplots/", paste0(compelem[3], collapse = "_"), "_BNN_mean_holdout_zoomed_RMSE.png"), width = 89, height = 89, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
  
  
  # Compare polygon sample averages to geochemical observations:
  r <- 1
  for(r in 1:nrow(bedrock_dt)){
    mapdat[polygon == r, c("K2O_sample_avg", "Fe2O3_sample_avg", "CaO_sample_avg") := as.list(colMeans(mapdat[train,][polygon == r, c("K2O", "Fe2O3", "CaO")])), ]
  }
  
  mapdat
  
  holdout_poly_samp <- na.omit(as.data.table(data.frame(Element = factor(rep(gsub("2.*|O.*", "", compelem), each = nrow(y_test)), levels = gsub("2.*|O.*", "", compelem)), 
                                                       obs = as.numeric(y_test), preds = as.numeric(unlist(mapdat[test, c("K2O_sample_avg", "Fe2O3_sample_avg", "CaO_sample_avg")])))))
  
  ggplot(holdout_poly_samp[sample(1:nrow(holdout_poly_samp), nrow(holdout_poly_samp))]) + theme_bw() + geom_abline(slope = 1, intercept = 0) + theme(aspect.ratio = 1) + 
    coord_fixed(ratio = 1, xlim = c(-4,3), ylim = c(-4,3)) +
    geom_point(aes(x = obs, y = preds, col = Element), shape = 16, size = 0.45, alpha = 0.25) +
    labs(x = "Observed (centred log-ratio)", y = "Predicted (centred log-ratio)",
         subtitle = paste0("R\U00B2 for ", gsub("2.*|O.*", "", compelem)[1], " = ", round(cor(holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[1], c(2:3)])[2]^2, 2), "   ", 
                           "R\U00B2 for ", gsub("2.*|O.*", "", compelem)[2], " = ", round(cor(holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[2], c(2:3)])[2]^2, 2), "   ", 
                           "R\U00B2 for ", gsub("2.*|O.*", "", compelem)[3], " = ", round(cor(holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[3], c(2:3)])[2]^2, 2), "\n",
                           "RMSE K = ", round(sqrt(mean((holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[1]]$preds-holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[1]]$obs)^2)), 2), "   ",
                           "RMSE Fe = ", round(sqrt(mean((holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[2]]$preds-holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[2]]$obs)^2)), 2), "   ",
                           "RMSE Ca = ", round(sqrt(mean((holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[3]]$preds-holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[3]]$obs)^2)), 2))) +
    theme(legend.position = c(0.05, 0.95), legend.justification = c(0,1), plot.subtitle = element_text(hjust = 0.5), plot.margin = margin(4,10,4,4))
  ggsave(paste0("paperplots/", paste0(compelem, collapse = "_"), "_sample_means_holdout_zoomed_RMSE.png"), width = 89, height = 89, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
  
  # Now for each element individually (as requested by reviewer):
    ggplot(holdout_poly_samp[sample(1:nrow(holdout_poly_samp), nrow(holdout_poly_samp))][Element == "K"]) + theme_bw() + geom_abline(slope = 1, intercept = 0) + theme(aspect.ratio = 1) + 
      coord_fixed(ratio = 1, xlim = c(-3,2), ylim = c(-3,2), expand = c(0,0)) +
      geom_point(aes(x = obs, y = preds), col = "#F8766D", shape = 16, size = 0.65, alpha = 0.5) +
      labs(x = "Observed (centred log-ratio)", y = "Predicted (centred log-ratio)",
           subtitle = paste0("R\U00B2 for ", gsub("2.*|O.*", "", compelem)[1], " = ", round(cor(holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[1], c(2:3)])[2]^2, 2), "   ", 
                             "RMSE K = ", round(sqrt(mean((holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[1]]$preds-holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[1]]$obs)^2)), 2))) +
      theme(legend.position = c(0.05, 0.95), legend.justification = c(0,1), plot.subtitle = element_text(hjust = 0.5), plot.margin = margin(4,10,4,4))
    ggsave(paste0("paperplots/", paste0(compelem[1], collapse = "_"), "_sample_means_holdout_poly_samp_zoomed_RMSE.png"), width = 89, height = 89, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
    
    ggplot(holdout_poly_samp[sample(1:nrow(holdout_poly_samp), nrow(holdout_poly_samp))][Element == "Fe"]) + theme_bw() + geom_abline(slope = 1, intercept = 0) + theme(aspect.ratio = 1) + 
      coord_fixed(ratio = 1, xlim = c(-2,3), ylim = c(-2,3), expand = c(0,0)) +
      geom_point(aes(x = obs, y = preds), col = "#00BA38", shape = 16, size = 0.65, alpha = 0.5) +
      labs(x = "Observed (centred log-ratio)", y = "Predicted (centred log-ratio)",
           subtitle = paste0("R\U00B2 for ", gsub("2.*|O.*", "", compelem)[2], " = ", round(cor(holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[2], c(2:3)])[2]^2, 2), "   ", 
                             "RMSE Fe = ", round(sqrt(mean((holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[2]]$preds-holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[2]]$obs)^2)), 2))) +
      theme(legend.position = c(0.05, 0.95), legend.justification = c(0,1), plot.subtitle = element_text(hjust = 0.5), plot.margin = margin(4,10,4,4))
    ggsave(paste0("paperplots/", paste0(compelem[2], collapse = "_"), "_sample_means_holdout_poly_samp_zoomed_RMSE.png"), width = 89, height = 89, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
    
    ggplot(holdout_poly_samp[sample(1:nrow(holdout_poly_samp), nrow(holdout_poly_samp))][Element == "Ca"]) + theme_bw() + geom_abline(slope = 1, intercept = 0) + theme(aspect.ratio = 1) + 
      coord_fixed(ratio = 1, xlim = c(-4,3), ylim = c(-4,3), expand = c(0,0)) +
      geom_point(aes(x = obs, y = preds), col = "#619CFF", shape = 16, size = 0.65, alpha = 0.5) +
      labs(x = "Observed (centred log-ratio)", y = "Predicted (centred log-ratio)",
           subtitle = paste0("R\U00B2 for ", gsub("2.*|O.*", "", compelem)[3], " = ", round(cor(holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[3], c(2:3)])[2]^2, 2), "   ", 
                             "RMSE Ca = ", round(sqrt(mean((holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[3]]$preds-holdout_poly_samp[Element == gsub("2.*|O.*", "", compelem)[3]]$obs)^2)), 2))) +
      theme(legend.position = c(0.05, 0.95), legend.justification = c(0,1), plot.subtitle = element_text(hjust = 0.5), plot.margin = margin(4,10,4,4))
    ggsave(paste0("paperplots/", paste0(compelem[3], collapse = "_"), "_sample_means_holdout_poly_samp_zoomed_RMSE.png"), width = 89, height = 89, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
  
  
  # Compare polygon BNN prediction averages to geochemical observations:
  r <- 1
  for(r in 1:nrow(bedrock_dt)){
    mapdat[polygon == r, c("K2O_bnn_avg", "Fe2O3_bnn_avg", "CaO_bnn_avg") := as.list(colMeans(preds[polygon == r, c("K2O", "Fe2O3", "CaO")])), ]
  }
  
  mapdat

  holdout_poly_bnn <- na.omit(as.data.table(data.frame(Element = factor(rep(gsub("2.*|O.*", "", compelem), each = nrow(y_test)), levels = gsub("2.*|O.*", "", compelem)), 
                                                       obs = as.numeric(y_test), preds = as.numeric(unlist(mapdat[test, c("K2O_bnn_avg", "Fe2O3_bnn_avg", "CaO_bnn_avg")])))))
  
  ggplot(holdout_poly_bnn[sample(1:nrow(holdout_poly_bnn), nrow(holdout_poly_bnn))]) + theme_bw() + geom_abline(slope = 1, intercept = 0) + theme(aspect.ratio = 1) + 
    coord_fixed(ratio = 1, xlim = c(-4,3), ylim = c(-4,3)) +
    geom_point(aes(x = obs, y = preds, col = Element), shape = 16, size = 0.45, alpha = 0.25) +
    labs(x = "Observed (centred log-ratio)", y = "Predicted (centred log-ratio)",
         subtitle = paste0("R\U00B2 for ", gsub("2.*|O.*", "", compelem)[1], " = ", round(cor(holdout_poly_bnn[Element == gsub("2.*|O.*", "", compelem)[1], c(2:3)])[2]^2, 2), "   ", 
                           "R\U00B2 for ", gsub("2.*|O.*", "", compelem)[2], " = ", round(cor(holdout_poly_bnn[Element == gsub("2.*|O.*", "", compelem)[2], c(2:3)])[2]^2, 2), "   ", 
                           "R\U00B2 for ", gsub("2.*|O.*", "", compelem)[3], " = ", round(cor(holdout_poly_bnn[Element == gsub("2.*|O.*", "", compelem)[3], c(2:3)])[2]^2, 2), "\n",
                           "RMSE K = ", round(sqrt(mean((holdout_poly_bnn[Element == gsub("2.*|O.*", "", compelem)[1]]$preds-holdout_poly_bnn[Element == gsub("2.*|O.*", "", compelem)[1]]$obs)^2)), 2), "   ",
                           "RMSE Fe = ", round(sqrt(mean((holdout_poly_bnn[Element == gsub("2.*|O.*", "", compelem)[2]]$preds-holdout_poly_bnn[Element == gsub("2.*|O.*", "", compelem)[2]]$obs)^2)), 2), "   ",
                           "RMSE Ca = ", round(sqrt(mean((holdout_poly_bnn[Element == gsub("2.*|O.*", "", compelem)[3]]$preds-holdout_poly_bnn[Element == gsub("2.*|O.*", "", compelem)[3]]$obs)^2)), 2))) +
    theme(legend.position = c(0.05, 0.95), legend.justification = c(0,1), plot.subtitle = element_text(hjust = 0.5), plot.margin = margin(4,10,4,4))
  ggsave(paste0("paperplots/", paste0(compelem, collapse = "_"), "_poly_bnn_holdout_zoomed_RMSE.png"), width = 89, height = 89, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
  
  # remotes::install_version("ggplot2", version = "3.4.4")
  
  # Ternary plot as legend ####
  library(ggplot2)
  library(ggtern)
  
  legendpix <- sample(1:nrow(preds), 1000000)
  ternplot <- ggtern(rgb_adj[legendpix,], aes(x = K2O, y = Fe2O3, z = CaO)) + geom_point(aes(col = simcol[legendpix]), size = 0.2, shape = 16) +
    scale_color_identity() + theme_bw(8) + theme_nolabels() + theme_noticks() + 
    labs(x = "K", xarrow = "relative enrichment in potassium", 
         y = "Fe", yarrow = "relative enrichment in iron", 
         z = "Ca", zarrow = "relative enrichment in calcium")
  ternplot
  ggsave(paste0("paperplots/", paste0(compelem, collapse = "_"), "_pred_map_ternary_legend.png"), width = 60, height = 60, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
  
  # Add this ternary plot as a legend to the map plots:
  
  # Convert to grob
  terngrob <- ggplotGrob(ternplot)

  # BNN map with legend:
  ggplot(preds) + geom_raster(aes(x = x, y = y, fill = simcol), interpolate = TRUE) + theme_bw() + coord_equal() +
    scale_fill_identity() + labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)") +
    theme(legend.justification = c(1, 1), legend.position = c(0.98, 0.90), legend.background=element_blank(),
          axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
    scale_y_continuous(breaks = seq(0,1000000, length.out = 6), expand = c(0,0), limits = c(0, 1225000)) +
    scale_x_continuous(breaks = seq(0,600000, length.out = 4), expand = c(0,0), limits = c(-200000, 664000)) +
    annotation_custom(grob = terngrob, xmin = -200000, xmax = 110000, ymin = 950000, ymax = 1225000)
  ggsave(paste0("paperplots/rgb_", paste0(compelem, collapse = "_"), "_bnn_map_clr_curve_legend.png"), width = 144, height = 200, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
  
  # Polygon map with legend:
  ggplot(data = bedrock_sf) + geom_sf(data = bedrock_sf, aes(fill = sample_avg_col), color = NA, lwd = 1) + theme_bw() +
    scale_fill_identity() + coord_sf(datum = sf::st_crs(27700)) + theme(legend.position="none") +
    labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)") +
    theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
    scale_y_continuous(breaks = seq(0,1000000, length.out = 6), expand = c(0,0), limits = c(0, 1225000)) +
    scale_x_continuous(breaks = seq(0,600000, length.out = 4), expand = c(0,0), limits = c(-200000, 664000)) +
    annotation_custom(grob = terngrob, xmin = -200000, xmax = 110000, ymin = 950000, ymax = 1225000)
  ggsave(paste0("paperplots/rgb_", paste0(compelem, collapse = "_"), "_polygon_sample_avgs_bnn_based_clr_curve_legend.png"), width = 144, height = 200, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
  
  # Raw observation map with legend:
  # First, legend needed (previous one used predictions not observations)
  
  mapdat_rgb_adj <- as.data.table(t(apply(mapdat[,4:6], 1, function(x) (x-predmapmean)/(pmax(predmapmax-predmapmean, predmapmean-predmapmin)*2)+0.5)))
  mapdat_rgb_adj[, col := mapdat$col, ]
  
  ternplot_obs <- ggtern(mapdat_rgb_adj, aes(x = K2O, y = Fe2O3, z = CaO)) + geom_point(aes(col = col), size = 0.2, shape = 16) +
    scale_color_identity() + theme_bw(8) + theme_nolabels() + theme_noticks() + 
    labs(x = "K", xarrow = "relative enrichment in potassium", 
         y = "Fe", yarrow = "relative enrichment in iron", 
         z = "Ca", zarrow = "relative enrichment in calcium")
  ternplot_obs
  
  terngrob_obs <- ggplotGrob(ternplot_obs)
  
  ggplot(data = mapdat[sample(1:nrow(mapdat)),]) + geom_point(aes(x = Easting_BNG, y = Northing_BNG, col = col), size = 0.5, shape = 16, alpha = 1) +
    geom_sf(data = bedrock_sf, fill = NA, col = alpha("black", 0.25)) + coord_sf(datum = sf::st_crs(27700)) +
    theme_bw() + scale_colour_identity() + labs(x = "Easting (metres BNG)", y = "Northing (metres BNG)") +
    theme(legend.justification = c(1, 1), legend.position = c(0.98, 0.90), legend.background=element_blank(),
          axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
    scale_y_continuous(breaks = seq(0,1000000, length.out = 6), expand = c(0,0), limits = c(0, 1225000)) +
    scale_x_continuous(breaks = seq(0,600000, length.out = 4), expand = c(0,0), limits = c(-200000, 664000)) +
    annotation_custom(grob = terngrob_obs, xmin = -200000, xmax = 110000, ymin = 950000, ymax = 1225000)
  ggsave(paste0("paperplots/rgb_", paste0(compelem, collapse = "_"), "_points_clr_curve_bedrock_lines_black_legend.png"), width = 144, height = 200, units = "mm", type = "cairo", dpi = 300, scale = 1.375)
