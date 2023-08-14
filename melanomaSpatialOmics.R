options(java.parameters = "-Xmx4g")
library(RBioFormats)
library("tiff")
library(tidyverse)
library(readr)
library(SpatialOmicsOverlay)
library(GeomxTools)
library(nanostringr)
library(NanoStringNCTools)
library(ggthemes)
library(ggiraph)
library(pheatmap)
library(cowplot)
library("viridis") 


#add image here
tifFile <- "/Users/isabella/SpatialOmics_Isabella/Input_files_analysis_june2023/HT-0692_Slide 6_TMA.ome.tiff"
#tifFile <- image_read("/Users/isabella/SpatialOmics_Isabella/Input_files_analysis_june2023/HT-0692_Slide 6_TMA.ome.tiff")
#tifFile <- image_rotate(tifFile, 90)
tifFile

# add lab worksheet here
melanomaTMA_LW <- "/Users/isabella/SpatialOmics_Isabella/Input_files_analysis_june2023/Downloads_from_DSP/P1001660009598A_LabWorksheet.txt"
# melanoma <- readLabWorksheet(melanomaTMA_LW, slideName = "HT-0692_Slide_6_TMA")
# generates coordinates from LW data
melanomaTMA <- readSpatialOverlay(ometiff = tifFile, annots = melanomaTMA_LW, 
                              slideName = "HT-0692_Slide_6_TMA", image = FALSE,  
                              saveFile = FALSE, outline = FALSE)
#full object
melanomaTMA

#sample names
head(sampNames(melanomaTMA))

#slide name
slideName(melanomaTMA)

#metadata of ROI overlays
#Height, Width, X, Y values are in pixels
head(meta(overlay(melanomaTMA)))

#coordinates of each ROI
head(coords(melanomaTMA))

#PLOTTING WITHOUT IMAGE
plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, legend = FALSE)


melanomaAnnots <- readLabWorksheet(lw = melanomaTMA_LW, slideName = "HT-0692_Slide_6_TMA")
# below is the sample tutorial code for extracting DCC file info
#muBrainGeomxSet <- readRDS(unzip(system.file("extdata", "muBrain_GxT.zip", 
 #                                           package = "SpatialOmicsOverlay")))
datadir <- "/Users/isabella/SpatialOmics_Isabella/Input_files_analysis_june2023/RCCfiles/20221116_P1001660009598A_RCC"
# datadir <- read_rcc(datadir)
# reading rcc file data only (no rlf or other files)
rccs <- dir(datadir, pattern = ".RCC", full.names = TRUE)
melanomaDataSet <- readNanoStringRccSet(rccs)
# melanomaDataSet <- rccData


initial_data_csv <- read_csv("/Users/isabella/SpatialOmics_Isabella/Input_files_analysis_june2023/RCCfiles/20221116_P1001660009598A_RCC/initial_dataset_6071_geomx_dsp_roi_collected_archivo_anotacion_final_060523.csv")
# load in initial dataset from annotation sheet
initial_data_dataframe <- as.data.frame(initial_data_csv)
# convert from csv to data.frame



# plot by segments
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = melanomaAnnots, 
                             plottingFactor = "segment")
# code stops here: cannot find target names in NanoString RCC Data Set
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                             plottingFactor = "CD56")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                             plottingFactor = "ROILabel")
melanomaTMA

head(plotFactors(melanomaTMA))

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD56", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD56 Expression in Melanoma TMA")


#lowest resolution = fastest speeds
checkValidRes(ometiff = tifFile)

# optimal res recommended by SpatialOmics
res <- 1
# add image back in to pot on top of
melanomaTMA <- addImageOmeTiff(overlay = melanomaTMA, ometiff = tifFile, res = res)
melanomaTMA

# only image
showImage(melanomaTMA)

# plots both image and plot now
# plotSpatialOverlay(overlay = melanomaTMA, colorBy = "segment", corner = "topcenter", 
#                    scaleBarWidth = 0.5, textDistance = 130, scaleBarColor = "cyan")
# 
plotSpatialOverlay(overlay = melanomaTMA, colorBy = "segment", corner = "topcenter", 
                    scaleBarWidth = 0.5, textDistance = 130, scaleBarColor = "cyan",
                   fluorLegend = TRUE)
# 
 gp <- plotSpatialOverlay(overlay = melanomaTMA, colorBy = "segment", 
                          corner = "bottomright")
 legend <- fluorLegend(melanomaTMA, nrow = 2, textSize = 4, 
                       boxColor = "grey85", alpha = 0.3)
 cowplot::ggdraw() +
   cowplot::draw_plot(gp) +
   cowplot::draw_plot(legend, scale = 0.105, x = 0.1, y = -0.25)
 ggsave("image_example.png", dpi = 1000)
 
melanomaTMA <- t(melanomaTMA)
plotSpatialOverlay(overlay = melanomaTMA, colorBy = "segment", scaleBar = FALSE)
# 
# melanomaTMA <- cropTissue(overlay = melanomaTMA, buffer = 0.05)
# plotSpatialOverlay(overlay = melanomaTMA, colorBy = "ROILabel", legend = FALSE, scaleBar = FALSE)+
#   viridis::scale_fill_viridis(option = "C")
# samps <- melanomaAnnots$Sample_ID[melanomaAnnots$segment == "Full ROI" & 
#                                    melanomaAnnots$slide.name == slideName(melanomaTMA)]
#melanomaTMACrop <- cropSamples(overlay = melanomaTMA, sampleIDs = samps, sampsOnly = TRUE)

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "FAP-alpha", scaleBar = TRUE, 
                   corner = "bottomleft", textDistance = 4, scaleBarWidth = 0.2, dpi = 300)+
  ggplot2::scale_fill_gradient2(low = "gray", high = "red", 
                                mid = "yellow", midpoint = 80)+
  ggplot2::labs(title = "FAP-alpha Expression in Melanoma TMA")



