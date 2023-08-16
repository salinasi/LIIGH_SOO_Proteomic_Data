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

tifFile

# add lab worksheet here
melanomaTMA_LW <- "/Users/isabella/SpatialOmics_Isabella/Input_files_analysis_june2023/Downloads_from_DSP/P1001660009598A_LabWorksheet.txt"
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
datadir <- "/Users/isabella/SpatialOmics_Isabella/Input_files_analysis_june2023/RCCfiles/20221116_P1001660009598A_RCC"
# reading rcc file data only (no rlf or other files)
rccs <- dir(datadir, pattern = ".RCC", full.names = TRUE)
melanomaDataSet <- readNanoStringRccSet(rccs)

## change data inputs from annotation sheet here
initial_data_csv <- read_csv("/Users/isabella/SpatialOmics_Isabella/Input_files_analysis_june2023/RCCfiles/20221116_P1001660009598A_RCC/initial_dataset_6071_geomx_dsp_roi_collected_archivo_anotacion_final_060523.csv")
# load in initial dataset from annotation sheet
initial_data_dataframe <- as.data.frame(initial_data_csv)
# convert from csv to data.frame

########### CD34 ###########
# plot by segments
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = melanomaAnnots, 
                             plottingFactor = "segment")
# code stops here: cannot find target names in NanoString RCC Data Set
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                             plottingFactor = "CD34")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                             plottingFactor = "ROILabel")
melanomaTMA

head(plotFactors(melanomaTMA))

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD34", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD34 Expression in Melanoma TMA")

#lowest resolution = fastest speeds
checkValidRes(ometiff = tifFile)

# optimal res recommended by SpatialOmics (resolution lowest to highest: 9 -> 1)
res <- 9
# add image back in to pot on top of
melanomaTMA <- addImageOmeTiff(overlay = melanomaTMA, ometiff = tifFile, res = res)
melanomaTMA

# only image
showImage(melanomaTMA)

# plotSpatialOverlay(overlay = melanomaTMA, colorBy = "segment", corner = "topcenter", 
#                     scaleBarWidth = 0.5, textDistance = 130, scaleBarColor = "cyan",
#                    fluorLegend = TRUE)
# # 
#  gp <- plotSpatialOverlay(overlay = melanomaTMA, colorBy = "segment", 
#                           corner = "bottomright")
#  legend <- fluorLegend(melanomaTMA, nrow = 2, textSize = 4, 
#                        boxColor = "grey85", alpha = 0.3)
#  cowplot::ggdraw() +
#    cowplot::draw_plot(gp) +
#    cowplot::draw_plot(legend, scale = 0.105, x = 0.1, y = -0.25)
#  ggsave("image_example.png", dpi = 1000)
#  
# melanomaTMA <- t(melanomaTMA)
# plotSpatialOverlay(overlay = melanomaTMA, colorBy = "segment", scaleBar = FALSE)

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD34", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 600)+
  ggplot2::labs(title = "CD34 Expression in Melanoma TMA")

ggsave("CD34_expression.png", dpi= 1000)


samps <- melanomaTMAAnnots$Sample_ID[melanomaTMAAnnots$QC_status == "PASSED" & 
                                   melanomeTMAAnnots$slide.name == slideName(melanomaTMA)]

melanomaTMACrop <- cropSamples(overlay = melanomaTMA, sampleIDs = samps, sampsOnly = TRUE)
plotSpatialOverlay(overlay = melanomaTMACrop, colorBy = "CD34", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 600)
########### FAP-alpha ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "FAP-alpha")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "FAP-alpha", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "FAP-alpha Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "FAP-alpha", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 80)+
  ggplot2::labs(title = "FAP-alpha Expression in Melanoma TMA")

ggsave("FAP_alpha_expression.png", dpi= 1000)



########### CD14 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD14")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD14", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD14 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD14", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 700)+
  ggplot2::labs(title = "CD14 Expression in Melanoma TMA")

ggsave("CD14_expression.png", dpi= 1000)


########### CD45RO ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD45RO")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD45RO", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD45RO Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD45RO", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 350)+
  ggplot2::labs(title = "CD45RO Expression in Melanoma TMA")

ggsave("CD45RO_expression.png", dpi= 1000)


########### CD163 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD163")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD163", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD163 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD163", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 200)+
  ggplot2::labs(title = "CD163 Expression in Melanoma TMA")

ggsave("CD163_expression.png", dpi= 1000)


########### FOXP3 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "FOXP3")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "FOXP3", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "FOXP3 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "FOXP3", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 22)+
  ggplot2::labs(title = "FOXP3 Expression in Melanoma TMA")

ggsave("FOXP3_expression.png", dpi= 1000)


########### CD66b ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD66b")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD66b", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD66b Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD66b", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 350)+
  ggplot2::labs(title = "CD66b Expression in Melanoma TMA")

ggsave("CD66b_expression.png", dpi= 1000)


########### CD40 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD40")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD40", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD40 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD40", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 125)+
  ggplot2::labs(title = "CD40 Expression in Melanoma TMA")

ggsave("CD40_expression.png", dpi= 1000)


########### CD127 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD127")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD127", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD127 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD127", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 1200)+
  ggplot2::labs(title = "CD127 Expression in Melanoma TMA")

ggsave("CD127_expression.png", dpi= 1000)


########### CD25 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD25")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD25", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD25 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD25", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 150)+
  ggplot2::labs(title = "CD25 Expression in Melanoma TMA")

ggsave("CD25_expression.png", dpi= 1000)


########### CD27 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD27")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD27", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD27 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD27", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 60)+
  ggplot2::labs(title = "CD27 Expression in Melanoma TMA")

ggsave("CD27_expression.png", dpi= 1000)


########### PD-L2 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "PD-L2")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "PD-L2", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "PD-L2 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "PD-L2", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 25)+
  ggplot2::labs(title = "PD-L2 Expression in Melanoma TMA")

ggsave("PD_L2_expression.png", dpi= 1000)


########### CD80 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD80")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD80", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD80 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD80", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 20)+
  ggplot2::labs(title = "CD80 Expression in Melanoma TMA")

ggsave("CD80_expression.png", dpi= 1000)


########### ICOS ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "ICOS")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "ICOS", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "ICOS Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "ICOS", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 125)+
  ggplot2::labs(title = "ICOS Expression in Melanoma TMA")

ggsave("ICOS_expression.png", dpi= 1000)


########### CD44 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD44")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD44", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD44 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD44", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 12000)+
  ggplot2::labs(title = "CD44 Expression in Melanoma TMA")

ggsave("CD44_expression.png", dpi= 1000)


########### VISTA ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "VISTA")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "VISTA", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "VISTA Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "VISTA", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 125)+
  ggplot2::labs(title = "VISTA Expression in Melanoma TMA")

ggsave("VISTA_expression.png", dpi= 1000)


########### IDO1 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "IDO1")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "IDO1", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "IDO1 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "IDO1", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 90)+
  ggplot2::labs(title = "IDO1 Expression in Melanoma TMA")

ggsave("IDO1_expression.png", dpi= 1000)


########### OX40L ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "OX40L")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "OX40L", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "OX40L Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "OX40L", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 700)+
  ggplot2::labs(title = "OX40L Expression in Melanoma TMA")

ggsave("OX40L_expression.png", dpi= 1000)


########### Tim-3 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "Tim-3")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "Tim-3", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "Tim-3 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "Tim-3", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 125)+
  ggplot2::labs(title = "Tim-3 Expression in Melanoma TMA")

ggsave("Tim_3_expression.png", dpi= 1000)


########### LAG3 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "LAG3")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "LAG3", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "LAG3 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "LAG3", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 20)+
  ggplot2::labs(title = "LAG3 Expression in Melanoma TMA")

ggsave("LAG3_expression.png", dpi= 1000)


########### STING ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "STING")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "STING", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "STING Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "STING", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 700)+
  ggplot2::labs(title = "STING Expression in Melanoma TMA")

ggsave("STING_expression.png", dpi= 1000)


########### 4-1BB ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "4-1BB")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "4-1BB", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "4-1BB Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "4-1BB", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 55)+
  ggplot2::labs(title = "4-1BB Expression in Melanoma TMA")

ggsave("4_1BB_expression.png", dpi= 1000)


########### GITR ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "GITR")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "GITR", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "GITR Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "GITR", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 30)+
  ggplot2::labs(title = "GITR Expression in Melanoma TMA")

ggsave("GITR_expression.png", dpi= 1000)


########### ARG1 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "ARG1")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "ARG1", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "ARG1 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "ARG1", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 330)+
  ggplot2::labs(title = "ARG1 Expression in Melanoma TMA")

ggsave("ARG1_expression.png", dpi= 1000)


########### B7-H3 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "B7-H3")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "B7-H3", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "B7-H3 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "B7-H3", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 500)+
  ggplot2::labs(title = "B7-H3 Expression in Melanoma TMA")

ggsave("B7_H3_expression.png", dpi= 1000)


########### CD68 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD68")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD68", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD68 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD68", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 700)+
  ggplot2::labs(title = "CD68 Expression in Melanoma TMA")

ggsave("CD68_expression.png", dpi= 1000)


########### Ms IgG2a ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "Ms IgG2a")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "Ms IgG2a", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "Ms IgG2a Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "Ms IgG2a", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 22)+
  ggplot2::labs(title = "Ms IgG2a Expression in Melanoma TMA")

ggsave("Ms_IgG2a_expression.png", dpi= 1000)


########### GAPDH ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "GAPDH")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "GAPDH", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "GAPDH Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "GAPDH", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 13500)+
  ggplot2::labs(title = "GAPDH Expression in Melanoma TMA")

ggsave("GAPDH_expression.png", dpi= 1000)


########### HYB-NEG ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "HYB-NEG")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "HYB-NEG", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "HYB-NEG Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "HYB-NEG", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 3)+
  ggplot2::labs(title = "HYB-NEG Expression in Melanoma TMA")

ggsave("HYB_NEG_expression.png", dpi= 1000)


########### PD-1 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "PD-1")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "PD-1", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "PD-1 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "PD-1", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 80)+
  ggplot2::labs(title = "PD-1 Expression in Melanoma TMA")

ggsave("PD_1_expression.png", dpi= 1000)


########### CD4 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD4")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD4", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD4 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD4", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 400)+
  ggplot2::labs(title = "CD4 Expression in Melanoma TMA")

ggsave("CD4_expression.png", dpi= 1000)


########### Histone H3 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "Histone H3")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "Histone H3", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "Histone H3 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "Histone H3", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 14000)+
  ggplot2::labs(title = "Histone H3 Expression in Melanoma TMA")

ggsave("Histone_H3_expression.png", dpi= 1000)


########### Beta-2-microglobulin ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "Beta-2-microglobulin")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "Beta-2-microglobulin", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "Beta-2-microglobulin Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "Beta-2-microglobulin", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 115)+
  ggplot2::labs(title = "Beta-2-microglobulin Expression in Melanoma TMA")

ggsave("Beta_2_microglobulin_expression.png", dpi= 1000)


########### Ki-67 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "Ki-67")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "Ki-67", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "Ki-67 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "Ki-67", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 115)+
  ggplot2::labs(title = "Ki-67 Expression in Melanoma TMA")

ggsave("Ki_67_expression.png", dpi= 1000)


########### CD20 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD20")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD20", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD20 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD20", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 65)+
  ggplot2::labs(title = "CD20 Expression in Melanoma TMA")

ggsave("CD20_expression.png", dpi= 1000)


########### Fibronectin ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "Fibronectin")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "Fibronectin", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "Fibronectin Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "Fibronectin", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 7000)+
  ggplot2::labs(title = "Fibronectin Expression in Melanoma TMA")

ggsave("Fibronectin_expression.png", dpi= 1000)


########### Rb IgG ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "Rb IgG")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "Rb IgG", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "Rb IgG Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "Rb IgG", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 80)+
  ggplot2::labs(title = "Rb IgG Expression in Melanoma TMA")

ggsave("Rb_IgG_expression.png", dpi= 1000)


########### CD8 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD8")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD8", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD8 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD8", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 1300)+
  ggplot2::labs(title = "CD8 Expression in Melanoma TMA")

ggsave("CD8_expression.png", dpi= 1000)


########### PD-L1 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "PD-L1")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "PD-L1", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "PD-L1 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "PD-L1", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 33)+
  ggplot2::labs(title = "PD-L1 Expression in Melanoma TMA")

ggsave("PD-L1_expression.png", dpi= 1000)


########### PanCk ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "PanCk")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "PanCk", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "PanCk Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "PanCk", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 7500)+
  ggplot2::labs(title = "PanCk Expression in Melanoma TMA")

ggsave("PanCk_expression.png", dpi= 1000)


########### CD3 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD3")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD3", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD3 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD3", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 1000)+
  ggplot2::labs(title = "CD3 Expression in Melanoma TMA")

ggsave("CD3_expression.png", dpi= 1000)


########### CD56 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD56")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD56", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD56 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD56", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 400)+
  ggplot2::labs(title = "CD56 Expression in Melanoma TMA")

ggsave("CD56_expression.png", dpi= 1000)


########### CD11c ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD11c")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD11c", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD11c Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD11c", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 1300)+
  ggplot2::labs(title = "CD11c Expression in Melanoma TMA")

ggsave("CD11c_expression.png", dpi= 1000)


########### SMA ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "SMA")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "SMA", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "SMA Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "SMA", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 12000)+
  ggplot2::labs(title = "SMA Expression in Melanoma TMA")

ggsave("SMA_expression.png", dpi= 1000)


########### S6 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "S6")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "S6", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "S6 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "S6", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 3000)+
  ggplot2::labs(title = "S6 Expression in Melanoma TMA")

ggsave("S6_expression.png", dpi= 1000)


########### Ms IgG1 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "Ms IgG1")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "Ms IgG1", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "Ms IgG1 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "Ms IgG1", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 60)+
  ggplot2::labs(title = "Ms IgG1 Expression in Melanoma TMA")

ggsave("Ms_IgG1_expression.png", dpi= 1000)


########### HLA-DR ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "HLA-DR")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "HLA-DR", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "HLA-DR Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "HLA-DR", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 1300)+
  ggplot2::labs(title = "HLA-DR Expression in Melanoma TMA")

ggsave("HLA_DR_expression.png", dpi= 1000)


########### CTLA4 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CTLA4")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CTLA4", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CTLA4 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CTLA4", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 14)+
  ggplot2::labs(title = "CTLA4 Expression in Melanoma TMA")

ggsave("CTLA4_expression.png", dpi= 1000)


########### GZMB ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "GZMB")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "GZMB", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "GZMB Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "GZMB", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 525)+
  ggplot2::labs(title = "GZMB Expression in Melanoma TMA")

ggsave("GZMB_expression.png", dpi= 1000)


########### HYB-POS ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "HYB-POS")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "HYB-POS", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "HYB-POS Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "HYB-POS", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 1600)+
  ggplot2::labs(title = "HYB-POS Expression in Melanoma TMA")

ggsave("HYB-POS_expression.png", dpi= 1000)


########### CD45 ###########
# plot by target expression level
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = initial_data_dataframe, 
                                 plottingFactor = "CD45")
# plot by adding in ROI labels
melanomaTMA <- addPlottingFactor(overlay = melanomaTMA, annots = 1:length(sampNames(melanomaTMA)), 
                                 plottingFactor = "ROILabel")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = FALSE, colorBy = "CD45", 
                   scaleBarWidth = 0.3, scaleBarColor = "green") +
  viridis::scale_color_viridis()+
  ggplot2::labs(title = "CD45 Expression in Melanoma TMA")

plotSpatialOverlay(overlay = melanomaTMA, hiRes = TRUE, colorBy = "CD45", scaleBar = TRUE, 
                   corner = "bottomright", textDistance = 4, scaleBarWidth = 0.05, dpi = 1000)+
  ggplot2::scale_fill_gradient2(low = "yellow", high = "red", 
                                mid = "orange", midpoint = 4000)+
  ggplot2::labs(title = "CD45 Expression in Melanoma TMA")

ggsave("CD45_expression.png", dpi= 1000)
