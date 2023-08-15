datadir <- "/Users/isabella/SpatialOmics_Isabella/Input_files_analysis_june2023/RCCfiles/20221116_P1001660009598A_RCC"
# datadir <- read_rcc(datadir)
rccs <- dir(datadir, pattern = ".RCC", full.names = TRUE)

# Just RCC data
solidTumorNoRlfPheno <- readNanoStringRccSet(rccs)
varLabels(solidTumorNoRlfPheno)
fvarLabels(solidTumorNoRlfPheno)

pheno <- file.path(datadir, "Initial_Dataset.csv")
solidTumor <-
  readNanoStringRccSet(rccs, rlfFile = NULL, phenoDataFile = pheno)



# read in RCC files
datadir <- "/Users/isabella/SpatialOmics_Isabella/Input_files_analysis_june2023/RCCfiles/20221116_P1001660009598A_RCC"
# datadir <- read_rcc(datadir)
# reading rcc file data only (no rlf or other files)
rccs <- dir(datadir, pattern = ".RCC", full.names = TRUE)
# read in RLF file
rlf_file <- file.path(datadir, "DSP_v1.0.rlf")
rlfInitialData <- readRlfFile(rlf_file)
# read in sample annotation
initial_data_annotation <- file.path(datadir, "initial_dataset_6071_geomx_dsp_roi_collected_archivo_anotacion_final_060523.csv")
# load all the input files into demoData (S4 object)
initialData <- readNanoStringRccSet(rccs, rlfFile = rlf_file, 
                                 phenoDataFile = initial_data_annotation)
solidTumorNoPheno <- readNanoStringRccSet(rccs, rlfFile = rlf_file)
