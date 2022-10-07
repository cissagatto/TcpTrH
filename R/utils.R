###############################################################################
# TEST COMMUNITIES PARTITIONS
# Copyright (C) 2021
#
# This code is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. This code is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
# Public License for more details.
#
# Elaine Cecilia Gatto | Prof. Dr. Ricardo Cerri | Prof. Dr. Mauri Ferrandin
# Federal University of Sao Carlos (UFSCar: https://www2.ufscar.br/) Campus
# Sao Carlos Computer Department (DC: https://site.dc.ufscar.br/)
# Program of Post Graduation in Computer Science
# (PPG-CC: http://ppgcc.dc.ufscar.br/)
# Bioinformatics and Machine Learning Group
# (BIOMAL: http://www.biomal.ufscar.br/)
#
###############################################################################


###############################################################################
# SET WORKSAPCE                                                               #
###############################################################################
FolderRoot = "~/TCP-TR-H-Clus"
FolderScripts = paste(FolderRoot, "/R", sep="")



#########################################################################
# FUNCTION DIRECTORIES
#   Objective:
#      Creates all the necessary folders for the project.
#   Parameters:
#      dataset_name: name of the dataset
#      folderResults: path to save process the algorithm. Example:
#                     "/dev/shm/birds", "/scratch/birds",
#                     "/home/usuario/birds", "/C:/Users/usuario/birds"
#   Return:
#      All path directories
######################################################################
directories <- function(dataset_name, folderResults, similarity){

  retorno = list()

  #############################################################################
  # RESULTS FOLDER:                                                           #
  # Parameter from command line. This folder will be delete at the end of the #
  # execution. Other folder is used to store definitely the results.          #
  # Example: "/dev/shm/j-GpositiveGO"                                         #
  #############################################################################
  retorno$folderResults = folderResults
  if(dir.exists(folderResults) == TRUE){
    setwd(folderResults)
    dir_folderResults = dir(folderResults)
    n_folderResults = length(dir_folderResults)
  } else {
    dir.create(folderResults)
    setwd(folderResults)
    dir_folderResults = dir(folderResults)
    n_folderResults = length(dir_folderResults)
  }
  #cat("\n\t", folderResults)


  #############################################################################
  # UTILS FOLDER:                                                             #
  # Get information about the files within folder utils that already exists   #
  # in the project. It's needed to run CLUS framework and convert CSV files   #
  # in ARFF files correctly.                                                  #
  # "~/TCP-KNN-H-Clus/utils"                                                  #
  #############################################################################
  folderUtils = paste(FolderRoot, "/utils", sep="")
  retorno$folderUtils = folderUtils
  if(dir.exists(folderUtils) == TRUE){
    setwd(folderUtils)
    dir_folderUtils = dir(folderUtils)
    n_folderUtils = length(dir_folderUtils)
  } else {
    dir.create(folderUtils)
    setwd(folderUtils)
    dir_folderUtils = dir(folderUtils)
    n_folderUtils = length(dir_folderUtils)
  }


  #############################################################################
  # DATASETS FOLDER:                                                          #
  # Get the information within DATASETS folder that already exists in the     #
  # project. This folder store the files from cross-validation and will be    #
  # use to get the label space to modeling the label correlations and         #
  # compute silhouete to choose the best hybrid partition.                    #
  # "/dev/shm/j-GpositiveGO/datasets"                                         #
  #############################################################################
  folderDatasets = paste(folderResults, "/datasets", sep="")
  retorno$folderDatasets = folderDatasets
  if(dir.exists(folderDatasets) == TRUE){
    setwd(folderDatasets)
    dir_folderDatasets = dir(folderDatasets)
    n_folderDatasets = length(dir_folderDatasets)
  } else {
    dir.create(folderDatasets)
    setwd(folderDatasets)
    dir_folderDatasets = dir(folderDatasets)
    n_folderDatasets = length(dir_folderDatasets)
  }

  #############################################################################
  # LABEL SPACE FOLDER:                                                       #
  # Path to the specific label space from the dataset that is runing.         #
  # This folder store the label space for each FOLD from the cross-validation #
  # which was computed in the Cross-Validation Multi-Label code.              #
  # In this way, we don't need to load the entire dataset into the running    #
  # "/dev/shm/j-GpositiveGO/datasets/LabelSpace"                              #
  #############################################################################
  folderLabelSpace = paste(folderDatasets, "/LabelSpace", sep="")
  retorno$folderLabelSpace = folderLabelSpace
  if(dir.exists(folderLabelSpace) == TRUE){
    setwd(folderLabelSpace)
    dir_folderLabelSpace = dir(folderLabelSpace)
    n_folderLabelSpace = length(dir_folderLabelSpace)
  } else {
    dir.create(folderLabelSpace)
    setwd(folderLabelSpace)
    dir_folderLabelSpace = dir(folderLabelSpace)
    n_folderLabelSpace = length(dir_folderLabelSpace)
  }


  #############################################################################
  # NAMES LABELS FOLDER:                                                      #
  # Get the names of the labels from this dataset. This will be used in the   #
  # code to create the groups for each partition. Is a way to guarantee the   #
  # use of the correct names labels.                                          #
  # "/dev/shm/j-GpositiveGO/datasets/NamesLabels"                             #
  #############################################################################
  folderNamesLabels = paste(folderDatasets, "/NamesLabels", sep="")
  retorno$folderNamesLabels = folderNamesLabels
  if(dir.exists(folderNamesLabels) == TRUE){
    setwd(folderNamesLabels)
    dir_folderNamesLabels = dir(folderNamesLabels)
    n_folderNamesLabels = length(dir_folderNamesLabels)
  } else {
    dir.create(folderNamesLabels)
    setwd(folderNamesLabels)
    dir_folderNamesLabels = dir(folderNamesLabels)
    n_folderNamesLabels = length(dir_folderNamesLabels)
  }


  #############################################################################
  # CROSS VALIDATION FOLDER:                                                  #
  # Path to the folders and files from cross-validation for the specific      #
  # dataset                                                                   #
  # "/dev/shm/j-GpositiveGO/datasets/CrossValidation"                         #
  #############################################################################
  folderCV = paste(folderDatasets, "/CrossValidation", sep="")
  retorno$folderCV = folderCV
  if(dir.exists(folderCV) == TRUE){
    setwd(folderCV)
    dir_folderCV = dir(folderCV)
    n_folderCV = length(dir_folderCV)
  } else {
    dir.create(folderCV)
    setwd(folderCV)
    dir_folderCV = dir(folderCV)
    n_folderCV = length(dir_folderCV)
  }


  #############################################################################
  # TRAIN CROSS VALIDATION FOLDER:                                            #
  # Path to the train files from cross-validation for the specific dataset    #                                                                   #
  # "/dev/shm/j-GpositiveGO/datasets/CrossValidation/Tr"                      #
  #############################################################################
  folderCVTR = paste(folderCV, "/Tr", sep="")
  retorno$folderCVTR = folderCVTR
  if(dir.exists(folderCVTR) == TRUE){
    setwd(folderCVTR)
    dir_folderCVTR = dir(folderCVTR)
    n_folderCVTR = length(dir_folderCVTR)
  } else {
    dir.create(folderCVTR)
    setwd(folderCVTR)
    dir_folderCVTR = dir(folderCVTR)
    n_folderCVTR = length(dir_folderCVTR)
  }


  #############################################################################
  # TEST CROSS VALIDATION FOLDER:                                             #
  # Path to the test files from cross-validation for the specific dataset     #                                                                   #
  # "/dev/shm/j-GpositiveGO/datasets/CrossValidation/Ts"                      #
  #############################################################################
  folderCVTS = paste(folderCV, "/Ts", sep="")
  retorno$folderCVTS = folderCVTS
  if(dir.exists(folderCVTS) == TRUE){
    setwd(folderCVTS)
    dir_folderCVTS = dir(folderCVTS)
    n_folderCVTS = length(dir_folderCVTS)
  } else {
    dir.create(folderCVTS)
    setwd(folderCVTS)
    dir_folderCVTS = dir(folderCVTS)
    n_folderCVTS = length(dir_folderCVTS)
  }


  #############################################################################
  # VALIDATION CROSS VALIDATION FOLDER:                                       #
  # Path to the validation files from cross-validation for the specific       #
  # dataset                                                                   #
  # "/dev/shm/j-GpositiveGO/datasets/CrossValidation/Vl"                      #
  #############################################################################
  folderCVVL = paste(folderCV, "/Vl", sep="")
  retorno$folderCVVL = folderCVVL
  if(dir.exists(folderCVVL) == TRUE){
    setwd(folderCVVL)
    dir_folderCVVL = dir(folderCVVL)
    n_folderCVVL = length(dir_folderCVVL)
  } else {
    dir.create(folderCVVL)
    setwd(folderCVVL)
    dir_folderCVVL = dir(folderCVVL)
    n_folderCVVL = length(dir_folderCVVL)
  }


  #############################################################################
  # FOLDER PARTITIONS                                                         #
  # "/dev/shm/j-GpositiveGO/Partitions"                                       #
  #############################################################################
  folderPartitions = paste(folderResults, "/Partitions", sep="")
  retorno$folderPartitions = folderPartitions
  if(dir.exists(folderPartitions) == TRUE){
    setwd(folderPartitions)
    dir_folderPartitions = dir(folderPartitions)
    n_folderPartitions = length(dir_folderPartitions)
  } else {
    dir.create(folderPartitions)
    setwd(folderPartitions)
    dir_folderPartitions = dir(folderPartitions)
    n_folderPartitions = length(dir_folderPartitions)
  }


  #############################################################################
  # FOLDER COMMUNITIES                                                        #
  # "/dev/shm/j-GpositiveGO/Communities"                                      #
  #############################################################################
  folderCommunities = paste(folderResults, "/Communities", sep="")
  retorno$folderCommunities = folderCommunities
  if(dir.exists(folderCommunities ) == TRUE){
    setwd(folderCommunities)
    dir_folderCommunities  = dir(folderCommunities )
    n_folderCommunities = length(dir_folderCommunities )
  } else {
    dir.create(folderCommunities )
    setwd(folderCommunities )
    dir_folderCommunities = dir(folderCommunities )
    n_folderCommunities  = length(dir_folderCommunities )
  }

  #############################################################################
  # FOLDER VALIDATION SILHOUETTE                                              #
  # "/dev/shm/j-GpositiveGO/Val-Silho"                                        #
  #############################################################################
  folderValSilho = paste(folderResults, "/Val-Silho", sep="")
  retorno$folderValSilho = folderValSilho
  if(dir.exists(folderValSilho) == TRUE){
    setwd(folderValSilho)
    dir_folderValSilho = dir(folderValSilho)
    n_folderValSilho = length(dir_folderValSilho)
  } else {
    dir.create(folderValSilho)
    setwd(folderValSilho)
    dir_folderValSilho = dir(folderValSilho)
    n_folderValSilho = length(dir_folderValSilho)
  }

  #############################################################################
  # FOLDER VALIDATION MACRO-F1                                                #
  #  "/dev/shm/j-GpositiveGO/Val-MaF1"                                        #
  #############################################################################
  folderValMaF1 = paste(folderResults, "/Val-MaF1", sep="")
  retorno$folderValMaF1 = folderValMaF1
  if(dir.exists(folderValMaF1) == TRUE){
    setwd(folderValMaF1)
    dir_folderValMaF1 = dir(folderValMaF1)
    n_folderValMaF1 = length(dir_folderValMaF1)
  } else {
    dir.create(folderValMaF1)
    setwd(folderValMaF1)
    dir_folderValMaF1 = dir(folderValMaF1)
    n_folderValMaF1 = length(dir_folderValMaF1)
  }

  #############################################################################
  # FOLDER VALIDATION MACRO-F1                                                #
  # "/dev/shm/j-GpositiveGO/Val-MiF1"                                         #
  #############################################################################
  folderValMiF1 = paste(folderResults, "/Val-MiF1", sep="")
  retorno$folderValMiF1 = folderValMiF1
  if(dir.exists(folderValMiF1) == TRUE){
    setwd(folderValMiF1)
    dir_folderValMiF1 = dir(folderValMiF1)
    n_folderValMiF1 = length(dir_folderValMiF1)
  } else {
    dir.create(folderValMiF1)
    setwd(folderValMiF1)
    dir_folderValMiF1 = dir(folderValMiF1)
    n_folderValMiF1 = length(dir_folderValMiF1)
  }

  #############################################################################
  # FOLDER TEST SILHOUETTE                                                    #
  # "/dev/shm/j-GpositiveGO/Test-Silho"                                       #
  #############################################################################
  folderTestSilho = paste(folderResults, "/Test-Silho", sep="")
  retorno$folderTestSilho = folderTestSilho
  if(dir.exists(folderTestSilho) == TRUE){
    setwd(folderTestSilho)
    dir_folderTestSilho = dir(folderTestSilho)
    n_folderTestSilho = length(dir_folderTestSilho)
  } else {
    dir.create(folderTestSilho)
    setwd(folderTestSilho)
    dir_folderTestSilho = dir(folderTestSilho)
    n_folderTestSilho = length(dir_folderTestSilho)
  }


  #############################################################################
  # FOLDER TEST MACRO F1                                                      #
  # "/dev/shm/j-GpositiveGO/Test-MaF1"                                        #
  #############################################################################
  folderTestMaF1 = paste(folderResults, "/Test-MaF1", sep="")
  retorno$folderTestMaF1 = folderTestMaF1
  if(dir.exists(folderTestMaF1) == TRUE){
    setwd(folderTestMaF1)
    dir_folderTestMaF1 = dir(folderTestMaF1)
    n_folderTestMaF1 = length(dir_folderTestMaF1)
  } else {
    dir.create(folderTestMaF1)
    setwd(folderTestMaF1)
    dir_folderTestMaF1 = dir(folderTestMaF1)
    n_folderTestMaF1 = length(dir_folderTestMaF1)
  }


  #############################################################################
  # FOLDER TEST MICRO F1                                                      #
  # "/dev/shm/j-GpositiveGO/Test-MiF1"                                        #
  #############################################################################
  folderTestMiF1 = paste(folderResults, "/Test-MiF1", sep="")
  retorno$folderTestMiF1 = folderTestMiF1
  if(dir.exists(folderTestMiF1) == TRUE){
    setwd(folderTestMiF1)
    dir_folderTestMiF1 = dir(folderTestMiF1)
    n_folderTestMiF1 = length(dir_folderTestMiF1)
  } else {
    dir.create(folderTestMiF1)
    setwd(folderTestMiF1)
    dir_folderTestMiF1 = dir(folderTestMiF1)
    n_folderTestMiF1 = length(dir_folderTestMiF1)
  }

  #############################################################################
  # FOLDER RESULTS REPORTS                                                    #
  # "/dev/shm/j-GpositiveGO/Reports"                                          #
  #############################################################################
  folderReports = paste(folderResults, "/Reports", sep="")
  retorno$folderReports = folderReports
  if(dir.exists(folderReports) == TRUE){
    setwd(folderReports)
    dir_folderReports = dir(folderReports)
    n_folderReports = length(dir_folderReports)
  } else {
    dir.create(folderReports)
    setwd(folderReports)
    dir_folderReports = dir(folderReports)
    n_folderReports = length(dir_folderReports)
  }

  #############################################################################
  # FOLDER ROOT REPORTS                                                       #
  # "~/TCP-KNN-H-Clus/Reports"                                                #
  #############################################################################
  folderRootReports = paste(FolderRoot, "/Reports", sep="")
  retorno$folderRootReports = folderRootReports
  if(dir.exists(folderRootReports) == TRUE){
    setwd(folderRootReports)
    dir_folderRootReports = dir(folderRootReports)
    n_folderRootReports = length(dir_folderRootReports)
  } else {
    dir.create(folderRootReports)
    setwd(folderRootReports)
    dir_folderRootReports = dir(folderRootReports)
    n_folderRootReports = length(dir_folderRootReports)
  }


  #############################################################################
  # FOLDER ROOT REPORTS SIMILARITY                                            #
  # "~/TCP-KNN-H-Clus/Reports/jaccard"                                        #
  #############################################################################
  folderRepSim = paste(folderRootReports, "/", similarity, sep="")
  retorno$folderRepSim = folderRepSim
  if(dir.exists(folderRepSim) == TRUE){
    setwd(folderRepSim)
    dir_folderRepSim = dir(folderRepSim)
    n_folderRepSim = length(dir_folderRepSim)
  } else {
    dir.create(folderRepSim)
    setwd(folderRepSim)
    dir_folderRepSim = dir(folderRepSim)
    n_folderRepSim = length(dir_folderRepSim)
  }

  #############################################################################
  # FOLDER ROOT REPORTS SIMILARITY MACRO F1                                   #
  # "~/TCP-KNN-H-Clus/Reports/jaccard/Silhouette"                             #
  #############################################################################
  folderRepSilho = paste(folderRepSim, "/Silhoutte", sep="")
  retorno$folderRepSilho = folderRepSilho
  if(dir.exists(folderRepSilho) == TRUE){
    setwd(folderRepSilho)
    dir_folderRepSilho = dir(folderRepSilho)
    n_folderRepSilho = length(dir_folderRepSilho)
  } else {
    dir.create(folderRepSilho)
    setwd(folderRepSilho)
    dir_folderRepSilho = dir(folderRepSilho)
    n_folderRepSilho = length(dir_folderRepSilho)
  }

  #############################################################################
  # FOLDER ROOT REPORTS SIMILARITY MACRO F1                                   #
  # "~/TCP-KNN-H-Clus/Reports/jaccard/MacroF1"                                #
  #############################################################################
  folderRepMaF1 = paste(folderRepSim, "/MacroF1", sep="")
  retorno$folderRepMaF1 = folderRepMaF1
  if(dir.exists(folderRepMaF1) == TRUE){
    setwd(folderRepMaF1)
    dir_folderRepMaF1 = dir(folderRepMaF1)
    n_folderRepMaF1 = length(dir_folderRepMaF1)
  } else {
    dir.create(folderRepMaF1)
    setwd(folderRepMaF1)
    dir_folderRepMaF1 = dir(folderRepMaF1)
    n_folderRepMaF1 = length(dir_folderRepMaF1)
  }

  #############################################################################
  # FOLDER ROOT REPORTS SIMILARITY MACRO F1                                   #
  # "~/TCP-KNN-H-Clus/Reports/jaccard/MicroF1"                                #
  #############################################################################
  folderRepMiF1 = paste(folderRepSim, "/MicroF1", sep="")
  retorno$folderRepMiF1 = folderRepMiF1
  if(dir.exists(folderRepMiF1) == TRUE){
    setwd(folderRepMiF1)
    dir_folderRepMiF1 = dir(folderRepMiF1)
    n_folderRepMiF1 = length(dir_folderRepMiF1)
  } else {
    dir.create(folderRepMiF1)
    setwd(folderRepMiF1)
    dir_folderRepMiF1 = dir(folderRepMiF1)
    n_folderRepMiF1 = length(dir_folderRepMiF1)
  }

  return(retorno)
  gc()
}



##################################################################################################
# FUNCTION LABEL SPACE                                                                           #
#   Objective                                                                                    #
#       Separates the label space from the rest of the data to be used as input for              #
#       calculating correlations                                                                 #
#   Parameters                                                                                   #
#       ds: specific dataset information                                                         #
#       dataset_name: dataset name. It is used to save files.                                    #
#       number_folds: number of folds created                                                    #
#       folderResults: folder where to save results                                              #
#   Return:                                                                                      #
#       Training set labels space                                                                #
##################################################################################################
labelSpace <- function(parameters){

  retorno = list()

  # return all fold label space
  classes = list()

  # from the first FOLD to the last
  k = 1
  while(k<=number_folds){

    # cat("\n\tFold: ", k)

    # enter folder train
    setwd(parameters$Folders$folderCVTR)

    # get the correct fold cross-validation
    nome_arquivo = paste(dataset_name, "-Split-Tr-", k, ".csv", sep="")

    # open the file
    arquivo = data.frame(read.csv(nome_arquivo))

    # split label space from input space
    classes[[k]] = arquivo[,ds$LabelStart:ds$LabelEnd]

    # get the names labels
    namesLabels = c(colnames(classes[[k]]))

    # increment FOLD
    k = k + 1

    # garbage collection
    gc()

  } # End While of the 10-folds

  # return results
  retorno$NamesLabels = namesLabels
  retorno$Classes = classes
  return(retorno)

  gc()
  cat("\n################################################################")
  cat("\n# FUNCTION LABEL SPACE: END                                    #")
  cat("\n################################################################")
  cat("\n\n\n\n")
}


###################################################################
# FUNCTION INFO DATA SET
#  Objective
#     Gets the information that is in the "datasets-hpmlk.csv" file.
#  Parameters
#     dataset: the specific dataset
#  Return
#     Everything in the "datasets-hpmlk.csv" file.
####################################################################
infoDataSet <- function(dataset){

  retorno = list()

  retorno$id = dataset$ID
  retorno$name = dataset$Name
  retorno$instances = dataset$Instances
  retorno$inputs = dataset$Inputs
  retorno$labels = dataset$Labels
  retorno$LabelsSets = dataset$LabelsSets
  retorno$single = dataset$Single
  retorno$maxfreq = dataset$MaxFreq
  retorno$card = dataset$Card
  retorno$dens = dataset$Dens
  retorno$mean = dataset$Mean
  retorno$scumble = dataset$Scumble
  retorno$tcs = dataset$TCS
  retorno$attStart = dataset$AttStart
  retorno$attEnd = dataset$AttEnd
  retorno$labStart = dataset$LabelStart
  retorno$labEnd = dataset$LabelEnd
  retorno$distinct = dataset$Distinct
  retorno$xn = dataset$xn
  retorno$yn = dataset$yn
  retorno$gridn = dataset$gridn
  retorno$xt = dataset$xt
  retorno$yt = dataset$yt
  retorno$gridt = dataset$gridt

  return(retorno)

  gc()
}


#############################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com
# Thank you very much!
#############################################################################
