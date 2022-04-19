##################################################################################################
# Test the Best Partition with CLUS                                                              #
# Copyright (C) 2021                                                                             #
#                                                                                                #
# This code is free software: you can redistribute it and/or modify it under the terms of the    #
# GNU General Public License as published by the Free Software Foundation, either version 3 of   #
# the License, or (at your option) any later version. This code is distributed in the hope       #
# that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    #
# more details.                                                                                  #
#                                                                                                #
# Elaine Cecilia Gatto | Prof. Dr. Ricardo Cerri | Prof. Dr. Mauri Ferrandin                     #
# Federal University of Sao Carlos (UFSCar: https://www2.ufscar.br/) Campus Sao Carlos           #
# Computer Department (DC: https://site.dc.ufscar.br/)                                           #
# Program of Post Graduation in Computer Science (PPG-CC: http://ppgcc.dc.ufscar.br/)            #
# Bioinformatics and Machine Learning Group (BIOMAL: http://www.biomal.ufscar.br/)               #
#                                                                                                #
##################################################################################################



FolderRoot = "~/TCP-TR-H/"
FolderScripts = paste(FolderRoot, "/R/", sep="")



##################################################################################################
# FUNCTION DIRECTORIES                                                                           #
#   Objective:                                                                                   #
#      Creates all the necessary folders for the project.                                        #
#   Parameters:                                                                                  #
#      dataset_name: name of the dataset                                                         #
#      folderResults: path to save process the algorithm. Example: "/dev/shm/birds",             #
#                     "/scratch/birds", "/home/usuario/birds", "/C:/Users/usuario/birds"         #
#   Return:                                                                                      #
#      All path directories                                                                      #
##################################################################################################
directories <- function(dataset_name, folderResults){

  retorno = list()

  #############################################################################
  # RESULTS FOLDER:                                                           #
  # Parameter from command line. This folder will be delete at the end of the #
  # execution. Other folder is used to store definitely the results.          #
  # Example: "/dev/shm/res"                                                   #
  #############################################################################
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
  # "/home/elaine/TestCommunitiesPartitions/utils"                            #
  #############################################################################
  folderUtils = paste(FolderRoot, "/utils", sep="")
  #cat("\n\t", folderUtils)
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
  # "/home/elaine/TestCommunitiesPartitions/datasets"                         #
  #############################################################################
  folderDatasets = paste(folderResults, "/datasets", sep="")
  #cat("\n\t", folderDatasets)
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
  # "/home/elaine/TestCommunitiesPartitions/datasets/GpositiveGO/LabelSpace"  #
  #############################################################################
  folderLabelSpace = paste(folderDatasets, "/LabelSpace", sep="")
  #cat("\n\t", folderLabelSpace)
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
  # "/home/elaine/TestCommunitiesPartitions/datasets/GpositiveGO/NamesLabels" #
  #############################################################################
  folderNamesLabels = paste(folderDatasets, "/NamesLabels", sep="")
  #cat("\n\t", folderNamesLabels)
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
  # "/home/elaine/TestCommunitiesPartitions/datasets/GpositiveGO/CrossValidation"
  #############################################################################
  folderCV = paste(folderDatasets, "/CrossValidation", sep="")
  #cat("\n\t", folderCV)
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
  # "/home/elaine/TestCommunitiesPartitions/datasets/GpositiveGO/CrossValidation/Tr"
  #############################################################################
  folderCVTR = paste(folderCV, "/Tr", sep="")
  #cat("\n\t", folderCVTR)
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
  # "/home/elaine/TestCommunitiesPartitions/datasets/GpositiveGO/CrossValidation/Ts"
  #############################################################################
  folderCVTS = paste(folderCV, "/Ts", sep="")
  #cat("\n\t", folderCVTS)
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
  # "/home/elaine/TestCommunitiesPartitions/datasets/GpositiveGO/CrossValidation/Vl"
  #############################################################################
  folderCVVL = paste(folderCV, "/Vl", sep="")
  #cat("\n\t", folderCVVL)
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
  # PARTITIONS FOLDER:                                                        #
  # Folder that has the partitions previously computed                        #
  # "/home/elaine/TestCommunitiesPartitions/Partitions"                       #
  #############################################################################
  folderPartitions = paste(folderResults, "/Partitions", sep="")
  #cat("\n\t", folderPartitions)
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
  # RESULTS PARTITIONS FOLDER:                                                #
  # Folder to store the results from partitioning the label correlations      #
  # "/dev/shm/t1/GpositiveGO/Partitions"                                      #
  #############################################################################
  folderValidate = paste(folderResults, "/Validation", sep="")
  #cat("\n\t", folderResPartitions)
  if(dir.exists(folderValidate) == TRUE){
    setwd(folderValidate)
    dir_folderValidate= dir(folderValidate)
    n_folderValidate = length(dir_folderValidate)
  } else {
    dir.create(folderValidate)
    setwd(folderValidate)
    dir_folderValidate= dir(folderValidate)
    n_folderValidate = length(dir_folderValidate)
  }

  #############################################################################
  # REPORTS FOLDER:                                                           #
  # Folder to store all results from this code in the ROOT FOLDER             #
  # The folder REPORTS created before can be anywhere you want, and it will   #
  # be delete in the end of the execution of this code. So, this folder here  #
  # is meant to store definitely the results to be used after to analyse.     #
  # "/home/elaine/TestCommunitiesPartitions/Reports"                          #
  #############################################################################
  folderTest = paste(folderResults, "/Test", sep="")
  #cat("\n\t", folderReports)
  if(dir.exists(folderTest) == TRUE){
    setwd(folderTest)
    dir_folderTest = dir(folderTest)
    n_folderTest = length(dir_folderTest)
  } else {
    dir.create(folderTest)
    setwd(folderTest)
    dir_folderTest = dir(folderTest)
    n_folderTest = length(dir_folderTest)
  }


  folderCommunities = paste(folderResults, "/Communities", sep="")
  #cat("\n\t", folderPartitions)
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


  folderReports = paste(folderResults, "/Reports", sep="")
  #cat("\n\t", folderPartitions)
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
  # RETURN ALL PATHS                                                          #
  #############################################################################
  retorno$folderResults = folderResults
  retorno$folderUtils = folderUtils
  retorno$folderDatasets = folderDatasets
  retorno$folderLabelSpace = folderLabelSpace
  retorno$folderNamesLabels = folderNamesLabels
  retorno$folderCV = folderCV
  retorno$folderCVTR = folderCVTR
  retorno$folderCVTS = folderCVTS
  retorno$folderCVVL = folderCVVL
  retorno$folderPartitions = folderPartitions
  retorno$folderValidate = folderValidate
  retorno$folderTest = folderTest
  retorno$folderCommunities = folderCommunities
  retorno$folderReports = folderReports


  #############################################################################
  # RETURN ALL DIRS                                                           #
  #############################################################################
  retorno$dir_folderResults = dir_folderResults
  retorno$dir_folderUtils = dir_folderUtils
  retorno$dir_folderDatasets = dir_folderDatasets
  retorno$dir_folderLabelSpace = dir_folderLabelSpace
  retorno$dir_folderNamesLabels = dir_folderNamesLabels
  retorno$dir_folderCV = dir_folderCV
  retorno$dir_folderCVTR = dir_folderCVTR
  retorno$dir_folderCVTS = dir_folderCVTS
  retorno$dir_folderCVVL = dir_folderCVVL
  retorno$dir_folderPartitions = dir_folderPartitions
  retorno$dir_folderValidate = dir_folderValidate
  retorno$dir_folderTest = dir_folderTest


  #############################################################################
  # RETURN ALL LENGHTS                                                        #
  #############################################################################
  retorno$n_folderResults = n_folderResults
  retorno$n_folderUtils = n_folderUtils
  retorno$n_folderDatasets = n_folderDatasets
  retorno$n_folderLabelSpace = n_folderLabelSpace
  retorno$n_folderNamesLabels = n_folderNamesLabels
  retorno$n_folderCV = n_folderCV
  retorno$n_folderCVTR = n_folderCVTR
  retorno$n_folderCVTS = n_folderCVTS
  retorno$n_folderCVVL = n_folderCVVL
  retorno$n_folderPartitions = n_folderPartitions
  retorno$n_folderValidate = n_folderValidate
  retorno$n_folderTest = n_folderTest

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
labelSpace <- function(ds, dataset_name, number_folds, folderResults){

  retorno = list()

  # return all fold label space
  classes = list()

  # get the directories
  diretorios = directories(dataset_name, folderResults)

  # from the first FOLD to the last
  k = 1
  while(k<=number_folds){

    # cat("\n\tFold: ", k)

    # enter folder train
    setwd(diretorios$folderCVTR)

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
  cat("\n##################################################################################################")
  cat("\n# FUNCTION LABEL SPACE: END                                                                      #")
  cat("\n##################################################################################################")
  cat("\n\n\n\n")
}


################################################################################################
# FUNCTION INFO DATA SET                                                                       #
#  Objective                                                                                   #
#     Gets the information that is in the "datasets-hpmlk.csv" file.                           #
#  Parameters                                                                                  #
#     dataset: the specific dataset                                                            #
#  Return                                                                                      #
#     Everything in the "datasets-hpmlk.csv" file.                                             #
################################################################################################
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


##################################################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com                                   #
# Thank you very much!                                                                           #
##################################################################################################
