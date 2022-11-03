rm(list=ls())

cat("\n\n##############################################################")
  cat("\n# START TCP-TR-H                                             #")
  cat("\n##############################################################\n\n")

###############################################################################
# TEST COMMUNITIES PARTITIONS
# Copyright (C) 2022
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


cat("\n\n##############################################################")
cat("\n# TCP-TR-H-Clus SET WORK SPACE                                  #")
cat("\n##############################################################\n\n")
FolderRoot = "~/TCP-TR-H-Clus"
FolderScripts = paste(FolderRoot, "/R", sep="")


cat("\n\n##############################################################")
cat("\n# TCP-TR-H-Clus LOAD SOURCES                                    #")
cat("\n##############################################################\n\n")
setwd(FolderScripts)
source("libraries.R")

setwd(FolderScripts)
source("utils.R")

setwd(FolderScripts)
source("run.R")


cat("\n\n##############################################################")
cat("\n# TCP-TR-H-Clus OPTIONS CONFIGURATIONS                          #")
cat("\n##############################################################\n\n")
options(java.parameters = "-Xmx64g")
options(show.error.messages = TRUE)
options(scipen=20)


cat("\n\n##############################################################")
cat("\n# TCP-TR-H-Clus READ DATASETS                                   #")
cat("\n##############################################################\n\n")
setwd(FolderRoot)
datasets <- data.frame(read.csv("datasets-original.csv"))


cat("\n\n##############################################################")
cat("\n# TCP-TR-H-Clus GET THE ARGUMENTS COMMAND LINE                  #")
cat("\n##############################################################\n\n")
args <- commandArgs(TRUE)



#############################################################################
# FIRST ARGUMENT: getting specific dataset information being processed      #
# from csv file                                                             #
#############################################################################


# config_file = "/home/cissa/TCP-TR-H-Clus/config-files-2/jaccard/Macro-F1/ctjma-GpositiveGO.csv"
# config_file = "/home/cissa/TCP-TR-H-Clus/config-files-2/rogers/Macro-F1/ctroma-GpositiveGO.csv"

# config_file = "/home/cissa/TCP-TR-H-Clus/config-files-2/jaccard/Silhouete/ctjs-GpositiveGO.csv"
# config_file = "/home/cissa/TCP-TR-H-Clus/config-files-2/rogers/Silhouette/ctros-GpositiveGO.csv"



config_file <- args[1]


if(file.exists(config_file)==FALSE){
  cat("\n################################################################")
  cat("\n# Missing Config File! Verify the following path:              #")
  cat("\n# ", config_file, "                                            #")
  cat("\n################################################################\n\n")
  break
} else {
  cat("\n########################################")
  cat("\n# Properly loaded configuration file!  #")
  cat("\n########################################\n\n")
}


cat("\n########################################")
cat("\n# PARAMETERS READ                    #\n")
config = data.frame(read.csv(config_file))
print(config)
cat("\n########################################\n\n")

parameters = list()

# DATASET_PATH
dataset_path = toString(config$Value[1])
dataset_path = str_remove(dataset_path, pattern = " ")
parameters$Path.Dataset = dataset_path

# TEMPORARTY_PATH
folderResults = toString(config$Value[2])
folderResults = str_remove(folderResults, pattern = " ")
parameters$Folder.Results = folderResults

# PARTITIONS_PATH
Partitions_Path = toString(config$Value[3])
Partitions_Path = str_remove(Partitions_Path, pattern = " ")
parameters$Path.Partitions = Partitions_Path

# VALIDATION
validation = as.numeric(config$Value[4])
parameters$Validation = validation

# SIMILARITY
similarity = toString(config$Value[5])
similarity = str_remove(similarity, pattern = " ")
parameters$Similarity = similarity

# DATASET_NAME
dataset_name = toString(config$Value[6])
dataset_name = str_remove(dataset_name, pattern = " ")
parameters$Dataset.Name = dataset_name

# DATASET_NAME
number_dataset = as.numeric(config$Value[7])
parameters$Number.Dataset = number_dataset

# NUMBER_FOLDS
number_folds = as.numeric(config$Value[8])
parameters$Number.Folds = number_folds

# NUMBER_CORES
number_cores = as.numeric(config$Value[9])
parameters$Number.Cores = number_cores

# DATASET_INFO
ds = datasets[number_dataset,]
parameters$Dataset.Info = ds


cat("\n################################################################\n")
print(ds)
cat("\n# DATASET PATH: \t", dataset_path)
cat("\n# TEMPORARY PATH: \t", folderResults)
cat("\n# PARTITIONS PATH: \t", Partitions_Path)
cat("\n# VALIDATION: \t", validation)
cat("\n# SIMILARITY:  \t", similarity)
cat("\n# DATASET NAME:  \t", dataset_name)
cat("\n# NUMBER DATASET: \t", number_dataset)
cat("\n# NUMBER X-FOLDS CROSS-VALIDATION: \t", number_folds)
cat("\n# NUMBER CORES: \t", number_cores)
cat("\n################################################################\n\n")


###############################################################################
# Creating temporary processing folder                                        #
###############################################################################
if(dir.exists(folderResults) == FALSE) {dir.create(folderResults)}



###############################################################################
# Creating all directories that will be needed for code processing            #
###############################################################################
cat("\n######################")
cat("\n# Get directories    #")
cat("\n######################\n")
diretorios <- directories(dataset_name, folderResults, similarity)
print(diretorios)
cat("\n\n")



###############################################################################
# Copying datasets from ROOT folder on server                                 #
###############################################################################

cat("\n####################################################################")
cat("\n# Checking the dataset tar.gz file                                 #")
cat("\n####################################################################\n\n")
str00 = paste(dataset_path, "/", ds$Name,".tar.gz", sep = "")
str00 = str_remove(str00, pattern = " ")

if(file.exists(str00)==FALSE){

  cat("\n######################################################################")
  cat("\n# The tar.gz file for the dataset to be processed does not exist!    #")
  cat("\n# Please pass the path of the tar.gz file in the configuration file! #")
  cat("\n# The path entered was: ", str00, "                                  #")
  cat("\n######################################################################\n\n")
  break

} else {

  # COPIANDO
  str01 = paste("cp ", str00, " ", diretorios$folderDatasets, sep = "")
  res = system(str01)
  if (res != 0) {
    cat("\nError: ", str01)
    break
  }

  # DESCOMPACTANDO
  str02 = paste("tar xzf ", diretorios$folderDatasets, "/", ds$Name,
                ".tar.gz -C ", diretorios$folderDatasets, sep = "")
  res = system(str02)
  if (res != 0) {
    cat("\nError: ", str02)
    break
  }

  str29 = paste("cp -r ", diretorios$folderDatasets, "/", ds$Name,
                "/CrossValidation/* ", diretorios$folderResults,
                "/datasets/CrossValidation/", sep="")
  res=system(str29)
  #if(res!=0){break}else{cat("\ncopiou")}

  str30 = paste("cp -r ",diretorios$folderDatasets, "/", ds$Name,
                "/LabelSpace/* ", diretorios$folderResults,
                "/datasets/LabelSpace/", sep="")
  res=system(str30)
  #if(res!=0){break}else{cat("\ncopiou")}

  str31 = paste("cp -r ", diretorios$folderDatasets, "/", ds$Name,
                "/NamesLabels/* ", diretorios$folderResults,
                "/datasets/NamesLabels/", sep="")
  res=system(str31)
  #if(res!=0){break}else{cat("\ncopiou")}

  str32 = paste("rm -r ", diretorios$folderResults,
                "/datasets/", ds$Name, sep="")
  print(system(str32))
  #if(res!=0){break}else{cat("\napagou")}

  #APAGANDO
  str03 = paste("rm ", diretorios$folderDatasets, "/", ds$Name,
                ".tar.gz", sep = "")
  res = system(str03)

  cat("\n####################################################################")
  cat("\n# tar.gz file of the DATASET loaded correctly!                     #")
  cat("\n####################################################################\n\n")


}



###############################################################################
# Copying PARTITIONS from ROOT folder on server                               #
###############################################################################

cat("\n####################################################################")
cat("\n# Checking the PARTITIONS tar.gz file                              #")
cat("\n####################################################################\n\n")
str00 = paste(Partitions_Path, "/", ds$Name,".tar.gz", sep = "")
str00 = str_remove(str00, pattern = " ")

if(file.exists(str00)==FALSE){

  cat("\n######################################################################")
  cat("\n# The tar.gz file for the dataset to be processed does not exist!    #")
  cat("\n# Please pass the path of the tar.gz file in the configuration file! #")
  cat("\n# The path entered was: ", str00, "                                  #")
  cat("\n######################################################################\n\n")
  break

} else {

  cat("\n####################################################################")
  cat("\n# tar.gz file of the PARTITION loaded correctly!                   #")
  cat("\n####################################################################\n\n")

  # COPIANDO
  str01 = paste("cp ", str00, " ", diretorios$folderPartitions, sep = "")
  res = system(str01)
  if (res != 0) {
    cat("\nError: ", str01)
    break
  }

  # DESCOMPACTANDO
  str02 = paste("tar xzf ", diretorios$folderPartitions, "/", ds$Name,
                ".tar.gz -C ", diretorios$folderPartitions, sep = "")
  res = system(str02)
  if (res != 0) {
    cat("\nError: ", str02)
    break
  }

  #APAGANDO
  str03 = paste("rm ", diretorios$folderPartitions, "/", ds$Name,
                ".tar.gz", sep = "")
  res = system(str03)
  if (res != 0) {
    cat("\nError: ", str03)
    break
  }

  str31 = paste("cp -r ", diretorios$folderPartitions, "/", ds$Name,
                "/Partitions/* ", diretorios$folderResults,
                "/Partitions/", sep="")
  res=system(str31)

  str31 = paste("cp -r ", diretorios$folderPartitions, "/", ds$Name,
                "/Communities/* ", diretorios$folderResults,
                "/Communities/", sep="")
  res=system(str31)

  system(paste("rm -r ", diretorios$folderPartitions, "/", ds$Name, sep=""))

}


#####################################
parameters$Folders = diretorios
#####################################


cat("\n\n###################################################################")
cat("\n# EXECUTE TCP-TR-H-Clus                                              #")
cat("\n#####################################################################\n\n")
timeTCP = system.time(res <- execute(parameters))


cat("\n\n###################################################################")
cat("\n# TCP-TR-H-Clus SAVE RUNTIME                                     #")
cat("\n#####################################################################\n\n")
result_set <- t(data.matrix(timeTCP))


if(parameters$Validation==1){

  setwd(diretorios$folderTestSilho)
  write.csv(result_set, "Runtime.csv")
  print(timeTCP)
  cat("\n")

  cat("\n\n#######################################################")
  cat("\n# COPY TEST TO GOOGLE DRIVE                           #")
  cat("\n#########################################################\n\n")
  origem1 = parameters$Folders$folderTestSilho
  destino1 = paste("nuvem:Clus/Communities/Test/",
                   similarity, "/Silhouette/", dataset_name,
                   "/Tr-H/Tested", sep="")
  comando1 = paste("rclone copy ", origem1, " ",
                   destino1, sep="")
  cat("\n\n\n", comando1, "\n\n\n")
  a = print(system(comando1))
  a = as.numeric(a)
  if(a != 0){
    stop("Erro RCLONE")
    quit("yes")
  }
  cat("\n\n")

  cat("\n\n#######################################################")
  cat("\n# COPY TEST TO GOOGLE DRIVE                           #")
  cat("\n#########################################################\n\n")
  origem1 = parameters$Folders$folderReports
  destino1 = paste("nuvem:Clus/Communities/Test/",
                   similarity, "/Silhouette/", dataset_name,
                   "/Tr-H/", sep="")
  comando1 = paste("rclone copy ", origem1, " ",
                   destino1, sep="")
  cat("\n\n\n", comando1, "\n\n\n")
  a = print(system(comando1))
  a = as.numeric(a)
  if(a != 0){
    stop("Erro RCLONE")
    quit("yes")
  }
  cat("\n\n")


  # cat("\n\n######")
  # cat("\n# COPY #")
  # cat("\n########\n\n")
  #
  # str2 = paste("cp -r ", diretorios$folderTestSilho ,
  #              " ", diretorios$folderRepSilho, sep="")
  # print(system(str2))
  #
  # str2 = paste("cp -r ", diretorios$folderValSilho ,
  #              " ", diretorios$folderRepSilho , sep="")
  # print(system(str2))
  #
  # str2 = paste("cp -r ", diretorios$folderReports ,
  #              "/* ", diretorios$folderRepSilho , sep="")
  # print(system(str2))



} else if(parameters$Validation==2){


  setwd(diretorios$folderTestMaF1)
  write.csv(result_set, "Runtime.csv")
  print(timeTCP)
  cat("\n")

  # cat("\n\n#######################################################")
  # cat("\n# COPY TEST TO GOOGLE DRIVE                           #")
  # cat("\n#########################################################\n\n")
  # origem1 = parameters$Folders$folderTestMaF1
  # destino1 = paste("nuvem:Clus/Communities/Test/",
  #                  similarity, "/Macro-F1/", dataset_name,
  #                  "/Tr-H/Tested", sep="")
  # comando1 = paste("rclone copy ", origem1, " ",
  #                  destino1, sep="")
  # cat("\n\n\n", comando1, "\n\n\n")
  # a = print(system(comando1))
  # a = as.numeric(a)
  # if(a != 0){
  #   stop("Erro RCLONE")
  #   quit("yes")
  # }
  # cat("\n\n")

  # cat("\n\n#######################################################")
  # cat("\n# COPY TEST TO GOOGLE DRIVE                           #")
  # cat("\n#########################################################\n\n")
  # origem1 = parameters$Folders$folderReports
  # destino1 = paste("nuvem:Clus/Communities/Test/",
  #                  similarity, "/Macro-F1/", dataset_name,
  #                  "/Tr-H/", sep="")
  # comando1 = paste("rclone copy ", origem1, " ",
  #                  destino1, sep="")
  # cat("\n\n\n", comando1, "\n\n\n")
  # a = print(system(comando1))
  # a = as.numeric(a)
  # if(a != 0){
  #   stop("Erro RCLONE")
  #   quit("yes")
  # }
  # cat("\n\n")


  str = paste(diretorios$folderRepSilho, "/",
              parameters$Dataset.Name, sep="")

  str2 = paste("cp -r ", diretorios$folderTestSiho, " ", str, sep="")
  print(system(str2))

  str2 = paste("cp -r ", diretorios$folderReports , "/* ", str , sep="")
  print(system(str2))



} else {


  setwd(diretorios$folderTestMiF1)
  write.csv(result_set, "Runtime.csv")
  print(timeTCP)
  cat("\n")

  # cat("\n\n#######################################################")
  # cat("\n# COPY TEST TO GOOGLE DRIVE                           #")
  # cat("\n#########################################################\n\n")
  # origem1 = parameters$Folders$folderTestMiF1
  # destino1 = paste("nuvem:Clus/Communities/Test/",
  #                  similarity, "/Micro-F1/", dataset_name,
  #                  "/Tr-H/Tested", sep="")
  # comando1 = paste("rclone copy ", origem1, " ",
  #                  destino1, sep="")
  # cat("\n\n\n", comando1, "\n\n\n")
  # a = print(system(comando1))
  # a = as.numeric(a)
  # if(a != 0){
  #   stop("Erro RCLONE")
  #   quit("yes")
  # }
  # cat("\n\n")

  # cat("\n\n#######################################################")
  # cat("\n# COPY TEST TO GOOGLE DRIVE                           #")
  # cat("\n#########################################################\n\n")
  # origem1 = parameters$Folders$folderReports
  # destino1 = paste("nuvem:Clus/Communities/Test/",
  #                  similarity, "/Micro-F1/", dataset_name,
  #                  "/Tr-H/", sep="")
  # comando1 = paste("rclone copy ", origem1, " ",
  #                  destino1, sep="")
  # cat("\n\n\n", comando1, "\n\n\n")
  # a = print(system(comando1))
  # a = as.numeric(a)
  # if(a != 0){
  #   stop("Erro RCLONE")
  #   quit("yes")
  # }
  # cat("\n\n")



  str = paste(diretorios$folderRepMaF1, "/",
              parameters$Dataset.Name, sep="")

  str2 = paste("cp -r ", diretorios$folderTestMaF1, " ", str, sep="")
  print(system(str2))

  str2 = paste("cp -r ", diretorios$folderReports , "/* ", str , sep="")
  print(system(str2))


}



cat("\n\n##############################################################")
  cat("\n# TCP-TR-H-Clus CLEAN                                            #")
  cat("\n##############################################################\n\n")
str2 = paste("rm -rf ", diretorios$folderResults, sep="")
print(system(str2))
rm(list = ls())
gc()

cat("\n\n##############################################################")
  cat("\n# TCP-TR-H-Clus END                                              #")
  cat("\n##############################################################\n")

rm(list = ls())

##################################################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com                                   #
# Thank you very much!                                                                           #
##################################################################################################
