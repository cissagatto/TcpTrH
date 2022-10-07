rm(list = ls())

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


cat("\n\n################################################################")
cat("\n# TCP-TR-H: SET WORK SPACE                                      #")
cat("\n##############################################################\n\n")
FolderRoot = "~/TCP-TR-H-Clus"
FolderScripts = paste(FolderRoot, "/R", sep="")


###############################################################################
# LOAD LIBRARY/PACKAGE                                                        #
###############################################################################
setwd(FolderScripts)
source("libraries.R")


###############################################################################
# READING DATASET INFORMATION FROM DATASETS-ORIGINAL.CSV                      #
###############################################################################
setwd(FolderRoot)
datasets = data.frame(read.csv("datasets-original.csv"))
n = nrow(datasets)


###############################################################################
# CREATING FOLDER TO SAVE CONFIG FILES                                        #
###############################################################################
# "~/TCP-KNN-H-Clus/tcp-clus-config-files"
FolderCF = paste(FolderRoot, "/config-files", sep="")
if(dir.exists(FolderCF)==FALSE){dir.create(FolderCF)}

similarity.name = c("jaccard", "rogers")
similarity.nick = c("j", "ro")

# 1 = silhouette
# 2 = macro f1
# 3 = micro f1
validation = c(1,2,3)
validation.name = c("Silhouette", "Macro-F1", "Micro-F1")
validation.nick = c("s", "ma","mi")

# da primeira similaridade para a ultima
s = 1
while(s<=length(similarity.name)){

  # "~/TCP-KNN-H-Clus/tcp-clus-config-files/jaccard"
  FolderSimilarity = paste(FolderCF, "/", similarity.name[s], sep="")
  if(dir.exists(FolderSimilarity)==FALSE){dir.create(FolderSimilarity)}

  # da primeira validação para a ultima
  v = 1
  while(v<=length(validation)){

    # "~/TCP-KNN-H-Clus/tcp-clus-config-files/jaccard/1"
    FolderValidation = paste(FolderSimilarity, "/", validation.name[v], sep="")
    if(dir.exists(FolderValidation)==FALSE){dir.create(FolderValidation)}

    d = 1
    while(d<=n){

      # specific dataset
      ds = datasets[d,]

      cat("\n\n===============================================")
      cat("\nSimilarity \t", similarity.name[s])
      cat("\nValidation \t", validation.name[v])
      cat("\nDataset \t", ds$Name)
      cat("\n===============================================")

      # Confi File Name
      # "~/TCP-KNN-H-Clus/tcp-clus-config-files/jaccard/Silhouette/
      # jaccard-3s-bbc1000.csv"
      file_name = paste(FolderValidation, "/tr", similarity.nick[s],
                        "-", ds$Name, ".csv", sep="")

      # Starts building the configuration file
      output.file <- file(file_name, "wb")

      # Config file table header
      write("Config, Value",
            file = output.file, append = TRUE)

      # Absolute path to the folder where the dataset's "tar.gz" is stored

      write("Dataset_Path, /home/u704616/Datasets",
            file = output.file, append = TRUE)

      # write("Dataset_Path, /home/cissa/Datasets",
      #      file = output.file, append = TRUE)

      # job name
      job_name = paste("tr", similarity.nick[s], "-", ds$Name, sep = "")

      fn = paste(similarity.nick[s], validation.nick[v], "tr-",
                 ds$Name, sep = "")

      # folder_name = paste("\"/scratch/", job_name, "\"", sep = "")
      # folder_name = paste("~/Exhaustive-MiF1-ECC/", job_name, sep = "")
      # folder_name = paste("~/tmp/", job_name, sep = "")
      # folder_name = paste("/dev/shm/", fn, sep = "")
      folder_name = paste("/scratch/", job_name, sep = "")

      # Absolute path to the folder where temporary processing will be done.
      # You should use "scratch", "tmp" or "/dev/shm", it will depend on the
      # cluster model where your experiment will be run.
      str1 = paste("Temporary_Path, ", folder_name, sep="")
      write(str1,file = output.file, append = TRUE)

      str = paste("/home/u704616/Partitions/CDM/",
                  similarity.name[s], sep="")

      # str = paste("/home/cissa/Partitions/CDM/",
      #            similarity.name[s], sep="")

      str2 = paste("Partitions_Path, ", str,  sep="")
      write(str2, file = output.file, append = TRUE)

      str4 = paste("Validation, ", validation[v], sep="")
      write(str4, file = output.file, append = TRUE)

      str4 = paste("Similarity, ", similarity.name[s], sep="")
      write(str4, file = output.file, append = TRUE)

      # dataset name
      str3 = paste("Dataset_Name, ", ds$Name, sep="")
      write(str3, file = output.file, append = TRUE)

      # Dataset number according to "datasets-original.csv" file
      str2 = paste("Number_Dataset, ", ds$Id, sep="")
      write(str2, file = output.file, append = TRUE)

      # Number used for X-Fold Cross-Validation
      write("Number_Folds, 10", file = output.file, append = TRUE)

      # Number of cores to use for parallel processing
      write("Number_Cores, 10", file = output.file, append = TRUE)

      # finish writing to the configuration file
      close(output.file)

      # increment
      d = d + 1
      gc()
    } # fim do dataset

    v = v + 1
    gc()
  } # fim da validação

  s = s + 1
  gc()
} # fim da similaridade

rm(list = ls())

###############################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com                #
# Thank you very much!                                                        #                                #
###############################################################################
