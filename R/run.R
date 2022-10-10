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


###############################################################################
# SET WORKSAPCE                                                               #
###############################################################################
FolderRoot = "~/TCP-TR-H-Clus"
FolderScripts = paste(FolderRoot, "/R", sep="")


##################################################################################################
#
##################################################################################################
execute <- function(parameters){

  FolderRoot = "~/TCP-TR-H-Clus"
  FolderScripts = paste(FolderRoot, "/R", sep="")

  if(parameters$Number.Cores  == 0){
    cat("\n\n##################################################################################################")
    cat("\n# Zero is a disallowed value for number_cores. Please choose a value greater than or equal to 1. #")
    cat("\n##################################################################################################\n\n")
  } else {
    cl <- parallel::makeCluster(parameters$Number.Cores)
    doParallel::registerDoParallel(cl)
    print(cl)

    if(parameters$Number.Cores==1){
      cat("\n\n###########################################################")
      cat("\n# RUN: Running Sequentially!                              #")
      cat("\n###########################################################\n\n")
    } else {
      cat("\n\n######################################################################")
      cat("\n# RUN: Running in parallel with ", parameters$Number.Cores, " cores! #")
      cat("\n######################################################################\n\n")
    }
  }

  retorno = list()

  setwd(FolderScripts)
  source("libraries.R")

  setwd(FolderScripts)
  source("utils.R")

  setwd(FolderScripts)
  source("validateSilho.R")

  setwd(FolderScripts)
  source("testSilho.R")

  setwd(FolderScripts)
  source("validateMaF1.R")

  setwd(FolderScripts)
  source("testMaF1.R")

  cat("\n\n################################################################")
  cat("\n# RUN: Get the label space                                     #")
  cat("\n################################################################\n\n")
  timeLabelSpace = system.time(resLS <- label.space(parameters))
  parameters$LabelSpace = resLS

  cat("\n\n################################################")
  cat("\n# CHOOSED                                        #")
  cat("\n#################################################\n\n")
  timeChoosed = system.time(resChoosed <- choosed(parameters))
  parameters$Choosed = resChoosed


  if(validation==1){


    cat("\n\n#############################################################")
    cat("\n# VALIDATION WITH SILHOUETTE                                #")
    cat("\n#############################################################\n\n")
    timeVal = system.time(resTHP <- silhouette(parameters))


    cat("\n\n#############################################################")
    cat("\n# VERIFYING THRESHOLDS                                      #")
    cat("\n#############################################################\n\n")
    timeVeri = system.time(resVeri <- silho.verifying.tresholds(parameters))
    parameters$valid_tr = resVeri$tr_valid


    cat("\n\n#############################################################")
    cat("\n# BEST SILHOUETTE                                             #")
    cat("\n#############################################################\n\n")
    timeBest = system.time(resBest <- silho.best.partitions(parameters))
    parameters$bests = resBest


    cat("\n\n#######################################################")
    cat("\n# COPY VALIDATION TO GOOGLE DRIVE                     #")
    cat("\n#########################################################\n\n")
    origem1 = parameters$Folders$folderValSilho
    destino1 = paste("nuvem:Clus/Communities/Test/",
                     similarity, "/Silhouette/", dataset_name,
                     "/Tr-H/Validation", sep="")
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
      cat("\n# DELETING VALIDATION                                 #")
      cat("\n#########################################################\n\n")
    system(paste("rm -r ", parameters$Folders$folderValSilho, sep=""))


    if(length(parameters$valid_tr)==0){


    } else {

      cat("\n\n############################################################")
      cat("\n# TEST WITH SIHOUETTE                                      #")
      cat("\n############################################################\n\n")
      timeTest = system.time(resTHP <- silho.test(parameters))


      cat("\n\n#######################################################")
      cat("\n# COPY TO GOOGLE DRIVE                                #")
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


      cat("\n\n#############################################################")
      cat("\n# RUN: Save Runtime                                         #")
      cat("\n#############################################################\n\n")
      Runtime = rbind(timeLabelSpace,
                      timeVeri,
                      timeVal,
                      timeBest,
                      timeTest)
      setwd(parameters$Folders$folderTestSilho)
      write.csv(Runtime, paste(dataset_name,
                               "-test-Runtime.csv", sep=""),
                row.names = FALSE)
    }


  } else if (validation==2){

    cat("\n\n#############################################################")
    cat("\n# VALIDATION WITH CLUS MACRO-F1                             #")
    cat("\n#############################################################\n\n")
    timeVal = system.time(resTHP <- maf1.validate(parameters))


    cat("\n\n#############################################################")
    cat("\n# BEST PARTITIONS MACRO-F1                                  #")
    cat("\n#############################################################\n\n")
    parameters$Best = 8
    timeBest = system.time(resTHP <- maf1.best.partitions(parameters))


    cat("\n\n#############################################################")
    cat("\n# RUN COPY VALIDATION TO GOOGLE DRIVE                       #")
    cat("\n#############################################################\n\n")
    origem1 = parameters$Folders$folderValMaF1
    destino1 = paste("nuvem:Clus/Communities/Test/",
                     similarity, "/Macro-F1/", dataset_name,
                     "/Tr-H/Validation", sep="")
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

    cat("\n\n#############################################################")
    cat("\n# DELETING VALIDATION DIR                                     #")
    cat("\n#############################################################\n\n")
    system(paste("rm -r ", parameters$Folders$folderValMaF1, sep=""))


    cat("\n\n#############################################################")
    cat("\n# TEST WITH CLUS MACRO-F1                                   #")
    cat("\n#############################################################\n\n")
    timeTest = system.time(resTHP <- maf1.test(parameters))


    cat("\n\n#############################################################")
    cat("\n# RUN: Save Runtime                                         #")
    cat("\n#############################################################\n\n")
    Runtime = rbind(timeLabelSpace,
                    timeChoosed,
                    timeVal,
                    timeBest,
                    timeTest)
    setwd(parameters$Folders$folderTestMaF1)
    write.csv(Runtime, paste(dataset_name,
                             "-test-Runtime.csv", sep=""),
              row.names = FALSE)

  } else {

    cat("\n\n############################################################")
    cat("\n# VALIDATION WITH CLUS MICRO-F1                            #")
    cat("\n############################################################\n\n")
    timeVal = system.time(resTHP <- validate(resChoosed, parameters))


    cat("\n\n#############################################################")
    cat("\n# BEST PARTITIONS MICRO-F1                                  #")
    cat("\n#############################################################\n\n")
    parameters$Best = 13
    timeBest = system.time(resTHP <- bestPartition(parameters))


    cat("\n\n#############################################################")
    cat("\n# RUN COPY VALIDATION TO GOOGLE DRIVE                       #")
    cat("\n#############################################################\n\n")
    origem1 = parameters$Folders$folderValMiF1
    destino1 = paste("nuvem:Clus/Communities/Test/",
                     similarity, "/Micro-F1/", dataset_name,
                     "/Tr-H/Validation", sep="")
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


    cat("\n\n#############################################################")
    cat("\n# TEST WITH CLUS MICRO-F1                                   #")
    cat("\n#############################################################\n\n")
    timeTest = system.time(resTHP <- test(parameters))


    cat("\n\n#############################################################")
    cat("\n# RUN: Save Runtime                                         #")
    cat("\n#############################################################\n\n")
    Runtime = rbind(timeLabelSpace,
                    timeChoosed,
                    timeVal,
                    timeBest,
                    timeTest)
    setwd(parameters$Folders$folderTestMiF1)
    write.csv(Runtime, paste(dataset_name,
                             "-test-Runtime.csv", sep=""),
              row.names = FALSE)

  }


  cat("\n\n#############################################################")
  cat("\n# RUN: Stop Parallel                                          #")
  cat("\n###############################################################\n\n")
  on.exit(stopCluster(cl))

  cat("\n\n#############################################################")
  cat("\n# RUN: END                                                    #")
  cat("\n###############################################################\n\n")

  gc()

}

##################################################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com                                   #
# Thank you very much!                                                                           #
##################################################################################################
