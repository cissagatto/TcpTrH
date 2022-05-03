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
# Runs for all datasets listed in the "datasets.csv" file                                        #
# n_dataset: number of the dataset in the "datasets.csv"                                         #
# number_cores: number of cores to paralell                                                      #
# number_folds: number of folds for cross validation                                             #
# delete: if you want, or not, to delete all folders and files generated                         #
##################################################################################################
execute <- function(ds, dataset_name, number_dataset,
                    number_folds, number_cores,
                    folderResults, diretorios){

  diretorios = diretorios

  if(number_cores == 0){
    cat("\n\n################################################################################################")
    cat("\n#Zero is a disallowed value for number_cores. Please choose a value greater than or equal to 1.#")
    cat("\n##################################################################################################\n\n")
  } else {
    cl <- parallel::makeCluster(number_cores)
    doParallel::registerDoParallel(cl)
    print(cl)

    if(number_cores==1){
      cat("\n\n###############################################################")
      cat("\n# ====> RUN: Running Sequentially!                              #")
      cat("\n#################################################################\n\n")
    } else {
      cat("\n\n###############################################################")
      cat("\n# ====> RUN: Running in parallel with ", number_cores, " cores! #")
      cat("\n#################################################################\n\n")
    }
  }

  retorno = list()


  setwd(FolderScripts)
  source("testClus.R")

  setwd(FolderScripts)
  source("validate.R")

  setwd(FolderScripts)
  source("libraries.R")
  

  cat("\n\n###################################################################")
  cat("\n# ====> RUN: Get dataset information: ", number_dataset, "          #")
  cat("\n#####################################################################\n\n")
  setwd(FolderRoot)
  datasets <- data.frame(read.csv("datasets-2022.csv"))
  ds = datasets[number_dataset,]
  info = infoDataSet(ds)
  dataset_name = toString(ds$Name)



  cat("\n\n###################################################################")
  cat("\n# ====> RUN: Get the label space                                    #")
  cat("\n#####################################################################\n\n")
  timeLabelSpace = system.time(resLS <- labelSpace3(ds,
                                                   dataset_name,
                                                   number_folds,
                                                   folderResults))


  cat("\n\n###################################################################")
  cat("\n# ====> RUN: Compute Silhouete                                      #")
  cat("\n#####################################################################\n\n")
  timeSilho = system.time(resTHP <- comuputeSilhouete3(ds,
                                                      resLS,
                                                      dataset_name,
                                                      number_folds,
                                                      folderResults))


  cat("\n\n###################################################################")
  cat("\n# ====> RUN: VERIFY ELEGIBLE PARTITIONS TO TEST                     #")
  cat("\n#####################################################################\n\n")
  timeVerify = system.time(resVT <- verifyingTresholds3(ds,
                                                       dataset_name,
                                                       number_dataset,
                                                       number_folds,
                                                       number_cores,
                                                       folderResults,
                                                       diretorios))


  cat("\n\n###################################################################")
  cat("\n# ====> TCP-TR-NH: COPY TO GOOGLE DRIVE                             #")
  cat("\n#####################################################################\n\n")
  origem1 = diretorios$folderValidate
  destino1 = paste("cloud:[2022]ResultadosExperimentos/Communities/Test/",
                   similarity, "/", dataset_name, "/TR-H/Reports/Validation", sep="")
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


  valid_tr = resVT$tr_valid

  if(length(valid_tr)==0){

    cat("\n NADA PARA TESTAR!")

  } else {

    cat("\n\n###################################################################")
    cat("\n# ====> RUN: build and test partitions                              #")
    cat("\n#####################################################################\n\n")
    timeBuild = system.time(resBT <- buildAndTest3(ds,
                                                   dataset_name,
                                                   number_dataset,
                                                   number_folds,
                                                   number_cores,
                                                   folderResults,
                                                   diretorios,
                                                   valid_tr))


    cat("\n\n###################################################################")
    cat("\n# ====> RUN: Matrix Confusion                                       #")
    cat("\n#####################################################################\n\n")
    timeSplit = system.time(resGather <- juntaResultados3(ds,
                                                         dataset_name,
                                                         number_dataset,
                                                         number_folds,
                                                         number_cores,
                                                         folderResults,
                                                         diretorios,
                                                         valid_tr))


    cat("\n\n###################################################################")
    cat("\n# ====> RUN: Evaluation Fold                                        #")
    cat("\n#####################################################################\n\n")
    timeAvalia = system.time(resEval <- avaliaTest3(ds,
                                                   dataset_name,
                                                   number_dataset,
                                                   number_folds,
                                                   number_cores,
                                                   folderResults,
                                                   diretorios,
                                                   valid_tr))


    cat("\n\n###################################################################")
    cat("\n# ====> RUN: Gather Evaluation                                      #")
    cat("\n#####################################################################\n\n")
    timeGather = system.time(resGE <- juntaAvaliacoes3(ds,
                                                      dataset_name,
                                                      number_dataset,
                                                      number_folds,
                                                      number_cores,
                                                      folderResults,
                                                      diretorios,
                                                      valid_tr))


    cat("\n\n#################################################################")
    cat("\n# ====> RUN: Organize Evaluation                                  #")
    cat("\n###################################################################\n\n")
    timeOrg = system.time(resGE <- organizaAvaliacoes3(ds,
                                                      dataset_name,
                                                      number_dataset,
                                                      number_folds,
                                                      number_cores,
                                                      folderResults,
                                                      diretorios,
                                                      valid_tr))


  cat("\n\n###################################################################")
  cat("\n# ====> TCP-TR-NH: COPY TO GOOGLE DRIVE                             #")
  cat("\n#####################################################################\n\n")
  origem1 = diretorios$folderTest
  destino1 = paste("cloud:[2022]ResultadosExperimentos/Communities/Test/",
                   similarity, "/", dataset_name, "/TR-H/Reports/Tested", sep="")
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


    cat("\n\n###################################################################")
    cat("\n# ====> RUN: Save Runtime                                           #")
    cat("\n#####################################################################\n\n")
    Runtime = rbind(timeLabelSpace,
                    timeSilho,
                    timeVerify,
                    timeBuild,
                    timeSplit,
                    timeAvalia,
                    timeGather,
                    timeOrg)
    setwd(diretorios$folderTest)
    write.csv(Runtime, paste(dataset_name,
                             "-testPartition-Runtime.csv", sep=""),
              row.names = FALSE)

  }


  cat("\n\n###################################################################")
  cat("\n# ====> RUN: Stop Parallel                                          #")
  cat("\n#####################################################################\n\n")
  on.exit(stopCluster(cl))

  cat("\n\n###################################################################")
  cat("\n# ====> RUN: END                                                    #")
  cat("\n#####################################################################\n\n")


  if(interactive()==TRUE){ flush.console() }
  gc()
}

##################################################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com                                   #
# Thank you very much!                                                                           #
##################################################################################################
