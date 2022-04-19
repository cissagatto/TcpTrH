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
labelSpace3 <- function(ds, dataset_name, number_folds, folderResults){

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




##################################################################################################
# FUNCTION COMPUTE SILHOUETE                                                                     #
#   Objective:                                                                                   #
#      Compute sillhouete for each partition/fold                                                #
#   Parameters:                                                                                  #
#      ds: information about the specific dataset                                                #
#      resLS: label space from the specific dataset                                              #
#      number_dataset: number of the specific dataset                                            #
#      number_cores: number of cores to process in paralel                                       #
#      number_folds: number of folds for the cross-validation                                    #
#      folderResults: folder to process                                                          #
#   Return:                                                                                      #
#      Silhouete Graphics                                                                        #
#      Silhouete Values                                                                          #
##################################################################################################
comuputeSilhouete3 <- function (ds, resLS, dataset_name, number_folds, folderResults){

  if(interactive()==TRUE){ flush.console() }

  diretorios = directories(dataset_name, folderResults)

  #cat("\nFrom 1 to 10 folds!")
  f = 1
  silhoueteParalel <- foreach(f = 1:number_folds) %dopar% {
  #while(f<=number_folds){

    cat("\nFold: ", f)

    ############################################################################################################
    FolderRoot = "~/TCP-TR-H/"
    FolderScripts = paste(FolderRoot, "/R/", sep="")

    setwd(FolderScripts)
    source("libraries.R")

    setwd(FolderScripts)
    source("utils.R")


    ##############################################################################
    constroiParticoes <- function(TotalParticoes){

      data <- data.frame(matrix(NA,    # Create empty data frame
                                nrow = TotalParticoes,
                                ncol = 2))

      names(data) = c("numberPartition", "numberGroup")

      i = 1
      a = 1
      while(i<=nrow(data)){
        data[i,1] = a + 1
        data[i,2] = a + 1
        i = i + 1
        a = a + 1
        gc()
      }

      return(data)

    }

    fold = c(0)
    tr = c(0)
    part = c(0)
    maximo = c(0)
    minimo = c(0)
    mediana = c(0)
    media = c(0)
    primeiroQuadrante = c(0)
    terceiroQuadrante = c(0)
    valueSilhouete = c(0)
    bestPartition = data.frame(fold, tr, part, maximo, minimo, mediana, media, primeiroQuadrante, terceiroQuadrante, valueSilhouete)

    ########################################################################
    FolderSplitVal = paste(diretorios$folderValidate, "/Split-", f, sep="")
    if(dir.exists(FolderSplitVal)==FALSE){dir.create(FolderSplitVal)}

    FolderPSplit = paste(diretorios$folderPartitions, "/Split-",f, sep="")
    FolderSplitComm = paste(diretorios$folderCommunities, "/Split-",f, sep="")

    ########################################################################
    # get the space label
    espacoDeRotulos = data.frame(resLS$Classes[f])
    espacoDeRotulos2 = data.frame(t(espacoDeRotulos))
    labels = rownames(espacoDeRotulos2)
    espacoDeRotulos2 = cbind(labels, espacoDeRotulos2)
    espacoDeRotulos2 = data.frame(espacoDeRotulos2[order(espacoDeRotulos2$labels, decreasing = FALSE),])

    ########################################################################
    cat("\n\nGrupos por particão")
    setwd(FolderPSplit)
    tr_H = data.frame(read.csv(paste("fold-", f, "-tr-h-choosed.csv", sep="")))
    total_tr_H = nrow(tr_H)

    # do primeiro tr ao ultimo
    u = 0
    while(u<total_tr_H){

      cat("\n#=========================================================")
      cat("\n#Threshold = ", u)
      cat("\n#=========================================================")

      diretorios = diretorios

      FolderPartTr = paste(FolderPSplit, "/Tr-", u, sep="")
      FolderPartComm = paste(FolderSplitComm, "/Tr-", u, sep="")

      FolderTrVal = paste(FolderSplitVal, "/Tr-", u, sep="")
      if(dir.exists(FolderTrVal)==FALSE){dir.create(FolderTrVal)}

      FolderPSplit = paste(diretorios$folderPartitions, "/Split-",f, sep="")
      setwd(FolderPSplit)
      tr_H = data.frame(read.csv(paste("fold-", f, "-tr-h-choosed.csv", sep="")))
      total_tr_H = nrow(tr_H)

      a = u + 1
      tr_H = tr_H[a,]

      if(tr_H$method=="none"){
        #cat("METHOD == NONE")

      } else {
        #cat("METHOD != NONE")

        setwd(FolderPartComm)
        particoes = data.frame(read.csv(paste("tr-", u, "-", tr_H$method,"-partitions-hierarchical.csv", sep="")))
        TotalParticoes = ncol(particoes)-1
        numPart = (ds$Labels-1)

        ########################################################################
        fold = c(0)
        part = c(0)
        tr = c (0)
        maximo = c(0)
        minimo = c(0)
        mediana = c(0)
        media = c(0)
        primeiroQuadrante = c(0)
        terceiroQuadrante = c(0)
        valueSilhouete = c(0)
        Silhouete = data.frame(fold, part, tr, maximo, minimo, mediana, media,
                               primeiroQuadrante, terceiroQuadrante, valueSilhouete)

        w = 2
        while(w<=numPart){

          cat("\nPartition ", w)

          FolderPartVal = paste(FolderTrVal, "/Partition-", w, sep="")
          if(dir.exists(FolderPartVal)==FALSE){dir.create(FolderPartVal)}

          ########################################################################
          # get the number of groups for this partition
          particao = particoes[,c(1,w)]
          names(particao) = c("labels", "groups")

          res = constroiParticoes(TotalParticoes)
          res = filter(res, numberPartition == w)
          numGroups = as.numeric(res$numberGroup)

          ########################################################################
          if(numGroups==1){
            #cat("\nOnly one group of labels (global partition)")

            fold = f
            part = w
            tr = u
            maximo = NA
            minimo = NA
            mediana = NA
            media = NA
            primeiroQuadrante = NA
            terceiroQuadrante = NA
            valueSilhouete = NA
            Silhouete = rbind(Silhouete, data.frame(fold, part, tr, maximo, minimo, mediana, media,
                                                    primeiroQuadrante, terceiroQuadrante, valueSilhouete))
            setwd(FolderTrVal)
            write.csv(Silhouete[-1,], paste("fold-", f, "-tr-", u, "-silho.csv", sep=""), row.names = FALSE)

          } else {
            #cat("\ntwo or more labels in the group")

            groups_label_space = cbind(particao, espacoDeRotulos2)
            groups_label_space = groups_label_space[,c(-1,-3)]
            a = dist(groups_label_space)
            b = as.dist(a)
            sil = silhouette(groups_label_space[,1], b)
            sil = sortSilhouette(sil)

            setwd(FolderPartVal)
            write.csv(sil, paste("silho-fold", f, "-tr-", u, "-part-", w, ".csv", sep=""), row.names = FALSE)

            if(all(is.na(sil))==TRUE){

              #cat("\nOne label per group (local partition)\n")
              fold = f
              part = w
              tr = u
              maximo = NA
              minimo = NA
              mediana = NA
              media = NA
              primeiroQuadrante = NA
              terceiroQuadrante = NA
              valueSilhouete = NA
              Silhouete = rbind(Silhouete, data.frame(fold, part, tr, maximo, minimo, mediana, media,
                                                      primeiroQuadrante, terceiroQuadrante, valueSilhouete))
              setwd(FolderTrVal)
              write.csv(Silhouete[-1,], paste("fold-", f, "-tr-", u, "-silho.csv", sep=""), row.names = FALSE)

            } else {

              #cat("\nMore than one label per group\n")

              setwd(FolderPartVal)
              pdf(paste("silho-fold-", f, "-tr-", u, "-part-", w, ".pdf", sep=""), width = 10, height = 8)
              print(plot(sil))
              dev.off()
              cat("\n")

              setwd(FolderPartVal)
              pdf(paste("fviz-silh-fold-", f, "-tr-", u, "-part-", w, ".pdf", sep=""), width = 10, height = 8)
              print(fviz_silhouette(sil))
              dev.off()
              cat("\n")

              # Summary of silhouette analysis
              si.sum = summary(sil)
              res.si.sum = unlist(si.sum)

              fold = f
              tr = u
              part = w
              maximo = res.si.sum$si.summary.Max.
              minimo = res.si.sum$si.summary.Min.
              mediana = res.si.sum$si.summary.Median
              media = res.si.sum$si.summary.Mean
              primeiroQuadrante = res.si.sum$`si.summary.1st Qu.`
              terceiroQuadrante = res.si.sum$`si.summary.3rd Qu.`
              valueSilhouete = res.si.sum$avg.width
              Silhouete = rbind(Silhouete, data.frame(fold, part, tr, maximo, minimo, mediana, media,
                                                      primeiroQuadrante, terceiroQuadrante, valueSilhouete))

              setwd(FolderTrVal)
              write.csv(Silhouete[-1,], paste("fold-", f, "-tr-", u, "-silho.csv", sep=""), row.names = FALSE)

            } # fim do if

          } # fim do if

          if(interactive()==TRUE){ flush.console() }
          w = w + 1
          gc()

        } # fim da partição

        Silhouete = Silhouete[-1,]
        indice = as.numeric(which.max(Silhouete$valueSilhouete))
        silhouete2 = Silhouete[indice,]
        bestPartition = rbind(bestPartition, silhouete2)
      }

      u = u + 1
      gc()

    } # fim do tr

    setwd(FolderSplitVal)
    write.csv(bestPartition[-1,], paste("fold-", f, "-best-silho.csv", sep=""), row.names = FALSE)

    #f = f + 1
    if(interactive()==TRUE){ flush.console() }
    gc()

  } # fim do fold

  if(interactive()==TRUE){ flush.console() }
  gc()
  cat("\n##################################################################################################")
  cat("\n# END COMPUTE SILHOUETE                                                                          #")
  cat("\n##################################################################################################")
  cat("\n\n\n\n")
}


##################################################################################################
# FUNCTION ASD                                                                                   #
#   Objective                                                                                    #
#       Compute statistics about the partitions                                                  #
#   Parameters                                                                                   #
#       ds: specific dataset information                                                         #
#       dataset_name: dataset name. It is used to save files.                                    #
#   Return                                                                                       #
#       Sum, mean, median, standart deviation, max and min partitions                            #
##################################################################################################
asd <- function(ds, dataset_name, diretorios, namesLabels){

  # function return
  retorno = list()

  # get the best partitions of the dataset
  pasta = paste(diretorios$folderDatasetReports, "/Validation", sep="")
  setwd(pasta)

  # computes statistics
  soma = apply(particoes, 2, sum)
  media = apply(particoes, 2, mean)
  mediana = apply(particoes, 2, median)
  desvioPadrao = apply(particoes, 2, sd)
  minimo = apply(particoes, 2, min)
  maximo = apply(particoes, 2, max)
  sumario = rbind(soma, media, mediana, desvioPadrao, minimo, maximo)

  # saves results in the RESULTS folder
  setwd(diretorios$folderResultsDataset)
  write.csv(sumario, paste(dataset_name, "-statistic-sumary-best-part.csv", sep=""))

  setwd(Folder)
  write.csv(sumario, paste(dataset_name, "-statistic-sumary-best-part.csv", sep=""))

  # function return
  retorno$sumario = sumario
  return(retorno)

  gc()
  cat("\n##################################################################################################")
  cat("\n# Statistics: END                                                                                #")
  cat("\n##################################################################################################")
  cat("\n\n\n\n")
}





##################################################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com                                   #
# Thank you very much!                                                                           #
##################################################################################################
