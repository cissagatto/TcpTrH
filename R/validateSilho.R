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



########################################################################
#
########################################################################
label.space <- function(parameters){

  retorno = list()

  # return all fold label space
  classes = list()

  # from the first FOLD to the last
  k = 1
  while(k<=parameters$Number.Folds){

    # cat("\n\tFold: ", k)

    # enter folder train
    setwd(parameters$Folders$folderCVTR)

    # get the correct fold cross-validation
    nome_arquivo = paste(dataset_name, "-Split-Tr-", k, ".csv", sep="")

    # open the file
    arquivo = data.frame(read.csv(nome_arquivo))

    # split label space from input space
    classes[[k]] = arquivo[,parameters$Dataset.Info$LabelStart:parameters$Dataset.Info$LabelEnd]

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



########################################################################
#
########################################################################
silhouette <- function(parameters){

  f = 1
  silhoueteParalel <- foreach(f = 1:parameters$Number.Folds ) %dopar% {
  #while(f<=number_folds){

    cat("\nFold: ", f)

    ############################################################
    FolderRoot = "~/TCP-TR-H-Clus"
    FolderScripts = paste(FolderRoot, "/R", sep="")

    setwd(FolderScripts)
    source("libraries.R")

    setwd(FolderScripts)
    source("utils.R")


    #########################################################################
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
    bestPartition = data.frame(fold, tr, part, maximo,
                               minimo, mediana, media, primeiroQuadrante,
                               terceiroQuadrante, valueSilhouete)


    ########################################################################
    FolderSplitVal = paste(parameters$Folders$folderValSilho,
                           "/Split-", f, sep="")
    if(dir.exists(FolderSplitVal)==FALSE){dir.create(FolderSplitVal)}

    FolderPSplit = paste(parameters$Folders$folderPartitions,
                         "/Split-",f, sep="")

    FolderSplitComm = paste(parameters$Folders$folderCommunities,
                            "/Split-",f, sep="")

    ########################################################################
    # get the space label
    espacoDeRotulos = data.frame(parameters$LabelSpace$Classes[f])
    espacoDeRotulos2 = data.frame(t(espacoDeRotulos))
    labels = rownames(espacoDeRotulos2)
    espacoDeRotulos2 = cbind(labels, espacoDeRotulos2)
    espacoDeRotulos2 = data.frame(espacoDeRotulos2[order(espacoDeRotulos2$labels, decreasing = FALSE),])

    ########################################################################
    cat("\n\nGrupos por particão")
    setwd(FolderPSplit)
    tr_H = data.frame(read.csv(paste("fold-", f,
                                     "-tr-h-choosed.csv", sep="")))
    total_tr_H = nrow(tr_H)

    # do primeiro tr ao ultimo
    u = 0
    while(u<total_tr_H){

      cat("\n#=========================================================")
      cat("\n#Threshold = ", u)
      cat("\n#=========================================================")

      FolderPartTr = paste(FolderPSplit, "/Tr-", u, sep="")
      FolderPartComm = paste(FolderSplitComm, "/Tr-", u, sep="")

      FolderTrVal = paste(FolderSplitVal, "/Tr-", u, sep="")
      if(dir.exists(FolderTrVal)==FALSE){dir.create(FolderTrVal)}

      FolderPSplit = paste(parameters$Folders$folderPartitions,
                           "/Split-",f, sep="")
      setwd(FolderPSplit)
      tr_H = data.frame(read.csv(paste("fold-", f,
                                       "-tr-h-choosed.csv", sep="")))
      total_tr_H = nrow(tr_H)

      a = u + 1
      tr_H = tr_H[a,]

      if(tr_H$method=="none"){
        #cat("METHOD == NONE")

      } else {
        #cat("METHOD != NONE")

        setwd(FolderPartComm)
        particoes = data.frame(read.csv(paste("tr-", u,
                                              "-", tr_H$method,
                                              "-partitions-hierarchical.csv",
                                              sep="")))
        TotalParticoes = ncol(particoes)-1
        numPart = (parameters$Dataset.Info$Labels-1)

        ######################################################################
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
        Silhouete = data.frame(fold, part, tr, maximo, minimo, mediana,
                               media, primeiroQuadrante, terceiroQuadrante,
                               valueSilhouete)

        w = 2
        while(w<=numPart){

          cat("\nPartition ", w)

          FolderPartVal = paste(FolderTrVal, "/Partition-", w, sep="")
          if(dir.exists(FolderPartVal)==FALSE){dir.create(FolderPartVal)}

          ####################################################################
          # get the number of groups for this partition
          particao = particoes[,c(1,w)]
          names(particao) = c("labels", "groups")

          res = constroiParticoes(TotalParticoes)
          res = filter(res, numberPartition == w)
          numGroups = as.numeric(res$numberGroup)

          ###################################################################
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
            Silhouete = rbind(Silhouete, data.frame(fold, part, tr, maximo,
                                                    minimo, mediana, media,
                                                    primeiroQuadrante,
                                                    terceiroQuadrante,
                                                    valueSilhouete))
            setwd(FolderTrVal)
            write.csv(Silhouete[-1,], paste("fold-", f, "-tr-", u,
                                            "-silho.csv", sep=""),
                      row.names = FALSE)

          } else {
            #cat("\ntwo or more labels in the group")

            groups_label_space = cbind(particao, espacoDeRotulos2)
            groups_label_space = groups_label_space[,c(-1,-3)]
            a = dist(groups_label_space)
            b = as.dist(a)
            sil = silhouette(groups_label_space[,1], b)
            sil = sortSilhouette(sil)

            setwd(FolderPartVal)
            write.csv(sil, paste("silho-fold", f, "-tr-", u, "-part-", w,
                                 ".csv", sep=""), row.names = FALSE)

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
              Silhouete = rbind(Silhouete, data.frame(fold, part, tr, maximo,
                                                      minimo, mediana, media,
                                                      primeiroQuadrante,
                                                      terceiroQuadrante,
                                                      valueSilhouete))
              setwd(FolderTrVal)
              write.csv(Silhouete[-1,], paste("fold-", f, "-tr-", u,
                                              "-silho.csv", sep=""),
                        row.names = FALSE)

            } else {

              #cat("\nMore than one label per group\n")

              setwd(FolderPartVal)
              pdf(paste("silho-fold-", f, "-tr-", u, "-part-", w,
                        ".pdf", sep=""), width = 10, height = 8)
              print(plot(sil))
              dev.off()
              cat("\n")

              setwd(FolderPartVal)
              pdf(paste("fviz-silh-fold-", f, "-tr-", u, "-part-", w,
                        ".pdf", sep=""), width = 10, height = 8)
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
              Silhouete = rbind(Silhouete, data.frame(fold, part, tr, maximo,
                                                      minimo, mediana, media,
                                                      primeiroQuadrante,
                                                      terceiroQuadrante,
                                                      valueSilhouete))

              setwd(FolderTrVal)
              write.csv(Silhouete[-1,], paste("fold-", f, "-tr-", u,
                                              "-silho.csv", sep=""), row.names = FALSE)

            } # fim do if

          } # fim do if

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
    write.csv(bestPartition[-1,], paste("fold-", f, "-best-silho.csv", sep=""),
              row.names = FALSE)

    #f = f + 1
    gc()

  } # fim do fold

  gc()
  cat("\n##############################################################")
  cat("\n# END COMPUTE SILHOUETE                                      #")
  cat("\n##############################################################")
  cat("\n\n\n\n")
}


########################################################################
#
########################################################################
silho.best.partitions <- function(parameters){

  retorno = list()
  all.silho = data.frame()

  f = 1
  while(f<=parameters$Number.Folds){
    pasta = paste(parameters$Folders$folderValSilho, "/Split-", f, sep="")
    nome = paste(pasta, "/fold-", f, "-best-silho.csv", sep="")
    arquivo = data.frame(read.csv(nome))
    all.silho  = rbind(all.silho , arquivo)
    f = f + 1
    gc()
  }

  setwd(parameters$Folders$folderReports)
  write.csv(all.silho , "all-best-silho.csv", row.names = FALSE)

  n = length(parameters$valid_tr)
  t = 0
  while(t<n){
    res.1 = data.frame(filter(all.silho, tr==t))
    setwd(parameters$Folders$folderReports)
    write.csv(res.1, paste("tr-",t, "-all-best-silho.csv", sep=""),
              row.names = FALSE)

    res.2 = data.frame(partition = res.1$part,
                       silhouette = res.1$valueSilhouete)

    # computes statistics
    soma = apply(res.2, 2, sum)
    media = apply(res.2, 2, mean)
    mediana = apply(res.2, 2, median)
    desvioPadrao = apply(res.2, 2, sd)
    minimo = apply(res.2, 2, min)
    maximo = apply(res.2, 2, max)
    sumario = rbind(soma, media, mediana, desvioPadrao, minimo, maximo)

    setwd(parameters$Folders$folderReports)
    write.csv(all.silho, paste("tr-",t,
                               "-summary-best-silho.csv", sep=""),
              row.names = FALSE)

    t = t + 1
    gc()
  }

  all.choosed.partitions = data.frame()
  all.choosed.methods = data.frame()
  all.eb.partitions = data.frame()
  all.fg.partitions = data.frame()
  all.wt.partitions = data.frame()

  f = 1
  while(f<=parameters$Number.Folds){
    FolderPartSplit = paste(parameters$Folders$folderPartitions,
                            "/Split-", f, sep="")
    setwd(FolderPartSplit)
    choosed = data.frame(read.csv(paste("fold-", f,
                                        "-tr-h-choosed.csv", sep="")))
    all.choosed.methods = rbind(all.choosed.methods, choosed)

    n = length(parameters$valid_tr)

    t=0
    while(t<n){

      Folder = paste(FolderPartSplit, "/Tr-", t, sep="")
      setwd(Folder)

      nome.eb = paste("tr-", t, "-eb-partitions-hierarchical.csv", sep="")
      choosed.eb = data.frame(read.csv(nome.eb))
      choosed.eb = cbind(fold = f, tr = t, choosed.eb)
      all.eb.partitions = rbind(all.eb.partitions , choosed.eb)

      nome.fg = paste("tr-", t, "-fg-partitions-hierarchical.csv", sep="")
      choosed.fg = data.frame(read.csv(nome.fg))
      choosed.fg = cbind(fold = f, tr = t, choosed.fg)
      all.fg.partitions = rbind(all.fg.partitions , choosed.fg)

      nome.wt = paste("tr-", t, "-wt-partitions-hierarchical.csv", sep="")
      choosed.wt = data.frame(read.csv(nome.wt))
      choosed.wt = cbind(fold = f, tr = t, choosed.wt)
      all.wt.partitions = rbind(all.wt.partitions , choosed.wt)

      all.eb.partitions.2 = cbind(method = "eb", all.eb.partitions)
      all.fg.partitions.2 = cbind(method = "fg", all.fg.partitions)
      all.wt.partitions.2 = cbind(method = "wt", all.wt.partitions)

      t = t + 1
      gc()
    }

    f = f + 1
    gc()
  }


  all.choosed.partitions = rbind(all.eb.partitions.2,
                                 all.fg.partitions.2,
                                 all.wt.partitions.2)

  setwd(parameters$Folders$folderReports)
  write.csv(all.choosed.methods, "all-choosed.csv", row.names = FALSE)
  write.csv(all.eb.partitions, "all-eb-partitions.csv", row.names = FALSE)
  write.csv(all.fg.partitions, "all-fg-partitions.csv", row.names = FALSE)
  write.csv(all.wt.partitions, "all-wt-partitions.csv", row.names = FALSE)
  write.csv(all.choosed.partitions, "all-choosed-partitions.csv", row.names = FALSE)

  retorno$all.silho = all.silho
  retorno$all.eb.partitions = all.eb.partitions
  retorno$all.fg.partitions = all.fg.partitions
  retorno$all.wt.partitions = all.wt.partitions
  retorno$all.choosed.methods = all.choosed.methods
  retorno$all.choosed.partitions = all.choosed.partitions

  return(retorno)

  gc()
  cat("\n###############################################################")
  cat("\n# Statistics: END                                             #")
  cat("\n###############################################################")
  cat("\n\n\n\n")
}



##################################################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com                                   #
# Thank you very much!                                                                           #
##################################################################################################
