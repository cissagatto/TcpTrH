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
choosed <- function(parameters){

  retorno = list()

  # generating partitions and groups
  partitions = seq(1,parameters$Dataset.Info$Labels,by=1)
  groups = seq(1,parameters$Dataset.Info$Labels,by=1)
  all.partitions = data.frame(partitions, groups)

  partitions = seq(2,(parameters$Dataset.Info$Labels-1),by=1)
  groups = seq(2,(parameters$Dataset.Info$Labels-1),by=1)
  all.hybrid.partitions = data.frame(partitions, groups)

  # saving information
  setwd(parameters$Folders$folderReports)
  write.csv(all.partitions, "all-partitions.csv", row.names = FALSE)

  # data frames
  todos.tr.h = data.frame()
  todos.eb = data.frame()
  todos.fg = data.frame()
  todos.wt = data.frame()
  todos = data.frame()

  # vectors
  total = c(0)
  nomes = c()

  # for all folds
  f = 1
  while(f<=parameters$Number.Folds){

    #  "/dev/shm/j-GpositiveGO/Partitions/Split-1"
    FolderSplit = paste(parameters$Folders$folderPartition,
                        "/Split-", f,
                        sep="")

    # "/dev/shm/ro-ma-knh-GpositiveGO/Communities/Split-1"
    FolderSC = paste(parameters$Folders$folderCommunities,
                     "/Split-", f,
                     sep="")

    # file name
    # "/dev/shm/j-GpositiveGO/Partitions/Split-1/fold-1-Tr-h-choosed.csv"
    tr.h = paste(FolderSplit, "/fold-", f, "-tr-h-choosed.csv", sep="")

    # open file with all methods chosen for FOLD=1
    a.tr.h = data.frame(read.csv(tr.h))

    num.spar = nrow(a.tr.h)

    # how many methods there are?
    total[f] = nrow(a.tr.h)

    # what is it folds?
    nomes[f] = paste("fold-",f,sep="")

    # save the data frame with that information
    todos.tr.h = rbind(todos.tr.h, a.tr.h)

    num.tr = length(parameters$Valid.TR)

    # through all Tr sparcification
    k = 0
    while(k<num.tr){

      cat("\n#===================================================#")
      cat("\n# FOLD \t", f, "                                      #")
      cat("\n# TR \t", k, "                                        #")
      cat("\n#===================================================#\n")

      # "/dev/shm/j-GpositiveGO/Partitions/Split-1/Tr-1"
      FolderTr = paste(FolderSplit, "/Tr-", k, sep="")

      #  "/dev/shm/ro-ma-knh-GpositiveGO/Communities/Split-1/Tr-1"
      FolderTC = paste(FolderSC, "/Tr-", k, sep="")

      # file names
      eb = paste(FolderTr, "/tr-", k,
                 "-eb-partitions-hierarchical.csv",
                 sep="")

      eb.c = paste(FolderTC, "/tr-", k,
                   "-eb-partitions-hierarchical.csv",
                   sep="")

      fg = paste(FolderTr, "/tr-", k,
                 "-fg-partitions-hierarchical.csv",
                 sep="")

      fg.c = paste(FolderTC, "/tr-", k,
                   "-fg-partitions-hierarchical.csv",
                   sep="")

      wt = paste(FolderTr, "/tr-", k,
                 "-wt-partitions-hierarchical.csv",
                 sep="")

      wt.c = paste(FolderTC, "/tr-", k,
                   "-wt-partitions-hierarchical.csv",
                   sep="")

      if(file.exists(eb)==TRUE){

        # open files
        a.eb = data.frame(read.csv(eb))
        a.fg = data.frame(read.csv(fg))
        a.wt = data.frame(read.csv(wt))

        # gather info with the specific fold and Tr
        d.eb = data.frame(fold = f, Tr = k, a.eb)
        d.fg = data.frame(fold = f, Tr = k, a.fg)
        d.wt = data.frame(fold = f, Tr = k, a.wt)

        todos.eb = rbind(todos.eb, d.eb)
        todos.fg = rbind(todos.fg, d.fg)
        todos.wt = rbind(todos.wt, d.wt)

        m.eb = data.frame(fold = f, Tr = k, method = "eb", a.eb)
        m.fg = data.frame(fold = f, Tr = k, method = "fg", a.fg)
        m.wt = data.frame(fold = f, Tr = k, method = "wt", a.wt)

        res = rbind(m.eb, m.fg, m.wt)
        todos = rbind(todos, res)

      } else {

        # open files
        a.eb = data.frame(read.csv(eb.c))
        a.fg = data.frame(read.csv(fg.c))
        a.wt = data.frame(read.csv(wt.c))

        # gather info with the specific fold and Tr
        d.eb = data.frame(fold = f, Tr = k, a.eb)
        d.fg = data.frame(fold = f, Tr = k, a.fg)
        d.wt = data.frame(fold = f, Tr = k, a.wt)

        todos.eb = rbind(todos.eb, d.eb)
        todos.fg = rbind(todos.fg, d.fg)
        todos.wt = rbind(todos.wt, d.wt)

        m.eb = data.frame(fold = f, Tr = k, method = "eb", a.eb)
        m.fg = data.frame(fold = f, Tr = k, method = "fg", a.fg)
        m.wt = data.frame(fold = f, Tr = k, method = "wt", a.wt)

        res = rbind(m.eb, m.fg, m.wt)
        todos = rbind(todos, res)

      }

      k = k + 1
      gc()
    } # fim do Tr

    f = f + 1
    gc()
  } # fim do fold

  setwd(parameters$Folders$folderReports)

  # saving sparcification
  sparcification = data.frame(nomes, total)
  n = mean(sparcification$total)
  sparcification.2 = data.frame(sparcification,
                                minimum = min(sparcification$total),
                                maximum = max(sparcification$total),
                                mean = trunc(n,1))
  names(sparcification.2)[c(1,2)]=c("folds", "tr")
  write.csv(sparcification.2, "sparcification.csv", row.names = FALSE)

  # saving all information
  write.csv(todos, "all-partitions.choosed.csv", row.names = FALSE)
  write.csv(todos.tr.h, "all-Tr-h-choosed.csv", row.names = FALSE)
  write.csv(todos.eb, "all-eb-partitions.csv", row.names = FALSE)
  write.csv(todos.fg, "all-fg-partitions.csv", row.names = FALSE)
  write.csv(todos.wt, "all-wt-partitions.csv", row.names = FALSE)

  # return
  retorno$all.partitions.choosed = todos
  retorno$all.methods.choosed = todos.tr.h
  retorno$all.eb.partitions = todos.eb
  retorno$all.fg.partitions = todos.fg
  retorno$all.wt.partitions = todos.wt
  retorno$sparcification = sparcification.2
  retorno$all.partitions = all.partitions
  retorno$all.hybrid.partitions = all.hybrid.partitions
  return(retorno)

  cat("\n\n##########################################################")
  cat("\n# FINISH CHOOSED                                           #")
  cat("\n############################################################\n\n")

}






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




###############################################################################
#
###############################################################################
verifying.tresholds <- function(parameters){

  retorno = list()

  todos = data.frame()
  total = c()

  f = 1
  while(f<=10){
    cat("\nF =",f)
    FolderPSplit = paste(parameters$Folders$folderPartitions,
                         "/Split-", f, sep="")
    setwd(FolderPSplit)
    escolhidos = data.frame(read.csv(paste("fold-",f,
                                           "-tr-h-choosed.csv",
                                           sep="")))
    total[f] = nrow(escolhidos)
    todos = rbind(todos, escolhidos)
    f = f + 1
    gc()
  }

  maximo = max(total)
  nomes = c()
  x = 1
  while(x<=maximo){
    nomes[x] = paste("tr-",x-1, sep="")
    x = x + 1
    gc()
  }

  res = list()
  totalFolds = c()
  escolhidoFinal = c()
  x = 1
  while(x<=maximo){
    res[[x]] = filter(todos, sparsification==nomes[x])
    a = nrow(res[[x]])
    totalFolds[x] = a

    if(a==0){
      #cat("\nNÃO TEM NADA")
    } else {
      setwd(parameters$Folders$folderPartitions)
      write.csv(res[[x]], paste(nomes[x],".csv", sep=""),
                row.names = FALSE)

      setwd(parameters$Folders$folderReports)
      write.csv(res[[x]], paste(nomes[x],".csv", sep=""),
                row.names = FALSE)
    }

    if(a==10){
      res2 = filter(res[[x]], method=="none")
      if(nrow(res2)==0){
        escolhidoFinal[x] = nomes[x]
      } else {
        #cat("\nnão dá pra usar!")
      }
    }

    x = x + 1
    gc()
  }

  setwd(parameters$Folders$folderPartitions)
  write.csv(data.frame(escolhidoFinal), "escolhidos.csv")

  setwd(parameters$Folders$folderReports)
  write.csv(data.frame(escolhidoFinal), "escolhidos.csv")

  cat("\n##############################################")
  cat("\n# END: verifying.tresholds             #")
  cat("\n##############################################")

  cat("\n")
  gc()
  cat("\n")

  retorno$tr_valid = escolhidoFinal
  return(retorno)

}
