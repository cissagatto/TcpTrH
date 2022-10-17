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

    num.tr = length(parameters$valid_tr)

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
maf1.best.partitions <- function(parameters){

  retorno = list()
  parameters = parameters

  if(parameters$Best == 8){

    cat("\n\n#=================================================#")
      cat("\n# MACRO F1                                        #")
      cat("\n#=================================================#\n\n")

    num.tr = length(parameters$valid_tr)

    k = 0
    while(k<num.tr){

      # "/dev/shm/j-GpositiveGO/Val-MaF1/Tr-1"
      MaF1.FolderTr = paste(parameters$Folders$folderValMaF1,
                             "/Tr-", k, sep="")

      apagar = c(0)
      MaF1.all = data.frame(apagar)
      MaF1.Folds.Avaliados.all = data.frame()

      nomes = c()
      measures = c()

      num.part = parameters$Dataset.Info$Labels-1

      p = 2
      a = 1
      while(p<=num.part){

        cat("\n#=================================================#")
        cat("\n# Tr ", k, " Partition ", p, "                   #")
        cat("\n#=================================================#\n")

        # "/dev/shm/j-GpositiveGO/Val-MaF1/Tr-1/Partition-2"
        MaF1.FolderPartition = paste(MaF1.FolderTr,
                                     "/Partition-", p, sep="")

        #  "/dev/shm/j-GpositiveGO/Val-MaF1/Partition-2/
        # Partition-2-Mean-10-folds-Validation.csv"
        MaF1.nome = paste(MaF1.FolderPartition, "/Partition-", p,
                          "-Mean-10-folds-Validation.csv", sep="")

        MaF1.avaliado = data.frame(read.csv(MaF1.nome))
        measures = MaF1.avaliado$measures

        MaF1.all = cbind(MaF1.all, MaF1.avaliado$Mean.10Folds)
        nomes[a] = paste("Partition-", p, sep="")

        # Partition-2-Evaluated-Validation.csv
        MaF1.nome.2 = paste(MaF1.FolderPartition, "/Partition-", p,
                            "-Evaluated-Validation.csv", sep="")

        MaF1.Folds.Avaliados = data.frame(read.csv(MaF1.nome.2))
        MaF1.Folds.Avaliados = MaF1.Folds.Avaliados[parameters$Best,]
        MaF1.Folds.Avaliados$measures = paste("partition-", p, sep="")
        MaF1.Folds.Avaliados.all = rbind(MaF1.Folds.Avaliados.all,
                                         MaF1.Folds.Avaliados)

        setwd(MaF1.FolderPartition)
        unlink(MaF1.nome)

        a = a + 1
        p = p + 1
        gc()
      }

      MaF1.all = MaF1.all[,-1]
      names(MaF1.all) = nomes
      MaF1.all = cbind(measures, MaF1.all)

      setwd(MaF1.FolderTr)
      write.csv(MaF1.all, "all-partitions-macrof1.csv",
                row.names = FALSE)

      # verificando qual particão é melhor em cada fold
      sparc = k
      fold = c(0)
      partition = c(0)
      value = c(0)
      result = data.frame(sparc, fold, partition, value)

      f = 1
      while(f<=parameters$Number.Folds){

        cat("\n#=================================================#")
        cat("\n# FOLD ", f, "                                    #")
        cat("\n#=================================================#\n")

        a = f + 1
        value_max = as.numeric(max(MaF1.Folds.Avaliados.all[,a]))
        index_max = which.max(MaF1.Folds.Avaliados.all[,a])

        fold = f
        partition = as.numeric(index_max+1)
        value = as.numeric(value_max)
        result.1 = data.frame(sparc, fold, partition, value)
        result = rbind(result, result.1)

        f = f + 1
        gc()
      }

      result = result[-1,]
      setwd(parameters$Folders$folderReports)
      write.csv(result,
                paste("tr-", k, "-Best-MacroF1.csv", sep=""),
                row.names = FALSE)

      frequency = data.frame(count(result, partition))
      names(frequency) = c("partition", "frequency")

      write.csv(frequency,
                paste("tr-", k, "-frequency-best-macroF1.csv", sep=""),
                row.names = FALSE)

      soma = apply(result, 2, sum)
      media = apply(result, 2, mean)
      mediana = apply(result, 2, median)
      desvioPadrao = apply(result, 2, sd)
      minimo = apply(result, 2, min)
      maximo = apply(result, 2, max)
      sumario = rbind(soma, media, mediana, desvioPadrao,
                      minimo, maximo)

      write.csv(sumario,
                paste("tr-", k, "-summary-best-macroF1.csv", sep=""),
                row.names = FALSE)

      k = k + 1
      gc()
    }

  } else {

    cat("\n\n#=================================================#")
    cat("\n# MICRO F1                                        #")
    cat("\n#=================================================#\n\n")


    all.partitions = data.frame(read.csv("all-partitions.csv"))
    ultimo = nrow(all.partitions)
    all.partitions = all.partitions[c(-1, -ultimo),]
    num.part = parameters$Dataset.Info$Labels-1

    num.tr = length(parameters$valid_tr)

    k = 1
    while(k<num.tr){

      # "/dev/shm/j-GpositiveGO/Val-MaF1/Tr-1"
      MiF1.FolderTr = paste(parameters$Folders$folderValMiF1,
                             "/Tr-", k, sep="")

      apagar = c(0)
      MiF1.all = data.frame(apagar)
      MiF1.Folds.Avaliados.all = data.frame()

      nomes = c()
      measures = c()

      p = 2
      a = 1
      while(p<=num.part){

        cat("\n#=================================================#")
        cat("\n# TR ", k, " Partition ", p, "                   #")
        cat("\n#=================================================#\n")

        # "/dev/shm/j-GpositiveGO/Val-MaF1/Tr-1/Partition-2"
        MiF1.FolderPartition = paste(MiF1.FolderTr,
                                     "/Partition-", p, sep="")

        #  "/dev/shm/j-GpositiveGO/Val-MaF1/Partition-2/
        # Partition-2-Mean-10-folds-Validation.csv"
        MiF1.nome = paste(MiF1.FolderPartition, "/Partition-", p,
                          "-Mean-10-folds-Validation.csv", sep="")

        MiF1.avaliado = data.frame(read.csv(MaF1.nome))
        measures = MiF1.avaliado$measures

        MiF1.all = cbind(MiF1.all, MiF1.avaliado$Mean.10Folds)
        nomes[a] = paste("Partition-", p, sep="")

        # Partition-2-Evaluated-Validation.csv
        MiF1.nome.2 = paste(MiF1.FolderPartition, "/Partition-", p,
                            "-Evaluated-Validation.csv", sep="")

        MiF1.Folds.Avaliados = data.frame(read.csv(MiF1.nome.2))
        MiF1.Folds.Avaliados = MiF1.Folds.Avaliados[parameters$Best,]
        MiF1.Folds.Avaliados$measures = paste("partition-", p, sep="")
        MiF1.Folds.Avaliados.all = rbind(MiF1.Folds.Avaliados.all,
                                         MiF1.Folds.Avaliados)

        setwd(MiF1.FolderPartition)
        unlink(MiF1.nome)

        a = a + 1
        p = p + 1
        gc()
      }

      MiF1.all = MiF1.all[,-1]
      names(MiF1.all) = nomes
      MiF1.all = cbind(measures, MiF1.all)

      setwd(MiF1.FolderTr)
      write.csv(MiF1.all, "all-partitions-microf1.csv",
                row.names = FALSE)

      # verificando qual particão é melhor em cada fold
      sparc = k
      fold = c(0)
      partition = c(0)
      value = c(0)
      result = data.frame(sparc, fold, partition, value)

      f = 1
      while(f<=parameters$Number.Folds){

        cat("\n#=================================================#")
        cat("\n# FOLD ", f, "                                    #")
        cat("\n#=================================================#\n")

        a = f + 1
        value_max = as.numeric(max(MiF1.Folds.Avaliados.all[,a]))
        index_max = which.max(MiF1.Folds.Avaliados.all[,a])

        fold = f
        partition = as.numeric(index_max+1)
        value = as.numeric(value_max)
        result.1 = data.frame(sparc, fold, partition, value)
        result = rbind(result, result.1)

        f = f + 1
        gc()
      }

      result = result[-1,]
      setwd(parameters$Folders$folderReports)
      write.csv(result,
                paste("Tr-", k, "-Best-MicroF1.csv", sep=""),
                row.names = FALSE)

      frequency = data.frame(count(result, partition))
      names(frequency) = c("partition", "frequency")

      write.csv(frequency,
                paste("Tr-", k, "-frequency-best-microF1.csv", sep=""),
                row.names = FALSE)

      soma = apply(result, 2, sum)
      media = apply(result, 2, mean)
      mediana = apply(result, 2, median)
      desvioPadrao = apply(result, 2, sd)
      minimo = apply(result, 2, min)
      maximo = apply(result, 2, max)
      sumario = rbind(soma, media, mediana, desvioPadrao,
                      minimo, maximo)

      write.csv(sumario,
                paste("Tr-", k, "-summary-best-microF1.csv", sep=""),
                row.names = FALSE)

      k = k + 1
      gc()
    }

  } # fim do IF/ELSE

  gc()
  cat("\n###############################################################")
  cat("\n# Best Partitions: END                                        #")
  cat("\n###############################################################")
  cat("\n\n\n\n")

}





######################################################################
#
######################################################################
maf1.validate.partitions <- function(parameters){

  parameters = parameters

  num.spar = c()
  num.fold = c()
  num.part = c()
  labels.distribution = data.frame()

  f = 1
  validateParalel <- foreach(f = 1:parameters$Number.Folds) %dopar%{
  #while(f<=parameters$Number.Folds){

    parameters = parameters

    ################################################################
    num.fold = f
    num.spar = parameters$k
    num.part = parameters$id_part

    ################################################################
    FolderRoot = "~/TCP-TR-H-Clus"
    FolderScripts = paste(FolderRoot, "/R", sep="")

    ################################################################
    setwd(FolderScripts)
    source("libraries.R")

    setwd(FolderScripts)
    source("utils.R")

    ##################################################################
    converteArff <- function(arg1, arg2, arg3){
      str = paste("java -jar ", parameters$Folders$folderUtils,
                  "/R_csv_2_arff.jar ", arg1, " ", arg2, " ", arg3, sep="")
      print(system(str))
      cat("\n")
    }

    ########################################################################
    # "/dev/shm/j-GpositiveGO/Val-MaF1/Tr-1/Partition-2/Split-1"
    FolderSplitVal = paste(parameters$FolderPartition, "/Split-", f, sep="")
    if(dir.exists(FolderSplitVal)==FALSE){dir.create(FolderSplitVal)}

    # "/dev/shm/j-GpositiveGO/Communities/Split-1"
    FolderSplitComm = paste(parameters$Folders$folderCommunities,
                            "/Split-", f, sep="")

    # "/dev/shm/j-GpositiveGO/Partitions/Split-1"
    FolderPartSplit = paste(parameters$Folders$folderPartitions,
                            "/Split-",f, sep="")

    #  "/dev/shm/j-GpositiveGO/Communities/Split-1/Tr-1"
    FolderTrComm = paste(FolderSplitComm, "/Tr-", parameters$k, sep="")

    ########################################################################
    escolhidos = filter(parameters$Choosed$all.methods.choosed, split == f)
    Tr = paste("tr-", parameters$k, sep="")
    escolhido = filter(escolhidos, sparsification == Tr)

    # /dev/shm/j-GpositiveGO/Communities/Split-1/Tr-1/
    # Tr-1-eb-partitions-hierarchical.csv
    nome = paste(FolderTrComm, "/tr-", parameters$k, "-",
                 toString(escolhido$method),
                 "-partitions-hierarchical.csv",
                 sep="")

    particoes = data.frame(read.csv(nome))
    labels = particoes$labels
    groups = particoes[,parameters$id_part]
    particao = data.frame(labels, groups)

    res = data.frame(filter(parameters$Choosed$all.hybrid.partitions,
                            partitions == parameters$id_part))

    num.labels = c()
    num.group = c()

    g = 1
    while(g<=as.numeric(res$groups)){

      ###############################################################
      nome_grupo = paste("grupo_", g, sep="")
      nome_grupo_2 = paste("Group-", g, sep="")

      cat("\n\n#=========================================================#")
        cat("\n# FOLD ", f, " GROUP ", g, "                              #")
        cat("\n#=========================================================#\n\n")

      #################################################################
      #cat("\n\nCreating folder group")
      #  "/dev/shm/j-GpositiveGO/Val-MaF1/Tr-1/Partition-2/Split-1/Group-1"
      FolderGroup = paste(FolderSplitVal, "/", nome_grupo_2, sep="")
      if(dir.exists(FolderGroup)==FALSE){dir.create(FolderGroup) }

      ################################################################
      #cat("\n\nGet the labels of this group")
      particao.especifica = particao %>% filter(., particao$groups == g)

      ################################################################
      num.labels[g] = nrow(particao.especifica)
      num.group[g] = g

      ###############################################################
      #cat("\n\nMount group")
      totalLabels = nrow(particao.especifica)

      #################################################################

      # GpositiveGO-Split-Tr-1.csv
      nomeTr = paste(parameters$Dataset.Name, "-Split-Tr-", f, ".csv", sep="")

      # GpositiveGO-Split-Vl-1.csv
      nomeVl = paste(parameters$Dataset.Name, "-Split-Vl-", f, ".csv", sep="")

      #################################################################

      # "/dev/shm/j-GpositiveGO/datasets/CrossValidation/Tr/
      # GpositiveGO-Split-Tr-1.csv"
      nomeTr2 = paste(parameters$Folders$folderCVTR, "/", nomeTr, sep="")

      # "/dev/shm/j-GpositiveGO/datasets/CrossValidation/Vl/
      # GpositiveGO-Split-Vl-1.csv"
      nomeVl2 = paste(parameters$Folders$folderCVVL, "/", nomeVl, sep="")

      ###################################################################
      #cat("\nTRAIN: MOUNT GROUP\n")
      arquivoTr = data.frame(read.csv(nomeTr2), stringsAsFactors = F)
      atributosTr = arquivoTr[parameters$Dataset.Info$AttStart:parameters$Dataset.Info$AttEnd]
      classesTr = select(arquivoTr, particao.especifica$labels)
      thisGroupTr = cbind(atributosTr, classesTr)
      ncols = ncol(thisGroupTr)

      ###############################################################
      #cat("\nTRAIN: Save CSV\n")
      nomeCsTr = paste("grupo_Tr_", g, ".csv", sep="")
      nomeArTr = paste("grupo_Tr_", g, ".arff", sep="")
      setwd(FolderGroup)
      write.csv(thisGroupTr, nomeCsTr, row.names = FALSE)

      #####################################################################
      #cat("\nTRAIN: Start End Targets\n")
      inicio = parameters$Dataset.Info$LabelStart
      fim = ncol(thisGroupTr)
      ifr = data.frame(inicio, fim)
      setwd(FolderGroup)
      write.csv(ifr, "inicioFimRotulos.csv", row.names = FALSE)

      #####################################################################
      #cat("\nTRAIN: Convert CSV to ARFF\n")
      setwd(FolderGroup)
      arg1 = nomeCsTr
      arg2 = nomeArTr
      arg3 = paste(inicio, "-", fim, sep="")
      converteArff(arg1, arg2, arg3)

      ################################################################
      #cat("\n\nTRAIN: Verify and correct {0} and {1}\n")
      arquivo = paste(FolderGroup, "/", nomeArTr, sep="")
      str0 = paste("sed -i 's/{0}/{0,1}/g;s/{1}/{0,1}/g' ", arquivo, sep="")
      print(system(str0))

      ###################################################################
      #cat("\n\nVALIDATION: MOUNT GROUP\n")
      arquivoVl = data.frame(read.csv(nomeVl2), stringsAsFactors = F)
      atributosVl = arquivoVl[parameters$Dataset.Info$AttStart:parameters$Dataset.Info$AttEnd]
      classesVl = select(arquivoVl, particao.especifica$labels)
      thisGroupVl = cbind(atributosVl, classesVl)

      #################################################################
      #cat("\n\nVALIDATION: Save CSV\n")
      nomeCsVl = paste("grupo_Vl_", g, ".csv", sep="")
      nomeArVl = paste("grupo_Vl_", g, ".arff", sep="")
      setwd(FolderGroup)
      write.csv(thisGroupVl, nomeCsVl, row.names = FALSE)

      ##################################################################
      #cat("\n\nVALIDATION: Convert CSV to ARFF\n")
      setwd(FolderGroup)
      arg1 = nomeCsVl
      arg2 = nomeArVl
      arg3 = paste(inicio, "-", fim, sep="")
      converteArff(arg1, arg2, arg3)

      #################################################################
      #cat("\n\nVALIDATION: Verify and correct {0} and {1}\n")
      arquivo = paste(FolderGroup, "/", nomeArVl, sep="")
      str0 = paste("sed -i 's/{0}/{0,1}/g;s/{1}/{0,1}/g' ", arquivo, sep="")
      system(str0)

      #################################################################
      if(totalLabels==1){

        cat("\n\n#======================================================#")
          cat("\n# SINGLE LABEL                                         #")
          cat("\n#======================================================#\n\n")

        setwd(FolderGroup)
        nome_config = paste("grupo_", g, ".s", sep="")
        sink(nome_config, type = "output")

        cat("[General]")
        cat("\nCompatibility = MLJ08")

        cat("\n")
        cat("\n[Data]")
        nome_arquivo_2 = paste("grupo_Tr_", g, ".arff", sep="")
        cat(paste("\nFile = ", nome_arquivo_2, sep=""))
        nome_arquivo_3 = paste("grupo_Vl_", g, ".arff", sep="")
        cat(paste("\nTestSet = ", nome_arquivo_3, sep=""))

        cat("\n")
        cat("\n[Attributes]")
        cat("\nReduceMemoryNominalAttrs = yes")

        cat("\n")
        cat("\n[Attributes]")
        cat(paste("\nTarget = ", inicio, sep=""))
        cat("\nWeights = 1")

        cat("\n")
        cat("\n[Tree]")
        cat("\nHeuristic = VarianceReduction")
        cat("\nFTest = [0.001,0.005,0.01,0.05,0.1,0.125]")

        cat("\n")
        cat("\n[Model]")
        cat("\nMinimalWeight = 5.0")

        cat("\n")
        cat("\n[Output]")
        cat("\nWritePredictions = {Test}")
        cat("\n")
        sink()

      } else {

        cat("\n\n#======================================================#")
          cat("\n# MULTI-LABEL                                          #")
          cat("\n#======================================================#\n\n")

        setwd(FolderGroup)
        nome_config = paste("grupo_", g, ".s", sep="")
        sink(nome_config, type = "output")

        cat("[General]")
        cat("\nCompatibility = MLJ08")

        cat("\n")
        cat("\n[Data]")
        nome_arquivo_2 = paste("grupo_Tr_", g, ".arff", sep="")
        cat(paste("\nFile = ", nome_arquivo_2, sep=""))
        nome_arquivo_3 = paste("grupo_Vl_", g, ".arff", sep="")
        cat(paste("\nTestSet = ", nome_arquivo_3, sep=""))

        cat("\n")
        cat("\n[Attributes]")
        cat("\nReduceMemoryNominalAttrs = yes")

        cat("\n")
        cat("\n[Attributes]")
        cat(paste("\nTarget = ", inicio, "-", fim, sep=""))
        cat("\nWeights = 1")

        cat("\n")
        cat("\n[Tree]")
        cat("\nHeuristic = VarianceReduction")
        cat("\nFTest = [0.001,0.005,0.01,0.05,0.1,0.125]")

        cat("\n")
        cat("\n[Model]")
        cat("\nMinimalWeight = 5.0")

        cat("\n")
        cat("\n[Output]")
        cat("\nWritePredictions = {Test}")
        cat("\n")
        sink()
      }

      ######################################################################
      #cat("\nExecute CLUS")
      nome_config2 = paste(FolderGroup, "/", nome_config, sep="")
      setwd(FolderGroup)
      str = paste("java -jar ", parameters$Folders$folderUtils,
                  "/Clus.jar ", nome_config2, sep="")
      print(system(str))

      ####################################################################
      #cat("\n\nOpen inicioFimRotulos.csv")
      targets = data.frame(read.csv("inicioFimRotulos.csv"))

      #################################################################
      #cat("\n\nOpen predictions")
      namae2 = paste(FolderGroup, "/", nome_grupo, ".test.pred.arff", sep="")
      predicoes = data.frame(foreign::read.arff(namae2), use_xml = FALSE)

      ####################################################################
      #cat("\nS\nPLIT PREDICTIS")

      if(targets$inicio == targets$fim){

        cat("\n\n#======================================================#")
          cat("\n# SINGLE LABEL                                         #")
          cat("\n#======================================================#\n\n")

        ###################################################################
        #cat("\n\nSave Y_true")
        setwd(FolderGroup)
        classes = data.frame(predicoes[,1])
        names(classes) = colnames(predicoes)[1]
        write.csv(classes, "y_true.csv", row.names = FALSE)

        ###################################################################
        #cat("\n\nSave Y_true")
        rot = paste("Pruned.p.", colnames(predicoes)[1], sep="")
        pred = data.frame(predicoes[,rot])
        names(pred) = colnames(predicoes)[1]
        setwd(FolderGroup)
        write.csv(pred, "y_predict.csv", row.names = FALSE)

        ################################################################
        rotulos = c(colnames(classes))
        n_r = length(rotulos)

        gc()

      } else {

        cat("\n\n#====================================================#")
          cat("\n# MULTI-LABEL                                        #")
          cat("\n#====================================================#\n\n")

        #################################################################
        comeco = 1+(targets$fim - targets$inicio)

        ####################################################################
        #cat("\n\nSave Y_true")
        classes = data.frame(predicoes[,1:comeco])
        setwd(FolderGroup)
        write.csv(classes, "y_true.csv", row.names = FALSE)

        ################################################################
        #cat("\n\nSave Y_true")
        rotulos = c(colnames(classes))
        n_r = length(rotulos)
        nomeColuna = c()
        t = 1
        while(t <= n_r){
          nomeColuna[t] = paste("Pruned.p.", rotulos[t], sep="")
          t = t + 1
          gc()
        }
        pred = data.frame(predicoes[nomeColuna])
        names(pred) = rotulos
        setwd(FolderGroup)
        write.csv(pred, "y_predict.csv", row.names = FALSE)

        gc()
      } # FIM DO ELSE

      # gather labels distribution in groups
      labels.res = data.frame(num.spar, num.fold,
                              num.part, num.group, num.labels)

      cat("\nDELETING UNECESSARY FOLDERS")
      nome1 = paste("grupo_Tr_", g, ".arff", sep="")
      nome2 = paste("grupo_Vl_", g, ".arff", sep="")
      nome3 = paste("grupo_Tr_", g, ".csv", sep="")
      nome4 = paste("grupo_Vl_", g, ".csv", sep="")
      nome5 = paste("grupo_", g, ".model", sep="")
      nome6 = paste("grupo_", g, ".s", sep="")
      nome7 = "Variance_RHE_1.csv"

      setwd(FolderGroup)
      unlink(nome1, recursive = TRUE)
      unlink(nome2, recursive = TRUE)
      unlink(nome3, recursive = TRUE)
      unlink(nome4, recursive = TRUE)
      unlink(nome5, recursive = TRUE)
      unlink(nome6, recursive = TRUE)
      unlink(nome7, recursive = TRUE)

      g = g + 1
      gc()
    } # fim do grupo

    # gather all labels distribution
    labels.distribution = rbind(labels.distribution, labels.res)

    # save label distribution
    setwd(parameters$FolderPartition)
    namae = paste("partition-", parameters$id_part,
                  "-distribution-labels.csv", sep="")
    write.csv(labels.distribution, namae , row.names = FALSE)

    setwd(parameters$Folders$folderReports)
    namae = paste("Tr-", parameters$k, "-partition-", parameters$id_part,
                  "-distribution-labels.csv", sep="")
    write.csv(labels.distribution, namae , row.names = FALSE)

    #f = f + 1
    gc()
  } # fim do foreach

  cat("\n\n##############################################################")
    cat("\n# BUILD AND VALIDATE PARTITION ", parameters$id_part ," END! #")
    cat("\n##############################################################\n\n")

  gc()
} # fim da função



################################################################
#
###############################################################
maf1.val.gather.predicts <- function(parameters){

  f = 1
  gatherParal <- foreach(f = 1:parameters$Number.Folds) %dopar%{
  #while(f<=parameters$Number.Folds){

    # data frame
    apagar = c(0)
    y_true = data.frame(apagar)
    y_pred = data.frame(apagar)

    cat("\n#======================================================#")
    cat("\n# FOLD ", f, "                                         #")
    cat("\n#======================================================#\n")

    ##########################################################################
    # "/dev/shm/j-GpositiveGO/Val-MaF1/Tr-1/Partition-2/Split-1"
    FolderSplit = paste(parameters$FolderPartition, "/Split-", f, sep="")

    ##########################################################################
    res = data.frame(filter(parameters$Choosed$all.hybrid.partitions, partitions == parameters$id_part))

    g = 1
    while(g<=as.numeric(res$groups)){

      cat("\n#======================================================#")
      cat("\n# GROUP ", g, "                                        #")
      cat("\n#======================================================#\n")

      FolderGroup = paste(FolderSplit, "/Group-", g, sep="")

      #cat("\n\nGather y_true ", g)
      setwd(FolderGroup)
      y_true_gr = data.frame(read.csv("y_true.csv"))
      y_true = cbind(y_true, y_true_gr)

      #cat("\n\nGather y_predict ", g)
      y_pred_gr = data.frame(read.csv("y_predict.csv"))
      y_pred = cbind(y_pred, y_pred_gr)

      # cat("\n\nDeleting files")
      unlink("y_true.csv", recursive = TRUE)
      unlink("y_predict.csv", recursive = TRUE)
      unlink("inicioFimRotulos.csv", recursive = TRUE)

      g = g + 1
      gc()
    }

    #cat("\n\nSave files ", g, "\n")
    setwd(FolderSplit)
    y_pred = y_pred[,-1]
    y_true = y_true[,-1]
    write.csv(y_pred, "y_predict.csv", row.names = FALSE)
    write.csv(y_true, "y_true.csv", row.names = FALSE)

    #f = f + 1

    gc()
  } # fim do foreach

  gc()
  cat("\n############################################################")
  cat("\n# Gather Predicts: END                                     #")
  cat("\n############################################################")
  cat("\n\n\n\n")

} # fim da função


####################################################################
#
#####################################################################
maf1.val.evaluate <- function(parameters){

  f = 1
  evalParal <- foreach(f = 1:parameters$Number.Folds) %dopar%{

    cat("\n#====================================================#")
    cat("\n# FOLD ", f, "                                       #")
    cat("\n#====================================================#\n")


    ################################################################
    FolderRoot = "~/TCP-TR-H-Clus"
    FolderScripts = paste(FolderRoot, "/R", sep="")

    ################################################################
    setwd(FolderScripts)
    source("libraries.R")

    setwd(FolderScripts)
    source("utils.R")


    # data frame
    apagar = c(0)
    confMatPartitions = data.frame(apagar)
    partitions = c()

    # specifyin folder for the fold
    FolderSplit = paste(parameters$FolderPartition, "/Split-", f, sep="")

    # get the true and predict lables
    setwd(FolderSplit)
    y_true = data.frame(read.csv("y_true.csv"))
    y_pred = data.frame(read.csv("y_predict.csv"))

    # compute measures multilabel
    y_true2 = data.frame(sapply(y_true, function(x) as.numeric(as.character(x))))
    y_true3 = mldr_from_dataframe(y_true2 ,
                                  labelIndices = seq(1,ncol(y_true2 )),
                                  name = "y_true2")
    y_pred2 = sapply(y_pred, function(x) as.numeric(as.character(x)))

    #cat("\n\t\tSave Confusion Matrix")
    setwd(FolderSplit)
    salva3 = paste("Conf-Mat-Fold-", f, ".txt", sep="")
    sink(file=salva3, type="output")
    confmat = multilabel_confusion_matrix(y_true3, y_pred2)
    print(confmat)
    sink()

    # creating a data frame
    confMatPart = multilabel_evaluate(confmat)
    confMatPart = data.frame(confMatPart)
    names(confMatPart) = paste("Fold-", f, sep="")
    namae = paste("Split-", f,"-Evaluated.csv", sep="")
    write.csv(confMatPart, namae)

    # delete files
    setwd(FolderSplit)
    unlink("y_true.csv", recursive = TRUE)
    unlink("y_predict.csv", recursive = TRUE)

    gc()
  } # end folds

  gc()
  cat("\n############################################################")
  cat("\n# Evaluation Folds: END                                    #")
  cat("\n############################################################")
  cat("\n\n\n\n")
}


####################################################################
#
###################################################################
maf1.val.gather.evaluated <- function(parameters){

  # vector with names
  measures = c("accuracy", "average-precision", "clp", "coverage", "F1",
               "hamming-loss", "macro-AUC", "macro-F1",
               "macro-precision", "macro-recall", "margin-loss",
               "micro-AUC","micro-F1", "micro-precision",
               "micro-recall", "mlp", "one-error", "precision",
               "ranking-loss", "recall", "subset-accuracy", "wlp")

  # data frame
  apagar = c(0)
  avaliado = data.frame(apagar)
  folds = c(0)
  nomesFolds = c(0)


  # from fold = 1 to number_folders
  f = 1
  while(f<=parameters$Number.Folds){

    cat("\n#======================================================#")
    cat("\n# FOLD ", f, "                                         #")
    cat("\n#======================================================#\n")

    # specifying folder for the fold
    FolderSplit = paste(parameters$FolderPartition, "/Split-", f, sep="")
    setwd(FolderSplit)
    str = paste("Split-", f, "-Evaluated.csv", sep="")
    avaliado.res = data.frame(read.csv(str))
    names(avaliado.res)[1] = "medidas"
    avaliado.res = data.frame(avaliado.res[order(avaliado.res$medidas, decreasing = FALSE),])
    avaliado.res = data.frame(avaliado.res[,-1])
    avaliado = cbind(avaliado, avaliado.res)
    #names(avaliado)[f+1] = paste("Fold-", f, sep="")
    nomesFolds[f] = paste("Fold-", f, sep="")

    setwd(FolderSplit)
    unlink(str, recursive = TRUE)

    f = f + 1
    gc()

  } # end folds

  #cat("\nSAVE MEASURES")
  avaliado$apagar = measures
  colnames(avaliado) = c("measures", nomesFolds)

  nome3 = paste(parameters$FolderPartition, "/Partition-",
                parameters$id_part, "-Evaluated-Validation.csv",
                sep="")
  write.csv(avaliado, nome3, row.names = FALSE)

  nome4 = paste(parameters$FolderPartition, "/Partition-",
                parameters$id_part,
                "-Mean-10-folds-Validation.csv", sep="")
  avaliado.2 = avaliado[,-1]
  avaliado.2 = data.frame(apply(avaliado.2, 1, mean))
  colnames(avaliado.2) = "Mean-10Folds"
  avaliado.2 = cbind(measures, avaliado.2)

  write.csv(avaliado.2, nome4, row.names = FALSE)

  gc()
  cat("\n############################################################")
  cat("\n# Evaluated Partition: END                                 #")
  cat("\n############################################################")
  cat("\n\n\n\n")
}



##########################################################################
#
##########################################################################
maf1.validate <- function(parameters){

  FolderRoot = "~/TCP-TR-H-Clus"
  FolderScripts = paste(FolderRoot, "/R", sep="")

  setwd(parameters$Folders$folderReports)
  spar = data.frame(read.csv("sparcification.csv"))

  setwd(parameters$Folders$folderNamesLabels)
  namae = paste(parameters$Dataset.Name, "-NamesLabels.csv", sep="")
  namesLabels = data.frame(read.csv(namae))
  names(namesLabels) = c("Index", "Name")

  parameters$namesLabes = namesLabels
  num.tr = length(parameters$valid_tr)

  k = 0
  while(k<num.tr){

    # "/dev/shm/j-GpositiveGO/Val-MaF1/Tr-1"
    FolderTr = paste(parameters$Folders$folderValMaF1,"/Tr-", k, sep="")
    if(dir.exists(FolderTr)==FALSE){dir.create(FolderTr)}

    # from partition 2 to last partition (n)
    count = 2
    id_part = 2
    while(id_part<parameters$Dataset.Info$Labels){

      cat("\n\n#=================================================#")
        cat("\n# TR ", parameters$k, " PARTITION ", id_part, "  #")
        cat("\n#=================================================#\n\n")

      # "/dev/shm/j-GpositiveGO/Val-MaF1/Tr-1/Partition-2"
      FolderPartition = paste(FolderTr, "/Partition-", id_part, sep="")
      if(dir.exists(FolderPartition)==FALSE){dir.create(FolderPartition)}

      parameters$k = k
      parameters$id_part = id_part
      parameters$FolderPartition = FolderPartition


      cat("\n\n#########################################################")
      cat("\n# VALIDATE: maf1.validate.partitions()                  #")
      cat("\n#########################################################\n\n")
      timeBuild = system.time(resVP <- maf1.validate.partitions(parameters))


      cat("\n\n##########################################################")
      cat("\n# VALIDATE: maf1.gather.predicts()                       #")
      cat("\n##########################################################\n\n")
      timeSplit = system.time(resGatherVal <- maf1.val.gather.predicts(parameters))


      cat("\n\n##########################################################")
      cat("\n# VALIDATE: maf1.evalute()                               #")
      cat("\n#######################################################3##\n\n")
      timeAvalia = system.time(resEVal <- maf1.val.evaluate(parameters))


      cat("\n\n##########################################################")
      cat("\n# VALIDATE: maf1.gather.evalution()                      #")
      cat("\n##########################################################\n\n")
      timeGather = system.time(resValEval <- maf1.val.gather.evaluated(parameters))


      cat("\n\n##########################################################")
      cat("\n# VALIDATE: Save Runtime                                   #")
      cat("\n############################################################\n\n")
      Runtime = rbind(timeBuild, timeSplit, timeAvalia, timeGather)
      setwd(parameters$FolderPartition)
      name2 = paste("Tr-", parameters$k,
                    "-partition-", parameters$id_part,
                    "-Runtime-Validation.csv", sep="")
      write.csv(Runtime, name2)

      id_part = id_part + 1
      count = count + 1
      gc()
    } # fim da partição

    k = k + 1
    gc()
  } # FIM DO Tr

  gc()
  cat("\n\n##########################################################")
  cat("\n# END VALIDATE                                             #")
  cat("\n############################################################\n\n")
}


#######################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com        #
# Thank you very much!                                                #
#######################################################################
