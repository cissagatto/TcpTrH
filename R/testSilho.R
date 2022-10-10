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



###############################################################################
#
###############################################################################
silho.build.data.frame <- function(){

  data <- data.frame(matrix(NA,    # Create empty data frame
                            nrow = 22,
                            ncol = 11))

  measures = c("accuracy","average-precision","clp","coverage","F1","hamming-loss","macro-AUC",
               "macro-F1","macro-precision","macro-recall","margin-loss","micro-AUC","micro-F1",
               "micro-precision","micro-recall","mlp","one-error","precision","ranking-loss",
               "recall","subset-accuracy","wlp")

  data$X1 = measures

  return(data)

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

      setwd(parameters$Folders$folderValSilho)
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

  setwd(parameters$Folders$folderValSilho)
  write.csv(data.frame(escolhidoFinal), "escolhidos.csv")

  cat("\n##############################################")
  cat("\n# END: silho.verifying.tresholds             #")
  cat("\n##############################################")

  cat("\n")
  gc()
  cat("\n")

  retorno$tr_valid = escolhidoFinal
  return(retorno)

}



###############################################################################
#
###############################################################################
silho.build.test <- function(parameters){

  f = 1
  buildParalel <- foreach(f = 1:parameters$Number.Folds ) %dopar%{
  #while(f<=parameters$Number.Folds){

    cat("\n#=========================================================")
    cat("\n#Fold: ", f)
    cat("\n#=========================================================")

    #####################################################################
    FolderRoot = "~/TCP-TR-H-Clus"
    FolderScripts = paste(FolderRoot, "/R", sep="")

    setwd(FolderScripts)
    source("utils.R")

    setwd(FolderScripts)
    source("libraries.R")


    #######################################################################
    #cat("\nLOAD FUNCTION CONVERT ARFF \n")
    converteArff <- function(arg1, arg2, arg3){
      str = paste("java -jar ", parameters$Folders$folderUtils,
                  "/R_csv_2_arff.jar ", arg1, " ", arg2, " ", arg3, sep="")
      print(system(str))
      cat("\n")
    }

    ##########################################################################
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

    #####################################################################
    FolderPSplit = paste(parameters$Folders$folderPartitions,
                         "/Split-", f, sep="")

    FolderSplitComm = paste(parameters$Folders$folderCommunities,
                            "/Split-", f, sep="")

    FolderTSplit = paste(parameters$Folders$folderTestSilho,
                         "/Split-", f, sep="")
    if(dir.create(FolderTSplit)==FALSE){dir.create(FolderTSplit)}


    ########################################################################################
    cat("\nOpen Train file ", f, "\n")
    setwd(parameters$Folders$folderCVTR)
    nome_arq_tr = paste(parameters$Dataset.Name, "-Split-Tr-", f, ".csv", sep="")
    cat("\n\t\t", nome_arq_tr)
    arquivo_tr = data.frame(read.csv(nome_arq_tr))

    ########################################################################################
    cat("\nOpen Validation file ", f, "\n")
    setwd(parameters$Folders$folderCVVL)
    nome_arq_vl = paste(parameters$Dataset.Name, "-Split-Vl-", f, ".csv", sep="")
    cat("\n\t\t", nome_arq_vl)
    arquivo_vl = data.frame(read.csv(nome_arq_vl))

    ########################################################################################
    cat("\nOpen Test file ", f, "\n")
    setwd(parameters$Folders$folderCVTS)
    nome_arq_ts = paste(parameters$Dataset.Name, "-Split-Ts-", f, ".csv", sep="")
    cat("\n\t\t", nome_arq_ts)
    arquivo_ts = data.frame(read.csv(nome_arq_ts))

    # juntando treino com validação
    arquivo_tr2 = rbind(arquivo_tr, arquivo_vl)

    u = 0
    while(u<length(parameters$valid_tr)){

        cat("\n#=========================================================")
        cat("\n# tr = ", u)
        cat("\n#=========================================================")

        FolderPartTR = paste(FolderPSplit, "/Tr-", u, sep="")
        FolderPartComm = paste(FolderSplitComm, "/Tr-", u, sep="")

        FolderTest = paste(FolderTSplit, "/Tr-", u, sep="")
        if(dir.create(FolderTest)==FALSE){dir.create(FolderTest)}

        setwd(parameters$Folders$folderReports)
        all.best = data.frame(read.csv("all-best-silho.csv"))
        all.best = filter(all.best, fold == f)
        all.best = filter(all.best, tr == u)
        num.part = all.best$part
        num.groups = all.best$part

        setwd(parameters$Folders$folderReports)
        all.choosed = data.frame(read.csv("all-choosed.csv"))
        all.choosed = data.frame(filter(all.choosed, split == f))
        escolhido = parameters$valid_tr[u+1]
        all.choosed = filter(all.choosed, all.choosed$sparsification == escolhido)
        method = all.choosed$method

        setwd(parameters$Folders$folderReports)
        nome = paste("all-", method, "-partitions.csv",sep="")
        all.partitions = data.frame(read.csv(nome))
        all.partitions = data.frame(filter(all.partitions, fold == f))
        all.partitions = data.frame(filter(all.partitions, tr == u))
        all.partitions = all.partitions[, c(-1,-2)]

        labels = all.partitions$labels
        groups = all.partitions[,num.groups]
        partition = data.frame(labels, groups)

        g = 1
        while(g<=num.groups){

          cat("\n#=========================================================")
          cat("\n#Group = ", g)
          cat("\n#=========================================================")

          FolderTestGroup = paste(FolderTest, "/Group-", g, sep="")
          if(dir.exists(FolderTestGroup)== FALSE){dir.create(FolderTestGroup) }

          ####################################################################
          #cat("\nSpecific Group: ", g, "\n")
          grupoEspecifico = filter(partition, groups == g)

          ####################################################################
          cat("\nTRAIN: Mount Group ", g, "\n")
          atributos_tr = arquivo_tr2[parameters$Dataset.Info$AttStart:parameters$Dataset.Info$AttEnd]
          n_a = ncol(atributos_tr)
          classes_tr = select(arquivo_tr2, grupoEspecifico$labels)
          n_c = ncol(classes_tr)
          grupo_tr = cbind(atributos_tr, classes_tr)
          fim_tr = ncol(grupo_tr)

          ####################################################################
          cat("\n\tTRAIN: Save Group", g, "\n")
          setwd(FolderTestGroup)
          nome_tr = paste(parameters$Dataset.Name, "-split-tr-", f, "-group-", g,
                          ".csv", sep="")
          write.csv(grupo_tr, nome_tr, row.names = FALSE)

          #################################################################
          cat("\n\tINICIO FIM TARGETS: ", g, "\n")
          inicio = parameters$Dataset.Info$LabelStart
          fim = fim_tr
          ifr = data.frame(inicio, fim)
          write.csv(ifr, "inicioFimRotulos.csv", row.names = FALSE)

          #################################################################
          cat("\n\tTRAIN: Convert Train CSV to ARFF ", g , "\n")
          nome_arquivo_2 = paste(parameters$Dataset.Name, "-split-tr-", f,
                                 "-group-", g, ".arff", sep="")
          arg1Tr = nome_tr
          arg2Tr = nome_arquivo_2
          arg3Tr = paste(inicio, "-", fim, sep="")
          str = paste("java -jar ", parameters$Folders$folderUtils,
                      "/R_csv_2_arff.jar ", arg1Tr, " ", arg2Tr, " ",
                      arg3Tr, sep="")
          print(system(str))

          ###################################################################
          cat("\n\tTRAIN: Verify and correct {0} and {1} ", g , "\n")
          arquivo = paste(FolderTestGroup, "/", arg2Tr, sep="")
          str0 = paste("sed -i 's/{0}/{0,1}/g;s/{1}/{0,1}/g' ", arquivo, sep="")
          print(system(str0))

          ###################################################################
          cat("\n\tTEST: Mount Group: ", g, "\n")
          atributos_ts = arquivo_ts[parameters$Dataset.Info$AttStart:parameters$Dataset.Info$AttEnd]
          classes_ts = select(arquivo_ts, grupoEspecifico$labels)
          grupo_ts = cbind(atributos_ts, classes_ts)
          fim_ts = ncol(grupo_ts)
          cat("\n\tTest Group Mounted: ", g, "\n")

          ################################################################
          cat("\n\tTEST: Save Group ", g, "\n")
          setwd(FolderTestGroup)
          nome_ts = paste(parameters$Dataset.Name, "-split-ts-", f, "-group-", g,
                          ".csv", sep="")
          write.csv(grupo_ts, nome_ts, row.names = FALSE)

          ##################################################################
          cat("\n\tTEST: Convert CSV to ARFF ", g , "\n")
          nome_arquivo_3 = paste(parameters$Dataset.Name, "-split-ts-", f,"-group-",
                                 g , ".arff", sep="")
          arg1Ts = nome_ts
          arg2Ts = nome_arquivo_3
          arg3Ts = paste(inicio, "-", fim, sep="")
          str = paste("java -jar ", parameters$Folders$folderUtils,
                      "/R_csv_2_arff.jar ", arg1Ts, " ", arg2Ts, " ", arg3Ts, sep="")
          print(system(str))

          ##################################################################
          cat("\n\tTEST: Verify and correct {0} and {1} ", g , "\n")
          arquivo = paste(FolderTestGroup, "/", arg2Ts, sep="")
          str0 = paste("sed -i 's/{0}/{0,1}/g;s/{1}/{0,1}/g' ", arquivo, sep="")
          cat("\n")
          print(system(str0))
          cat("\n")

          #####################################################################
          cat("\nCreating .s file for clus")
          if(inicio == fim){

            nome_config = paste(parameters$Dataset.Name, "-split-", f, "-group-", g,
                                ".s", sep="")
            sink(nome_config, type = "output")

            cat("[General]")
            cat("\nCompatibility = MLJ08")

            cat("\n\n[Data]")
            cat(paste("\nFile = ", nome_arquivo_2, sep=""))
            cat(paste("\nTestSet = ", nome_arquivo_3, sep=""))

            cat("\n\n[Attributes]")
            cat("\nReduceMemoryNominalAttrs = yes")

            cat("\n\n[Attributes]")
            cat(paste("\nTarget = ", fim, sep=""))
            cat("\nWeights = 1")

            cat("\n")
            cat("\n[Tree]")
            cat("\nHeuristic = VarianceReduction")
            cat("\nFTest = [0.001,0.005,0.01,0.05,0.1,0.125]")

            cat("\n\n[Model]")
            cat("\nMinimalWeight = 5.0")

            cat("\n\n[Output]")
            cat("\nWritePredictions = {Test}")
            cat("\n")
            sink()

            ################################################################

            cat("\nExecute CLUS: ", g , "\n")
            nome_config2 = paste(FolderTestGroup, "/", nome_config, sep="")
            str = paste("java -jar ", parameters$Folders$folderUtils,
                        "/Clus.jar ", nome_config2, sep="")
            print(system(str))

          } else {

            nome_config = paste(parameters$Dataset.Name, "-split-", f, "-group-", g,
                                ".s", sep="")
            sink(nome_config, type = "output")

            cat("[General]")
            cat("\nCompatibility = MLJ08")

            cat("\n\n[Data]")
            cat(paste("\nFile = ", nome_arquivo_2, sep=""))
            cat(paste("\nTestSet = ", nome_arquivo_3, sep=""))

            cat("\n\n[Attributes]")
            cat("\nReduceMemoryNominalAttrs = yes")

            cat("\n\n[Attributes]")
            cat(paste("\nTarget = ", inicio, "-", fim, sep=""))
            cat("\nWeights = 1")

            cat("\n")
            cat("\n[Tree]")
            cat("\nHeuristic = VarianceReduction")
            cat("\nFTest = [0.001,0.005,0.01,0.05,0.1,0.125]")

            cat("\n\n[Model]")
            cat("\nMinimalWeight = 5.0")

            cat("\n\n[Output]")
            cat("\nWritePredictions = {Test}")
            cat("\n")
            sink()

            cat("\nExecute CLUS: ", g , "\n")
            nome_config2 = paste(FolderTestGroup, "/", nome_config, sep="")
            str = paste("java -jar ", parameters$Folders$folderUtils,
                        "/Clus.jar ", nome_config2, sep="")
            print(system(str))

          }

          ################################################################
          cat("\n\nOpen predictions")
          nomeDoArquivo = paste(FolderTestGroup, "/", parameters$Dataset.Name,
                                "-split-", f,"-group-", g, ".test.pred.arff", sep="")
          predicoes = data.frame(foreign::read.arff(nomeDoArquivo))

          ################################################################
          cat("\nS\nPLIT PREDICTIS")
          if(inicio == fim){
            #cat("\n\nOnly one label in this group")

            ################################################################
            cat("\n\nSave Y_true")
            setwd(FolderTestGroup)
            classes = data.frame(predicoes[,1])
            names(classes) = colnames(predicoes)[1]
            write.csv(classes, "y_true.csv", row.names = FALSE)

            ################################################################
            cat("\n\nSave Y_true")
            rot = paste("Pruned.p.", colnames(predicoes)[1], sep="")
            pred = data.frame(predicoes[,rot])
            names(pred) = colnames(predicoes)[1]
            setwd(FolderTestGroup)
            write.csv(pred, "y_predict.csv", row.names = FALSE)

            ################################################################
            rotulos = c(colnames(classes))
            n_r = length(rotulos)
            gc()

          } else {

            ##################################################################
            cat("\n\nMore than one label in this group")
            comeco = 1+(fim - inicio)

            ##################################################################
            cat("\n\nSave Y_true")
            classes = data.frame(predicoes[,1:comeco])
            setwd(FolderTestGroup)
            write.csv(classes, "y_true.csv", row.names = FALSE)


            ##################################################################
            cat("\n\nSave Y_true")
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
            setwd(FolderTestGroup)
            write.csv(pred, "y_predict.csv", row.names = FALSE)
            gc()
          } # FIM DO ELSE

          # deleting files
          um = paste(parameters$Dataset.Name, "-split-", f, "-group-", g, ".model", sep="")
          dois = paste(parameters$Dataset.Name, "-split-", f, "-group-", g, ".s", sep="")
          tres = paste(parameters$Dataset.Name, "-split-tr-", f, "-group-", g, ".arff", sep="")
          quatro = paste(parameters$Dataset.Name, "-split-ts-", f, "-group-", g, ".arff", sep="")
          cinco = paste(parameters$Dataset.Name, "-split-tr-", f, "-group-", g, ".csv", sep="")
          seis = paste(parameters$Dataset.Name, "-split-ts-", f, "-group-", g, ".csv", sep="")
          sete = paste(parameters$Dataset.Name, "-split-", f, "-group-", g, ".out", sep="")
          oito = "Variance_RHE_1.csv"

          setwd(FolderTestGroup)
          unlink(um, recursive = TRUE)
          unlink(dois, recursive = TRUE)
          unlink(tres, recursive = TRUE)
          unlink(quatro, recursive = TRUE)
          unlink(cinco, recursive = TRUE)
          unlink(seis, recursive = TRUE)
          unlink(sete, recursive = TRUE)
          unlink(oito, recursive = TRUE)

          g = g + 1
          gc()
        } # fim do grupo

        u = u + 1
        gc()
      } # end tr

    # f = f + 1
    gc()
  } # fim do for each

  cat("\n##############################################")
  cat("\n# END: silho.build.test                      #")
  cat("\n##############################################")

  cat("\n")
  gc()
  cat("\n")
}


###############################################################################
#
###############################################################################
silho.gather.predicts <- function(parameters){

  #cat("\nFrom 1 to 10 folds!")
  # start build partitions
  # do fold 1 até o último fold
  f = 1
  gatherR <- foreach(f = 1:parameters$Number.Folds) %dopar%{
  # while(f<=parameters$Number.Folds){

    ############################################################################################################
    FolderRoot = "~/TCP-TR-H-Clus"
    FolderScripts = paste(FolderRoot, "/R", sep="")

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

    ########################################################################################
    #cat("\nSelect Partition for", f, "\n")
    FolderPartSplit = paste(parameters$Folders$folderPartitions,
                            "/Split-", f, sep="")

    # "/dev/shm/j-ma-trh-GpositiveGO/Test-Silho/Split-1"
    FolderTSplit = paste(parameters$Folders$folderTestSilho,
                         "/Split-", f, sep="")

    # tr
    k = 0
    while(k<length(parameters$valid_tr)){

      #  "/dev/shm/j-ma-trh-GpositiveGO/Partitions/Split-1/Tr-0"
      FolderPartTR = paste(FolderPartSplit, "/Tr-", k, sep="")

      # "/dev/shm/j-ma-trh-GpositiveGO/Test-Silho/Split-1/Tr-0"
      FolderTP = paste(FolderTSplit, "/Tr-", k, sep="")

      ########################################################################################
      apagar = c(0)
      y_true = data.frame(apagar)
      y_pred = data.frame(apagar)

      ######################################################################################################################
      bests = data.frame(parameters$bests$all.silho)
      bests = data.frame(filter(bests, bests$fold == f))
      bests = data.frame(filter(bests, bests$tr == k))

      methods = data.frame(parameters$bests$all.choosed.methods)
      methods = data.frame(filter(methods, methods$split==f))
      nome = paste("tr-", k, sep="")
      methods = data.frame(filter(methods, methods$sparsification==nome))

      partitions = data.frame(parameters$bests$all.choosed.partitions)
      partitions = data.frame(filter(partitions, partitions$fold==f))
      partitions = data.frame(filter(partitions, partitions$tr==k))
      partitions = data.frame(filter(partitions,
                                     partitions$method == methods$method))

      partitions = partitions[,c(-1:-3)]
      labels = partitions$labels
      groups = partitions[,bests$part]
      partition = data.frame(labels, groups)

      # GROUP
      g = 1
      while(g<=bests$part){

        cat("\n\n#==============================")
        cat("\n# Threshold \t", k)
        cat("\n# Fold \t\t", f)
        cat("\n# Groups \t", g)
        cat("\n#==============================")

          FolderTGroup = paste(FolderTP, "/Group-", g, sep="")

          #cat("\n\nGather y_true ", g)
          setwd(FolderTGroup)
          y_true_gr = data.frame(read.csv("y_true.csv"))
          y_true = cbind(y_true, y_true_gr)

          setwd(FolderTGroup)
          #cat("\n\nGather y_predict ", g)
          y_pred_gr = data.frame(read.csv("y_predict.csv"))
          y_pred = cbind(y_pred, y_pred_gr)

          cat("\n\nDeleting files")
          unlink("y_true.csv", recursive = TRUE)
          unlink("y_predict.csv", recursive = TRUE)
          unlink("inicioFimRotulos.csv", recursive = TRUE)
          unlink(paste(parameters$Dataset.Name, "-split-",f, "-group-",g,".test.pred.arff", sep=""),recursive=TRUE)

          g = g + 1
          gc()
        } # FIM DO GRUPO

      #cat("\n\nSave files ", g, "\n")
      setwd(FolderTP)
      y_pred = y_pred[,-1]
      y_true = y_true[,-1]
      write.csv(y_pred, "y_predict.csv", row.names = FALSE)
      write.csv(y_true, "y_true.csv", row.names = FALSE)

      k = k + 1
      gc()
    } # FIM DO tr

    # f = f + 1
    gc()
  } # fim do foreach


  cat("\n##############################################")
  cat("\n# END: silho.gather.predicts                 #")
  cat("\n##############################################")

  cat("\n")
  gc()
  cat("\n")


} # fim da função


###############################################################################
#
###############################################################################
silho.evaluate.test <- function(parameters){

  f = 1
  avalParal <- foreach(f = 1:parameters$Number.Folds) %dopar%{
  #while(f<=number_folds){

    cat("\n#=========================================================")
    cat("\n#Fold: ", f)
    cat("\n#=========================================================")

    folders = list()

    FolderRoot = "~/TCP-TR-H/"
    FolderScripts = paste(FolderRoot, "/R/", sep="")

    ########################################################################################
    #cat("\nSelect Partition for", f, "\n")
    FolderPartSplit = paste(parameters$Folders$folderPartitions, "/Split-", f, sep="")

    # "/dev/shm/j-ma-trh-GpositiveGO/Test-Silho/Split-1"
    FolderTSplit = paste(parameters$Folders$folderTestSilho, "/Split-", f, sep="")

    # "/dev/shm/j-ma-trh-GpositiveGO/Communities/Split-1"
    FolderSplitComm = paste(parameters$Folders$folderCommunities, "/Split-", f, sep="")

    # tr
    k = 0
    while(k<length(parameters$valid_tr)){

      cat("\n#=========================================================")
      cat("\n#Threshold = ", k)
      cat("\n#=========================================================")

      # "/dev/shm/j-ma-trh-GpositiveGO/Partitions/Split-1/Tr-0"
      FolderPartTr = paste(FolderPartSplit, "/Tr-", k, sep="")

      # "/dev/shm/j-ma-trh-GpositiveGO/Test-Silho/Split-1/Tr-0"
      FolderTP = paste(FolderTSplit, "/Tr-", k, sep="")

      # "/dev/shm/j-ma-trh-GpositiveGO/Communities/Split-1/Tr-0"
      FolderPartComm = paste(FolderSplitComm, "/Tr-", k, sep="")

      #cat("\nData frame")
      apagar = c(0)
      confMatPartitions = data.frame(apagar)
      partitions = c()

      #cat("\nGet the true and predict lables")
      setwd(FolderTP)
      y_true = data.frame(read.csv("y_true.csv"))
      y_pred = data.frame(read.csv("y_predict.csv"))

      #cat("\nCompute measures multilabel")
      y_true2 = data.frame(sapply(y_true, function(x) as.numeric(as.character(x))))
      y_true3 = mldr_from_dataframe(y_true2 , labelIndices = seq(1,ncol(y_true2 )), name = "y_true2")
      y_pred2 = sapply(y_pred, function(x) as.numeric(as.character(x)))

      #cat("\nSave Confusion Matrix")
      setwd(FolderTP)
      salva3 = paste("Conf-Mat-Fold-", f, "-tr-", k, ".txt", sep="")
      sink(file=salva3, type="output")
      confmat = multilabel_confusion_matrix(y_true3, y_pred2)
      print(confmat)
      sink()

      #cat("\nCreating a data frame")
      confMatPart = multilabel_evaluate(confmat)
      confMatPart = data.frame(confMatPart)
      names(confMatPart) = paste("Fold-", f, "-tr-", k, sep="")
      namae = paste("Split-", f, "-tr-", k,"-Evaluated.csv", sep="")
      setwd(FolderTP)
      write.csv(confMatPart, namae)

      #cat("\nDelete files")
      setwd(FolderTP)
      unlink("y_true.csv", recursive = TRUE)
      unlink("y_predict.csv", recursive = TRUE)

      k = k + 1
      gc()
    } # FIM DO tr

    #f = f + 1
    gc()
  } # fim do for each


  cat("\n##############################################")
  cat("\n# END: silho.evaluate.test                   #")
  cat("\n##############################################")

  cat("\n")
  gc()
  cat("\n")

}



###############################################################################
#
###############################################################################
silho.gather.evaluated <- function(parameters){

  # vector with names
  measures = c("accuracy","average-precision","clp","coverage","F1","hamming-loss","macro-AUC",
               "macro-F1","macro-precision","macro-recall","margin-loss","micro-AUC","micro-F1",
               "micro-precision","micro-recall","mlp","one-error","precision","ranking-loss",
               "recall","subset-accuracy","wlp")

  TRS = data.frame()

  # from fold = 1 to number_folders
  f = 1
  while(f<=parameters$Number.Folds){

    # data frame
    apagar = c(0)
    avaliadoFinal = data.frame(apagar, measures)
    avaliadoTr = data.frame(apagar, measures)
    folds = c(0)
    threshold = c(0)
    nomesThreshold = c(0)
    nomesFolds = c(0)

    ########################################################################################
    FolderPartSplit = paste(parameters$Folders$folderPartitions, "/Split-", f, sep="")
    FolderTSplit = paste(parameters$Folders$folderPartitions, "/Split-", f, sep="")
    FolderSplit = paste(parameters$Folders$folderTestSilho, "/Split-", f, sep="")

    numberTRs = 0
    # tr
    k = 0
    while(k<length(parameters$valid_tr)){

      cat("\n#==========================================")
      cat("\n#Threshold ", k)
      cat("\n#Fold ", f)
      cat("\n#===========================================")

      FolderTestT = paste(FolderSplit, "/Tr-", k, sep="")
      FolderPart = paste(FolderTSplit, "/Tr-", k, sep="")

      ######################################################################################################################
      setwd(FolderTestT)
      #setwd(Foldertr)
      str = paste("Split-", f, "-tr-", k, "-Evaluated.csv", sep="")
      avaliado = data.frame(read.csv(str))
      names(avaliado)[1] = "medidas"
      avaliadoTr = cbind(avaliadoTr, avaliado[,2])
      a = k + 1
      nomesThreshold[a] = paste("Fold-", f, "-tr-", k, sep="")
      names(avaliadoTr)[a+2] = nomesThreshold[a]
      unlink(str, recursive = TRUE)

      numberTRs = numberTRs + 1

      k = k + 1
      gc()
    } # FIM DO tr

    fold = f
    TRS = rbind(TRS, data.frame(fold, numberTRs))

    avaliadoTr = avaliadoTr[,-1]
    setwd(FolderSplit)
    write.csv(avaliadoTr, paste("Evaluated-Fold-", f, ".csv", sep=""), row.names = FALSE)

    f = f + 1
    gc()

  } # end folds

  setwd(parameters$Folders$folderTestSilho)
  write.csv(TRS, "Number-TRs.csv", row.names = FALSE)


  cat("\n##############################################")
  cat("\n# END: silho.gather.evaluated                #")
  cat("\n##############################################")

  cat("\n")
  gc()
  cat("\n")

}


###############################################################################
#
###############################################################################
silho.organize.evaluation <- function(parameters){

  dfs = list()
  dfs2 = list()

  x = 1
  while(x<=parameters$Number.Folds){
    dfs[[x]] = silho.build.data.frame()
    x = x + 1
    gc()
  }

   # from fold = 1 to number_folders
  f = 1
  while(f<=parameters$Number.Folds){
    cat("\nFold: ", f)

    ########################################################################################
    FolderPartSplit = paste(parameters$Folders$folderPartitions, "/Split-", f, sep="")
    setwd(FolderPartSplit)
    tr_H = data.frame(read.csv(paste("fold-", f, "-tr-h-choosed.csv", sep="")))
    total_tr_H = nrow(tr_H)

    FolderTest = paste(parameters$Folders$folderTestSilho, "/Split-", f, sep="")

    ######################################################################################################################
    #setwd(FolderTemptr)
    setwd(FolderTest)
    str = paste("Evaluated-Fold-", f, ".csv", sep="")
    dfs2[[f]] = data.frame(read.csv(str))
    numcol = ncol(dfs2[[f]])-1

    unlink(str, recursive = TRUE)

    f = f + 1
    gc()

  } # end folds

  numCol = ncol(dfs2[[1]])-1

  # vector with names
  measures = c("accuracy","average-precision","clp","coverage","F1","hamming-loss","macro-AUC",
               "macro-F1","macro-precision","macro-recall","margin-loss","micro-AUC","micro-F1",
               "micro-precision","micro-recall","mlp","one-error","precision","ranking-loss",
               "recall","subset-accuracy","wlp")
  apagar = c(0)
  nomestr = c()
  nomes = c()

  setwd(parameters$Folders$folderTestSilho)
  trs = data.frame(read.csv("Number-TRs.csv"))
  minimo = data.frame(apply(trs, 2, min))
  names(minimo) = "minimo"
  minimo = as.numeric(minimo[2,])

  k = 0
  while(k<length(parameters$valid_tr)){
    cat("\n\nK: ", k)

    resultado = data.frame(measures, apagar)
    nomesFolds = c()
    nometr1 = paste("Evaluated-10Folds-tr-", k, ".csv", sep="")
    nometr2 = paste("Mean-10Folds-tr-", k, ".csv", sep="")
    nometr3 = paste("Median-10Folds-tr-", k, ".csv", sep="")
    nometr4 = paste("SD-10Folds-tr-", k, ".csv", sep="")

    f = 1
    while(f<=number_folds){
      cat("\n\tF: ", f)

      # pegando apenas o fold especifico
      res = data.frame(dfs2[[f]])

      if(ncol(res)==2){
        #cat("\nApenas um threshold")
        resultado = cbind(resultado, res)
        nomes[f] = paste("Fold-",f,"-tr-", k, sep="")

      } else {
        #cat("\nMais de um threshold")

        res2 = res[,-1]
        nomesColunas = colnames(res2)

        # pegando a partir da segunda coluna
        a = k + 1
        res3 = res2[,a]

        resultado = cbind(resultado, res3)
        b = ncol(resultado)
        names(resultado)[b] = nomesColunas[a]

        nomes[f] = paste("Fold-",f,"-tr-", k, sep="")

      }

      f = f + 1
      gc()
    } # fim do fold


    resultado = data.frame(resultado[,-2])
    colnames(resultado) = c("measures", nomes)
    setwd(parameters$Folders$folderTestSilho)
    write.csv(resultado, nometr1, row.names = FALSE)

    # calculando a média dos 10 folds para cada medida
    media = data.frame(apply(resultado[,-1], 1, mean))
    media = cbind(measures, media)
    names(media) = c("Measures", "Mean10Folds")
    write.csv(media, nometr2, row.names = FALSE)

    mediana = data.frame(apply(resultado[,-1], 1, median))
    mediana = cbind(measures, mediana)
    names(mediana) = c("Measures", "Median10Folds")
    write.csv(mediana, nometr3, row.names = FALSE)

    dp = data.frame(apply(resultado[,-1], 1, sd))
    dp = cbind(measures, dp)
    names(dp) = c("Measures", "SD10Folds")
    write.csv(dp, nometr4, row.names = FALSE)

    k = k + 1
    gc()
  } # fim do k


  cat("\n##############################################")
  cat("\n# END: silho.organized.evaluation           #")
  cat("\n##############################################")

  cat("\n")
  gc()
  cat("\n")

}

silho.test <- function(parameters){

  cat("\n\n#########################################################")
    cat("\n# TEST: silho.build.test()                              #")
    cat("\n#########################################################\n\n")
  timeBuild = system.time(resBT <- silho.build.test(parameters))


  cat("\n\n#################################################")
    cat("\n# TEST: silho.gather.predicts()                 #")
    cat("\n#################################################\n\n")
  timeSplit = system.time(resGather <- silho.gather.predicts(parameters))


  cat("\n\n##############################################")
    cat("\n# TEST: silho.evaluate.test()                #")
    cat("\n##############################################\n\n")
  timeAvalia = system.time(resEval <- silho.evaluate.test(parameters))


  cat("\n\n####################################################")
    cat("\n# TEST: silho.gather.evaluation()                  #")
    cat("\n####################################################\n\n")
  timeGather = system.time(resGE <- silho.gather.evaluated(parameters))


  cat("\n\n##########################################################")
    cat("\n# TEST: silho.organize.evaluation()                      #")
    cat("\n##########################################################\n\n")
  timeOrg = system.time(resGE <- silho.organize.evaluation(parameters))


  cat("\n\n###################################################################")
  cat("\n# TEST: Save Runtime                                           #")
  cat("\n#####################################################################\n\n")
  Runtime = rbind(timeBuild,
                  timeSplit,
                  timeAvalia,
                  timeGather,
                  timeOrg)
  setwd(parameters$Folders$folderTestSilho)
  write.csv(Runtime, paste(parameters$Dataset.Name,
                           "-testPartition-Runtime.csv", sep=""),
            row.names = FALSE)


}


##################################################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com                                   #
# Thank you very much!                                                                           #
##################################################################################################
