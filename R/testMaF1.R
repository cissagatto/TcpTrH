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
FolderRoot = "~/TCP-Tr-H-Clus"
FolderScripts = paste(FolderRoot, "/R", sep="")



#################################################################
#
#################################################################
maf1.build.data.frame <- function(){

  data <- data.frame(matrix(NA,    # Create empty data frame
                            nrow = 22,
                            ncol = 11))

  measures = c("accuracy", "average-precision", "clp", "coverage", "F1",
               "hamming-loss", "macro-AUC", "macro-F1", "macro-precision",
               "macro-recall", "margin-loss", "micro-AUC", "micro-F1",
               "micro-precision", "micro-recall", "mlp", "one-error",
               "precision", "ranking-loss", "recall", "subset-accuracy", "wlp")

  data$X1 = measures

  return(data)

}


##############################################################################
# FUNCTION BUILD AND TEST SELECTED HYBRID PARTITION                          #
#   Objective                                                                #
#   Parameters                                                               #
##############################################################################
maf1.build.test <- function(parameters){

  f = 1
  #buildParalel <- foreach(f = 1:parameters$Number.Folds) %dopar%{
  while(f<=parameters$Number.Folds){

    cat("\n#=========================================================")
    cat("\n# Fold: ", f)
    cat("\n#=========================================================")

    parameters = parameters

    ########################################################################
    FolderRoot = "~/TCP-TR-H-Clus"
    FolderScripts = paste(FolderRoot, "/R", sep="")

    setwd(FolderScripts)
    source("utils.R")

    setwd(FolderScripts)
    source("libraries.R")

    setwd(FolderScripts)
    source("misc.R")


    #########################################################################
    converteArff <- function(arg1, arg2, arg3){
      str = paste("java -jar ", parameters$FoldersfolderUtils,
                  "/R_csv_2_arff.jar ", arg1, " ", arg2, " ", arg3, sep="")
      print(system(str))
      cat("\n")
    }

    ###########################################################################

    # "/dev/shm/j-GpositiveGO/Partitions/Split-1"
    FolderPSplit = paste(parameters$Folders$folderPartitions,
                         "/Split-", f, sep="")

    # "/dev/shm/j-GpositiveGO/Communities/Split-1"
    FolderSplitComm = paste(parameters$Folders$folderCommunities,
                            "/Split-", f, sep="")

    # "/dev/shm/j-GpositiveGO/Test-MaF1/Split-1"
    FolderTSplit = paste(parameters$Folders$folderTestMaF1,
                         "/Split-", f, sep="")
    if(dir.create(FolderTSplit)==FALSE){dir.create(FolderTSplit)}

    ########################################################################
    #cat("\nOpen Train file ", f, "\n")
    setwd(parameters$Folders$folderCVTR)
    nome_arq_tr = paste(parameters$Dataset.Name, "-Split-Tr-", f,
                        ".csv", sep="")
    arquivo_tr = data.frame(read.csv(nome_arq_tr))

    #####################################################################
    #cat("\nOpen Validation file ", f, "\n")
    setwd(parameters$Folders$folderCVVL)
    nome_arq_vl = paste(parameters$Dataset.Name, "-Split-Vl-", f,
                        ".csv", sep="")
    arquivo_vl = data.frame(read.csv(nome_arq_vl))

    ########################################################################
    #cat("\nOpen Test file ", f, "\n")
    setwd(parameters$Folders$folderCVTS)
    nome_arq_ts = paste(parameters$Dataset.Name, "-Split-Ts-", f,
                        ".csv", sep="")
    arquivo_ts = data.frame(read.csv(nome_arq_ts))

    # juntando treino com validaçao
    arquivo_tr2 = rbind(arquivo_tr, arquivo_vl)

    setwd(parameters$Folders$folderReports)
    sparcification = data.frame(read.csv("sparcification.csv"))
    sparcification = sparcification[f,]

    num.tr = length(parameters$Valid.TR)

    k = 0
    while(k<num.tr){

      cat("\n#=========================================================")
      cat("\n# TR = ", k)
      cat("\n#=========================================================")

      # "/dev/shm/j-GpositiveGO/Partitions/Split-1/Tr-1"
      FolderPartTr = paste(FolderPSplit, "/Tr-", k, sep="")

      # "/dev/shm/j-GpositiveGO/Communities/Split-1/Tr-1"
      FolderPartComm = paste(FolderSplitComm, "/Tr-", k, sep="")

      # "/dev/shm/j-GpositiveGO/Val-MaF1/Split-1/Tr-1"
      # FolderVTr = paste(FolderVSplit, "/Tr-", k, sep="")

      # "/dev/shm/j-GpositiveGO/Test-MaF1/Split-1/Tr-1"
      FolderTr = paste(FolderTSplit, "/Tr-", k, sep="")
      if(dir.create(FolderTr)==FALSE){dir.create(FolderTr)}

      # QUEM É A MELHOR PARTIÇÃO?
      setwd(parameters$Folders$folderReports)
      str = paste("tr-", k, "-Best-MacroF1.csv", sep="")
      best = data.frame(read.csv(str))
      best = best[f,]
      num.part = best$partition
      num.groups = best$partition

      # QUAL É O MÉTODO ESCOLHIDO?
      methods = filter(parameters$Choosed$all.methods.choosed, split == f)
      Tr = paste("tr-", k, sep="")
      methods = filter(methods, sparsification == Tr)


      # QUEM É A PARTIÇÃO?
      partitions = filter(parameters$Choosed$all.partitions.choosed, fold == f)
      partitions = filter(partitions, partitions$Tr == k)
      partitions = filter(partitions, partitions$method == methods$method)

      # ORGANIZANDO A PARTIÇÃO
      partition = partitions[,c(-1,-2,-3)]
      labels = partition$labels
      groups = partition[,num.part]
      partition = data.frame(labels, groups)

      g = 1
      while(g<=num.groups){

        cat("\n#=========================================================")
        cat("\n# Group = ", g)
        cat("\n#=========================================================")

        ###############################################################
        FolderTestGroup = paste(FolderTr, "/Group-", g, sep="")
        if(dir.exists(FolderTestGroup)== FALSE){dir.create(FolderTestGroup) }


        ###############################################################
        #cat("\nSpecific Group: ", g, "\n")
        grupoEspecifico = filter(partition, groups == g)


        ################################################################
        cat("\nTRAIN: Mount Group ", g, "\n")
        atributos_tr = arquivo_tr2[parameters$Dataset.Info$AttStart:parameters$Dataset.Info$AttEnd]
        n_a = ncol(atributos_tr)
        classes_tr = select(arquivo_tr2, grupoEspecifico$label)
        n_c = ncol(classes_tr)
        grupo_tr = cbind(atributos_tr, classes_tr)
        fim_tr = ncol(grupo_tr)


        #####################################################################
        #cat("\n\tTRAIN: Save Group", g, "\n")
        setwd(FolderTestGroup)
        nome_tr = paste(parameters$Dataset.Name, "-split-tr-", f,
                        "-group-", g, ".csv", sep="")
        write.csv(grupo_tr, nome_tr, row.names = FALSE)


        #####################################################################
        #cat("\n\tINICIO FIM TARGETS: ", g, "\n")
        inicio = parameters$Dataset.Info$LabelStart
        fim = fim_tr
        ifr = data.frame(inicio, fim)
        write.csv(ifr, "inicioFimRotulos.csv", row.names = FALSE)


        ####################################################################
        #cat("\n\tTRAIN: Convert Train CSV to ARFF ", g , "\n")
        nome_arquivo_2 = paste(parameters$Dataset.Name, "-split-tr-", f,
                               "-group-", g, ".arff", sep="")
        arg1Tr = nome_tr
        arg2Tr = nome_arquivo_2
        arg3Tr = paste(inicio, "-", fim, sep="")
        str = paste("java -jar ", parameters$Folders$folderUtils,
                    "/R_csv_2_arff.jar ", arg1Tr, " ", arg2Tr, " ",
                    arg3Tr, sep="")
        print(system(str))


        ##################################################################
        #cat("\n\tTRAIN: Verify and correct {0} and {1} ", g , "\n")
        arquivo = paste(FolderTestGroup, "/", arg2Tr, sep="")
        str0 = paste("sed -i 's/{0}/{0,1}/g;s/{1}/{0,1}/g' ", arquivo, sep="")
        print(system(str0))


        ######################################################################
        #cat("\n\tTEST: Mount Group: ", g, "\n")
        atributos_ts = arquivo_ts[parameters$Dataset.Info$AttStart:parameters$Dataset.Info$AttEnd]
        classes_ts = select(arquivo_ts, grupoEspecifico$label)
        grupo_ts = cbind(atributos_ts, classes_ts)
        fim_ts = ncol(grupo_ts)


        ######################################################################
        #cat("\n\tTEST: Save Group ", g, "\n")
        setwd(FolderTestGroup)
        nome_ts = paste(parameters$Dataset.Name, "-split-ts-", f, "-group-",
                        g, ".csv", sep="")
        write.csv(grupo_ts, nome_ts, row.names = FALSE)


        ####################################################################
        #cat("\n\tTEST: Convert CSV to ARFF ", g , "\n")
        nome_arquivo_3 = paste(parameters$Dataset.Name, "-split-ts-", f,
                               "-group-", g , ".arff", sep="")
        arg1Ts = nome_ts
        arg2Ts = nome_arquivo_3
        arg3Ts = paste(inicio, "-", fim, sep="")
        str = paste("java -jar ", parameters$Folders$folderUtils,
                    "/R_csv_2_arff.jar ", arg1Ts, " ", arg2Ts, " ",
                    arg3Ts, sep="")
        system(str)


        #####################################################################
        #cat("\n\tTEST: Verify and correct {0} and {1} ", g , "\n")
        arquivo = paste(FolderTestGroup, "/", arg2Ts, sep="")
        str0 = paste("sed -i 's/{0}/{0,1}/g;s/{1}/{0,1}/g' ", arquivo, sep="")
        cat("\n")
        print(system(str0))
        cat("\n")


        #####################################################################
        #cat("\nCreating .s file for clus")
        if(inicio == fim){

          nome_config = paste(parameters$Dataset.Name, "-split-", f, "-group-",
                              g, ".s", sep="")
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

          ###################################################################
          cat("\nExecute CLUS: ", g , "\n")
          nome_config2 = paste(FolderTestGroup, "/", nome_config, sep="")
          str = paste("java -jar ", parameters$Folders$folderUtils,
                      "/Clus.jar ", nome_config2, sep="")
          print(system(str))
          cat("\n")

        } else {

          nome_config = paste(parameters$Dataset.Name, "-split-", f,
                              "-group-", g, ".s", sep="")
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
          cat("\n")

        }

        ##################################################################
        #cat("\n\nOpen predictions")
        nomeDoArquivo = paste(FolderTestGroup, "/", parameters$Dataset.Name,
                              "-split-", f,"-group-", g,
                              ".test.pred.arff", sep="")
        predicoes = data.frame(foreign::read.arff(nomeDoArquivo))


        #####################################################################
        #cat("\nS\nPLIT PREDICTIS")
        if(inicio == fim){
          #cat("\n\nOnly one label in this group")

          ###################################################################
          #cat("\n\nSave Y_true")
          setwd(FolderTestGroup)
          classes = data.frame(predicoes[,1])
          names(classes) = colnames(predicoes)[1]
          write.csv(classes, "y_true.csv", row.names = FALSE)

          #################################################################
          #cat("\n\nSave Y_true")
          rot = paste("Pruned.p.", colnames(predicoes)[1], sep="")
          pred = data.frame(predicoes[,rot])
          names(pred) = colnames(predicoes)[1]
          setwd(FolderTestGroup)
          write.csv(pred, "y_predict.csv", row.names = FALSE)

          ####################################################################
          rotulos = c(colnames(classes))
          n_r = length(rotulos)
          gc()

        } else {

          ##############################################################
          #cat("\n\nMore than one label in this group")
          comeco = 1+(fim - inicio)


          ####################################################################
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
        oito = paste("Variance_RHE_1.csv")

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

      k = k + 1
      gc()
    } # end Tr

    f = f + 1
    gc()
  } # fim do for each

  gc()
  cat("\n###########################################################")
  cat("\n# TEST MACRO F1: Build and Test Hybrid Partitions End     #")
  cat("\n###########################################################")
  cat("\n\n\n\n")
}


##############################################################################
# FUNCTION GATHER PREDICTIONS - BUILD CONFUSION MATRIX                       #
#   Objective                                                                #
#   Parameters                                                               #
##############################################################################
maf1.gather.predicts <- function(parameters){

  f = 1
  gatherR <- foreach(f = 1:parameters$Number.Folds) %dopar%{
  # while(f<=parameters$Number.Folds){


    ########################################################################
    FolderRoot = "~/TCP-TR-H-Clus"
    FolderScripts = paste(FolderRoot, "/R", sep="")

    setwd(FolderScripts)
    source("utils.R")

    setwd(FolderScripts)
    source("libraries.R")

    ##########################################################################

    # "/dev/shm/j-GpositiveGO/Partitions/Split-1"
    FolderPSplit = paste(parameters$Folders$folderPartitions,
                         "/Split-", f, sep="")

    # "/dev/shm/j-GpositiveGO/Communities/Split-1"
    FolderSplitComm = paste(parameters$Folders$folderCommunities,
                            "/Split-", f, sep="")

    # "/dev/shm/j-GpositiveGO/Test-MaF1/Split-1"
    FolderTSplit = paste(parameters$Folders$folderTestMaF1,
                         "/Split-", f, sep="")

    setwd(parameters$Folders$folderReports)
    sparcification = data.frame(read.csv("sparcification.csv"))

    num.tr = length(parameters$Valid.TR)

    k = 0
    while(k<num.tr){

      # "/dev/shm/j-GpositiveGO/Partitions/Split-1/Tr-1"
      FolderPartTr = paste(FolderPSplit, "/Tr-", k, sep="")

      # "/dev/shm/j-GpositiveGO/Communities/Split-1/Tr-1"
      FolderPartComm = paste(FolderSplitComm, "/Tr-", k, sep="")

      # "/dev/shm/j-GpositiveGO/Test-MaF1/Split-1/Tr-1"
      FolderTr = paste(FolderTSplit, "/Tr-", k, sep="")

      # QUEM É A MELHOR PARTIÇÃO?
      setwd(parameters$Folders$folderReports)
      str = paste("tr-", k, "-Best-MacroF1.csv", sep="")
      best = data.frame(read.csv(str))
      best = best[f,]
      num.part = best$partition
      num.groups = best$partition

      # QUAL É O MÉTODO ESCOLHIDO?
      methods = filter(parameters$Choosed$all.methods.choosed, split == f)
      Tr = paste("tr-", k, sep="")
      methods = filter(methods, sparsification == Tr)

      # QUEM É A PARTIÇÃO?
      partitions = filter(parameters$Choosed$all.partitions.choosed, fold == f)
      partitions = filter(partitions, partitions$Tr == k)
      partitions = filter(partitions, partitions$method == methods$method)

      # ORGANIZANDO A PARTIÇÃO
      partition = partitions[,c(-1,-2,-3)]
      labels = partition$labels
      groups = partition[,num.part]
      partition = data.frame(labels, groups)

      ################################################################
      apagar = c(0)
      y_true = data.frame(apagar)
      y_pred = data.frame(apagar)

      # GROUP
      g = 1
      while(g<=num.groups){

        cat("\n\n#=========================================================")
        cat("\n# Fold \t\t", f)
        cat("\n# Tr \t\t", k)
        cat("\n# Group \t", g)
        cat("\n#=========================================================")

        FolderTestGroup = paste(FolderTr, "/Group-", g, sep="")

        #######################################################
        #cat("\nSpecific Group: ", g, "\n")
        grupoEspecifico = filter(partition, groups == g)

        #cat("\n\nGather y_true ", g)
        setwd(FolderTestGroup)
        #setwd(FolderTG)
        y_true_gr = data.frame(read.csv("y_true.csv"))
        y_true = cbind(y_true, y_true_gr)

        setwd(FolderTestGroup)
        #setwd(FolderTG)
        #cat("\n\nGather y_predict ", g)
        y_pred_gr = data.frame(read.csv("y_predict.csv"))
        y_pred = cbind(y_pred, y_pred_gr)

        #cat("\n\nDeleting files")
        #unlink("y_true.csv", recursive = TRUE)
        #unlink("y_predict.csv", recursive = TRUE)
        #unlink("inicioFimRotulos.csv", recursive = TRUE)

        g = g + 1
        gc()
      } # FIM DO GRUPO

      #cat("\n\nSave files ", g, "\n")
      setwd(FolderTr)
      y_pred = y_pred[,-1]
      y_true = y_true[,-1]
      write.csv(y_pred, "y_predict.csv", row.names = FALSE)
      write.csv(y_true, "y_true.csv", row.names = FALSE)

      k = k + 1
      gc()
    } # FIM DO Tr

    #f = f + 1
    gc()
  } # fim do foreach

  gc()
  cat("\n###############################################################")
  cat("\n# Gather Predicts: END                                        #")
  cat("\n###############################################################")
  cat("\n\n\n\n")

} # fim da função


##############################################################################
# FUNCTION EVALUATE TESTED HYBRID PARTITIONS                                 #
#   Objective                                                                #
#   Parameters                                                               #
##############################################################################
maf1.evaluate.test <- function(parameters){

  f = 1
  avalParal <- foreach(f = 1:parameters$Number.Folds) %dopar%{
    #while(f<=parameters$Number.Folds){

    cat("\n#=========================================================")
    cat("\n# Fold: ", f)
    cat("\n#=========================================================")

    ########################################################################
    FolderRoot = "~/TCP-TR-H-Clus"
    FolderScripts = paste(FolderRoot, "/R", sep="")

    setwd(FolderScripts)
    source("utils.R")

    setwd(FolderScripts)
    source("libraries.R")

    ###########################################################################

    # "/dev/shm/j-GpositiveGO/Partitions/Split-1"
    FolderPSplit = paste(parameters$Folders$folderPartitions,
                         "/Split-", f, sep="")

    # "/dev/shm/j-GpositiveGO/Communities/Split-1"
    FolderSplitComm = paste(parameters$Folders$folderCommunities,
                            "/Split-", f, sep="")


    # "/dev/shm/j-GpositiveGO/Test-MaF1/Split-1"
    FolderTSplit = paste(parameters$Folders$folderTestMaF1,
                         "/Split-", f, sep="")


    setwd(parameters$Folders$folderReports)
    sparcification = data.frame(read.csv("sparcification.csv"))

    num.tr = length(parameters$Valid.TR)

    # Tr
    k = 0
    while(k<num.tr){

      cat("\n#=========================================================")
      cat("\n# Tr = ", k)
      cat("\n#=========================================================")

      # "/dev/shm/j-GpositiveGO/Partitions/Split-1/Tr-1"
      FolderPartTr = paste(FolderPSplit, "/Tr-", k, sep="")

      # "/dev/shm/j-GpositiveGO/Communities/Split-1/Tr-1"
      FolderPartComm = paste(FolderSplitComm, "/Tr-", k, sep="")

      # "/dev/shm/j-GpositiveGO/Test-MaF1/Split-1/Tr-1"
      FolderTr = paste(FolderTSplit, "/Tr-", k, sep="")

      #cat("\nData frame")
      apagar = c(0)
      confMatPartitions = data.frame(apagar)
      partitions = c()

      #cat("\nGet the true and predict lables")
      setwd(FolderTr)
      #setwd(FolderTr)
      y_true = data.frame(read.csv("y_true.csv"))
      y_pred = data.frame(read.csv("y_predict.csv"))

      #cat("\nCompute measures multilabel")
      y_true2 = data.frame(sapply(y_true, function(x) as.numeric(as.character(x))))
      y_true3 = mldr_from_dataframe(y_true2 , labelIndices = seq(1,ncol(y_true2 )), name = "y_true2")
      y_pred2 = sapply(y_pred, function(x) as.numeric(as.character(x)))

      #cat("\nSave Confusion Matrix")
      setwd(FolderTr)
      #setwd(FolderTr)
      salva3 = paste("Conf-Mat-Fold-", f, "-Tr-", k, ".txt", sep="")
      sink(file=salva3, type="output")
      confmat = multilabel_confusion_matrix(y_true3, y_pred2)
      print(confmat)
      sink()

      #cat("\nCreating a data frame")
      confMatPart = multilabel_evaluate(confmat)
      confMatPart = data.frame(confMatPart)
      names(confMatPart) = paste("Fold-", f, "-Tr-", k, sep="")
      namae = paste("Split-", f, "-Tr-", k,"-Evaluated.csv", sep="")
      setwd(FolderTr)
      write.csv(confMatPart, namae)

      #cat("\nDelete files")
      setwd(FolderTr)
      #unlink("y_true.csv", recursive = TRUE)
      #unlink("y_predict.csv", recursive = TRUE)

      k = k + 1
      gc()
    } # FIM DO Tr

    #f = f + 1
    gc()
  } # fim do for each

  gc()
  cat("\n############################################################")
  cat("\n# TEST MACRO F1: Evaluation Folds END                      #")
  cat("\n############################################################")
  cat("\n\n\n\n")
}



##############################################################################
# FUNCTION GATHER EVALUATION                                                 #
#   Objective                                                                #
#   Parameters                                                               #
##############################################################################
maf1.gather.evaluation <- function(parameters){

  # vector with names
  measures = c("accuracy","average-precision","clp","coverage","F1",
               "hamming-loss","macro-AUC", "macro-F1","macro-precision",
               "macro-recall","margin-loss","micro-AUC","micro-F1",
               "micro-precision","micro-recall","mlp","one-error",
               "precision","ranking-loss", "recall","subset-accuracy","wlp")

  # from fold = 1 to number_folders
  f = 1
  while(f<=parameters$Number.Folds){

    ########################################################################
    FolderRoot = "~/TCP-TR-H-Clus"
    FolderScripts = paste(FolderRoot, "/R", sep="")

    setwd(FolderScripts)
    source("utils.R")

    setwd(FolderScripts)
    source("libraries.R")

    ###########################################################################

    # "/dev/shm/j-GpositiveGO/Partitions/Split-1"
    FolderPSplit = paste(parameters$Folders$folderPartitions,
                         "/Split-", f, sep="")

    # "/dev/shm/j-GpositiveGO/Communities/Split-1"
    FolderSplitComm = paste(parameters$Folders$folderCommunities,
                            "/Split-", f, sep="")

    # "/dev/shm/j-GpositiveGO/Test-MaF1/Split-1"
    FolderTSplit = paste(parameters$Folders$folderTestMaF1,
                         "/Split-", f, sep="")

    setwd(parameters$Folders$folderReports)
    sparcification = data.frame(read.csv("sparcification.csv"))
    n = mean(sparcification$minimum)

    #rm(avaliadoTr)

    apagar = c(0)
    avaliadoTr = data.frame(apagar)
    nomesThreshold = c("")

    num.tr = length(parameters$Valid.TR)

    k = 0
    while(k<num.tr){

      cat("\n#==========================")
      cat("\n# Fold \t", f)
      cat("\n# Tr \t", k)
      cat("\n#==========================")

      # "/dev/shm/j-GpositiveGO/Partitions/Split-1/Tr-1"
      FolderPartTr = paste(FolderPSplit, "/Tr-", k, sep="")

      # "/dev/shm/j-GpositiveGO/Communities/Split-1/Tr-1"
      FolderPartComm = paste(FolderSplitComm, "/Tr-", k, sep="")

      # "/dev/shm/j-GpositiveGO/Test-MaF1/Split-1/Tr-1"
      FolderTr = paste(FolderTSplit, "/Tr-", k, sep="")

      ######################################################################
      setwd(FolderTr)
      str = paste("Split-", f, "-Tr-", k, "-Evaluated.csv", sep="")
      avaliado = data.frame(read.csv(str))
      avaliadoTr = cbind(avaliadoTr, avaliado[,2])
      nomesThreshold[k] = paste("Fold-", f, "-Tr-", k, sep="")
      #names(avaliadoTr)[k+2] = nomesThreshold[k]
      #unlink(str, recursive = TRUE)

      k = k + 1
      gc()
    } # FIM DO Tr

    avaliadoTr = avaliadoTr[,-1]
    names(avaliadoTr) = nomesThreshold
    avaliadoTr = cbind(measures, avaliadoTr)
    setwd(FolderTSplit)
    write.csv(avaliadoTr, paste("Evaluated-Fold-", f, ".csv", sep=""),
              row.names = FALSE)

    f = f + 1
    gc()

  } # end folds


  gc()
  cat("\n###############################################################")
  cat("\n# TEST MACRO F1: Gather Evaluations End                       #")
  cat("\n###############################################################")
  cat("\n\n\n\n")
}



##############################################################################
# FUNCTION ORGANIZE EVALUATION                                               #
#   Objective                                                                #
#   Parameters                                                               #
##############################################################################
maf1.organize.evaluation <- function(parameters){

  dfs = list()
  dfs2 = list()

  x = 1
  while(x<=parameters$Number.Folds){
    dfs[[x]] = maf1.build.data.frame()
    x = x + 1
    gc()
  }

  # from fold = 1 to number_folders
  f = 1
  while(f<=parameters$Number.Folds){

    cat("\n#=========================================================")
    cat("\n# Fold: ", f)
    cat("\n#=========================================================")


    ########################################################################
    FolderRoot = "~/TCP-TR-H-Clus"
    FolderScripts = paste(FolderRoot, "/R", sep="")

    ###########################################################################
    # "/dev/shm/j-GpositiveGO/Test-MaF1/Split-1"
    FolderTSplit = paste(parameters$Folders$folderTestMaF1,
                         "/Split-", f, sep="")

    #########################################################
    #setwd(FolderTempTr)
    setwd(FolderTSplit)
    str = paste("Evaluated-Fold-", f, ".csv", sep="")
    dfs2[[f]] = data.frame(read.csv(str))

    #unlink(str, recursive = TRUE)

    f = f + 1
    gc()

  } # end folds

  numCol = ncol(dfs2[[1]])-1

  # vector with names
  measures = c("accuracy", "average-precision", "clp", "coverage",
               "F1", "hamming-loss", "macro-AUC", "macro-F1",
               "macro-precision", "macro-recall", "margin-loss",
               "micro-AUC", "micro-F1", "micro-precision", "micro-recall",
               "mlp", "one-error", "precision", "ranking-loss",
               "recall", "subset-accuracy", "wlp")
  apagar = c(0)
  nomesTr = c()
  nomes = c()

  k = 1
  while(k<=numCol){

    resultado = data.frame(measures, apagar)
    nomesFolds = c()
    a = k - 1
    nomeTr1 = paste("Evaluated-10Folds-Tr-", a, ".csv", sep="")
    nomeTr2 = paste("Mean-10Folds-Tr-", a, ".csv", sep="")
    nomeTr3 = paste("Median-10Folds-Tr-", a, ".csv", sep="")
    nomeTr4 = paste("SD-10Folds-Tr-", a, ".csv", sep="")

    f = 1
    while(f<=parameters$Number.Folds){

      cat("\n\n#====================================================")
      cat("\n# FOLD: ", f)
      cat("\n# Tr: ", k)
      cat("\n#====================================================")

      # pegando apenas o fold especifico
      res = data.frame(dfs2[[f]])
      nomesColunas = colnames(res)

      # pegando a partir da segunda coluna
      a = k + 1
      res = res[,a]

      resultado = cbind(resultado, res)
      b = ncol(resultado)
      names(resultado)[b] = nomesColunas[a]

      nomes[f] = paste("Fold-",f,"-Tr-", k, sep="")

      f = f + 1
      gc()
    } # fim do fold

    resultado = data.frame(resultado[,-2])
    setwd(parameters$Folders$folderTestMaF1)
    write.csv(resultado, nomeTr1, row.names = FALSE)

    # calculando a média dos 10 folds para cada medida
    media = data.frame(apply(resultado[,-1], 1, mean))
    media = cbind(measures, media)
    names(media) = c("Measures", "Mean10Folds")
    write.csv(media, nomeTr2, row.names = FALSE)

    mediana = data.frame(apply(resultado[,-1], 1, median))
    mediana = cbind(measures, mediana)
    names(mediana) = c("Measures", "Median10Folds")
    write.csv(mediana, nomeTr3, row.names = FALSE)

    dp = data.frame(apply(resultado[,-1], 1, sd))
    dp = cbind(measures, dp)
    names(dp) = c("Measures", "SD10Folds")
    write.csv(dp, nomeTr4, row.names = FALSE)

    k = k + 1
    gc()
  } # fim do k

  gc()
  cat("\n################################################################")
  cat("\n# TEST MACRO F1: Organize End                                  #")
  cat("\n################################################################")
  cat("\n\n\n\n")
}



##################################################################
#
##################################################################
maf1.test <- function(parameters){


  cat("\n\n#############################################################")
    cat("\n# TEST: maf1.build.test()                                   #")
    cat("\n#############################################################\n\n")
  timeBuild = system.time(resBT <- maf1.build.test(parameters))


  cat("\n\n############################################################")
    cat("\n# TEST: maf1.gather.predicts()                             #")
    cat("\n############################################################\n\n")
  timeSplit = system.time(resGather <- maf1.gather.predicts(parameters))


  cat("\n\n#############################################################")
    cat("\n# TEST: maf1.evaluate.test()                                #")
    cat("\n#############################################################\n\n")
  timeAvalia = system.time(resEval <- maf1.evaluate.test(parameters))


  cat("\n\n#############################################################")
    cat("\n# TEST: maf1.gather.evaluation()                            #")
     cat("\n############################################################\n\n")
  timeGather = system.time(resGE <- maf1.gather.evaluation(parameters))


  cat("\n\n#####################################################")
    cat("\n# TEST: maf1.organize.evaluation()                  #")
    cat("\n#####################################################\n\n")
  timeOrg = system.time(resOA <- maf1.organize.evaluation(parameters))


  cat("\n\n############################################################")
    cat("\n# RUN TEST MACRO F1: Save Runtime                          #")
    cat("\n############################################################\n\n")
  Runtime = rbind(timeBuild,
                  timeSplit,
                  timeAvalia,
                  timeGather,
                  timeOrg)
  setwd(parameters$Folders$folderTestSilho)
  write.csv(Runtime, paste(parameters$Dataset.Name,
                           "-test-maf1-runtime-1.csv", sep=""),
            row.names = FALSE)
}


#########################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com
# Thank you very much!
#########################################################################
