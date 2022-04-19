cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: HEAD                                             #")
getwd()
cat("\n#####################################################################\n\n")

rm(list=ls())

##################################################################################################
# TEST COMMUNITIES PARTITIONS                                                                    #
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


cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: SET WORK SPACE                                   #")
cat("\n#####################################################################\n\n")
FolderRoot = "~/TCP-TR-H/"
FolderScripts = paste(FolderRoot, "/R/", sep="")


cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: OPTIONS CONFIGURATIONS                           #")
cat("\n#####################################################################\n\n")
options(java.parameters = "-Xmx32g")
options(show.error.messages = TRUE)
options(scipen=20)


cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: READ DATASETS-2022                               #")
cat("\n#####################################################################\n\n")
setwd(FolderRoot)
datasets <- data.frame(read.csv("datasets-2022.csv"))


cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: GET THE ARGUMENTS COMMAND LINE                   #")
cat("\n#####################################################################\n\n")
args <- commandArgs(TRUE)


cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: DATASET INFORMATION                              #")
cat("\n#####################################################################\n\n")
ds <- datasets[args[1],]

number_dataset <- as.numeric(args[1])
cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: NUMBER DATASET: ", number_dataset, "             #")
cat("\n#####################################################################\n\n")


number_cores <- as.numeric(args[2])
cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: GET THE NUMBER CORES: ", number_cores, "         #")
cat("\n#####################################################################\n\n")


number_folds <- as.numeric(args[3])
cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: GET THE NUMBER FOLDS: ", number_folds, "         #")
cat("\n#####################################################################\n\n")


similarity <- toString(args[4])
cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: GET THE STRING SIMILARITY: ", similarity, "      #")
cat("\n#####################################################################\n\n")


folderResults <- toString(args[5])
cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: GET THE STRING FOLDER RESULTS: ", folderResults, " #")
cat("\n#####################################################################\n\n")


dataset_name <- toString(ds$Name)
cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: GET THE DATASET NAME: ", dataset_name, "         #")
cat("\n#####################################################################\n\n")


#ds <- datasets[22,]
#number_dataset = 22
#number_cores = 10
#number_folds = 10
#similarity = "Jaccard"
#folderResults = "/dev/shm/res"
#dataset_name = ds$Name


cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: CREATING FOLDER RESULTS TEMP                     #")
cat("\n#####################################################################\n\n")
if(dir.exists(folderResults)==FALSE){ dir.create(folderResults)}
cat("\n")

cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H:  LISTINGS                                       #")
cat("\n#####################################################################\n\n")
dir(folderResults)
cat("\n")


cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: LOAD SOURCES R                                   #")
cat("\n#####################################################################\n\n")

setwd(FolderScripts)
source("run.R")

setwd(FolderScripts)
source("utils.R")


cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: CREATING DIRECTORIES                             #")
cat("\n#####################################################################\n\n")
diretorios = directories(dataset_name, folderResults)


cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: ORGANIZANDO AS COISAS                            #")
cat("\n#####################################################################\n\n")

cat("\nCOPIANDO PARTIÇÕES")
str20 = paste("cp ~/TCP-TR-H/Partitions/", similarity ,"/", ds$Name,
              ".tar.gz ", diretorios$folderResults, sep="")
res = system(str20)
if(res!=0){break}else{cat("\ncopiou")}


cat("\nDESCOMPACTANDO PARTIÇÕES")
str21 = paste("tar xzf ", diretorios$folderResults, "/", ds$Name,
              ".tar.gz -C ", diretorios$folderResults, sep="")
res=system(str21)
if(res!=0){break}else{cat("\ndescompactou")}


cat("\nAPAGANDO TAR")
str22 = paste("rm ", diretorios$folderResults, "/",
              ds$Name, ".tar.gz", sep="")
res=system(str22)
if(res!=0){break}else{cat("\napagou")}


cat("\nCOPIANDO COMUNIDADES")
str23 = paste("cp -r ", diretorios$folderResults, "/", ds$Name,
              "/Communities/* ",  diretorios$folderResults,
              "/Communities/", sep="")
res=system(str23)
if(res!=0){break}else{cat("\ncopiou")}


cat("\nCOPIANDO PARTITIONS")
str24 = paste("cp -r ", diretorios$folderResults, "/", ds$Name,
              "/Partitions/* ", diretorios$folderResults,
              "/Partitions/",sep="")
res=system(str24)
if(res!=0){break}else{cat("\ncopiou")}


cat("\nAPAGANDO")
str25 = paste("rm -r ", diretorios$folderResults, "/", ds$Name, sep="")
print(system(str25))
if(res!=0){break}else{cat("\ncopiou")}


cat("\nCOPIANDO DATASETS")
str26 = paste("cp ~/TCP-TR-H/Datasets/", ds$Name, ".tar.gz ",
              diretorios$folderResults, "/datasets/", sep="")
res=system(str26)
if(res!=0){break}else{cat("\ncopiou")}


cat("\nDESCOMPACTANDO DATASETS")
str27 = paste("tar xzf ", diretorios$folderResults ,
              "/datasets/", ds$Name, ".tar.gz -C ",
              diretorios$folderResults, "/datasets", sep="")
res=system(str27)
if(res!=0){break}else{cat("\ndescompactou")}


cat("\n APAGANDO TAR")
str28 = paste("rm ", diretorios$folderResults, "/datasets/",
              ds$Name, ".tar.gz", sep="")
res=system(str28)
if(res!=0){break}else{cat("\napagou")}


cat("\nCOPIANDO CROSS VALIDATION")
str29 = paste("cp -r ", diretorios$folderDatasets, "/", ds$Name,
              "/CrossValidation/* ", diretorios$folderResults,
              "/datasets/CrossValidation/", sep="")
res=system(str29)
if(res!=0){break}else{cat("\ncopiou")}


cat("\nCOPIANDO CROSS LABEL SPACE")
str30 = paste("cp -r ",diretorios$folderDatasets, "/", ds$Name,
              "/LabelSpace/* ", diretorios$folderResults,
              "/datasets/LabelSpace/", sep="")
res=system(str30)
if(res!=0){break}else{cat("\ncopiou")}


cat("\nCOPIANDO CROSS NAMES LABEL")
str31 = paste("cp -r ", diretorios$folderDatasets, "/", ds$Name,
              "/NamesLabels/* ", diretorios$folderResults,
              "/datasets/NamesLabels/", sep="")
res=system(str31)
if(res!=0){break}else{cat("\ncopiou")}


cat("\nAPAGANDO PASTA")
str32 = paste("rm -r ", diretorios$folderResults,
              "/datasets/", ds$Name, sep="")
print(system(str32))
if(res!=0){break}else{cat("\napagou")}



cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: START                                            #")
cat("\n#####################################################################\n\n")
timeTCP = system.time(res <- execute(ds, dataset_name, number_dataset,
                                     number_folds, number_cores,
                                     folderResults, diretorios))


cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: SAVE RUNTIME                                     #")
cat("\n#####################################################################\n\n")
result_set <- t(data.matrix(timeTCP))
setwd(diretorios$folderTest)
write.csv(result_set, "Runtime.csv")
print(timeTCP)
cat("\n")

str2 = paste("rm -rf ", diretorios$folderDatasets, sep="")
print(system(str2))



cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: COMPRESS TEST RESULTS                             #")
cat("\n#####################################################################\n\n")
str = paste(" tar -zcf ", diretorios$folderTest ,"/test.tar.gz ",
            diretorios$folderReports, sep="")
print(system(str))


str2 = paste("cp ", diretorios$folderTest ,"/test.tar.gz ", diretorios$folderReports, sep="")
print(system(str2))

str2 = paste("rm -rf ", diretorios$folderTest, sep="")
print(system(str2))



cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: COMPRESS VALIDATE RESULTS                             #")
cat("\n#####################################################################\n\n")
str = paste(" tar -zcf ", diretorios$folderValidate ,"/test.tar.gz ",
            diretorios$folderReports, sep="")
print(system(str))

str2 = paste("cp ", diretorios$folderValidate ,"/test.tar.gz ",
             diretorios$folderReports, sep="")
print(system(str2))

str2 = paste("rm -rf ", diretorios$folderValidate, sep="")
print(system(str2))


#################################################################
str2 = paste("rm -rf ", diretorios$folderCommunities, sep="")
print(system(str2))

str2 = paste("rm -rf ", diretorios$folderPartitions, sep="")
print(system(str2))
################################################################



cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-NH: COPY TO HOME                                     #")
cat("\n#####################################################################\n\n")

str0 = "~/TCP-TR-H/Reports"
if(dir.exists(str0)==FALSE){dir.create(str0)}

str1 = paste(str0, "/", similarity, sep="")
if(dir.exists(str1)==FALSE){dir.create(str1)}

str2 = paste(str1, "/TR-H", sep="")
if(dir.exists(str2)==FALSE){dir.create(str2)}

str3 = paste(str2, "/", dataset_name, sep="")
if(dir.exists(str3)==FALSE){dir.create(str3)}

str4 = paste("cp -r ", diretorios$folderReports, "/* ", str3, sep="")
print(system(str4))

cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: CLEAN                                            #")
cat("\n#####################################################################\n\n")
str2 = paste("rm -rf ", diretorios$folderResults, sep="")
print(system(str2))
rm(list = ls())
gc()

cat("\n\n###################################################################")
cat("\n# ====> TCP-TR-H: END                                              #")
cat("\n#####################################################################\n")


##################################################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com                                   #
# Thank you very much!                                                                           #
##################################################################################################
