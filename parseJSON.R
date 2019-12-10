#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("jsonlite")
library("methods")

result_json <- fromJSON(args)

########### Node Annotation ###########

dfID <- as.data.frame(result_json$elements$nodes$data$id)
colnames(dfID) <- "IDname"
dfx <- as.data.frame(result_json$elements$nodes$position$x)
colnames(dfx) <- ("positionX")
dfy <- as.data.frame(result_json$elements$nodes$position$y)
colnames(dfy) <- ("positionY")
dfname <- as.data.frame(result_json$elements$nodes$data$name)
colnames(dfname) <- "Name"
df <- data.frame(dfID, dfx, dfy,dfname)

file1 <- unlist(strsplit(args, split = ".cyjs"))
file2 <- "_Nodes.csv"
filename  <- paste(file1,file2,sep = "",collapse = NULL)

write.csv(x = df,file = filename)

############ Edge Annotation #########

Source <- as.data.frame(result_json$elements$edges$data$source)
colnames(Source) <- "Source"
target <- as.data.frame(result_json$elements$edges$data$target)
colnames(target) <- "target"
InterA <- as.data.frame(result_json$elements$edges$data$name)
colnames(InterA) <- "Interction"

dfEDGE <- data.frame(Source,target,InterA)

file3 <- "_Edges.csv"
filenameE <- paste(file1,file3,sep = "",collapse = NULL)

write.csv(x = dfEDGE,file = filenameE)
