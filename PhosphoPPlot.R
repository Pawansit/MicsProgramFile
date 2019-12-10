library(corrplot)
A <- read.csv("Phopho-S-InD.csv", header = F)
B <- read.csv("Phopho-S-InD_Relative.csv",header = F)
C <- read.csv("Phopho-T-InD.csv",header = F)
D <- read.csv("Phopho-T-InD_Relative.csv",header = F)


############### Array Plot (Serine) #############################

E <- read.csv("Phopho-S-Array.csv",header = T)
G <- read.csv("Phopho-S-Array_Relative.csv",header = T)
E <- as.matrix(E)
G <- as.matrix(G)
tiff("Sp-array.tiff",width = 30,height = 25,units = "cm", res = 600)
par(mfrow = c(2, 1)) 
barplot(E[1,1:20],col = rainbow(20))
barplot(G[1,1:20],col = rainbow(20))
dev.off()


################## Array Plot (Theonine) ##########################

I <- read.csv("Phopho-T-Array.csv",header = T)
J <- read.csv("Phopho-T-Array_Relative.csv",header = T)
I <- as.matrix(I)
J <- as.matrix(J)
tiff("Tp-array.tiff",width = 30,height = 25,units = "cm", res = 600)
par(mfrow = c(2, 1)) 
barplot(I[1,1:20],col = rainbow(20))
barplot(J[1,1:20],col = rainbow(20))
dev.off()




###########
col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white",
                           "cyan", "#007FFF", "blue","#00007F"))

col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))

col3 <- colorRampPalette(c("red", "white", "blue"))

col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F",
                           "cyan", "#007FFF", "blue","#00007F"))
col5 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F", 
                           "cyan", "#007FFF", "blue","#00007F")) 
col6 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white", 
                           "cyan", "#007FFF", "blue","#00007F"))


#############

corrplotMatrix <- function (A,Name)
{
  A <- as.matrix(A)
  As <- apply(A[2:21,2:14], 2,as.numeric)
  Ass <- as.data.frame(As)
  row.names(Ass)<- A[2:21,1]
  colnames(Ass) <- c(-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6)
  tiff(Name,width = 30,height = 25,units = "cm", res = 600)
  corrplot(as.matrix(Ass),is.cor = FALSE, col = col4(7))
  dev.off()
}


corrplotMatrix(A,"Sp-Original.tiff")
corrplotMatrix(B,"Sp-Relative.tiff")
corrplotMatrix(C,"Tp-Original.tiff")
corrplotMatrix(D,"Tp-Relative.tiff")


