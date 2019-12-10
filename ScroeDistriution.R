library(ggplot2)
require(gridExtra)


################# Distribution of the Confedence Score in the String database for Plasmodium, Human and Mouse 

A <- read.table("Human/ScoreHuman.txt",header = T)
B <- read.table("Plasmodium/Score_Plas.txt",header = T)
C <- read.table("Mouse/Scroemouse.txt",header = T)



Hu <- hist(A$combined_score)
Pl <- hist(B$combined_score)
Mu <- hist(C$combined_score)

BCount <- c(Hu$counts[1],Pl$counts[1],Mu$counts[1],Hu$counts[2],Pl$counts[2],Mu$counts[2],Hu$counts[3],Pl$counts[3],Mu$counts[3],Hu$counts[4],Pl$counts[4],Mu$counts[4],Hu$counts[5],Pl$counts[5],Mu$counts[5],Hu$counts[6],Pl$counts[6],Mu$counts[6],Hu$counts[7],Pl$counts[7],Mu$counts[7],Hu$counts[8],Pl$counts[8],Mu$counts[8],Hu$counts[9],Pl$counts[9],Mu$counts[9],Hu$counts[10],Pl$counts[10],Mu$counts[10])

df <- data.frame(Probe = rep(c("100","200","300","400","500","600","700","800","900","990"),each=3), S2=c("Human","Plasmodium","Mouse","Human","Plasmodium","Mouse","Human","Plasmodium","Mouse","Human","Plasmodium","Mouse","Human","Plasmodium","Mouse","Human","Plasmodium","Mouse","Human","Plasmodium","Mouse","Human","Plasmodium","Mouse","Human","Plasmodium","Mouse","Human","Plasmodium","Mouse"), Count= BCount)


tiff("Score_Distribution.tiff", width = 6, height = 5, units = "in", res = 300)

ggplot(data=df, aes(x=Probe, y=Count, fill=S2)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()

dev.off()

############### Distribution of the Interaction Type in Plasmodium, Human and Mouse ####################

D <- read.table("Plasmodium/5833.protein.links.full.v11.0.txt",sep = ' ',header = T)
E <- read.table("Human/9606.protein.links.full.v11.0.txt",sep = ' ',header = T)
G <- read.table("Mouse/10090.protein.links.full.v11.0.txt",sep = ' ',header = T)

Df_Plas <- Distribution(D,100)
Df_Hum  <- Distribution(E,100) 
Df_Mpu  <- Distribution(G,100)


Percentage <-  c(Df_Plas$Count[1],Df_Hum$Count[1],Df_Mpu$Count[1],Df_Plas$Count[2],Df_Hum$Count[2],Df_Mpu$Count[2],Df_Plas$Count[3],Df_Hum$Count[3],Df_Mpu$Count[3],Df_Plas$Count[4],Df_Hum$Count[4],Df_Mpu$Count[4],Df_Plas$Count[5],Df_Hum$Count[5],Df_Mpu$Count[5],Df_Plas$Count[6],Df_Hum$Count[6],Df_Mpu$Count[6],Df_Plas$Count[7],Df_Hum$Count[7],Df_Mpu$Count[7],Df_Plas$Count[8],Df_Hum$Count[8],Df_Mpu$Count[8],Df_Plas$Count[9],Df_Hum$Count[9],Df_Mpu$Count[9],Df_Plas$Count[10],Df_Hum$Count[10],Df_Mpu$Count[10],Df_Plas$Count[11],Df_Hum$Count[11],Df_Mpu$Count[11],Df_Plas$Count[12],Df_Hum$Count[12],Df_Mpu$Count[12],Df_Plas$Count[13],Df_Hum$Count[13],Df_Mpu$Count[13])
Property <- rep(c("neighborhood","neighborhood_T","fusion","cooccurence","homology","coexpression","coexpression_T","experiments","experiments_T","database","database_T","textmining","textmining_T"),each=3)
Species <- c("Plasmodium","Human","Mouse","Plasmodium","Human","Mouse","Plasmodium","Human","Mouse","Plasmodium","Human","Mouse","Plasmodium","Human","Mouse","Plasmodium","Human","Mouse","Plasmodium","Human","Mouse","Plasmodium","Human","Mouse","Plasmodium","Human","Mouse","Plasmodium","Human","Mouse","Plasmodium","Human","Mouse","Plasmodium","Human","Mouse","Plasmodium","Human","Mouse")

df_new <- data.frame(Property,Species,Percentage)

tiff("Property_Distribution.tiff", width = 10, height = 5, units = "in", res = 300)
ggplot(df_new,aes(x=Property, y=Percentage, fill=Species, label =Percentage)) + 
  geom_bar(stat="identity", position=position_dodge())+
  geom_col(position  = position_dodge2(width = 0.9, preserve  = "single"), show.legend = F) +
  geom_text(position = position_dodge2(width = 0.9, preserve = "single"), angle = 90, vjust=0.25, hjust=0) +
  ylim(0,10)+
  theme(text = element_text(size=8))
dev.off()
  

Distribution <- function(record,cutoff)
{
  S   <- dim(record) 
  Nb  <- record$neighborhood[D$neighborhood >= cutoff]
  NbT <- record$neighborhood_transferred[D$neighborhood_transferred >= cutoff]
  Fu  <- record$fusion[D$fusion >= cutoff]
  Co  <- record$cooccurence[D$cooccurence >= cutoff]
  Ho  <- record$homology[D$homology >= cutoff]
  Ce  <- record$coexpression[D$coexpression >= cutoff ]
  CeT <- record$coexpression_transferred[D$coexpression_transferred >=cutoff]
  Ex  <- record$experiments[D$experiments >= cutoff]
  ExT <- record$experiments_transferred[D$experiments_transferred >= cutoff]
  Da  <- record$database[D$database >= cutoff]
  DaT <- record$database_transferred[D$database_transferred >= cutoff]
  Tx  <- record$textmining[D$textmining >= cutoff]
  TxT <- record$textmining_transferred[D$textmining_transferred >= cutoff]
  CoB <- record$combined_score[D$combined_score >= cutoff]
  
  Aa1  <- round(((length(Nb)/S[1])*100),1)
  Aa2  <- round(((length(NbT)/S[1])*100),1)
  Aa3  <- round(((length(Fu)/S[1])*100),1)
  Aa4  <- round(((length(Co)/S[1])*100),1)
  Aa5  <- round(((length(Ho)/S[1])*100),1)
  Aa6  <- round(((length(Ce)/S[1])*100),1)
  Aa7  <- round(((length(CeT)/S[1])*100),1)
  Aa8  <- round(((length(Ex)/S[1])*100),1)
  Aa9  <- round(((length(ExT)/S[1])*100),1)
  Aa10 <- round(((length(Da)/S[1])*100),1)
  Aa11 <- round(((length(DaT)/S[1])*100),1)
  Aa12 <- round(((length(Tx)/S[1])*100),1)
  Aa13 <- round(((length(TxT)/S[1])*100),1)
  Aa14 <- round(((length(CoB)/S[1])*100),1)
  
  Property  <- c("neighborhood","neighborhood_transferred","fusion","cooccurence","homology","coexpression","coexpression_transferred","experiments","experiments_transferred","database","database_transferred","textmining","textmining_transferred","combined_score")
  Count     <- c(Aa1,Aa2,Aa3,Aa4,Aa5,Aa6,Aa7,Aa8,Aa9,Aa10,Aa11,Aa12,Aa13,Aa14)
  df <- data.frame(Property, Count)
  return(df)
}











