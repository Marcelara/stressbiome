#PLSDA gene expression and phytohormones

library(ropls)
library(readxl)
library(PERMANOVA)
library(santaR)

library(plsdepot)
library(ggplot2)
library(vegan)

Statistics_gene_expression_and_phytohormones <- read_excel("./Statistics gene expression and phytohormones.xlsx", 
                                                           sheet = "BSF Shoot")
View(Statistics_gene_expression_and_phytohormones)


####BSF shoot####
data<-Statistics_gene_expression_and_phytohormones
View(data)
matrixdata<-data[,-1:-4]#remove first 4 columns (sample names)
Metadata<-data[,2:4]
MetaData<-as.matrix(Metadata)


#PLS-DA may give better separation
matrixdataPLSDA<-opls(matrixdata,MetaData)#R2X=0.687,R2Y=0.0862,Q2=0.0764,RMSEE=5.72 pre=1, ort=0

matrixdataPLSDA<-opls(matrixdata,MetaData,scaleC = "pareto")#R2X=1,R2Y=0.14,Q2=0.125,RMSEE=5.94 pre=1, ort=0

matrixdataPLSDA<-opls(matrixdata,MetaData,scaleC = "center")#R2X=1,R2Y=0.159,Q2=0.143,RMSEE=5.94 pre=1, ort=0
#all single component models

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata,MetaData,predI=2)#R2X=0.766,R2Y=0.149,Q2=0.0749, RMSEE=5.29, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata,MetaData,predI=2,scaleC = "pareto")#R2X=1,R2Y=0.238,Q2=0.161, RMSEE=5.6, pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata,MetaData,predI=2,scaleC = "center")#R2X=1,R2Y=0.194,Q2=0.163, RMSEE=5.86, pre=2, ort=0
#matrixdataPLSDA gives best split of samples

scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata$Herbivore)  # Add Herbivore as a factor
scores_df$Time_point <- as.numeric(Metadata$Time_point)  # Add Time_point as numeric

#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore), size = factor(Time_point))) +
  geom_point() +
  labs(x = "t1 (69%)", y = "t2 (8%)") +
  scale_color_manual(values = c("black", "#B4C7E7", "#2F5597")) +  # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,16,17)) +  # Customize shape for Herbivore
  scale_size_manual(values = c(2, 6)) +  # Customize size for Time_point
  xlim(-10,10) + #set x-axis
  ylim(-5,5) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#nice, it works :D now add loading arrows

loadings<-matrixdataPLSDA@loadingMN
loadings_df<-as.data.frame(loadings)

loadings_df$p1 <- as.numeric(loadings_df$p1)
loadings_df$p2 <- as.numeric(loadings_df$p2)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*6),
                                                       yend = (p2*6)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*6), y = (p2*6), label = rownames(loadings_df),
                angle = 0.45, hjust = -.1, vjust = .2),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)



#extract VIP scores and CV-value
VIPmatrixdataPLSDA<-getVipVn(matrixdataPLSDA, orthoL = FALSE)
write.csv(VIPmatrixdataPLSDA)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata$Treatment<-as.factor(Metadata$Treatment)
Metadata$Herbivore<-as.factor(Metadata$Herbivore)
Metadata$Time_point<-as.factor(Metadata$Time_point)
adonis2(DistanceMatrix~Treatment*Herbivore*Time_point,data=Metadata,permutations = 999, method="bray")
#Error in terms.formula(formula, data = data) : invalid model formula in ExtractVars

##subset 2h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "BSF Shoot")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data2h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==2)
View(data2h)
data2h<-na.omit(data2h)
matrixdata2h<-data2h[,-1:-4]#remove first 4 columns (sample names)
Metadata2h<-data2h[,2:3]
Metadata2hMatrix<-as.matrix(Metadata2h)

#PLS-DA may give better separation
matrixdataPLSDA2h<-opls(matrixdata2h,Metadata2hMatrix)#no significant axis

matrixdataPLSDA2h<-opls(matrixdata2h,Metadata2hMatrix,scaleC = "pareto")#no significant axis

matrixdataPLSDA2h<-opls(matrixdata2h,Metadata2hMatrix,scaleC = "center")#no significant axis
#all no component models

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata2h,Metadata2hMatrix,predI=2)#R2X=0.726,R2Y=0.149,Q2=-0.00171, RMSEE=0.784, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata2h,Metadata2hMatrix,predI=2,scaleC = "pareto")#R2X=0.997,R2Y=0.07,Q2=-0.0275, RMSEE=0.82, pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata2h,Metadata2hMatrix,predI=2,scaleC = "center")#R2X=1,R2Y=0.0622,Q2=-0.033, RMSEE=0.823, pre=2, ort=0
#matrixdataPLSDA gives best split of samples

scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata2h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata2h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (58%)", y = "t2 (15%)") +
  scale_color_manual(values = c("black", "#B4C7E7", "#2F5597")) +  # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-6.5,6.5) + #set x-axis
  ylim(-2,2) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#looks good :) separate cluster of Plutella infested frass and exuviae treated plants to the rigth. 
#One big cluster to the left with Plutella control on top, uninfested frass and control in the middle and uninfested 
#control with most Delia infested plants in the bottom. Except for Delia exuviae which has 2 in the middle and 2 on top

#add loading arrows

loadings<-matrixdataPLSDA@loadingMN
loadings_df<-as.data.frame(loadings)


Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*2),
                                                       yend = (p2*2)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*2), y = (p2*2), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = -.1, vjust = .2),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata2h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata2h$Treatment<-as.factor(Metadata2h$Treatment)
Metadata2h$Herbivore<-as.factor(Metadata2h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata2h,permutations = 999, method="bray")

##subset 24h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "BSF Shoot")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data24h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==24)
View(data24h)
data24h<-na.omit(data24h)
matrixdata24h<-data24h[,-1:-4]#remove first 4 columns (sample names)
Metadata24h<-data24h[,2:3]
Metadata24hMatrix<-as.matrix(Metadata24h)

#PLS-DA may give better separation
matrixdataPLSDA2h<-opls(matrixdata24h,Metadata24hMatrix)#no significant axis

matrixdataPLSDA2h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "pareto")#no significant axis

matrixdataPLSDA2h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "center")#no significant axis
#all no component models

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2)#R2X=0.84,R2Y=0.285,Q2=0.0596, RMSEE=0.724, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2,scaleC = "pareto")#R2X=1,R2Y=0.105,Q2=-0.1, RMSEE=0.81, pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2,scaleC = "center")#R2X=1,R2Y=0.0729,Q2=-0.0305, RMSEE=0.825, pre=2, ort=0
#matrixdataPLSDA gives best split of samples


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata24h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata24h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (52%)", y = "t2 (32%)") +
  scale_color_manual(values = c("black", "#B4C7E7", "#2F5597")) +  # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-4,4) + #set x-axis
  ylim(-4,4) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#Nice separation of all plutella infested plants from the others and separate for soil treatments. Also separation for 
#uninfested plants per soil treatment. No clear pattern for Delia infested plants

#add loading arrows
loadings<-matrixdataPLSDA@loadingMN
loadings_df<-as.data.frame(loadings)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*5),
                                                       yend = (p2*5)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*5), y = (p2*5), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = 1, vjust = -.5),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata24h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata24h$Treatment<-as.factor(Metadata24h$Treatment)
Metadata24h$Herbivore<-as.factor(Metadata24h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata24h,permutations = 999, method="bray")

####BSF roots####
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "BSF Root")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data<-Copy_of_Statistics_gene_expression_and_phytohormones
View(data)
data<-na.omit(data)
matrixdata<-data[,-1:-4]#remove first 4 columns (sample names)
Metadata<-data[,2:4]
MetadataMatrix<-as.matrix(Metadata)

#PLS-DA may give better separation
matrixdataPLSDA<-opls(matrixdata,MetadataMatrix)#R2X=0.687,R2Y=0.394,Q2=0.345,RMSEE=4.6 pre=2, ort=0

matrixdataparetroPLSDA<-opls(matrixdata,MetadataMatrix,scaleC = "pareto")#no significant axis

matrixdatacenterPLSDA<-opls(matrixdata,MetadataMatrix,scaleC = "center")#no significant axis
#matrixdataPLSDA gives significant model, others are no component models


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata$Herbivore)  # Add Herbivore as a factor
scores_df$Time_point <- as.numeric(Metadata$Time_point)  # Add Time_point as a factor

#Plot the scores using ggplot2
ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore), size = factor(Time_point))) +
  geom_point() +
  labs(x = "t1 (56%)", y = "t2 (13%)") +
  scale_color_manual(values = c("black", "#B4C7E7", "#2F5597")) +  # Customize color scheme for Treatment
  scale_shape_manual(values = c(16, 17,18)) +  # Customize shape for Herbivore
  scale_size_manual(values = c(2, 4)) +  # Customize size for Time_point
  xlim(-6,6) + #set x-axis
  ylim(-3,3) + #set x-axis
  theme_minimal()
#Exuviae splits from the other soil amendments. No clear herbivore or time effect

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata$Treatment<-as.factor(Metadata$Treatment)
Metadata$Herbivore<-as.factor(Metadata$Herbivore)
Metadata$"Time-Point"<-as.factor(Metadata$Time_point)
adonis2(DistanceMatrix~Treatment*Herbivore*Time_point,data=Metadata,permutations = 999, method="bray")


##subset 2h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "BSF Root")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data2h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==2)
View(data2h)
data2h<-na.omit(data2h)
matrixdata2h<-data2h[,-1:-4]#remove first 4 columns (sample names)
Metadata2h<-data2h[,2:3]
Metadata2hMatrix<-as.matrix(Metadata2h)

#PLS-DA may give better separation
matrixdataPLSDA2h<-opls(matrixdata2h,Metadata2hMatrix)#R2X=0.555,R2Y=0.676,Q2=0.515,RMSEE=0.476 pre=2, ort=0

matrixdataPLSDA2hpareto<-opls(matrixdata2h,Metadata2hMatrix,scaleC = "pareto")#no significant axis

matrixdataPLSDA2hcenter<-opls(matrixdata2h,Metadata2hMatrix,scaleC = "center")#no significant axis
#matrixdataPLSDA gives significant model, others are no component models



scores<-matrixdataPLSDA2h@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata2h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata2h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (33%)", y = "t2 (22%)") +
  scale_color_manual(values = c("black", "#B4C7E7", "#2F5597")) +  # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-4.5,4.5) + #set x-axis
  ylim(-3,3) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#separation between exuviae and the other amendments. Herbivores unclear due to few replicates

#add loading arrows
loadings<-matrixdataPLSDA2h@loadingMN
loadings_df<-as.data.frame(loadings)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*5),
                                                       yend = (p2*5)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*5), y = (p2*5), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = 1, vjust = .5),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata2h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata2h$Treatment<-as.factor(Metadata2h$Treatment)
Metadata2h$Herbivore<-as.factor(Metadata2h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata2h,permutations = 999, method="bray")


##subset 24h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "BSF Root")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data24h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==24)
View(data24h)
data24h<-na.omit(data24h)
matrixdata24h<-data24h[,-1:-4]#remove first 4 columns (sample names)
Metadata24h<-data24h[,2:3]
Metadata24hMatrix<-as.matrix(Metadata24h)

#PLS-DA may give better separation
matrixdataPLSDA24h<-opls(matrixdata24h,Metadata24hMatrix)#R2X=0.636,R2Y=0.366,Q2=0.324,RMSEE=0.656 pre=1, ort=0

matrixdataPLSDAparetro24h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "pareto")#R2X=0.999,R2Y=0.349,Q2=0.271,RMSEE=0.672 pre=2, ort=0

matrixdataPLSDAcenter24h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "center")#R2X=1,R2Y=0.232,Q2=0.195,RMSEE=0.712 pre=1, ort=0
#Pareto scaling most significant axis, but is data really pareto distributed?

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2)#R2X=0.736,R2Y=0.484,Q2=0.336, RMSEE=0.602, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2,scaleC = "pareto")#R2X=0.999,R2Y=0.349,Q2=0.271,RMSEE=0.672 pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2,scaleC = "center")#R2X=1,R2Y=0.307,Q2=0.22, RMSEE=0.69, pre=2, ort=0
#matrixdataPLSDA gives best split of samples


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata24h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata24h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (64%)", y = "t2 (10%)") +
  scale_color_manual(values = c("black", "#B4C7E7", "#2F5597")) +  # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-5,5) + #set x-axis
  ylim(-3,3) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#clear separation of exuviae amendmend from frass and control. No clear herbivore effect
#add loading arrows
loadings<-matrixdataPLSDA@loadingMN
loadings_df<-as.data.frame(loadings)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*3),
                                                       yend = (p2*3)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*3), y = (p2*3), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = -.1, vjust = .2),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)


#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata24h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata24h$Treatment<-as.factor(Metadata24h$Treatment)
Metadata24h$Herbivore<-as.factor(Metadata24h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata24h,permutations = 999, method="bray")

####HC shoots####
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "HC shoot")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data<-Copy_of_Statistics_gene_expression_and_phytohormones
View(data)
data<-na.omit(data)
matrixdata<-data[,-1:-4]#remove first 4 columns (sample names)
Metadata<-data[,2:4]
MetadataMatrix<-as.matrix(Metadata)

#PLS-DA may give better separation
matrixdataPLSDA<-opls(matrixdata,MetadataMatrix)#R2X=0.798,R2Y=0.195,Q2=0.154,RMSEE=4.7 pre=2, ort=0

matrixdataparetroPLSDA<-opls(matrixdata,MetadataMatrix,scaleC = "pareto")#R2X=1,R2Y=0.439,Q2=0.294,RMSEE=4.69 pre=3, ort=0

matrixdatacenterPLSDA<-opls(matrixdata,MetadataMatrix,scaleC = "center")#R2X=1,R2Y=0.212,Q2=0.151,RMSEE=5.8 pre=2, ort=0
#matrixdataPLSDA gives best spread of data. But VSP2 has an enormous pull in all three models


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata$Herbivore)  # Add Herbivore as a factor
scores_df$Time_point <- as.numeric(Metadata$Time_point)  # Add Time_point as a factor

#Plot the scores using ggplot2
ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore), size = factor(Time_point))) +
  geom_point() +
  labs(x = "t1 (68%)", y = "t2 (12%)") +
  scale_color_manual(values = c("black", "#FDD8A1", "#FFC000")) +  # Customize color scheme for Treatment
  scale_shape_manual(values = c(16, 17,18)) +  # Customize shape for Herbivore
  scale_size_manual(values = c(2, 4)) +  # Customize size for Time_point
  #xlim(-6,6) + #set x-axis
  #ylim(-3,3) + #set x-axis
  theme_minimal()
#clear separation of Plutella 24 hpi, rest is all over the place

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata$Treatment<-as.factor(Metadata$Treatment)
Metadata$Herbivore<-as.factor(Metadata$Herbivore)
Metadata$"Time-Point"<-as.factor(Metadata$Time_point)
adonis2(DistanceMatrix~Treatment*Herbivore*Time_point,data=Metadata,permutations = 999, method="bray")



##subset 2h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "HC shoot")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data2h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==2)
View(data2h)
data2h<-na.omit(data2h)
matrixdata2h<-data2h[,-1:-4]#remove first 4 columns (sample names)
Metadata2h<-data2h[,2:3]
Metadata2hMatrix<-as.matrix(Metadata2h)

#PLS-DA may give better separation
matrixdataPLSDA2h<-opls(matrixdata2h,Metadata2hMatrix)#no significant axis

matrixdataPLSDA2hpareto<-opls(matrixdata2h,Metadata2hMatrix,scaleC = "pareto")#no significant axis

matrixdataPLSDA2hcenter<-opls(matrixdata2h,Metadata2hMatrix,scaleC = "center")#no significant axis
#all no component models

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata2h,Metadata2hMatrix,predI=2)#R2X=0.613,R2Y=0.112,Q2=-0.34, RMSEE=0.792, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata2h,Metadata2hMatrix,predI=2,scaleC = "pareto")#R2X=0.973,R2Y=0.0291,Q2=-0.2,RMSEE=0.827 pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata2h,Metadata2hMatrix,predI=2,scaleC = "center")#R2X=0.999,R2Y=0.0211,Q2=-0.191, RMSEE=0.831, pre=2, ort=0
#matrixdataPLSDA gives best split of samples


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata2h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata2h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (18%)", y = "t2 (43%)") +
  scale_color_manual(values = c("black", "#FDD8A1", "#FFC000")) + # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-4,4) + #set x-axis
  ylim(-4,4) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#no clear separation or grouping

#add loading arrows
loadings<-matrixdataPLSDA@loadingMN
loadings_df<-as.data.frame(loadings)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*6),
                                                       yend = (p2*6)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*6), y = (p2*6), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = .2, vjust = 1.4),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata2h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata2h$Treatment<-as.factor(Metadata2h$Treatment)
Metadata2h$Herbivore<-as.factor(Metadata2h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata2h,permutations = 999, method="bray")

##subset 24h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "HC shoot")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data24h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==24)
View(data24h)
data24h<-na.omit(data24h)
data24h<-data24h[-24,]
matrixdata24h<-data24h[,-1:-4]#remove first 4 columns (sample names)
Metadata24h<-data24h[,2:3]
Metadata24hMatrix<-as.matrix(Metadata24h)

#PLS-DA may give better separation
matrixdataPLSDA24h<-opls(matrixdata24h,Metadata24hMatrix)#no significant axis

matrixdataPLSDAparetro24h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "pareto")#no significant axis

matrixdataPLSDAcenter24h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "center")#no significant axis
#all no component models

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2)#R2X=0.889,R2Y=0.099,Q2=-0.0236, RMSEE=0.794, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2,scaleC = "pareto")#R2X=1,R2Y=0.0274,Q2=-0.109,RMSEE=0.824 pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2,scaleC = "center")#R2X=1,R2Y=0.0157,Q2=-0.083, RMSEE=0.829, pre=2, ort=0
#matrixdataPLSDA gives best split of samples


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata24h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata24h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (60%)", y = "t2 (29%)") +
  scale_color_manual(values = c("black", "#FDD8A1", "#FFC000")) +   # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-4.5,4.5) + #set x-axis
  ylim(-5.5,5.5) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#clear separation of Plutella infested plants, but soil amendment is mixed

#add loading arrows
loadings<-matrixdataPLSDA@loadingMN
loadings_df<-as.data.frame(loadings)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*4),
                                                       yend = (p2*4)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*4), y = (p2*4), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = .4, vjust = 1.3),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata24h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata24h$Treatment<-as.factor(Metadata24h$Treatment)
Metadata24h$Herbivore<-as.factor(Metadata24h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata24h,permutations = 999, method="bray")


####HC roots####
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "HC root")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data<-Copy_of_Statistics_gene_expression_and_phytohormones
View(data)
data<-na.omit(data)
matrixdata<-data[,-1:-4]#remove first 4 columns (sample names)
Metadata<-data[,2:4]
MetadataMatrix<-as.matrix(Metadata)

#PLS-DA may give better separation
matrixdataPLSDA<-opls(matrixdata,MetadataMatrix)#R2X=0.619,R2Y=0.29,Q2=0.21,RMSEE=4.8 pre=2, ort=0

matrixdataparetroPLSDA<-opls(matrixdata,MetadataMatrix,scaleC = "pareto")#no component model

matrixdatacenterPLSDA<-opls(matrixdata,MetadataMatrix,scaleC = "center")#no component model
#matrixdataPLSDA gives best spread of data.others are no component models


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata$Herbivore)  # Add Herbivore as a factor
scores_df$Time_point <- as.numeric(Metadata$Time_point)  # Add Time_point as a factor

#Plot the scores using ggplot2
ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore), size = factor(Time_point))) +
  geom_point() +
  labs(x = "t1 (51%)", y = "t2 (11%)") +
  scale_color_manual(values = c("black", "#FDD8A1", "#FFC000")) +  # Customize color scheme for Treatment
  scale_shape_manual(values = c(16, 17,18)) +  # Customize shape for Herbivore
  scale_size_manual(values = c(2, 4)) +  # Customize size for Time_point
  xlim(-7.5,7.5) + #set x-axis
  ylim(-3,3) + #set x-axis
  theme_minimal()
#clear separation of Delia 24 hpi, rest is all over the place

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata$Treatment<-as.factor(Metadata$Treatment)
Metadata$Herbivore<-as.factor(Metadata$Herbivore)
Metadata$"Time-Point"<-as.factor(Metadata$Time_point)
adonis2(DistanceMatrix~Treatment*Herbivore*Time_point,data=Metadata,permutations = 999, method="bray")


##subset 2h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "HC root")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data2h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==2)
View(data2h)
data2h<-na.omit(data2h)
matrixdata2h<-data2h[,-1:-4]#remove first 4 columns (sample names)
Metadata2h<-data2h[,2:3]
Metadata2hMatrix<-as.matrix(Metadata2h)

#PLS-DA may give better separation
matrixdataPLSDA2h<-opls(matrixdata2h,Metadata2hMatrix)#no significant axis

matrixdataPLSDA2hpareto<-opls(matrixdata2h,Metadata2hMatrix,scaleC = "pareto")#no significant axis

matrixdataPLSDA2hcenter<-opls(matrixdata2h,Metadata2hMatrix,scaleC = "center")#no significant axis
#all no component models

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata2h,Metadata2hMatrix,predI=2)#R2X=0.386,R2Y=0.265,Q2=-0.21, RMSEE=0.737, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata2h,Metadata2hMatrix,predI=2,scaleC = "pareto")#R2X=0.966,R2Y=0.121,Q2=-0.315,RMSEE=0.807 pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata2h,Metadata2hMatrix,predI=2,scaleC = "center")#R2X=0.999,R2Y=0.0981,Q2=-0.187, RMSEE=0.817, pre=2, ort=0
#matrixdataPLSDA gives best split of samples. OH_JA is removed from analyses because all samples have zero OH_JA


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata2h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata2h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (18%)", y = "t2 (21%)") +
  scale_color_manual(values = c("black", "#FDD8A1", "#FFC000")) + # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-3,3) + #set x-axis
  ylim(-5,5) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#no clear separation or grouping
#add loading arrows
loadings<-matrixdataPLSDA@loadingMN
loadings_df<-as.data.frame(loadings)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*3),
                                                       yend = (p2*3)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*3), y = (p2*3), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = -.1, vjust = .2),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata2h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata2h$Treatment<-as.factor(Metadata2h$Treatment)
Metadata2h$Herbivore<-as.factor(Metadata2h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata2h,permutations = 999, method="bray")

##subset 24h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "HC root")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data24h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==24)
View(data24h)
data24h<-na.omit(data24h)
data24h<-data24h[-24,]
matrixdata24h<-data24h[,-1:-4]#remove first 4 columns (sample names)
Metadata24h<-data24h[,2:3]
Metadata24hMatrix<-as.matrix(Metadata24h)

#PLS-DA may give better separation
matrixdataPLSDA24h<-opls(matrixdata24h,Metadata24hMatrix)#R2X=0.503,R2Y=0.361,Q2=0.296, RMSEE=0.675, pre=1, ort=0

matrixdataPLSDAparetro24h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "pareto")#R2X=0.973,R2Y=0.136,Q2=0.0809, RMSEE=0.79, pre=1, ort=0

matrixdataPLSDAcenter24h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "center")#R2X=0.999,R2Y=0.134,Q2=0.0717, RMSEE=0.791, pre=1, ort=0
#all single component models

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2)#R2X=0.676,R2Y=0.421,Q2=0.278, RMSEE=0.654, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2,scaleC = "pareto")#R2X=0.993,R2Y=0.187,Q2=0.0383,RMSEE=0.778 pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2,scaleC = "center")#R2X=1,R2Y=0.157,Q2=0.00972, RMSEE=0.794, pre=2, ort=0
#matrixdataPLSDA gives best split of samples


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata24h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata24h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (50%)", y = "t2 (17%)") +
  scale_color_manual(values = c("black", "#FDD8A1", "#FFC000")) +   # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-5,5) + #set x-axis
  ylim(-4.5,4.5) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#clear separation of Delia infested plants. Exuviae amendment forms a tight cluster, control to the bottom and exuviae more to the top

#add loading arrows
loadings<-matrixdataPLSDA@loadingMN
loadings_df<-as.data.frame(loadings)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*5),
                                                       yend = (p2*5)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*5), y = (p2*5), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = 1.1, vjust = .2),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata24h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata24h$Treatment<-as.factor(Metadata24h$Treatment)
Metadata24h$Herbivore<-as.factor(Metadata24h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata24h,permutations = 999, method="bray")


####MW shoots####
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "MW shoot")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data<-Copy_of_Statistics_gene_expression_and_phytohormones
View(data)
data<-na.omit(data)
matrixdata<-data[,-1:-4]#remove first 4 columns (sample names)
Metadata<-data[,2:4]
MetadataMatrix<-as.matrix(Metadata)

#PLS-DA may give better separation
matrixdataPLSDA<-opls(matrixdata,MetadataMatrix)#no significant axis

matrixdataparetroPLSDA<-opls(matrixdata,MetadataMatrix,scaleC = "pareto")#R2X=1,R2Y=0.127,Q2=0.109,RMSEE=6 pre=1, ort=0

matrixdatacenterPLSDA<-opls(matrixdata,MetadataMatrix,scaleC = "center")#R2X=1,R2Y=0.142,Q2=0.122,RMSEE=6 pre=1, ort=0
#no or single component models

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata,MetadataMatrix,predI=2)#R2X=0.647,R2Y=0.14,Q2=0.0598, RMSEE=5.98, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata,MetadataMatrix,predI=2,scaleC = "pareto")#R2X=1,R2Y=0.177,Q2=-0.0102,RMSEE=5.85 pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata,MetadataMatrix,predI=2,scaleC = "center")#R2X=1,R2Y=0.155,Q2=0.00435, RMSEE=6, pre=2, ort=0
#matrixdataPLSDA gives best split of samples

scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata$Herbivore)  # Add Herbivore as a factor
scores_df$Time_point <- as.numeric(Metadata$Time_point)  # Add Time_point as a factor

#Plot the scores using ggplot2
ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore), size = factor(Time_point))) +
  geom_point() +
  labs(x = "t1 (55%)", y = "t2 (10%)") +
  scale_color_manual(values = c("black", "#C5E0B4", "#548235")) +  # Customize color scheme for Treatment
  scale_shape_manual(values = c(16, 17,18)) +  # Customize shape for Herbivore
  scale_size_manual(values = c(2, 4)) +  # Customize size for Time_point
  xlim(-10,10) + #set x-axis
  ylim(-3.5,3.5) + #set x-axis
  theme_minimal()
#clear separation of Plutella 24 hpi, frass clearly separate, exuviae as well but closer to control

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata$Treatment<-as.factor(Metadata$Treatment)
Metadata$Herbivore<-as.factor(Metadata$Herbivore)
Metadata$"Time-Point"<-as.factor(Metadata$Time_point)
adonis2(DistanceMatrix~Treatment*Herbivore*Time_point,data=Metadata,permutations = 999, method="bray")

##subset 2h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "MW shoot")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data2h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==2)
View(data2h)
data2h<-na.omit(data2h)
matrixdata2h<-data2h[,-1:-4]#remove first 4 columns (sample names)
Metadata2h<-data2h[,2:3]
Metadata2hMatrix<-as.matrix(Metadata2h)

#PLS-DA may give better separation
matrixdataPLSDA2h<-opls(matrixdata2h,Metadata2hMatrix)#no significant axis

matrixdataPLSDA2hpareto<-opls(matrixdata2h,Metadata2hMatrix,scaleC = "pareto")#no significant axis

matrixdataPLSDA2hcenter<-opls(matrixdata2h,Metadata2hMatrix,scaleC = "center")#no significant axis
#all no component models

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata2h,Metadata2hMatrix,predI=2)#R2X=0.552,R2Y=0.227,Q2=-0.249, RMSEE=0.745, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata2h,Metadata2hMatrix,predI=2,scaleC = "pareto")#R2X=0.99,R2Y=0.168,Q2=-0.0449,RMSEE=0.773 pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata2h,Metadata2hMatrix,predI=2,scaleC = "center")#R2X=1,R2Y=0.116,Q2=-0.0453, RMSEE=0.796, pre=2, ort=0
#matrixdataPLSDA gives best split of samples


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata2h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata2h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (16%)", y = "t2 (36%)") +
  scale_color_manual(values = c("black", "#C5E0B4", "#548235")) + # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-4.5,4.5) + #set x-axis
  ylim(-4,4) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#Plutella exuviae to the top, but with one control Plutella

#add loading arrows
loadings<-matrixdataPLSDA@loadingMN
loadings_df<-as.data.frame(loadings)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*5),
                                                       yend = (p2*5)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*5), y = (p2*5), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = .2, vjust = -.6),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)

#PerMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata2h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata2h$Treatment<-as.factor(Metadata2h$Treatment)
Metadata2h$Herbivore<-as.factor(Metadata2h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata2h,permutations = 999, method="bray")

##subset 24h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "MW shoot")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data24h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==24)
View(data24h)
data24h<-na.omit(data24h)
data24h<-data24h[-24,]
matrixdata24h<-data24h[,-1:-4]#remove first 4 columns (sample names)
Metadata24h<-data24h[,2:3]
Metadata24hMatrix<-as.matrix(Metadata24h)

#PLS-DA may give better separation
matrixdataPLSDA24h<-opls(matrixdata24h,Metadata24hMatrix)#R2X=0.648,R2Y=0.0964,Q2=0.0665, RMSEE=0.82, pre=1, ort=0

matrixdataPLSDAparetro24h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "pareto")#no significant axis

matrixdataPLSDAcenter24h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "center")#no significant axis
#all no component models

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2)#R2X=0.788,R2Y=0.2,Q2=0.0776, RMSEE=0.789, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2,scaleC = "pareto")#R2X=1,R2Y=0.123,Q2=-0.0321,RMSEE=0.821 pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2,scaleC = "center")#R2X=1,R2Y=0.0819,Q2=-0.0471, RMSEE=0.838, pre=2, ort=0
#matrixdataPLSDA gives best split of samples


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata24h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata24h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (65%)", y = "t2 (14%)") +
  scale_color_manual(values = c("black", "#C5E0B4", "#548235")) +   # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-7.5,7.5) + #set x-axis
  ylim(-3.5,3.5) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#clear separation of Plutella infested plants, but soil amendment is mixed

#add loading arrows
loadings<-matrixdataPLSDA@loadingMN
loadings_df<-as.data.frame(loadings)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*4),
                                                       yend = (p2*4)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*4), y = (p2*4), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = .2, vjust = 1),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata24h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata24h$Treatment<-as.factor(Metadata24h$Treatment)
Metadata24h$Herbivore<-as.factor(Metadata24h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata24h,permutations = 999, method="bray")


####MW roots####
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "MW root")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data<-Copy_of_Statistics_gene_expression_and_phytohormones
View(data)
data<-na.omit(data)
matrixdata<-data[,-1:-4]#remove first 4 columns (sample names)
Metadata<-data[,2:4]
MetadataMatrix<-as.matrix(Metadata)

#PLS-DA may give better separation
matrixdataPLSDA<-opls(matrixdata,MetadataMatrix)#R2X=0.696,R2Y=0.213,Q2=0.151,RMSEE=5.7 pre=2, ort=0

matrixdataparetroPLSDA<-opls(matrixdata,MetadataMatrix,scaleC = "pareto")#no component model

matrixdatacenterPLSDA<-opls(matrixdata,MetadataMatrix,scaleC = "center")#no component model
#matrixdataPLSDA gives best spread of data.others are no component models


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata$Herbivore)  # Add Herbivore as a factor
scores_df$Time_point <- as.numeric(Metadata$Time_point)  # Add Time_point as a factor

#Plot the scores using ggplot2
ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore), size = factor(Time_point))) +
  geom_point() +
  labs(x = "t1 (56%)", y = "t2 (13%)") +
  scale_color_manual(values = c("black", "#C5E0B4", "#548235")) + # Customize color scheme for Treatment
  scale_shape_manual(values = c(16, 17,18)) +  # Customize shape for Herbivore
  scale_size_manual(values = c(2, 4)) +  # Customize size for Time_point
  xlim(-10,10) + #set x-axis
  ylim(-6,6) + #set x-axis
  theme_minimal()
#clear separation of Delia 24 hpi, rest is all over the place

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata$Treatment<-as.factor(Metadata$Treatment)
Metadata$Herbivore<-as.factor(Metadata$Herbivore)
Metadata$Time_point<-as.factor(Metadata$Time_point)
adonis2(DistanceMatrix~Treatment*Herbivore*Time_point,data=Metadata,permutations = 999, method="bray")

##subset 2h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "MW root")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data2h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==2)
View(data2h)
data2h<-na.omit(data2h)
matrixdata2h<-data2h[,-1:-4]#remove first 4 columns (sample names)
Metadata2h<-data2h[,2:3]
Metadata2hMatrix<-as.matrix(Metadata2h)

#PLS-DA may give better separation
matrixdataPLSDA2h<-opls(matrixdata2h,Metadata2hMatrix)#no significant axis

matrixdataPLSDA2hpareto<-opls(matrixdata2h,Metadata2hMatrix,scaleC = "pareto")#R2X=0.317,R2Y=0.163,Q2=0.0602,RMSEE=0.793 pre=1, ort=0

matrixdataPLSDA2hcenter<-opls(matrixdata2h,Metadata2hMatrix,scaleC = "center")#no significant axis
#all no or single component models. COOH_JA_Ile is zero in all samples, so removed from analysis

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata2h,Metadata2hMatrix,predI=2)#R2X=0.299,R2Y=0.539,Q2=0.194, RMSEE=0.6, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata2h,Metadata2hMatrix,predI=2,scaleC = "pareto")#R2X=0.599,R2Y=0.227,Q2=-0.0998,RMSEE=0.778 pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata2h,Metadata2hMatrix,predI=2,scaleC = "center")#R2X=0.683,R2Y=0.151,Q2=-0.0508, RMSEE=0.817, pre=2, ort=0
#matrixdataPLSDA gives best split of samples. OH_JA is removed from analyses because all samples have zero OH_JA


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata2h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata2h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (15%)", y = "t2 (15%)") +
  scale_color_manual(values = c("black", "#C5E0B4", "#548235")) + # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-3,3) + #set x-axis
  ylim(-3.5,3.5) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#no clear separation or grouping

#add loading arrows
loadings<-matrixdataPLSDA@loadingMN
loadings_df<-as.data.frame(loadings)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*4),
                                                       yend = (p2*4)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*4), y = (p2*4), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = .1, vjust = -.5),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)

#PerMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata2h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata2h$Treatment<-as.factor(Metadata2h$Treatment)
Metadata2h$Herbivore<-as.factor(Metadata2h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata2h,permutations = 999, method="bray")

##subset 24h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "MW root")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data24h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==24)
View(data24h)
data24h<-na.omit(data24h)
matrixdata24h<-data24h[,-1:-4]#remove first 4 columns (sample names)
Metadata24h<-data24h[,2:3]
Metadata24hMatrix<-as.matrix(Metadata24h)

#PLS-DA may give better separation
matrixdataPLSDA24h<-opls(matrixdata24h,Metadata24hMatrix)#R2X=0.626,R2Y=0.285,Q2=0.21, RMSEE=0.711, pre=1, ort=0

matrixdataPLSDAparetro24h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "pareto")#no significant axis

matrixdataPLSDAcenter24h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "center")#no significant axis
#all single or no component models

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2)#R2X=0.728,R2Y=0.34,Q2=-0.0686, RMSEE=0.695, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2,scaleC = "pareto")#R2X=0.997,R2Y=0.239,Q2=0.0619,RMSEE=0.744 pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2,scaleC = "center")#R2X=1,R2Y=0.232,Q2=0.0654, RMSEE=0.745, pre=2, ort=0
#matrixdataPLSDA gives best split of samples


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata24h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata24h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (63%)", y = "t2 (10%)") +
  scale_color_manual(values = c("black", "#C5E0B4", "#548235")) +   # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-7.5,7.5) + #set x-axis
  ylim(-2,2) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#clear separation of Delia infested plants. 

#add loading arrows
loadings<-matrixdataPLSDA@loadingMN
loadings_df<-as.data.frame(loadings)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*2),
                                                       yend = (p2*2)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*2), y = (p2*2), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = 1, vjust = -.3),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata24h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata24h$Treatment<-as.factor(Metadata24h$Treatment)
Metadata24h$Herbivore<-as.factor(Metadata24h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata24h,permutations = 999, method="bray")



####BSF shoot####
Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                           sheet = "BSF SA")
View(Statistics_gene_expression_and_phytohormones)
data<-Statistics_gene_expression_and_phytohormones
View(data)
data<-na.omit(data)
matrixdata<-data[,-1:-4]#remove first 4 columns (sample names)
Metadata<-data[,2:4]
MetaData<-as.matrix(Metadata)


#PLS-DA may give better separation
matrixdataPLSDA<-opls(matrixdata,MetaData)#R2X=0.179,R2Y=0.207,Q2=0.155,RMSEE=9.98 pre=1, ort=0

matrixdataPartetoPLSDA<-opls(matrixdata,MetaData,scaleC = "pareto")#no significant axis

matrixdataCenterPLSDA<-opls(matrixdata,MetaData,scaleC = "center")#no significant axis
#all single component models

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata,MetaData,predI=2)#R2X=0.331,R2Y=0.288,Q2=0.155, RMSEE=8.14, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata,MetaData,predI=2,scaleC = "pareto")#R2X=0.548,R2Y=0.465,Q2=0.015, RMSEE=10.2, pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata,MetaData,predI=2,scaleC = "center")#R2X=0.834,R2Y=0.171,Q2=-0.123, RMSEE=12.9, pre=2, ort=0
#matrixdataPLSDA gives best split of samples

scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata$Herbivore)  # Add Herbivore as a factor
scores_df$Time_point <- as.numeric(Metadata$Time_point)  # Add Time_point as numeric

#Plot the scores using ggplot2
ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore), size = factor(Time_point))) +
  geom_point() +
  labs(x = "t1 (18%)", y = "t2 (15%)") +
  scale_color_manual(values = c("black", "#B4C7E7", "#2F5597","#203864")) +  # Customize color scheme for Treatment
  scale_shape_manual(values = c(16, 17)) +  # Customize shape for Herbivore
  scale_size_manual(values = c(2, 6)) +  # Customize size for Time_point
  xlim(-4.5,4.5) + #set x-axis
  ylim(-5,5) + #set x-axis
  theme_minimal()
#nice, it works :D now fix the axis and maybe split per time-point. I think there is no need to change the metadata to numbers, but need to check


#extract VIP scores and CV-value
VIPmatrixdataPLSDA<-getVipVn(matrixdataPLSDA, orthoL = FALSE)
write.csv(VIPmatrixdataPLSDA)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata$Treatment<-as.factor(Metadata$Treatment)
Metadata$Herbivore<-as.factor(Metadata$Herbivore)
Metadata$Time_point<-as.factor(Metadata$Time_point)
adonis2(DistanceMatrix~Treatment*Herbivore*Time_point,data=Metadata,permutations = 999, method="bray")


##subset 24h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "BSF SA")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data24h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==24)
View(data24h)
data24h<-na.omit(data24h)
matrixdata24h<-data24h[,-1:-4]#remove first 4 columns (sample names)
Metadata24h<-data24h[,2:3]
Metadata24hMatrix<-as.matrix(Metadata24h)

#PLS-DA may give better separation
matrixdataPLSDA24h<-opls(matrixdata24h,Metadata24hMatrix)#no significant axis

matrixdataParetoPLSDA24h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "pareto")#no significant axis

matrixdataCenterPLSDA24h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "center")#no significant axis
#all no component models

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2)#R2X=0.342,R2Y=0.313,Q2=-0.181, RMSEE=0.774, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2,scaleC = "pareto")#R2X=0.894,R2Y=0.0761,Q2=-0.208, RMSEE=0.867, pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2,scaleC = "center")#R2X=0.932,R2Y=0.0633,Q2=-0.204, RMSEE=0.872, pre=2, ort=0
#matrixdataPLSDA gives best split of samples

scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata24h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata24h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (15%)", y = "t2 (19%)") +
  scale_color_manual(values = c("black", "#B4C7E7", "#2F5597","#203864")) +  # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-4,4) + #set x-axis
  ylim(-6,6) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#no clear separation

#add loading arrows
loadings<-matrixdataPLSDA@loadingMN
loadings_df<-as.data.frame(loadings)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*6),
                                                       yend = (p2*6)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*6), y = (p2*6), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = .6, vjust = 1.2),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata24h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata24h$Treatment<-as.factor(Metadata24h$Treatment)
Metadata24h$Herbivore<-as.factor(Metadata24h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata24h,permutations = 999, method="bray")

##subset 72h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "BSF SA")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data72h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==72)
View(data72h)
data72h<-na.omit(data72h)
matrixdata72h<-data72h[,-1:-4]#remove first 4 columns (sample names)
Metadata72h<-data72h[,2:3]
Metadata72hMatrix<-as.matrix(Metadata72h)

#PLS-DA may give better separation
matrixdataPLSDA2h<-opls(matrixdata72h,Metadata72hMatrix)#R2X=0.297,R2Y=0.225,Q2=0.112, RMSEE=0.791, pre=1, ort=0

matrixdataPLSDA2h<-opls(matrixdata72h,Metadata72hMatrix,scaleC = "pareto")#no significant axis

matrixdataPLSDA2h<-opls(matrixdata72h,Metadata72hMatrix,scaleC = "center")#no significant axis
#all single or no component models

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata72h,Metadata72hMatrix,predI=2)#R2X=0.413,R2Y=0.4,Q2=0.0717, RMSEE=0.694, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata72h,Metadata72hMatrix,predI=2,scaleC = "pareto")#R2X=0.798,R2Y=0.26,Q2=-0.037, RMSEE=0.743, pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata72h,Metadata72hMatrix,predI=2,scaleC = "center")#R2X=0.9,R2Y=0.236,Q2=0.0316, RMSEE=0.752, pre=2, ort=0
#matrixdataPLSDA gives best split of samples


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata72h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata72h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (30%)", y = "t2 (12%)") +
  scale_color_manual(values = c("black", "#B4C7E7", "#2F5597","#203864")) + # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-6,6) + #set x-axis
  ylim(-3,3) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#no clear patterns

#add loading arrows
loadings<-matrixdataPLSDA@loadingMN
loadings_df<-as.data.frame(loadings)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*4),
                                                       yend = (p2*4)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*4), y = (p2*4), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = .5, vjust = -.2),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata72h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata72h$Treatment<-as.factor(Metadata72h$Treatment)
Metadata72h$Herbivore<-as.factor(Metadata72h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata72h,permutations = 999, method="bray")

####HC shoot####
Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                           sheet = "HC SA")
View(Statistics_gene_expression_and_phytohormones)
data<-Statistics_gene_expression_and_phytohormones
View(data)
data<-na.omit(data)
matrixdata<-data[,-1:-4]#remove first 4 columns (sample names)
Metadata<-data[,2:4]
MetaData<-as.matrix(Metadata)


#PLS-DA may give better separation
matrixdataPLSDA<-opls(matrixdata,MetaData)#R2X=0.466,R2Y=0.477,Q2=0.34,RMSEE=7.83 pre=3, ort=0

matrixdataPartetoPLSDA<-opls(matrixdata,MetaData,scaleC = "pareto")#R2X=0.43,R2Y=0.629,Q2=0.574,RMSEE=8.13 pre=1, ort=0

matrixdataCenterPLSDA<-opls(matrixdata,MetaData,scaleC = "center")#R2X=0.643,R2Y=0.653,Q2=0.604,RMSEE=8.31 pre=1, ort=0
#matrixdataPLSDA best model


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata$Herbivore)  # Add Herbivore as a factor
scores_df$Time_point <- as.numeric(Metadata$Time_point)  # Add Time_point as numeric

#Plot the scores using ggplot2
ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore), size = factor(Time_point))) +
  geom_point() +
  labs(x = "t1 (21%)", y = "t2 (13%)") +
  scale_color_manual(values = c("black", "#767171", "#FBDC85","#FFC000")) +  # Customize color scheme for Treatment
  scale_shape_manual(values = c(16, 17)) +  # Customize shape for Herbivore
  scale_size_manual(values = c(2, 6)) +  # Customize size for Time_point
  xlim(-6,6) + #set x-axis
  ylim(-3.5,3.5) + #set x-axis
  theme_minimal()
#nice, it works :D now fix the axis and maybe split per time-point. I think there is no need to change the metadata to numbers, but need to check


#extract VIP scores and CV-value
VIPmatrixdataPLSDA<-getVipVn(matrixdataPLSDA, orthoL = FALSE)
write.csv(VIPmatrixdataPLSDA)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata$Treatment<-as.factor(Metadata$Treatment)
Metadata$Herbivore<-as.factor(Metadata$Herbivore)
Metadata$Time_point<-as.factor(Metadata$Time_point)
adonis2(DistanceMatrix~Treatment*Herbivore*Time_point,data=Metadata,permutations = 999, method="bray")


##subset 24h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "HC SA")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data24h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==24)
View(data24h)
data24h<-na.omit(data24h)
matrixdata24h<-data24h[,-1:-4]#remove first 4 columns (sample names)
Metadata24h<-data24h[,2:3]
Metadata24hMatrix<-as.matrix(Metadata24h)

#PLS-DA may give better separation
matrixdataPLSDA24h<-opls(matrixdata24h,Metadata24hMatrix)#R2X=0.266,R2Y=0.335,Q2=0.273, RMSEE=0.592, pre=1, ort=0

matrixdataParetoPLSDA24h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "pareto")#R2X=0.608,R2Y=0.443,Q2=0.281, RMSEE=0.625, pre=2, ort=0

matrixdataCenterPLSDA24h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "center")#R2X=0.737,R2Y=0.361,Q2=0.203, RMSEE=0.705, pre=2, ort=0


#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata24h,Metadata24hMatrix,predI=2)#R2X=0.432,R2Y=0.442,Q2=0.289, RMSEE=0.569, pre=2, ort=0
#matrixdataPLSDA gives best split of samples

scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata24h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata24h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (27%)", y = "t2 (16%)") +
  scale_color_manual(values = c("black", "#767171", "#FBDC85","#FFC000")) + # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-6,6) + #set x-axis
  ylim(-4,4) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#no clear separation

#add loading arrows
loadings<-matrixdataPLSDA@loadingMN
loadings_df<-as.data.frame(loadings)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*5),
                                                       yend = (p2*5)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*5), y = (p2*5), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = .5, vjust = -.3),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata24h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata24h$Treatment<-as.factor(Metadata24h$Treatment)
Metadata24h$Herbivore<-as.factor(Metadata24h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata24h,permutations = 999, method="bray")

##subset 72h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "HC SA")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data72h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==72)
View(data72h)
data72h<-na.omit(data72h)
matrixdata72h<-data72h[,-1:-4]#remove first 4 columns (sample names)
Metadata72h<-data72h[,2:3]
Metadata72hMatrix<-as.matrix(Metadata72h)

#PLS-DA may give better separation
matrixdataPLSDA72h<-opls(matrixdata72h,Metadata72hMatrix)#no significant axis

matrixdataPLSDA72h<-opls(matrixdata72h,Metadata72hMatrix,scaleC = "pareto")#no significant axis

matrixdataPLSDA72h<-opls(matrixdata72h,Metadata72hMatrix,scaleC = "center")#no significant axis
#all no component models

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata72h,Metadata72hMatrix,predI=2)#R2X=0.345,R2Y=0.456,Q2=0.0661, RMSEE=0.703, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata72h,Metadata72hMatrix,predI=2,scaleC = "pareto")#R2X=0.696,R2Y=0.292,Q2=-0.0202, RMSEE=0.785, pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata72h,Metadata72hMatrix,predI=2,scaleC = "center")#R2X=0.912,R2Y=0.171,Q2=-0.0401, RMSEE=0.819, pre=2, ort=0
#matrixdataPLSDA gives best split of samples


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata72h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata72h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (17%)", y = "t2 (17%)") +
  scale_color_manual(values = c("black", "#767171", "#FBDC85","#FFC000")) + # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-5.5,5.5) + #set x-axis
  ylim(-4,4) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#no clear patterns

#add loading arrows
loadings<-matrixdataPLSDA@loadingMN
loadings_df<-as.data.frame(loadings)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*4),
                                                       yend = (p2*4)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*4), y = (p2*4), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = 0, vjust = 1),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata72h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata72h$Treatment<-as.factor(Metadata72h$Treatment)
Metadata72h$Herbivore<-as.factor(Metadata72h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata72h,permutations = 999, method="bray")

####MW shoot####
Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                           sheet = "MW SA")
View(Statistics_gene_expression_and_phytohormones)
data<-Statistics_gene_expression_and_phytohormones
View(data)
data<-na.omit(data)
matrixdata<-data[,-1:-4]#remove first 4 columns (sample names)
Metadata<-data[,2:4]
MetaData<-as.matrix(Metadata)


#PLS-DA may give better separation
matrixdataPLSDA<-opls(matrixdata,MetaData)#R2X=0.378,R2Y=0.376,Q2=0.294,RMSEE=6.74 pre=2, ort=0

matrixdataPartetoPLSDA<-opls(matrixdata,MetaData,scaleC = "pareto")#R2X=0.548,R2Y=0.376,Q2=0.345,RMSEE=11 pre=2, ort=0

matrixdataCenterPLSDA<-opls(matrixdata,MetaData,scaleC = "center")#no significant axis
#matrixdataPLSDA best model


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata$Herbivore)  # Add Herbivore as a factor
scores_df$Time_point <- as.numeric(Metadata$Time_point)  # Add Time_point as numeric

#Plot the scores using ggplot2
ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore), size = factor(Time_point))) +
  geom_point() +
  labs(x = "t1 (20%)", y = "t2 (18%)") +
  scale_color_manual(values = c("black", "#AFABAB", "#C5E0B4","#548235")) +  # Customize color scheme for Treatment
  scale_shape_manual(values = c(16, 17)) +  # Customize shape for Herbivore
  scale_size_manual(values = c(2, 6)) +  # Customize size for Time_point
  xlim(-5,5) + #set x-axis
  ylim(-6.5,6.5) + #set x-axis
  theme_minimal()
#nice, it works :D now fix the axis and maybe split per time-point. I think there is no need to change the metadata to numbers, but need to check


#extract VIP scores and CV-value
VIPmatrixdataPLSDA<-getVipVn(matrixdataPLSDA, orthoL = FALSE)
write.csv(VIPmatrixdataPLSDA)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata$Treatment<-as.factor(Metadata$Treatment)
Metadata$Herbivore<-as.factor(Metadata$Herbivore)
Metadata$Time_point<-as.factor(Metadata$Time_point)
adonis2(DistanceMatrix~Treatment*Herbivore*Time_point,data=Metadata,permutations = 999, method="bray")


##subset 24h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "MW SA")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data24h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==24)
View(data24h)
data24h<-na.omit(data24h)
matrixdata24h<-data24h[,-1:-4]#remove first 4 columns (sample names)
Metadata24h<-data24h[,2:3]
Metadata24hMatrix<-as.matrix(Metadata24h)

#PLS-DA may give better separation
matrixdataPLSDA24h<-opls(matrixdata24h,Metadata24hMatrix)#R2X=0.416,R2Y=0.367,Q2=0.197, RMSEE=0.645, pre=2, ort=0

matrixdataParetoPLSDA24h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "pareto")#No significant axis

matrixdataCenterPLSDA24h<-opls(matrixdata24h,Metadata24hMatrix,scaleC = "center")#No significant axis
#matrixdataPLSDA gives best split of samples. COOH_JA_Ile droped from analysis, since all values are 0

scores<-matrixdataPLSDA24h@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata24h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata24h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (28%)", y = "t2 (14%)") +
  scale_color_manual(values = c("black", "#AFABAB", "#C5E0B4","#548235")) +  # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-6,6) + #set x-axis
  ylim(-4,4) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#no clear separation

#add loading arrows
loadings<-matrixdataPLSDA24h@loadingMN
loadings_df<-as.data.frame(loadings)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*5),
                                                       yend = (p2*5)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*5), y = (p2*5), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = 1, vjust = -.3),
            data = loadings_df,
            colour = "black",inherit.aes = FALSE)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata24h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata24h$Treatment<-as.factor(Metadata24h$Treatment)
Metadata24h$Herbivore<-as.factor(Metadata24h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata24h,permutations = 999, method="bray")

##subset 72h
Copy_of_Statistics_gene_expression_and_phytohormones <- read_excel("aaWerk/jaar 5/Herbivore performance gene expression and phytohormone alaysis/Copy of Statistics gene expression and phytohormones.xlsx", 
                                                                   sheet = "MW SA")
View(Copy_of_Statistics_gene_expression_and_phytohormones)
data72h<-subset(Copy_of_Statistics_gene_expression_and_phytohormones,Time_point==72)
View(data72h)
data72h<-na.omit(data72h)
matrixdata72h<-data72h[,-1:-4]#remove first 4 columns (sample names)
Metadata72h<-data72h[,2:3]
Metadata72hMatrix<-as.matrix(Metadata72h)

#PLS-DA may give better separation
matrixdataPLSDA72h<-opls(matrixdata72h,Metadata72hMatrix)#R2X=0.257,R2Y=0.355,Q2=0.318, RMSEE=0.802, pre=1, ort=0

matrixdataPLSDA72h<-opls(matrixdata72h,Metadata72hMatrix,scaleC = "pareto")#R2X=0.467,R2Y=0.135,Q2=0.0946, RMSEE=0.814, pre=1, ort=0

matrixdataPLSDA72h<-opls(matrixdata72h,Metadata72hMatrix,scaleC = "center")#no significant axis
#all no component models

#request 2 predictive axis, even when not significant
matrixdataPLSDA<-opls(matrixdata72h,Metadata72hMatrix,predI=2)#R2X=0.375,R2Y=0.531,Q2=0.297, RMSEE=0.677, pre=2, ort=0

matrixdataparetoPLSDA<-opls(matrixdata72h,Metadata72hMatrix,predI=2,scaleC = "pareto")#R2X=0.807,R2Y=0.227,Q2=0.0507, RMSEE=0.81, pre=2, ort=0

matrixdatacenterPLSDA<-opls(matrixdata72h,Metadata72hMatrix,predI=2,scaleC = "center")#R2X=0.907,R2Y=0.161,Q2=-0.241, RMSEE=0.825, pre=2, ort=0
#matrixdataPLSDA gives best split of samples


scores<-matrixdataPLSDA@scoreMN
scores_df<-as.data.frame(scores)

scores_df$p1 <- as.numeric(scores_df$p1)
scores_df$p2 <- as.numeric(scores_df$p2)
scores_df$Treatment <- as.numeric(Metadata72h$Treatment)  # Add Treatment as a factor
scores_df$Herbivore <- as.numeric(Metadata72h$Herbivore)  # Add Herbivore as a factor


#Plot the scores using ggplot2
scoreplot<-ggplot(scores_df, aes(x = p1, y = p2, color = factor(Treatment), shape = factor(Herbivore))) +
  geom_point(size=4) +
  labs(x = "t1 (26%)", y = "t2 (12%)") +
  scale_color_manual(values = c("black", "#AFABAB", "#C5E0B4","#548235")) +  # Customize color scheme for Treatment
  scale_shape_manual(values = c(15,19,17)) +  # Customize shape for Herbivore
  xlim(-4,4) + #set x-axis
  ylim(-2.5,2.5) + #set x-axis
  theme_minimal()+
  theme(axis.text.x = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),
        axis.text.y = element_text( size = 12, angle = 0, hjust = 0, vjust = 0, face = "plain"),  
        axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = 1, face = "plain"))
#no clear patterns

#add loading arrows
loadings<-matrixdataPLSDA@loadingMN
loadings_df<-as.data.frame(loadings)

Biplot<-scoreplot+geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = (p1*3),
                                                       yend = (p2*3)), arrow = arrow(length = unit(1, "picas")),
                               colour = "black",inherit.aes = FALSE)+
  geom_text(aes(x = (p1*3), y = (p2*3), label = rownames(loadings_df),
                size = 12, angle = 0.45, hjust = -.08, vjust = -0.1),
            data = loadings_df,
            colour = "gray32",inherit.aes = FALSE)

#PERMANOVA
UVmatrixdata<-santaR:::scaling_UV(matrixdata72h)
DistanceMatrix<-dist(UVmatrixdata,method="euclidean")
Metadata72h$Treatment<-as.factor(Metadata72h$Treatment)
Metadata72h$Herbivore<-as.factor(Metadata72h$Herbivore)
adonis2(DistanceMatrix~Treatment*Herbivore,Metadata72h,permutations = 999, method="bray")
