#
##
### Supporting R code for publication in Ecology and Evolution
### "Which fruit traits correlate reliably with sugar content, 
###  and does this vary with dispersal syndrome?"
### Authors: Amanda C. Luzardo, Alice C. Poirier, Laís A. A. Moreira, 
### Allegra N. DePasquale, Linh M. N. Nguyen, Omer Nevo & Amanda D. Melin
##
#

# Libraries (here using R 4.4.1)
library(devtools)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(dplyr)
library(lme4)
library(lme4)
library(DHARMa)
library(corrplot)
library(qgraph)


# for colour package and accessing jcolors 
# install.packages("devtools")
devtools::install_github("jaredhuling/jcolors")
library(jcolors)
jcolors('default')




### 1) Datasets ----
# Banana dataset: Luzardo-et-al_EcolEvol_Dataset_Bananas_Oct2024.csv
ba.dat<- read.csv(file.choose(), header=T, na.strings="NA", sep=",")
ba.dat[,c(1:5)]<- lapply(ba.dat[,c(1:5)], factor)
ba.dat<- mutate(ba.dat, Shannon.index.std=Shannon.index/Weight)
str(ba.dat) 

# Mango dataset: Luzardo-et-al_EcolEvol_Dataset_Mangoes_Oct2024.csv
ma.dat<- read.csv(file.choose(), header=T, na.strings="NA", sep=",")
ma.dat[,c(1:5)]<- lapply(ma.dat[,c(1:5)], factor)
ma.dat<- mutate(ma.dat, Shannon.index.std=Shannon.index/Weight)
str(ma.dat) 

# Currant dataset: Luzardo-et-al_EcolEvol_Dataset_Currants_Oct2024.csv
cu.dat<- read.csv(file.choose(), header=T, na.strings="NA", sep=",")
cu.dat[,c(1:6)]<- lapply(cu.dat[,c(1:6)], factor)
cu.dat<- mutate(cu.dat, Shannon.index.std=Shannon.index/Sum.weight)
str(cu.dat) 



# Remove spoiled fruits for subsequent analyses (only bananas)
###
# mean and median sugar content per day
ba.mean<- data.frame(ba.dat[which(ba.dat$Sugar.content!="NA"),] %>% group_by(Day) %>%
                       summarise(mean.sug=mean(Sugar.content), med.sug=median(Sugar.content)))
ba.mean

max(ba.mean$mean.sug)   #21.79 = Day 4
0.85*max(ba.mean$mean.sug) #18.52
0.95*max(ba.mean$mean.sug) #20.70

# so remove everything after Day 10
ba.dat.good<- filter(ba.dat, as.numeric(Day)<11)




### 2) Correlation plots ----
## Figure 1 ----
ba.day.sugar<-
  ggplot(ba.dat, aes(x=as.numeric(Day), y=Sugar.content))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  scale_x_discrete(name="Day", limits=c("Day01","Day02","Day03","Day04","Day05","Day06","Day07","Day08","Day09","Day10","Day11","Day12","Day13","Day14","Day15","Day16"),
                   labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"))+
  ylab("Sugar content (% Brix)")+
  labs(title=expression(paste(bold("Bananas"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="right")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))

#ma.day.sugar<-
  ggplot(ma.dat, aes(x=as.numeric(Day), y=Sugar.content))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  scale_x_discrete(name="Day", limits=c("Day01","Day02","Day03","Day04","Day05","Day06","Day07","Day08","Day09","Day10","Day11","Day12","Day13","Day14","Day15","Day16"),
                   labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"))+
  ylab("Sugar content (% Brix)")+
  labs(title=expression(paste(bold("Mangoes"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="right")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black")) 

cu.day.sugar<-
  ggplot(cu.dat, aes(x=as.numeric(Day), y=Sugar.content))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  scale_x_discrete(name="Day", limits=c("Day01","Day02","Day03","Day04","Day05","Day06","Day07","Day08","Day09","Day10","Day11","Day12","Day13","Day14","Day15","Day16"),
                   labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"))+
  ylab("Sugar content (% Brix)")+
  labs(title=expression(paste(bold("Currants"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="right")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black")) 



## Suppl. Figs S4-S6 ----
# Correlation fruit colour variables and days: y = Primate.blue-yellow.contrast
# Replace y variable with Primate.red-green.contrast, Primate.luminance, Passerine.blue-yellow.contrast, Passerine.red-green.contrast, or Passerine.luminance, as needed.
ba.day.BY<-
  ggplot(ba.dat.good, aes(x=as.numeric(Day), y=Primate.blue.yellow.contrast))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  scale_x_discrete(name="", limits=c("Day01","Day02","Day03","Day04","Day05","Day06","Day07","Day08","Day09","Day10"),
                   labels=c("1","2","3","4","5","6","7","8","9","10"))+
  xlab("Day")+
  ylab("Blue-yellow contrast (primate)")+
  labs(title=expression(paste(bold("Bananas"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))

ma.day.BY<-
  ggplot(ma.dat, aes(x=as.numeric(Day), y=Primate.blue.yellow.contrast))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  scale_x_discrete(name="Day", limits=c("Day01","Day02","Day03","Day04","Day05","Day06","Day07","Day08","Day09","Day10","Day11","Day12","Day13","Day14","Day15","Day16"),
                   labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"))+
  xlab("Day")+
  ylab("Blue-yellow contrast (primate)")+
  labs(title=expression(paste(bold("Mangoes"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black")) 

cu.day.BY<-
  ggplot(cu.dat, aes(x=as.numeric(Day), y=Primate.blue.yellow.contrast))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  scale_x_discrete(name="Day", limits=c("Day01","Day02","Day03","Day04","Day05","Day06","Day07","Day08","Day09","Day10","Day11","Day12","Day13","Day14","Day15","Day16"),
                   labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"))+
  xlab("Day")+
  ylab("Blue-yellow contrast (primate)")+
  labs(title=expression(paste(bold("Currants"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black")) 


# Correlation fruit colour variables and sugar content: y = Primate.blue-yellow.contrast
ba.sug.BY<-
  ggplot(ba.dat.good, aes(x=Sugar.content, y=Primate.blue.yellow.contrast))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  ylab("Blue-yellow contrast (primate)")+
  xlab("Sugar content (% Brix)")+
  labs(title=expression(paste(bold("Bananas"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))

ma.sug.BY<-
  ggplot(ma.dat, aes(x=Sugar.content, y=Primate.blue.yellow.contrast))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  ylab("Blue-yellow contrast (primate)")+
  xlab("Sugar content (% Brix)")+
  labs(title=expression(paste(bold("Mangoes"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black")) 

cu.sug.BY<-
  ggplot(cu.dat, aes(x=Sugar.content, y=Primate.blue.yellow.contrast))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  ylab("Blue-yellow contrast (primate)")+
  xlab("Sugar content (% Brix)")+
  labs(title=expression(paste(bold("Currants"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black")) 



## Suppl. Fig S7 ----
# Correlation fruit odour complexity and days
ba.day.odo<-
  ggplot(ba.dat.good, aes(x=as.numeric(Day), y=Shannon.index.std))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  scale_x_discrete(name="", limits=c("Day01","Day02","Day03","Day04","Day05","Day06","Day07","Day08","Day09","Day10"),
                   labels=c("1","2","3","4","5","6","7","8","9","10"))+
  ylab("Odour complexity (Shannon index)")+
  labs(title=expression(paste(bold("Bananas"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))

ma.day.odo<-
  ggplot(ma.dat, aes(x=as.numeric(Day), y=Shannon.index.std))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  scale_x_discrete(name="Day", limits=c("Day01","Day02","Day03","Day04","Day05","Day06","Day07","Day08","Day09","Day10","Day11","Day12","Day13","Day14","Day15","Day16"),
                   labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"))+
  ylab("Odour complexity (Shannon index)")+
  labs(title=expression(paste(bold("Mangoes"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black")) 

cu.day.odo<-
  ggplot(cu.dat, aes(x=as.numeric(Day), y=Shannon.index.std))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  scale_x_discrete(name="Day", limits=c("Day01","Day02","Day03","Day04","Day05","Day06","Day07","Day08","Day09","Day10","Day11","Day12","Day13","Day14","Day15","Day16"),
                   labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"))+
  ylab("Odour complexity (Shannon index)")+
  labs(title=expression(paste(bold("Currants"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black")) 


# Correlation fruit odour complexity and sugar content
ba.sug.odo<-
  ggplot(ba.dat.good, aes(x=Sugar.content, y=Shannon.index.std))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  ylab("Odour complexity (Shannon index)")+
  xlab("Sugar content (% Brix)")+
  labs(title=expression(paste(bold("Bananas"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))

ma.sug.odo<-
  ggplot(ma.dat, aes(x=Sugar.content, y=Shannon.index.std))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  ylab("Odour complexity (Shannon index)")+
  xlab("Sugar content (% Brix)")+
  labs(title=expression(paste(bold("Mangoes"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black")) 

cu.sug.odo<-
  ggplot(cu.dat, aes(x=Sugar.content, y=Shannon.index.std))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  ylab("Odour complexity (Shannon index)")+
  xlab("Sugar content (% Brix)")+
  labs(title=expression(paste(bold("Currants"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black")) 



## Suppl. Fig S8 ----
# Correlation fruit hardness and days
ba.day.instant<-
  ggplot(ba.dat.good, aes(x=as.numeric(Day), y=Instant.modulus.test))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  scale_x_discrete(name="", limits=c("Day01","Day02","Day03","Day04","Day05","Day06","Day07","Day08","Day09","Day10"),
                   labels=c("1","2","3","4","5","6","7","8","9","10"))+
  ylab("Hardness (instant modulus test MPa)")+
  labs(title=expression(paste(bold("Bananas"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))

ma.day.instant<-
  ggplot(ma.dat, aes(x=as.numeric(Day), y=Instant.modulus.test))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  scale_x_discrete(name="Day", limits=c("Day01","Day02","Day03","Day04","Day05","Day06","Day07","Day08","Day09","Day10","Day11","Day12","Day13","Day14","Day15","Day16"),
                   labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"))+
  ylab("Hardness (instant modulus test MPa)")+
  labs(title=expression(paste(bold("Mangoes"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black")) 

cu.day.instant<-
  ggplot(cu.dat, aes(x=as.numeric(Day), y=Instant.modulus.test))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  scale_x_discrete(name="Day", limits=c("Day01","Day02","Day03","Day04","Day05","Day06","Day07","Day08","Day09","Day10","Day11","Day12","Day13","Day14","Day15","Day16"),
                   labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"))+
  ylab("Hardness (instant modulus test MPa)")+
  labs(title=expression(paste(bold("Currants"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
       axis.text.y=element_text(colour="black")) 
  

# Correlation fruit hardness and sugar content
ba.sug.instant<-
  ggplot(ba.dat.good, aes(x=Sugar.content, y=Instant.modulus.test))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  ylab("Hardness (instant modulus test MPa)")+
  xlab("Sugar content (% Brix)")+
  labs(title=expression(paste(bold("Bananas"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))

ma.sug.instant<-
  ggplot(ma.dat, aes(x=Sugar.content, y=Instant.modulus.test))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  ylab("Hardness (instant modulus test MPa)")+
  xlab("Sugar content (% Brix)")+
  labs(title=expression(paste(bold("Mangoes"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black")) 

cu.sug.instant<-
  ggplot(cu.dat, aes(x=Sugar.content, y=Instant.modulus.test))+
  geom_point(aes(color=Ripeness.stage))+
  geom_smooth(method="loess", se=T, colour="black", na.rm=TRUE)+
  ylab("Hardness (instant modulus test MPa)")+
  xlab("Sugar content (% Brix)")+
  labs(title=expression(paste(bold("Currants"))))+
  theme_classic(base_size=10)+ 
  scale_color_manual(values=c("#339900","#FFFF00","#FF9933","#663300"), name="Ripeness stage",
                     labels=c("Day1-4", "Day5-8", "Day9-10", "Day13-16"))+ 
  theme(legend.position="none")+
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black")) 

  
  
### 3) Suppl. Table S3: Spearman correlation tests ----
# by Day ----
# colour: blue-yellow, primate visual system
cor.test(as.numeric(ba.dat.good$Day), ba.dat.good$Primate.blue.yellow.contrast, method = "spearman")
cor.test(as.numeric(ma.dat$Day), ma.dat$Primate.blue.yellow.contrast, method = "spearman")
cor.test(as.numeric(cu.dat$Day), cu.dat$Primate.blue.yellow.contrast, method = "spearman")
# colour: blue-yellow, passerine visual system
cor.test(as.numeric(ba.dat.good$Day), ba.dat.good$Passerine.blue.yellow.contrast, method = "spearman")
cor.test(as.numeric(ma.dat$Day), ma.dat$Passerine.blue.yellow.contrast, method = "spearman")
cor.test(as.numeric(cu.dat$Day), cu.dat$Passerine.blue.yellow.contrast, method = "spearman")

# colour: red-green, primate visual system 
cor.test(as.numeric(ba.dat.good$Day), ba.dat.good$Primate.red.green.contrast, method = "spearman")
cor.test(as.numeric(ma.dat$Day), ma.dat$Primate.red.green.contrast, method = "spearman")
cor.test(as.numeric(cu.dat$Day), cu.dat$Primate.red.green.contrast, method = "spearman")
# colour: red-green, passerine visual system
cor.test(as.numeric(ba.dat.good$Day), ba.dat.good$Passerine.red.green.contrast, method = "spearman")
cor.test(as.numeric(ma.dat$Day), ma.dat$Passerine.red.green.contrast, method = "spearman")
cor.test(as.numeric(cu.dat$Day), cu.dat$Passerine.red.green.contrast, method = "spearman")

# colour: luminance, primate visual system 
cor.test(as.numeric(ba.dat.good$Day), ba.dat.good$Primate.luminance, method = "spearman")
cor.test(as.numeric(ma.dat$Day), ma.dat$Primate.luminance, method = "spearman")
cor.test(as.numeric(cu.dat$Day), cu.dat$Primate.luminance, method = "spearman")
# colour: luminance, passerine visual system 
cor.test(as.numeric(ba.dat.good$Day), ba.dat.good$Passerine.luminance, method = "spearman")
cor.test(as.numeric(ma.dat$Day), ma.dat$Passerine.luminance, method = "spearman")
cor.test(as.numeric(cu.dat$Day), cu.dat$Passerine.luminance, method = "spearman")

# odour complexity: Shannon index
cor.test(as.numeric(ba.dat.good$Day), ba.dat.good$Shannon.index.std, method = "spearman")
cor.test(as.numeric(ma.dat$Day), ma.dat$Shannon.index.std, method = "spearman")
cor.test(as.numeric(cu.dat$Day), cu.dat$Shannon.index.std, method = "spearman")

# hardness: instant modulus test
cor.test(as.numeric(ba.dat.good$Day), ba.dat.good$Instant.modulus.test, method = "spearman")
cor.test(as.numeric(ma.dat$Day), ma.dat$Instant.modulus.test, method = "spearman")
cor.test(as.numeric(cu.dat$Day), cu.dat$Instant.modulus.test, method = "spearman")


# by Sugar content ----
# colour: blue-yellow, primate visual system
cor.test(ba.dat.good$Sugar.content, ba.dat.good$Primate.blue.yellow.contrast, method = "spearman")
cor.test(ma.dat$Sugar.content, ma.dat$Primate.blue.yellow.contrast, method = "spearman")
cor.test(cu.dat$Sugar.content, cu.dat$Primate.blue.yellow.contrast, method = "spearman")
# colour: blue-yellow, passerine visual system
cor.test(ba.dat.good$Sugar.content, ba.dat.good$Passerine.blue.yellow.contrast, method = "spearman")
cor.test(ma.dat$Sugar.content, ma.dat$Passerine.blue.yellow.contrast, method = "spearman")
cor.test(cu.dat$Sugar.content, cu.dat$Passerine.blue.yellow.contrast, method = "spearman")

# colour: red-green, primate visual system 
cor.test(ba.dat.good$Sugar.content, ba.dat.good$Primate.red.green.contrast, method = "spearman")
cor.test(ma.dat$Sugar.content, ma.dat$Primate.red.green.contrast, method = "spearman")
cor.test(cu.dat$Sugar.content, cu.dat$Primate.red.green.contrast, method = "spearman")
# colour: red-green, passerine visual system
cor.test(ba.dat.good$Sugar.content, ba.dat.good$Passerine.red.green.contrast, method = "spearman")
cor.test(ma.dat$Sugar.content, ma.dat$Passerine.red.green.contrast, method = "spearman")
cor.test(cu.dat$Sugar.content, cu.dat$Passerine.red.green.contrast, method = "spearman")

# colour: luminance, primate visual system 
cor.test(ba.dat.good$Sugar.content, ba.dat.good$Primate.luminance, method = "spearman")
cor.test(ma.dat$Sugar.content, ma.dat$Primate.luminance, method = "spearman")
cor.test(cu.dat$Sugar.content, cu.dat$Primate.luminance, method = "spearman")
# colour: luminance, passerine visual system 
cor.test(ba.dat.good$Sugar.content, ba.dat.good$Passerine.luminance, method = "spearman")
cor.test(ma.dat$Sugar.content, ma.dat$Passerine.luminance, method = "spearman")
cor.test(cu.dat$Sugar.content, cu.dat$Passerine.luminance, method = "spearman")

# odour complexity: Shannon index
cor.test(ba.dat.good$Sugar.content, ba.dat.good$Shannon.index.std, method = "spearman")
cor.test(ma.dat$Sugar.content, ma.dat$Shannon.index.std, method = "spearman")
cor.test(cu.dat$Sugar.content, cu.dat$Shannon.index.std, method = "spearman")

# hardness: instant modulus test
cor.test(ba.dat.good$Sugar.content, ba.dat.good$Instant.modulus.test, method = "spearman")
cor.test(ma.dat$Sugar.content, ma.dat$Instant.modulus.test, method = "spearman")
cor.test(cu.dat$Sugar.content, cu.dat$Instant.modulus.test, method = "spearman")



### 4) LMMs ----
# Table 1: Primate visual system ----
ba.mod.pri<- lmer(Sugar.content~Primate.blue.yellow.contrast+Primate.red.green.contrast+Primate.luminance+Shannon.index.std+Instant.modulus.test+(1|Fruit.batch), data=ba.dat.good)
summary(ba.mod.pri)

ma.mod.pri<- lmer(Sugar.content~Primate.blue.yellow.contrast+Primate.red.green.contrast+Primate.luminance+Shannon.index.std+Instant.modulus.test+(1|Fruit.batch), data=ma.dat)
summary(ma.mod.pri)

cu.mod.pri<- lmer(Sugar.content~Primate.blue.yellow.contrast+Primate.red.green.contrast+Primate.luminance+Shannon.index.std+Instant.modulus.test+(1|Fruit.batch), data=cu.dat)
summary(cu.mod.pri)


# Table 2: Passerine visual system ----
ba.mod.pas<- lmer(Sugar.content~Passerine.blue.yellow.contrast+Passerine.red.green.contrast+Passerine.luminance+Shannon.index.std+Instant.modulus.test+(1|Fruit.batch), data=ba.dat.good)
summary(ba.mod.pas)

ma.mod.pas<- lmer(Sugar.content~Passerine.blue.yellow.contrast+Passerine.red.green.contrast+Passerine.luminance+Shannon.index.std+Instant.modulus.test+(1|Fruit.batch), data=ma.dat)
summary(ma.mod.pas)

cu.mod.pas<- lmer(Sugar.content~Passerine.blue.yellow.contrast+Passerine.red.green.contrast+Passerine.luminance+Shannon.index.std+Instant.modulus.test+(1|Fruit.batch), data=cu.dat)
summary(cu.mod.pas)

# Checking model residuals (package DHARMa) ----
plot(simulateResiduals(fittedModel=ba.mod.pri, n=1000, refit=F, use.u=F)) #ok
plot(simulateResiduals(fittedModel=ma.mod.pri, n=1000, refit=F, use.u=F)) #ok
plot(simulateResiduals(fittedModel=cu.mod.pri, n=1000, refit=F, use.u=F)) #ok

plot(simulateResiduals(fittedModel=ba.mod.pas, n=1000, refit=F, use.u=F)) #ok
plot(simulateResiduals(fittedModel=ma.mod.pas, n=1000, refit=F, use.u=F)) #ok
plot(simulateResiduals(fittedModel=cu.mod.pas, n=1000, refit=F, use.u=F)) #ok



### 5) Pairwise correlations across all fruit traits ----
# Pairwise Spearman's ⍴ correlation tests ----
# Primate visual system
ba.part.pri = subset(ba.dat.good, select=c("Sugar.content", "Primate.blue.yellow.contrast" , "Primate.red.green.contrast", "Primate.luminance","Shannon.index.std","Instant.modulus.test"))
cor(ba.part.pri, method = "spearman")
ba.part.pri[is.na(ba.part.pri)] = 0
cor.ba.pri<-cor(ba.part.pri, method = "spearman")
cor.ba.pri

ma.part.pri = subset(ma.dat, select=c("Sugar.content", "Primate.blue.yellow.contrast" , "Primate.red.green.contrast", "Primate.luminance","Shannon.index.std","Instant.modulus.test"))
cor(ma.part.pri, method = "spearman")
ma.part.pri[is.na(ma.part.pri)] = 0
cor.ma.pri<-cor(ma.part.pri, method = "spearman")
cor.ma.pri

cu.part.pri = subset(cu.dat, select=c("Sugar.content", "Primate.blue.yellow.contrast" , "Primate.red.green.contrast", "Primate.luminance","Shannon.index.std","Instant.modulus.test"))
cor(cu.part.pri, method = "spearman")
cu.part.pri[is.na(cu.part.pri)] = 0
cor.cu.pri<-cor(cu.part.pri, method = "spearman")
cor.cu.pri


# Passerine visual system
ba.part.pas = subset(ba.dat.good, select=c("Sugar.content", "Passerine.blue.yellow.contrast" , "Passerine.red.green.contrast", "Passerine.luminance","Shannon.index.std","Instant.modulus.test"))
cor(ba.part.pas, method = "spearman")
ba.part.pas[is.na(ba.part.pas)] = 0
cor.ba.pas<-cor(ba.part.pas, method = "spearman")
cor.ba.pas

ma.part.pas = subset(ma.dat, select=c("Sugar.content", "Passerine.blue.yellow.contrast" , "Passerine.red.green.contrast", "Passerine.luminance","Shannon.index.std","Instant.modulus.test"))
cor(ma.part.pas, method = "spearman")
ma.part.pas[is.na(ma.part.pas)] = 0
cor.ma.pas<-cor(ma.part.pas, method = "spearman")
cor.ma.pas

cu.part.pas = subset(cu.dat, select=c("Sugar.content", "Passerine.blue.yellow.contrast" , "Passerine.red.green.contrast", "Passerine.luminance","Shannon.index.std","Instant.modulus.test"))
cor(cu.part.pas, method = "spearman")
cu.part.pas[is.na(cu.part.pas)] = 0
cor.cu.pas<-cor(cu.part.pas, method = "spearman")
cor.cu.pas


# Figure 3 ----
# Primate visual system
#cor.ba.pri.qg <-
qgraph(cor(cor.ba.pri), shape= "circle" , posCol= "darkorange" , negCol= "darkblue" , layout= "circle" , vsize=10, legend= F, labels = c("S", "H", "BY" ,"RG", "L", "O"), title= "A) Bananas (Primate visual system)", title.cex=1.2, label.cex=1.5)

#cor.ma.pri.qg <-
qgraph(cor(cor.ma.pri), shape= "circle" , posCol= "darkorange" , negCol= "darkblue" , layout= "circle" , vsize=10, legend= F, labels = c("S", "H", "BY" ,"RG", "L", "O"), title= "B) Mangoes (Primate visual system)", title.cex=1.2, label.cex=1.5)

#cor.cu.pri.qg <-
qgraph(cor(cor.cu.pri), shape= "circle" , posCol= "darkorange" , negCol= "darkblue" , layout= "circle" , vsize=10, legend= F, labels = c("S", "H", "BY" ,"RG", "L", "O"), title= "C) Currants (Primate visual system)", title.cex=1.2, label.cex=1.5)


# Passerine visual system
#cor.ba.pas.qg <-
qgraph(cor(cor.ba.pas), shape= "circle" , posCol= "darkorange" , negCol= "darkblue" , layout= "circle" , vsize=10, legend= F, labels = c("S", "H", "BY" ,"RG", "L", "O"), title= "D) Bananas (Passerine visual system)", title.cex=1.2, label.cex=1.5)

#cor.ma.pas.qg <-
qgraph(cor(cor.ma.pas), shape= "circle" , posCol= "darkorange" , negCol= "darkblue" , layout= "circle" , vsize=10, legend= F, labels = c("S", "H", "BY" ,"RG", "L", "O"), title= "E) Mangoes (Passerine visual system)", title.cex=1.2, label.cex=1.5)

#cor.cu.pas.qg <-
qgraph(cor(cor.cu.pas), shape= "circle" , posCol= "darkorange" , negCol= "darkblue" , layout= "circle" , vsize=10, legend= F, labels = c("S", "H", "BY" ,"RG", "L", "O"), title= "F) Currants (Passerine visual system)", title.cex=1.2, label.cex=1.5)


