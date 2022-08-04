# R Code to accompany our publication in *Behavioral Neuroscience* 
<sub>Frankot, M. A., O'Hearn, C. M., Blancke, A. M., Rodriguez, B., Pechacek, K. M., Gandhi, J., Hu, G., Martens, K. M., & Vonder Haar, C. (2022). Acute gut microbiome changes after traumatic brain injury are associated with chronic deficits in decision-making and impulsivity in male rats. Behavioral Neuroscience. Advance online publication. https://doi.org/10.1037/bne0000532 </sub>

### This code can be used to analyze and visualize alpha diversity, beta diversity, and taxa abundance for 16s sequencing of the gut microbiome. These data were collected in the [Vonder Haar lab](https://github.com/vonderhaarlab), and the full publication can be found [here](https://psycnet.apa.org/record/2022-85473-001).

In brief, the goal of the project was to determine the effects of traumatic brain injury (TBI) and a high-fat diet (HFD) on the gut microbiome across several collection timepoints. We found that TBI shifted the gut microbiome, and some acute changes were associated with behavioral outcomes.

1. Begin by loading libraries and reading in the data, which can be found [here](https://github.com/mfrankz/BNE_Publication/blob/main/ps.rds).
```
#set up
library(phyloseq)
library(ggplot2)
setwd("path")
ps <- readRDS("path/ps.rds")
sample_data(ps)$Group <- factor(sample_data(ps)$Group, 
                                levels = c("LFD Sham", "HFD Sham", "LFD TBI", "HFD TBI"))
```
```
#set plotting theme
my_theme<-theme(
  plot.title = element_text(size=45, face="bold"),
  axis.title.x = element_text(size=40, face="bold"),
  axis.title.y = element_text(size=40, face="bold"),
  axis.text.y = element_text(size=35, face="bold", color="black"),
  axis.text.x = element_text(size=25, angle=0, hjust = 0.4, face="bold", color="black"),
  legend.title = element_text(size = 40, face="bold"),
  legend.text = element_text(size = 35, face="bold"),
  strip.text.x = element_text(size = 25, face="bold"), 
  strip.background = element_rect(color="white", fill="white"),
  panel.background = element_rect(fill="white", colour="white"),
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.key=element_blank(),
  legend.key.height = unit(1.5, "cm"),
  legend.key.width = unit(1, "cm"))
  ```
  
2. We will now inspect alpha diversity, a measure of within-sample richness/abundance. Alpha diversity was calculated using the Shannon Index.
```
#extract alpha diversity values into a dataframe
library(plotrix)
alpha_df <- estimate_richness(ps, split = TRUE, measure = "Shannon")
alpha_df <- merge(alpha_df, sample_data(ps), by=0)
colnames(alpha_df)[colnames(alpha_df)=="Shannon"] <- "alpha_diversity"
StdEr<-aggregate(alpha_diversity ~ ID +  Group + Time, alpha_df, mean)
StdEr<-aggregate(alpha_diversity ~ Group + Time, StdEr, std.error)
colnames(StdEr)[colnames(StdEr)=="alpha_diversity"] <- "Alpha.Error"
temp<-aggregate(alpha_diversity~Group+Time, FUN=mean, data=alpha_df)
alpha_plot<-merge(temp, StdEr, by=c("Group", "Time"), all=T)
```
```
#plot alpha diversity
ggplot(data=alpha_plot, aes(x=Time, y= alpha_diversity, color=Group, group=Group, shape=Group))+
  geom_errorbar(aes(ymin=alpha_diversity, ymax=alpha_diversity+Alpha.Error), 
                width=.15, size=1.5)+
  geom_line(size=2.5)+
  geom_point(aes(fill=Group), color="black",size=8,stroke=1)+
  scale_fill_manual(values = c("#8BC589", "#2DA05A", "#32648C", "#234664"))+
  scale_color_manual(values = c("#8BC589", "#2DA05A", "#32648C", "#234664"))+
  scale_shape_manual(values=c(24, 21, 24, 21))+
  ylab("Shannon Score")+
  xlab("Days Post-Injury")+   
  ggtitle("Alpha Diversity")+
  coord_cartesian(ylim=c(3,5))+
  my_theme
ggsave("Alpha.png", width = 28, height = 20, units = "cm")
```
<img src="https://github.com/mfrankz/BNE_Publication/blob/main/Alpha.png" width="600">

```
#analyze alpha diversity using linear mixed-effects regression

#Create baseline variable 
temp<-aggregate(cbind(alpha_diversity)~ID, data=subset(alpha_df, Time == "Pre"), mean)
colnames(temp)[colnames(temp)=="alpha_diversity"] <- "Baseline"
alpha_df<-merge(alpha_df,temp, by=c("ID"), all=T)

#perform analysis
library(lme4)
library(lmerTest)
alpha_df$Diet<-as.factor(alpha_df$Diet)
alpha_df$Injury<-as.factor(alpha_df$Injury)
alpha_df$Time<-as.factor(alpha_df$Time)
alpha_df$ID<-as.factor(alpha_df$ID)
alpha_df$Diet<-relevel(alpha_df$Diet, ref="LFD")
alpha<-subset(alpha_df, alpha_df$Time != "Pre")
alpha$Time<-droplevels(alpha$Time)
LMERalpha<-lmer(scale(alpha_diversity)~Injury*Diet*Time + 
                scale(Baseline) +(1|ID), 
                data=alpha)
summary(LMERalpha) #note that there are many methods for significance testing with mixed models

#note: post-publication, we found that adding baseline alpha diversity causes overfitting. Below, is a modified analysis to account for this:
#significance via model comparison without the baseline variable
m1<-lmer(scale(alpha_diversity)~Time +(1|ID), 
                data=subset(alpha_df,alpha_df$Time != "Pre"))
m2<-lmer(scale(alpha_diversity)~Time+Diet +(1|ID), 
         data=subset(alpha_df,alpha_df$Time != "Pre"))
m3<-lmer(scale(alpha_diversity)~Time+Diet+Injury +(1|ID), 
         data=subset(alpha_df,alpha_df$Time != "Pre"))
anova(m1, m2, m3)
```

2. We will now inspect beta diversity, a measure of dissimilarity between samples. Beta diversity was calculated using Weighted Unifrac Distance and reduced using PCoA for data visualization.
```
#plot beta diversity
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.pcoa = ordinate(ps.prop, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(ps, ord.pcoa, color="Group", shape="Group")+
  facet_wrap(~Time)+
  geom_point(aes(fill=Group), color="black",size=6)+
  scale_fill_manual(values = c("#8BC589", "#2DA05A", "#32648C", "#234664"))+
  scale_color_manual(values = c("#8BC589", "#2DA05A", "#32648C", "#234664"))+
  scale_shape_manual(values=c(24, 21, 24, 21))+
  ggtitle("PCoA on Weighted-UniFrac Distance")+
  xlim(-0.52,0.52)+
  stat_ellipse(type = "norm", linetype = 2, size=1.5)+
  my_theme
ggsave("Beta.png", width = 28, height = 20, units = "cm")
```
<img src="https://github.com/mfrankz/BNE_Publication/blob/main/Beta.png" width="600">

```
#Analyze beta diversity using PERMANOVA
#isolate beta diversity scores in dataframe
library(reshape2)
beta = phyloseq::distance(ps, "wunifrac") #option for weighted UniFrac values 
beta_df = melt(as.matrix(beta))

#format data for PERMANOVA
beta_df <- reshape(beta_df, 
                   timevar = "Var2",
                   idvar = c("Var1"),
                   direction = "wide")
names(beta_df)[names(beta_df) == "Var1"] <- "SampleID"
temp<-as.data.frame(sample_data(ps))
temp$SampleID<-as.factor(temp$SampleID)
beta_df<-cbind(beta_df, temp) 
beta_df<-subset(beta_df, select = -c(1))

#create baseline variable
beta_df$Avg <- (rowSums(beta_df[,1:96]))/96
temp<-aggregate(Avg~ID, data=subset(beta_df, Time == "Pre"), mean)
colnames(temp)[colnames(temp)=="Avg"] <- "Baseline"
beta_df<-merge(beta_df,temp, by=c("ID"), all=T)

#prep data for analysis
beta_df<-subset(beta_df, Time!="Pre")
beta_df$Time<-droplevels(beta_df$Time)
IVs <- subset(beta_df, select = c("SampleID", "ID", "Time", "Injury", "Group", "Diet"))
DVs <- subset(beta_df, select = -c(SampleID, ID, Time, Injury, Group, Diet))
set.seed(50)
permanova <- vegan::adonis2(DVs ~ Injury*Diet*Time, data=IVs, permutations=999)
permanova
#note: permanova draws random permutations from the data. Here, we have set a random seed for reproducibility. The seed was not set for the manuscript analyses, so your results might look slightly different
```

3. We will now inspect the abundance of specific taxa at the phylum level. 
```
#isolate phylum level data
filt <- genefilter_sample(ps, filterfun_sample(function(x) x > 3))
ps1 <- prune_taxa(filt, ps)
ps1 <- transform_sample_counts(ps1, function(x) x/sum(x))
phy <- tax_glom(ps1, taxrank="Phylum") #note: can change taxrank to family, genus, etc.
phy <- prune_taxa(taxa_sums(phy)>=0.01, phy)
phy_df <- psmelt(phy)
```
```
#visualize TBI subject
phy_TBI<- subset(phy_df, phy_df$Injury == "TBI")
ggplot(phy_TBI, aes(x = Time, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0))+ 
  facet_wrap(~ID)+  	
  ggtitle("TBI Subjects") +
  xlab("Days Post-Injury")+
  my_theme+
  theme(axis.text.x=element_text(size=14))
ggsave("TBI_phylum.png", width = 34, height = 20, units = "cm")
```
<img src="https://github.com/mfrankz/BNE_Publication/blob/main/TBI_phylum.png" width="600">
```
#visualize sham subjects
phy_SHAM<- subset(phy_df, phy_df$Injury == "Sham")
ggplot(phy_SHAM, aes(x = Time, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0))+ 
  facet_wrap(~ID)+  	
  ggtitle("Sham Subjects") +
  xlab("Days Post-Injury")+
  my_theme+
  theme(axis.text.x=element_text(size=14))
ggsave("SHAM_phylum.png", width = 34, height = 20, units = "cm")
```
<img src="https://github.com/mfrankz/BNE_Publication/blob/main/SHAM_phylum.png" width="600">
