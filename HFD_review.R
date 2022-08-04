### HFD review edits 

#set up
library(ggplot2)
library(phyloseq)
setwd("~/Desktop/School/Lab/WetLab/Microbiome/HFD/Final")
ps <- readRDS("~/Desktop/School/Lab/WetLab/Microbiome/HFD/Final/ps.pruned.rds")
my_theme<-    theme(
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
  legend.key.width = unit(1, "cm")
) 




#relevel group variable
sample_data(ps)$Group <- factor(sample_data(ps)$Group, 
                                levels = c("LFD Sham", "HFD Sham", "LFD TBI", "HFD TBI"))



#analyze class-level data within specific phyla
filt <- genefilter_sample(ps, filterfun_sample(function(x) x > 3))
ps1 <- prune_taxa(filt, ps)
ps1 <- transform_sample_counts(ps1, function(x) x/sum(x))
class <- tax_glom(ps1, taxrank="Class") 
class <- prune_taxa(taxa_sums(class)>=0.01, class)
otu <- as.data.frame(otu_table(class))
TAX1 <- as(tax_table(class), "matrix")
tax <- as.data.frame(TAX1)
tax <- data.frame(lapply(tax, as.character), stringsAsFactors=FALSE)
colnames(otu) <- tax$Class
otu$Other <- (1-rowSums(otu))
tax[nrow(tax)+1,"Class"] <- "Other"
rownames(tax) <- tax$Class
class@otu_table <- otu_table(as.matrix(otu), taxa_are_rows=FALSE)
class@tax_table <- tax_table(as.matrix(tax))
class_df <- psmelt(class)

#select phylum to analyze
classes<-subset(class_df, Phylum=="Firmicutes")
classes<-subset(class_df, Phylum=="Bacteroidetes")


#prep data for plotting
classes$Time <- as.factor(classes$Time)
classes$Abundance <- as.numeric(classes$Abundance)
classes$Group <- as.factor(classes$Group)
classes$Phylum <- as.factor(classes$Phylum)
classes$Time <- relevel(classes$Time, "Pre")
classes[is.na(classes)] <- 0
library(plotrix)
StdEr<-aggregate(Abundance ~ ID + Diet + Injury + Time + Class, classes, mean)
StdEr<-aggregate(Abundance ~ Diet + Injury + Time + Class, StdEr, std.error)
colnames(StdEr)[colnames(StdEr)=="Abundance"] <- "Class.Error"

library(tidyr)
class_avg <- aggregate(classes[,"Abundance"],
                     by = list(classes[,"Time"], 
                               classes[,"Injury"], 
                               classes[,"Diet"],
                               classes[,"Class"]),
                     FUN = function(x) mean = mean(x))
class_avg <- do.call(data.frame, class_avg)
colnames(class_avg) <- c("Time", "Injury", "Diet", "Class", "Abundance")
class_avg$Group <- paste(class_avg$Diet, class_avg$Injury)

class_avg<-merge(class_avg, StdEr, by=c("Diet", "Injury", "Time", "Class"), all=T)
class_avg$Class<-as.factor(class_avg$Class)
class_avg$Class<-factor(class_avg$Class, levels = c("Bacilli", "Erysipelotrichia", "Clostridia"))

ggplot(data=class_avg, aes(Time, Abundance, color=Group, group=Group, shape=Group))+ 
  geom_line(size=2.5)+
  geom_point(aes(fill=Group), color="black", size=7)+
  facet_wrap(~Class)+
  ggtitle("Class Abundance Within Firmicutes")+   
  scale_fill_manual(values = c("#8BC589", "#2DA05A", "#32648C", "#234664"))+
  scale_color_manual(values = c("#8BC589", "#2DA05A", "#32648C", "#234664"))+
  scale_shape_manual(values=c(24, 21, 24, 21))+
  scale_x_discrete(breaks=c("Pre","D3","D30", "D60"),labels=c("Base", "3", "30", "60"))+  
  xlab("Days Post-Injury")+
  ylab("Relative Abundance")+
  geom_errorbar(aes(ymin=Abundance, ymax=Abundance+Class.Error), width=.1, size=2)+
  my_theme
ggsave("class_firm.png", width = 38, height = 22, units = "cm")


#run PERMANOVA
temp <- subset(classes, Time != "Pre")
temp <- subset(temp, select = -c(OTU, Sample, SampleID, Group, Kingdom,Phylum))
wide <- reshape(temp, 
                       timevar = "Class",
                       idvar = c("ID", "Time", "Injury", "Diet"),
                       direction = "wide")
IVs <- subset(wide, select = c(ID, Time, Injury, Diet))
DVs <- subset(wide, select = -c(ID, Time, Injury, Diet))
DVs[is.na(DVs)] <- 0
library(vegan)
permanova <- adonis(DVs ~ Injury*Diet*Time, data=IVs, permutations=999)
permanova

#injury x time interaction 
temp<-subset(wide, Time=="D3")
wilcox.test(temp$Abundance.Clostridia~ temp$Injury)
wilcox.test(temp$Abundance.Bacilli~ temp$Injury)
wilcox.test(temp$Abundance.Erysipelotrichia~ temp$Injury)

#diet effect
wilcox.test(wide$Abundance.Clostridia~ wide$Diet)
wilcox.test(wide$Abundance.Bacilli~ wide$Diet)
wilcox.test(wide$Abundance.Erysipelotrichia~ wide$Diet)

#time effects
library(FSA)
dunnTest(Abundance.Clostridia ~ Time,
               data=wide,
               method="bh")
dunnTest(Abundance.Bacilli ~ Time,
         data=wide,
         method="bh")
dunnTest(Abundance.Erysipelotrichia ~ Time,
         data=wide,
         method="bh")







#analyze order-level data 
filt <- genefilter_sample(ps, filterfun_sample(function(x) x > 3))
ps1 <- prune_taxa(filt, ps)
ps1 <- transform_sample_counts(ps1, function(x) x/sum(x))
order <- tax_glom(ps1, taxrank="Order") 
order <- prune_taxa(taxa_sums(order)>=0.01, order)
otu <- as.data.frame(otu_table(order))
TAX1 <- as(tax_table(order), "matrix")
tax <- as.data.frame(TAX1)
tax <- data.frame(lapply(tax, as.character), stringsAsFactors=FALSE)
colnames(otu) <- tax$Order
otu$Other <- (1-rowSums(otu))
tax[nrow(tax)+1,"Order"] <- "Other"
rownames(tax) <- tax$Order
order@otu_table <- otu_table(as.matrix(otu), taxa_are_rows=FALSE)
order@tax_table <- tax_table(as.matrix(tax))
order_df <- psmelt(order)

#select class to analyze
order<-subset(order_df, Class=="Clostridia")
order<-subset(order_df, Class=="Bacteroidia")

#prep data for plotting
order$Time <- as.factor(order$Time)
order$Abundance <- as.numeric(order$Abundance)
order$Group <- as.factor(order$Group)
order$Phylum <- as.factor(order$Phylum)
order$Time <- relevel(order$Time, "Pre")
order[is.na(order)] <- 0
library(plotrix)
StdEr<-aggregate(Abundance ~ ID + Diet + Injury + Time + Order, order, mean)
StdEr<-aggregate(Abundance ~ Diet + Injury + Time + Order, StdEr, std.error)
colnames(StdEr)[colnames(StdEr)=="Abundance"] <- "Order.Error"

library(tidyr)
order_avg <- aggregate(order[,"Abundance"],
                       by = list(order[,"Time"], 
                                 order[,"Injury"], 
                                 order[,"Diet"],
                                 order[,"Order"]),
                       FUN = function(x) mean = mean(x))
order_avg <- do.call(data.frame, order_avg)
colnames(order_avg) <- c("Time", "Injury", "Diet", "Order", "Abundance")
order_avg$Group <- paste(order_avg$Diet, order_avg$Injury)
order_avg<-merge(order_avg, StdEr, by=c("Diet", "Injury", "Time", "Order"), all=T)

#plot order data
ggplot(data=order_avg, aes(Time, Abundance, color=Group, group=Group, shape=Group))+ 
  geom_line(size=2.5)+
  geom_point(aes(fill=Group), color="black", size=7)+
  facet_wrap(~Order)+
  ggtitle("Order Abundance Within Bacteroidia")+   
  scale_fill_manual(values = c("#8BC589", "#2DA05A", "#32648C", "#234664"))+
  scale_color_manual(values = c("#8BC589", "#2DA05A", "#32648C", "#234664"))+
  scale_shape_manual(values=c(24, 21, 24, 21))+
  scale_x_discrete(breaks=c("Pre","D3","D30", "D60"),labels=c("Base", "3", "30", "60"))+  
  xlab("Days Post-Injury")+
  ylab("Relative Abundance")+
  geom_errorbar(aes(ymin=Abundance, ymax=Abundance+Order.Error), width=.1, size=2)+
  my_theme+
  theme(strip.text.x = element_text(size=30))
ggsave("order_bac.png", width = 30, height = 22, units = "cm")


#permanova despite having 1 DV (friedman doesn't work for multiple IVs)
temp <- subset(order, Time != "Pre")
temp<-temp[,c(3, 5:8,13)]
wide <- reshape(temp, 
                timevar = "Order",
                idvar = c("ID", "Time", "Injury", "Diet"),
                direction = "wide")
IVs <- subset(wide, select = c(ID, Time, Injury, Diet))
DVs <- subset(wide, select = -c(ID, Time, Injury, Diet))
DVs[is.na(DVs)] <- 0
permanova <- adonis(DVs ~ Injury*Diet*Time, data=IVs, permutations=999)
permanova

#post hoc for time
dunnTest(Abundance.Clostridiales ~ Time,
         data=wide,
         method="bh")
temp<-subset(wide, Time=="D60")
wilcox.test(temp$Abundance.Clostridiales~ temp$Injury)









#additional picrust analysis 
PiCrustPWsD60 <- read_csv("~/Desktop/School/Lab/WetLab/Microbiome/PiCrust2/R/Aldex2 Analyses/PiCrustPWsD60.csv")
aldex2_PW1<-PiCrustPWsD60 
tyr<-subset(aldex2_PW1, pathway=="PWY-6630"|pathway=="TYRFUMCAT-PWY"| pathway=="PWY-6628")
tryp<-subset(aldex2_PW1, pathway=="PWY-6629"|pathway=="TRPSYN-PWY")
chor<-subset(aldex2_PW1, pathway=="ALL-CHORISMATE-PWY"|pathway=="PWY-6163"|pathway=="ARO-PWY"|pathway=="PWY-6165")

#plot
ggplot(data=aldex2_PW1)+
  geom_point(aes(x=diff.btw,y=we.ep),color="black",size=6, alpha=0.8)+
  geom_point(data=tyr,aes(x=diff.btw,y=we.ep), color="blue",size=12)+
  geom_point(data=tryp,aes(x=diff.btw,y=we.ep), color="red", size=12)+
  geom_point(data=chor,aes(x=diff.btw,y=we.ep),color="purple", size=12)+
  xlab("Median Difference")+
  ylab("p-value")+
  ggtitle("C. D60 Post-Injury Pathways")+ 
  coord_cartesian(ylim=c(0.00001,1), xlim=c(-6,4))+
  geom_hline(yintercept=0.05,size=5, color="grey",linetype="dashed")+scale_y_continuous(trans="log10")+scale_colour_manual(values=c("#000000","#0000FF", "#ff0000","#8F00FF"),labels = c("All Paths", "Tyrosine Synthesis", "Tryptophan Synthesis", "Chorismate Synthesis"))+
  theme(
    plot.title = element_text(size=35, face="bold"),
    axis.title.x = element_text(size=30, face="bold"),
    axis.title.y = element_text(size=30, face="bold"),
    axis.text.y = element_text(size=26, face="bold", color="black"),
    axis.text.x = element_text(size=26, angle=0, hjust = 0.4, face="bold", color="black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face="bold"),
    legend.key=element_blank(),
    strip.text.x = element_text(size = 30, face="bold"),  #strip is the facet
    strip.background = element_rect(color="white", fill="white"), #removes background color from facet
    panel.background = element_rect(fill="white", colour="white"), #makes plot background white
    panel.border = element_rect(colour = "black", fill=NA, size=2), #adds black border around graph
    legend.position="bottom"
  ) + guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("D60_PWs.png", width = 38, height = 22, units = "cm")





#model comparison 
Behavior_Analysis <- read_csv("Behavior.Analysis.csv")
BA<-Behavior_Analysis
BA$Injury<-as.factor(BA$Injury)
BA$Diet<-as.factor(BA$Diet)
BA$Time<-as.factor(BA$Time)
relevel(BA$Diet, ref="LFD")
BA$Diet<-relevel(BA$Diet, ref="LFD")
BA$Week <- as.numeric(BA$Week) 
BA <- subset(BA, Week>0)

library(lme4) 


### Pre models ###

#choice~alpha diversity (Pre) 
m1<-lmer(scale(t.PctChoice) ~ 
           Injury*scale(c.Week) 	
         + scale(t.ChoiceBase) + (1|ID), data=subset(BA, Time=="Pre"))
m1Boot <- bootMer(m1, AIC, nsim = 100, re.form = NA)
confint(m1Boot) 
median(m1Boot$t)




#Prematures ~ Alpha Diversity 
vals=data.frame(Time=double(), AIC1=double(), AIC2=double(), AIC3=double(), AIC4=double(),
                BIC1=double(), BIC2=double(), BIC3=double(), BIC4=double(),
                A1_lower=double(), A1_upper=double(), B1_lower=double(), B1_upper=double(),
                A2_lower=double(), A2_upper=double(),B2_lower=double(), B2_upper=double(),
                A3_lower=double(), A3_upper=double(),B3_lower=double(), B3_upper=double(),
                A4_lower=double(), A4_upper=double(),B4_lower=double(), B4_upper=double())
for(x in levels(BA$Time)){
  m2<-lmer(scale(t.PctPre) ~ 
             scale(alpha_diversity)*scale(c.Week) 
           + scale(t.PreBase) + (1|ID), data=subset(BA, Time==x))	
  m1<-lmer(scale(t.PctPre) ~ 
             Injury*scale(c.Week) 	
           + scale(t.PreBase) + (1|ID), data=subset(BA, Time==x))
  m3<-lmer(scale(t.PctPre) ~ 
             scale(alpha_diversity)*scale(c.Week) 
           + Injury*scale(c.Week) 
           + scale(t.PreBase) + (1|ID), data=subset(BA, Time==x))
  m4<-lmer(scale(t.PctPre) ~ 
             scale(alpha_diversity)*scale(c.Week)*Diet*Injury
           + scale(t.PreBase) + (1|ID), data=subset(BA, Time==x))
  m2Boot <- bootMer(m2, AIC, nsim = 100, re.form = NA);m2BootB <- bootMer(m2, BIC, nsim = 100, re.form = NA)
  m1Boot <- bootMer(m1, AIC, nsim = 100, re.form = NA);m1BootB <- bootMer(m1, BIC, nsim = 100, re.form = NA)
  m3Boot <- bootMer(m3, AIC, nsim = 100, re.form = NA);m3BootB <- bootMer(m3, BIC, nsim = 100, re.form = NA)
  m4Boot <- bootMer(m4, AIC, nsim = 100, re.form = NA);m4BootB <- bootMer(m4, BIC, nsim = 100, re.form = NA)
  AIC1<-median(m1Boot$t);BIC1<-median(m1BootB$t)
  AIC2<-median(m2Boot$t);BIC2<-median(m2BootB$t)
  AIC3<-median(m3Boot$t);BIC3<-median(m3BootB$t)
  AIC4<-median(m4Boot$t);BIC4<-median(m4BootB$t)
  temp1_AIC<-as.data.frame(confint(m1Boot)); temp1_BIC<-as.data.frame(confint(m1BootB))
  colnames(temp1_AIC)<-c("low", "up");colnames(temp1_BIC)<-c("low", "up")
  temp2_AIC<-as.data.frame(confint(m2Boot)); temp2_BIC<-as.data.frame(temp2_BIC<-confint(m2BootB))
  colnames(temp2_AIC)<-c("low", "up");colnames(temp2_BIC)<-c("low", "up")
  temp3_AIC<-as.data.frame(confint(m3Boot)); temp3_BIC<-as.data.frame(confint(m3BootB))
  colnames(temp3_AIC)<-c("low", "up");colnames(temp3_BIC)<-c("low", "up")
  temp4_AIC<-as.data.frame(confint(m4Boot)); temp4_BIC<-as.data.frame(confint(m4BootB))
  colnames(temp4_AIC)<-c("low", "up");colnames(temp4_BIC)<-c("low", "up")
  vals[nrow(vals)+1,]<-c(x, AIC1, AIC2, AIC3, AIC4, BIC1, BIC2, BIC3, BIC4,
                         temp1_AIC$low,temp1_AIC$up, temp1_BIC$low, temp1_BIC$up,
                         temp2_AIC$low,temp2_AIC$up, temp2_BIC$low, temp2_BIC$up, 
                         temp3_AIC$low,temp3_AIC$up, temp3_BIC$low, temp3_BIC$up, 
                         temp4_AIC$low,temp4_AIC$up, temp4_BIC$low, temp4_BIC$up)
}
#combine CIs into single variable 
vals$Time<-as.factor(vals$Time)
vals<-vals %>% dplyr::mutate_if(is.character,as.numeric)
vals<-vals %>% dplyr::mutate_if(is.numeric,round,digits=3)
vals$AIC1_conf<-paste(vals$A1_lower, vals$A1_upper, sep=",")
vals$AIC2_conf<-paste(vals$A2_lower, vals$A2_upper, sep=",")
vals$AIC3_conf<-paste(vals$A3_lower, vals$A3_upper, sep=",")
vals$AIC4_conf<-paste(vals$A4_lower, vals$A4_upper, sep=",")
vals$BIC1_conf<-paste(vals$B1_lower, vals$B1_upper, sep=",")
vals$BIC2_conf<-paste(vals$B2_lower, vals$B2_upper, sep=",")
vals$BIC3_conf<-paste(vals$B3_lower, vals$B3_upper, sep=",")
vals$BIC4_conf<-paste(vals$B4_lower, vals$B4_upper, sep=",")
vals<-subset(vals,select=-c(A1_lower, A1_upper, A2_lower, A2_upper,
                            A3_lower, A3_upper, A4_lower, A4_upper,
                            B1_lower, B1_upper, B2_lower, B2_upper,
                            B3_lower, B3_upper, B4_lower, B4_upper))

write.table(vals, "Pre.alpha.boot.csv", append=T, sep=",", col.names=NA)





#Prematures ~ Beta Diversity 
vals=data.frame(Time=double(), AIC1=double(), AIC2=double(), AIC3=double(), AIC4=double(),
                BIC1=double(), BIC2=double(), BIC3=double(), BIC4=double(),
                A1_lower=double(), A1_upper=double(), B1_lower=double(), B1_upper=double(),
                A2_lower=double(), A2_upper=double(),B2_lower=double(), B2_upper=double(),
                A3_lower=double(), A3_upper=double(),B3_lower=double(), B3_upper=double(),
                A4_lower=double(), A4_upper=double(),B4_lower=double(), B4_upper=double())
for(x in levels(BA$Time)){
  m2<-lmer(scale(t.PctPre) ~ 
             scale(beta_diversity)*scale(c.Week) 
           + scale(t.PreBase) + (1|ID), data=subset(BA, Time==x))	
  m1<-lmer(scale(t.PctPre) ~ 
             Injury*scale(c.Week) 	
           + scale(t.PreBase) + (1|ID), data=subset(BA, Time==x))
  m3<-lmer(scale(t.PctPre) ~ 
             scale(beta_diversity)*scale(c.Week) 
           + Injury*scale(c.Week) 
           + scale(t.PreBase) + (1|ID), data=subset(BA, Time==x))
  m4<-lmer(scale(t.PctPre) ~ 
             scale(beta_diversity)*scale(c.Week)*Diet*Injury
           + scale(t.PreBase) + (1|ID), data=subset(BA, Time==x))
  m2Boot <- bootMer(m2, AIC, nsim = 100, re.form = NA);m2BootB <- bootMer(m2, BIC, nsim = 100, re.form = NA)
  m1Boot <- bootMer(m1, AIC, nsim = 100, re.form = NA);m1BootB <- bootMer(m1, BIC, nsim = 100, re.form = NA)
  m3Boot <- bootMer(m3, AIC, nsim = 100, re.form = NA);m3BootB <- bootMer(m3, BIC, nsim = 100, re.form = NA)
  m4Boot <- bootMer(m4, AIC, nsim = 100, re.form = NA);m4BootB <- bootMer(m4, BIC, nsim = 100, re.form = NA)
  AIC1<-median(m1Boot$t);BIC1<-median(m1BootB$t)
  AIC2<-median(m2Boot$t);BIC2<-median(m2BootB$t)
  AIC3<-median(m3Boot$t);BIC3<-median(m3BootB$t)
  AIC4<-median(m4Boot$t);BIC4<-median(m4BootB$t)
  temp1_AIC<-as.data.frame(confint(m1Boot)); temp1_BIC<-as.data.frame(confint(m1BootB))
  colnames(temp1_AIC)<-c("low", "up");colnames(temp1_BIC)<-c("low", "up")
  temp2_AIC<-as.data.frame(confint(m2Boot)); temp2_BIC<-as.data.frame(temp2_BIC<-confint(m2BootB))
  colnames(temp2_AIC)<-c("low", "up");colnames(temp2_BIC)<-c("low", "up")
  temp3_AIC<-as.data.frame(confint(m3Boot)); temp3_BIC<-as.data.frame(confint(m3BootB))
  colnames(temp3_AIC)<-c("low", "up");colnames(temp3_BIC)<-c("low", "up")
  temp4_AIC<-as.data.frame(confint(m4Boot)); temp4_BIC<-as.data.frame(confint(m4BootB))
  colnames(temp4_AIC)<-c("low", "up");colnames(temp4_BIC)<-c("low", "up")
  vals[nrow(vals)+1,]<-c(x, AIC1, AIC2, AIC3, AIC4, BIC1, BIC2, BIC3, BIC4,
                         temp1_AIC$low,temp1_AIC$up, temp1_BIC$low, temp1_BIC$up,
                         temp2_AIC$low,temp2_AIC$up, temp2_BIC$low, temp2_BIC$up, 
                         temp3_AIC$low,temp3_AIC$up, temp3_BIC$low, temp3_BIC$up, 
                         temp4_AIC$low,temp4_AIC$up, temp4_BIC$low, temp4_BIC$up)
}
#combine CIs into single variable 
vals$Time<-as.factor(vals$Time)
vals<-vals %>% dplyr::mutate_if(is.character,as.numeric)
vals<-vals %>% dplyr::mutate_if(is.numeric,round,digits=3)
vals$AIC1_conf<-paste(vals$A1_lower, vals$A1_upper, sep=",")
vals$AIC2_conf<-paste(vals$A2_lower, vals$A2_upper, sep=",")
vals$AIC3_conf<-paste(vals$A3_lower, vals$A3_upper, sep=",")
vals$AIC4_conf<-paste(vals$A4_lower, vals$A4_upper, sep=",")
vals$BIC1_conf<-paste(vals$B1_lower, vals$B1_upper, sep=",")
vals$BIC2_conf<-paste(vals$B2_lower, vals$B2_upper, sep=",")
vals$BIC3_conf<-paste(vals$B3_lower, vals$B3_upper, sep=",")
vals$BIC4_conf<-paste(vals$B4_lower, vals$B4_upper, sep=",")
vals<-subset(vals,select=-c(A1_lower, A1_upper, A2_lower, A2_upper,
                            A3_lower, A3_upper, A4_lower, A4_upper,
                            B1_lower, B1_upper, B2_lower, B2_upper,
                            B3_lower, B3_upper, B4_lower, B4_upper))

write.table(vals, "Pre.beta.boot.csv", append=T, sep=",", col.names=NA)








#P2 ~ Alpha Diversity 
vals=data.frame(Time=double(), AIC1=double(), AIC2=double(), AIC3=double(), AIC4=double(),
                BIC1=double(), BIC2=double(), BIC3=double(), BIC4=double(),
                A1_lower=double(), A1_upper=double(), B1_lower=double(), B1_upper=double(),
                A2_lower=double(), A2_upper=double(),B2_lower=double(), B2_upper=double(),
                A3_lower=double(), A3_upper=double(),B3_lower=double(), B3_upper=double(),
                A4_lower=double(), A4_upper=double(),B4_lower=double(), B4_upper=double())
for(x in levels(BA$Time)){
  m2<-lmer(scale(t.PctChoice) ~ 
             scale(alpha_diversity)*scale(c.Week) 
           + scale(t.ChoiceBase) + (1|ID), data=subset(BA, Time==x))	
  m1<-lmer(scale(t.PctChoice) ~ 
             Injury*scale(c.Week) 	
           + scale(t.ChoiceBase) + (1|ID), data=subset(BA, Time==x))
  m3<-lmer(scale(t.PctChoice) ~ 
             scale(alpha_diversity)*scale(c.Week) 
           + Injury*scale(c.Week) 
           + scale(t.ChoiceBase) + (1|ID), data=subset(BA, Time==x))
  m4<-lmer(scale(t.PctChoice) ~ 
             scale(alpha_diversity)*scale(c.Week)*Diet*Injury
           + scale(t.ChoiceBase) + (1|ID), data=subset(BA, Time==x))
  m2Boot <- bootMer(m2, AIC, nsim = 100, re.form = NA);m2BootB <- bootMer(m2, BIC, nsim = 100, re.form = NA)
  m1Boot <- bootMer(m1, AIC, nsim = 100, re.form = NA);m1BootB <- bootMer(m1, BIC, nsim = 100, re.form = NA)
  m3Boot <- bootMer(m3, AIC, nsim = 100, re.form = NA);m3BootB <- bootMer(m3, BIC, nsim = 100, re.form = NA)
  m4Boot <- bootMer(m4, AIC, nsim = 100, re.form = NA);m4BootB <- bootMer(m4, BIC, nsim = 100, re.form = NA)
  AIC1<-median(m1Boot$t);BIC1<-median(m1BootB$t)
  AIC2<-median(m2Boot$t);BIC2<-median(m2BootB$t)
  AIC3<-median(m3Boot$t);BIC3<-median(m3BootB$t)
  AIC4<-median(m4Boot$t);BIC4<-median(m4BootB$t)
  temp1_AIC<-as.data.frame(confint(m1Boot)); temp1_BIC<-as.data.frame(confint(m1BootB))
  colnames(temp1_AIC)<-c("low", "up");colnames(temp1_BIC)<-c("low", "up")
  temp2_AIC<-as.data.frame(confint(m2Boot)); temp2_BIC<-as.data.frame(temp2_BIC<-confint(m2BootB))
  colnames(temp2_AIC)<-c("low", "up");colnames(temp2_BIC)<-c("low", "up")
  temp3_AIC<-as.data.frame(confint(m3Boot)); temp3_BIC<-as.data.frame(confint(m3BootB))
  colnames(temp3_AIC)<-c("low", "up");colnames(temp3_BIC)<-c("low", "up")
  temp4_AIC<-as.data.frame(confint(m4Boot)); temp4_BIC<-as.data.frame(confint(m4BootB))
  colnames(temp4_AIC)<-c("low", "up");colnames(temp4_BIC)<-c("low", "up")
  vals[nrow(vals)+1,]<-c(x, AIC1, AIC2, AIC3, AIC4, BIC1, BIC2, BIC3, BIC4,
                         temp1_AIC$low,temp1_AIC$up, temp1_BIC$low, temp1_BIC$up,
                         temp2_AIC$low,temp2_AIC$up, temp2_BIC$low, temp2_BIC$up, 
                         temp3_AIC$low,temp3_AIC$up, temp3_BIC$low, temp3_BIC$up, 
                         temp4_AIC$low,temp4_AIC$up, temp4_BIC$low, temp4_BIC$up)
}
#combine CIs into single variable 
vals$Time<-as.factor(vals$Time)
vals<-vals %>% dplyr::mutate_if(is.character,as.numeric)
vals<-vals %>% dplyr::mutate_if(is.numeric,round,digits=3)
vals$AIC1_conf<-paste(vals$A1_lower, vals$A1_upper, sep=",")
vals$AIC2_conf<-paste(vals$A2_lower, vals$A2_upper, sep=",")
vals$AIC3_conf<-paste(vals$A3_lower, vals$A3_upper, sep=",")
vals$AIC4_conf<-paste(vals$A4_lower, vals$A4_upper, sep=",")
vals$BIC1_conf<-paste(vals$B1_lower, vals$B1_upper, sep=",")
vals$BIC2_conf<-paste(vals$B2_lower, vals$B2_upper, sep=",")
vals$BIC3_conf<-paste(vals$B3_lower, vals$B3_upper, sep=",")
vals$BIC4_conf<-paste(vals$B4_lower, vals$B4_upper, sep=",")
vals<-subset(vals,select=-c(A1_lower, A1_upper, A2_lower, A2_upper,
                            A3_lower, A3_upper, A4_lower, A4_upper,
                            B1_lower, B1_upper, B2_lower, B2_upper,
                            B3_lower, B3_upper, B4_lower, B4_upper))

write.table(vals, "P2.alpha.boot.csv", append=T, sep=",", col.names=NA)






#P2 ~ Beta Diversity 
vals=data.frame(Time=double(), AIC1=double(), AIC2=double(), AIC3=double(), AIC4=double(),
                BIC1=double(), BIC2=double(), BIC3=double(), BIC4=double(),
                A1_lower=double(), A1_upper=double(), B1_lower=double(), B1_upper=double(),
                A2_lower=double(), A2_upper=double(),B2_lower=double(), B2_upper=double(),
                A3_lower=double(), A3_upper=double(),B3_lower=double(), B3_upper=double(),
                A4_lower=double(), A4_upper=double(),B4_lower=double(), B4_upper=double())
for(x in levels(BA$Time)){
  m2<-lmer(scale(t.PctChoice) ~ 
             scale(beta_diversity)*scale(c.Week) 
           + scale(t.ChoiceBase) + (1|ID), data=subset(BA, Time==x))	
  m1<-lmer(scale(t.PctChoice) ~ 
             Injury*scale(c.Week) 	
           + scale(t.ChoiceBase) + (1|ID), data=subset(BA, Time==x))
  m3<-lmer(scale(t.PctChoice) ~ 
             scale(beta_diversity)*scale(c.Week) 
           + Injury*scale(c.Week) 
           + scale(t.ChoiceBase) + (1|ID), data=subset(BA, Time==x))
  m4<-lmer(scale(t.PctChoice) ~ 
             scale(beta_diversity)*scale(c.Week)*Diet*Injury
           + scale(t.ChoiceBase) + (1|ID), data=subset(BA, Time==x))
  m2Boot <- bootMer(m2, AIC, nsim = 100, re.form = NA);m2BootB <- bootMer(m2, BIC, nsim = 100, re.form = NA)
  m1Boot <- bootMer(m1, AIC, nsim = 100, re.form = NA);m1BootB <- bootMer(m1, BIC, nsim = 100, re.form = NA)
  m3Boot <- bootMer(m3, AIC, nsim = 100, re.form = NA);m3BootB <- bootMer(m3, BIC, nsim = 100, re.form = NA)
  m4Boot <- bootMer(m4, AIC, nsim = 100, re.form = NA);m4BootB <- bootMer(m4, BIC, nsim = 100, re.form = NA)
  AIC1<-median(m1Boot$t);BIC1<-median(m1BootB$t)
  AIC2<-median(m2Boot$t);BIC2<-median(m2BootB$t)
  AIC3<-median(m3Boot$t);BIC3<-median(m3BootB$t)
  AIC4<-median(m4Boot$t);BIC4<-median(m4BootB$t)
  temp1_AIC<-as.data.frame(confint(m1Boot)); temp1_BIC<-as.data.frame(confint(m1BootB))
  colnames(temp1_AIC)<-c("low", "up");colnames(temp1_BIC)<-c("low", "up")
  temp2_AIC<-as.data.frame(confint(m2Boot)); temp2_BIC<-as.data.frame(temp2_BIC<-confint(m2BootB))
  colnames(temp2_AIC)<-c("low", "up");colnames(temp2_BIC)<-c("low", "up")
  temp3_AIC<-as.data.frame(confint(m3Boot)); temp3_BIC<-as.data.frame(confint(m3BootB))
  colnames(temp3_AIC)<-c("low", "up");colnames(temp3_BIC)<-c("low", "up")
  temp4_AIC<-as.data.frame(confint(m4Boot)); temp4_BIC<-as.data.frame(confint(m4BootB))
  colnames(temp4_AIC)<-c("low", "up");colnames(temp4_BIC)<-c("low", "up")
  vals[nrow(vals)+1,]<-c(x, AIC1, AIC2, AIC3, AIC4, BIC1, BIC2, BIC3, BIC4,
                         temp1_AIC$low,temp1_AIC$up, temp1_BIC$low, temp1_BIC$up,
                         temp2_AIC$low,temp2_AIC$up, temp2_BIC$low, temp2_BIC$up, 
                         temp3_AIC$low,temp3_AIC$up, temp3_BIC$low, temp3_BIC$up, 
                         temp4_AIC$low,temp4_AIC$up, temp4_BIC$low, temp4_BIC$up)
}
#combine CIs into single variable 
vals$Time<-as.factor(vals$Time)
vals<-vals %>% dplyr::mutate_if(is.character,as.numeric)
vals<-vals %>% dplyr::mutate_if(is.numeric,round,digits=3)
vals$AIC1_conf<-paste(vals$A1_lower, vals$A1_upper, sep=",")
vals$AIC2_conf<-paste(vals$A2_lower, vals$A2_upper, sep=",")
vals$AIC3_conf<-paste(vals$A3_lower, vals$A3_upper, sep=",")
vals$AIC4_conf<-paste(vals$A4_lower, vals$A4_upper, sep=",")
vals$BIC1_conf<-paste(vals$B1_lower, vals$B1_upper, sep=",")
vals$BIC2_conf<-paste(vals$B2_lower, vals$B2_upper, sep=",")
vals$BIC3_conf<-paste(vals$B3_lower, vals$B3_upper, sep=",")
vals$BIC4_conf<-paste(vals$B4_lower, vals$B4_upper, sep=",")
vals<-subset(vals,select=-c(A1_lower, A1_upper, A2_lower, A2_upper,
                            A3_lower, A3_upper, A4_lower, A4_upper,
                            B1_lower, B1_upper, B2_lower, B2_upper,
                            B3_lower, B3_upper, B4_lower, B4_upper))

write.table(vals, "P2.beta.boot.csv", append=T, sep=",", col.names=NA)
