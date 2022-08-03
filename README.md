# Code to accompany our publication in *Behavioral Neuroscience* 

### This code can be used to analyze and visualize alpha diversity, beta diversity, and taxa abundance for 16s sequencing of the gut microbiome. These data were collected in the [Vonder Haar lab](https://github.com/vonderhaarlab), and the full publication can be found [here](https://psycnet.apa.org/record/2022-85473-001).

In brief, the goal of the project was to determine the effects of traumatic brain injury (TBI) and a high-fat diet (HFD) on the gut microbiome across several collection timepoints. We found that TBI shifted the gut microbiome, and some acute changes were associated with behavioral outcomes.

1. Begin by loading libraries and reading in the data, which can be found [here]().
```
#set up
library(phyloseq)
library(ggplot2)
setwd("path")
ps <- readRDS("path/ps.rds")
sample_data(ps)$Group <- factor(sample_data(ps)$Group, 
                                levels = c("LFD Sham", "HFD Sham", "LFD TBI", "HFD TBI"))

```
