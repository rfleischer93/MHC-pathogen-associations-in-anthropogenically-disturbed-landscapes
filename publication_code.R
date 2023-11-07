
# R script for the manuscript "Hidden zoonotic risk: immunogenetic-pathogen networks shrink in a generalist rodent inhabiting disturbed landscapes"
# by 

# Ramona Fleischer
# ramona.fleischer@uni-ulm.de


### Load libraries 

packages <-
  c(
    "plyr", "tidyverse", "dplyr", "reshape2", "expss",
    "ggplot2", "RColorBrewer", "sjPlot", "sjmisc", "ggsignif", "ggbiplot",
    "rstatix", "car", "vegan", "ggrepel", "DT", "ggpubr",   "remotes", 
    "effects", "glmmTMB","lme4", "lmerTest", "boot", "ggeffects", "cooccur", "nlme", "brglm2", "usdm", "MuMIn") 

lapply(packages, function(y) {
  # # check if installed, if not install
  # if (!y %in% installed.packages()[, "Package"])
  #   install.packages(y)
  
  # load package
  try(require(y, character.only = T), silent = T)
})


### Load data

meta <- read.csv("prse_meta_full.csv")




################################################################################################################
############################################# Pathogens overview ###############################################
################################################################################################################


### Pathogen diversity across landscapes ###

# make landscape names short and easy
levels(meta$landscape) <- list(forest ="continuous.forest" ,fragment="forest.fragment",island="island") 

# define landscape pairs to compare
my_comparisons <- list(c("forest", "fragment"),c("forest", "island"),c("fragment", "island")) # landscapes

# number of nematodes (NNI)

landscape.nni <- meta %>%
  mutate(landscape = fct_relevel(landscape,"forest","island","fragment")) %>%
  ggplot(aes(x = landscape, y = NNI_helm, fill=landscape))+
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_point(aes(fill = landscape), size=1)+
  scale_fill_manual(values=c("chartreuse4","#54bed6", "orange2"))+  
  theme_bw(base_size = 13)+xlab("")+ylab("NNI")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  scale_x_discrete(labels=c("forest" = "C", "island" = "I", "fragment" = "A"))+
  ylim(-1,11)+
  ggsignif::geom_signif(comparisons = list(c("forest", "fragment"),c("forest", "island"),c("fragment", "island")),
                        y_position= c(7,8,9),
                        test = "wilcox.test",
                        map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05))

# number of viruses (NVI)

landscape.nvi <- meta %>%
  mutate(landscape = fct_relevel(landscape,"forest","island","fragment")) %>%
  ggplot(aes(x = landscape, y = NVI, fill=landscape))+
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_point(aes(fill = landscape), size=1)+
  scale_fill_manual(values=c("chartreuse4","#54bed6", "orange2"))+  
  theme_bw(base_size = 13)+xlab("")+ylab("NVI")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  scale_x_discrete(labels=c("forest" = "C", "island" = "I", "fragment" = "A"))+
  ylim(-1,11)+
  ggsignif::geom_signif(comparisons = list(c("forest", "fragment"),c("forest", "island"),c("fragment", "island")),
                        y_position= c(7,8,9),
                        test = "wilcox.test",
                        map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05))

# scaled nematode egg counts (FEC)

landscape.fec <- meta %>%
  mutate(landscape = fct_relevel(landscape,"forest","island","fragment")) %>%
  ggplot(aes(x = landscape, y = ZFEC_per_g, fill=landscape))+
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_point(aes(fill = landscape), size=1)+
  scale_fill_manual(values=c("chartreuse4","#54bed6", "orange2"))+  
  theme_bw(base_size = 13)+xlab("")+ylab("FEC")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  scale_x_discrete(labels=c("forest" = "C", "island" = "I", "fragment" = "A"))+
  ylim(-1,11)+
  ggsignif::geom_signif(comparisons = list(c("forest", "fragment"),c("forest", "island"),c("fragment", "island")),
                        y_position= c(7,8,9),
                        test = "wilcox.test",
                        map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05))

# Trypanosoma infection

meta_tryp <- subset(meta, Try =="1") # subset only Trypanosoma infected individuals

landscape.try <- 
  meta_tryp %>%
  mutate(landscape = fct_relevel(landscape,"forest","island","fragment")) %>%
  ggplot(aes(x=Try, fill = landscape)) + 
  geom_bar(stat="count", color = "black",  position="dodge") + 
  scale_fill_manual(values=c("chartreuse4","#54bed6", "orange2"), labels =c("C", "I", "A"))+ 
  ylab("Trypanosoma infected") + 
  xlab("                    
                       ")+
  theme_bw(base_size=13)+
  ylim(-1,11)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.position = "none")

# test
pairwise.wilcox.test(meta$Try, meta$landscape, p.adjust.method = "BH")


### Figure 1A ###

pathogen_div_plot <- ggarrange(landscape.nni, landscape.nvi, landscape.fec, landscape.try, common.legend = TRUE, legend = "none", ncol=4)
pathogen_div_plot  




### Pathogen Multivariate analyses ###

# subset pathogens
meta_pathogen <- subset(meta, select=c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11","N12","N13","Try","HdV","HpV","PbV","PcV"))

# make pca
pca_pathogen <- prcomp(meta_pathogen)

# overview
str(pca_pathogen)
print(pca_pathogen)

plot(pca_pathogen)
summary(pca_pathogen)

# simple plot
ggbiplot(pca_pathogen)


# Extract axis 1+2 into new df to plot pca in ggplot
pca_pathogen_scores <- pca_pathogen$x
pca_pathogen_scores <- data.frame(pca_pathogen_scores[,1:2])

# add field id and landscape info
pca_pathogen_scores$field_id <- meta$field_id
pca_pathogen_scores$landscape <- meta$landscape


# plot with ggplot2

pathogen_comp_mds <- ggplot(pca_pathogen_scores)+
  geom_point(aes(x=PC1,y=PC2, color=landscape), size=3) + 
  stat_ellipse(aes(x=PC1,y=PC2,color=landscape),type = "norm", size = 1) +
  theme_bw() + 
  scale_color_manual(values=c("chartreuse4", "orange2","#54bed6"))+ 
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) +
  labs(x = "Axis 1 [16.8 %]", y ="Axis 2 [13.4 %]")+
  theme(axis.text = element_text( size=12), axis.title = element_text( size=12))+
  annotate(geom="text", x=-0.7, y=1.5, label="R²= 0.07, p < 0.001", size = 4)


# calculate distance matrix for permutation test using jaccard (for presence / absence data)

pathogen_dist <- vegdist(meta_pathogen, method = "jaccard")

pathogen_adonis <- adonis2(pathogen_dist ~ landscape, data = pca_pathogen_scores) 
pathogen_adonis


################
### Figure 1 ###
################  

ggarrange(pathogen_div_plot,pathogen_comp_mds, labels = c("A", "B"), ncol=2, legend = "none", widths = c(1.5,1))




################################################################################################################
################################################ MHC overview ##################################################
################################################################################################################

# Rename alleles correctly

# Alleles

names(meta)[names(meta) == "PrseDRB.003"] <- "PrseDRB*003"
names(meta)[names(meta) == "PrseDRB.004"] <- "PrseDRB*004"
names(meta)[names(meta) == "PrseDRB.005"] <- "PrseDRB*005"
names(meta)[names(meta) == "PrseDRB.006"] <- "PrseDRB*006"
names(meta)[names(meta) == "PrseDRB.007"] <- "PrseDRB*007"
names(meta)[names(meta) == "PrseDRB.008"] <- "PrseDRB*008"
names(meta)[names(meta) == "PrseDRB.009"] <- "PrseDRB*009"
names(meta)[names(meta) == "PrseDRB.011"] <- "PrseDRB*011"
names(meta)[names(meta) == "PrseDRB.012"] <- "PrseDRB*012"
names(meta)[names(meta) == "PrseDRB.013"] <- "PrseDRB*013"
names(meta)[names(meta) == "PrseDRB.014"] <- "PrseDRB*014"
names(meta)[names(meta) == "PrseDRB.015"] <- "PrseDRB*015"
names(meta)[names(meta) == "PrseDRB.016"] <- "PrseDRB*016"
names(meta)[names(meta) == "PrseDRB.017"] <- "PrseDRB*017"
names(meta)[names(meta) == "PrseDRB.018"] <- "PrseDRB*018"
names(meta)[names(meta) == "PrseDRB.019"] <- "PrseDRB*019"
names(meta)[names(meta) == "PrseDRB.021"] <- "PrseDRB*021"
names(meta)[names(meta) == "PrseDRB.023"] <- "PrseDRB*023"
names(meta)[names(meta) == "PrseDRB.024"] <- "PrseDRB*024"
names(meta)[names(meta) == "PrseDRB.025"] <- "PrseDRB*025"
names(meta)[names(meta) == "PrseDRB.026"] <- "PrseDRB*026"
names(meta)[names(meta) == "PrseDRB.027"] <- "PrseDRB*027"
names(meta)[names(meta) == "PrseDRB.028"] <- "PrseDRB*028"
names(meta)[names(meta) == "PrseDRB.029"] <- "PrseDRB*029"
names(meta)[names(meta) == "PrseDRB.030"] <- "PrseDRB*030"
names(meta)[names(meta) == "PrseDRB.031"] <- "PrseDRB*031"
names(meta)[names(meta) == "PrseDRB.032"] <- "PrseDRB*032"
names(meta)[names(meta) == "PrseDRB.033"] <- "PrseDRB*033"
names(meta)[names(meta) == "PrseDRB.036"] <- "PrseDRB*036"
names(meta)[names(meta) == "PrseDRB.037"] <- "PrseDRB*037"
names(meta)[names(meta) == "PrseDRB.038"] <- "PrseDRB*038"
names(meta)[names(meta) == "PrseDRB.039"] <- "PrseDRB*039"
names(meta)[names(meta) == "PrseDRB.040"] <- "PrseDRB*040"
names(meta)[names(meta) == "PrseDRB.041"] <- "PrseDRB*041"
names(meta)[names(meta) == "PrseDRB.043"] <- "PrseDRB*043"
names(meta)[names(meta) == "PrseDRB.044"] <- "PrseDRB*044"
names(meta)[names(meta) == "PrseDRB.045"] <- "PrseDRB*045"
names(meta)[names(meta) == "PrseDRB.046"] <- "PrseDRB*046"
names(meta)[names(meta) == "PrseDRB.047"] <- "PrseDRB*047"
names(meta)[names(meta) == "PrseDRB.048"] <- "PrseDRB*048"
names(meta)[names(meta) == "PrseDRB.049"] <- "PrseDRB*049"
names(meta)[names(meta) == "PrseDRB.050"] <- "PrseDRB*050"
names(meta)[names(meta) == "PrseDRB.051"] <- "PrseDRB*051"
names(meta)[names(meta) == "PrseDRB.052"] <- "PrseDRB*052"
names(meta)[names(meta) == "PrseDRB.053"] <- "PrseDRB*053"
names(meta)[names(meta) == "PrseDRB.054"] <- "PrseDRB*054"
names(meta)[names(meta) == "PrseDRB.055"] <- "PrseDRB*055"
names(meta)[names(meta) == "PrseDRB.056"] <- "PrseDRB*056"
names(meta)[names(meta) == "PrseDRB.058"] <- "PrseDRB*058"
names(meta)[names(meta) == "PrseDRB.061"] <- "PrseDRB*061"
names(meta)[names(meta) == "PrseDRB.062"] <- "PrseDRB*062"
names(meta)[names(meta) == "PrseDRB.063"] <- "PrseDRB*063"
names(meta)[names(meta) == "PrseDRB.064"] <- "PrseDRB*064"
names(meta)[names(meta) == "PrseDRB.065"] <- "PrseDRB*065"
names(meta)[names(meta) == "PrseDRB.066"] <- "PrseDRB*066"
names(meta)[names(meta) == "PrseDRB.068"] <- "PrseDRB*068"
names(meta)[names(meta) == "PrseDRB.070"] <- "PrseDRB*070"
names(meta)[names(meta) == "PrseDRB.071"] <- "PrseDRB*071"
names(meta)[names(meta) == "PrseDRB.073"] <- "PrseDRB*073"
names(meta)[names(meta) == "PrseDRB.078"] <- "PrseDRB*078"
names(meta)[names(meta) == "PrseDRB.079"] <- "PrseDRB*079"
names(meta)[names(meta) == "PrseDRB.080"] <- "PrseDRB*080"
names(meta)[names(meta) == "PrseDRB.001"] <- "PrseDRB*001"
names(meta)[names(meta) == "PrseDRB.002"] <- "PrseDRB*002"
names(meta)[names(meta) == "PrseDRB.010"] <- "PrseDRB*010"
names(meta)[names(meta) == "PrseDRB.020"] <- "PrseDRB*020"
names(meta)[names(meta) == "PrseDRB.022"] <- "PrseDRB*022"
names(meta)[names(meta) == "PrseDRB.035"] <- "PrseDRB*035"
names(meta)[names(meta) == "PrseDRB.042"] <- "PrseDRB*042"



### MHC constitution across landscapes ###


meta_c <- subset(meta, landscape == "forest")
meta_i <- subset(meta, landscape == "island")
meta_a <- subset(meta, landscape == "fragment")

meta_ci <- subset(meta, landscape != "fragment") # since ST7 is absent in forested fragments, we'll create this subset, too


# Overview, number of individuals with/without each ST per landscape

table(meta$ST1,meta$landscape)

table(meta$ST3,meta$landscape)

table(meta$ST4,meta$landscape)

table(meta$ST5,meta$landscape)

table(meta$ST6,meta$landscape)

table(meta$ST7,meta$landscape)

table(meta$ST8,meta$landscape)

table(meta$ST10,meta$landscape)

table(meta$ST11,meta$landscape)

table(meta$ST12,meta$landscape)

table(meta$ST13,meta$landscape)

table(meta$ST14,meta$landscape)

table(meta$ST15,meta$landscape)




### MHC allele frequency overall ###

# Mean NrAlleles in the population + standard deviation
mean(meta$NoDifAA_alleles)
sd(meta$NoDifAA_alleles)

# make metadata into long format
prse_meta_long <- gather(meta, allele, freq, "PrseDRB*001":"PrseDRB*079", factor_key=TRUE) # for alleles
prse_meta_long$freq <- factor(prse_meta_long$freq)


# Plot presence / absence
prse_meta_long1<-subset(prse_meta_long, freq != 0) # Several alleles not present in this data subset 

# Line at frequency = 5 (alleles present in >5 individuals are above this line)

ggplot(prse_meta_long1, aes(x = allele, fill=freq))+
  geom_histogram(stat="count")+
  labs(x=" ",y="Allele frequency")+
  scale_fill_manual(values=c("grey45"))+ 
  geom_hline(aes(yintercept=5),linetype="dashed",size=0.5)+
  ylab("relative frequency")+
  xlab("")+
  theme_bw(base_size=12)+
  theme(legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_text(size = 10))



## MHC allele frequency per landscape

prse_meta_long <- gather(meta, allele, freq, "PrseDRB*001":"PrseDRB*079", factor_key=TRUE) # for alleles

MHCsummary<-ddply(prse_meta_long, c("landscape","allele"), summarise, 
                  frequency=sum(freq)/length(allele), 
                  length=length(allele),
                  number=sum(freq))



### Supp Figure 4a ###

MHCsummary$allele <- gsub("[PrseDRB]","", MHCsummary$allele) # shorten allele names

suppfigure_4a <- MHCsummary %>%
  mutate(landscape = fct_relevel(landscape,"forest","island","fragment")) %>%
  ggplot(aes(allele, frequency, fill=landscape))+
  geom_bar(position="dodge", stat="identity") +
  theme_classic(base_size=14)+
  ylab("Relative frequency")+
  xlab("")+
  scale_fill_manual(labels = c("C", "I", "A"), values=c("chartreuse4", "#54bed6", "orange2"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size= 10))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")



# subset only alleles that show co-occurrence with a pathogen in at least one of the landscapes

MHCsummary_sub <- subset(MHCsummary, allele =="*025"| 
                           allele =="*018"|
                           allele =="*015"|
                           allele =="*014"|
                           allele =="*013"|
                           allele =="*012"|
                           allele =="*011"|
                           allele =="*010"|
                           allele =="*009"|
                           allele =="*008"|
                           allele =="*007"|
                           allele =="*006"|
                           allele =="*005"|
                           allele =="*004"|
                           allele =="*003"|
                           allele =="*002"|
                           allele =="*001")



### Figure 1b ###

figure1b_allele <- MHCsummary_sub %>%
  mutate(allele = fct_relevel(allele,"*025", "*018", "*015", "*014", "*013", "*012", "*011", "*010", "*009", "*008", "*007", "*006", "*005", "*004", "*003", "*002", "*001")) %>%
  mutate(landscape = fct_relevel(landscape,"forest","island","fragment")) %>%
  ggplot(aes(x = allele, y=frequency, fill=landscape)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_classic(base_size = 14)+
  ylab(" ")+
  xlab("")+
  scale_fill_manual(labels = c("C", "I", "A"), values=c("chartreuse4", "#54bed6", "orange2"))+  
  theme(axis.text.x = element_text(size= 12),axis.text.y = element_text(size= 12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none")+ 
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  scale_y_continuous(breaks=seq(0,1,0.2))




### Glms - do alleles differ in frequency across landscapes? ###

# control for season and sex

meta$landscape <- factor(meta$landscape, levels = c("forest", "fragment", "island")) # compare with continuous forests
meta$landscape <- factor(meta$landscape, levels = c("island", "forest", "fragment")) # compare with islands

names(meta)[names(meta) == "PrseDRB*001"] <- "PrseDRB_001"

m=glm(PrseDRB_001 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta)[names(meta) == "PrseDRB_001"] <- "PrseDRB*001"

names(meta)[names(meta) == "PrseDRB*002"] <- "PrseDRB_002"

m=glm(PrseDRB_002 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta)[names(meta) == "PrseDRB_002"] <- "PrseDRB*002"

names(meta)[names(meta) == "PrseDRB*003"] <- "PrseDRB_003"

m=glm(PrseDRB_003 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta)[names(meta) == "PrseDRB_003"] <- "PrseDRB*003"

names(meta)[names(meta) == "PrseDRB*004"] <- "PrseDRB_004"

m=glm(PrseDRB_004 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta)[names(meta) == "PrseDRB_004"] <- "PrseDRB*004"

names(meta)[names(meta) == "PrseDRB*005"] <- "PrseDRB_005"

m=glm(PrseDRB_005 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta)[names(meta) == "PrseDRB_005"] <- "PrseDRB*005"

names(meta)[names(meta) == "PrseDRB*006"] <- "PrseDRB_006"

m=glm(PrseDRB_006 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta)[names(meta) == "PrseDRB_006"] <- "PrseDRB*006"

names(meta)[names(meta) == "PrseDRB*007"] <- "PrseDRB_007"

m=glm(PrseDRB_007 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta)[names(meta) == "PrseDRB_007"] <- "PrseDRB*007"

names(meta)[names(meta) == "PrseDRB*008"] <- "PrseDRB_008"

m=glm(PrseDRB_008 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta)[names(meta) == "PrseDRB_008"] <- "PrseDRB*008"

names(meta)[names(meta) == "PrseDRB*009"] <- "PrseDRB_009"

m=glm(PrseDRB_009 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta)[names(meta) == "PrseDRB_009"] <- "PrseDRB*009"

names(meta)[names(meta) == "PrseDRB*010"] <- "PrseDRB_010"

m=glm(PrseDRB_010 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta)[names(meta) == "PrseDRB_010"] <- "PrseDRB*010"

names(meta)[names(meta) == "PrseDRB*011"] <- "PrseDRB_011"

m=glm(PrseDRB_011 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta)[names(meta) == "PrseDRB_011"] <- "PrseDRB*011"

names(meta)[names(meta) == "PrseDRB*012"] <- "PrseDRB_012"

m=glm(PrseDRB_012 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta)[names(meta) == "PrseDRB_012"] <- "PrseDRB*012"

names(meta)[names(meta) == "PrseDRB*013"] <- "PrseDRB_013"

m=glm(PrseDRB_013 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta)[names(meta) == "PrseDRB_013"] <- "PrseDRB*013"

names(meta)[names(meta) == "PrseDRB*014"] <- "PrseDRB_014"

m=glm(PrseDRB_014 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta_ci)[names(meta_ci) == "PrseDRB_014"] <- "PrseDRB*014"

names(meta_ci)[names(meta_ci) == "PrseDRB*015"] <- "PrseDRB_015"

table(meta$"PrseDRB*015",meta$landscape)   # only in C, I

m=glm(PrseDRB_015 ~ landscape + season + sex, family=binomial, data=meta_ci, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta_ci)[names(meta_ci) == "PrseDRB_015"] <- "PrseDRB*015"

names(meta_ci)[names(meta_ci) == "PrseDRB*016"] <- "PrseDRB_016"

table(meta$"PrseDRB*016",meta$landscape) # only in C, I

m=glm(PrseDRB_016 ~ landscape + season + sex, family=binomial, data=meta_ci, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta_ci)[names(meta_ci) == "PrseDRB_016"] <- "PrseDRB*016"

names(meta_ci)[names(meta_ci) == "PrseDRB*017"] <- "PrseDRB_017"

table(meta$"PrseDRB*017",meta$landscape)   # only in C, I

m=glm(PrseDRB_017 ~ landscape + season + sex, family=binomial, data=meta_ci, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta_ci)[names(meta_ci) == "PrseDRB_017"] <- "PrseDRB*017"


table(meta$"PrseDRB*018",meta$landscape) # only in the forest
table(meta$"PrseDRB*019",meta$landscape) # only in the forest

table(meta$"PrseDRB*020",meta$landscape) 

table(meta$"PrseDRB*021",meta$landscape) # only in the forest

table(meta$"PrseDRB*022",meta$landscape) 
table(meta$"PrseDRB*023",meta$landscape) 
table(meta$"PrseDRB*025",meta$landscape) 

names(meta_ci)[names(meta_ci) == "PrseDRB*025"] <- "PrseDRB_025"

m=glm(PrseDRB_025 ~ landscape + season + sex, family=binomial, data=meta_ci, na.action=na.fail)
summary(m)
plot(allEffects(m))

names(meta_ci)[names(meta_ci) == "PrseDRB_025"] <- "PrseDRB*025"

table(meta$"PrseDRB*026",meta$landscape) # only on islands
table(meta$"PrseDRB*027",meta$landscape) # only on islands
table(meta$"PrseDRB*028",meta$landscape) # only in fragments
table(meta$"PrseDRB*029",meta$landscape) # only in the forest
table(meta$"PrseDRB*031",meta$landscape) # only in the forest
table(meta$"PrseDRB*032",meta$landscape) # only in the forest

table(meta$"PrseDRB*038",meta$landscape) 
table(meta$"PrseDRB*039",meta$landscape)

table(meta$"PrseDRB*040",meta$landscape) # only in the forest
table(meta$"PrseDRB*043",meta$landscape) # only on islands
table(meta$"PrseDRB*044",meta$landscape) # only in the forest
table(meta$"PrseDRB*045",meta$landscape) # only in the forest
table(meta$"PrseDRB*047",meta$landscape) # only in fragments
table(meta$"PrseDRB*048",meta$landscape) # only on islands
table(meta$"PrseDRB*049",meta$landscape) # only on islands
table(meta$"PrseDRB*052",meta$landscape) # only on islands
table(meta$"PrseDRB*054",meta$landscape) # only in the forest
table(meta$"PrseDRB*055",meta$landscape) # only on islands
table(meta$"PrseDRB*056",meta$landscape) # only in the forest
table(meta$"PrseDRB*058",meta$landscape) # only in the forest
table(meta$"PrseDRB*070",meta$landscape) # only in fragments
table(meta$"PrseDRB*071",meta$landscape) # only on islands
table(meta$"PrseDRB*078",meta$landscape) # only on islands
table(meta$"PrseDRB*079",meta$landscape) # only on islands





### MHC ST frequency overall ###


# Mean NrST in the population
mean(meta$NoOfSTs)  
sd(meta$NoOfSTs)   

# make metadata into long format
prse_st_long <- gather(meta, st, freq, ST1:ST15, factor_key=TRUE) # for STs

# Line at frequency = 10 (STs present in >10 individuals are above this line)

prse_st_long1<-subset(prse_st_long, freq != 0) 

ggplot(prse_st_long1, aes(x = st, fill=as.factor(freq)))+
  geom_histogram(stat="count")+
  labs(x=" ",y="ST frequency")+
  scale_fill_manual(values=c("grey45"))+ 
  geom_hline(aes(yintercept=10),linetype="dashed",size=0.5)+
  theme_bw(base_size=14)+
  theme(legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_text(size = 11))



# relative ST frequency in the population

prse_meta_long_st <- gather(meta, st, freq, ST1:ST15, factor_key=TRUE) # for STs

MHCsummary.st<-ddply(prse_meta_long_st, c("landscape","st"), summarise, 
                     frequency=sum(freq)/length(st), 
                     length=length(st),
                     number=sum(freq))


## Supp Figure 4b ##

suppfigure_4b <- MHCsummary.st %>%
  mutate(landscape = fct_relevel(landscape,"forest","island","fragment")) %>%
  mutate(st = fct_relevel(st,"ST1","ST3","ST4","ST5","ST6","ST7","ST8","ST10","ST11","ST12","ST13","ST14","ST15")) %>%
  ggplot(aes(x = st , y= frequency, fill=landscape)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_classic(base_size = 14)+
  ylab("Relative frequency")+
  xlab("")+
  scale_fill_manual(labels = c("C", "I", "A"), values=c("chartreuse4", "#54bed6", "orange2"))+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")



## Figure 1b ##

MHCsummary.st_sub <- subset(MHCsummary.st, st =="ST15" |
                              st =="ST14" |
                              st =="ST13" |
                              st =="ST12" |
                              st =="ST11" |
                              st =="ST8" |
                              st =="ST7" |
                              st =="ST6" |
                              st =="ST5")


figure1b_st <- MHCsummary.st_sub %>%
  mutate(landscape = fct_relevel(landscape,"forest","island","fragment")) %>%
  mutate(st = fct_relevel(st,"ST15","ST14","ST13","ST12","ST11","ST8","ST7","ST6","ST5")) %>%
  ggplot(aes(x = st , y=frequency, fill=landscape)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_classic(base_size = 14)+
  ylab("relative frequency")+
  xlab("")+
  scale_fill_manual(labels = c("C", "I", "A"), values=c("chartreuse4", "#54bed6", "orange2"))+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  
  theme(legend.position = "none")+ 
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  scale_y_continuous(breaks=seq(0,1,0.2))



#####################
### Supp Figure 4 ###
#####################

ggarrange(suppfigure_4a, suppfigure_4b, nrow=2, labels = c("A", "B"))



###########################
### Part of Figure 2B/C ###
###########################

ggarrange(figure1b_allele + coord_flip(), figure1b_st + coord_flip(), nrow=2, legend = "none", widths = c(1.3, 1))





## Glms - do STs differ in frequency across landscapes?

# control for season and sex
# run once for both orders to include all comparisons

meta$landscape <- factor(meta$landscape, levels = c("forest", "fragment", "island"))
meta$landscape <- factor(meta$landscape, levels = c("island", "forest", "fragment")) 

# ST1 is present in all individuals, thus not tested

m=glm(ST3 ~ landscape + season + sex, family=binomial, data=meta_ci, na.action=na.fail)
summary(m)
plot(allEffects(m))

m=glm(ST5 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

m=glm(ST6 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

m=glm(ST7 ~ landscape + season + sex, family=binomial, data=meta_ci, na.action=na.fail)
summary(m)
plot(allEffects(m))

m=glm(ST8 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

m=glm(ST11 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

m=glm(ST12 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

m=glm(ST13 ~ landscape + season + sex, family=binomial, data=meta_ci, na.action=na.fail)
summary(m)
plot(allEffects(m))

m=glm(ST14 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))

m=glm(ST15 ~ landscape + season + sex, family=binomial, data=meta, na.action=na.fail)
summary(m)
plot(allEffects(m))




## Neutral and MHC diversity indices across landscapes

# Add neutral genetic diversity 

# Data from: Not all effects are direct: How habitat disturbance affects multiple components of biodiversity. Nina Isabell Schwensow, Alexander Christoph Heni,
# Julian Schmid, B. Karina Montero, Stefan Dominik Brändel, Tanja Katharina Halczok, Gerd Mayer, Gloria Fackelmann, Kerstin Wilhelm, Simone Sommer

# load file

SEM_dataset <- read.csv("D:\\Ramona_F\\2PhD Project\\Prse_Ramona\\SEM_Dataset.csv")

# remove the plantation individuals

SEM_dataset <- subset(SEM_dataset, landscape != "plantation")

# rename the landscape variable similar to the other df

levels(SEM_dataset$landscape) <- list(forest ="continuous.forest" ,fragment="forest.fragment",island="island") 


# Now plot genetic diversity

my_comparisons <- list(c("forest", "fragment"),c("forest", "island"),c("fragment", "island"))



# GenDiv
neutral <- SEM_dataset %>%
  mutate(landscape = fct_relevel(landscape,"forest","island","fragment")) %>%
  ggplot(aes(x = landscape, y = GenDiv, fill = landscape))+
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(fill = landscape), size=1, colour = "black", width =0.01)+
  ggsignif::geom_signif(comparisons = list(c("forest", "fragment"),c("forest", "island"),c("fragment", "island")),
                        y_position= c(1,1.15,1.25),
                        test = "wilcox.test",
                        map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05))+
  scale_fill_manual(values=c("chartreuse4","#54bed6","orange2"))+
  theme_bw(base_size=13)+xlab("")+
  ylab("GenDiv")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  scale_x_discrete(labels=c("forest" = "C", "island" = "I", "fragment" = "A"))+
  ylim(0, 1.4)+
  theme(axis.title.y = element_text(size = 14), axis.ticks.x = element_blank(),axis.text.x = element_blank())


# NrAlleles
p.diff.allele <-  meta %>%
  mutate(landscape = fct_relevel(landscape,"forest","island","fragment")) %>%
  ggplot(aes(x = landscape, y = NoDifAA_alleles, fill = landscape))+
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(fill = landscape), size=1, colour = "black", width =0.01)+
  ggsignif::geom_signif(comparisons = list(c("forest", "fragment"),c("forest", "island"),c("fragment", "island")),
                        y_position= c(9,10,10.5),
                        test = "wilcox.test",
                        map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05))+
  scale_fill_manual(values=c("chartreuse4", "#54bed6","orange2"))+
  theme_bw(base_size=13)+xlab("")+
  ylab(expression(Nr[ Alleles]))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  scale_x_discrete(labels=c("forest" = "C", "island" = "I", "fragment" = "A"))+
  ylim(0, 11.5)+
  theme(axis.title.y = element_text(size = 16), axis.ticks.x = element_blank(),axis.text.x = element_blank())


# NrST
p.st <-  meta %>%
  mutate(landscape = fct_relevel(landscape,"forest","island","fragment")) %>%
  ggplot(aes(x = landscape, y = NoOfSTs, fill = landscape))+
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(fill = landscape), size=1, colour = "black", width =0.01)+
  ggsignif::geom_signif(comparisons = list(c("forest", "fragment"),c("forest", "island"),c("fragment", "island")),
                        y_position= c(9,10,10.5),
                        test = "wilcox.test",
                        map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05))+
  scale_fill_manual(values=c("chartreuse4", "#54bed6","orange2"))+
  theme_bw(base_size=13)+xlab("")+
  ylab(expression(Nr[ ST]))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  scale_x_discrete(labels=c("forest" = "C", "island" = "I", "fragment" = "A"))+
  ylim(0, 11.5)+
  theme(axis.title.y = element_text(size = 16), axis.ticks.x = element_blank(),axis.text.x = element_blank())


# PdistAllele
p.pdist.sum <-  meta %>%
  mutate(landscape = fct_relevel(landscape,"forest","island","fragment")) %>%
  ggplot(aes(x = landscape, y = AA_pdist_sum, fill = landscape))+
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(fill = landscape), size=1, colour = "black", width =0.01)+
  ggsignif::geom_signif(comparisons = list(c("forest", "fragment"),c("forest", "island"),c("fragment", "island")),
                        y_position= c(9,10,10.5),
                        test = "wilcox.test",
                        map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05))+
  scale_fill_manual(values=c("chartreuse4", "#54bed6","orange2"))+
  theme_bw(base_size=13)+xlab("")+
  ylab(expression(Pdist[ Allele]))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  scale_x_discrete(labels=c("forest" = "C", "island" = "I", "fragment" = "A"))+
  ylim(0, 11.5)+
  theme(axis.title.y = element_text(size = 16), axis.ticks.x = element_blank(),axis.text.x = element_blank())


# PdistPSS
p.pdistPSS.sum <-  meta %>%
  mutate(landscape = fct_relevel(landscape,"forest","island","fragment")) %>%
  ggplot(aes(x = landscape, y = PSS_pdist_sum, fill = landscape))+
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(fill = landscape), size=1, colour = "black", width =0.01)+
  ggsignif::geom_signif(comparisons = list(c("forest", "fragment"),c("forest", "island"),c("fragment", "island")),
                        y_position= c(23,26,28),
                        test = "wilcox.test",
                        map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05))+
  scale_fill_manual(values=c("chartreuse4", "#54bed6","orange2"))+
  theme_bw(base_size=13)+xlab("")+
  ylab(expression(Pdist[ PSS]))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  scale_x_discrete(labels=c("forest" = "C", "island" = "I", "fragment" = "A"))+
  ylim(0, 30)+
  theme(axis.title.y = element_text(size = 16), axis.ticks.x = element_blank(),axis.text.x = element_blank())


# for legend
meta %>%
  mutate(landscape = fct_relevel(landscape,"forest","island","fragment")) %>%
  ggplot(aes(x = landscape, y = PSS_pdist_sum, fill = landscape))+
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(fill = landscape), size=1, colour = "black", width =0.01)+
  ggsignif::geom_signif(comparisons = list(c("forest", "fragment"),c("forest", "island"),c("fragment", "island")),
                        y_position= c(23,26,28),
                        test = "wilcox.test",
                        map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05))+
  scale_fill_manual(labels = c("C", "I","A"),values=c("chartreuse4", "#54bed6","orange2"))+
  theme_bw()+xlab("")+
  theme_bw(base_size=15)+
  ylab(expression(Pdist[ PSS]))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  scale_x_discrete(labels=c("forest" = "C", "island" = "I", "fragment" = "A"))+ylim(0, 32)



### Figure 2a ###

figure2a <- ggarrange(neutral, p.diff.allele, p.st, p.pdist.sum, p.pdistPSS.sum, legend ="none", nrow=1)
figure2a






# SNPs vs. MHC

# merge snp and mhc datasets 

merge_test <- merge(meta, SEM_dataset, by="field_id")
# we have joined data for SNP and MHC diversity of 95 individuals 


# Is neutral genetic diversity correlated with MHC diversity estimates in spiny rats?

snp.mhc <- ggplot(data = merge_test, mapping = aes(x = GenDiv, y = NoDifAA_alleles)) + 
  geom_point(size = 3) +
  geom_smooth()+
  ylim(0,10) +
  stat_cor(method = "pearson", y.position=10)+
  ylab(expression(Nr[ Alleles]))+
  theme_bw(base_size=15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

snp.mhc1 <- ggplot(data = merge_test, mapping = aes(x = GenDiv, y = NoOfSTs)) + 
  geom_point(size = 3) +
  geom_smooth()+
  ylim(0,10) +
  stat_cor(method = "pearson", y.position=10)+
  ylab(expression(Nr[ ST]))+
  theme_bw(base_size=15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

snp.mhc2 <- ggplot(data = merge_test, mapping = aes(x = GenDiv, y = AA_pdist_sum)) + 
  geom_point(size = 3) +
  geom_smooth()+
  ylim(0,10) +
  stat_cor(method = "pearson", y.position=10)+
  ylab(expression(Pdist[ Allele]))+
  theme_bw(base_size=15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

snp.mhc3 <- ggplot(data = merge_test, mapping = aes(x = GenDiv, y = PSS_pdist_sum)) + 
  geom_point(size = 3) +
  geom_smooth()+
  ylim(0,25) +
  stat_cor(method = "pearson", y.position=20)+
  ylab(expression(Pdist[ PSS]))+
  theme_bw(base_size=15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


######################
### Supp Figure S3 ###
######################

ggarrange(snp.mhc, snp.mhc1, snp.mhc2, snp.mhc3)









## Glms ##

# first check for colinearitry among continuous variables: Variance Inflation Factor and test for multicollinearity
# MHC variables

# calculate vif for a subset of mhc diversity variables 
# only distances among amino acids, not on nucleotide level included

m <- lm(NNI_helm ~ NoDifAA_alleles+NoOfSTs+AA_pdist_sum+AA_jttdist_sum+PSS_pdist_sum, data=meta)
car::vif(m)

var <- meta[,c("NoDifAA_alleles", "NoOfSTs", "AA_pdist_sum", "AA_jttdist_sum", "PSS_pdist_sum")]
vif(var)

var <- meta[,c("NoDifAA_alleles", "NoOfSTs", "AA_pdist_sum")]
vif(var)
#         Variables       VIF
# 1 NoDifAA_alleles  9.547281
# 2         NoOfSTs  3.034448
# 3    AA_pdist_sum 10.039582


# the following are methods which exclude highly colinear variables using distinct strategies
# vifcor, first find a pair of variables which has the maximum linear correlation (greater than th), and exclude one of them which has greater VIF. 
#The procedure is repeated untill no variable with a high corrrelation coefficient (grater than threshold) with other variables remains. 
v1 <- vifcor(var, th=0.9) # identify collinear variables that should be excluded
v1
3 variables from the 5 input variables have collinearity problem: 
  
  AA_pdist_sum PSS_pdist_sum AA_jttdist_sum 

After excluding the collinear variables, the linear correlation coefficients ranges between: 
  min correlation ( NoOfSTs ~ NoDifAA_alleles ):  0.801258 
max correlation ( NoOfSTs ~ NoDifAA_alleles ):  0.801258 

---------- VIFs of the remained variables -------- 
  Variables      VIF
1 NoDifAA_alleles 2.793408
2         NoOfSTs 2.793408


# vifstep calculate VIF for all variables, exclude one with highest VIF (greater than threshold), repeat the procedure untill no variables with VIF greater than th remains.
v2 <- vifstep(var, th=10) # identify collinear variables that should be excluded
v2
2 variables from the 5 input variables have collinearity problem: 
  
  AA_pdist_sum PSS_pdist_sum 

After excluding the collinear variables, the linear correlation coefficients ranges between: 
  min correlation ( NoOfSTs ~ NoDifAA_alleles ):  0.801258 
max correlation ( AA_jttdist_sum ~ NoDifAA_alleles ):  0.9418751 

---------- VIFs of the remained variables -------- 
  Variables      VIF
1 NoDifAA_alleles 9.356420
2         NoOfSTs 2.954063
3  AA_jttdist_sum 9.369180



## draw a correlation matrix 

var <- meta[,c("NoDifAA_alleles", "NoOfSTs", "AA_pdist_sum", "PSS_pdist_sum")]
vif(var)

mhc.cor = cor(var, method = "pearson")
mhc.cor

#install.packages("corrplot")
library(corrplot)
dev.off()


######################
### Supp Figure S6 ###
######################

corrplot(mhc.cor)


#install.packages("sjPlot")
library(sjPlot)

# Now model



##########################################################################################
###################### Does MHC diversity explain pathogen richness ######################
##########################################################################################

# Model per landscape

######################## C #####################

table(meta_c$estimated_habitat_size_ha)

hist(meta_c$NNI_helm)
m1 <- glm(NNI_helm ~ season + sex + NoDifAA_alleles, meta_c, family = Gamma(link = "inverse"))
summary(m1)
plot(allEffects(m1))
m.tab <- sjPlot::plot_model(m1)
sjPlot::tab_model(m1)
AIC(m1) # 196.9353


meta_c$NNI_helm <- as.numeric(meta_c$NNI_helm)

m1 <- glm(NNI_helm ~ season + sex + NoOfSTs , meta_c, family = Gamma(link = "inverse"))
summary(m1)
plot(allEffects(m1))
m.tab1 <- sjPlot::plot_model(m1)
sjPlot::tab_model(m1)
AIC(m1)

m1 <- glm(NNI_helm ~ season + sex + AA_pdist_sum, meta_c, family = Gamma(link = "inverse"))
summary(m1)
plot(allEffects(m1))
m.tab2 <- sjPlot::plot_model(m1)
sjPlot::tab_model(m1)
AIC(m1)

m1 <- glm(NNI_helm ~ season + sex + PSS_pdist_sum, meta_c, family = Gamma(link = "inverse"))
summary(m1)
plot(allEffects(m1))
sjPlot::plot_model(m1)
sjPlot::tab_model(m1)
AIC(m1)

m1 <- glm(NNI_helm ~ season + sex + AA_jttdist_sum, meta_c, family = Gamma(link = "inverse"))
summary(m1)
plot(allEffects(m1))
sjPlot::plot_model(m1)
sjPlot::tab_model(m1)
AIC(m1)



hist(meta_c$NVI)
m2 <- glm(NVI ~ season + sex + NoDifAA_alleles, meta_c, family = gaussian(link = "identity"))
summary(m2)
plot(allEffects(m2))
sjPlot::plot_model(m2)
sjPlot::tab_model(m2)
AIC(m2)

m2 <- glm(NVI ~ season + sex + NoOfSTs, meta_c, family = gaussian(link = "identity"))
summary(m2)
plot(allEffects(m2))
sjPlot::plot_model(m2)
sjPlot::tab_model(m2)
AIC(m2)

m2 <- glm(NVI ~ season + sex + AA_pdist_sum, meta_c, family = gaussian(link = "identity"))
summary(m2)
plot(allEffects(m2))
sjPlot::plot_model(m2)
sjPlot::tab_model(m2)
AIC(m2)

m2 <- glm(NVI ~ season + sex + PSS_pdist_sum, meta_c, family = gaussian(link = "identity"))
summary(m2)
plot(allEffects(m2))
sjPlot::plot_model(m2)
sjPlot::tab_model(m2)
AIC(m2)

m2 <- glm(NVI ~ season + sex + AA_jttdist_sum, meta_c, family = gaussian(link = "identity"))
summary(m2)
plot(allEffects(m2))
sjPlot::plot_model(m2)
sjPlot::tab_model(m2)
AIC(m2)





######################## I #####################

meta_i_red <- subset(meta_i, NNI_helm != "0")

table(meta_i_red$NNI_helm) 

hist(meta_i_red$NNI_helm)
m1 <- glm(NNI_helm ~ season + sex + NoDifAA_alleles, meta_i_red, family = Gamma(link = "inverse"))
summary(m1)
plot(allEffects(m1))
sjPlot::plot_model(m1)
sjPlot::tab_model(m1)
AIC(m1)

m1 <- glm(NNI_helm ~ season + sex + NoOfSTs , meta_i_red, family = Gamma(link = "inverse"))
summary(m1)
plot(allEffects(m1))
sjPlot::plot_model(m1)
sjPlot::tab_model(m1)
AIC(m1)

m1 <- glm(NNI_helm ~ season + sex + AA_pdist_sum, meta_i_red, family = Gamma(link = "inverse"))
summary(m1)
plot(allEffects(m1))
sjPlot::plot_model(m1)
sjPlot::tab_model(m1)
AIC(m1)

m1 <- glm(NNI_helm ~ season + sex + PSS_pdist_sum, meta_i_red, family = Gamma(link = "inverse"))
summary(m1)
plot(allEffects(m1))
sjPlot::plot_model(m1)
sjPlot::tab_model(m1)
AIC(m1)




table(meta_i_red$NVI, meta_i_red$season)

hist(meta_i_red$NVI)
m2 <- glm(NVI ~ season + sex + NoDifAA_alleles, meta_i_red, family = gaussian(link = "identity"))
summary(m2)
plot(allEffects(m2))
sjPlot::plot_model(m2)
sjPlot::tab_model(m2)
AIC(m2)

m2 <- glm(NVI ~ season + sex + NoOfSTs, meta_i_red, family = gaussian(link = "identity"))
summary(m2)
plot(allEffects(m2))
sjPlot::plot_model(m2)
sjPlot::tab_model(m2)
AIC(m2)

m2 <- glm(NVI ~ season + sex + AA_pdist_sum, meta_i_red, family = gaussian(link = "identity"))
summary(m2)
plot(allEffects(m2))
sjPlot::plot_model(m2)
sjPlot::tab_model(m2)
AIC(m2)

m2 <- glm(NVI ~ season + sex + PSS_pdist_sum, meta_i_red, family = gaussian(link = "identity"))
summary(m2)
plot(allEffects(m2))
sjPlot::plot_model(m2)
sjPlot::tab_model(m2)
AIC(m2)




######################## A #####################

meta_a <- subset(meta_a, NNI_helm != "0")

hist(meta_a$NNI_helm)
m1 <- glm(NNI_helm ~ season + sex + NoDifAA_alleles, meta_a, family = Gamma(link = "inverse"))
summary(m1)
plot(allEffects(m1))
sjPlot::plot_model(m1)
sjPlot::tab_model(m1)
AIC(m1)

m1 <- glm(NNI_helm ~ season + sex + NoOfSTs , meta_a, family = Gamma(link = "inverse"))
summary(m1)
plot(allEffects(m1))
sjPlot::plot_model(m1)
sjPlot::tab_model(m1)
AIC(m1)

m1 <- glm(NNI_helm ~ season + sex + AA_pdist_sum, meta_a, family = Gamma(link = "inverse"))
summary(m1)
plot(allEffects(m1))
sjPlot::plot_model(m1)
sjPlot::tab_model(m1)
AIC(m1)

m1 <- glm(NNI_helm ~ season + sex + PSS_pdist_sum, meta_a, family = Gamma(link = "inverse"))
summary(m1)
plot(allEffects(m1))
sjPlot::plot_model(m1)
sjPlot::tab_model(m1)
AIC(m1)




hist(meta_a$NVI)
m2 <- glm(NVI ~ season + sex + NoDifAA_alleles, meta_a, family = gaussian(link = "identity"))
summary(m2)
plot(allEffects(m2))
sjPlot::plot_model(m2)
sjPlot::tab_model(m2)
AIC(m2)

m2 <- glm(NVI ~ season + sex + NoOfSTs, meta_a, family = gaussian(link = "identity"))
summary(m2)
plot(allEffects(m2))
sjPlot::plot_model(m2)
sjPlot::tab_model(m2)
AIC(m2)

m2 <- glm(NVI ~ season + sex + AA_pdist_sum, meta_a, family = gaussian(link = "identity"))
summary(m2)
plot(allEffects(m2))
sjPlot::plot_model(m2)
sjPlot::tab_model(m2)
AIC(m2)

m2 <- glm(NVI ~ season + sex + PSS_pdist_sum, meta_a, family = gaussian(link = "identity"))
summary(m2)
plot(allEffects(m2))
sjPlot::plot_model(m2)
sjPlot::tab_model(m2)
AIC(m2)






#####################################################################################################################################################
########################################################## SINGLE PATHOGENS #########################################################################
#####################################################################################################################################################


# arrange data for plotting
pathogens_long <- gather(meta, pathogen, Count, N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, N13, Try, HdV, HpV, PbV, PcV, factor_key = TRUE)
table(pathogens_long$pathogen, pathogens_long$landscape)

pathogens_long <- subset(pathogens_long, Count != "0") # remove zero counts

table(pathogens_long$pathogen, pathogens_long$landscape)
table(pathogens_long$Count, pathogens_long$landscape)

colorlist <- c("#63ab00",
               "#5f4ed4",
               "#f0b907",
               "#4d78ff",
               "#69de78",
               "#a034bb",
               "#4a5d00",
               "#ff80f0",
               "#e5531a",
               "#789aff",
               "#a80429",
               "#384f94",
               "#e2c28c",
               "#7d3b69",
               "#ff7b89",
               "#e1b6f9",
               "#813f37",
               "#d295bb")

# plot
ggplot(pathogens_long, aes(x =landscape, fill=pathogen))+
  geom_bar(position="stack", color = "black")+
  labs(x=" ",y="Number of individuals")+
  theme_classic2(base_size=14)+
  xlab("")+
  scale_fill_manual(values = colorlist)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.x = element_blank(), axis.text.x=element_text(size=14))





pathogens_long2 <- gather(meta, pathogen, Count, N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, N13, Try, HdV, HpV, PbV, PcV, factor_key = TRUE)

ggplot(pathogens_long, aes(x =pathogen, fill=landscape))+
  geom_bar(position="identity")+
  scale_fill_manual(values=c("chartreuse4","#54bed6", "orange2"))+ 
  labs(x=" ",y="Number of individuals")+
  theme_classic2(base_size=14)+
  xlab("")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.x = element_blank(), axis.text.x=element_text(size=14))

pathogens_long$landscape <- relevel(pathogens_long$landscape, "forest")

ggplot(pathogens_long, aes(x =pathogen, fill=landscape))+
  geom_bar(position="identity", color="black")+
  scale_fill_manual(values=c("chartreuse4","#54bed6", "orange2"))+ 
  labs(x=" ",y="Number of infected individuals")+
  theme_bw(base_size=14)+
  xlab("")+
  facet_wrap(~landscape, nrow=3)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.x = element_blank(), axis.text.x=element_text(size=11))+
  theme(strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"))




ggplot(pathogens_long, aes(landscape, group = pathogen)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") + 
  scale_y_continuous(labels=scales::percent) +
  ylab("Prevalence [%]") +
  facet_grid(~pathogen)



ggplot(pathogens_long, aes(pathogen, group = landscape, fill = landscape)) + 
  geom_bar(aes(y = ..prop..),  stat="count",position="dodge", color = "black") + 
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(values=c("chartreuse4","#54bed6", "orange2"))+ 
  ylab("Prevalence [%]") +
  theme_classic(base_size=14)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.x = element_blank(), axis.text.x=element_text(size=11))


library(plyr)

pathogens_long2 <- gather(meta, pathogen, Count, N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, N13, Try, HdV, HpV, PbV, PcV, factor_key = TRUE)

pathogen_summary <-ddply(pathogens_long2, c("landscape","pathogen"), summarise, 
                         frequency=sum(Count)/length(landscape),
                         length=length(landscape), 
                         number=sum(Count))

######################
### Supp Figure S2 ###
######################

pathogen_summary %>%
  mutate(landscape = fct_relevel(landscape,"forest","island","fragment")) %>%
  ggplot(aes(pathogen, frequency, fill = landscape)) + 
  geom_bar(stat="identity", color = "black", position="dodge")+
  scale_fill_manual(values=c("chartreuse4","#54bed6", "orange2"), labels =c("C", "I", "A"))+ 
  ylab("Prevalence [%]") +
  xlab("")+
  theme_classic(base_size=13)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks.x = element_blank(), axis.text.x=element_text(size=11))







################################################################################################################
##################################### Pathogen pressure across landscapes ######################################
################################################################################################################


table(meta$N1, meta$landscape)


# 1) Hdv

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(HdV~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(HdV~landscape,meta, family = binomial(link=logit))
summary(m1)



# 2) Hepacivirus

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(HpV~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(HpV~landscape,meta, family = binomial(link=logit))
summary(m1)



# 3) Picobirnavirus

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(PbV~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(PbV~landscape,meta, family = binomial(link=logit))
summary(m1)




# 4) Picornavirus

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(PcV~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(PcV~landscape,meta, family = binomial(link=logit))
summary(m1)




# 5) Trypanosoma

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(Try~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(Try~landscape,meta, family = binomial(link=logit))
summary(m1)




# 6) N1

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(N1~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(N1~landscape,meta, family = binomial(link=logit))
summary(m1)




# 7) N2 

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(N2~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(N2~landscape,meta, family = binomial(link=logit))
summary(m1)




# 8) N3

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(N3~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(N3~landscape,meta, family = binomial(link=logit))
summary(m1)




# 9) N4

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(N4~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(N4~landscape,meta, family = binomial(link=logit))
summary(m1)




# 10) N5

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(N5~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(N5~landscape,meta, family = binomial(link=logit))
summary(m1)




# 11) N6

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(N6~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(N6~landscape,meta, family = binomial(link=logit))
summary(m1)




# 12) N7

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(N7~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(N7~landscape,meta, family = binomial(link=logit))
summary(m1)




# 13) N8 

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(N8~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(N8~landscape,meta, family = binomial(link=logit))
summary(m1)



# 14) N9

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(N9~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(N9~landscape,meta, family = binomial(link=logit))
summary(m1)




# 15) N10

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(N10~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(N10~landscape,meta, family = binomial(link=logit))
summary(m1)




# 16) N11

table(meta$landscape, meta$N11)

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(N11~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(N11~landscape,meta, family = binomial(link=logit))
summary(m1)




# 17) N12

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(N12~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(N12~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "island")




# 18) N13

meta$landscape <- relevel(meta$landscape, "forest")

m1 <- glm(N13~landscape,meta, family = binomial(link=logit))
summary(m1)

meta$landscape <- relevel(meta$landscape, "fragment")

m1 <- glm(N13~landscape,meta, family = binomial(link=logit))
summary(m1)





####################################################################################################################################
########################################################### Co-occurrence ##########################################################
####################################################################################################################################


## Subset data

# Continuous forest
prse_meta_C <- subset(meta, landscape == "continuous.forest")

# Island
prse_meta_I  <- subset(meta, landscape == "island")

# Forest fragments in argicultural matrix
prse_meta_A <- subset(meta, landscape == "forest.fragment")


##############################################################################################
########################################## Alleles ###########################################
##############################################################################################


####################################### Continuous Forest ####################################


input_alleles_c <- subset(prse_meta_C, select=c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11","N12","N13","Try","HdV","HpV","PbV","PcV",
                                                "PrseDRB*001","PrseDRB*002","PrseDRB*003","PrseDRB*004","PrseDRB*005","PrseDRB*006","PrseDRB*007","PrseDRB*008","PrseDRB*009","PrseDRB*010","PrseDRB*011","PrseDRB*012",
                                                "PrseDRB*013","PrseDRB*014","PrseDRB*015","PrseDRB*016","PrseDRB*017","PrseDRB*018","PrseDRB*019","PrseDRB*020","PrseDRB*021","PrseDRB*022","PrseDRB*023","PrseDRB*025",
                                                "PrseDRB*029","PrseDRB*031","PrseDRB*032","PrseDRB*039","PrseDRB*040","PrseDRB*044","PrseDRB*045","PrseDRB*054","PrseDRB*056","PrseDRB*058"))

# loop to automate co-occurences by random subset

dataALL <- NULL
input <- input_alleles_c

start <- 1
end <- 1000

cutoff_min <- 4
cutoff_max <- 46

i = start

set.seed(123)

for(i in start:end){
  data_subset <-  sample_n(input, 50)
  transposed_subset <- as.data.frame(t(data_subset[, (-1)]))
  transposed_subset$count <- rowSums(transposed_subset)
  filtered_subset <- transposed_subset %>%
    filter(count > cutoff_min) 
  filtered_subset <- filtered_subset %>%
    filter(count < cutoff_max) 
  Input_matrix <- filtered_subset[, c(1:ncol(filtered_subset)-1)]
  Co <- cooccur(mat = as.matrix(Input_matrix), 
                type = "spp_site", 
                thresh = FALSE, 
                spp_names = TRUE)
  effectsize_Co  <- effect.sizes(Co, standardized = TRUE, matrix = FALSE)     # standardized effect sizes
  names(effectsize_Co) <- c("sp1_name", "sp2_name", "effects")
  
  Complete_Co<-prob.table(Co) %>%
    filter(sp1_name =="N1" | sp1_name=="N2"  | sp1_name=="N3" | sp1_name=="N4" | sp1_name=="N5" | sp1_name=="N6" | sp1_name=="N7" | sp1_name=="N8" | sp1_name=="N9" | sp1_name=="N10" | sp1_name=="N11" | sp1_name=="N12" | sp1_name=="N13" | sp1_name=="Try" | sp1_name=="HdV" | sp1_name=="HpV" | sp1_name=="PbV" | sp1_name=="PcV") %>%
    left_join(effectsize_Co, by = c("sp1_name", "sp2_name"))
  
  Co_left <- Complete_Co %>%
    filter(p_lt < 0.1)
  Co_right <- Complete_Co %>%
    filter(p_gt < 0.1)
  
  Final_Co=rbind(Co_left,Co_right)
  
  Final_Co$obs.v.exp <- round((Final_Co$exp_cooccur-Final_Co$obs_cooccur)*-1, digits = 3)
  Final_Co$NR_IDs <- ncol(Input_matrix)
  Final_Co$round <- i
  dataALL=rbind(dataALL,Final_Co)
  
  i  = as.numeric (i + 1)
  
}


# Now filter only significant results

Sign_left <- dataALL %>%
  filter(p_lt < 0.05)
Sign_right <- dataALL %>%
  filter(p_gt < 0.05)

dataFilt=rbind(Sign_left,Sign_right)

res_Co_alleles <- dataFilt

#write.csv(res_Co_alleles, "res_alleles_C1000.csv")
#res_Co_alleles <- read.csv("res_alleles_C1000.csv")



# Draw heatmap of the results


heatmap.dat.allele <- res_Co_alleles %>%
  mutate(Asterisks = ifelse(p_gt  <= 0.001 |p_lt  <= 0.001, "***",
                            ifelse(p_gt <= 0.01 | p_lt  <= 0.01, "**",
                                   ifelse(p_gt  <= 0.05 | p_lt  <= 0.05,  "*", NA))))

head(heatmap.dat.allele)

table(heatmap.dat.allele$effects)

ggplot(heatmap.dat.allele, aes(sp2_name, sp1_name, fill=effects)) + 
  geom_tile()+
  #scale_fill_gradient2(high="#D53E4F", mid="white", low="#3288BD", limits = c(-0.1, 0.9))+
  scale_fill_gradient2(high="#D53E4F", mid="white", low="#3288BD")+
  facet_wrap(~sp1_name) + 
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(angle = 90))+
  geom_text(data = heatmap.dat.allele, aes(x = sp2_name, y = sp1_name, label = Asterisks), size = 4)



# Barplot that shows number of associations found across 


tbl <- data.frame(table(res_Co_alleles$round))

ggplot(tbl, aes(factor(Var1), Freq, fill = Var1)) +     
  geom_bar(stat="identity", position = "dodge")+
  ylab("Nr of associations")+
  xlab("Co-occurrence round")+
  theme_classic(base_size=14)+
  theme(legend.position = "NONE")

mean(tbl$Freq)
median(tbl$Freq)

table(tbl$Freq)



# How stable are the associations, how often are they found?

head(res_Co_alleles)

res_Co_alleles$association <- paste(res_Co_alleles$sp1_name, res_Co_alleles$sp2_name)

tbl2 <- data.frame(table(res_Co_alleles$association))

# exclude spurious associations (present in < 5% of co-occurrence rounds)

tbl2 %>%
  filter(Freq > 50)

ggplot(tbl2, aes(factor(Var1), Freq, fill = Var1)) +     
  geom_bar(stat="identity", position = "dodge")+
  ylab("Nr of times the association was found in 1000 repeated rounds")+
  xlab("Association pair")+
  theme_classic(base_size=14)+
  theme(legend.position = "NONE") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#write.csv(data.frame(tbl2), "robustness_associations_C1000.csv")



####################################### Forested islands ####################################

input_alleles_i <- subset(prse_meta_I, select=c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11","N12","N13","Try","HdV","HpV","PbV","PcV",
                                                "PrseDRB*001","PrseDRB*002","PrseDRB*003","PrseDRB*004","PrseDRB*005","PrseDRB*006","PrseDRB*007","PrseDRB*008","PrseDRB*009","PrseDRB*010","PrseDRB*011","PrseDRB*012",
                                                "PrseDRB*013","PrseDRB*014","PrseDRB*015","PrseDRB*016","PrseDRB*017","PrseDRB*018","PrseDRB*019","PrseDRB*020","PrseDRB*021","PrseDRB*022","PrseDRB*023","PrseDRB*025",
                                                "PrseDRB*029","PrseDRB*031","PrseDRB*032","PrseDRB*039","PrseDRB*040","PrseDRB*044","PrseDRB*045","PrseDRB*054","PrseDRB*056","PrseDRB*058"))


dataALL <- NULL
input <- input_alleles_i

start <- 1
end <- 1000

cutoff_min <- 4
cutoff_max <- 46

i = start

set.seed(123)

for(i in start:end){
  data_subset <-  sample_n(input, 50)
  transposed_subset <- as.data.frame(t(data_subset[, (-1)]))
  transposed_subset$count <- rowSums(transposed_subset)
  filtered_subset <- transposed_subset %>%
    filter(count > cutoff_min) 
  filtered_subset <- filtered_subset %>%
    filter(count < cutoff_max) 
  Input_matrix <- filtered_subset[, c(1:ncol(filtered_subset)-1)]
  Co <- cooccur(mat = as.matrix(Input_matrix), 
                type = "spp_site", 
                thresh = FALSE, 
                spp_names = TRUE)
  effectsize_Co  <- effect.sizes(Co, standardized = TRUE, matrix = FALSE)     # standardized effect sizes
  names(effectsize_Co) <- c("sp1_name", "sp2_name", "effects")
  
  Complete_Co<-prob.table(Co) %>%
    filter(sp1_name =="N1" | sp1_name=="N2"  | sp1_name=="N3" | sp1_name=="N4" | sp1_name=="N5" | sp1_name=="N6" | sp1_name=="N7" | sp1_name=="N8" | sp1_name=="N9" | sp1_name=="N10" | sp1_name=="N11" | sp1_name=="N12" | sp1_name=="N13" | sp1_name=="Try" | sp1_name=="HdV" | sp1_name=="HpV" | sp1_name=="PbV" | sp1_name=="PcV") %>%
    left_join(effectsize_Co, by = c("sp1_name", "sp2_name"))
  
  Co_left <- Complete_Co %>%
    filter(p_lt < 0.1)
  Co_right <- Complete_Co %>%
    filter(p_gt < 0.1)
  
  Final_Co=rbind(Co_left,Co_right)
  
  Final_Co$obs.v.exp <- round((Final_Co$exp_cooccur-Final_Co$obs_cooccur)*-1, digits = 3)
  Final_Co$NR_IDs <- ncol(Input_matrix)
  Final_Co$round <- i
  dataALL=rbind(dataALL,Final_Co)
  
  i  = as.numeric (i + 1)
  
}


# Now filter only significant results

Sign_left <- dataALL %>%
  filter(p_lt < 0.05)
Sign_right <- dataALL %>%
  filter(p_gt < 0.05)

dataFilt=rbind(Sign_left,Sign_right)

res_Co_alleles <- dataFilt

#write.csv(res_Co_alleles, "res_alleles_i1000.csv")
#res_Co_alleles <- read.csv("res_alleles_i1000.csv")


# Draw heatmap of the results

heatmap.dat.allele <- res_Co_alleles %>%
  mutate(Asterisks = ifelse(p_gt  <= 0.001 |p_lt  <= 0.001, "***",
                            ifelse(p_gt <= 0.01 | p_lt  <= 0.01, "**",
                                   ifelse(p_gt  <= 0.05 | p_lt  <= 0.05,  "*", NA))))

head(heatmap.dat.allele)

table(heatmap.dat.allele$effects)

ggplot(heatmap.dat.allele, aes(sp2_name, sp1_name, fill=effects)) + 
  geom_tile()+
  #scale_fill_gradient2(high="#D53E4F", mid="white", low="#3288BD", limits = c(-0.1, 0.9))+
  scale_fill_gradient2(high="#D53E4F", mid="white", low="#3288BD")+
  facet_wrap(~sp1_name) + 
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(angle = 90))+
  geom_text(data = heatmap.dat.allele, aes(x = sp2_name, y = sp1_name, label = Asterisks), size = 4)



# Barplot that shows number of associations found across 

tbl <- data.frame(table(res_Co_alleles$round))

ggplot(tbl, aes(factor(Var1), Freq, fill = Var1)) +     
  geom_bar(stat="identity", position = "dodge")+
  ylab("Nr of associations")+
  xlab("Co-occurrence round")+
  theme_classic(base_size=14)+
  theme(legend.position = "NONE")

mean(tbl$Freq)
median(tbl$Freq)

table(tbl$Freq)


# How stable are the associations, how often are they found?

head(res_Co_alleles)

res_Co_alleles$association <- paste(res_Co_alleles$sp1_name, res_Co_alleles$sp2_name)


tbl2 <- data.frame(table(res_Co_alleles$association))

# exclude spurious associations (present in < 5% of co-occurrence rounds)

tbl2 %>%
  filter(Freq > 50)

ggplot(tbl2, aes(factor(Var1), Freq, fill = Var1)) +     
  geom_bar(stat="identity", position = "dodge")+
  ylab("Nr of times the association was found in 1000 repeated rounds")+
  xlab("Association pair")+
  theme_classic(base_size=14)+
  theme(legend.position = "NONE") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#write.csv(data.frame(tbl2), "robustness_associations_I1000.csv")




############################# Forest fragments in agricultural matrix ##############################

input_alleles_a <- subset(prse_meta_A, select=c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11","N12","N13","Try","HdV","HpV","PbV","PcV",
                                                "PrseDRB*001","PrseDRB*002","PrseDRB*003","PrseDRB*004","PrseDRB*005","PrseDRB*006","PrseDRB*007","PrseDRB*008","PrseDRB*009","PrseDRB*010","PrseDRB*011","PrseDRB*012",
                                                "PrseDRB*013","PrseDRB*014","PrseDRB*015","PrseDRB*016","PrseDRB*017","PrseDRB*018","PrseDRB*019","PrseDRB*020","PrseDRB*021","PrseDRB*022","PrseDRB*023","PrseDRB*025",
                                                "PrseDRB*029","PrseDRB*031","PrseDRB*032","PrseDRB*039","PrseDRB*040","PrseDRB*044","PrseDRB*045","PrseDRB*054","PrseDRB*056","PrseDRB*058"))


dataALL <- NULL
input <- input_alleles_a

start <- 1
end <- 1000

cutoff_min <- 4
cutoff_max <- 46

i = start

set.seed(123)

for(i in start:end){
  data_subset <-  sample_n(input, 50)
  transposed_subset <- as.data.frame(t(data_subset[, (-1)]))
  transposed_subset$count <- rowSums(transposed_subset)
  filtered_subset <- transposed_subset %>%
    filter(count > cutoff_min) 
  filtered_subset <- filtered_subset %>%
    filter(count < cutoff_max) 
  Input_matrix <- filtered_subset[, c(1:ncol(filtered_subset)-1)]
  Co <- cooccur(mat = as.matrix(Input_matrix), 
                type = "spp_site", 
                thresh = FALSE, 
                spp_names = TRUE)
  effectsize_Co  <- effect.sizes(Co, standardized = TRUE, matrix = FALSE)     # standardized effect sizes
  names(effectsize_Co) <- c("sp1_name", "sp2_name", "effects")
  
  Complete_Co<-prob.table(Co) %>%
    filter(sp1_name =="N1" | sp1_name=="N2"  | sp1_name=="N3" | sp1_name=="N4" | sp1_name=="N5" | sp1_name=="N6" | sp1_name=="N7" | sp1_name=="N8" | sp1_name=="N9" | sp1_name=="N10" | sp1_name=="N11" | sp1_name=="N12" | sp1_name=="N13" | sp1_name=="Try" | sp1_name=="HdV" | sp1_name=="HpV" | sp1_name=="PbV" | sp1_name=="PcV") %>%
    left_join(effectsize_Co, by = c("sp1_name", "sp2_name"))
  
  Co_left <- Complete_Co %>%
    filter(p_lt < 0.1)
  Co_right <- Complete_Co %>%
    filter(p_gt < 0.1)
  
  Final_Co=rbind(Co_left,Co_right)
  
  Final_Co$obs.v.exp <- round((Final_Co$exp_cooccur-Final_Co$obs_cooccur)*-1, digits = 3)
  Final_Co$NR_IDs <- ncol(Input_matrix)
  Final_Co$round <- i
  dataALL=rbind(dataALL,Final_Co)
  
  i  = as.numeric (i + 1)
  
}




# Now filter only significant results

Sign_left <- dataALL %>%
  filter(p_lt < 0.05)
Sign_right <- dataALL %>%
  filter(p_gt < 0.05)

dataFilt=rbind(Sign_left,Sign_right)

res_Co_alleles <- dataFilt

#write.csv(res_Co_alleles, "res_alleles_a1000.csv")
#res_Co_alleles <- read.csv("res_alleles_a1000.csv")


# Draw heatmap of the results


heatmap.dat.allele <- res_Co_alleles %>%
  mutate(Asterisks = ifelse(p_gt  <= 0.001 |p_lt  <= 0.001, "***",
                            ifelse(p_gt <= 0.01 | p_lt  <= 0.01, "**",
                                   ifelse(p_gt  <= 0.05 | p_lt  <= 0.05,  "*", NA))))

head(heatmap.dat.allele)

table(heatmap.dat.allele$effects)

ggplot(heatmap.dat.allele, aes(sp2_name, sp1_name, fill=effects)) + 
  geom_tile()+
  #scale_fill_gradient2(high="#D53E4F", mid="white", low="#3288BD", limits = c(-0.1, 0.9))+
  scale_fill_gradient2(high="#D53E4F", mid="white", low="#3288BD")+
  facet_wrap(~sp1_name) + 
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(angle = 90))+
  geom_text(data = heatmap.dat.allele, aes(x = sp2_name, y = sp1_name, label = Asterisks), size = 4)



# Barplot that shows number of associations found across 


tbl <- data.frame(table(res_Co_alleles$round))

ggplot(tbl, aes(factor(Var1), Freq, fill = Var1)) +     
  geom_bar(stat="identity", position = "dodge")+
  ylab("Nr of associations")+
  xlab("Co-occurrence round")+
  theme_classic(base_size=14)+
  theme(legend.position = "NONE")

mean(tbl$Freq)



# How stable are the associations, how often are they found?

head(res_Co_alleles)

res_Co_alleles$association <- paste(res_Co_alleles$sp1_name, res_Co_alleles$sp2_name)


tbl2 <- data.frame(table(res_Co_alleles$association))

# exclude spurious associations (present in < 5% of co-occurrence rounds)

tbl2 %>%
  filter(Freq > 50)

ggplot(tbl2, aes(factor(Var1), Freq, fill = Var1)) +     
  geom_bar(stat="identity", position = "dodge")+
  ylab("Nr of times the association was found in 1000 repeated rounds")+
  xlab("Association pair")+
  theme_classic(base_size=14)+
  theme(legend.position = "NONE") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#write.csv(data.frame(tbl2), "robustness_associations_A1000.csv")





##############################################################################################
############################################ STs #############################################
##############################################################################################


####################################### Continuous Forest ####################################

input_STs_c <- subset(prse_meta_C, select=c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11","N12","N13","Try","HdV","HpV","PbV","PcV",
                                            "ST1","ST2","ST3","ST4","ST5","ST6","ST7","ST8","ST9","ST10","ST11","ST12","ST13","ST14","ST15"))


dataALL <- NULL
input <- input_STs_c

start <- 1
end <- 1000

cutoff_min <- 4
cutoff_max <- 46

cutoff_min_st <- 9
cutoff_max_st <- 41


i = start

set.seed(123)

for(i in start:end){
  data_subset <-  sample_n(input, 50)
  transposed_subset <- as.data.frame(t(data_subset[, (-1)]))
  transposed_subset$count <- rowSums(transposed_subset)
  
  pathogen_sub=transposed_subset[c(1:17),]
  
  filtered_p_subset <- pathogen_sub %>%
    filter(count > cutoff_min) 
  filtered_p_subset <- filtered_p_subset %>%
    filter(count < cutoff_max) 
  
  st_sub=transposed_subset[c(18:51),]
  
  filtered_st_sub <- st_sub %>%
    filter(count > cutoff_min_st) 
  filtered_st_sub <- filtered_st_sub %>%
    filter(count < cutoff_max_st) 
  
  filtered_subset=rbind(filtered_p_subset,filtered_st_sub)
  
  Input_matrix <- filtered_subset[, c(1:ncol(filtered_subset)-1)]
  Co <- cooccur(mat = as.matrix(Input_matrix), 
                type = "spp_site", 
                thresh = FALSE, 
                spp_names = TRUE)
  effectsize_Co  <- effect.sizes(Co, standardized = TRUE, matrix = FALSE)     # standardized effect sizes
  names(effectsize_Co) <- c("sp1_name", "sp2_name", "effects")
  
  Complete_Co<-prob.table(Co) %>%
    filter(sp1_name =="N1" | sp1_name=="N2"  | sp1_name=="N3" | sp1_name=="N4" | sp1_name=="N5" | sp1_name=="N6" | sp1_name=="N7" | sp1_name=="N8" | sp1_name=="N9" | sp1_name=="N10" | sp1_name=="N11" | sp1_name=="N12" | sp1_name=="N13" | sp1_name=="Try" | sp1_name=="HdV" | sp1_name=="HpV" | sp1_name=="PbV" | sp1_name=="PcV") %>%
    left_join(effectsize_Co, by = c("sp1_name", "sp2_name"))
  
  Co_left <- Complete_Co %>%
    filter(p_lt < 0.2)
  Co_right <- Complete_Co %>%
    filter(p_gt < 0.2)
  
  Final_Co=rbind(Co_left,Co_right)
  
  Final_Co$obs.v.exp <- round((Final_Co$exp_cooccur-Final_Co$obs_cooccur)*-1, digits = 3)
  Final_Co$NR_IDs <- ncol(Input_matrix)
  Final_Co$round <- i
  dataALL=rbind(dataALL,Final_Co)
  
  i  = as.numeric (i + 1)
  
}



# Now filter only significant results

Sign_left <- dataALL %>%
  filter(p_lt < 0.05)
Sign_right <- dataALL %>%
  filter(p_gt < 0.05)

dataFilt=rbind(Sign_left,Sign_right)

res_Co_sts_c <- dataFilt

#write.csv(res_Co_sts_c, "res_sts_C1000.csv")
#res_Co_sts_c <- read.csv("res_sts_C1000.csv")


# Draw heatmap of the results


heatmap.dat.st <- res_Co_sts_c %>%
  mutate(Asterisks = ifelse(p_gt  <= 0.001 |p_lt  <= 0.001, "***",
                            ifelse(p_gt <= 0.01 | p_lt  <= 0.01, "**",
                                   ifelse(p_gt  <= 0.05 | p_lt  <= 0.05,  "*", NA))))

head(heatmap.dat.st)

ggplot(heatmap.dat.st, aes(sp2_name, sp1_name, fill=effects)) + 
  geom_tile()+
  #scale_fill_gradient2(high="#D53E4F", mid="white", low="#3288BD", limits = c(-0.1, 0.9))+
  scale_fill_gradient2(high="#D53E4F", mid="white", low="#3288BD")+
  facet_wrap(~sp1_name) + 
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(angle = 90))+
  geom_text(data = heatmap.dat.st, aes(x = sp2_name, y = sp1_name, label = Asterisks), size = 4)



# Barplot that shows number of associations found across 

tbl <- data.frame(table(res_Co_sts_c$round))

ggplot(tbl, aes(factor(Var1), Freq, fill = Var1)) +     
  geom_bar(stat="identity", position = "dodge")+
  ylab("Nr of associations")+
  xlab("Co-occurrence round")+
  theme_classic(base_size=14)+
  theme(legend.position = "NONE")

mean(tbl$Freq)



# How stable are the associations, how often are they found?

res_Co_sts_c$association <- paste(res_Co_sts_c$sp1_name, res_Co_sts_c$sp2_name)

tbl2 <- data.frame(table(res_Co_sts_c$association))

table(tbl2$Freq, tbl2$Var1)

# exclude spurious associations (present in < 5% of co-occurrence rounds)

tbl2 %>%
  filter(Freq > 50)

ggplot(tbl2, aes(factor(Var1), Freq, fill = Var1)) +     
  geom_bar(stat="identity", position = "dodge")+
  ylab("Nr of times the association was found in 100 repeated rounds")+
  xlab("Association pair")+
  theme_classic(base_size=14)+
  theme(legend.position = "NONE") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




####################################### Forested islands ####################################

input_STs_i <- subset(prse_meta_I, select=c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11","N12","N13","Try","HdV","HpV","PbV","PcV",
                                            "ST1","ST2","ST3","ST4","ST5","ST6","ST7","ST8","ST9","ST10","ST11","ST12","ST13","ST14","ST15"))


dataALL <- NULL
input <- input_STs_i

start <- 1
end <- 1000

cutoff_min <- 4
cutoff_max <- 46

cutoff_min_st <- 9
cutoff_max_st <- 41


i = start

set.seed(123)

for(i in start:end){
  data_subset <-  sample_n(input, 50)
  transposed_subset <- as.data.frame(t(data_subset[, (-1)]))
  transposed_subset$count <- rowSums(transposed_subset)
  
  pathogen_sub=transposed_subset[c(1:17),]
  
  filtered_p_subset <- pathogen_sub %>%
    filter(count > cutoff_min) 
  filtered_p_subset <- filtered_p_subset %>%
    filter(count < cutoff_max) 
  
  st_sub=transposed_subset[c(18:32),]
  
  filtered_st_sub <- st_sub %>%
    filter(count > cutoff_min_st) 
  filtered_st_sub <- filtered_st_sub %>%
    filter(count < cutoff_max_st) 
  
  filtered_subset=rbind(filtered_p_subset,filtered_st_sub)
  
  Input_matrix <- filtered_subset[, c(1:ncol(filtered_subset)-1)]
  Co <- cooccur(mat = as.matrix(Input_matrix), 
                type = "spp_site", 
                thresh = FALSE, 
                spp_names = TRUE)
  effectsize_Co  <- effect.sizes(Co, standardized = TRUE, matrix = FALSE)     # standardized effect sizes
  names(effectsize_Co) <- c("sp1_name", "sp2_name", "effects")
  
  Complete_Co<-prob.table(Co) %>%
    filter(sp1_name =="N1" | sp1_name=="N2"  | sp1_name=="N3" | sp1_name=="N4" | sp1_name=="N5" | sp1_name=="N6" | sp1_name=="N7" | sp1_name=="N8" | sp1_name=="N9" | sp1_name=="N10" | sp1_name=="N11" | sp1_name=="N12" | sp1_name=="N13" | sp1_name=="Try" | sp1_name=="HdV" | sp1_name=="HpV" | sp1_name=="PbV" | sp1_name=="PcV") %>%
    left_join(effectsize_Co, by = c("sp1_name", "sp2_name"))
  
  Co_left <- Complete_Co %>%
    filter(p_lt < 0.1)
  Co_right <- Complete_Co %>%
    filter(p_gt < 0.1)
  
  Final_Co=rbind(Co_left,Co_right)
  
  Final_Co$obs.v.exp <- round((Final_Co$exp_cooccur-Final_Co$obs_cooccur)*-1, digits = 3)
  Final_Co$NR_IDs <- ncol(Input_matrix)
  Final_Co$round <- i
  dataALL=rbind(dataALL,Final_Co)
  
  i  = as.numeric (i + 1)
  
}



# Now filter only significant results

Sign_left <- dataALL %>%
  filter(p_lt < 0.05)
Sign_right <- dataALL %>%
  filter(p_gt < 0.05)

dataFilt=rbind(Sign_left,Sign_right)

res_Co_sts_i <- dataFilt

#write.csv(res_Co_sts_i, "res_sts_i1000.csv")
#res_Co_sts_i <- read.csv("res_sts_i1000.csv")


# Draw heatmap of the results


heatmap.dat.st <- res_Co_sts_i %>%
  mutate(Asterisks = ifelse(p_gt  <= 0.001 |p_lt  <= 0.001, "***",
                            ifelse(p_gt <= 0.01 | p_lt  <= 0.01, "**",
                                   ifelse(p_gt  <= 0.05 | p_lt  <= 0.05,  "*", NA))))

head(heatmap.dat.st)

ggplot(heatmap.dat.st, aes(sp2_name, sp1_name, fill=effects)) + 
  geom_tile()+
  #scale_fill_gradient2(high="#D53E4F", mid="white", low="#3288BD", limits = c(-0.1, 0.9))+
  scale_fill_gradient2(high="#D53E4F", mid="white", low="#3288BD")+
  facet_wrap(~sp1_name) + 
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(angle = 90))+
  geom_text(data = heatmap.dat.st, aes(x = sp2_name, y = sp1_name, label = Asterisks), size = 4)



# Barplot that shows number of associations found across 

tbl <- data.frame(table(res_Co_sts_i$round))

ggplot(tbl, aes(factor(Var1), Freq, fill = Var1)) +     
  geom_bar(stat="identity", position = "dodge")+
  ylab("Nr of associations")+
  xlab("Co-occurrence round")+
  theme_classic(base_size=14)+
  theme(legend.position = "NONE")

mean(tbl$Freq)



# How stable are the associations, how often are they found?

head(res_Co_sts_i)

res_Co_sts_i$association <- paste(res_Co_sts_i$sp1_name, res_Co_sts_i$sp2_name)


tbl2 <- data.frame(table(res_Co_sts_i$association))

table(tbl2$Freq, tbl2$Var1)

# exclude spurious associations (present in < 5% of co-occurrence rounds)

tbl2 %>%
  filter(Freq > 50)

ggplot(tbl2, aes(factor(Var1), Freq, fill = Var1)) +     
  geom_bar(stat="identity", position = "dodge")+
  ylab("Nr of times the association was found in 100 repeated rounds")+
  xlab("Association pair")+
  theme_classic(base_size=14)+
  theme(legend.position = "NONE") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





############################### Forest fragments in agricultural matrix ##############################

input_STs_a <- subset(prse_meta_A, select=c("N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11","N12","N13","Try","HdV","HpV","PbV","PcV",
                                            "ST1","ST2","ST3","ST4","ST5","ST6","ST7","ST8","ST9","ST10","ST11","ST12","ST13","ST14","ST15"))


dataALL <- NULL
input <- input_STs_a

start <- 1
end <- 1000

cutoff_min <- 4
cutoff_max <- 46

cutoff_min_st <- 9
cutoff_max_st <- 41


i = start

set.seed(123)

for(i in start:end){
  data_subset <-  sample_n(input, 50)
  transposed_subset <- as.data.frame(t(data_subset[, (-1)]))
  transposed_subset$count <- rowSums(transposed_subset)
  
  pathogen_sub=transposed_subset[c(1:17),]
  
  filtered_p_subset <- pathogen_sub %>%
    filter(count > cutoff_min) 
  filtered_p_subset <- filtered_p_subset %>%
    filter(count < cutoff_max) 
  
  st_sub=transposed_subset[c(18:32),]
  
  filtered_st_sub <- st_sub %>%
    filter(count > cutoff_min_st) 
  filtered_st_sub <- filtered_st_sub %>%
    filter(count < cutoff_max_st) 
  
  filtered_subset=rbind(filtered_p_subset,filtered_st_sub)
  
  Input_matrix <- filtered_subset[, c(1:ncol(filtered_subset)-1)]
  Co <- cooccur(mat = as.matrix(Input_matrix), 
                type = "spp_site", 
                thresh = FALSE, 
                spp_names = TRUE)
  effectsize_Co  <- effect.sizes(Co, standardized = TRUE, matrix = FALSE)     # standardized effect sizes
  names(effectsize_Co) <- c("sp1_name", "sp2_name", "effects")
  
  Complete_Co<-prob.table(Co) %>%
    filter(sp1_name =="N1" | sp1_name=="N2"  | sp1_name=="N3" | sp1_name=="N4" | sp1_name=="N5" | sp1_name=="N6" | sp1_name=="N7" | sp1_name=="N8" | sp1_name=="N9" | sp1_name=="N10" | sp1_name=="N11" | sp1_name=="N12" | sp1_name=="N13" | sp1_name=="Try" | sp1_name=="HdV" | sp1_name=="HpV" | sp1_name=="PbV" | sp1_name=="PcV") %>%
    left_join(effectsize_Co, by = c("sp1_name", "sp2_name"))
  
  Co_left <- Complete_Co %>%
    filter(p_lt < 0.1)
  Co_right <- Complete_Co %>%
    filter(p_gt < 0.1)
  
  Final_Co=rbind(Co_left,Co_right)
  
  Final_Co$obs.v.exp <- round((Final_Co$exp_cooccur-Final_Co$obs_cooccur)*-1, digits = 3)
  Final_Co$NR_IDs <- ncol(Input_matrix)
  Final_Co$round <- i
  dataALL=rbind(dataALL,Final_Co)
  
  i  = as.numeric (i + 1)
  
}



# Now filter only significant results

Sign_left <- dataALL %>%
  filter(p_lt < 0.05)
Sign_right <- dataALL %>%
  filter(p_gt < 0.05)

dataFilt=rbind(Sign_left,Sign_right)

res_Co_sts_a <- dataFilt

#write.csv(res_Co_sts_a, "res_sts_a1000.csv")
res_Co_sts_a <- read.csv("res_sts_a1000.csv")



# Draw heatmap of the results


heatmap.dat.st <- res_Co_sts_a %>%
  mutate(Asterisks = ifelse(p_gt  <= 0.001 |p_lt  <= 0.001, "***",
                            ifelse(p_gt <= 0.01 | p_lt  <= 0.01, "**",
                                   ifelse(p_gt  <= 0.05 | p_lt  <= 0.05,  "*", NA))))

head(heatmap.dat.st)

ggplot(heatmap.dat.st, aes(sp2_name, sp1_name, fill=effects)) + 
  geom_tile()+
  #scale_fill_gradient2(high="#D53E4F", mid="white", low="#3288BD", limits = c(-0.1, 0.9))+
  scale_fill_gradient2(high="#D53E4F", mid="white", low="#3288BD")+
  facet_wrap(~sp1_name) + 
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(angle = 90))+
  geom_text(data = heatmap.dat.st, aes(x = sp2_name, y = sp1_name, label = Asterisks), size = 4)



# Barplot that shows number of associations found across 


tbl <- data.frame(table(res_Co_sts_a$round))

ggplot(tbl, aes(factor(Var1), Freq, fill = Var1)) +     
  geom_bar(stat="identity", position = "dodge")+
  ylab("Nr of associations")+
  xlab("Co-occurrence round")+
  theme_classic(base_size=14)+
  theme(legend.position = "NONE")

mean(tbl$Freq)
median(tbl$Freq)

table(tbl$Freq)



# How stable are the associations, how often are they found?

head(res_Co_sts_a)

res_Co_sts_a$association <- paste(res_Co_sts_a$sp1_name, res_Co_sts_a$sp2_name)

tbl2 <- data.frame(table(res_Co_sts_a$association))

# exclude spurious associations (present in < 5% of co-occurrence rounds)

tbl2 %>%
  filter(Freq > 50)

ggplot(tbl2, aes(factor(Var1), Freq, fill = Var1)) +     
  geom_bar(stat="identity", position = "dodge")+
  ylab("Nr of times the association was found in 100 repeated rounds")+
  xlab("Association pair")+
  theme_classic(base_size=14)+
  theme(legend.position = "NONE") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




##########################################################################################################################################################################################################
### MHC and pathogen dissimilarity
##########################################################################################################################################################################################################

### Explore pathogen - MHC coevolution ###

# load library

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)



# first subset pathogens, ST and alleles

pathogens <- subset(meta, select=c(N1:N13,Try:PcV))
str(pathogens)

supertypes <- subset(meta, select=c(ST1:ST15))
str(supertypes)

alleles <- subset(meta, select=c(PrseDRB.001:PrseDRB.079))
str(alleles)




####################################################################################################################################
#################################################### Mantel tests per landscape ####################################################
####################################################################################################################################


##############################################################################################
####################################### Continuous Forest ####################################
##############################################################################################

# Subset C
meta_C <- subset(meta,landscape=="continuous.forest")

alleles_C <- subset(meta_C, select=c(PrseDRB.001:PrseDRB.079))
pathogens_C <- subset(meta_C, select=c(N1:N13,Try:PcV))
supertypes_C <- subset(meta_C, select=c(ST1:ST15))




############################ pathogens & STs #################################

# Order data frames 
pathogens_C <-pathogens_C[order(rownames(pathogens_C)),]
supertypes_C <- supertypes_C[order(rownames(supertypes_C)),]

# Distance matrices
pathogens.dist_C <- vegdist(pathogens_C, method='jaccard', binary=F)
#pathogens.dist_C1<-pathogens.dist_C/sd(pathogens.dist_C)

supertypes.dist_C <- vegdist(supertypes_C, method='jaccard', binary=F)
#supertypes.dist_C1<-supertypes.dist_C/sd(supertypes.dist_C)

# Check that order of distance matrices match
labels(pathogens.dist_C) == labels(supertypes.dist_C)
#labels(pathogens.dist_C1) == labels(supertypes.dist_C1)

# Mantel test
set.seed(123)
vegan::mantel(pathogens.dist_C,supertypes.dist_C, permutations = 9999)
#vegan::mantel(pathogens.dist_C1,supertypes.dist_C1, permutations = 9999)

# Correlation figure
corrdat_C <- data.frame(supertypes.vector=as.vector(supertypes.dist_C),pathogens.vector=as.vector(pathogens.dist_C))

corr.C1 <- ggplot(corrdat_C, aes(x=supertypes.vector, y=pathogens.vector)) + 
  geom_point(colour='blue', size=1) + 
  geom_smooth(method=lm)+theme_bw()+
  ggtitle("Continuous forest")+
  labs(x="MHC ST dissimilarity", y="Pathogen dissimilarity")+
  geom_text(x=0.3, y=0.95, label="r = 0.03, p = 0.306", fontface="italic")+
  ylim(0,1) +
  xlim(0,1)

# Correlation figure normalized matrices
#corrdat1 <- data.frame(supertypes.vector=as.vector(supertypes.dist_C1),pathogens.vector=as.vector(pathogens.dist_C1))
#ggplot(corrdat1, aes(x=supertypes.vector, y=pathogens.vector))



########################## pathogens & alleles ###############################

### Correlate ecological and genetic networks
alleles_C <- alleles_C[order(rownames(alleles_C)),]

# Distance matrices
alleles.dist_C <- vegdist(alleles_C, method='jaccard', binary=F)

labels(pathogens.dist_C) == labels(alleles.dist_C)

# Mantel test
set.seed(123)
vegan::mantel(pathogens.dist_C, alleles.dist_C, permutations = 9999)

# Correlation figure
corrdat_C <- data.frame(alleles.vector=as.vector(alleles.dist_C),pathogens.vector=as.vector(pathogens.dist_C))

corr.C <- ggplot(corrdat_C, aes(x=alleles.vector, y=pathogens.vector)) + 
  geom_point(colour='blue', size=1) + 
  geom_smooth(method=lm)+theme_bw() +
  ggtitle("Continuous forest")+
  labs(x="MHC allele dissimilarity", y="Pathogen dissimilarity")+
  geom_text(x=0.3, y=0.95, label="r = 0.13, p = 0.003", fontface="italic")+
  ylim(0,1) +
  xlim(0,1)




#############################################################################################
####################################### Forested Islands ####################################
#############################################################################################

# Subset I
meta_I <- subset(meta,landscape=="island")

alleles_I <- subset(meta_I, select=c(PrseDRB.001:PrseDRB.079))
pathogens_I <- subset(meta_I, select=c(N1:N13,Try:PcV))
supertypes_I <- subset(meta_I, select=c(ST1:ST15))



############################ pathogens & STs #################################

# Order data frames 
pathogens_I <-pathogens_I[order(rownames(pathogens_I)),]
supertypes_I <- supertypes_I[order(rownames(supertypes_I)),]

# Distance matrices
pathogens.dist_I <- vegdist(pathogens_I, method='jaccard', binary=F)
#pathogens.dist_I1<-pathogens.dist_I/sd(pathogens.dist_I)

supertypes.dist_I <- vegdist(supertypes_I, method='jaccard', binary=F)
#supertypes.dist_I1<-supertypes.dist_I/sd(supertypes.dist_I)

# Check that order of distance matrices match
labels(pathogens.dist_I) == labels(supertypes.dist_I)
#labels(pathogens.dist_I1) == labels(supertypes.dist_I1)

# Mantel test
set.seed(123)
vegan::mantel(pathogens.dist_I,supertypes.dist_I, permutations = 9999)
#vegan::mantel(pathogens.dist_I1,supertypes.dist_I1, permutations = 9999)

# Correlation figure
corrdat_I <- data.frame(supertypes.vector=as.vector(supertypes.dist_I),pathogens.vector=as.vector(pathogens.dist_I))

corr.I1 <- ggplot(corrdat_I, aes(x=supertypes.vector, y=pathogens.vector)) + 
  geom_point(colour='blue', size=1) +
  geom_smooth(method=lm)+theme_bw()+
  ggtitle("Forested islands")+
  labs(x="MHC ST dissimilarity", y="Pathogen dissimilarity")+
  geom_text(x=0.3, y=0.95, label="r = 0.004, p = 0.463", fontface="italic")+
  ylim(0,1) +
  xlim(0,1)



########################## pathogens & alleles ###############################

### Correlate ecological and genetic networks
alleles_I <- alleles_I[order(rownames(alleles_I)),]

# Distance matrices
alleles.dist_I <- vegdist(alleles_I, method='jaccard', binary=F)

labels(pathogens.dist_I) == labels(alleles.dist_I)

# Mantel test
set.seed(123)
vegan::mantel(pathogens.dist_I, alleles.dist_I, permutations = 9999)

# Correlation figure
corrdat_I <- data.frame(alleles.vector=as.vector(alleles.dist_I),pathogens.vector=as.vector(pathogens.dist_I))

corr.I <- ggplot(corrdat_I, aes(x=alleles.vector, y=pathogens.vector)) + 
  geom_point(colour='blue', size=1) + 
  geom_smooth(method=lm)+theme_bw()+
  ggtitle("Forested islands")+
  labs(x="MHC allele dissimilarity", y="Pathogen dissimilarity")+
  geom_text(x=0.3, y=0.95, label="r = -0.03, p = 0.739", fontface="italic")+
  ylim(0,1) +
  xlim(0,1)




##############################################################################################
########################## Forest fragments in agricultural matrix ###########################
##############################################################################################

# Subset A
meta_A <- subset(meta,landscape=="forest.fragment")

alleles_A <- subset(meta_A, select=c(PrseDRB.001:PrseDRB.079))
supertypes_A <- subset(meta_A, select=c(ST1:ST15))
pathogens_A <- subset(meta_A, select=c(N1:N13,Try:PcV))



############################ pathogens & ST #################################

# Order data frames 
pathogens_A <-pathogens_A[order(rownames(pathogens_A)),]
supertypes_A <- supertypes_A[order(rownames(supertypes_A)),]

# Distance matrices
pathogens.dist_A <- vegdist(pathogens_A, method='jaccard', binary=F)
#pathogens.dist_A1<-pathogens.dist_A/sd(pathogens.dist_A)

supertypes.dist_A <- vegdist(supertypes_A, method='jaccard', binary=F)
#supertypes.dist_A1<-supertypes.dist_A/sd(supertypes.dist_A)

# Check that order of distance matrices match
labels(pathogens.dist_A) == labels(supertypes.dist_A)
#labels(pathogens.dist_A1) == labels(supertypes.dist_A1)

# Mantel test
set.seed(123)
vegan::mantel(pathogens.dist_A,supertypes.dist_A, permutations = 9999)
#vegan::mantel(pathogens.dist_A1,supertypes.dist_A1, permutations = 9999)

# Correlation figure
corrdat_A <- data.frame(supertypes.vector=as.vector(supertypes.dist_A),pathogens.vector=as.vector(pathogens.dist_A))

corr.A1 <- ggplot(corrdat_A, aes(x=supertypes.vector, y=pathogens.vector)) + 
  geom_point(colour='blue', size=1) + 
  geom_smooth(method=lm)+theme_bw() + 
  ggtitle("Forest in agricultural matrix")+
  labs(x="MHC ST dissimilarity", y="Pathogen dissimilarity")+
  geom_text(x=0.3, y=0.95, label="r = 0.03, p = 0.270", fontface="italic")+
  ylim(0,1) +
  xlim(0,1)



########################## pathogens & alleles ###############################

### Correlate ecological and genetic networks
alleles_A <- alleles_A[order(rownames(alleles_A)),]

# Distance matrices
alleles.dist_A <- vegdist(alleles_A, method='jaccard', binary=F)

matrix_a <- as.matrix(alleles.dist_A)
#write.csv(matrix_a, "matrix_a.csv")

labels(pathogens.dist_A) == labels(alleles.dist_A)

# Mantel test
set.seed(123)
vegan::mantel(pathogens.dist_A, alleles.dist_A, permutations = 9999)

# Correlation figure
corrdat_A <- data.frame(alleles.vector=as.vector(alleles.dist_A),pathogens.vector=as.vector(pathogens.dist_A))

corr.A <- ggplot(corrdat_A, aes(x=alleles.vector, y=pathogens.vector)) + 
  geom_point(colour='blue', size=1) + 
  geom_smooth(method=lm)+theme_bw() +
  ggtitle("Forest in agricultural matrix")+
  labs(x="MHC allele dissimilarity", y="Pathogen dissimilarity")+
  geom_text(x=0.3, y=0.95, label="r = 0.04, p = 0.134", fontface="italic")+
  ylim(0,1) +
  xlim(0,1)


# Plot comparison of correlations

#####################
### Supp Figure 5 ###
#####################

correlation.plots <- ggarrange(corr.C,corr.I,corr.A, ncol=3)
correlation.plots1 <- ggarrange(corr.C1,corr.I1,corr.A1, ncol=3)

ggarrange(correlation.plots, correlation.plots1, nrow=2, labels = c("a)", "b)"))

