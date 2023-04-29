setwd("~/Library/CloudStorage/GoogleDrive-nn33@hawaii.edu/Shared drives/Rhizobia of Hawai'i/Analyses/R analyses")
library(phyloseq)
library(agricolae)
library(ggplot2)
library(Rmisc)
library(ggpubr)

##===================================================================================================
## Read in 16S data
##===================================================================================================

#Read in metadata file and merge other metrics with the metadata file
metadata_table <- read.csv(file="metadata-site.csv", header=TRUE) %>% 
  as_tibble()

#Read in sample x genus table
otu_table <- read.csv(file="table-genus-proportion.csv", header=TRUE, row.names=1) %>%
  t() %>% 
  as_tibble()

#Format data tables for phyloseq
otu <- otu_table(otu_table, taxa_are_rows=FALSE)
metadata <- sample_data(metadata_table)

#Merge data tables into a phyloseq object
merged_tables <- phyloseq(otu, metadata)

##===================================================================================================
##  STATISTICS: COMPARISONS BETWEEN genus and site
##===================================================================================================
#PERMANOVA comparisons among community members
#Make distance matrices
braydist <- phyloseq::distance(merged_tables, method="bray")
eucldist <- phyloseq::distance(merged_tables, method="euclidean")

#Run PERMANOVA for 16S data
set.seed(2023)
adonis2(braydist ~ LandUse, data=metadata_table, permutations=999) #NS

set.seed(2023)
adonis2(eucldist ~ LandUse, data=metadata_table, permutations=999) #NS

##===================================================================================================
## STATISTICS: Comparison between effective inocula
##===================================================================================================

#Read in data
rhizobia<-read.csv("Rhizobia_inoculum.csv", header=TRUE)

##Testing ANOVA assumptions; perform ANOVA
# Using the Shapiro Test -- if p-value not significant = data is normal = OK for ANOVA
shapiro.test(rhizobia$weight)#NS
shapiro.test(rhizobia$height_cm)#NS
shapiro.test(rhizobia$volume)#NS
shapiro.test(rhizobia$nodule)#NS

# Test for homogeneity of variances -- if p-value not significant = data is homogeneous = OK for ANOVA
bartlett.test(rhizobia$weight~rhizobia$treatment)#S
bartlett.test(rhizobia$height_cm~rhizobia$treatment)#NS
bartlett.test(rhizobia$volume~rhizobia$treatment)#NS
bartlett.test(rhizobia$nodule~rhizobia$treatment)#S

#-----------
#Plant weight (ANOVA assumptions failed, using a non-parametric test)
#Kruskal-Wallis test
kruskal_rhizobiaw <- kruskal(rhizobia$weight, rhizobia$treatment, group=TRUE, p.adj="BH")$groups

#create new data frame
rhizobiaw1<-data.frame(rhizobia$treatment,rhizobia$weight);rhizobiaw1

#get standard error
colnames(rhizobiaw1)<-c("treatment","weight")
rhizobiaw1SE<-summarySE(rhizobiaw1,measurevar="weight",groupvars="treatment");rhizobiaw1SE

#add significant letters
letters1<-c("a","b","c")
rhizobiaw1SE1<-cbind(rhizobiaw1SE,letters1);rhizobiaw1SE1

#Plot
(plot1 <- ggplot(rhizobiaw1SE,aes(x=treatment,y=weight)) +
  geom_bar(stat="identity",fill="grey50", width=0.55) +
  geom_errorbar(aes(ymin=weight-ci,ymax=weight+ci),width=0.1) +
  geom_text(aes(label=letters1,y=weight+(ci+1.5)))+
  labs(x="", y="biomass (g)") +
  theme_light()
)

#-----------
#Plant weight
kruskal_rhizobiah <- kruskal(rhizobia$height_mm, rhizobia$treatment, group=TRUE, p.adj="BH")$groups

#create new data frame
rhizobiah1<-data.frame(rhizobia$treatment,rhizobia$height_mm);rhizobiah1

#get standard error
colnames(rhizobiah1)<-c("treatment","height_mm")
rhizobiah1SE<-summarySE(rhizobiah1,measurevar="height_mm",groupvars="treatment");rhizobiah1SE

#add significant letters
letters2<-c("a","a","a")
rhizobiah1SE2<-cbind(rhizobiah1SE,letters2);rhizobiah1SE2

#Plot
(plot2 <- ggplot(rhizobiah1SE,aes(x=treatment,y=height_mm))+
  geom_bar(stat="identity",fill="grey50", width=0.55)+
  geom_errorbar(aes(ymin=height_mm-ci,ymax=height_mm+ci),width=0.1)+
  geom_text(aes(label=letters2,y=height_mm+(ci+30)))+
    labs(x="", y="stem height (mm)") +
    theme_light()
)

#-----------
#Plant volume
kruskal_rhizobiav <- kruskal(rhizobia$volume, rhizobia$treatment, group=TRUE, p.adj="BH")$groups

#create new data frame
rhizobiav1<-data.frame(rhizobia$treatment,rhizobia$volume);rhizobiav1

#get standard error
colnames(rhizobiav1)<-c("treatment","volume")
rhizobiav1SE<-summarySE(rhizobiav1,measurevar="volume",groupvars="treatment");rhizobiav1SE

#add significant letters
letters3<-c("a","a","a")
rhizobiav1SE3<-cbind(rhizobiav1SE,letters3);rhizobiav1SE3

#Plot
(plot3 <- ggplot(rhizobiav1SE,aes(x=treatment,y=volume)) +
  geom_bar(stat="identity",fill="grey50", width=0.55) +
  geom_errorbar(aes(ymin=volume-ci,ymax=volume+ci),width=0.1) +
  geom_text(aes(label=letters3, y=volume+(ci+700))) +
  labs(x="", y=bquote('stem volume '(mm^3))) +
  theme_light()
)

#-----------
#Plant nodules
kruskal_rhizobiah <- kruskal(rhizobia$nodules, rhizobia$treatment, group=TRUE, p.adj="BH")$groups

#create new data frame
rhizobian1<-data.frame(rhizobia$treatment,rhizobia$nodules);rhizobian1

#get standard error
colnames(rhizobian1)<-c("treatment","nodules")
rhizobian1SE<-summarySE(rhizobian1,measurevar="nodules",groupvars="treatment");rhizobian1SE

#add significant letters
letters4<-c("a","b","b")
rhizobian1SE4<-cbind(rhizobian1SE,letters4);rhizobian1SE4

#Plot
(plot4 <- ggplot(rhizobian1SE,aes(x=treatment,y=nodules))+
  geom_bar(stat="identity",fill="grey50", width=0.55)+
  geom_errorbar(aes(ymin=nodules-ci,ymax=nodules+ci),width=0.1)+
  geom_text(aes(label=letters4,y=nodules+(ci+1.5)))+
  labs(x="", y="number of nodules") +
  theme_light()
)

#Plot together
figure<-ggarrange(plot1,plot2,plot3,plot4,labels=c("A","B","C","D"),ncol=2,nrow=2);figure

#Save the plot
ggsave("Figure3-rhizobia-inoculation.pdf", device="pdf", width=7, height=5.5)
