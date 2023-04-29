setwd("~/Library/CloudStorage/GoogleDrive-nn33@hawaii.edu/Shared drives/Rhizobia of Hawai'i/Analyses/R analyses")
library(tidyverse)
library(ggplot2)

#Read in data
rhizo <- read.csv("table-genus-proportion.csv")
rhizo[rhizo==0]=NA #convert 0's to NAs so that they don't plot

#graph
(rhizobubble <- rhizo %>% 
# rownames_to_column(var="Genus") %>% #use if rownames is not a column
  gather(location, abundance, -Genus) %>% #make tidy
  ggplot(aes(location, Genus, y=reorder(Genus,desc(Genus)))) +
  geom_point(aes(size=abundance), color="#555555", fill="#d65d5f", shape=21, alpha=0.8) +
  labs(x="", y="", size="Abundance (%)") +
  scale_x_discrete(position="bottom", 
                   limit=c("A.WaimeaRidge", "B.NorthOahu2", "C.Poamoho", "D.Lyon", "E.Waimanalo", "F.NorthOahu1", "G.West_Oahu", "H.WaimeaValley", "I.WaimeaBG", "J.UHM", "K.KCBG", "L.WR"), 
                   labels = c("Waimea Ridge", "North Oahu Farm 2", "Poamoho Research Station", "Lyon Arboretum", "Waimanalo Research Station", "North Oahu Farm 1", "West Oahu Farm", "Waimea Valley", "Waimea Botanical Garden", "University of Hawai'i campus", "Koko Crater Botancical Garden", "Wa'ahila Ridge")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, color="#111111", size=11)) +
  theme(axis.text.y = element_text(face="italic", color="#111111", size=11))
)

ggsave("Figure2-2-rhizobubble.pdf", device="pdf", width=7.3, height=7.3)
 #color="#555555", fill="#d65d5f"
#nice blue color="#50587e", fill="#436fb6"

---------------------------

#hierachical cluster
#Read in data
rhizo2 <- read.csv("table-genus-proportion1.csv", row.names="Genus") %>% 
  t()

#Scale the data to make them comparable
rhizo2 <- scale(rhizo2)

# Make dissimilarity matrix
set.seed(2023)
rhizodist <- dist(rhizo2, method="euclidean")

# Hierarchical clustering using Complete Linkage
rhizocluster <- hclust(rhizodist, method ="complete")

# Plot the obtained dendrogram
plot(rhizocluster, hang=-1)

ggsave("Figure2-dendrogram.pdf", device="pdf", width=4.25, height=2.8)
