setwd("~/Library/CloudStorage/GoogleDrive-nn33@hawaii.edu/Shared drives/Rhizobia of Hawai'i/Analyses/R analyses")

library(ggmap)
# library(usmap) #This package has nice maps but coordinates have to be connected
#   plot_usmap(include = c("HI"))

#read in data
data <- read.csv("site-coordinates.csv")

#fetch map
mapWorld <-  #world map, x-coordinate=longitude, y-coordinate=latitude
  borders(database="world", colour="grey80", fill="grey95")

#map simple coordinate points, zoom into Oahu
(bubbleMap <-
  ggplot(data) +
  mapWorld +
  geom_point(aes(x=Long, y=Lat, size=samples), color="#7570b3", alpha=0.8) +
  coord_quickmap(xlim=c(-158.3, -157.6), ylim=c(21.25, 21.725)) +
  labs(x="Longitude", y="Latitude") +
  theme_bw()
  #+ theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14))
)

ggsave("Figure1-Map.pdf", device="pdf", width=6.4, height=3)

# #map pie plot with library(scatterpie). I have not changed these parameters to reflect data for this project yet.
# (pie_plot <-
#   ggplot(cluster_sum_wide) +
#   mapWorld +
#   geom_scatterpie(aes(x=mean_long, y=mean_lat, r=total), alpha=0.8, data=cluster_sum_wide, cols=c("complete", "incomplete")) +
#   geom_scatterpie_legend(radius = c(1,10,20), x=-160, y=-55) +
#   coord_quickmap()
# )
