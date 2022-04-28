## Make a map of TAM_05 with source trees highlighted and Gentry #46088 featured

library(tidyverse)
library(ggplot2)
library(reshape)
library(ggstar)

# Load plot data from "specimen cleaning and tree matching.R"
tam <- Live_trees %>%
  filter(Plot == "TAM_05")
tam$source <- "false"

tam.source <- source_trees %>%
  filter(Plot == "TAM_05")
tam.source$source <- "true"

tam.all <- full_join(tam, tam.source)

tam.hymen <- tam.source %>%
  filter(CollectionNumber == 46088)

# Make the map
cols <- c("black","#377EB8")
map <- ggplot(tam.all, aes(y = Y, x = X))+
  labs(title = "Tambopata plot three", x = "Meters", y = "Meters", color = "") +
  scale_x_continuous(expand = c(0,0), limits = c(0,100), labels=c(seq(0,100,20)), breaks=c(seq(0,100,20)))+ #assigns the x axis limits of 0 and 100, with 20 meter lines
  scale_y_continuous(expand = c(0,0), limits = c(0,100), labels=c(seq(0,100,20)), breaks=c(seq(0,100,20))) + #assigns the y axis limits of 0 and 100
  theme_classic() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = c(1.3,0.95),
        plot.title = element_text(size=22),
        text = element_text(size=20)) +
  geom_point(shape=16, aes(color=source, size = DBH1)) +
  scale_color_manual(values = cols, labels = c("Living trees", "Living source trees with \n herbarium specimens")) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  scale_size(guide = "none") +
  geom_star(data = tam.hymen, show.legend=FALSE, fill="#E41A1C", size=7) +
  coord_fixed()
map
ggsave("TAM_05_map.png")
