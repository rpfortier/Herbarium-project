### Data analysis for individual tree trait changes through time ###

#set working directory

# Load libraries
library(dplyr)
library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(MuMIn)
library(stringr)
library(ggpubr)
library(PairedData)
library(raster)
library(terra)
library(gridExtra)
library(zoo)
library(rstatix) 

###summarize and clean trait data
leaves_modern <- read.csv("leaf_traits_data.csv")

#remove stratum 2 (middle of crown)
leaves_modern <- leaves_modern[which(leaves_modern$Stratum != 2),]

#calculate elw
leaves_modern$elw_fresh <- leaves_modern$circle_perimeter_fresh/pi
leaves_modern$elw <- leaves_modern$circl_perim/pi

#average the three dry leaf thickness measurements
leaves_modern$thickness <- apply(leaves_modern[,c("thickness1", "thickness2", "thickness3")], 1, mean)
leaves_modern$thickness <- leaves_modern$thickness/1000 #convert from micrometers to millimeters
leaves_modern$thickness_fresh <- leaves_modern$thickness_fresh/1000 #convert from micrometers to millimeters

#convert mass to grams
leaves_modern$mass <- leaves_modern$mass/1000 #convert from mg to g

#calculate leaf mass per area (LMA)
leaves_modern$lma_fresh <- (leaves_modern$mass)/leaves_modern$leaf_area_fresh*1000 
leaves_modern$lma <- (leaves_modern$mass)/leaves_modern$area*1000 

#add the species for each leaf
trees2023 <- read.csv("trees_collected.csv")

x <- left_join(leaves_modern, trees2023, by = c("Plot", "Tag"))
leaves_modern$species <- x$Species
rm(x)

#summarize stomatal traits
#start with size
stomata_dry <- read_excel("C:/Users/rfortier/Dropbox/Peru - Individual Tree Project/Stomatal peels/Stomatal analysis.xlsx")
stomata_dry <- stomata_dry %>%
  group_by(Image_name) %>%
  summarise("Pore_length" = mean(Pore_length),
            "Guard_cell_width" = mean(Guard_cell_width))

stomata_fresh <- read_excel("C:/Users/rfortier/Dropbox/Peru - Individual Tree Project/Stomatal peels/Stomatal analysis.xlsx", 
                                sheet = "FRESH")
stomata_fresh$Species <- NULL
stomata_fresh$Image_name <- sub(".tif", ".jpg", stomata_fresh$Image_name)
#summarize size traits
stomata_fresh <- stomata_fresh %>%
  group_by(Image_name) %>%
  summarise("Pore_length_fresh" = mean(Pore_length_fresh),
            "Guard_cell_width_fresh" = mean(Guard_cell_width_fresh))

#combine dry and fresh stomatal data
stomata_size <- left_join(stomata_dry, stomata_fresh, by = "Image_name")

#parse out plot, tag, strata, and leaf info
stomata_size$Image_name <- substr(stomata_size$Image_name,1,nchar(stomata_size$Image_name)-4)
stomata_size$leaf <- str_sub(stomata_size$Image_name, start = -1)
stomata_size$leaf <- toupper(stomata_size$leaf)
stomata_size$Stratum <- str_sub(stomata_size$Image_name, start = -2, end = -2)
stomata_size$Plot <- ifelse(startsWith(stomata_size$Image_name, "S"),str_sub(stomata_size$Image_name, start = 1, end = 3), str_sub(stomata_size$Image_name, start = 1, end = 6))
stomata_size$Tag <- sub("\\-.*", "", stomata_size$Image_name)
stomata_size$Tag <- ifelse(startsWith(stomata_size$Image_name, "S"),str_sub(stomata_size$Tag, start = 5),str_sub(stomata_size$Tag, start = 8))
stomata_size$Stratum <- as.numeric(stomata_size$Stratum)
stomata_size$Plot <- ifelse(stomata_size$Plot == "SEN", "Sendero", stomata_size$Plot)
stomata_size$Image_name <- NULL

#now do density
stomata_dry_dens <- read_excel("C:/Users/rfortier/Dropbox/Peru - Individual Tree Project/Stomatal peels/Stomatal analysis.xlsx", 
                               sheet = "DRY_Density")
stomata_fresh_dens <- read_excel("C:/Users/rfortier/Dropbox/Peru - Individual Tree Project/Stomatal peels/Stomatal analysis.xlsx", 
                               sheet = "FRESH_Density")
stomata_fresh_dens$Image_name <- sub(".tif", ".jpg", stomata_fresh_dens$Image_name)
stomata_dens <- left_join(stomata_dry_dens, stomata_fresh_dens, by = "Image_name")

#convert to number per mm^2
stomata_dens$Density <- stomata_dens$Density*10
stomata_dens$Density_fresh <- stomata_dens$Density_fresh*10

#parse out plot, tag, strata, and leaf info
stomata_dens$Image_name <- substr(stomata_dens$Image_name,1,nchar(stomata_dens$Image_name)-4)
stomata_dens$leaf <- str_sub(stomata_dens$Image_name, start = -1)
stomata_dens$leaf <- toupper(stomata_dens$leaf)
stomata_dens$Stratum <- str_sub(stomata_dens$Image_name, start = -2, end = -2)
stomata_dens$Plot <- ifelse(startsWith(stomata_dens$Image_name, "S"),str_sub(stomata_dens$Image_name, start = 1, end = 3), str_sub(stomata_dens$Image_name, start = 1, end = 6))
stomata_dens$Tag <- sub("\\-.*", "", stomata_dens$Image_name)
stomata_dens$Tag <- ifelse(startsWith(stomata_dens$Image_name, "S"),str_sub(stomata_dens$Tag, start = 5),str_sub(stomata_dens$Tag, start = 8))
stomata_dens$Stratum <- as.numeric(stomata_dens$Stratum)
stomata_dens$Plot <- ifelse(stomata_dens$Plot == "SEN", "Sendero", stomata_dens$Plot)
stomata_dens$Image_name <- NULL

#combine by plot, tag, stratum, and leaf
stomata <- left_join(stomata_size, stomata_dens, by = c("Plot", "Tag", "Stratum", "leaf"))

#combine with leaves_modern data
leaves_modern <- left_join(leaves_modern, stomata, by = c("Plot", "Tag", "Stratum", "leaf"))
  
#make function to calculate gmax
calculate_gmax <- function(pore_length, guard_cell_width, stomatal_density){
  ###calculate maximum stomatal conductance
  #pore_length = length of the stomatal pore in um
  #guard_cell_width = width of the guard cells in um
  #stomatal_density = number of stomata per m^2
  
  #constants
  d <- 0.0000249 # diffusivity of water vapor at 25 deg C in m^2 s^-1
  v <- 0.0224 # molar volume of air in m^3 mol^-1
  
  #convert stomatal density to m^2
  sd <- stomatal_density * 1e6
  
  #convert pore length to m
  pl <- pore_length * 1e-6
  
  #convert guard cell width to m
  gcw <- guard_cell_width * 1e-6
  
  #calculate the area of the stomatal pore in m^2
  amax <- pi * (pl/2) * gcw #maximum stomatal pore area
  
  #calculate stomatal depth, equal to guard cell width
  pd <- gcw
  
  #calculate the maximum stomatal conductance in mol m^-2 s^-1
  gmax <- ((d/v)*sd*amax) / 
    (pd + (pi/2)*sqrt(amax/pi))
  
  return(gmax)
}

leaves_modern$gmax <- calculate_gmax(leaves_modern$Pore_length, leaves_modern$Guard_cell_width, leaves_modern$Density)
leaves_modern$gmax_fresh <- calculate_gmax(leaves_modern$Pore_length_fresh, leaves_modern$Guard_cell_width_fresh, leaves_modern$Density_fresh)

###analysis of dry and fresh modern leaf measurements
#add tree identifier to include as random effect
leaves_modern$treeID <- paste0(leaves_modern$Plot, leaves_modern$Tag)

#thickness
sf1 <- ggplot(data = leaves_modern, aes(x = thickness, y = thickness_fresh)) +
  geom_point(size = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  xlab("mm") +
  ylab("mm") +
  #add title
  ggtitle("Leaf thickness") +
  theme_bw()
sf1

m <- lmer(thickness_fresh ~ thickness + (1|species/treeID), data=leaves_modern)
summary(m)
r.squaredGLMM(m)

#lma
sf2 <- ggplot(data = leaves_modern, aes(x = lma, y = lma_fresh)) + 
  geom_point(size = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  #geom_text(hjust = 0, vjust = 0) +
  xlab(bquote("g"~m^-2)) +
  ylab(bquote("g"~m^-2)) +
  ggtitle("LMA") +
  theme_bw() 
sf2

m <- lmer(lma_fresh ~ lma + (1|species/treeID), data=leaves_modern) 
summary(m)
r.squaredGLMM(m)

#elw
sf3 <- ggplot(data = leaves_modern, aes(x = elw, y = elw_fresh)) +
  geom_point(size = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  #geom_text(hjust = 0, vjust = 0) +
  xlab("cm") +
  ylab("cm") +
  #add title
  ggtitle("ELW") +
  theme_bw()
sf3
m <- lmer(elw_fresh ~ elw + (1|species/treeID), data=leaves_modern) 
summary(m)
r.squaredGLMM(m)
m <- summary(m)
#assign slope and intercept to a variable to incorporate into leaf temp models later
elw_s <- coef(m)[2] 
elw_i <- coef(m)[1]

#pore length
sf4 <- ggplot(data = leaves_modern, aes(x = Pore_length, y = Pore_length_fresh)) +
  geom_point(size = 0.7) +
  geom_smooth(method = "lm", se = FALSE, na.rm = TRUE) +
  #geom_text(hjust = 0, vjust = 0) +
  xlab("um") +
  ylab("um") +
  #add title
  ggtitle("Pore length") +
  theme_bw()
sf4
m <- lmer(Pore_length_fresh ~ Pore_length + (1|species/treeID), data=leaves_modern)
summary(m)
r.squaredGLMM(m)

#guard cell width
sf5 <- ggplot(data = leaves_modern, aes(x = Guard_cell_width, y = Guard_cell_width_fresh)) +
  geom_point(size = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  #label each point with the leaf code
  geom_text(aes(label = Tag), hjust = 0, vjust = 0) +
  xlab("um") +
  ylab("um") +
  #add title
  ggtitle("Guard cell width") +
  theme_bw()
sf5
m <- lmer(Guard_cell_width_fresh ~ Guard_cell_width + (1|species/treeID), data=leaves_modern)
summary(m)
r.squaredGLMM(m)

#stomatal size
leaves_modern$stomatal_size <- (pi * (leaves_modern$Pore_length/2) * leaves_modern$Guard_cell_width)
leaves_modern$stomatal_size_fresh <- (pi * (leaves_modern$Pore_length_fresh/2) * leaves_modern$Guard_cell_width_fresh)

sf6 <- ggplot(data = filter(leaves_modern, !is.na(stomatal_size_fresh)), aes(x = stomatal_size, y = stomatal_size_fresh)) +
  geom_point(size = 0.7, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, na.rm = TRUE) +
  xlab(bquote(""~um^2)) +
  ylab(bquote(""~um^2)) +
  #add title
  ggtitle("Stomatal size") +
  theme_bw()
sf6
m <- lmer(stomatal_size_fresh ~ stomatal_size + (1|species/treeID), data=leaves_modern)
summary(m)
r.squaredGLMM(m)

#density
sf7 <- ggplot(data = filter(leaves_modern, !is.na(Density_fresh)), aes(x = Density, y = Density_fresh)) +
  geom_point(size = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  xlab(bquote(""~mm^-2)) +
  ylab(bquote(""~mm^-2)) +
  #add title
  ggtitle("Stomatal density") +
  theme_bw()
sf7
m <- lmer(Density_fresh ~ Density + (1|species/treeID), data=leaves_modern)
summary(m)
r.squaredGLMM(m)

#gmax
sf8 <- ggplot(data = filter(leaves_modern, !is.na(gmax_fresh)), aes(x = gmax, y = gmax_fresh)) +
  geom_point(size = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  xlab(bquote("mol m-2 s-1")) +
  ylab(bquote("mol m-2 s-1")) +
  #add title
  ggtitle("gsmax") +
  theme_bw()
sf8
m <- lmer(gmax_fresh ~ gmax + (1|species/treeID), data=leaves_modern)
summary(m)
r.squaredGLMM(m)
m <- summary(m)
#assign slope and intercept to a variable to incorporate into leaf temp models later
gmax_s <- coef(m)[2] 
gmax_i <- coef(m)[1]

#combine sf plots in one figure
grid.arrange(sf1, sf2, sf3, sf6, sf7, sf8, ncol = 3, bottom = "Dry measurements", left = "Fresh measurements")

###analysis of average tree measurements
#average values for each trait at each stratum for each tree
trees_traits <- leaves_modern %>%
  group_by(Plot, Tag, Stratum) %>%
  dplyr::summarize("thickness" = mean(thickness),
            "thickness_fresh" = mean(thickness_fresh),
            "elw" = mean(elw),
            "elw_fresh" = mean(elw_fresh),
            "area" = mean(area),
            "area_fresh" = mean(leaf_area_fresh),
            "lma" = mean(lma),
            "lma_fresh" = mean(lma_fresh),
            "mass" = mean(mass),
            "guard_cell_width" = mean(Guard_cell_width, na.rm=TRUE),
            "guard_cell_width_fresh" = mean(Guard_cell_width_fresh, na.rm=TRUE),
            "pore_length" = mean(Pore_length, na.rm=TRUE),
            "pore_length_fresh" = mean(Pore_length_fresh, na.rm=TRUE),
            "stomatal_density" = mean(Density, na.rm=TRUE),
            "stomatal_density_fresh" = mean(Density_fresh, na.rm=TRUE),
            "gmax" = mean(gmax, na.rm=TRUE),
            "gmax_fresh" = mean(gmax_fresh, na.rm=TRUE),
            "year" = 2023,
            "month" = 9)

trees_traits <- left_join(trees_traits, trees2023, by = c("Plot" = "Plot", "Tag" = "Tag"))
trees_traits$treeID <- paste(trees_traits$Plot, trees_traits$Tag, sep = "")
#change stratum values to "bottom" and "top"
trees_traits$Stratum <- ifelse(trees_traits$Stratum == 1, "bottom", "top")
#calculate stomatal size
trees_traits$stomatal_size <- (pi * (trees_traits$pore_length/2) * trees_traits$guard_cell_width)

#which trees have data for pore_length but not stomatal_density?
x <- trees_traits %>%
  filter(is.na(stomatal_density) & !is.na(pore_length))
#what about the opposite?
y <- trees_traits %>%
  filter(is.na(pore_length) & !is.na(stomatal_density))

#what are the most common species?
trees_traits %>%
  group_by(Species) %>%
  dplyr::summarize("n" = n()) %>%
  arrange(desc(n)) %>%
  head(10)

####analyze if crown position affects leaf traits####
#thickness
x <- trees_traits[, c("treeID", "Stratum", "thickness")] %>%
  pivot_wider(names_from = Stratum,
              values_from = thickness)
t.test(x$bottom, x$top, paired = TRUE)
a <- ggplot(data = trees_traits, aes(x = Stratum, y = thickness)) +
  geom_boxplot(aes(group = Stratum)) +
  theme_bw() +
  ylab("mm") +
  xlab("") +
  ggtitle("Thickness") +
  geom_signif(comparisons = list(c("bottom", "top")), 
              map_signif_level=TRUE,
              test = "t.test",
              test.args = list(paired = TRUE)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  theme_bw() 

#lma
x <- trees_traits[, c("treeID", "Stratum", "lma")] %>%
  pivot_wider(names_from = Stratum,
              values_from = lma)
t.test(x$bottom, x$top, paired = TRUE)
b <- ggplot(data = trees_traits, aes(x = Stratum, y = lma)) +
  geom_boxplot(aes(group = Stratum)) +
  theme_bw() +
  ylab(bquote("g"~m^-2)) +
  xlab("") +
  ggtitle("LMA") +
  geom_signif(comparisons = list(c("bottom", "top")), 
              map_signif_level=TRUE,
              test = "t.test",
              test.args = list(paired = TRUE)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  theme_bw() 

#elw
x <- trees_traits[, c("treeID", "Stratum", "elw")] %>%
  pivot_wider(names_from = Stratum,
              values_from = elw)
t.test(x$bottom, x$top, paired = TRUE)
c <- ggplot(data = trees_traits, aes(x = Stratum, y = elw)) +
  geom_boxplot(aes(group = Stratum)) +
  theme_bw() +
  ylab("cm") +
  xlab("") +
  ggtitle("ELW") +
  geom_signif(comparisons = list(c("bottom", "top")), 
              map_signif_level=TRUE,
              test = "t.test",
              test.args = list(paired = TRUE)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  theme_bw() 

#stomatal size
x <- trees_traits[, c("treeID", "Stratum", "stomatal_size")] %>%
  pivot_wider(names_from = Stratum,
              values_from = stomatal_size)
t.test(x$bottom, x$top, paired = TRUE)
d <- ggplot(data = trees_traits, aes(x = Stratum, y = stomatal_size)) +
  geom_boxplot(aes(group = Stratum)) +
  theme_bw() +
  ylab(bquote(""~um^2)) +
  xlab("") +
  ggtitle("Stomatal size") +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  theme_bw() 

#stomatal density
x <- trees_traits[, c("treeID", "Stratum", "stomatal_density")] %>%
  pivot_wider(names_from = Stratum,
              values_from = stomatal_density)
t.test(x$bottom, x$top, paired = TRUE)
e <- ggplot(data = trees_traits, aes(x = Stratum, y = stomatal_density)) +
  geom_boxplot(aes(group = Stratum)) +
  theme_bw() +
  ylab(bquote(""~mm^-2)) +
  xlab("") +
  ggtitle("Stomatal density") +
  geom_signif(comparisons = list(c("bottom", "top")), 
              map_signif_level=TRUE,
              test = "t.test",
              test.args = list(paired = TRUE)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  theme_bw() 

#gmax
x <- trees_traits[, c("treeID", "Stratum", "gmax")] %>%
  pivot_wider(names_from = Stratum,
              values_from = gmax)
t.test(x$bottom, x$top, paired = TRUE)
f <- ggplot(data = trees_traits, aes(x = Stratum, y = gmax)) +
  geom_boxplot(aes(group = Stratum)) +
  theme_bw() +
  ylab("mol m-2 s-1") +
  xlab("") +
  ggtitle("gmax") +
  geom_signif(comparisons = list(c("bottom", "top")), 
              map_signif_level=TRUE,
              test = "t.test",
              test.args = list(paired = TRUE)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  theme_bw() 
f
grid.arrange(a, b, c, d, e, f, ncol = 3, bottom = "Crown position")

####analyze if DBH affects traits####
#thickness
m <- lmer(thickness ~ DBH2023 + (1|Species) + (1|Stratum) , data = trees_traits)
summary (m) #DBH positively affects thickness
r.squaredGLMM(m)
a <- ggplot(data = trees_traits, aes(x = DBH2023, y = thickness, color = Stratum)) +
  geom_point(size = .9) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("blue2", "red2")) +
  theme_bw() +
  ylab("mm") +
  xlab("") +
  ggtitle("Thickness")+
  theme(legend.position = "none")
a

#lma
m <- lmer(lma ~ DBH2023 + (1|Species) + (1|Stratum), trees_traits)
summary(m) #DBH positively affects lma
r.squaredGLMM(m)
b <- ggplot(data = trees_traits, aes(x = DBH2023, y = lma, color = Stratum)) +
  geom_point(size = 0.9) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("blue2", "red2")) +
  theme_bw() +
  ylab(bquote("g"~m^-2)) +
  xlab("") +
  ggtitle("LMA")+
  theme(legend.position = "none")
b

#elw
m <- lmer(elw ~ DBH2023 + (1|Species) + (1|Stratum), trees_traits)
summary(m) #DBH negatively affects elw
r.squaredGLMM(m)
c <- ggplot(data = trees_traits, aes(x = DBH2023, y = elw, color = Stratum)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("blue2", "red2")) +
  theme_bw() +
  ylab("cm") +
  xlab("") +
  ggtitle("ELW")+
  theme(legend.position = "none")
c

#stomatal size  
m <- lmer(stomatal_size ~ DBH2023 + (1|Species) + (1|Stratum), trees_traits)
summary(m) #DBH does not effect stomatal size
r.squaredGLMM(m)
d <- ggplot(data = trees_traits, aes(x = DBH2023, y = stomatal_size, color = Stratum)) +
  geom_point() +
  scale_color_manual(values = c("blue2", "red2")) +
  theme_bw() +
  ylab(bquote(""~um^2)) +
  xlab("") +
  ggtitle("Stomatal size")+
  theme(legend.position = "none")
d

#stomatal density
m <- lmer(stomatal_density ~ DBH2023 + (1|Species) + (1|Stratum), trees_traits)
summary(m) #DBH does not effect density
r.squaredGLMM(m)
e <- ggplot(data = trees_traits, aes(x = DBH2023, y = stomatal_density, color = Stratum)) +
  geom_point() +
  scale_color_manual(values = c("blue2", "red2")) +
  theme_bw() +
  ylab(bquote("#/"~mm^2)) +
  xlab("") +
  ggtitle("Stomatal density") +
  theme(legend.position = "none")
e

#gmax
m <- lmer(gmax ~ DBH2023 + (1|Species) + (1|Stratum), trees_traits)
summary(m) #DBH does not effect gmax
r.squaredGLMM(m)
f <- ggplot(data = trees_traits, aes(x = DBH2023, y = gmax, color = Stratum)) +
  geom_point() +
  scale_color_manual(values = c("blue2", "red2")) +
  theme_bw() +
  ylab("mol m-2 s-1") +
  xlab("") +
  ggtitle("gmax")+
  theme(legend.position = "none")
f

grid.arrange(a, b, c, d, e, f, ncol = 3, bottom = "DBH (cm)")

#write leaves_modern for additional spectral analysis
#write.csv(leaves_modern, "leaves_modern_cleaned.csv", row.names = FALSE)

########################
##### historic leaves data
#read in historic leaves data and get the Plot and tag 
leaves_historic <- read_excel("leaves_historic.xlsx")
leaves_historic$CollectionNumber <- as.character(leaves_historic$CollectionNumber)
source_trees <- read.csv("source_trees_collected.csv")
leaves_historic <- left_join(leaves_historic, source_trees, by = "CollectionNumber")
leaves_historic <- leaves_historic[,-c(18,19,22)]
#remove rows with "new leaf?" in the leaf_notes column
leaves_historic <- leaves_historic[!grepl("new leaf?", leaves_historic$leaf_notes),]
#remove rows with collection number 45941a, 57786, or 57770
leaves_historic <- leaves_historic[leaves_historic$CollectionNumber != "45941a" & leaves_historic$CollectionNumber != "57786" & leaves_historic$CollectionNumber != "57770",]
leaves_historic$sterile_collection <- ifelse(leaves_historic$sterile_collection == 1, "bottom", "top")
#rename column
colnames(leaves_historic)[colnames(leaves_historic) == "sterile_collection"] <- "Stratum"

#average the three thickness measurements
leaves_historic$thickness <- apply(leaves_historic[,c("thickness1", "thickness2", "thickness3")], 1, mean)
leaves_historic$thickness <- leaves_historic$thickness/1000 #convert from um to mm

#calculate elw
leaves_historic$elw <- leaves_historic$circle_perimeter/pi

#convert mass to grams
leaves_historic$mass <- leaves_historic$mass/1000 #convert from mg to g
leaves_historic$mass_w_glue <- leaves_historic$mass_w_glue/1000 #convert from mg to g

#LMA
leaves_historic$lma <- ifelse(!is.na(leaves_historic$mass),  #if the mass without glue is missing, use the mass with glue
                                  (leaves_historic$mass)/leaves_historic$leaf_area*1000,
                                  (leaves_historic$mass_w_glue)/leaves_historic$leaf_area*1000) #divide by 1000 to convert mg to g

#historical stomata data 
stomata_hist <- read_excel("C:/Users/rfortier/Dropbox/Peru - Individual Tree Project/Stomatal peels/Stomatal analysis.xlsx", 
                           sheet = "HISTORIC")
stomata_hist <- stomata_hist %>%
  group_by(Image_name) %>%
  summarise("Pore_length" = mean(Pore_length),
            "Guard_cell_width" = mean(Guard_cell_width))
stomata_hist_dens <- read_excel("C:/Users/rfortier/Dropbox/Peru - Individual Tree Project/Stomatal peels/Stomatal analysis.xlsx", 
                                sheet = "HISTORIC_Density")
stomata_hist_dens$Density <- stomata_hist_dens$Density*10 #convert to mm^2
stomata_hist <- left_join(stomata_hist, stomata_hist_dens, by = "Image_name")

#parse out collector, coll no, and leaf info
stomata_hist$Image_name <- substr(stomata_hist$Image_name,1,nchar(stomata_hist$Image_name)-4)
stomata_hist$leaf <- str_sub(stomata_hist$Image_name, start = -1)
stomata_hist$leaf <- toupper(stomata_hist$leaf)
#extract everything before the underscore in Image name column
stomata_hist$Collector <- str_extract(stomata_hist$Image_name, "^[^_]+")
#extract number in Image name column
stomata_hist$CollectionNumber <- str_extract(stomata_hist$Image_name, "[0-9]+")
stomata_hist$Image_name <- NULL

#combine with other historical data
leaves_historic <- left_join(leaves_historic, stomata_hist, by = c("Collector", "CollectionNumber", "leaf"), relationship = "many-to-many")

#write leaves_historic for additional spectral analysis
#write.csv(leaves_historic, "leaves_historic_cleaned.csv", row.names = FALSE)

#make supp table 1
supp_table1 <- leaves_historic %>%
  group_by(Collector, CollectionNumber, Month, year, Species.x, Plot, Tag) %>%
  dplyr::summarize("Herbarium" = paste(unique(Herbarium), collapse=", "),
            "n" = length(leaf))
#write.csv(supp_table1, "supp_table1.csv")

#There are CollectionNumber values that are present in source_trees and not in supp_table1. This is due to having collected source trees in the field and subsequently not finding the associated specimen in the herbaria.
#how many unique combinations of Plot and Tag are there in the historic leaves data?
length(unique(paste(supp_table1$Plot, supp_table1$Tag, sep = ""))) #150. This is the number of source trees.

#how many unique combinations of Plot and Tag have multiple collection numbers?
supp_table1 %>%
  group_by(Plot, Tag) %>%
  dplyr::summarize("n" = n()) %>%
  filter(n > 1) %>%
  nrow() #12. this is the number of trees associated with 2 historical specimens

#how many species are there?
length(unique(supp_table1$Species.x)) #40

stomata_spp <- c("Pouteria torta",
"Protium altissima",
"Brosimum lactescens",
"Eschweilera coriacea",
"Pseudolmedia macrophylla",
"Ocotea bofo",
"Licania heteromorpha",
"Pseudolmedia laevigata",
"Protium stevensonii",
"Helicostylis tomentosa",
"clarisia racemosa",
"Meliosma herbertii",
"Matisia ochrocalyx",
"Oxandra xylopioides",
"Ocotea indet",
"Siparuna decipiens",
"Micropholis guianensis")

#summarize averages of stomatal traits
stomata_hist <- stomata_hist %>%
  group_by(Collector, CollectionNumber) %>%
  dplyr::summarize("guard_cell_width" = mean(Guard_cell_width),
                   "pore_length" = mean(Pore_length),
                   "stomatal_density" = mean(Density))

###summarize other traits
trees_traits_hist <- leaves_historic %>%
  group_by(Collector, CollectionNumber, Plot, Tag, year, Month, Stratum) %>%
  dplyr::summarize("thickness" = mean(thickness),
            "elw" = mean(elw),
            "area" = mean(leaf_area),
            "lma" = mean(lma))

colnames(trees_traits_hist)[colnames(trees_traits_hist) == "Month"] <- "month"
trees_traits_hist$month <- as.numeric(trees_traits_hist$month)

#combine leaf traits and stomatal traits
trees_traits_hist <- left_join(trees_traits_hist, stomata_hist, by = c("Collector", "CollectionNumber"))

trees_traits_hist$Tag <- as.character(trees_traits_hist$Tag)

#make one long dataframe for data analysis
trees_traits_all <- bind_rows(trees_traits,trees_traits_hist)
trees_traits_all$treeID <- paste(trees_traits_all$Plot, trees_traits_all$Tag, sep = "")
trees_traits_all$group <- ""
trees_traits_all$group[trees_traits_all$year <= 1990] <- "1980s"
trees_traits_all$group <- ifelse(trees_traits_all$year > 2000 & trees_traits_all$year <= 2010, "2000s", trees_traits_all$group)
trees_traits_all$group <- ifelse(trees_traits_all$year > 2010 & trees_traits_all$year <= 2020, "2010s", trees_traits_all$group)
trees_traits_all$group[trees_traits_all$year > 2020] <- "2023"
trees_traits_all <- drop_na(trees_traits_all,Plot)

trees_traits_all$gmax <- calculate_gmax(trees_traits_all$pore_length, trees_traits_all$guard_cell_width, trees_traits_all$stomatal_density)

#add the species for each tree
trees2023 <- read.csv("trees_collected.csv")

x <- left_join(trees_traits_all, trees2023, by = c("Plot", "Tag"))
trees_traits_all$Species <- x$Species.y
rm(x)

historic_trees <- trees_traits_all %>%
  filter(year <= 2020)

historic_trees %>%
  group_by(Species) %>%
  dplyr::summarize("n" = n()) %>%
  arrange(desc(n)) %>%
  head(10)

#how many trees are in the historic trees that are a species in the stomata_spp vector?
length(historic_trees$Species[historic_trees$Species %in% stomata_spp])

##############
#make a data frame with only source trees, removing 2023 conspecifics. 
#if the historical collection is sterile, compare with modern leaves from bottom of crown.
#if the historical collection is not sterile, compare with modern leaves from top of crown.

#first remove trees that don't have corresponding historical collections
ind.trees <- trees_traits_all %>%
  group_by(treeID) %>%
  filter(any(year < 2023)) %>%
  ungroup()

#duplicate ind.trees before filtering out strata and add a second group column
ind.trees.supp <- ind.trees 
ind.trees.supp$group2 <- ifelse(ind.trees.supp$year < 2023, "historical", paste("2023", ind.trees.supp$Stratum, sep = "_"))
#reorder so historical shows up before 2023 bottom, then 2023 top
ind.trees.supp$group2 <- factor(ind.trees.supp$group2, levels = c("historical", "2023_bottom", "2023_top"))

#now filter out strata from 2023 that don't match the same tree's stratum from historical years
ind.trees <- ind.trees %>%
  group_by(treeID) %>% # Group by tree ID
  filter(!(year == 2023 & !Stratum %in% Stratum[year < 2023])) %>%
  ungroup()

#change all NaN to NA
ind.trees <- ind.trees %>% replace_na(list(guard_cell_width = NA, guard_cell_width_fresh = NA, pore_length = NA, pore_length_fresh = NA, stomatal_density = NA, stomatal_density_fresh = NA, gmax = NA))

#make a table of which traits were measured for each species
traits_table <- ind.trees %>%
  group_by(Species) %>%
  dplyr::summarize("thickness" = sum(!is.na(thickness)),
            "elw" = sum(!is.na(elw)),
            "lma" = sum(!is.na(lma)),
            "guard_cell_width" = sum(!is.na(guard_cell_width)),
            "pore_length" = sum(!is.na(pore_length)),
            "stomatal_density" = sum(!is.na(stomatal_density)),
            "gmax" = sum(!is.na(gmax)))

#analyze thickness by year
m <- lmer(thickness*100 ~ year + (1|Species/treeID), data=ind.trees) #year negatively affects thickness
summary(m)
r.squaredGLMM(m)

a <- ggplot(data = ind.trees, aes(x = year, y = thickness)) +
  geom_line(data = ind.trees, (aes(group = treeID)), color = "#82C0AB", alpha = 0.5) +
  geom_point(size = 1) +
  geom_smooth(method = "lm", se = T, color = "#2044BD") +
  theme_bw() +
  ylab("Thickness (mm)") +
  xlab("") 

#Use paired t-test to analyze if there is a difference between historical and modern leaves (both top and bottom)
x <- ind.trees.supp[, c("treeID", "group2", "thickness")] %>%
  pivot_wider(names_from = group2,
              values_from = thickness,
              values_fn = list(thickness = mean))
m1 <- t.test(x$historical, x$`2023_bottom`, paired = TRUE)
p1 <- ifelse(m1$p.value < 0.001, "***", 
             ifelse(m1$p.value < 0.01, "**", 
                    ifelse(m1$p.value < 0.05, "*", "ns")))
m2 <- t.test(x$historical, x$`2023_top`, paired = TRUE)
p2 <- ifelse(m2$p.value < 0.001, "***", 
             ifelse(m2$p.value < 0.01, "**", 
                    ifelse(m2$p.value < 0.05, "*", "ns")))

#make boxplots of historical, top, and bottom thickness
a1 <- ggplot(data = ind.trees.supp, aes(x = group2, y = thickness*100)) +
  geom_boxplot() +
  #show signicantly different groups from the paired t tests
  geom_signif(annotations = c(p1, p2),  # replace with your actual significance
              y_position = c(30, 35),      # adjust these to be above your boxplots
              xmin = c(1, 1), 
              xmax = c(2, 3),
              textsize = 4, 
              size = 0.5) +
  theme_bw() +
  ylab("mm") +
  xlab("") +
  ggtitle("Thickness") +
  scale_x_discrete(labels = c("historical" = "Historical", "2023_bottom" = "2023 Bottom", "2023_top" = "2023 Top")) 

#analyze lma by year
m <- lmer(lma ~ year + (1|Species/treeID), data=ind.trees) #year negatively affects lma
summary(m)
r.squaredGLMM(m)

b <- ggplot(data = ind.trees, aes(x = year, y = lma)) +
  geom_line(data = ind.trees, (aes(group = treeID)), color = "#82C0AB", alpha = 0.5) +
  geom_point(size = 1) +
  geom_smooth(method = "lm", se = T, color = "#2044BD") +
  theme_bw() +
  ylab(expression(LMA~(g~m^-3))) +
  xlab("") 
b

#Use paired t-test to analyze if there is a difference between historical and modern leaves (both top and bottom)
x <- ind.trees.supp[, c("treeID", "group2", "lma")] %>%
  pivot_wider(names_from = group2,
              values_from = lma,
              values_fn = list(lma = mean))
m1 <- t.test(x$historical, x$`2023_bottom`, paired = TRUE)
p1 <- ifelse(m1$p.value < 0.001, "***", 
             ifelse(m1$p.value < 0.01, "**", 
                    ifelse(m1$p.value < 0.05, "*", "ns")))
m2 <- t.test(x$historical, x$`2023_top`, paired = TRUE)
p2 <- ifelse(m2$p.value < 0.001, "***", 
             ifelse(m2$p.value < 0.01, "**", 
                    ifelse(m2$p.value < 0.05, "*", "ns")))

#make boxplots of historical, top, and bottom lma
b1 <- ggplot(data = ind.trees.supp, aes(x = group2, y = lma)) +
  geom_boxplot() +
  #show signicantly different groups from the paired t tests
  geom_signif(annotations = c(p1, p2),  # replace with your actual significance
              y_position = c(250, 270),      # adjust these to be above your boxplots
              xmin = c(1, 1), 
              xmax = c(2, 3),
              textsize = 4, 
              size = 0.5) +
  theme_bw() +
  ylab("g m-2") +
  xlab("") +
  ggtitle("LMA") +
  scale_x_discrete(labels = c("historical" = "Historical", "2023_bottom" = "2023 Bottom", "2023_top" = "2023 Top")) 

#analyze elw by year
m <- lmer(elw ~ year + (1|Species/treeID), data=ind.trees) #no year effect on elw
summary(m)
r.squaredGLMM(m)

c <- ggplot(data = ind.trees, aes(x = year, y = elw)) +
  geom_line(data = ind.trees, (aes(group = treeID)), color = "#82C0AB", alpha = 0.5) +
  geom_point(size = 1) +
  theme_bw() +
  ylab("ELW (cm)") +
  xlab("") 

#Use paired t-test to analyze if there is a difference between historical and modern leaves (both top and bottom)
x <- ind.trees.supp[, c("treeID", "group2", "elw")] %>%
  pivot_wider(names_from = group2,
              values_from = elw,
              values_fn = list(elw = mean))
m1 <- t.test(x$historical, x$`2023_bottom`, paired = TRUE)
p1 <- ifelse(m1$p.value < 0.001, "***", 
             ifelse(m1$p.value < 0.01, "**", 
                    ifelse(m1$p.value < 0.05, "*", "ns")))
m2 <- t.test(x$historical, x$`2023_top`, paired = TRUE)
p2 <- ifelse(m2$p.value < 0.001, "***", 
             ifelse(m2$p.value < 0.01, "**", 
                    ifelse(m2$p.value < 0.05, "*", "ns")))

#make boxplots of historical, top, and bottom elw
c1 <- ggplot(data = ind.trees.supp, aes(x = group2, y = elw)) +
  geom_boxplot() +
  #show signicantly different groups from the paired t tests
  geom_signif(annotations = c(p1, p2),  # replace with your actual significance
              y_position = c(15, 16),      # adjust these to be above your boxplots
              xmin = c(1, 1), 
              xmax = c(2, 3),
              textsize = 4, 
              size = 0.5) +
  theme_bw() +
  ylab("cm") +
  xlab("") +
  ggtitle("ELW") +
  scale_x_discrete(labels = c("historical" = "Historical", "2023_bottom" = "2023 Bottom", "2023_top" = "2023 Top")) 

#remove records that only have pore_length for a single year
x <- ind.trees %>%
  group_by(treeID) %>%
  filter(sum(!is.na(pore_length)) > 1)

#analyze stomatal size by year
m <- lmer((pi * (pore_length/2) * guard_cell_width) ~ year + (1|Species/treeID), data=x)
summary(m)
r.squaredGLMM(m)

d <- ggplot(data = x, aes(x = year, y = (pi * (pore_length/2) * guard_cell_width))) +
  geom_line(data = x, (aes(group = treeID)), color = "#82C0AB", alpha = 0.5) +
  geom_point(size = 1) +
  geom_smooth(method = "lm", se = T, color = "#2044BD") +
  theme_bw() +
  ylab(expression(Stomatal~size~(um^2))) +
  xlab("") 
d

#Use paired t-test to analyze if there is a difference between historical and modern leaves (both top and bottom)
x <- ind.trees.supp
x <- x %>%
  group_by(treeID) %>%
  filter(sum(!is.na(pore_length)) > 1)
x$stomatal_size <- pi * (x$pore_length/2) * x$guard_cell_width
x2 <- x[, c("treeID", "group2", "stomatal_size")] %>%
  pivot_wider(names_from = group2,
              values_from = stomatal_size,
              values_fn = list(stomatal_size = mean))
m1 <- t.test(x2$historical, x2$`2023_bottom`, paired = TRUE)
p1 <- ifelse(m1$p.value < 0.001, "***", 
             ifelse(m1$p.value < 0.01, "**", 
                    ifelse(m1$p.value < 0.05, "*", "ns")))
m2 <- t.test(x2$historical, x2$`2023_top`, paired = TRUE)
p2 <- ifelse(m2$p.value < 0.001, "***", 
             ifelse(m2$p.value < 0.01, "**", 
                    ifelse(m2$p.value < 0.05, "*", "ns")))

#make boxplots of historical, top, and bottom 
d1 <- ggplot(data = x, aes(x = group2, y = stomatal_size)) +
  geom_boxplot() +
  #show signicantly different groups from the paired t tests
  geom_signif(annotations = c(p1, p2),  # replace with your actual significance
              y_position = c(55, 60),      # adjust these to be above your boxplots
              xmin = c(1, 1), 
              xmax = c(2, 3),
              textsize = 4, 
              size = 0.5) +
  theme_bw() +
  ylab("um2") +
  xlab("") +
  ggtitle("Stomatal size") +
  scale_x_discrete(labels = c("historical" = "Historical", "2023_bottom" = "2023 Bottom", "2023_top" = "2023 Top")) 

#analyze stomatal density by year
x <- ind.trees %>%
  group_by(treeID) %>%
  filter(sum(!is.na(stomatal_density)) > 1) 

m <- lmer(stomatal_density ~ year + (1|Species/treeID), data=x)
summary(m)
r.squaredGLMM(m)

e <- ggplot(data = x, aes(x = year, y = stomatal_density)) +
  geom_line(data = ind.trees, (aes(group = treeID)), color = "#82C0AB", alpha = 0.5) +
  geom_point(size = 1) +
  theme_bw() +
  ylab(expression(Stomatal~density~(mm^-2))) +
  xlab("") 
e

#Use paired t-test to analyze if there is a difference between historical and modern leaves (both top and bottom)
x <- ind.trees.supp %>%
  group_by(treeID) %>%
  filter(sum(!is.na(stomatal_density)) > 1)
x2 <- x[, c("treeID", "group2", "stomatal_density")] %>%
  pivot_wider(names_from = group2,
              values_from = stomatal_density,
              values_fn = list(stomatal_density = mean))
m1 <- t.test(x2$historical, x2$`2023_bottom`, paired = TRUE)
p1 <- ifelse(m1$p.value < 0.001, "***", 
             ifelse(m1$p.value < 0.01, "**", 
                    ifelse(m1$p.value < 0.05, "*", "ns")))
m2 <- t.test(x2$historical, x2$`2023_top`, paired = TRUE)
p2 <- ifelse(m2$p.value < 0.001, "***", 
             ifelse(m2$p.value < 0.01, "**", 
                    ifelse(m2$p.value < 0.05, "*", "ns")))

#make boxplots of historical, top, and bottom elw
e1 <- ggplot(data = ind.trees.supp, aes(x = group2, y = stomatal_density)) +
  geom_boxplot() +
  #show signicantly different groups from the paired t tests
  geom_signif(annotations = c(p1, p2),  # replace with your actual significance
              y_position = c(2200, 2400),      # adjust these to be above your boxplots
              xmin = c(1, 1), 
              xmax = c(2, 3),
              textsize = 4, 
              size = 0.5) +
  theme_bw() +
  ylab("mm-2") +
  xlab("") +
  ggtitle("Stomatal density") +
  scale_x_discrete(labels = c("historical" = "Historical", "2023_bottom" = "2023 Bottom", "2023_top" = "2023 Top")) 

#analyze gsmax by year
x <- ind.trees %>%
  group_by(treeID) %>%
  filter(sum(!is.na(gmax)) > 1) 

m <- lmer(gmax ~ year + (1|Species/treeID), data=x)
summary(m)
r.squaredGLMM(m)

f <- ggplot(data = x, aes(x = year, y = gmax)) +
  geom_line(data = ind.trees, (aes(group = treeID)), color = "#82C0AB", alpha = 0.5) +
  geom_point(size = 1) +
  geom_smooth(method = "lm", se = T, color = "#2044BD") +
  theme_bw() +
  ylab(expression(g[smax]~(mol~m^-2~s^-1))) +
  xlab("") 
f

#Use paired t-test to analyze if there is a difference between historical and modern leaves (both top and bottom)
x <- ind.trees.supp %>%
  group_by(treeID) %>%
  filter(sum(!is.na(gmax)) > 1)
x2 <- x[, c("treeID", "group2", "gmax")] %>%
  pivot_wider(names_from = group2,
              values_from = gmax,
              values_fn = list(gmax = mean))
m1 <- t.test(x2$historical, x2$`2023_bottom`, paired = TRUE)
p1 <- ifelse(m1$p.value < 0.001, "***", 
             ifelse(m1$p.value < 0.01, "**", 
                    ifelse(m1$p.value < 0.05, "*", "ns")))
m2 <- t.test(x2$historical, x2$`2023_top`, paired = TRUE)
p2 <- ifelse(m2$p.value < 0.001, "***", 
             ifelse(m2$p.value < 0.01, "**", 
                    ifelse(m2$p.value < 0.05, "*", "ns")))

#make boxplots of historical, top, and bottom elw
f1 <- ggplot(data = ind.trees.supp, aes(x = group2, y = gmax)) +
  geom_boxplot() +
  #show signicantly different groups from the paired t tests
  geom_signif(annotations = c(p1, p2),  # replace with your actual significance
              y_position = c(6.5, 7),      # adjust these to be above your boxplots
              xmin = c(1, 1), 
              xmax = c(2, 3),
              textsize = 4, 
              size = 0.5) +
  theme_bw() +
  ylab("mol m-2 s-1") +
  xlab("") +
  ggtitle("Gsmax") +
  scale_x_discrete(labels = c("historical" = "Historical", "2023_bottom" = "2023 Bottom", "2023_top" = "2023 Top")) 

grid.arrange(a, d, b, e, c, f, ncol = 2, bottom = "Year") #portrait
grid.arrange(a, b, c, d, e, f, ncol = 3, bottom = "Year") #landscape
grid.arrange(a1, b1, c1, d1, e1, f1, ncol = 3, bottom = "Group") #landscape


#list all trees with density data but missing pore_length data
x <- ind.trees %>%
  filter(!is.na(pore_length) & is.na(stomatal_density))
#the opposite
y <- ind.trees %>%
  filter(is.na(pore_length) & !is.na(stomatal_density))


#########################
## leaf temperatures
library(tealeaves)

#input climate data by month
climate <- read.csv("cru_monthly_climate.csv")

#calculate average leaf temps for trees
#write a for loop to calculate leaf temperature for each tree in trees_traits_all
T_leaf <- rep(NA, length = nrow(trees_traits_all))
enviro_par <- make_enviropar()
leaf_par <- make_leafpar()
constants <- make_constants()

#since we have dry leaf measurements, we need to model fresh leaf measurements for ELW and gmax
#we model fresh leaf measurements in the loop

for (i in 1:nrow(trees_traits_all)){
  
  #leaf parameters. input leaf trait averages for each tree
  leaf_par$g_sw <- ((trees_traits_all$gmax[i]*gmax_s + gmax_i)*1000000) / (101.3246*1000) #estimate fresh (using slope and intercept of fresh ~ dry model), optimal (multiply by 0.2 [Dow et al. 2014 New Phytologist]) stomatal conductance and convert to umol m-2 s-1 P-1   
  leaf_par$leafsize <- (trees_traits_all$elw[i]*elw_s + elw_i)/100 #estimate fresh width and convert cm to m and 
  leaf_par$abs_s <- 0.95 #absorptance of a leaf to shortwave radiation 
  leaf_par$logit_sr <- 0.001 #stomatal ratio set at 0.001
  
  #environmental parameters. input climate data from collection month and year
  y <- trees_traits_all$year[i]
  m <- trees_traits_all$month[i]
  clim <- climate[climate$Year == y & climate$month == m,]
  
  enviro_par$T_air <- clim$tmax + 273.15 #convert max temp from degrees C to K
  enviro_par$RH <- clim$RH/100 #convert relative humidity from % to fraction 
  
  #calculate and extract leaf temperature if both g_sw and leafsize are not missing
  
  if(!is.na(leaf_par$g_sw) & !is.na(leaf_par$leafsize)) {
         temp <- tleaf(leaf_par, enviro_par, constants, quiet= TRUE, set_units = TRUE)
         T_leaf[[i]] <- temp$T_leaf
         }
  if(is.na(leaf_par$g_sw) | is.na(leaf_par$leafsize)) {
    T_leaf[[i]] <- NA
  }
  
}

trees_traits_all$T_leaf <- T_leaf
#convert T_leaf to celsius
trees_traits_all$T_leaf <- trees_traits_all$T_leaf - 273.15

#remake ind.trees data frame from trees_traits_all now that we have leaf temperatures
ind.trees <- trees_traits_all %>%
  group_by(treeID) %>%
  filter(any(year < 2023)) %>%
  ungroup()

#duplicate ind.trees before filtering out strata and add a second group column
ind.trees.supp <- ind.trees 
ind.trees.supp$group2 <- ifelse(ind.trees.supp$year < 2023, "historical", paste("2023", ind.trees.supp$Stratum, sep = "_"))
#reorder so historical shows up before 2023 bottom, then 2023 top
ind.trees.supp$group2 <- factor(ind.trees.supp$group2, levels = c("historical", "2023_bottom", "2023_top"))

#now filter out strata from 2023 that don't match the same tree's stratum from historical years
ind.trees <- ind.trees %>%
  group_by(treeID) %>% # Group by tree ID
  filter(!(year == 2023 & !Stratum %in% Stratum[year < 2023])) %>%
  ungroup()

#remove records with a treeID that doesn't show up in multiple years
ind.trees <- ind.trees %>%
  group_by(treeID) %>%
  filter(n() > 1)

#change all NaN to NA
ind.trees <- ind.trees %>% replace_na(list(guard_cell_width = NA, guard_cell_width_fresh = NA, pore_length = NA, pore_length_fresh = NA, stomatal_density = NA, stomatal_density_fresh = NA, gmax = NA))

#remove historical trees that don't have 2023 data for T_leaf
x <- ind.trees %>%
  group_by(treeID) %>%
  filter(sum(!is.na(T_leaf)) > 1)

#do a model
m <- lmer(T_leaf ~ year + (1|Species/treeID), data = x)
summary(m)
r.squaredGLMM(m)

#plot leaf temperature by year
ggplot(data = x, aes(x = year, y = T_leaf)) +
  theme_bw() +
  geom_line(aes(group = treeID), color = "#82C0AB", alpha = 1, linewidth = 0.2) +
  geom_point(size = 1.5) +
  geom_abline(slope = 0.036, intercept = -35.74, color = "#C42021", linewidth = 1, linetype = "dashed") +
  geom_smooth(method = "lm", color = "#2044BD", se = T, linewidth = 1.3) +
  ylab("Leaf Temperature (°C)") +
  #add a second y axis for air temperature
  scale_y_continuous(sec.axis = sec_axis(~ . -10, name = "Air Temperature (°C)")) +
  xlab("Year") +
  #make axis text larger
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) 

#Use paired t-test to analyze if there is a difference between historical and modern leaves (both top and bottom)
x <- ind.trees.supp %>%
  group_by(treeID) %>%
  filter(sum(!is.na(T_leaf)) > 1)
x2 <- x[, c("treeID", "group2", "T_leaf")] %>%
  pivot_wider(names_from = group2,
              values_from = T_leaf,
              values_fn = list(T_leaf = mean))
m1 <- t.test(x2$historical, x2$`2023_bottom`, paired = TRUE)
p1 <- ifelse(m1$p.value < 0.001, "***", 
             ifelse(m1$p.value < 0.01, "**", 
                    ifelse(m1$p.value < 0.05, "*", "ns")))
m2 <- t.test(x2$historical, x2$`2023_top`, paired = TRUE)
p2 <- ifelse(m2$p.value < 0.001, "***", 
             ifelse(m2$p.value < 0.01, "**", 
                    ifelse(m2$p.value < 0.05, "*", "ns")))

#make boxplot of historical, top, and bottom 
ggplot(data = ind.trees.supp, aes(x = group2, y = T_leaf)) +
  geom_boxplot() +
  #show signicantly different groups from the paired t tests
  geom_signif(annotations = c(p1, p2),  
              y_position = c(40, 41),      # adjust these to be above your boxplots
              xmin = c(1, 1), 
              xmax = c(2, 3),
              textsize = 4, 
              size = 0.5) +
  theme_bw() +
  ylab("Leaf Temperature (°C)") +
  xlab("Group") +
  ggtitle("") +
  scale_x_discrete(labels = c("historical" = "Historical", "2023_bottom" = "2023 Bottom", "2023_top" = "2023 Top")) 

####make a table of each species with the mean and standard deviation of each trait
traits <- c("thickness", "elw", "lma",  "stomatal_size", "stomatal_density", "gmax", "T_leaf")
trait_table <- trees_traits_all %>%
  group_by(Species) %>%
  summarise(across(all_of(traits), list(mean = mean, sd = sd), na.rm = TRUE))

#change stomatal_density to NA if stomatal_size is NA
trait_table$stomatal_density_mean[is.na(trait_table$stomatal_size_mean)] <- NA
trait_table$stomatal_density_sd[is.na(trait_table$stomatal_size_sd)] <- NA
#same for gmax
trait_table$gmax_mean[is.na(trait_table$stomatal_size_mean)] <- NA
trait_table$gmax_sd[is.na(trait_table$stomatal_size_sd)] <- NA
#same for T_leaf
trait_table$T_leaf_mean[is.na(trait_table$stomatal_size_mean)] <- NA
trait_table$T_leaf_sd[is.na(trait_table$stomatal_size_sd)] <- NA

#make list of species with number of conspecifics and source trees
species_list <- trees_traits_all %>%
  group_by(Species) %>%
  summarise(n = n_distinct(Tag))

n_hist <- ind.trees %>%
  filter(year < 2023) %>%
  group_by(Species) %>%
  summarise(n = n())

#combine
species_list <- merge(species_list, n_hist, by = "Species", all = TRUE)
species_list <- merge(species_list, trait_table, by = "Species", all = TRUE)

#rename columns
species_list <- species_list %>%
  rename(n_cont = n.x, n_hist = n.y)
#replace all NA with "-"
species_list[is.na(species_list)] <- "-"

#write.csv(species_list, "supp_table2.csv")

#####now model leaf temps as if all leaves were collected at the same time
#first identify the hottest and coolest six month periods
#create a rolling sum (or average) of 6 months of max temperatures
climate$rolling_6mo <- rollapply(climate$temp, width = 6, FUN = mean, align = "right", fill = NA)
#find the index of the max and min rolling average
max_index <- which.max(climate$rolling_6mo)
min_index <- which.min(climate$rolling_6mo)

#get the corresponding 6-month period
start_date_hot <- climate$X[max_index - 5]
end_date_hot <- climate$X[max_index]
start_date_cold <- climate$X[min_index - 5]
end_date_cold <- climate$X[min_index]

#save the results
hot_season <- climate[max_index:(max_index - 5),]
cold_season <- climate[min_index:(min_index - 5),]

#initiate leaf temp parameters
T_leaf <- rep(NA, length = nrow(trees_traits_all))
enviro_par <- make_enviropar()
leaf_par <- make_leafpar()
constants <- make_constants()

for (i in 1:nrow(trees_traits_all)){
  
  #leaf parameters. input leaf trait averages for each tree
  leaf_par$g_sw <- ((trees_traits_all$gmax[i]*gmax_s + gmax_i)*1000000) / (101.3246*1000) #estimate fresh stomatal conductance and convert to umol m-2 s-1 Pa-1
  leaf_par$leafsize <- (trees_traits_all$elw[i]*elw_s + elw_i)/100 #estimate fresh width and convert cm to m and 
  leaf_par$abs_s <- 0.85 #absorptance of a leaf to shortwave radiation 
  leaf_par$logit_sr <- 0.001 #stomatal ratio set at 0.001
  
  #environmental parameters. make a new data frame from climate with just one row with averages for each variable
  clim <- climate %>%
    dplyr::summarize(tmax = mean(tmax), 
              RH = mean(RH))
  
  enviro_par$T_air <- clim$tmax + 273.15 #convert max temp from degrees C to K
  enviro_par$RH <- clim$RH/100 #convert relative humidity from % to fraction 
  
  #calculate and extract leaf temperature if both g_sw and leafsize are not missing
  
  if(!is.na(leaf_par$g_sw) & !is.na(leaf_par$leafsize)) {
    temp <- tleaf(leaf_par, enviro_par, constants, quiet= TRUE, set_units = TRUE)
    T_leaf[[i]] <- temp$T_leaf
  }
  if(is.na(leaf_par$g_sw) | is.na(leaf_par$leafsize)) {
    T_leaf[[i]] <- NA
  }
  
}
trees_traits_all$T_leaf_mean <- T_leaf

#same but using the modern hot season climate data
for (i in 1:nrow(trees_traits_all)){
  
  #leaf parameters. input leaf trait averages for each tree
  leaf_par$g_sw <- ((trees_traits_all$gmax[i]*gmax_s + gmax_i)*1000000) / (101.3246*1000) #estimate fresh stomatal conductance and convert to umol m-2 s-1 Pa-1
  leaf_par$leafsize <- (trees_traits_all$elw[i]*elw_s + elw_i)/100 #estimate fresh width and convert cm to m and 
  leaf_par$abs_s <- 0.85 #absorptance of a leaf to shortwave radiation 
  leaf_par$logit_sr <- 0.001 #stomatal ratio set at 0.001
  
  #environmental parameters. use hot climate
  clim <- hot_season %>%
    dplyr::summarize(tmax = mean(tmax), 
              RH = mean(RH))
  
  enviro_par$T_air <- clim$tmax + 273.15 #convert max temp from degrees C to K
  enviro_par$RH <- clim$RH/100 #convert relative humidity from % to fraction 
  
  #calculate and extract leaf temperature if both g_sw and leafsize are not missing
  
  if(!is.na(leaf_par$g_sw) & !is.na(leaf_par$leafsize)) {
    temp <- tleaf(leaf_par, enviro_par, constants, quiet= TRUE, set_units = TRUE)
    T_leaf[[i]] <- temp$T_leaf
  }
  if(is.na(leaf_par$g_sw) | is.na(leaf_par$leafsize)) {
    T_leaf[[i]] <- NA
  }
  
}
trees_traits_all$T_leaf_hot <- T_leaf

#same but using the historical wet season climate data
for (i in 1:nrow(trees_traits_all)){
  
  #leaf parameters. input leaf trait averages for each tree
  leaf_par$g_sw <- ((trees_traits_all$gmax[i]*gmax_s + gmax_i)*1000000) / (101.3246*1000) #estimate fresh stomatal conductance and convert to umol m-2 s-1 Pa-1
  leaf_par$leafsize <- (trees_traits_all$elw[i]*elw_s + elw_i)/100 #estimate fresh width and convert cm to m and 
  leaf_par$abs_s <- 0.85 #absorptance of a leaf to shortwave radiation 
  leaf_par$logit_sr <- 0.001 #stomatal ratio set at 0.001
  
  #environmental parameters. use cold climate
  clim <- cold_season %>%
    dplyr::summarize(tmax = mean(tmax), 
              RH = mean(RH))
  
  enviro_par$T_air <- clim$tmax + 273.15 #convert max temp from degrees C to K
  enviro_par$RH <- clim$RH/100 #convert relative humidity from % to fraction 
  
  #calculate and extract leaf temperature if both g_sw and leafsize are not missing
  
  if(!is.na(leaf_par$g_sw) & !is.na(leaf_par$leafsize)) {
    temp <- tleaf(leaf_par, enviro_par, constants, quiet= TRUE, set_units = TRUE)
    T_leaf[[i]] <- temp$T_leaf
  }
  if(is.na(leaf_par$g_sw) | is.na(leaf_par$leafsize)) {
    T_leaf[[i]] <- NA
  }
  
}
trees_traits_all$T_leaf_cold <- T_leaf

#convert T_leaf to celsius
trees_traits_all$T_leaf_mean <- trees_traits_all$T_leaf_mean - 273.15
trees_traits_all$T_leaf_hot <- trees_traits_all$T_leaf_hot - 273.15
trees_traits_all$T_leaf_cold <- trees_traits_all$T_leaf_cold - 273.15

#remake ind.trees data frame from trees_traits_all now that we have leaf temperatures
ind.trees <- trees_traits_all %>%
  group_by(treeID) %>%
  filter(any(year < 2023)) %>%
  ungroup()

#now filter out strata from 2023 that don't match the same tree's stratum from historical years
ind.trees <- ind.trees %>%
  group_by(treeID) %>% # Group by tree ID
  filter(!(year == 2023 & !Stratum %in% Stratum[year < 2023])) %>%
  ungroup()

#remove records with a treeID that doesn't show up in multiple years
ind.trees <- ind.trees %>%
  group_by(treeID) %>%
  filter(n() > 1)

#change all NaN to NA
ind.trees <- ind.trees %>% replace_na(list(guard_cell_width = NA, guard_cell_width_fresh = NA, pore_length = NA, pore_length_fresh = NA, stomatal_density = NA, stomatal_density_fresh = NA, gmax = NA))

#remove historical trees that don't have 2023 data for T_leaf
#x <- ind.trees %>%
#  group_by(treeID) %>%
#  filter(sum(!is.na(T_leaf)) > 1)

############
#####Run this chunk of code changing the T-leaf column where indicated to analyze historical (cool), mean, and contemporary (hot) climates
#do a t test to see if historical and modern leaves are different
x <- ind.trees %>% 
  filter(group != "2000s" & group != "2010s" & sum(!is.na(T_leaf)) >1) #Don't change T_leaf here. This removes trees that don't have leaf temps for both collections
x <- x[, c("treeID", "group", "T_leaf_hot")] %>% #change T_leaf column to test the cool, mean, and hot climates 
  pivot_wider(names_from = group,
              values_from = T_leaf_hot, #change T_leaf column to test the cool, mean, and hot climates 
              values_fn = mean)
x <- x %>%
  rename("Historical" = `1980s`,
         "Contemporary" = `2023`)

t.test(x$Contemporary, x$Historical, paired = TRUE)

x <- x %>%
  pivot_longer(cols = c("Historical", "Contemporary"),
               names_to = "group",
               values_to = "T_leaf_hot") #change T_leaf column to test the cool, mean, and hot climates

x$group <- factor(x$group, levels = c("Historical", "Contemporary"))

ggplot(data = x, aes(x = group, y = T_leaf_hot)) + #change T_leaf column to visualize the cool, mean, and hot climates
  geom_boxplot(aes(group = group), fill = "gray95", outliers = FALSE) +
  #geom_line(aes(group = treeID), color = "gray50", alpha = 0.4) +
  geom_jitter(width = 0.05, height = 0, size = 1.5, color = "black", alpha = 0.4) +
  theme_bw() +
  ylab("Leaf Temperature (°C)") +
  xlab("") +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

#####

#############sensitivity analysis using Latin hypercubes
library(lhs)

# define parameter ranges
n_samples <- 1000  # Number of Monte Carlo samples
param_bounds <- data.frame(
  tmax = runif(n_samples, 27, 35),  # Air temp: 27-35. This is the range of tmax in climate df
  RH = runif(n_samples, 0.59, 0.82),  # Relative humidity: 0.59-0.82. This is the range of RH in climate df
  gmax = runif(n_samples, min((ind.trees$gmax*gmax_s + gmax_i)*1000000 / (101.3246*1000), na.rm = TRUE), 
                          max((ind.trees$gmax*gmax_s + gmax_i)*1000000 / (101.3246*1000), na.rm = TRUE)),  # Stomatal conductance. This is the range of gmax in umol/m2/s
  elw = runif(n_samples, min((ind.trees$elw*elw_s + elw_i)/100, na.rm = TRUE), max((ind.trees$elw*elw_s + elw_i)/100, na.rm = TRUE))  # ELW (m) 
)

#param_bounds <- data.frame(
#  min = c(27, 0.59, min((ind.trees$gmax*gmax_s + gmax_i)*1000000 / (101.3246*1000), na.rm = TRUE), min((ind.trees$elw*elw_s + elw_i)/100, na.rm = TRUE)),   # Minimum values
#  max = c(35, 0.82, max((ind.trees$gmax*gmax_s + gmax_i)*1000000 / (101.3246*1000), na.rm = TRUE), max((ind.trees$elw*elw_s + elw_i)/100, na.rm = TRUE))        # Maximum values
#)

colnames(param_bounds) <- c("tmax", "RH", "gmax", "elw")

# Step 2: Generate Latin Hypercube Samples
lhs_samples <- randomLHS(n_samples, ncol(param_bounds))  # LHS values between 0 and 1
lhs_samples <- as.data.frame(lhs_samples)

# rescale samples to the parameter ranges
scaled_samples <- as.data.frame(matrix(0, ncol=4, nrow=n_samples))
scaled_samples$V1 <- (lhs_samples$V1 * (max(param_bounds$tmax) - min(param_bounds$tmax)) + min(param_bounds$tmax))
scaled_samples$V2 <- (lhs_samples$V2 * (max(param_bounds$RH) - min(param_bounds$RH)) + min(param_bounds$RH))
scaled_samples$V3 <- (lhs_samples$V3 * (max(param_bounds$gmax) - min(param_bounds$gmax)) + min(param_bounds$gmax))
scaled_samples$V4 <- (lhs_samples$V4 * (max(param_bounds$elw) - min(param_bounds$elw)) + min(param_bounds$elw))

# Convert samples into a dataframe
param_samples <- as.data.frame(scaled_samples)
colnames(param_samples) <- colnames(param_bounds)  # Use original parameter names

# Step 4: Run Monte Carlo simulations
leaf_temps <- numeric(n_samples)  # Preallocate for speed

#Make parameters
enviro_par <- make_enviropar()
leaf_par <- make_leafpar()
constants <- make_constants()

for (i in 1:n_samples) {

    enviro_par$T_air = param_samples$tmax[i] + 273.15
    enviro_par$RH = param_samples$RH[i]/100
    leaf_par$g_sw = param_samples$gmax[i]
    leaf_par$leafsize = param_samples$elw[i]
  
  # Compute leaf temperature
  model <- tleaf(leaf_par, enviro_par, constants, quiet= TRUE, set_units = TRUE)
  leaf_temps[i] <- model$T_leaf
}

# Step 5: Store results in a dataframe
sensitivity_df <- cbind(param_samples, Leaf_Temp = leaf_temps)
sensitivity_df$Leaf_Temp <- sensitivity_df$Leaf_Temp-273.15 # Convert from Kelvin to Celsius

# Step 6: Perform Sensitivity Analysis using Linear Regression
lm_model <- lm(Leaf_Temp ~ ., data = sensitivity_df)
summary(lm_model)  # Check which variables are most influential

####visualize effect sizes
# Standardize parameters and leaf temperature
sensitivity_scaled <- as.data.frame(scale(sensitivity_df))

# Fit a standardized regression model
lm_model_scaled <- lm(Leaf_Temp ~ ., data = sensitivity_scaled)
summary (lm_model_scaled)

# Extract standardized coefficients
coefficients_df <- as.data.frame(summary(lm_model_scaled)$coefficients)
coefficients_df$Parameter <- rownames(coefficients_df)
colnames(coefficients_df) <- c("Estimate", "Std_Error", "t_value", "p_value", "Parameter")

# Remove the intercept and rename columns
coefficients_df <- coefficients_df %>% filter(Parameter != "(Intercept)")
coefficients_df$Parameter <- recode(coefficients_df$Parameter,
                                    "tmax" = "Maximum air temperature",
                                    "gmax" = "Stomatal conductance",
                                    "RH" = "Relative humidity",
                                    "elw" = "Effective leaf width"
)
order <- c("Stomatal conductance", "Effective leaf width", "Relative humidity", "Maximum air temperature")
coefficients_df$Parameter <- factor(coefficients_df$Parameter, levels = order)
coefficients_df$Category <- ifelse(coefficients_df$Parameter %in% c("Stomatal conductance", "Effective leaf width"), 
                                   "Leaf", "Environment")


# Plot effect sizes
ggplot(coefficients_df, aes(x = Parameter, y = Estimate, fill = Category)) +
  geom_abline(slope = 0, intercept = 0.2, color = "gray", linewidth = 0.5, linetype = "solid") +
  geom_abline(slope = 0, intercept = 0.4, color = "gray", linewidth = 0.5, linetype = "solid") +
  geom_abline(slope = 0, intercept = 0.6, color = "gray", linewidth = 0.5, linetype = "solid") +
  geom_abline(slope = 0, intercept = 0.8, color = "gray", linewidth = 0.5, linetype = "solid") +
  geom_abline(slope = 0, intercept = -0.2, color = "gray", linewidth = 0.5, linetype = "solid") +
  geom_abline(slope = 0, intercept = -0.4, color = "gray", linewidth = 0.5, linetype = "solid") +
  geom_abline(slope = 0, intercept = -0.6, color = "gray", linewidth = 0.5, linetype = "solid") +
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  labs(x = "", y = "Effect Size") +
  scale_y_continuous(breaks = seq(-1, 1, 0.2)) +
  scale_fill_manual(values = c("Environment" = "#FED789FF", "Leaf" = "#72874EFF")) +
  geom_abline(slope = 0, intercept = 0, color = "black", linewidth = 1, linetype = "dashed") +
  theme_classic() +
  #add the estimates as text on top of the bars
  geom_text(aes(label = round(Estimate, 2)), 
            position = position_stack(vjust = 0.5), 
            size = 4, 
            color = "black") +
  theme(axis.line.x = element_line(color = "black", size = 1, linetype = "solid"),
        axis.ticks.length.x = unit(.3, "cm"),
        axis.ticks.x = element_line(color = "black", size = 1, linetype = "solid"),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.line.y = element_blank(), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))


 ####
