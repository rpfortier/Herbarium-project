### Find MDD plot herbarium specimens for individual tree project
# set working directory
setwd("C:/Users/rxf568/OneDrive - University of Miami/Peru - Individual Tree Project")
library(dplyr)
library(stringr)
library(xml2)
library(rvest)   ##very important
library(tibble)
library(tidyverse)
library(ggplot2)
library(readxl)

# load Gentry's MDD collections
Gentry_MDD <- read_excel("Gentry_MadredeDios_Collections.xlsx", 
                         col_types = c("numeric", "text", "text", 
                         "text", "text", "text", "text", "text", 
                         "text", "text", "text", "text", "text", 
                         "numeric", "numeric", "date", "text", 
                         "numeric", "text", "numeric", "numeric", 
                         "text", "text"))

# Filter out everything that doesn't contain tree tag info in the description
Gentry_MDD <- filter(Gentry_MDD, grepl('Tree No|Tree no|tree no|Tree #|tree #|Tag|tag', Description))

#Optional, filter out everything that doesn't contain "plot" in the locality
Gentry_MDD <- filter(Gentry_MDD, grepl('plot|Plot|parcela|Parcela', Locality))

# Extract the tag numbers from the description
Gentry_MDD$Tag_no <- gsub("Tree No|Tree No.|Tree no|tree no|Tree #|tree #|Tag|tag|Tree No |Tree No. |Tree no |tree no |Tree # |tree # |Tag |tag ", "", Gentry_MDD$Description)
Gentry_MDD <- Gentry_MDD %>%
  separate(Tag_no, c("Tag", "post", sep = ","))
Gentry_MDD$post <- NULL
Gentry_MDD$`,` <- NULL

  # Fix a few missing tags manually
Gentry_MDD$Tag <- ifelse(Gentry_MDD$CollectionNumber == 51210, 6148, Gentry_MDD$Tag)
Gentry_MDD$Tag <- ifelse(Gentry_MDD$CollectionNumber == 68809, 1457, Gentry_MDD$Tag)
Gentry_MDD$Tag <- ifelse(Gentry_MDD$CollectionNumber == 45742, 412, Gentry_MDD$Tag)
Gentry_MDD$Tag <- ifelse(Gentry_MDD$CollectionNumber == 46003, 3106, Gentry_MDD$Tag)
Gentry_MDD$Tag <- ifelse(Gentry_MDD$CollectionNumber == 46124, 3357, Gentry_MDD$Tag)
Gentry_MDD$Tag <- ifelse(Gentry_MDD$CollectionNumber == 43693, 661, Gentry_MDD$Tag)
Gentry_MDD$Tag <- ifelse(Gentry_MDD$CollectionNumber == 51195, 6155, Gentry_MDD$Tag)
Gentry_MDD$Tag <- ifelse(Gentry_MDD$CollectionNumber == 45767, 4437, Gentry_MDD$Tag)

## Figure out which plots the specimens were collected from?
#Gentry_MDD <- arrange(Gentry_MDD, `*Coordinate Display`)

#ggplot(Gentry_MDD, aes(x=`Date Display`)) + geom_histogram()

### Make list of living trees across all MDD plots
# Load data
Live_trees <- read_excel("Live_plot_trees.xlsx")

# Remove dead trees and change tag variable name to match Gentry_MDD
Live_trees <- Live_trees %>%
  filter(F2 == "1")
Live_trees <- Live_trees %>%
  rename("Tag" = "Tag No",
         "Pv_tag" = "Pv. Tag No")

# Collate Gentry collections and live trees by tag (or previous tag) and family
source_trees <- inner_join( Live_trees, Gentry_MDD, by = c("Tag" = "Tag", "Family" = "FamilyName"))
source_trees2 <- inner_join( Live_trees, Gentry_MDD, by = c("Pv_tag" = "Tag", "Family" = "FamilyName"))

source_trees <- rbind(source_trees, source_trees2)
rm(source_trees2)
source_trees <- as.data.frame(unique((source_trees)))

# Remove some irrelevant columns
source_trees$Habitat <- NULL
source_trees$IsCultivated <- NULL
source_trees$CountryName <- NULL
source_trees$Vegetation <- NULL

# Find duplicate collection numbers and clean them 
n_occur <- data.frame(table(source_trees$CollectionNumber))
# There are many collection IDs that are duplicated. Best to extract them and clean them in Excel, then combine with the non-duplicates.

source_trees_duplicated <- source_trees[duplicated(source_trees$CollectionNumber) | duplicated(source_trees$CollectionNumber, fromLast = TRUE), ]

# Make another df with the non-duplicated numbers
source_trees_unique <- source_trees[!(duplicated(source_trees$CollectionNumber) | duplicated(source_trees$CollectionNumber, fromLast = TRUE)), ]
  
# Export dataframe as .csv
#write.csv(source_trees_duplicated, "Duplicated_Gentry_plot_specimens.csv")

# Reload cleaned dataframe
source_trees_duplicated_Cleaned <- read_csv("Duplicated_Gentry_plot_specimens_Cleaned.csv")

# Combine the unique collection numbers with the cleaned duplicates
source_trees <- rbind(source_trees_duplicated_Cleaned, source_trees_unique)
rm(source_trees_duplicated)
rm(source_trees_duplicated_Cleaned)
rm(source_trees_unique)

# Clean more based on plot tree and specimen locality
source_trees$plot.place <- ifelse(source_trees$Plot == "MNU_01"|source_trees$Plot == "MNU_03"|source_trees$Plot == "MNU_05"|source_trees$Plot == "MNU_06", "Manu", ifelse(source_trees$Plot== "CUZ_01"|source_trees$Plot == "CUZ_02"|source_trees$Plot == "CUZ_03"|source_trees$Plot == "CUZ_04", "Cuzco", "Tambopata"))

source_trees$spec.place <- ifelse(word(source_trees$Locality) == "Cocha", "Manu", ifelse(word(source_trees$Locality) == "Cusco", "Cuzco", "Tambopata"))

source_trees <- subset(source_trees, plot.place == spec.place)
source_trees$plot.place <- NULL
source_trees$spec.place <- NULL

source_trees$plot.genus <- word(source_trees$Species)
source_trees$spec.genus <- word(source_trees$`*FullName`)

subset.1 <- subset(source_trees, plot.genus != spec.genus)
subset.2 <- subset(source_trees, plot.genus == spec.genus)

# Export subset.1
#write.csv(subset.1, "subset1.csv")

# Reread clean subset.1
subset1_clean <- read_csv("subset1_clean.csv")

# Combine subsets to make final list
source_trees <- rbind(subset.2,subset1_clean)
rm(subset.1,subset.2,subset1_clean)

write.csv(source_trees, "source_trees.csv")


