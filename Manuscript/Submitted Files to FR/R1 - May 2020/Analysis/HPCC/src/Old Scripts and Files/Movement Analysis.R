#What about movement in the Okamoto Data?

library(geosphere) #distGeo()
library(data.table) #uniqueN()

# How many fish move reporting grids?
# How far do they move?
# Do multiple recaptures move?
# Does movement likelihood coorospond to DAL?
# Geodesic movement

##### Script Set Up #####
#### Clearing Worksapce
rm(list = ls())

#### Setting a Script Timer
script_timer = proc.time()

#### Declaring Directory Path
src_dir = file.path(getwd(), 'src')
data_dir = file.path(getwd(), 'data')
results_dir = file.path(getwd(), 'results')
figure_dir = file.path(getwd(), 'figure')

#### Installing Principle Dependencies
library('gdata') # read.xls()
library('FSA') # vbFuns()
library('nlstools') # 
library('notifyR') # send_push()
library('minpack.lm') # nlsLM
# install.packages('fishmethods')
library('fishmethods')

##### Data Setup: Loading and Cleaning Data Files #####
#### Mark Recapture Data
## Note: This data set is a fucking mess. Each line is a fish with location of capture, tag date, tagging depth, species, tag id, fork length, remarks, and then duplicate data for all of that each time it was recaptured up to 4 recaptures
mark_recapture_data = read.xls(file.path(data_dir, 'HO Mstr, temp (version 1).xlsx'), stringsAsFactors =  FALSE)

### Loading in reporting grid data
island_areas = read.csv("/Users/stephenscherrer/Google Drive/Weng Lab/Personal_Folders/Steve/dissertation work/BACIP Analysis/data/fish_zone_centroid_coords.txt", col.names = c('fid', 'area', 'perimeter', 'area_id', 'area_a', 'type', 'chart', 'reef_name', 'notes', 'island', 'orig_fid', 'lon', 'lat'))

### Making up approximate coordinates for a couple areas
## Area 351
new_area = island_areas[1, ]
new_area$area_id = 351; new_area$lat = 20.7556; new_area$lon = -157.548
island_areas = rbind(island_areas,new_area)
## Area 452
new_area = island_areas[1, ]
new_area$area_id = 452; new_area$lat = 20.803; new_area$lon = -157.834
island_areas = rbind(island_areas,new_area)

## A getting reporting grid areas
island_areas$island[island_areas$area_id >= 100 & island_areas$area_id < 200] = "hawaii"
island_areas$island[island_areas$area_id >= 300 & island_areas$area_id < 400] = "mkml"
island_areas$island[island_areas$area_id >= 400 & island_areas$area_id < 500] = "oahu"
island_areas$island[island_areas$area_id >= 500 & island_areas$area_id < 600] = "kauai"

## How many total fish do we have in the data set?
dim(mark_recapture_data)[1] # 4245!

### Renaming data columns
colnames(mark_recapture_data) = c('tag_date', 'location', 'station', 'depth_f', 'species', 'previously_tagged', 'tag_id','fork_length_in', 'remarks', 'recapture_1_date', 'recapture_1_location', 'recapture_1_station', 'recapture_1_depth_f', 'recapture_1_fork_length_in', 'weight_1_lbs', 'days_1_free', 'growth_1_in', 'distance_1_miles','retagged_1',
                                  'recapture_2_date', 'recapture_2_location', 'recapture_2_station', 'recapture_2_depth_f', 'recapture_2_fork_length_in', 'weight_2_lbs', 'days_2_free', 'growth_2_in', 'distance_2_miles', 'retagged_2',
                                  'recapture_3_date', 'recapture_3_location', 'recapture_3_station', 'recapture_3_depth_f', 'recapture_3_fork_length_in', 'weight_3_lbs', 'days_3_free', 'growth_3_in', 'distance_3_miles', 'retagged_3',
                                  'recapture_4_date', 'recapture_4_location', 'recapture_4_station', 'recapture_4_depth_f', 'recapture_4_fork_length_in', 'weight_4_lbs', 'days_4_free', 'growth_4_in', 'distance_4_miles', 'x_retagged')

#### Adusting Data Classes
### Formatting Dates (Converting Characters to POSIXct)
mark_recapture_data$tag_date = as.POSIXct(mark_recapture_data$tag_date, format = "%Y-%m-%d")
mark_recapture_data$recapture_1_date = as.POSIXct(mark_recapture_data$recapture_1_date, format = "%Y-%m-%d")
mark_recapture_data$recapture_2_date = as.POSIXct(mark_recapture_data$recapture_2_date, format = "%Y-%m-%d")
mark_recapture_data$recapture_3_date = as.POSIXct(mark_recapture_data$recapture_3_date, format = "%Y-%m-%d")
mark_recapture_data$recapture_4_date = as.POSIXct(mark_recapture_data$recapture_4_date, format = "%Y-%m-%d")

### Formatting fork lengths 
## Note: There are a couple fork lengths that have ?, *, or have no lengths recorded. 
## I have no idea what these are but they're qualifiers and so I'm going to let them go to NA and get dropped from analysis
in_to_cm = 2.54
mark_recapture_data$fork_length_cm = as.numeric(mark_recapture_data$fork_length_in) * in_to_cm
mark_recapture_data$recapture_1_fork_length_cm = as.numeric(mark_recapture_data$recapture_1_fork_length_in) * in_to_cm
mark_recapture_data$recapture_2_fork_length_cm = as.numeric(mark_recapture_data$recapture_2_fork_length_in) * in_to_cm
mark_recapture_data$recapture_3_fork_length_cm = as.numeric(mark_recapture_data$recapture_3_fork_length_in) * in_to_cm
mark_recapture_data$recapture_4_fork_length_cm = as.numeric(mark_recapture_data$recapture_4_fork_length_in) * in_to_cm

### Subsetting out only Opakapaka with tag IDs - That is, fish that were marked
mark_recapture_data = mark_recapture_data[mark_recapture_data$species == '1' & mark_recapture_data$tag_id != '', ]
dim(mark_recapture_data)[1] # This gets you to the previously published 4179 tagged paka number from Kobayashi, Okamoto, & Oishi . for some reason doesn't exclude fish marked 'died'





##### Analysis: #####
paka_moves = mark_recapture_data



#### How many fish move reporting grids?
head(paka_moves)
unique(paka_moves$location)
unique(paka_moves$recapture_1_location)

### Cleaning recapture locations
"Waianae side of Kaena pt." = 403 | 423
"Makapuu" = 408 | 428
"Kahoolawe" = 308 | 309
"Ilio Pt" = 312 | 311 | 332
"Ewa bch" = 401 | 421
"P.B. 1st. F" = 331
"kaena pt." = 404 | 403 | 424 | 423
"Kalaupapa" = 312 | 332
"Niihau-Kaula" = 508 | 528
"Kaula Rck." = 508 | 528
"Hnl. side banks" = 331
"Ewa Bch" = 401 | 421

move_df = data.frame("tag_id" = paka_moves$tag_id, 
           "liberty" = 0,
           "location_0" = paka_moves$location,
           "location_1" = paka_moves$recapture_1_location, 
           "location_2" = paka_moves$recapture_2_location, 
           "location_3" = paka_moves$recapture_3_location, 
           "location_4" = paka_moves$recapture_4_location,
           "moved" = NA, 
           stringsAsFactors = FALSE)

for(i in 1:length(move_df$tag_id)){
  tag_date = paka_moves$tag_date
    if(is.na(paka_moves$recapture_4_date[i]) == FALSE){
      move_df$liberty[i] = as.numeric(difftime(paka_moves$recapture_4_date[i], paka_moves$tag_date[i], units = "days"))
    } else if(is.na(paka_moves$recapture_3_date[i]) == FALSE){
      move_df$liberty[i] = as.numeric(difftime(paka_moves$recapture_3_date[i], paka_moves$tag_date[i], units = "days"))
    } else if(is.na(paka_moves$recapture_2_date[i]) == FALSE){
      move_df$liberty[i] = as.numeric(difftime(paka_moves$recapture_2_date[i], paka_moves$tag_date[i], units = "days"))
    } else if(is.na(paka_moves$recapture_1_date[i]) == FALSE){
      move_df$liberty[i] = as.numeric(difftime(paka_moves$recapture_1_date[i], paka_moves$tag_date[i], units = "days"))
    } else(move_df$liberty[i] = NA)
    }


paka_moves = mark_recapture_data[! is.na(mark_recapture_data$recapture_1_date), ]
  # dim(paka_moves)[1]
    # 487

## Converting recapture locations to numerics - Note, NAs are produced for any data that doesnt have a valid statistical grid location
paka_moves$recapture_1_location = as.numeric(paka_moves$recapture_1_location)
paka_moves$recapture_2_location = as.numeric(paka_moves$recapture_2_location)
paka_moves$recapture_3_location = as.numeric(paka_moves$recapture_3_location)
paka_moves$recapture_4_location = as.numeric(paka_moves$recapture_4_location)
move = data.frame("tag_id" = paka_moves$tag_id, "moved1" = NA, 'moved2' = NA, 'moved3' = NA, 'moved4' = 0, stringsAsFactors = FALSE)

## Getting number of moements for each fish and the total distance traveled as computed from the centeroid of the tagging location
moves_made = data.frame("tag_id" = paka_moves$tag_id, "movements" = 0, "distance" = 0)
move_hist = cbind(paka_moves$location, paka_moves$recapture_1_location, paka_moves$recapture_2_location, paka_moves$recapture_3_location, paka_moves$recapture_4_location)
for(i in 1:length(move_hist[ ,1])){
  print(i)
  if (!is.na(move_hist[i,5]) & move_hist[i, 5] != move_hist[i, 4]){
    moves_made$movements[i] = moves_made$movements[i] + 1
    moves_made$distance[i] = moves_made$distance[i] + distGeo(p1 = island_areas[which(island_areas$area_id %in% move_hist[i,5]) ,c("lon", "lat")], p2 = island_areas[which(island_areas$area_id %in% move_hist[i,4]) ,c("lon", "lat")])/1000
  } 
  if (!is.na(move_hist[i,4]) & move_hist[i, 4] != move_hist[i, 3]) {
    moves_made$movements[i] =  moves_made$movements[i] + 1
    moves_made$distance[i] = moves_made$distance[i] + distGeo(p1 = island_areas[which(island_areas$area_id %in% move_hist[i,4]) ,c("lon", "lat")], p2 = island_areas[which(island_areas$area_id %in% move_hist[i,3]) ,c("lon", "lat")])/1000
  } 
  if (!is.na(move_hist[i,3]) & move_hist[i, 3] != move_hist[i, 2]) {
    moves_made$movements[i] =  moves_made$movements[i] + 1
    moves_made$distance[i] = moves_made$distance[i] + distGeo(p1 = island_areas[which(island_areas$area_id %in% move_hist[i,3]) ,c("lon", "lat")], p2 = island_areas[which(island_areas$area_id %in% move_hist[i,2]) ,c("lon", "lat")])/1000
  } 
  if (!is.na(move_hist[i,2]) & move_hist[i, 2] != move_hist[i, 1]) {
    moves_made$movements[i] =  moves_made$movements[i] + 1
    if(dim(island_areas[which(island_areas$area_id %in% move_hist[i,2]) ,c("lon", "lat")])[1] > 0){
      moves_made$distance[i] = moves_made$distance[i] + distGeo(p1 = island_areas[which(island_areas$area_id %in% move_hist[i,2]) ,c("lon", "lat")], p2 = island_areas[which(island_areas$area_id %in% move_hist[i,1]) ,c("lon", "lat")])/1000
    }
}
}

unique(move_hist[which(move_hist[,1] %in% island_areas$area_id == FALSE),1])
unique(move_hist[which(move_hist[,2] %in% island_areas$area_id == FALSE),2])
unique(move_hist[which(move_hist[,3] %in% island_areas$area_id == FALSE),3])
unique(move_hist[which(move_hist[,4] %in% island_areas$area_id == FALSE),4])
unique(move_hist[which(move_hist[,5] %in% island_areas$area_id == FALSE),5])



## How many fish were recaptured?
length(moves_made$movements)
  # 487
## How many of those fish did not move?
length(which(moves_made$movements == 0))
  # 352
## How many made 1 movement?
length(which(moves_made$movements == 1))
  # 132
## How many moved twice?
length(which(moves_made$movements == 2))
  # 3
## What is the percentage of fish that moved
  round(length(which(moves_made$movements > 0)) / length(moves_made$movements), digits = 3)
    # 0.277
  
## Classifying movements as adjacent, intra_island, or inter_island

# How far do they move?
movement_index = which(moves_made$movements > 0)
moves_made$classification = as.numeric(NA)
for(i in movement_index){
  print(move_hist[i, ])
  moves_made$classification[i] = readline("(1) adjacent, (2) intra_island, (3) inter_island")
}


moves_made$classification[moves_made$classification == 1] = "adjacent"
moves_made$classification[moves_made$classification == 2] = "intra-island"
moves_made$classification[moves_made$classification == 3] = "inter-island"
moves_made$classification[is.na(moves_made$classification)] = "none"
write.csv(moves_made, file.path(data_dir, 'movement_classification.csv'), row.names = FALSE)
moves_made = read.csv(file.path(data_dir, 'movement_classification.csv'))

movement_summary = aggregate(moves_made$tag_id, by = list(moves_made$classification), FUN = uniqueN)
  colnames(movement_summary) = c('classification', 'frequency')
  ## Reordering before making barchart
  movement_summary = movement_summary[c(4, 1, 3, 2), ]
  ## Converting frequency to percentage
  movement_summary$percentage = movement_summary$frequency / sum(movement_summary$frequency) * 100
  
png("Movement of Recaptured Fish - Barchart.png",  width = 826, height = 216)
  barplot(height = movement_summary$percentage, names.arg = movement_summary$classification, ylim = c(0,100), main = "Movement of Recaptured Fish - Okamoto", ylab = '% of Recaptures')
  text(.7, 8 + movement_summary$percentage[movement_summary$classification == "none"], labels = movement_summary$frequency[movement_summary$classification == "none"])
  text(1.9, 8 + movement_summary$percentage[movement_summary$classification == "adjacent"], labels =  movement_summary$frequency[movement_summary$classification == "adjacent"])
  text(3.1, 8 + movement_summary$percentage[movement_summary$classification == "intra-island"], labels = movement_summary$frequency[movement_summary$classification == "intra-island"])
  text(4.3, 8 + movement_summary$percentage[movement_summary$classification == "inter-island"], labels = movement_summary$frequency[movement_summary$classification == "inter-island"])
dev.off()

# Do multiple recaptures move?
# Does movement likelihood coorospond to DAL?
# Geodesic movement


