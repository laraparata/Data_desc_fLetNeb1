#make global lethirnus nebulosus distribution map

install.packages("leaflet")
install.packages('rgbif')


library(leaflet)
library(rgbif)
library(dplyr)
library(robis)


# Retrieve occurrence data from GBIF
species_name <- "Lethrinus nebulosus"
occurrences <- occ_search(scientificName = species_name)  # Adjust limit as needed

# Extract the data
occ_data <- occurrences$data

# Keep only the records with latitude and longitude information
occ_data <- occ_data %>%
  filter(!is.na(occ_data$decimalLatitude),
         !is.na(occ_data$decimalLongitude))

# Retrieve occurrence data from OBIS
obis_occurrences <- occurrence(scientificname = species_name)

# Keep only the records with latitude and longitude information
obis_data <- obis_occurrences %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude))

# Create the interactive map using `leaflet`
leaflet() %>%
  addTiles() %>%
  addCircleMarkers(
    data = occ_data,
    ~ occ_data$decimalLongitude,
    ~ occ_data$decimalLatitude,
    popup = ~ paste(
      "<strong>Scientific Name:</strong>",
      scientificName,
      "<br>",
      "<strong>Country:</strong>",
      country,
      "<br>",
      "<strong>Event Date:</strong>",
      eventDate
    ),
    radius = 4,
    stroke = FALSE,
    fillOpacity = 0.5,
    color = 'red'
  ) %>%
  addCircleMarkers(
    data = obis_data,
    ~ decimalLongitude,
    ~ decimalLatitude,
    popup = ~ paste(
      "<strong>Scientific Name:</strong>",
      scientificName,
      "<br>",
      "<strong>Country:</strong>",
      country,
      "<br>",
      "<strong>Event Date:</strong>",
      eventDate
    ),
    radius = 4,
    stroke = FALSE,
    fillOpacity = 0.5,
    color = 'blue'
  )
