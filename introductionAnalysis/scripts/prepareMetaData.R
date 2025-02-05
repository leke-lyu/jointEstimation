library(magrittr)
library(dplyr)
library(zipcodeR)
library(stringr)

load("clusters.RData")
load("data/houstonDeltaMeta.RData")
zipCodeTable <- search_state('TX') %>% select(zipcode, county)
zipCodeTable[zipCodeTable$zipcode == "77328", "county"] <- "Montgomery County"

tips <- tipList %>% unlist(.)
filteredHoustonDeltaMeta <- houstonDeltaMeta %>% filter(GISAID_name %in% tips) 
filteredHoustonDeltaMeta <- merge(filteredHoustonDeltaMeta, zipCodeTable, by.x = "zip", by.y = "zipcode", all.x = TRUE)

filteredHoustonDeltaMeta <- filteredHoustonDeltaMeta %>%
  mutate(new_age_groups = case_when(
    is.na(age) ~ 'Unknown',
    age >= 0 & age <= 12 ~ 'Infants and Children',
    age >= 13 & age <= 19 ~ 'Teenagers',
    age >= 20 & age <= 35 ~ 'Young Adults',
    age >= 36 & age <= 55 ~ 'Middle-aged Adults',
    age >= 56 ~ 'Seniors',
    TRUE ~ 'Unknown'
  ))
  
filteredHoustonDeltaMeta <- filteredHoustonDeltaMeta %>%
  mutate(sex = ifelse(is.na(sex), "Unknown", sex))

filteredHoustonDeltaMeta$new_age_groups %>% table()
filteredHoustonDeltaMeta$sex %>% table()
filteredHoustonDeltaMeta$county %>% table()

save(filteredHoustonDeltaMeta, file = "filteredHoustonDeltaMeta.RData")





