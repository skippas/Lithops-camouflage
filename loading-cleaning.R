library(tidyverse)
library(broom)

# cleaning #
#* pfowl ####
pfow<- read.csv('data//Nikon D7100 CoastalOpt 105mm D65 to Peafowl D65.csv')
pfow<- pfow[, ! names(pfow) %in% c('X', 'lumMean')]
pfow <- separate(pfow, Label, into = c("pop_spp", "label"), sep = "/")
pfow <- separate(pfow, label, into = c("vis_image","uv_image", "substrate"), sep = "_")
pfow <- unite(pfow, "mspec", c("vis_image","uv_image"))
pfow$roi <- pfow$substrate
pfow$substrate <- substr(pfow$substrate, 0,1)
unique(pfow$substrate) # empty space substrates?
pfow<- pfow[pfow$substrate %in% c("a","b","c"),]
pfow$substrate[pfow$substrate == 'a']<- 'l'
pfow$substrate[pfow$substrate == 'b']<- 'r'
pfow$substrate[pfow$substrate == 'c']<- 's'
pfow$group <- "bg"
pfow$group[pfow$substrate == "l"] <- "l"

pfow[pfow$pop_spp == "namies_olivacea_x","pop_spp"] <- "namies_olivacea_XXON"
pfow[pfow$pop_spp == "pofadder_olivacea_x","pop_spp"] <- "pofadder_olivacea_XXOP"
pfow[pfow$pop_spp == "Liebeberg_hookeri_x","pop_spp"] <- "Liebeberg_hookeri_XXHL"
pfow[pfow$pop_spp == "Hopetown_aucampiae_x","pop_spp"] <- "Hopetown_aucampiae_XXAH"
pfow[pfow$pop_spp == "Haasriver2_Oztenia_JKL","pop_spp"] <- "Haasriver2_otzeniana_JKL"
pfow[pfow$pop_spp == "Haasriver1_Oztenia_IJK","pop_spp"] <- "Haasriver1_otzeniana_IJK"

lut<- read.csv('data//lut.csv')
pfow<- merge(pfow, lut, by = 'pop_spp')

pfow[pfow$species == "lesleii.venterii", "species"] <- "leslei.venteri"

cols <- sapply(pfow, is.numeric)
pfow<- pfow[! rowSums(pfow[cols] < 0) > 0, ]
rm(cols)

pfow$pop_spp <- as.factor(pfow$pop_spp)
pfow$abbrevs <- as.factor(pfow$abbrevs)

# Pops with low sample sizes
pfow <- pfow[!pfow$pop_spp %in% c("kangnas_marmorata_LM","khoeries_fulleri_KL"),]
# drop unused levels 
pfow[] <- lapply(pfow, function(x) if(is.factor(x)) factor(x) else x)

# Remove colour oultiers

# Lump subspecies / variants
pfow <- pfow %>% 
  mutate(species = fct_relevel(species, sort)) %>%
  group_by(species) %>%
  mutate(collapseSpp = case_when(
    n_distinct(abbrevs) >1 ~ as.character(species),
    n_distinct(abbrevs) <2 ~ "singlePopSpp")) %>%
  mutate(collapseSpp = as_factor(collapseSpp),
         collapseSpp = fct_recode(collapseSpp,
                                  bromf. = 'bromfieldiiSpp',
                                  compt. = 'comptoniiSpp',
                                  verruc. = 'verruculosa',
                                  dint. = 'dinteriSpp',
                                  singlePops = 'singlePopSpp'))

# create unique ID (necessary for accurately indexing original df with dist df)
pfow<- rownames_to_column(pfow, var = 'id')

# put in a table with abbrev, pop, species, sampling of subs..
# TO DO: add in the mspec counts. Difficult because summarise causes to lose. 

#* human ####
hum<- read.csv('data//Nikon D7100 CoastalOpt 105mm D65 to Human D65_all.csv')
hum<- hum[, ! names(hum) %in% 'X']
hum <- separate(hum, Label, into = c("pop_spp", "label"), sep = "/")
hum <- separate(hum, label, into = c("vis_image","uv_image", "substrate"), sep = "_")
hum <- unite(hum, "mspec", c("vis_image","uv_image"))
hum$roi <- hum$substrate
hum$substrate <- substr(hum$substrate, 0,1)
unique(hum$substrate) # empty space substrates?
hum<- hum[hum$substrate %in% c("a","b","c"),]
hum$substrate[hum$substrate == 'a']<- 'l'
hum$substrate[hum$substrate == 'b']<- 'r'
hum$substrate[hum$substrate == 'c']<- 's'
hum$group <- "bg"
hum$group[hum$substrate == "l"] <- "l"

hum[hum$pop_spp == "namies_olivacea_x","pop_spp"] <- "namies_olivacea_XXON"
hum[hum$pop_spp == "pofadder_olivacea_x","pop_spp"] <- "pofadder_olivacea_XXOP"
hum[hum$pop_spp == "Liebeberg_hookeri_x","pop_spp"] <- "Liebeberg_hookeri_XXHL"
hum[hum$pop_spp == "Hopetown_aucampiae_x","pop_spp"] <- "Hopetown_aucampiae_XXAH"
hum[hum$pop_spp == "Haasriver2_Oztenia_JKL","pop_spp"] <- "Haasriver2_otzeniana_JKL"
hum[hum$pop_spp == "Haasriver1_Oztenia_IJK","pop_spp"] <- "Haasriver1_otzeniana_IJK"

lut<- read.csv('data//lut.csv')
hum<- merge(hum, lut, by = 'pop_spp')

hum[hum$species == "lesleii.venterii", "species"] <- "leslei.venteri"

cols <- sapply(hum, is.numeric)
hum<- hum[! rowSums(hum[cols] < 0) > 0, ]
rm(cols)
hum %>% filter(lumMean < 0)
hum$pop_spp <- as.factor(hum$pop_spp)
hum$abbrevs <- as.factor(hum$abbrevs)

# Pops with low sample sizes
hum <- hum[!hum$pop_spp %in% c("kangnas_marmorata_LM","khoeries_fulleri_KL"),]
# drop unused levels 
hum[] <- lapply(hum, function(x) if(is.factor(x)) factor(x) else x)

# Remove colour oultiers

# Lump subspecies / variants
hum <- hum %>% 
  mutate(species = fct_relevel(species, sort)) %>%
  group_by(species) %>%
  mutate(collapseSpp = case_when(
    n_distinct(abbrevs) >1 ~ as.character(species),
    n_distinct(abbrevs) <2 ~ "singlePopSpp")) %>%
  mutate(collapseSpp = as_factor(collapseSpp),
         collapseSpp = fct_recode(collapseSpp,
                                  bromf. = 'bromfieldiiSpp',
                                  compt. = 'comptoniiSpp',
                                  verruc. = 'verruculosa',
                                  dint. = 'dinteriSpp',
                                  singlePops = 'singlePopSpp'))

# put in a table with abbrev, pop, species, sampling of subs..
# TO DO: add in the mspec counts. Difficult because summarise causes to lose. 

# create id
hum<- rownames_to_column(hum, var = 'id')

#* ferret ####
fer<- read.csv('data//ferretdataedited_peafowlremoved.csv')
fer<- fer[, ! names(fer) %in% 'X']
fer <- separate(fer, Label, into = c("pop_spp", "label"), sep = "/")
fer <- separate(fer, label, into = c("vis_image","uv_image", "substrate"), sep = "_")
fer <- unite(fer, "mspec", c("vis_image","uv_image"))
fer$roi <- fer$substrate
fer$substrate <- substr(fer$substrate, 0,1)
unique(fer$substrate) # empty space substrates?
fer<- fer[fer$substrate %in% c("a","b","c"),]
fer$substrate[fer$substrate == 'a']<- 'l'
fer$substrate[fer$substrate == 'b']<- 'r'
fer$substrate[fer$substrate == 'c']<- 's'
fer$group <- "bg"
fer$group[fer$substrate == "l"] <- "l"

fer[fer$pop_spp == "namies_olivacea_x","pop_spp"] <- "namies_olivacea_XXON"
fer[fer$pop_spp == "pofadder_olivacea_x","pop_spp"] <- "pofadder_olivacea_XXOP"
fer[fer$pop_spp == "Liebeberg_hookeri_x","pop_spp"] <- "Liebeberg_hookeri_XXHL"
fer[fer$pop_spp == "Hopetown_aucampiae_x","pop_spp"] <- "Hopetown_aucampiae_XXAH"
fer[fer$pop_spp == "Haasriver2_Oztenia_JKL","pop_spp"] <- "Haasriver2_otzeniana_JKL"
fer[fer$pop_spp == "Haasriver1_Oztenia_IJK","pop_spp"] <- "Haasriver1_otzeniana_IJK"

lut<- read.csv('data//lut.csv')
fer<- merge(fer, lut, by = 'pop_spp')

fer[fer$species == "lesleii.venterii", "species"] <- "leslei.venteri"

cols <- sapply(fer, is.numeric)
fer<- fer[! rowSums(fer[cols] < 0) > 0, ]
rm(cols)

fer$pop_spp <- as.factor(fer$pop_spp)
fer$abbrevs <- as.factor(fer$abbrevs)

# create id
fer<- rownames_to_column(fer, var = 'id')