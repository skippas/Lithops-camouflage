pfow<- pfow %>% unite('abbrevs__mspec', abbrevs, mspec, sep = '__') 
test<- pfow[pfow$pop_spp == as.character(unique(pfow$pop_spp)[1]),]
test<- test %>% filter(substrate == 'l')
ls<- split(test, test$abbrevs__mspec)
ls<- lapply(ls, rbind, pfow[pfow$substrate == 'r' & 
                              pfow$pop_spp == unique(pfow$pop_spp)[1] ,])


formatPavoColdist <- function(df){
  rowNamez <- paste(df$abbrevs__mspec, df$substrate, sep = '__')
  df <- as.data.frame(df)
  df <- df[, c("vsMean", "swMean", "mwMean", "lwMean","lumMean")]
  .rowNamesDF(df, make.names = T) <- rowNamez # makes unique names for dups
  colnames(df) <- c("vs", "s", "m", "l", "lum") 
  return(df)
}
ls<- lapply(ls, formatPavoColdist)
ls<- ls[sapply(ls, nrow) >= 2]

library(pavo)
dists<- lapply(ls, coldist, qcatch = 'Qi', n = c(1, 1.9, 2.2, 2.1),
               weber = .05, achromatic = T) 
dists<- bind_rows(dists, .id =  'abbrevs_mspec') %>% 
  separate(patch1, into = c(NA, NA, 'comparison'), sep = '__', remove = T) %>%
  separate(patch2, into = c(NA, NA, 'comparison2'), sep = '__', remove = T) %>%
  separate(abbrevs_mspec, into = c('abbrevs', 'mspec'), sep = '__') %>%
  unite(col = 'comparison', comparison, comparison2, remove = T, sep = '') 
