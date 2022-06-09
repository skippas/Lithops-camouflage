nearestDistances<- function(df, visSys){

df<- df[order(df$substrate),]
# choose 1 lithops
# split rock and soils
# split into list of dfs by mspec
# filter dfs (mspecs) that dont have at least 1 lithops and 1 substrate
set.seed(8)
lithSample<- df %>% filter(substrate == "l") %>% group_by(abbrevs, mspec) %>% 
  sample_n(size = 1) 

lr<- rbind(lithSample, df[df$substrate == "r",])

lr<- lr %>% unite('abbrevs__mspec', abbrevs, mspec, sep = '__') 
lr<- split(lr, lr$abbrevs__mspec) 

ls<- rbind(lithSample, df[df$substrate != "l",])
ls<- ls[ls$substrate != 'r',]
subset<- ls %>% group_by(abbrevs, mspec) %>%
  summarise(n_distinct(substrate)) %>%
  filter(`n_distinct(substrate)` <2) %>% ungroup %>% select(mspec)
ls<- ls[! ls$mspec %in% subset$mspec,]
ls<- ls %>% unite('abbrevs__mspec', abbrevs, mspec, sep = '__') 
ls<- split(ls, ls$abbrevs__mspec) 

cold4matPf <- function(df){
  rowNamez <- paste(df$abbrevs__mspec, df$substrate, df$roi, sep = '__')
  df <- as.data.frame(df)
  df <- df[, c("vsMean", "swMean", "mwMean", "lwMean","lumMean")]
  .rowNamesDF(df, make.names = T) <- rowNamez # makes unique names for dups
  colnames(df) <- c("vs", "s", "m", "l", "lum") 
  return(df)
}
cold4matHu <- function(df){
  rowNamez <- paste(df$abbrevs__mspec, df$substrate, df$roi, sep = '__')
  df <- as.data.frame(df)
  df <- df[, c("swMean", "mwMean", "lwMean","lumMean")]
  .rowNamesDF(df, make.names = T) <- rowNamez # makes unique names for dups
  colnames(df) <- c( "s", "m", "l", "lum") 
  return(df)
}
cold4matFe <- function(df){
  rowNamez <- paste(df$abbrevs__mspec, df$substrate, df$roi, sep = '__')
  df <- as.data.frame(df)
  df <- df[, c("swMean", "lwMean","lumMean")]
  .rowNamesDF(df, make.names = T) <- rowNamez # makes unique names for dups
  colnames(df) <- c( "s", "l", "lum") 
  return(df)
}

if(visSys == 'pfowl'){
  lr<- lapply(lr, cold4matPf)
  lr<- lapply(lr, coldist, qcatch = 'Qi', n = c(1, 1.9, 2.2, 2.1),
               weber = .05, achromatic = T, subset = c('__l', '__r'))
  ls<- lapply(ls, cold4matPf)
  ls<- lapply(ls, coldist, qcatch = 'Qi', n = c(1, 1.9, 2.2, 2.1),
              weber = .05, achromatic = T, subset = c('__l', '__s'))
  }

if(visSys == 'human'){
  lr<- lapply(lr, cold4matHu)
  lr<- lapply(lr, coldist, qcatch = 'Qi', n = c(1, 5.49, 10.99),
              weber = .05, achromatic = T, subset = c('__l', '__r'))
  ls<- lapply(ls, cold4matHu)
  ls<- lapply(ls, coldist, qcatch = 'Qi', n = c(1, 5.49, 10.99),
              weber = .05, achromatic = T, subset = c('__l', '__s'))
}

if(visSys == 'ferret'){
  lr<- lapply(lr, cold4matFe)
  lr<- lapply(lr, coldist, qcatch = 'Qi', n = c(1, 14),
              weber = .05, achromatic = T, subset = c('__l', '__r'))
  ls<- lapply(ls, cold4matFe)
  ls<- lapply(ls, coldist, qcatch = 'Qi', n = c(1, 14),
              weber = .05, achromatic = T, subset = c('__l', '__s'))
}

dists<- c(ls, lr)
dists<- bind_rows(dists, .id =  'abbrevs_mspec') %>% 
  separate(patch1, into = c(NA, NA, 'comparison','roi.l'), sep = '__', remove = T) %>%
  separate(patch2, into = c(NA, NA, 'comparison2', 'roi.bg'), sep = '__', remove = T) %>%
  separate(abbrevs_mspec, into = c('abbrevs', 'mspec'), sep = '__') %>%
  unite(col = 'comparison', comparison, comparison2, remove = T, sep = '') 

  dists$comparison<- gsub("[0-9]+","",dists$comparison)
   dists$comparison<- gsub("\\.","", dists$comparison)
  return(dists)
  
  
}
