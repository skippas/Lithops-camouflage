# IS THE ORDER OF THE COLUMNS GOING INTO COLDIST MATCHED TO THE CONE PROPORTION ORDER?
# IS THE LUM COLUMN ASSUMED TO BE LAST?

##### Image level dists ##### 
#stop emphasising geometric means in variables, its just a method. focus is on image-level dists vs overall. or nearest vs avg

# all data
gmean <- function(x) {
  exp(mean(log(x)))
}

pfowGm<- pfow %>% group_by(abbrevs, mspec, substrate) %>%
  summarise(across(contains('Mean'), gmean))

formatPavoColdist <- function(df, visSystem){
  df<- df[order(df$substrate),]
  df<- df %>% unite('abbrevs__mspec', abbrevs, mspec, sep = '__') %>%
    as.data.frame()
  rowNamez<- paste(df$abbrevs__mspec, df$substrate, sep = '__')
  .rowNamesDF(df, make.names = T) <- rowNamez # makes unique names for dups
  if(visSystem == 'pfowl'){
    df<- df %>% rename('vs' = vsMean, 's' = swMean, 'm' = mwMean, 
                       'l' = lwMean, 'lum' = dblMean)
  }
  if(visSystem == 'human'){
    df<- df %>% rename('s' = swMean, 'm' = mwMean, 
                       'l' = lwMean, 'lum' = lumMean)
  }
  if(visSystem == 'ferret'){
    df<- df %>% rename('s' = swMean,
                       'l' = lwMean, 'lum' = lumMean)
  }
  dfList<- split(df, df$abbrevs__mspec) 
  dfList<- lapply(dfList,
                  function(x) {x <- x %>% select(!c(abbrevs__mspec, substrate))})
  return(dfList)
}

pfowGm<- formatPavoColdist(pfowGm, visSystem = 'pfowl')
pfowGm<- pfowGm[sapply(pfowGm, nrow) >= 2]

library(pavo)
conePf<- c(1, 1.9, 2.2, 2.1)
coneHu<- c(1, 5.49, 10.99)
coneFe<- c(1, 14)

pfowGmD<- lapply(pfowGm, coldist, qcatch = 'Qi', n = conePf,
               weber = .05, achromatic = T) 

bindList<- function(list, visSystem){
  bind_rows(list, .id =  'abbrevs_mspec') %>%
    separate(patch1, into = c(NA, NA, 'comparison'), sep = '__', remove = T) %>%
    separate(patch2, into = c(NA, NA, 'comparison2'), sep = '__', remove = T) %>%
    separate(abbrevs_mspec, into = c('abbrevs', 'mspec'), sep = '__') %>%
    unite(col = 'comparison', comparison, comparison2, remove = T, sep = '') %>%
    mutate(visSys = visSystem)
}
pfowGmD<- bindList(pfowGmD, visSystem = 'pfowl')

rm(pfowGm)

humGm<- hum %>% group_by(abbrevs, mspec, substrate) %>%
  summarise(across(contains('Mean'), gmean))

humGm<- formatPavoColdist(humGm, visSystem = 'human')
humGm<- humGm[sapply(humGm, nrow) >= 2]

library(pavo)
humGmD<- lapply(humGm, coldist, qcatch = 'Qi', n = coneHu,
                  weber = .05, achromatic = T) 
humGmD<- bindList(humGmD, visSystem = 'human')

rm(humGm)

ferGm<- fer %>% group_by(abbrevs, mspec, substrate) %>%
  summarise(across(contains('Mean'), gmean))

ferGm<- formatPavoColdist(ferGm, visSystem = 'ferret')
ferGm<- ferGm[sapply(ferGm, nrow) >= 2] # much more elegent to get rid of low comparisons than I used in nearest dists

library(pavo)
ferGmD<- lapply(ferGm, coldist, qcatch = 'Qi', n = coneFe,
                  weber = .05, achromatic = T) 
ferGmD<- bindList(ferGmD, visSystem = 'ferret')

rm(ferGm)


# nearest dists

set.seed(8) # should this be inside function to ensure same sample every time?
lithSample<- function(df){
  df<- df %>% filter(substrate == 'l') %>%
    group_by(abbrevs, mspec) %>% slice_sample(n = 1)
}

pfND<- pfow %>% ungroup() %>% 
  filter(substrate != 'l') %>% rbind(., lithSample(pfow)) %>%
  select(contains('mean'),'substrate', 'abbrevs', 'mspec') 
  
pfND<- formatPavoColdist(pfND, visSystem = 'pfowl')
pfND<- pfND[sapply(pfND, nrow) >= 2] # hope which mspecs dont have soil can still be tracked

pfND<- lapply(pfND, coldist, qcatch = 'Qi', n = conePf,
                weber = .05, achromatic = T, subset = '__l')

pfND<- bindList(pfND, visSystem = 'pfowl')

cleanComparison<- function(df){
  df$comparison<- gsub("[0-9]+","",df$comparison)
  df$comparison<- gsub("\\.","", df$comparison)
  return(df)
}
pfND<- cleanComparison(pfND)

pfND<- pfND %>% pivot_longer(cols = c('dS', 'dL'), names_to = 'visInfo',
                                   values_to = 'distance') %>%
  group_by(abbrevs, mspec, visInfo, comparison) %>%
  slice_min(distance) 


huND<- hum %>% ungroup() %>% filter(substrate != 'l') %>% 
  rbind(., lithSample(hum)) %>%
  select(contains('mean'), 'substrate', 'abbrevs', 'mspec')

huND<- formatPavoColdist(huND, visSystem = 'human')
huND<- huND[sapply(huND, nrow) >= 2]

huND<- lapply(huND, coldist, qcatch = 'Qi', n = coneHu,
              weber = .05, achromatic = T, subset = '__l')

huND<- bindList(huND, visSystem = 'human')

huND<- cleanComparison(huND) %>% 
  pivot_longer(cols = c('dS', 'dL'), names_to = 'visInfo',
               values_to = 'distance') %>%
  group_by(abbrevs, mspec, visInfo, comparison) %>%
  slice_min(distance) 

feND<- fer %>% ungroup() %>% 
  filter(substrate != 'l') %>% rbind(., lithSample(fer)) %>%
  select(contains('mean'), 'substrate', 'abbrevs', 'mspec')
feND<- fer %>% ungroup() %>% select(contains('mean'), 'substrate', 'abbrevs', 'mspec')
feND<- formatPavoColdist(feND, visSystem = 'ferret')
feND<- feND[sapply(feND, nrow) >= 2]

feND<- lapply(feND, coldist, qcatch = 'Qi', n = coneFe,
              weber = .05, achromatic = T, subset = '__l')

feND<- bindList(feND, visSystem = 'ferret')

feND<- cleanComparison(feND) %>% 
  pivot_longer(cols = c('dS', 'dL'), names_to = 'visInfo',
               values_to = 'distance') %>%
  group_by(abbrevs, mspec, visInfo, comparison) %>%
  slice_min(distance) 

#### Overall distances ####
# avg data
pfOv<- pfow %>% group_by(abbrevs, mspec, substrate) %>%
  summarise(across(contains('Mean'), mean)) %>% # mean instead of gmean initially
  group_by(abbrevs, substrate) %>%
  summarise(across(contains('Mean'), gmean))

formatPavoColdist <- function(df, visSystem){
  if(visSystem == 'pfowl'){
    df<- df %>% rename('vs' = vsMean, 's' = swMean, 'm' = mwMean, 
                       'l' = lwMean, 'lum' = dblMean)
  }
  if(visSystem == 'human'){
    df<- df %>% rename('s' = swMean, 'm' = mwMean, 
                       'l' = lwMean, 'lum' = lumMean)
  }
  if(visSystem == 'ferret'){
    df<- df %>% rename('s' = swMean,
                       'l' = lwMean, 'lum' = lumMean)
  }
  df<- as.data.frame(df)
  df<- df[order(df$substrate),]
  rowNamez<- paste(df$abbrevs, df$substrate, sep = '__')
  .rowNamesDF(df) <- rowNamez # makes unique names for dups
  dfList<- split(df, df$abbrevs) 
  dfList<- lapply(dfList,
                  function(x) {x <- x %>% ungroup() %>%
                    select(!c(abbrevs, substrate))})
  return(dfList)
}

pfOv<- formatPavoColdist(pfOv, 'pfowl')
pfOv<- pfOv[sapply(pfOv, nrow) >= 2]

pfOv<- lapply(pfOv, coldist, qcatch = 'Qi', n = conePf,
              weber = .05, achromatic = T)

bindList<- function(list, visSystem){
  bind_rows(list, .id =  'abbrevs') %>%
    separate(patch1, into = c(NA, 'comparison'), sep = '__', remove = T) %>%
    separate(patch2, into = c(NA, 'comparison2'), sep = '__', remove = T) %>%
    separate(abbrevs, into = c('abbrevs'), sep = '__') %>%
    unite(col = 'comparison', comparison, comparison2, remove = T, sep = '') %>%
    mutate(visSys = visSystem)
}
pfOv<- bindList(pfOv, visSystem = 'pfowl')

huOv<- hum %>% group_by(abbrevs, mspec, substrate) %>% 
  summarise(across(contains('Mean'), mean)) %>% 
  group_by(abbrevs, substrate) %>%
  summarise(across(contains('Mean'), gmean))

huOv<- formatPavoColdist(huOv, 'human')
huOv<- huOv[sapply(huOv, nrow) >= 2]

huOv<- lapply(huOv, coldist, n = coneHu, qcatch = 'Qi', weber = .05,
              achromatic = T)
huOv<- bindList(huOv, 'human')

feOv<- fer %>% group_by(abbrevs, mspec, substrate) %>% 
  summarise(across(contains('Mean'), mean)) %>% 
  group_by(abbrevs, substrate) %>%
  summarise(across(contains('Mean'), gmean))

feOv<- formatPavoColdist(feOv, 'ferret')
feOv<- feOv[sapply(feOv, nrow) >= 2]

feOv<- lapply(feOv, coldist, n = coneFe, qcatch = 'Qi', weber = .05,
       achromatic = T)
feOv<- bindList(feOv, 'ferret')

# nearest data, quite a complex procedure, many potential pitfalls
# eg. filtering out duplicate lithops when going from coldist format to
# normal format
set.seed(8) # should this be inside function to ensure same sample every time?
lithSample<- function(df){
  df<- df %>% filter(substrate == 'l') %>%
    group_by(abbrevs, mspec) %>% slice_sample(n = 1)
}

pfOvN<- pfow %>% ungroup() %>% 
  filter(substrate != 'l') %>%
  rbind(., lithSample(pfow)) %>%
  select('vsMean', 'swMean', 'mwMean', 'lwMean', 'dblMean',
         'substrate', 'abbrevs', 'mspec', 'id') 

formatPavoColdist <- function(df, visSystem){
  df<- df[order(df$substrate),]
  df<- df %>% unite('abbrevs__mspec', abbrevs, mspec, sep = '__') %>%
    as.data.frame()
  rowNamez<- paste(df$abbrevs__mspec, df$id, df$substrate, sep = '__')
  .rowNamesDF(df, make.names = T) <- rowNamez 
  if(visSystem == 'pfowl'){
    df<- df %>% rename('vs' = vsMean, 's' = swMean, 'm' = mwMean, 
                       'l' = lwMean, 'lum' = dblMean)
  }
  if(visSystem == 'human'){
    df<- df %>% rename('s' = swMean, 'm' = mwMean, 
                       'l' = lwMean, 'lum' = lumMean)
  }
  if(visSystem == 'ferret'){
    df<- df %>% rename('s' = swMean,
                       'l' = lwMean, 'lum' = lumMean)
  }
  dfList<- split(df, df$abbrevs__mspec) 
  dfList<- lapply(dfList,
                  function(x) {x <- x %>% 
                    select(!c(abbrevs__mspec, substrate, id))})
  return(dfList)
}
pfOvN<- formatPavoColdist(pfOvN, 'pfowl')
pfOvN<- pfOvN[sapply(pfOvN, nrow) >= 2]

pfOvN<- lapply(pfOvN, coldist, qcatch = 'Qi', n = conePf,
              weber = .05, achromatic = T, subset = '__l')

bindList<- function(list, visSystem){
  bind_rows(list, .id =  'abbrevs__mspec') %>%
    separate(abbrevs__mspec, into = c('abbrevs', 'mspec'), sep = '__') %>%
    separate(patch1, into = c(NA, NA, 'id1', 'substrate1'), 
             sep = '__', remove = T) %>%
    separate(patch2, into = c(NA, NA, 'id2', 'substrate2'),
             sep = '__', remove = T) %>%
    unite(col = 'comparison', substrate1, substrate2, remove = T, sep = '') %>% 
    mutate(visSys = visSystem)
}
pfOvN<- bindList(pfOvN, 'pfowl')

pfOvN<- pfOvN %>% pivot_longer(cols = c('dS', 'dL'), names_to = 'visInfo',
                             values_to = 'distance') %>%
  group_by(abbrevs, mspec, visInfo, comparison) %>%
  slice_min(distance) 

indexN<- function(df, visSystem){
  if(visSystem == 'pfowl'){nrstIndex<- pfow}
  if(visSystem == 'human'){nrstIndex<- hum}
  if(visSystem == 'ferret'){nrstIndex<- fer}
  nrstIndex<- nrstIndex  %>% ungroup %>% 
    select(c(contains('mean'), abbrevs, mspec, id, substrate))
  nrstIndexL<- merge(nrstIndex, 
                     df[, c('id1', 'visInfo')],
                     by.x = 'id',
                     by.y = 'id1') %>%
    distinct()
  
  nrstIndexRS<- merge(nrstIndex,
                      df[, c('id2', 'visInfo')],
                      by.x = 'id',
                      by.y = 'id2')
  nrstIndexLRS<- rbind(nrstIndexL, nrstIndexRS)

}
pfOvN<- indexN(pfOvN, visSystem =  'pfowl')

pfOvN<- pfOvN %>% group_by(visInfo, abbrevs, substrate) %>%
  summarise(across(contains('Mean'), gmean))

formatAfterGmean <- function(df, visSystem){
  df<- df %>% relocate(any_of(c('vsMean', 'swMean', 'mwMean',
                       'lwMean', 'lumMean', 'dblMean')))
  if(visSystem == 'pfowl'){
    df<- df %>% rename('vs' = vsMean, 's' = swMean, 'm' = mwMean, 
                       'l' = lwMean, 'lum' = dblMean)
  }
  if(visSystem == 'human'){
    df<- df %>% rename('s' = swMean, 'm' = mwMean, 
                       'l' = lwMean, 'lum' = lumMean)
  }
  if(visSystem == 'ferret'){
    df<- df %>% rename('s' = swMean,
                       'l' = lwMean, 'lum' = lumMean)
  }
  df<- as.data.frame(df)
  df<- df[order(df$substrate),]
  dfList<- split(df, list(df$abbrevs, df$visInfo), sep = '__')
  dfList<- lapply(dfList, remove_rownames)
  dfList<- lapply(dfList, column_to_rownames, 'substrate')
  dfList<- lapply(dfList,
                  function(x) {x <- x %>% ungroup() %>%
                    select(!c(abbrevs, visInfo))})

}

pfOvN<- formatAfterGmean(pfOvN, 'pfowl')
pfOvN<- pfOvN[sapply(pfOvN, nrow) >= 2]

pfOvN<- lapply(pfOvN, coldist, qcatch = 'Qi', n = conePf,
              weber = .05, achromatic = T)

bindAfterGmean<- function(list, visSystem){
  bind_rows(list, .id =  'abbrevs__visInfo') %>%
    separate(abbrevs__visInfo, into = c('abbrevs', 'visInfo'), sep = '__') %>%
    unite(col = 'comparison', patch1, patch2, remove = T, sep = '') %>%
    mutate(visSys = visSystem)
}

pfOvN<- bindAfterGmean(pfOvN, 'pfowl')

# human

huOvN<- hum %>% ungroup() %>% 
  filter(substrate != 'l') %>%
  rbind(., lithSample(hum)) %>%
  select('swMean', 'mwMean', "lwMean", 'lumMean',
         'substrate', 'abbrevs', 'mspec', 'id') 

huOvN<- formatPavoColdist(huOvN, visSystem = 'human')
huOvN<- lapply(huOvN, coldist, n= coneHu, qcatch = 'Qi',
               weber = .05, achromatic = T, subset = '__l')
huOvN<- bindList(huOvN, 'human')

huOvN<- huOvN %>% pivot_longer(cols = c(dS, dL), names_to = 'visInfo',
                               values_to = 'distance') %>%
  group_by(abbrevs, mspec, visInfo, comparison) %>%
  slice_min(distance) 

huOvN<- indexN(huOvN, visSystem = 'human')
huOvN<- huOvN %>% group_by(abbrevs, visInfo, substrate) %>%
  summarise(across(contains('Mean'), gmean))

huOvN<- formatAfterGmean(huOvN, 'human')
huOvN<- lapply(huOvN, coldist, n = coneHu, qcatch = 'Qi',  
       weber = .05, achromatic = T) 
huOvN<- bindAfterGmean(huOvN, 'human')

# ferret
feOvN<- fer %>% ungroup %>%
  filter(substrate != 'l') %>%
  rbind(., lithSample(fer)) %>%
  select('swMean', 'lwMean', 'lumMean',
         'substrate', 'abbrevs', 'mspec', 'id')

feOvN<- formatPavoColdist(feOvN, 'ferret')
feOvN<- lapply(feOvN, coldist, n = coneFe, qcatch = 'Qi',
               weber = .05, achromatic = T, subset = '__l' )

feOvN<- bindList(feOvN, 'ferret')

feOvN<- feOvN %>% pivot_longer(cols = c(dS, dL), values_to = 'distance',
                               names_to = 'visInfo') %>%
  group_by(abbrevs, mspec, visInfo, comparison) %>%
  slice_min(distance)

feOvN<- indexN(feOvN, 'ferret')

feOvN<- feOvN %>% group_by(abbrevs, visInfo, substrate) %>%
  summarise(across(contains('Mean'), gmean))

feOvN<- formatAfterGmean(feOvN, 'ferret')
feOvN<- lapply(feOvN, coldist, n = coneFe, qcatch = 'Qi',
               weber = .05, achromatic = T, subset = 'l' )

feOvN<- bindAfterGmean(feOvN, visSystem = 'ferret')



