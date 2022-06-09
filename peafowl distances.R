source('loading-cleaning.R')
rm(lut)
source('generate-distances.R')

#### geometric means ####
#* pfowl ####
gmean <- function(x) {
  exp(mean(log(x)))
}
pfowGm<- pfow %>% group_by(abbrevs, mspec, substrate) %>%
  summarise(across(contains('Mean'), gmean))

pfowGm<- pfowGm %>% unite('abbrevs__mspec', abbrevs, mspec, sep = '__') 
pfowGm<- split(pfowGm, pfowGm$abbrevs__mspec) 

formatPavoColdist <- function(df){
  rowNamez <- paste(df$abbrevs__mspec, df$substrate, sep = '__')
  df <- as.data.frame(df)
  df <- df[, c("vsMean", "swMean", "mwMean", "lwMean","lumMean")]
  .rowNamesDF(df, make.names = T) <- rowNamez # makes unique names for dups
  colnames(df) <- c("vs", "s", "m", "l", "lum") 
  return(df)
}
pfowGm<- lapply(pfowGm, formatPavoColdist)
pfowGm<- pfowGm[sapply(pfowGm, nrow) >= 2]

library(pavo)
dists<- lapply(pfowGm, coldist, qcatch = 'Qi', n = c(1, 1.9, 2.2, 2.1),
               weber = .05, achromatic = T) 
dists<- bind_rows(dists, .id =  'abbrevs_mspec') %>% 
  separate(patch1, into = c(NA, NA, 'comparison'), sep = '__', remove = T) %>%
  separate(patch2, into = c(NA, NA, 'comparison2'), sep = '__', remove = T) %>%
  separate(abbrevs_mspec, into = c('abbrevs', 'mspec'), sep = '__') %>%
  unite(col = 'comparison', comparison, comparison2, remove = T, sep = '') 

pfowlDists<- dists
pfowlDists$visSys<- 'pfowl'
rm(pfowGm)

#* human ####
humGm<- hum %>% group_by(abbrevs, mspec, substrate) %>%
  summarise(across(contains('Mean'), gmean))

humGm<- humGm %>% unite('abbrevs__mspec', abbrevs, mspec, sep = '__') 
humGm<- split(humGm, humGm$abbrevs__mspec) 

formatPavoColdist <- function(df){
  rowNamez <- paste(df$abbrevs__mspec, df$substrate, sep = '__')
  df <- as.data.frame(df)
  df <- df[, c("swMean", "mwMean", "lwMean","lumMean")]
  .rowNamesDF(df, make.names = T) <- rowNamez # makes unique names for dups
  colnames(df) <- c("s", "m", "l", "lum") 
  return(df)
}
humGm<- lapply(humGm, formatPavoColdist)
humGm<- humGm[sapply(humGm, nrow) >= 2]

library(pavo)
humDists<- lapply(humGm, coldist, qcatch = 'Qi', n = c(1, 5.49, 10.99),
               weber = .05, achromatic = T) 
humDists<- bind_rows(humDists, .id =  'abbrevs_mspec') %>% 
  separate(patch1, into = c(NA, NA, 'comparison'), sep = '__', remove = T) %>%
  separate(patch2, into = c(NA, NA, 'comparison2'), sep = '__', remove = T) %>%
  separate(abbrevs_mspec, into = c('abbrevs', 'mspec'), sep = '__') %>%
  unite(col = 'comparison', comparison, comparison2, remove = T, sep = '') 
rm(humGm)

humDists$visSys<- 'human'
#* ferret ####
ferGm<- fer %>% group_by(abbrevs, mspec, substrate) %>%
  summarise(across(contains('Mean'), gmean))

ferGm<- ferGm %>% unite('abbrevs__mspec', abbrevs, mspec, sep = '__') 
ferGm<- split(ferGm, ferGm$abbrevs__mspec) 

formatPavoColdist <- function(df){
  rowNamez <- paste(df$abbrevs__mspec, df$substrate, sep = '__')
  df <- as.data.frame(df)
  df <- df[, c("swMean", "lwMean","lumMean")]
  .rowNamesDF(df, make.names = T) <- rowNamez # makes unique names for dups
  colnames(df) <- c("s", "l", "lum") 
  return(df)
}
ferGm<- lapply(ferGm, formatPavoColdist)
ferGm<- ferGm[sapply(ferGm, nrow) >= 2] # much more elegent to get rid of low comparisons than I used in nearest dists

library(pavo)
ferDists<- lapply(ferGm, coldist, qcatch = 'Qi', n = c(1, 14),
                  weber = .05, achromatic = T) 
ferDists<- bind_rows(ferDists, .id =  'abbrevs_mspec') %>% 
  separate(patch1, into = c(NA, NA, 'comparison'), sep = '__', remove = T) %>%
  separate(patch2, into = c(NA, NA, 'comparison2'), sep = '__', remove = T) %>%
  separate(abbrevs_mspec, into = c('abbrevs', 'mspec'), sep = '__') %>%
  unite(col = 'comparison', comparison, comparison2, remove = T, sep = '') 
rm(ferGm)

ferDists$visSys<- 'ferret'
#### nearest dists ####
source('functions\\nearestDistances.R')
distsPf<- nearestDistances(pfow, visSys = 'pfowl')
distsPf<- distsPf %>% pivot_longer(cols = c('dS', 'dL'), names_to = 'visInfo',
                       values_to = 'distance') %>%
  group_by(abbrevs, mspec, visInfo) %>%
  mutate(qtile = ntile(distance, 10), visSys = 'pfowl') %>%
  slice_min(qtile)

distsHu<- nearestDistances(hum, visSys = 'human')
distsHu<- distsHu %>% pivot_longer(cols = c('dS', 'dL'), names_to = 'visInfo',
                       values_to = 'distance') %>%
  group_by(abbrevs, mspec, visInfo) %>%
  mutate(qtile = ntile(distance, 10), visSys = 'human') %>%
  slice_min(qtile)

distsFe<- nearestDistances(fer, visSys = 'ferret')
distsFe<- distsFe %>% pivot_longer(cols = c('dS', 'dL'), names_to = 'visInfo',
                       values_to = 'distance') %>%
  group_by(abbrevs, mspec, visInfo) %>%
  mutate(qtile = ntile(distance, 10), visSys = 'ferret') %>%
  slice_min(qtile)

ndists<- rbind(distsPf, distsFe, distsHu)
ndists<- merge(ndists, unique(fer[,c('abbrevs', 'collapseSppShrt')]), all.x = T, 
      by = 'abbrevs') %>%
  mutate(dataset = 'nrst') %>%
  select(!qtile)
##### Stats 
dists<- rbind(pfowlDists, humDists, ferDists)
dists<- merge(dists, unique(fer[,c('abbrevs', 'collapseSppShrt')]), all.x = T, 
          by = 'abbrevs') %>% 
  filter(comparison != 'rs') %>%
  mutate(dataset = 'all') %>% 
  pivot_longer(cols = c('dS', 'dL'), names_to = 'visInfo', 
               values_to = 'distance')

nrstAll<- rbind(dists, ndists)

nrstAll %>% filter(visInfo == 'dL', visSys == 'human') %>%
ggplot()+
  geom_boxplot(aes(fill = comparison, x = dataset,  y = distance))+
  facet_grid(vars(visSys))

popDists<- dists %>% group_by(abbrevs, comparison, visSys) %>%
  summarise(across(c('dS', 'dL'), mean, names = 'mean_{.col}'))

png("output//scatterp-ferr-pfowl-vs-hum.png",
    type = 'cairo', units = "mm", res = 300,
    width = 215, height = 260)
popDists %>% pivot_wider(names_from = visSys, values_from = c('dS', 'dL')) %>% 
  pivot_longer(cols = c('dS_pfowl', 'dS_ferret', 'dL_pfowl', 'dL_ferret'), 
               names_to = 'visSys') %>% 
  separate(visSys, into = c('visInfo', 'otherVisSys')) %>%
  pivot_wider(names_from = 'visInfo', values_from = 'value',
              names_glue = '{visInfo}_other') %>%
  ggplot(aes(dS_human, dS_other))+
  geom_point(aes(colour = otherVisSys))+
  geom_abline(lty=2)+
  facet_wrap(vars(comparison))+
  theme_grey()+
  coord_equal()+
  labs(colour = 'compared visual system')+
  theme(legend.position = 'top')+
  ylab('ferret / peafowl chroma distance')+
  xlab('human chroma distance')
dev.off()

# Shift in rock vs soil match between visual systems
# should be summarising by mspec FIRST and then by abbrevs

# does the ranking change?
order<- popDists %>% pivot_wider(names_from = c('visSys', 'comparison'),
                              values_from = c('dS', 'dL'))
popDists$abbrevs<- factor(popDists$abbrevs, 
                       levels = levels(fct_reorder(order$abbrevs, 
                                                   order$dS_ferret_lr,.desc = T)))
p1<- popDists %>% filter(comparison != 'ls') %>%
  group_by(abbrevs, comparison, visSys) %>% 
  summarise(dS = mean(dS, na.rm = T)) %>%
  ggplot(aes(abbrevs, dS))+
  geom_point(aes(colour = visSys))+
  geom_hline(yintercept = 1, lty = 2)+
  coord_flip()+
  theme_bw()+
  scale_colour_viridis_d()+
  ylab('lithops-rock chromatic distance')+
  xlab('populations')+
  labs(colour = 'compared visual system')+
  theme(legend.position = 'top')
  
popDists$abbrevs<- factor(popDists$abbrevs, 
                       levels = levels(fct_reorder(order$abbrevs, 
                                                   order$dS_ferret_ls,.desc = T)))
p2<- popDists %>% filter(comparison != 'lr') %>%
  group_by(abbrevs, comparison, visSys) %>% 
  summarise(dS = mean(dS, na.rm = T)) %>%
  ggplot(aes(abbrevs, dS))+
  geom_point(aes(colour = visSys))+
  geom_hline(yintercept = 1, lty = 2)+
  coord_flip()+
  theme_bw()+
  scale_colour_viridis_d()+
  ylab('lithops-soil chromatic distance')+
  xlab('populations')+
  labs(colour = 'compared visual system')+
  theme(legend.position = 'top')

popDists$abbrevs<- factor(popDists$abbrevs, 
                       levels = levels(fct_reorder(order$abbrevs, 
                                                   order$dL_ferret_lr,.desc = T)))
p3<- popDists %>% filter(comparison != 'ls') %>%
  group_by(abbrevs, comparison, visSys) %>% 
  summarise(dL = mean(dL, na.rm = T)) %>%
  ggplot(aes(abbrevs, dL))+
  geom_point(aes(colour = visSys))+
  geom_hline(yintercept = 1, lty = 2)+
  coord_flip()+
  theme_bw()+
  scale_colour_viridis_d()+
  ylab('lithops-rock luminance distance')+
  xlab('populations')+
  labs(colour = 'compared visual system')+
  theme(legend.position = 'none')

popDists$abbrevs<- factor(popDists$abbrevs, 
                       levels = levels(fct_reorder(order$abbrevs, 
                                                   order$dL_ferret_ls,.desc = T)))
p4<- popDists %>% filter(comparison != 'lr') %>%
  group_by(abbrevs, comparison, visSys) %>% 
  summarise(dL = mean(dL, na.rm = T)) %>%
  ggplot(aes(abbrevs, dL))+
  geom_point(aes(colour = visSys))+
  geom_hline(yintercept = 1, lty = 2)+
  coord_flip()+
  theme_bw()+
  scale_colour_viridis_d()+
  ylab('lithops-soil luminance distance')+
  xlab('populations')+
  labs(colour = 'compared visual system')+
  theme(legend.position = 'none')

png("output//dotplot-lr&ls-chroma-distances-per-pop-&-vissystem.png",
    type = 'cairo', units = "mm", res = 300, width = 215, height = 215)
gridExtra::grid.arrange(p1, p2, ncol =2)
dev.off()

png("output//dotplot-lr&ls-lum-distances-per-pop-&-vissystem.png",
    type = 'cairo', units = "mm", res = 300, width = 215, height = 215)
gridExtra::grid.arrange(p3, p4, ncol =2)
dev.off()
rm(p1,p2,p3,p4)

# facet wrap by info type
# statistics 
# proportion of significant differences in each direction across vis systems
###** Paired T-test, anno df, boxplots
t_test <- function(df, mu = 0, alt = "two.sided", paired = T,
                   conf.level = .95,var.equal = F){
  tidy(t.test(df$lr, df$ls,
              mu = mu, alt = alt,
              conf.level = conf.level,
              paired = paired, var.equal = var.equal))
}


pairedTests <- dists %>%
  pivot_longer(cols = c('dS', 'dL'),
               values_to = 'distance', names_to = 'visInfo') %>%
  pivot_wider(names_from = 'comparison', values_from = 'distance') %>%
  group_by(visSys,collapseSppShrt, abbrevs, visInfo) %>%
  nest() %>% 
  mutate(ttest = map(data, t_test)) %>%
  unnest(ttest) %>%
  group_by(visSys, abbrevs, visInfo) %>%
  mutate(p.adj = p.adjust(p.value, n = 56, method = 'holm')) %>%
  mutate(sigSym = case_when(p.adj < 0.05 ~ "*"),
         distance =
           case_when(visInfo == 'dS'~ 40,
                     visInfo == 'dL'~ 20),
         sigDir = case_when(estimate > 0 ~ 'Soil< Rock',
                                           TRUE ~ 'Rock< Soil'),
         facetVar = case_when(visInfo == 'dS' ~ 'chroma',
                              T ~ 'luminance'),
         facetVar = paste(facetVar, visSys)) 

boxpLRLSppop<- dists %>% 
  pivot_longer(cols = c('dS', 'dL'),
                       values_to = 'distance', names_to = 'visInfo') %>%
  ggplot()+ 
  geom_boxplot(aes(x = abbrevs, y = distance, fill = comparison))+
  geom_text(data = pairedTests, key_glyph = 'point',
            aes(label = sigSym, colour = sigDir,
                x = abbrevs, y = distance))+
  coord_flip()+
  facet_grid(collapseSppShrt~facetVar,
             space = 'free_y', scales = 'free', switch = 'y')+
  scale_fill_manual(values = c("lr" = "#00BFC4", "ls" = "#F8766D"),
                    labels = c("lithops-rock", "lithops-soil"))+
  scale_colour_manual(values = c("Rock< Soil" = "#00BFC4",
                                 "Soil< Rock" = "#F8766D"))+
  guides(colour = guide_legend(
    override.aes = list(size = 5, shape = c(utf8ToInt("*"), utf8ToInt("*")))))+
  labs(fill = "Substrate contrast",
       colour = 'Significance direction')+
  theme(legend.position = 'top',
        strip.placement = 'outside',
        strip.text.y = element_text(face = 'italic'),
        strip.background = element_rect(fill = 'grey91'),
        panel.spacing.y = unit(0,'lines'))
png("output//boxpLRLSppop.png",
    type = 'cairo', units = "mm", res = 300, width = 215, height = 280)
boxpLRLSppop
dev.off()
rm(boxpLRLSppop, pairedTests)

lrvslsSig<- pairedTests %>% 
  filter(p.adj < 0.05) %>%
  group_by(visSys, visInfo) %>%
  summarise(n = n(),
            'lr' = sum(statistic < 0),
            'ls' = sum(statistic > 0),
            'prop_lr' = lr / n*100,
            'prop_ls' =  ls / n*100) %>%
  select(-(c(n, prop_lr, prop_ls))) %>%
  pivot_longer(cols = c('lr', 'ls'), names_to = 'comparison', 
               values_to = 'significance') %>% # all this faffing to get columns ordered
  pivot_wider(names_from = c(visInfo, visSys, comparison), 
              values_from = c(significance),
              names_sort = T)

# summary stats
lrvsls<- popDists %>% group_by(abbrevs, comparison, visSys) %>%
  summarise(across(c('dS', 'dL'), mean, na.rm = T)) %>%
  pivot_longer(names_to = 'visInfo', cols = c('dL', 'dS')) %>%
  group_by(visSys, abbrevs, visInfo) %>%
  slice(which.min(value)) %>%
  group_by(visSys, visInfo) %>%
  count(comparison) %>%
  pivot_wider(names_from = c(visInfo, visSys, comparison), values_from = n,
              names_sort = T)
write.table(rbind(lrvsls, lrvslsSig),
            file = "output//lrvslsTab.txt", sep = ",", quote = FALSE, row.names = F)
rm(lrvsls, lrvslsSig)

# populations below threshold in each vis sys?
belowThresh<- popDists %>%
  pivot_longer(names_to = 'visInfo', cols = c('dS', 'dL'),
               values_to = 'distance') %>%
  group_by(comparison, visInfo, visSys) %>%
  summarise(`n<1` = sum(distance < 1)) %>%
  pivot_wider(names_from = c(visInfo, visSys, comparison), values_from = `n<1`,
              names_sort = T)
write.table(belowThresh, file = "output//belowThreshTab.txt", sep = ",", quote = FALSE, row.names = F)
rm(belowThresh)

# mean and standard error of LR and LS dists across visual sysems
meanDistTable <- popDists %>% 
  group_by(abbrevs, comparison, visSys) %>%
  summarise(across(c('dS', 'dL'), mean, na.rm = T)) %>%
  pivot_longer(names_to = 'visInfo', cols = c('dL', 'dS')) %>%
  group_by(visSys, visInfo, comparison) %>%
  summarise(distance = mean(value), sd = sd(value)) %>%
  mutate(sd = round(sd, digits = 2),
         distance = round(distance, digits = 2),
         sd = paste0('(', sd, ')'),) %>%
  unite('distance', c(distance, sd), sep = ' ') %>%
  pivot_wider(names_from = c(visInfo, visSys, comparison), 
              values_from = c(distance), names_sort = T)
write.table(meanDistTable, file = "output//distTab.txt", sep = ",", quote = FALSE, row.names = F)
rm(meanDistTable)

# overall tests
overallTtest <- popDists %>% # should be tukey post hoc?
  pivot_longer(cols = c('dS', 'dL'),
               values_to = 'distance', names_to = 'visInfo') %>%
  pivot_wider(names_from = 'comparison', values_from = 'distance') %>%
  group_by(visInfo, visSys) %>%
  nest() %>%
  mutate(ttest = map(data, t_test)) %>%
  unnest(ttest) %>% ungroup() %>%
  mutate(labelz = 
           case_when(p.value < 0.001 ~ '***',
                     p.value < 0.01 ~ '**',
                     p.value < 0.05 ~ '*'))

lrlsBoxpOverall<- popDists %>%
  group_by(abbrevs, comparison, visSys) %>% 
  summarise(across(c('dS', 'dL'),
                   mean, na.rm = T)) %>%
  pivot_longer(names_to = 'visInfo', cols = c('dS', 'dL'),
               values_to = 'distance') %>%
  ggplot(aes(visSys, distance))+
  geom_boxplot(aes(fill = comparison))+
  geom_text(data = overallTtest,
            aes(y = Inf,label = labelz), vjust = 1,size = 5)+
  geom_hline(yintercept = 1, lty = 2)+
  scale_fill_viridis_d()+
  facet_wrap(vars(visInfo), 
             labeller = labeller(visInfo = c(dS = 'chroma', dL = 'luminance')))
png("output//dotplot-lr&ls-chroma-distances-per-pop-&-vissystem.png",
    type = 'cairo', units = "mm", res = 300, width = 215, height = 150)
lrlsBoxpOverall
dev.off()
rm(lrlsBoxpOverall, overallTtest)

# models chroma
mod<- lm(dists$dS~dists$visSys+dists$comparison+dists$visSys*dists$comparison)
summary(mod) # p values = different from 0

modred<- lm(dists$dS~dists$visSys)
summary(modred)

anova(modred,mod) # testing whether comparison is sig

modred2<- lm(dists$dS~dists$visSys+dists$comparison)

anova(modred2, mod) # testing whether interaction is significant

# luminance
mod<- lm(dists$dL~dists$visSys+dists$comparison+dists$visSys*dists$comparison)
summary(mod) # p values = different from 0


# exploring interaction
popDists %>%
  pivot_wider(names_from = c('visSys'),
              values_from = c('dS', 'dL')) %>%
  mutate(aviHumRatio = dS_pfowl / dS_human) %>%
  group_by(abbrevs, comparison) %>%
  summarise(aviHumRatio = mean(aviHumRatio, na.rm = T)) %>%
  ggplot(aes(abbrevs, aviHumRatio))+
  geom_point(aes(colour = comparison))+
  coord_flip()

popDists %>% 
  pivot_wider(names_from = c('visSys'),
              values_from = c('dS', 'dL')) %>%
  mutate(aviHumRatio = dS_pfowl / dS_human) %>%
  group_by(abbrevs, comparison) %>%
  summarise(aviHumRatio = mean(aviHumRatio, na.rm = T)) %>%
  ggplot(aes(comparison, aviHumRatio))+
  geom_boxplot(aes(fill = comparison))+
  theme(text = element_text(family = 'serif'))




