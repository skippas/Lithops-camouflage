sppLabeller<- c(olivacea = 'L. oliv.', localis = 'L. loc.',
                marmorataSpp = 'L. marmorata', divergensSpp = 'L. divergens',
                halliiSpp = 'L. hallii', leslieiSpp = 'L. lesliei',
                otzeniana = 'L. otz.', compt. = 'L. compt.',
                fulleri = 'L. fulleri', bromf. = 'L. bromf.',
                hookeriSpp = 'L. hookeri', meyeri = 'L. me.', 
                dint. = 'L. din.', verruc. = 'L. verru.', 
                singlePops = 'singlePops')

lut<- read.csv('data//lookup_table.csv')

names(sppLabeller)
levels(pfow$collapseSpp) %in% names(sppLabeller)
names(sppLabeller) %in% levels(pfow$collapseSpp) 

lut<- merge(lut, unique(pfow[, c('abbrevs', 'collapseSpp')]), by = 'abbrevs')
sppLabeller<- stack(sppLabeller)
names(sppLabeller) <- c('collapseSppShrt', 'collapseSpp')
lut<- merge(lut, sppLabeller)

write.csv(lut,
            file = "data//lut.csv", sep = ",", quote = FALSE, row.names = F)


