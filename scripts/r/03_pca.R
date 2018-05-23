library(plyr)
library(gridExtra)
library(openxlsx)
library(tidyverse)
library(randomForest)
library(magrittr)

setwd('~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/')

chrom = 'chr20'
chrom = 'onep'
chrom = 'allchr'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get files from bucket
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gsutil = '/Users/robert/google-cloud-sdk/bin/gsutil'
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/pca/',chrom,'/pca_values.tsv ~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/out/', chrom, '/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/pca/',chrom,'/pca_scores.tsv ~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/out/', chrom, '/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/pheno/final_pheno_2018_01_04.tsv ~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/out/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/',chrom,'/sample_qc_info.txt ~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/out/', chrom, '/'))

# per population
# for(p in c('afr', 'amr', 'asj', 'fin', 'nfe', 'oth')) {
#   system(paste0(gsutil, paste0(' cp gs://wgspd-wgs-v1_temp/qc_measures/allchr/pca_scores_',p,'.tsv ~/Dropbox/postdoc/projects/WGSPD/out/')))
# }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ibd_exc = read.table(paste0('out/ibd_greater_0884_',chrom,'.txt'), stringsAsFactors = F)[,1]
relpairs = read.table(paste0('out/ibd_greater_0884_pairs_',chrom,'.txt'), stringsAsFactors = F, h=T)

values = read.table(paste0('out/', chrom, '/pca_values.tsv'), stringsAsFactors = F)
scores = read.table(paste0('out/', chrom, '/pca_scores.tsv'), stringsAsFactors = F, h=T)
names(scores) = tolower(gsub('pcaScores\\.', '', names(scores)))
#pops = read.table(paste0('out/', chrom, '/sample_qc_info.txt'), stringsAsFactors = F, h=T, comment.char = '', sep='\t')
#pheno = read.table('out/final_pheno_2018_01_04.tsv', stringsAsFactors = F, h=T, comment.char = '', sep='\t')
popvec = c('CHB'='EAS', 'JPT'='EAS', 'CHS'='EAS', 'CDX'='EAS', 'KHV'='EAS',
           'CEU'='EUR', 'TSI'='EUR', 'FIN'='EUR', 'GBR'='EUR', 'IBS'='EUR',
           'YRI'='AFR', 'LWK'='AFR', 'GWD'='AFR', 'MSL'='AFR', 'ESN'='AFR', 'ASW'='AFR', 'ACB'='AFR',
           'MXL'='AMR', 'PUR'='AMR', 'CLM'='AMR', 'PEL'='AMR',
           'GIH'='SAS', 'PJL'='SAS', 'BEB'='SAS', 'STU'='SAS', 'ITU'='SAS')
onekgpops = read.xlsx('~/Downloads/20130606_sample_info.xlsx')[,c(1,3)]
onekgpops[,1] = paste0('1KG_', onekgpops[,1])
onekgpops %<>% mutate(super=popvec[Population])
write.table(onekgpops, '~/Downloads/onekgpops.tsv', sep='\t', col.names = T, row.names=F, quote=F)
system(paste0(gsutil, paste0(' cp ~/Downloads/onekgpops.tsv gs://ccdg-qc-multi/data/1000genomes/')))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot pcs against populations in 1KG
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m = right_join(onekgpops, scores, by=c('Sample'='s'))
m = data.frame(onekg = grepl('1KG_', m$Sample), imp=FALSE, m)

source('~/Dropbox/postdoc/projects/ccdg_qc/r/pop_colors.R')
cols = pop_color[,2]; names(cols) = pop_color[,1]

dat = m; dat$Population = factor(dat$Population, levels=unique(dat$Population[order(dat$super)]))
ggplot(dat, aes(pc1, pc2, col=Population)) + geom_point() + scale_colour_manual(values=cols)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# assign PC scores to relatives
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# impute_missing = function(scores, relpairs) {
#   # returns scores with filled in PCs
#   # scores: s, pc1... pc10
#   # relpairs: ID1, ID2
#   left_join(relpairs, scores, by=c('ID2'='s')) %>% ddply('ID1', function(x) apply(x[,-1:-2], 2, function(y) mean(y, na.rm=T))) %>% dplyr::rename(s=ID1)
# }

impute_missing = function(data, relpairs, missing) {
  # returns list of lenght 2: [[1]]: data with filled in PCs; [[2]]: vector of matched IDs
  # data: s, pc1... pc10
  # relpairs: ID1, ID2, Kinship
  # missing: IDs which need to be filled
  
  # find matches
  out = data.frame()
  mat = c()
  for(id in missing) {
    pot = relpairs %>% filter((ID1 == id & !ID2 %in% missing) | (ID2 == id & !ID1 %in% missing))
    relord = 1
    while(nrow(pot) == 0) {
      print(c(id, relord))
      pot = relpairs %>% filter((ID1 == id & ID2 %in% missing) | (ID2 == id & ID1 %in% missing))
      pot = pot[order(pot$Kinship, decreasing=TRUE),]
      closrel = ifelse(id == pot$ID1[relord], pot$ID2[relord], pot$ID1[relord])
      pot = relpairs %>% filter(((ID1 == closrel & !ID2 %in% missing) | (ID2 == closrel & !ID1 %in% missing)))
      relord = relord + 1
    }
    pot = pot[order(pot$Kinship, decreasing=TRUE),]
    p = c(pot$ID2, pot$ID1); p = p[p!=id & ! p %in% missing & p %in% data$s]
    mat[id] = p[1]
    #print(mat[id])
    out = rbind(out, data.frame(s=id, data[data$s == mat[id], -1]))
  }
  out[,1] = as.character(out[,1])
  list(out, mat)
}

imputed_list = impute_missing(scores, relpairs, setdiff(relpairs$ID1, scores[,1]))
imputed = imputed_list[[1]]
matched_relatives = imputed_list[[2]]
m2 = bind_rows(m, data.frame(onekg=F, imp=T, Sample=imputed$s, Population=NA, super=NA, imputed[,-1], stringsAsFactors = F))
m2 %<>% mutate(mypops = ifelse(Population %in% c('FIN', 'PUR'), Population, super))

matched_relatives_df = data.frame(names(matched_relatives), matched_relatives, stringsAsFactors = F)
names(matched_relatives_df) = c('excluded_relative', 'matched_relative')
write.table(matched_relatives_df, '~/Downloads/matched_relatives.tsv', col.names=T, row.names=F, quote=F, sep='\t')
system(paste0(gsutil, paste0(' cp ~/Downloads/matched_relatives.tsv gs://ccdg-qc-multi/qc_measures/pca/',chrom,'/matched_relatives.tsv')))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# predict populations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pop_forest = function(training_data, data, ntree=100, seed=42, pcs=1:6) {
  set.seed(seed)
  require(randomForest)
  form = formula(paste('as.factor(known_pop) ~', paste0('pc', pcs, collapse = ' + ')))
  forest = randomForest(form,
                        data = training_data,
                        importance = T,
                        ntree = ntree)
  print(forest)
  fit_data = data.frame(predict(forest, data, type='prob'), combined_sample = data$combined_sample, stringsAsFactors = F)
  fit_data %>%
    gather(predicted_pop, probability, -combined_sample) %>%
    group_by(combined_sample) %>%
    slice(which.max(probability))
}

trdat = m2 %>% filter(onekg) %>% select(known_pop=mypops, pc1:pc10)
tedat = m2 %>% filter(!onekg) %>% select(combined_sample=Sample, pc1:pc10)
predicted = as.data.frame(pop_forest(training_data = trdat, data = tedat))

m3 = m2 %>% left_join(predicted, by=c('Sample'='combined_sample')) %>% mutate(Population = ifelse(is.na(Population), ifelse(probability < .9, 'other', predicted_pop), Population), super=popvec[Population])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot pcs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d3 = m3[!m3$onekg,]
# nrow(d3) should be number of samples to keep (22524)
cols = c(gg_color_hue(length(unique(d3$Population))-1), '#D3D3D3')

d3$Population = factor(d3$Population, levels=c(sort(setdiff(unique(d3$Population), 'other')), 'other'))
th = theme(panel.background = element_blank(), legend.position = 'none', axis.text=element_text(size=12), axis.title = element_text(size=30))
co = scale_colour_manual(values=cols)
p12 = ggplot(d3, aes(pc1, pc2, col=Population)) + geom_point(size=.1) + th + co
p23 = ggplot(d3, aes(pc3, pc2, col=Population)) + geom_point(size=.1) + th + co
p34 = ggplot(d3, aes(pc3, pc4, col=Population)) + geom_point(size=.1) + th + co
p45 = ggplot(d3, aes(pc5, pc4, col=Population)) + geom_point(size=.1) + th + co
p56 = ggplot(d3, aes(pc5, pc6, col=Population)) + geom_point(size=.1) + th + co
p67 = ggplot(d3, aes(pc7, pc6, col=Population)) + geom_point(size=.1) + th + co
p0 = g_legend(ggplot(d3, aes(pc1, pc2, col=Population)) + geom_point() + theme(panel.background = element_blank(), legend.title = element_blank(), legend.key.size = unit(20, 'mm'), legend.text = element_text(size=24)) + co + guides(colour = guide_legend(override.aes = list(size=10))))

png(paste0('plots/',chrom,'/03_pca1to7_population.png'), width=1600, height=1000)
grid.arrange(p12, p23, p34, p45, p56, p0, nrow=2)
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# upload predicted populations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nam = paste0('~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/pca/',chrom,'/predicted_populations.tsv')
m4 = m3 %>% filter(!onekg) %>% select(-onekg, -mypops)
write.table(m4, nam, quote=F, col.names=T, row.names=F, sep='\t')
system(paste0(gsutil, ' cp ',nam,' gs://ccdg-qc-multi/qc_measures/pca/',chrom,'/'))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# combine with self-reported data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m4 = left_join(m3, pheno, by=c('Sample'='id'))
as.data.frame(table(m4$super, m4$race)) %>% dplyr::rename(predicted=Var1, reported=Var2) %>%
  ggplot(aes(predicted, Freq, fill=predicted)) + geom_bar(stat='identity') + facet_grid(~ reported) + theme(panel.background = element_blank())

#racestat = gather(as.data.frame(table(m4$super, m4$race)), self_reported, genetic, )


a = m4 %>% select(Sample, Population, super, race) %>% filter(!is.na(race)) %>% arrange(factor(super, levels=levels(as.factor(super))[c(1,3,2,5,4)]), Population, Sample) %>% mutate(o1=1:n()) %>% arrange(race, Population, Sample) %>% mutate(o2=1:n()) %>% select(Population, super, race, Sample, o1, o2) 

td1 = a %>% group_by(super) %>% summarise(t='o1', mid = mean(o1)) %>% mutate(Sample=NA)
td2 = a %>% group_by(Population) %>% summarise(t='o1', mid = mean(o1)) %>% mutate(Sample=NA, super=NA)
td3 = a %>% group_by(race) %>% summarise(t='o2', mid = mean(o2)) %>% mutate(Sample=NA, super=NA)

a %<>% gather(key='type', value='y', o1:o2)
ggplot(a, aes(type, y, group=Sample, col=super)) + geom_line(alpha=.1) + theme(panel.background = element_blank(), legend.position = 'none', text=element_text(size=16)) + geom_point(aes(type, y, col=super)) + geom_text(data=td1, aes(t, mid, label=super), nudge_x=-.2, size=5) + geom_text(data=td2, aes(t, mid, label=Population), nudge_x=-.1, size=5) + geom_text(data=td3, aes(t, mid, label=race), nudge_x=.04, size=5, hjust=0) + xlab('') + ylab('') + scale_x_discrete(labels=c('predicted', 'self-reported'))
ggsave(paste0('plots/',chrom,'/03_pca_selfreported_predicted.pdf'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot PCs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ggplot(m3[m3$onekg,], aes(pc1, pc2)) + theme(panel.background = element_blank(), legend.position = 'none') + scale_fill_manual(values=cols) + geom_point(data=m3[!m3$onekg,], aes(pc1, pc2, col=Population), size=1) + geom_density2d(aes(group=Population), col='black', bins=2) + facet_wrap( ~ Population, scales='free') + scale_colour_manual(values=cols)
ggsave(paste0('plots/',chrom,'/03_pca_1KG_density_by_population.pdf'))

ggplot(m3[m3$onekg,], aes(pc1, pc2)) + theme(panel.background = element_blank(), legend.position = 'none') + scale_fill_manual(values=cols) + geom_point(data=m3[!m3$onekg,], aes(pc1, pc2, col=probability > .6), size=1) + geom_density2d(aes(group=Population), col='black', bins=2) + facet_wrap( ~ Population, scales='free') #+ scale_colour_manual(values=cols)
ggsave(paste0('plots/',chrom,'/03_pca_1KG_density_by_population_prob_greater_03.pdf'))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# predict super-population only
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

trdat = m2 %>% filter(onekg) %>% select(known_pop=super, pc1:pc10)
tedat = m2 %>% filter(!onekg) %>% select(combined_sample=Sample, pc1:pc10)
predicted = as.data.frame(pop_forest(training_data = trdat, data = tedat))

m3 = m2 %>% mutate(super = ifelse(is.na(super), predicted$predicted_pop[match(Sample, predicted$combined_sample)], super)) %>% mutate(probability=predicted$probability[match(Sample, predicted$combined_sample)])



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot PCs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ggplot(m3[m3$onekg,], aes(pc1, pc2)) + theme(panel.background = element_blank(), legend.position = 'none') + scale_fill_manual(values=cols) + geom_point(data=m3[!m3$onekg,], aes(pc1, pc2, col=super), size=1) + geom_density2d(aes(group=super), col='black', bins=2) + facet_wrap( ~ super, scales='free') #+ scale_colour_manual(values=cols)

ggplot(m3[m3$onekg,], aes(pc1, pc2)) + theme(panel.background = element_blank(), legend.position = 'none') + scale_fill_manual(values=cols) + geom_point(data=m3[!m3$onekg,], aes(pc1, pc2, col=probability > .99), size=1) + geom_density2d(aes(group=super), col='black', bins=2) + facet_wrap( ~ super, scales='free') #+ scale_colour_manual(values=cols)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define strict super populations based on 2 mads only
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mads = 2
maddat = m4 %>% filter(onekg) %>% group_by(super) %>% summarise(median_pc1 = median(pc1), mad_pc1=mad(pc1),
                                                            median_pc2 = median(pc2), mad_pc2=mad(pc2),
                                                            median_pc3 = median(pc3), mad_pc3=mad(pc3),
                                                            median_pc4 = median(pc4), mad_pc4=mad(pc4),
                                                            median_pc5 = median(pc5), mad_pc5=mad(pc5),
                                                            median_pc6 = median(pc6), mad_pc6=mad(pc6)) %>% as.data.frame()
m5 = m4 %>% left_join(maddat, by='super') %>% mutate(within2mad = pc1 < median_pc1 + mads*mad_pc1 & pc1 > median_pc1 - mads*mad_pc1 &
                                                             pc2 < median_pc2 + mads*mad_pc2 & pc2 > median_pc2 - mads*mad_pc2 &
                                                             pc3 < median_pc3 + mads*mad_pc3 & pc3 > median_pc3 - mads*mad_pc3 &
                                                             pc4 < median_pc4 + mads*mad_pc4 & pc4 > median_pc4 - mads*mad_pc4 &
                                                             pc5 < median_pc5 + mads*mad_pc5 & pc5 > median_pc5 - mads*mad_pc5 &
                                                             pc6 < median_pc6 + mads*mad_pc6 & pc6 > median_pc6 - mads*mad_pc6)

ggplot(m5[m5$onekg,], aes(pc1, pc2)) + theme(panel.background = element_blank(), legend.position = 'none') + scale_fill_manual(values=cols) + geom_point(data=m5[!m5$onekg,], aes(pc1, pc2, col=within2mad), size=1) + geom_density2d(aes(group=super), col='black', bins=2) + facet_wrap( ~ super, scales='free')
ggsave(paste0('plots/',chrom,'/03_pca_1KG_density_by_super_2mad.pdf'))

th = theme(panel.background = element_blank(), legend.position = 'none')
p12 = ggplot(m5[!m5$onekg,], aes(pc1, pc2, col=super)) + geom_point(size=.1) + th
p23 = ggplot(m5[!m5$onekg,], aes(pc3, pc2, col=super)) + geom_point(size=.1) + th
p34 = ggplot(m5[!m5$onekg,], aes(pc3, pc4, col=super)) + geom_point(size=.1) + th
p45 = ggplot(m5[!m5$onekg,], aes(pc5, pc4, col=super)) + geom_point(size=.1) + th
p56 = ggplot(m5[!m5$onekg,], aes(pc5, pc6, col=super)) + geom_point(size=.1) + th
p67 = ggplot(m5[!m5$onekg,], aes(pc7, pc6, col=super)) + geom_point(size=.1) + th
p0 = g_legend(ggplot(m5[!m5$onekg,], aes(pc1, pc2, col=super)) + geom_point() + theme(panel.background = element_blank(), legend.title = element_blank()))

pdf(paste0('plots/',chrom,'/03_pca1to7_super.pdf'), width=14, height=8)
grid.arrange(p12, p23, p34, p45, p56, p0, nrow=2)
dev.off()

th = theme(panel.background = element_blank(), legend.position = 'none')
al = scale_alpha_manual(values=c(0,1))
gp = geom_point(size=.1, aes(alpha=within2mad))
p12 = ggplot(m5[!m5$onekg,], aes(pc1, pc2, col=super)) + gp + al + th
p23 = ggplot(m5[!m5$onekg,], aes(pc3, pc2, col=super)) + gp + al + th
p34 = ggplot(m5[!m5$onekg,], aes(pc3, pc4, col=super)) + gp + al + th
p45 = ggplot(m5[!m5$onekg,], aes(pc5, pc4, col=super)) + gp + al + th
p56 = ggplot(m5[!m5$onekg,], aes(pc5, pc6, col=super)) + gp + al + th
p67 = ggplot(m5[!m5$onekg,], aes(pc7, pc6, col=super)) + gp + al + th
p0 = g_legend(ggplot(m5[!m5$onekg & m5$within2mad,], aes(pc1, pc2, col=super)) + geom_point() + theme(panel.background = element_blank(), legend.title = element_blank()))

pdf(paste0('plots/',chrom,'/03_pca1to7_super_within2mad.pdf'), width=14, height=8)
grid.arrange(p12, p23, p34, p45, p56, p0, nrow=2)
dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add population specific PCs to predicted_pop file
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/pca/',chrom,'/predicted_populations.tsv ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/pca/',chrom,'/matched_relatives.tsv ~/Downloads/'))

predicted_pop = read.table('~/Downloads/predicted_populations.tsv', stringsAsFactors = F, h=T)
matched_relatives = read.table('~/Downloads/matched_relatives.tsv', stringsAsFactors = F, h=T)
# check if predicted pop contains all samples in samples_to_keep (22524)
  
popspecpcs = list()
for(p in setdiff(unique(predicted_pop$Population), 'other')) {
  #gsutil ls -lh gs://ccdg-qc-multi/qc_measures/pca/onep/perpop/pca_scores_*tsv
  system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/pca/',chrom,'/perpop/pca_scores_',p,'.tsv ~/Downloads/'))
  dat = read.table(paste0('~/Downloads/pca_scores_',p,'.tsv'), h=T, stringsAsFactors = F)
  names(dat)[-1] = paste0('pop_', names(dat)[-1])
  popspecpcs[[p]] = dat
}
dat = do.call(rbind, popspecpcs)

pp2 = left_join(predicted_pop, dat, by=c('Sample'='s'))
misids = pp2$Sample[is.na(pp2$pop_PC1)]
fillids = matched_relatives$matched_relative[match(misids, matched_relatives$excluded_relative)]
pp2[misids, paste0('pop_PC', 1:10)] = pp2[match(pp2$Sample, ), paste0('pop_PC', 1:10)]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add predicted populations to pheno file
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

predicted_pop = read.table('~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/pca/allchr/predicted_populations.tsv', h=T, stringsAsFactors = F)
pheno = read.table('~/Dropbox/postdoc/projects/ccdg_qc/data/phenotypes/combined_phenotypes_feb2018_plus_taichi.txt', h=T, stringsAsFactors = F, sep='\t', quote='')
pheno %<>% left_join(select(predicted_pop, Sample, Population, probability), by=c('subject_id'='Sample')) %>% filter(study_nickname != 'ccdg_broad_taichi_ds-irb-col-mds')
write.table(pheno, '~/Dropbox/postdoc/projects/ccdg_qc/data/phenotypes/combined_phenotypes_feb2018_notaichi_pop.txt', quote=F, col.names=T, row.names=F, sep='\t')

# a = pheno %>% select(subject_id, Population, race_ethnicity) %>% filter(!is.na(Population)) %>% arrange(factor(Population, levels=levels(as.factor(Population))[c(1,3,2,5,4,6,7)]), subject_id) %>% mutate(o1=1:n()) %>% arrange(race_ethnicity, Population, subject_id) %>% mutate(o2=1:n()) %>% select(Population, race_ethnicity, subject_id, o1, o2)
# 
# td2 = a %>% group_by(Population) %>% summarise(t='o1', mid = mean(o1)) %>% mutate(subject_id=NA)
# td3 = a %>% group_by(race_ethnicity) %>% summarise(t='o2', mid = mean(o2)) %>% mutate(subject_id=NA, Population=NA)
# 
# a %<>% gather(key='type', value='y', o1:o2)
# 
# ggplot(a, aes(type, y, group=subject_id, col=Population)) + geom_line(alpha=.1) + theme(panel.background = element_blank(), legend.position = 'none', text=element_text(size=16)) + geom_point(aes(type, y, col=Population)) + geom_text(data=td2, aes(t, mid, label=Population), nudge_x=-.1, size=5) + geom_text(data=td3, aes(t, mid, label=race_ethnicity), nudge_x=.04, size=5, hjust=0) + xlab('') + ylab('') + scale_x_discrete(labels=c('predicted', 'self-reported'))







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot predicted ancestry composition of cohorts
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pop_by_cohort = left_join(predicted_pop, manifest, by=c('Sample'='Sample.ID'))
pop_by_cohort %<>% group_by(Center, Filename_Prefix, Population) %>% summarize(count = n()) #%>% mutate(Filename_Prefix2 = as.factor(Filename_Prefix))

ggplot(as.data.frame(pop_by_cohort), aes(Filename_Prefix, count, fill=Population)) + geom_bar(stat='identity') + facet_grid(. ~ Center, scales='free_x', space='free_x', drop=T) + xlab('') + ylab('') + theme(panel.background = element_blank(), text=element_text(size=14), axis.text=element_text(angle=90, hjust=1))
ggsave(paste0('plots/',chrom,'/03_populations_by_cohort_and_center.pdf'), width=12, height=6)


# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # add PCs to phenotype file
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# pheno2 = left_join(pheno, bind_rows(dplyr::rename(scores, ID='s'), dplyr::rename(imputed, ID='ID1')))
# names(pheno2) = gsub('pcaScores\\.', '', names(pheno2))
# pheno2$random = rnorm(nrow(pheno2))
# pheno2$random = sample(0:1, nrow(pheno2), replace=T)
# pheno2$jitter = jitter(pheno2$scz_dich,50)
# 
# write.table(pheno2, 'out/final_pheno_2018_01_04_PCs.tsv', col.names=T, row.names=F, quote=F, sep='\t')
# system(paste0(gsutil, ' cp ~/Dropbox/postdoc/projects/WGSPD/out/final_pheno_2018_01_04_PCs.tsv gs://wgspd-wgs-v1_temp/pheno/'))
# 
# 
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #
# # Pop PCs
# #
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# # read scores
# 
# scorelist = list()
# for(p in c('afr', 'amr', 'asj', 'fin', 'nfe', 'oth')) {
#   scorelist[[p]] = read.table(paste0('out/pca_scores_',p,'.tsv'), stringsAsFactors = F, h=T)
# }
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # assign PC scores to relatives
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# implist = list()
# for(i in 1:length(scorelist)) {
#   sc = scorelist[[i]]
#   implist[[i]] = left_join(select(relpairs, 'ID1', s='ID2'), sc) %>% ddply('ID1', function(x) apply(x[,-1:-2], 2, function(y) mean(y, na.rm=T)))
# }
# 
# aa = abind(lapply(implist, `[`, -1), along=3)
# imputed2 = data.frame(implist[[1]][,1], apply(aa, 1:2, mean, na.rm=T), stringsAsFactors = F)
# names(imputed2) = paste0(names(imputed2), '_pop')
# names(imputed2)[1] = 'ID1'
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # add PCs to phenotype file
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# pheno2 = left_join(pheno, bind_rows(dplyr::rename(scores, ID='s'), dplyr::rename(imputed, ID='ID1')))
# x1 = dplyr::rename(imputed2, ID='ID1')
# x2 = dplyr::rename(do.call(rbind, scorelist), ID='s'); names(x2) = names(x1)
# pheno2 = left_join(pheno2, bind_rows(x2, x1))
# names(pheno2) = gsub('pcaScores\\.', '', names(pheno2))
# 
# write.table(pheno2, 'out/final_pheno_2018_01_04_PCs.tsv', col.names=T, row.names=F, quote=F, sep='\t')
# system(paste0(gsutil, ' cp ~/Dropbox/postdoc/projects/WGSPD/out/final_pheno_2018_01_04_PCs.tsv gs://wgspd-wgs-v1_temp/pheno/'))

