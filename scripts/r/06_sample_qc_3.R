#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# this script evaluates sample QC after variant QC
# one the same samples as the previous sample QC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library(plyr)
library(tidyverse)
library(magrittr)
library(viridis)

setwd('~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/')

chrom = 'chr20'
chrom = 'onep'
chrom = 'allchr'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define thresholds
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mincoverage = 15
mincallrate = 0.96

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# copy files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gsutil = '/Users/robert/google-cloud-sdk/bin/gsutil'
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/',chrom,'/sample_qc_info.txt ~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/',chrom,'/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/',chrom,'/sample_qc_info_postqc.txt ~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/',chrom,'/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/',chrom,'/01_sample_qc_keep.txt ~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/',chrom,'/'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sampleqc = read.table(paste0('~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/',chrom,'/sample_qc_info_postqc.txt'), h=T, stringsAsFactors = F, sep='\t')
sampleqc_pre = read.table(paste0('~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/',chrom,'/sample_qc_info.txt'), h=T, stringsAsFactors = F, sep='\t')
names(sampleqc_pre) = paste0(names(sampleqc_pre), '_pre')
samples_to_keep = read.table(paste0('~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/',chrom,'/01_sample_qc_keep.txt'), stringsAsFactors = F)[,1]
manifest = read.table('~/Dropbox/postdoc/projects/ccdg_qc/CCDG_Freeze_1_Subsetting_Manifest_2017-12-06.txt', sep='\t', h=T, stringsAsFactors = F, quote='')
manifest %<>% filter(!Remove)
sex_mismatch = read.table(paste0('~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/',chrom,'/sex_mismatch.txt'), stringsAsFactors = F)
pheno = read.table('~/Dropbox/postdoc/projects/ccdg_qc/data/phenotypes/combined_phenotypes_broad.txt', stringsAsFactors = F, h=T, comment.char = '', sep='\t')


# import ancestry
predicted_populations = read.table(paste0('~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/pca/',chrom,'/predicted_populations.tsv'), stringsAsFactors = F, h=T)

labs = as.data.frame(do.call(rbind, strsplit(manifest$Filename_Prefix, '_')))
dup = names(which(colSums(table(labs[,2], labs[,3]) > 0) > 1))
labs = paste0(ifelse(labs[,3] %in% dup, paste0(labs[,2], '_'), ''), labs[,3])

sampleqc2 = left_join(manifest, sampleqc, by=c('Sample.ID' = 's')) %>%
  mutate(label = labs) %>%
  filter(Sample.ID %in% samples_to_keep) %>%
  #left_join(predicted_populations %>% select(Sample.ID=Sample, Population), by='Sample.ID') %>%
  arrange(Center, Project, Population) %>%
  mutate(borders = (1:n()) %in% cumsum(rle(Project)$lengths)) %>%
  mutate(odd = as.numeric(factor(Project, levels=unique(Project)))%%2==0) %>%
  left_join(sampleqc_pre, by=c('Sample.ID'='s_pre'))


nums = rowSums(table(sampleqc2$Center, sampleqc2$Project) > 0)


popvec = c('CHB'='EAS', 'JPT'='EAS', 'CHS'='EAS', 'CDX'='EAS', 'KHV'='EAS',
           'CEU'='EUR', 'TSI'='EUR', 'FIN'='EUR', 'GBR'='EUR', 'IBS'='EUR',
           'YRI'='AFR', 'LWK'='AFR', 'GWD'='AFR', 'MSL'='AFR', 'ESN'='AFR', 'ASW'='AFR', 'ACB'='AFR',
           'MXL'='AMR', 'PUR'='AMR', 'CLM'='AMR', 'PEL'='AMR',
           'GIH'='SAS', 'PJL'='SAS', 'BEB'='SAS', 'STU'='SAS', 'ITU'='SAS')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot Ti/Tv, Het/Hom, Ins/Del
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

wdt = 20 
hgh = 5
custompops = c(sort(setdiff(unique(predicted_populations$Population), 'other')), 'other')
#source('~/Dropbox/postdoc/projects/ccdg_qc/r/pop_colors.R')
#cols = pop_color[,2]; names(cols) = pop_color[,1]
cols = c(gg_color_hue(length(custompops)-1), '#D3D3D3')

sampleqc3 = sampleqc2 %>%
  group_by_at(vars(Population, Project)) %>%
  mutate(index=row_number(),
         median_rTiTv=median(sample_qc.r_ti_tv_pre, na.rm=T),
         mad_rTiTv=mad(sample_qc.r_ti_tv_pre, na.rm=T, constant=1),
         median_rHetHomVar=median(sample_qc.r_het_hom_var_pre, na.rm=T),
         mad_rHetHomVar=mad(sample_qc.r_het_hom_var_pre, na.rm=T, constant=1),
         median_rInsertionDeletion=median(sample_qc.r_insertion_deletion_pre, na.rm=T),
         mad_rInsertionDeletion=mad(sample_qc.r_insertion_deletion_pre, na.rm=T, constant=1),
         Population2 = factor(Population, levels=custompops)) %>%
  arrange(Center, Project, Population)

borderdat = sampleqc3 %>%
  ungroup %>%
  slice(which(borders)) %>% 
  mutate(borderpos=which(sampleqc3$borders)) %>%
  select(Center, Project, Population, Population2, borderpos, label, odd, median_rTiTv, mad_rTiTv, median_rHetHomVar, mad_rHetHomVar, median_rInsertionDeletion, mad_rInsertionDeletion) %>%
  mutate(start=c(1,borderpos[-length(borderpos)])) %>%
  mutate(mid=(start+borderpos)/2)

# Ti/Tv
p = ggplot(sampleqc3, aes(x=1:nrow(sampleqc3), y=sample_qc.r_ti_tv, colour=Population2)) +
  geom_point(size=.5) +
  #geom_point(aes(size=Sample.ID %in% badids[,1]))  + scale_size_manual(values=c(.5, 3)) +
  ylab('Ti/Tv ratio') +
  xlab('') +
  facet_grid(~ Center, scale='free_x', space = 'free_x') +
  theme(legend.position = 'bottom', legend.title = element_blank(), panel.background = element_blank(), text=element_text(size=16)) +
  guides(col=guide_legend(nrow=1, override.aes = list(size=6))) +
  geom_point(aes(x=1:nrow(sampleqc3), y=median_rTiTv + 4*mad_rTiTv), col='black', shape='_') +
  geom_point(aes(x=1:nrow(sampleqc3), y=median_rTiTv - 4*mad_rTiTv), col='black', shape='_') +
  geom_vline(data=borderdat, aes(xintercept = borderpos), col='lightgrey') +
  scale_colour_manual(values=cols)
#geom_segment(data=borderdat, aes(x=start, xend=borderpos, y=median_rTiTv + 4*mad_rTiTv, yend=median_rTiTv + 4*mad_rTiTv, col=race)) +
#geom_segment(data=borderdat, aes(x=start, xend=borderpos, y=median_rTiTv - 4*mad_rTiTv, yend=median_rTiTv - 4*mad_rTiTv, col=race))
lb = geom_label(data=borderdat, aes(x=mid, y=2.10-(1:nrow(borderdat)%%4)*.01, label=label), col='black', angle=90, hjust=.5, alpha=.7)

ggsave(paste0('plots/',chrom,'/06_sample_qc_3_',chrom,'_titv.pdf'), p, width=wdt, height=hgh)
ggsave(paste0('plots/',chrom,'/06_sample_qc_3_',chrom,'_titv_labs.pdf'), p + lb, width=wdt, height=hgh)

# Ti/Tv arranged by population, coloured by center
sampleqc3_v2 = sampleqc3 %>% arrange(Population2, Center, Project)
p = ggplot(sampleqc3_v2, aes(x=1:nrow(sampleqc3_v2), y=sample_qc.r_ti_tv, colour=Center)) +
  geom_point(size=.5) +
  #geom_point(aes(size=Sample.ID %in% badids[,1]))  + scale_size_manual(values=c(.5, 3)) +
  ylab('Ti/Tv ratio') +
  xlab('') +
  facet_grid(~ Population2, scale='free_x', space = 'free_x') +
  theme(legend.position = 'bottom', legend.title = element_blank(), panel.background = element_blank(), text=element_text(size=16), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(col=guide_legend(nrow=1, override.aes = list(size=6))) +
  geom_point(aes(x=1:nrow(sampleqc3_v2), y=median_rTiTv + 4*mad_rTiTv, col=Center), shape='_') +
  geom_point(aes(x=1:nrow(sampleqc3_v2), y=median_rTiTv - 4*mad_rTiTv, col=Center), shape='_')

ggsave(paste0('plots/',chrom,'/06_sample_qc_3_',chrom,'_titv_bypop.pdf'), p, width=wdt, height=hgh)


wdt = 20 
hgh = 5
source('~/Dropbox/postdoc/projects/ccdg_qc/r/pop_colors.R')
cols = pop_color[,2]; names(cols) = pop_color[,1]

sampleqc3 = sampleqc2 %>%
  group_by_at(vars(Population, Project)) %>%
  mutate(index=row_number(),
         median_rTiTv=median(sa.qc.rTiTv_pre, na.rm=T),
         mad_rTiTv=mad(sa.qc.rTiTv_pre, na.rm=T, constant=1),
         median_rHetHomVar=median(sa.qc.rHetHomVar_pre, na.rm=T),
         mad_rHetHomVar=mad(sa.qc.rHetHomVar_pre, na.rm=T, constant=1),
         median_rInsertionDeletion=median(sa.qc.rInsertionDeletion_pre, na.rm=T),
         mad_rInsertionDeletion=mad(sa.qc.rInsertionDeletion_pre, na.rm=T, constant=1),
         Population2 = factor(Population, levels=names(popvec))) %>%
  arrange(Center, Project, Population)

borderdat = sampleqc3 %>%
  ungroup %>%
  slice(which(borders)) %>% 
  mutate(borderpos=which(sampleqc3$borders)) %>%
  select(Center, Project, Population, Population2, borderpos, label, odd, median_rTiTv, mad_rTiTv, median_rHetHomVar, mad_rHetHomVar, median_rInsertionDeletion, mad_rInsertionDeletion) %>%
  mutate(start=c(1,borderpos[-length(borderpos)])) %>%
  mutate(mid=(start+borderpos)/2)

# Ti/Tv
p = ggplot(sampleqc3, aes(x=1:nrow(sampleqc3), y=sa.qc.rTiTv, colour=Population2)) +
  geom_point(size=.5)  +
  ylab('Ti/Tv ratio') +
  xlab('') +
  facet_grid(~ Center, scale='free_x', space = 'free_x') +
  theme(legend.position = 'bottom', legend.title = element_blank(), , panel.background = element_blank(), text=element_text(size=16)) +
  guides(col=guide_legend(nrow=1)) +
  geom_point(aes(x=1:nrow(sampleqc3), y=median_rTiTv + 4*mad_rTiTv, col=Population2), shape='_') +
  geom_point(aes(x=1:nrow(sampleqc3), y=median_rTiTv - 4*mad_rTiTv, col=Population2), shape='_') +
  geom_vline(data=borderdat, aes(xintercept = borderpos), col='lightgrey') +
  scale_colour_manual(values=cols)
  #geom_segment(data=borderdat, aes(x=start, xend=borderpos, y=median_rTiTv + 4*mad_rTiTv, yend=median_rTiTv + 4*mad_rTiTv, col=race)) +
  #geom_segment(data=borderdat, aes(x=start, xend=borderpos, y=median_rTiTv - 4*mad_rTiTv, yend=median_rTiTv - 4*mad_rTiTv, col=race))
lb = geom_label(data=borderdat, aes(x=mid, y=2.10-(1:nrow(borderdat)%%4)*.01, label=label), col='black', angle=90, hjust=.5, alpha=.7)
ggsave(paste0('plots/',chrom,'/06_sample_qc_3_',chrom,'_titv.pdf'), p, width=wdt, height=hgh)
ggsave(paste0('plots/',chrom,'/06_sample_qc_3_',chrom,'_titv_labs.pdf'), p + lb, width=wdt, height=hgh)



# Het/Hom
p = ggplot(sampleqc3, aes(x=1:nrow(sampleqc3), y=sa.qc.rHetHomVar, colour=Population2)) +
  geom_point(size=.5)  +
  ylab('Het/Hom ratio') +
  xlab('') +
  facet_grid(~ Center, scale='free_x', space = 'free_x') +
  theme(legend.position = 'none', panel.background = element_blank(), text=element_text(size=16)) +
  guides(col=guide_legend(nrow=1)) +
  geom_point(aes(x=1:nrow(sampleqc3), y=median_rHetHomVar + 4*mad_rHetHomVar, col=Population2), shape='_', col='black') +
  geom_point(aes(x=1:nrow(sampleqc3), y=median_rHetHomVar - 4*mad_rHetHomVar, col=Population2), shape='_', col='black') +
  geom_vline(data=borderdat, aes(xintercept = borderpos), col='lightgrey') +
  #geom_segment(data=borderdat, aes(x=start, xend=borderpos, y=median_rHetHomVar + 4*mad_rHetHomVar, yend=median_rHetHomVar + 4*mad_rHetHomVar, col=Population2)) +
  #geom_segment(data=borderdat, aes(x=start, xend=borderpos, y=median_rHetHomVar - 4*mad_rHetHomVar, yend=median_rHetHomVar - 4*mad_rHetHomVar, col=Population2)) +
  scale_colour_manual(values=cols)
lb = geom_label(data=borderdat, aes(x=mid, y=3.2-(1:nrow(borderdat)%%4)*.2, label=label), col='black', angle=90, hjust=.5, alpha=.7)
ggsave(paste0('plots/',chrom,'/06_sample_qc_3_',chrom,'_hethom.pdf'), p, width=wdt, height=hgh)
ggsave(paste0('plots/',chrom,'/06_sample_qc_3_',chrom,'_hethom_labs.pdf'), p + lb, width=wdt, height=hgh)

#ggplot(sampleqc3, aes(sa.qc.rHetHomVar, fill=odd)) + geom_density() + facet_wrap( ~ label) + scale_fill_manual(values=cols) + xlab('Het/Hom ratio') + theme(legend.position = 'none', panel.background = element_blank())
#ggsave(paste0('plots/',chrom,'/01_sample_qc_',chrom,'_hethom_density.pdf'))

# ggplot(sampleqc3, aes(x=label, y=sa.qc.rHetHomVar, fill=odd)) +
#   geom_violin() +
#   facet_grid( ~ Center, scale='free_x', space = 'free_x'  ) +
#   scale_fill_manual(values=cols) +
#   theme(legend.position = 'none', panel.background = element_blank(), text=element_text(size=16), axis.text.x = element_text(angle=45, hjust=1)) +
#   geom_point(data=borderdat, aes(x=label, y=median_rHetHomVar - 4*mad_rHetHomVar), col='black', shape=95, size=5) +
#   geom_point(data=borderdat, aes(x=label, y=median_rHetHomVar + 4*mad_rHetHomVar), col='black', shape=95, size=5) +
#   xlab('') +
#   ylab('Het/Hom ratio')
# ggsave(paste0('plots/',chrom,'/01_sample_qc_',chrom,'_hethom_violin.pdf'), width=wdt, height=hgh)




# Ins/Del
p = ggplot(sampleqc3, aes(x=1:nrow(sampleqc3), y=sa.qc.rInsertionDeletion, colour=Population2)) +
  geom_point(size=.5)  +
  ylab('Insertion/Deletion ratio') +
  xlab('') +
  facet_grid(~ Center, scale='free_x', space = 'free_x') +
  theme(legend.position = 'none', panel.background = element_blank(), text=element_text(size=16)) +
  guides(col=guide_legend(nrow=1)) +
  geom_point(aes(x=1:nrow(sampleqc3), y=median_rInsertionDeletion + 4*mad_rInsertionDeletion, col=Population2), shape='_') +
  geom_point(aes(x=1:nrow(sampleqc3), y=median_rInsertionDeletion - 4*mad_rInsertionDeletion, col=Population2), shape='_') +
  geom_vline(data=borderdat, aes(xintercept = borderpos), col='lightgrey') +
  #geom_segment(data=borderdat, aes(x=start, xend=borderpos, y=median_rHetHomVar + 4*mad_rHetHomVar, yend=median_rHetHomVar + 4*mad_rHetHomVar, col=Population2)) +
  #geom_segment(data=borderdat, aes(x=start, xend=borderpos, y=median_rHetHomVar - 4*mad_rHetHomVar, yend=median_rHetHomVar - 4*mad_rHetHomVar, col=Population2)) +
  scale_colour_manual(values=cols)
lb = geom_label(data=borderdat, aes(x=mid, y=1.14-(1:nrow(borderdat)%%4)*.01, label=label), col='black', angle=90, hjust=.5, alpha=.7)
ggsave(paste0('plots/',chrom,'/06_sample_qc_3_',chrom,'_insdel.pdf'), p, width=wdt, height=hgh)
ggsave(paste0('plots/',chrom,'/06_sample_qc_3_',chrom,'_insdel_labs.pdf'), p + lb, width=wdt, height=hgh)

# ggplot(sampleqc3, aes(sa.qc.rHetHomVar, fill=odd)) + geom_density() + facet_wrap( ~ label) + scale_fill_manual(values=cols) + xlab('Insertion/Deletion ratio') + theme(legend.position = 'none', panel.background = element_blank())
# ggsave(paste0('plots/',chrom,'/01_sample_qc_',chrom,'_insdel_density.pdf'))
# 
# ggplot(sampleqc3, aes(x=label, y=sa.qc.rHetHomVar, fill=odd)) +
#   geom_violin() +
#   facet_grid( ~ Center, scale='free_x', space = 'free_x'  ) +
#   scale_fill_manual(values=cols) +
#   theme(legend.position = 'none', panel.background = element_blank(), text=element_text(size=16), axis.text.x = element_text(angle=45, hjust=1)) +
#   geom_point(data=borderdat, aes(x=label, y=median_rHetHomVar - 4*mad_rHetHomVar), col='black', shape=95, size=5) +
#   geom_point(data=borderdat, aes(x=label, y=median_rHetHomVar + 4*mad_rHetHomVar), col='black', shape=95, size=5) +
#   xlab('') +
#   ylab('Insertion/Deletion ratio')
# ggsave(paste0('plots/',chrom,'/01_sample_qc_',chrom,'_insdel_violin.pdf'), width=wdt, height=hgh)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# count and define samples to be excluded
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sampleqc4 = sampleqc3 %>%
  #select(Center, Project, Sample.ID, label) %>%
  ungroup() %>%
  transmute(Center, Project, Sample.ID, label,
            exc_callrate = sa.qc.callRate < mincallrate,
            exc_coverage = sa.qc.dpMean < mindepth,
            exc_titv = sa.qc.rTiTv > median_rTiTv + 4*mad_rTiTv |
              sa.qc.rTiTv < median_rTiTv - 4*mad_rTiTv,
            exc_hethom = sa.qc.rHetHomVar > median_rHetHomVar + 4*mad_rHetHomVar |
              sa.qc.rHetHomVar < median_rHetHomVar - 4*mad_rHetHomVar,
            exc_insdel = sa.qc.rInsertionDeletion > median_rInsertionDeletion + 4*mad_rInsertionDeletion |
              sa.qc.rInsertionDeletion > median_rInsertionDeletion + 4*mad_rInsertionDeletion,
            exc_any = exc_callrate | exc_coverage | exc_titv | exc_hethom | exc_insdel)

# summary project
exctab = sampleqc4 %>% group_by(label) %>%
  summarize(Center=Center[1],
            callrate = mean(exc_callrate),
            coverage = mean(exc_coverage),
            titv = mean(exc_titv),
            hethom = mean(exc_hethom),
            insdel = mean(exc_insdel),
            any = mean(exc_any))

exctab %>% gather(reason, exc_rate=callrate:any) %>%
  ggplot(aes(reason, label)) + geom_tile(aes(fill=value)) + scale_fill_viridis() + xlab('Reason for sample exclusion') + ylab('') + theme(text=element_text(size=16), panel.background = element_blank(), axis.text.x = element_text(angle=45, hjust=1))
ggsave(paste0('plots/',chrom,'/06_sample_qc_3_',chrom,'_summary_project.pdf'))

# summary center
exctab = sampleqc4 %>% group_by(Center) %>%
  summarize(
    callrate = mean(exc_callrate),
    coverage = mean(exc_coverage),
    titv = mean(exc_titv),
    hethom = mean(exc_hethom),
    insdel = mean(exc_insdel),
    any = mean(exc_any))

# exctab %>% gather(reason, exc_rate=callrate:any) %>%
#   ggplot(aes(reason, Center)) + geom_tile(aes(fill=value)) + scale_fill_viridis() + xlab('Reason for sample exclusion') + ylab('') + theme(text=element_text(size=16), panel.background = element_blank(), axis.text.x = element_text(angle=45, hjust=1))
# ggsave(paste0('plots/',chrom,'/01_sample_qc_',chrom,'_summary_project.pdf'))


sampleqc4 %>% filter(!exc_any) %>% select(Sample.ID) %>%
  write.table(paste0('~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/',chrom,'/06_sample_qc_3_keep.txt'), quote=F, col.names=F, row.names=F)

system(paste0(gsutil, ' cp ~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/',chrom,'/06_sample_qc_3_keep.txt gs://ccdg-qc-multi/qc_measures/',chrom,'/'))
