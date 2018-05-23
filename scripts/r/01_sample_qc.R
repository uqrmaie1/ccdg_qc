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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rename8 = read.table('../data/rename_eight.txt', stringsAsFactors = F, h=T)
sampleqc = read.table(paste0('~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/',chrom,'/sample_qc_info.txt'), h=T, stringsAsFactors = F, sep='\t')
manifest = read.table('~/Dropbox/postdoc/projects/ccdg_qc/CCDG_Freeze_1_Subsetting_Manifest_2017-12-06.txt', sep='\t', h=T, stringsAsFactors = F, quote='')
manifest %<>% filter(!Remove)

if('NA12892' %in% sampleqc[,1]) sampleqc[match(rename8[,1], sampleqc[,1]),1] = rename8[,2]

labs = as.data.frame(do.call(rbind, strsplit(manifest$Filename_Prefix, '_')))
dup = names(which(colSums(table(labs[,2], labs[,3]) > 0) > 1))
labs = paste0(ifelse(labs[,3] %in% dup, paste0(labs[,2], '_'), ''), labs[,3])

sampleqc2 = left_join(manifest, sampleqc, by=c('Sample.ID' = 's')) %>%
  mutate(borders = (1:n()) %in% cumsum(rle(Project)$lengths)) %>%
  mutate(label = labs) %>%
  arrange(Center, Project) %>%
  mutate(odd = as.numeric(factor(Project, levels=unique(Project)))%%2==0)
borderdat = sampleqc2 %>% slice(which(borders)) %>%
  mutate(borderpos=which(sampleqc2$borders)) %>%
  select(Center, Project, borderpos, label, odd) %>%
  mutate(start=c(1,borderpos[-length(borderpos)])) %>%
  mutate(mid=(start+borderpos)/2)
nums = rowSums(table(sampleqc2$Center, sampleqc2$Project) > 0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot call rate, coverage
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cols = c('steelblue3', 'seagreen4')
wdt = 20 
hgh = 5

# call rate
p = ggplot(sampleqc2, aes(x=1:nrow(sampleqc2), y=sample_qc.call_rate, colour=odd)) +
  geom_point()  + scale_size_manual(values=c(.5, 13)) +
  #geom_point(size=.5) +
  facet_grid( ~ Center, scale='free_x', space = "free_x") +
  theme(legend.position = 'none', panel.background = element_blank(), text=element_text(size=16)) + geom_hline(yintercept=mincallrate, linetype=1, size=.3) +
  scale_color_manual(values=cols) +
  ylab('Call rate') +
  xlab('') +
  geom_vline(data=borderdat, aes(xintercept = borderpos), col='lightgrey')
lb = geom_label(data=borderdat, aes(x=mid, y=1.0-(1:nrow(borderdat)%%4)*.004, label=label), col='black', angle=90, hjust=.5, alpha=.7)
ggsave(paste0('plots/',chrom,'/01_sample_qc_',chrom,'_callrate.pdf'), p, width=wdt, height=hgh)
ggsave(paste0('plots/',chrom,'/01_sample_qc_',chrom,'_callrate_labs.pdf'), p + lb, width=wdt, height=hgh)

# coverage
p = ggplot(sampleqc2, aes(x=1:nrow(sampleqc2), y=sample_qc.dp_mean, colour=odd)) +
  geom_point()  + scale_size_manual(values=c(.5, 13)) +
  #geom_point(size=.5) +
  facet_grid( ~ Center, scale='free_x', space = 'free_x') +
  theme(legend.position = 'none', panel.background = element_blank(), text=element_text(size=16)) +
  geom_hline(yintercept=mincoverage, linetype=1, size=.3) +
  scale_color_manual(values=cols) +
  ylab('Coverage') +
  xlab('') +
  geom_vline(data=borderdat, aes(xintercept = borderpos), col='lightgrey')
lb = geom_label(data=borderdat, aes(x=mid, y=60-(1:nrow(borderdat)%%4)*2, label=label), col='black', angle=90, hjust=1, alpha=.7)
ggsave(paste0('plots/',chrom,'/01_sample_qc_',chrom,'_coverage.pdf'), p, width=wdt, height=hgh)
ggsave(paste0('plots/',chrom,'/01_sample_qc_',chrom,'_coverage_labs.pdf'), p + lb, width=wdt, height=hgh)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# count and define samples to be excluded
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sampleqc4 = sampleqc2 %>%
  #select(Center, Project, Sample.ID, label) %>%
  ungroup() %>%
  transmute(Center, Project, Sample.ID, label,
            exc_callrate = sample_qc.call_rate < mincallrate,
            exc_coverage = sample_qc.dp_mean < mincoverage,
            exc_any = exc_callrate | exc_coverage)

# summary project
exctab = sampleqc4 %>% group_by(label) %>%
  summarize(Center=Center[1],
            callrate = mean(exc_callrate),
            coverage = mean(exc_coverage),
            any = mean(exc_any))

exctab %>% gather(reason, exc_rate=callrate:any) %>%
  ggplot(aes(reason, label)) + geom_tile(aes(fill=value)) + scale_fill_viridis() + xlab('Reason for sample exclusion') + ylab('') + theme(text=element_text(size=16), panel.background = element_blank(), axis.text.x = element_text(angle=45, hjust=1))
ggsave(paste0('plots/',chrom,'/01_sample_qc_',chrom,'_summary_project.pdf'))

# summary center
exctab = sampleqc4 %>% group_by(Center) %>%
  summarize(
    callrate = mean(exc_callrate),
    coverage = mean(exc_coverage),
    any = mean(exc_any))

exctab %>% gather(reason, exc_rate=callrate:any) %>%
  ggplot(aes(reason, Center)) + geom_tile(aes(fill=value)) + scale_fill_viridis() + xlab('Reason for sample exclusion') + ylab('') + theme(text=element_text(size=16), panel.background = element_blank(), axis.text.x = element_text(angle=45, hjust=1))
ggsave(paste0('plots/',chrom,'/01_sample_qc_',chrom,'_summary_project.pdf'))


sampleqc4 %>% filter(!exc_any) %>% select(Sample.ID) %>%
  write.table(paste0('~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/',chrom,'/01_sample_qc_keep.txt'), quote=F, col.names=F, row.names=F)

system(paste0(gsutil, ' cp ~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/',chrom,'/01_sample_qc_keep.txt gs://ccdg-qc-multi/qc_measures/',chrom,'/01_sample_qc_keep.txt'))


