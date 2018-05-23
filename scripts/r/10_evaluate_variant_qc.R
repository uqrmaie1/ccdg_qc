# !diagnostics off

library(plyr)
library(tidyverse)
library(magrittr)
library(viridis)

setwd('~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/')

chrom = 'chr20'
chrom = 'onep'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# copy files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gsutil = '/Users/robert/google-cloud-sdk/bin/gsutil'
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/',chrom,'/variants_table_pre.tsv.gz ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/',chrom,'/variants_table_post.tsv.gz ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/mendel/',chrom,'/all1.txt ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/mendel/',chrom,'/all2.txt ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/mendel/',chrom,'/perfam1.txt ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/mendel/',chrom,'/perfam2.txt ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/mendel/',chrom,'/persample1.txt ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/mendel/',chrom,'/persample2.txt ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/mendel/',chrom,'/pervariant1.txt ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/mendel/',chrom,'/pervariant2.txt ~/Downloads/'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pre = as.data.frame(read_table2('~/Downloads/variants_table_pre.tsv.gz'))
post = as.data.frame(read_table2('~/Downloads/variants_table_post.tsv.gz'))

comb = bind_rows(data.frame(stage='pre', pre, stringsAsFactors = F),
                 data.frame(stage='post', post, stringsAsFactors = F))
comblong = comb %>% gather('metric', 'value', qual,qc.callRate:qc.pHWE)


all1 = read.table('~/Downloads/all1.txt', sep='\t', h=T, stringsAsFactors = F)
all2 = read.table('~/Downloads/all2.txt', sep='\t', h=T, stringsAsFactors = F)
persample1 = as.data.frame(read_table2('~/Downloads/persample1.txt'))
persample2 = as.data.frame(read_table2('~/Downloads/persample2.txt'))
pervariant1 = as.data.frame(read_table2('~/Downloads/pervariant1.txt'))
pervariant2 = as.data.frame(read_table2('~/Downloads/pervariant2.txt'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot qc measures pre and post
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

vsplit = do.call(rbind, strsplit(comb$v, ':'))
varlen = pmax(nchar(vsplit[,3]), nchar(vsplit[,4]))
comb$varlen = varlen
comb$issnp = varlen == 1

comblong %>% sample_n(1e4) %>% ggplot(aes(value, fill=stage)) + geom_density(alpha=.5) + facet_wrap(~ metric, scales='free') + scale_x_log10()
numsamp = 1e5

# AC
p = bind_rows(comb %>% filter(stage=='pre') %>% sample_n(numsamp), comb %>% filter(stage=='post') %>% sample_n(numsamp)) %>% ggplot(aes(qc.AC, fill=stage)) + geom_histogram(position='dodge') + scale_x_log10(breaks=10^(0:4)) + theme(panel.background = element_blank(), axis.line = element_line(), text=element_text(size=18), legend.title = element_blank()) + xlab('allele count')
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_ac_histogram.pdf'), p)

# AF
p = comb %>% sample_n(numsamp) %>% ggplot(aes(qc.AF, col=stage)) + stat_ecdf(geom='step') + scale_x_log10(breaks=c(.0001, .001, .01, .1)) + theme(panel.background = element_blank(), axis.line = element_line(), text=element_text(size=18), legend.title = element_blank()) + xlab('minor allele frequency') + ylab('cumulative fraction of variants')
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_maf_ecdf.pdf'), p)

# hethomvar ratio
p = comb %>% sample_n(numsamp) %>% ggplot(aes(qc.nHet/qc.nHomVar, fill=stage)) + geom_density(alpha=.4, adjust=1, trim=T) + scale_x_log10() + theme(panel.background = element_blank(), axis.line = element_line(), text=element_text(size=18), legend.title = element_blank())
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_hethomvar_density.pdf'), p)

# dpMean
p = comb %>% sample_n(numsamp) %>% ggplot(aes(qc.dpMean, fill=stage)) + geom_density(alpha=.4) + scale_x_continuous(limits=c(10,40)) + theme(panel.background = element_blank(), axis.line = element_line(), text=element_text(size=18), legend.title = element_blank()) + xlab('dpMean')
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_dpMean_density.pdf'), p)

# dpStDev
p = comb %>% sample_n(numsamp) %>% ggplot(aes(qc.dpStDev, fill=stage)) + geom_density(alpha=.4,) + scale_x_continuous(limits=c(5,15)) + theme(panel.background = element_blank(), axis.line = element_line(), text=element_text(size=18), legend.title = element_blank()) + xlab('dpStDev')
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_dpStDev_density.pdf'), p)

# gqMean
p = comb %>% sample_n(numsamp) %>% ggplot(aes(qc.gqMean, fill=stage)) + geom_density(alpha=.4) + scale_x_continuous() + theme(panel.background = element_blank(), axis.line = element_line(), text=element_text(size=18), legend.title = element_blank()) + xlab('gqMean')
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_gqMean_density.pdf'), p)

# gqStDev
p = comb %>% sample_n(numsamp) %>% ggplot(aes(qc.gqStDev, fill=stage)) + geom_density(alpha=.4) + scale_x_continuous(limits=c(10,40)) + theme(panel.background = element_blank(), axis.line = element_line(), text=element_text(size=18), legend.title = element_blank()) + xlab('gqStDev')
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_gqStDev_density.pdf'), p)

# pHWE
p = comb %>% sample_n(numsamp) %>% ggplot(aes(qc.pHWE, fill=stage)) + geom_density(alpha=.4) + scale_x_continuous() + theme(panel.background = element_blank(), axis.line = element_line(), text=element_text(size=18), legend.title = element_blank()) + xlab('pHWE')
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_pHWE_density.pdf'), p)

comb %>% filter(qc.nHomVar > 0) %>% sample_n(numsamp) %>% ggplot(aes(qc.pHWE, fill=stage)) + geom_density(alpha=.4) + scale_x_continuous() + theme(panel.background = element_blank(), axis.line = element_line(), text=element_text(size=18), legend.title = element_blank()) + xlab('pHWE')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot mendel errors
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


prem = left_join(pre, pervariant1, by='v') %>% mutate(nError = replace(nError, is.na(nError), 0))
postm = left_join(post, pervariant2, by='v') %>% mutate(nError = replace(nError, is.na(nError), 0))

comb = bind_rows(data.frame(stage='pre', prem, stringsAsFactors = F),
                 data.frame(stage='post', postm, stringsAsFactors = F))

breaks = rev(10^(seq(0,log10(unique(sort(post$qc.AF))[2]),len=100)))
comb %<>% mutate(afbin = cut(qc.AF, breaks))

mafbinsum = comb %>% group_by(stage, afbin) %>% summarise(meanmend = mean(nError), count=n())

cols = c('orange', 'red', 'green', 'blue')
x = c(4, 21, 43, 64)
txtdat = data.frame(x=x, lab=c('0.01%' ,'0.1%', '1%', '10%'), afbin=NA, pre=NA, post=NA, stage=NA, stringsAsFactors = F)
p = ggplot(mafbinsum, aes(afbin, meanmend, shape=stage)) + geom_point() + theme(axis.text.x=element_text(angle=90, size=8), panel.background = element_blank(), axis.line=element_line(), text=element_text(size=18), legend.title = element_blank()) + xlab('allele frequency bin') + ylab('mean mendel errors per variant') + geom_vline(xintercept=x, col=cols, linetype=2) + geom_text(data=txtdat, aes(x, .5, label=lab, col=as.factor(x)), hjust=-.01, size=6, show.legend = F) + scale_colour_manual(values=cols)
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_mendelerror_maf1.pdf'), p)



dat = spread(mafbinsum %>% select(-count), stage, meanmend) %>% filter(post != 0)
p = ggplot(dat, aes(afbin, pre/post)) + geom_point() + theme(axis.text.x=element_text(angle=90, size=8), panel.background = element_blank(), axis.line=element_line(), text=element_text(size=18), legend.position='none') + xlab('allele frequency bin') + ylab('reduction of mendel errors') + scale_y_log10() + geom_vline(xintercept=x, col=cols, linetype=2) + geom_text(data=txtdat, aes(x, 5, label=lab, col=as.factor(x)), hjust=-.01, size=6) + scale_colour_manual(values=cols)
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_mendelerror_maf2.pdf'), p)



dat = spread(mafbinsum %>% select(-meanmend), stage, count) #%>% filter(post != 0)
p = ggplot(dat, aes(afbin, pre/post)) + geom_point() + theme(axis.text.x=element_text(angle=90, size=8), panel.background = element_blank(), axis.line=element_line(), text=element_text(size=18), legend.position='none') + xlab('allele frequency bin') + ylab('reduction of variants') + geom_vline(xintercept=x, col=cols, linetype=2) + geom_text(data=txtdat, aes(x, .5, label=lab, col=as.factor(x)), hjust=-.01, size=6) + scale_colour_manual(values=cols) + scale_y_continuous(limits=c(0,2))
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_maf_prepost.pdf'), p)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mendel error types
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

breaks = rev(10^(seq(0,log10(unique(sort(pre$qc.AF))[2]),len=100)))
transitions = c("A/A x A/A -> A/G", "T/T x T/T -> T/C", "C/C x C/C -> C/T", "G/G x G/G -> G/A")
transversions = c("G/G x G/G -> G/C", "C/C x C/C -> C/G", "T/T x T/T -> T/G", "C/C x C/C -> C/A", "T/T x T/T -> T/A", "A/A x A/A -> A/T", "A/A x A/A -> A/C", "G/G x G/G -> G/T")

mendel = bind_rows(all1 %>% select(v, code, error) %>% left_join(pre %>% select(v, qc.AF), by='v') %>% mutate(stage='pre'),
                   all2 %>% select(v, code, error) %>% left_join(pre %>% select(v, qc.AF), by='v') %>% mutate(stage='post')) %>%
  mutate(afbin = cut(qc.AF, breaks), trans=ifelse(error %in% transitions, 'transition', ifelse(error %in% transversions, 'transversion', NA)))

summ = mendel %>% group_by(afbin, stage, code) %>% summarise(meanmend = n())
ggplot(summ, aes(afbin, meanmend, col=stage)) + geom_point() + facet_grid(stage ~ as.factor(code), scales='free') + theme(axis.text.x = element_blank())
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_mendel_errors_by_type.pdf'))


a = all1 %>% filter(code == 2 & nchar(error) == 16)

# transition transversion rate pre and post qc
mendel %>% filter(!duplicated(mendel) & code==2 & as.numeric(afbin)==1) %>% group_by(stage) %>% summarise(is = mean(trans == 'transition', na.rm=T), ve = mean(trans=='transversion', na.rm=T))
# A tibble: 2 x 3
# stage    is    ve
# <chr> <dbl> <dbl>
#   1 post  0.638 0.362
# 2 pre   0.458 0.542
# > 0.638/0.362
# [1] 1.762431

#singletons:
# > 0.570/0.430
# [1] 1.325581

# annotate by dp
predp = read_table2('~/Downloads/gtable_pre.txt.gz')
postdp = read_table2('~/Downloads/gtable_post.txt.gz')

ab = left_join(all2, perfam2[,1:3], by='fid') %>% left_join(predp, by=c('v'='v', 'father'='s')) %>% left_join(predp, by=c('v'='v', 'mother'='s')) %>% mutate(dpmin=pmin(DP.x, DP.y)) %>% left_join(select(prem, v, qc.AF)) %>% mutate(afbin = cut(qc.AF, breaks), dpminbin=cut(dpmin, quantile(ab$dpmin, seq(0,1,len=6))))

summ = ab %>% group_by(afbin, dpminbin) %>% summarise(numerr = n())
ggplot(summ, aes(afbin, numerr, col=dpminbin)) + geom_point() + theme(panel.background = element_blank(), axis.text.x = element_blank()) + ylab('# of MEs after QC') + xlab('AF bin')

ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_mendel_errors_by_maf_and_dpmin.pdf'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot mendel errors v3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

vvcfids = left_join(all2, manifest, by=c('s'='Sample.ID')) %>% filter(Filename_Prefix == 'ccdg_washu_vccontrols_gru') %$% unique(fid) 
all1novvc = all1 %>% filter(!fid %in% vvcfids)
all2novvc = all2 %>% filter(!fid %in% vvcfids)

errormode_pre = as.data.frame(table(all1novvc$v, all1novvc$code), stringsAsFactors = F) %>% spread(Var2, Freq) %>% select(Var1, as.character(1:7))
errormode_post = as.data.frame(table(all2novvc$v, all2novvc$code), stringsAsFactors = F) %>% spread(Var2, Freq) %>% select(Var1, as.character(1:7))

names(errormode_pre) = c('v', paste0('me', names(errormode_pre)[-1]))
names(errormode_post) = c('v', paste0('me', names(errormode_post)[-1]))

prem = left_join(pre, errormode_pre, by='v'); for(i in 1:7) {prem[is.na(prem[,paste0('me', i)]),paste0('me', i)] = 0}
postm = left_join(post, errormode_post, by='v'); for(i in c(1:7)) {postm[is.na(postm[,paste0('me', i)]),paste0('me', i)] = 0}


vsplit = do.call(rbind, strsplit(prem$v, ':'))
varlen = pmax(nchar(vsplit[,3]), nchar(vsplit[,4]))
prem$varlen = varlen
prem$issnp = varlen == 1
vsplit = do.call(rbind, strsplit(postm$v, ':'))
varlen = pmax(nchar(vsplit[,3]), nchar(vsplit[,4]))
postm$varlen = varlen
postm$issnp = varlen == 1



comb = bind_rows(data.frame(stage='pre', prem, stringsAsFactors = F),
                 data.frame(stage='post', postm, stringsAsFactors = F))


breaks = rev(10^(seq(0,log10(unique(sort(pre$qc.AF))[2]*.99),len=100)))
comb %<>% filter(qc.AF > 0)
comb %<>% mutate(afbin = cut(qc.AF, breaks))
comb %<>% gather(error_mode, mecount, me1:me7)


mafbinsum = comb %>% filter(issnp) %>%  group_by(stage, afbin, error_mode) %>% summarise(meanmend = mean(mecount), summend=sum(mecount), count=n())
mafbinsum$stage2 = factor(mafbinsum$stage, levels=levels(as.factor(mafbinsum$stage))[2:1])
mafbinsum %<>% mutate(pertriopergenome = summend*100/(179-24))

namemap = c('me1' = 'HV/HV -> Het',
            'me2' = 'HR/HR -> Het',
            'me3' = 'HR/!HR -> HV',
            'me4' = '!HR/HR -> HV',
            'me5' = 'HR/HR -> HV',
            'me6' = '!HV/HV -> HR',
            'me7' = 'HV/!HV -> HR')

mafbinsum$error_mode2 = namemap[match(mafbinsum$error_mode, names(namemap))]

txtdat = mafbinsum %>% group_by(stage, error_mode2) %>% summarise(sum=sum(pertriopergenome)) %>% filter(stage=='post')
txtdat %<>% mutate(stage2 = factor('post', levels=levels(mafbinsum$stage2)))
#levels(txtdat$stage2) = rev(levels(mafbinsum$stage2))

p = ggplot(mafbinsum, aes(as.numeric(afbin), pertriopergenome, shape=stage2)) + geom_point() + theme(axis.text.x=element_blank(), axis.line=element_line(), text=element_text(size=18), legend.title = element_blank(), axis.ticks.x = element_blank()) + xlab('allele frequency bin') + ylab('Mendel errors per trio')  + scale_colour_manual(values=cols) + facet_grid(stage2 ~ error_mode2, scale='free_y') + geom_text(data=txtdat, aes(50, 10, label=round(sum)), size=8)
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_mendelerror_maf4_notsnps.pdf'), p, width=12, height=8)

badfams = perfam2 %>% filter(nErrors > 20) %>% select(fid)
badvars = all2 %>% filter(!fid %in% badfams$fid) %>% select(v) %$% unique(v)
badids = badfams %>% left_join(all2, by='fid') %>% select(s) 


cols = c('orange', 'red', 'green', 'blue', 'purple')
x = c(4, 21, 43, 64, 93)
txtdat = data.frame(x=x, lab=c('0.01%' ,'0.1%', '1%', '10%', '50%'), cols=cols, afbin=NA, pre=NA, post=NA, stage=NA, error_mode=NA, stringsAsFactors = F)
ggplot(mafbinsum %>% filter(error_mode=='me2'), aes(afbin, meanmend, shape=stage2)) + geom_point() + theme(axis.text.x=element_blank(), panel.background = element_blank(), axis.line=element_line(), text=element_text(size=18), legend.title = element_blank()) + xlab('allele frequency bin') + ylab('mean mendel errors per variant')  + scale_colour_manual(values=cols) + facet_grid(stage2 ~ error_mode, scale='free_y') #+ geom_vline(data=a, aes(xintercept=x, col=cols), linetype=2) #+ geom_text(data=txtdat, aes(x, .5, label=lab, col=as.factor(x)), hjust=-.01, size=6, show.legend = F)



ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_mendelerror_maf1.pdf'), p)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot mendel errors v4
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

left_join(all2, manifest, by=c('s'='Sample.ID')) %>%
  group_by(Filename_Prefix, s, code) %>%
  summarise(n=n()) %>%
  left_join(predicted_populations[,2:4], by=c('s'='Sample')) %>%
  filter(code != 5 & code != 10) %>% mutate(code2=paste(code, namemap[paste0('me', code)])) %>%
  ggplot(aes(Filename_Prefix, n, col=Population)) + geom_jitter(width=.2) + facet_wrap(~ code2, scales='free_x') + scale_colour_manual(values=cols) + xlab('') + ylab('errors per trio') + theme(axis.text.x=element_text(angle=90), panel.background = element_blank(), panel.grid.major.y = element_line(colour='black')) + coord_flip()
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_mendelerror_maf5.pdf'))



left_join(all2, manifest, by=c('s'='Sample.ID')) %>%
  group_by(Filename_Prefix, s) %>%
  summarise(n=n()) %>%
  left_join(predicted_populations[,2:4], by=c('s'='Sample')) %>%
  ggplot(aes(Filename_Prefix, n, col=Population)) + geom_jitter(width=.2) + scale_colour_manual(values=cols) + xlab('') + ylab('errors per trio') + theme(axis.text.x=element_text(angle=90), panel.background = element_blank(), panel.grid.major.y = element_line(colour='black')) + coord_flip()
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_mendelerror_maf6.pdf'))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# all variants
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# copy files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

chrom = 'trios'

gsutil = '/Users/robert/google-cloud-sdk/bin/gsutil'
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/',chrom,'/variants_table_pre.tsv.gz ~/Downloads/mendel_trios/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/',chrom,'/variants_table_post.tsv.gz ~/Downloads/mendel_trios/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/mendel/',chrom,'/all1.txt ~/Downloads/mendel_trios/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/mendel/',chrom,'/all2.txt ~/Downloads/mendel_trios/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/mendel/',chrom,'/perfam1.txt ~/Downloads/mendel_trios/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/mendel/',chrom,'/perfam2.txt ~/Downloads/mendel_trios/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/mendel/',chrom,'/persample1.txt ~/Downloads/mendel_trios/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/mendel/',chrom,'/persample2.txt ~/Downloads/mendel_trios/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/mendel/',chrom,'/pervariant1.txt ~/Downloads/mendel_trios/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/mendel/',chrom,'/pervariant2.txt ~/Downloads/mendel_trios/'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


namemap2 = c('HomVar, HomVar -> Het',
            'HomRef, HomRef -> Het',
            'HomRef, !HomRef -> HomVar',
            '!HomRef, HomRef -> HomVar',
            'HomRef, HomRef -> HomVar',
            '!HomVar, HomVar -> HomRef',
            'HomVar, !HomVar -> HomRef',
            'HomVar, HomVar -> HomRef')

namemap3 = c('HomVar, HomVar -> Het',
             'HomRef, HomRef -> Het',
             'HomRef -> HomVar',
             'HomRef -> HomVar',
             'HomRef -> HomVar',
             'HomVar -> HomRef',
             'HomVar -> HomRef',
             'HomVar -> HomRef')


pre = as.data.frame(read_table2('~/Downloads/mendel_trios/variants_table_pre_freqonly.tsv'))
post = as.data.frame(read_table2('~/Downloads/mendel_trios/variants_table_post_freqonly.tsv'))

comb = bind_rows(data.frame(stage='pre', pre, stringsAsFactors = F),
                 data.frame(stage='post', post, stringsAsFactors = F))

all1 = read.table('~/Downloads/mendel_trios/all1.txt', sep='\t', h=T, stringsAsFactors = F)
all2 = read.table('~/Downloads/mendel_trios/all2.txt', sep='\t', h=T, stringsAsFactors = F)
perfam1 = as.data.frame(read_table2('~/Downloads/mendel_trios/perfam1.txt'))
perfam2 = as.data.frame(read_table2('~/Downloads/mendel_trios/perfam2.txt'))
persample1 = as.data.frame(read_table2('~/Downloads/mendel_trios/persample1.txt'))
persample2 = as.data.frame(read_table2('~/Downloads/mendel_trios/persample2.txt'))
pervariant1 = as.data.frame(read_table2('~/Downloads/pervariant1.txt'))
pervariant2 = as.data.frame(read_table2('~/Downloads/pervariant2.txt'))


vsplit = do.call(rbind, strsplit(all1$v, ':'))
varlen = pmax(nchar(vsplit[,3]), nchar(vsplit[,4]))

mall1 = left_join(all1, pre, by='v') %>%
  mutate(issnp = varlen==1) %>%
  #filter(fam_id %in% (filter(perfam1, children==1) %$% fam_id)) %>%
  left_join(select(predicted_populations, Sample, Population, super), by=c('s'='Sample')) %>%
  left_join(select(manifest, Sample.ID, Filename_Prefix), by=c('s'='Sample.ID'))

vsplit = do.call(rbind, strsplit(all2$v, ':'))
varlen = pmax(nchar(vsplit[,3]), nchar(vsplit[,4]))

mall2 = left_join(all2, post, by='v') %>%
  mutate(issnp = varlen==1) %>%
  #filter(fam_id %in% (filter(perfam2, children==1) %$% fam_id)) %>%
  left_join(select(predicted_populations, Sample, Population, super), by=c('s'='Sample')) %>%
  left_join(select(manifest, Sample.ID, Filename_Prefix), by=c('s'='Sample.ID'))


mall = mall1 %>% mutate(stage='pre') %>% bind_rows(mutate(mall2, stage='post'))
breaks = rev(10^(seq(0,log10(min(mall$qc.AF)),len=100)))
mall %<>% mutate(afbin = cut(qc.AF, breaks))
mall %<>% mutate(code3=namemap3[code])
#occursinoneperson = colSums(table(mall$s, mall$v)); occursinoneperson = names(occursinoneperson)[occursinoneperson == 1]
mall %<>% mutate(singleton = qc.AC == 1)

# collapse variants
errpersample = mall %>% filter(code <= 8) %>% group_by(stage, s, code, code3, singleton, issnp, Population, super, Filename_Prefix) %>% summarise(nerrors = n())
errpersample %<>% ungroup %>% complete(nesting(Filename_Prefix, s, super, Population), nesting(code, code3), stage, singleton, issnp, fill = list(nerrors = 0))

errpersample %>% filter(singleton, issnp, Filename_Prefix=='ccdg_baylor_sol_hmb-col') %>% head() 

# collapse error codes
errpersample2 = errpersample %>% group_by(stage, s, code3, Population, super, Filename_Prefix, singleton, issnp) %>% summarise(nerrors = sum(nerrors))

# collapse error codes and cohorts
errpersample3 = errpersample2 %>% group_by(s, singleton, issnp, stage, Population) %>% summarise(nerrors = sum(nerrors))
errpersample3 %<>% mutate(singleton2 = ifelse(singleton, ' Singleton', 'not Singleton'),
                          issnp2 = ifelse(issnp, ' SNV', 'other'))

ggplot(errpersample3, aes(stage, nerrors)) + geom_violin() + geom_jitter(alpha=.3, size=.5, width=.2, aes(col=Population)) + facet_wrap(issnp2 ~ singleton2, scales='free') + coord_flip() + scale_colour_manual(values=cols) + theme(panel.background = element_blank()) + xlab('') + ylab('Mendelian errors per trio')
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_mendelerror_7.pdf'), height=5, width=7)

# filter post, singleton, SNVs; collapse cohorts
errpersample4 = errpersample2 %>% filter(issnp, stage=='post') %>% group_by(s, singleton, Population, code3) %>% summarise(nerrors = sum(nerrors)) %>% mutate(singleton2 = ifelse(singleton, ' Singleton', 'not Singleton'))

ggplot(errpersample4, aes('a', nerrors)) + geom_violin() + geom_jitter(alpha=.3, size=.5, width=.2, aes(col=Population)) + facet_wrap(singleton2 ~ code3, nrow=4, dir='v' ) + coord_flip() + scale_colour_manual(values=cols) + theme(panel.background = element_blank(), axis.text.y = element_blank()) + xlab('') + ylab('Mendelian errors per trio')
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_mendelerror_8.pdf'), height=5, width=7)


# filter post, singleton, SNVs, HomRef, HomRef -> Het
errpersample5 = errpersample2 %>% filter(issnp, singleton, stage=='post', code3=='HomRef, HomRef -> Het') %>%  mutate(singleton2 = ifelse(singleton, ' Singleton', 'not Singleton'))
ggplot(errpersample5, aes(Filename_Prefix, nerrors)) + geom_violin() + geom_jitter(alpha=.3, size=.5, width=.2, aes(col=Population)) + coord_flip() + scale_colour_manual(values=cols) + theme(panel.background = element_blank()) + xlab('') + ylab('Mendelian errors per trio') #+ geom_hline(yintercept=c(60,80), linetype=2)
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_mendelerror_9.pdf'), height=5, width=7)

ggplot(errpersample5, aes(Filename_Prefix, nerrors)) + geom_violin() + geom_jitter(alpha=.3, size=.5, width=.2, aes(col=Population)) + coord_flip(ylim=c(0,150)) + scale_colour_manual(values=cols) + theme(panel.background = element_blank()) + xlab('') + ylab('Mendelian errors per trio')  + geom_hline(yintercept=c(60,80), linetype=2)
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_mendelerror_10.pdf'), height=5, width=7)

errpersample35 = errpersample2 %>% filter(Filename_Prefix != 'ccdg_broad_mcgovern_gru' & Filename_Prefix != 'ccdg_washu_vccontrols_gru') %>% group_by(s, singleton, issnp, stage, Population) %>% summarise(nerrors = sum(nerrors))
errpersample35 %<>% mutate(singleton2 = ifelse(singleton, ' Singleton', 'not Singleton'),
                          issnp2 = ifelse(issnp, ' SNV', 'other'))
ggplot(errpersample35, aes(stage, nerrors)) + geom_violin() + geom_jitter(alpha=.3, size=.5, width=.2, aes(col=Population)) + facet_wrap(issnp2 ~ singleton2, scales='free') + coord_flip() + scale_colour_manual(values=cols) + theme(panel.background = element_blank()) + xlab('') + ylab('Mendelian errors per trio')
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_mendelerror_11.pdf'), height=5, width=7)






a = left_join(all2, post, by='v') %>%
  #filter(fam_id %in% (filter(perfam2, children==1) %$% fam_id)) %>%
  #left_join(select(predicted_populations, Sample, Population, super), by=c('s'='Sample')) %>%
  left_join(select(manifest, Sample.ID, Filename_Prefix), by=c('s'='Sample.ID')) %>% group_by(s, Filename_Prefix) %>% summarise(n=n())

mall %>% select(s, Filename_Prefix) %>% group_by(s, Filename_Prefix) %>% summarise(n=n()) %$% table(Filename_Prefix)

#%$% table(Filename_Prefix)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Evaluate TDT
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'gsutil cp gs://ccdg-qc-multi/qc_measures/tdt/onep/tdt_pre_table.txt ~/Downloads/'
'gsutil cp gs://ccdg-qc-multi/qc_measures/tdt/onep/tdt_post_table.txt ~/Downloads/'
mincount = 10

pre = read.table('~/Downloads/tdt_pre_table.txt', h=T, stringsAsFactors = F)
post = read.table('~/Downloads/tdt_post_table.txt', h=T, stringsAsFactors = F)

dat = bind_rows(data.frame(stage='pre', pre, stringsAsFactors = F), data.frame(stage='post', post, stringsAsFactors = F))
dat %<>% mutate(ratio = t/(t+u), snv=nchar(alleles) == 9)


dat %>% filter(u > mincount | t > mincount) %>% ggplot(aes(ratio)) + geom_density(aes(fill = stage), alpha=.5)

ppre = dat %>% filter(u > mincount | t > mincount, stage=='pre', snv) %$% ggqq(p_value, reduce=0, title='pre QC SNV')
ppost = dat %>% filter(u > mincount | t > mincount, stage=='post', snv) %$% ggqq(p_value, reduce=0, title='post QC SNV')
ppre2 = dat %>% filter(u > mincount | t > mincount, stage=='pre', !snv) %$% ggqq(p_value, reduce=0, title='pre QC other')
ppost2 = dat %>% filter(u > mincount | t > mincount, stage=='post', !snv) %$% ggqq(p_value, reduce=0, title='post QC other')

pdf(paste0('plots/',chrom,'/10_evaluate_variant_qc_tdt_qq.pdf'), width=10, height=8)
grid.arrange(ppre, ppost, ppre2, ppost2, nrow=2)
dev.off()






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Evaluate batch GWAS
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

chrom = 'onep'
minac = 10

system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/',chrom,'/variant_qc_batch_gwas.txt.gz ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/',chrom,'/variant_qc_batch_sample_files/variant_qc_batch_gwas_*tsv ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/', chrom, '/variant_qc_table.txt.gz ~/Downloads/'))

batchgwas = as.data.frame(read_table2('~/Downloads/variant_qc_batch_gwas.txt.gz'))
qc_table = as.data.frame(read_table2('~/Downloads/variant_qc_table.txt.gz'))

batch_sample_files = list()
for(f in list.files(path='~/Downloads/', pattern='variant_qc_batch_gwas_.+tsv')) {
  batch_sample_files[[f]] = read.table(paste0('~/Downloads/', f), h=T, stringsAsFactors = F)
}

bg = batchgwas %>% select(locus, alleles, info.AF, grep('linreg.+\\.AC|\\.beta|\\.y_transpose_x|\\.p_value', names(batchgwas)))
bg2 = bg %>% gather('name', 'value', -locus, -alleles, -info.AF)
bg2 = bg2 %>% mutate(value = as.numeric(gsub('[]"[]', '', value)),
                   pop=gsub('_.+', '', gsub('linreg_', '', name)),
                   coh = gsub('.+_ccdg', 'ccdg', gsub('\\..+', '', name)),
                   metric = gsub('.+\\.', '', name),
                   popcoh = paste(pop, coh, sep='_')) %>% select(-name)
bg3 = bg2 %>% spread(metric, value) %>% filter(!is.na(p_value))
bg3 = bg3 %>% mutate(trustpvalue = y_transpose_x > minac | AC - y_transpose_x > minac)

#bg4 = bg3 %>% mutate(pass = paste(locus, alleles) %in% paste(post$locus, post$alleles))
bg4 = bg3 %>% left_join(qc_table, by=c('locus', 'alleles')) %>% mutate_at(vars(qc.step0:qcpass), as.logical)

bg5 = bg4 %>% filter(trustpvalue)
bg5$qcpass = bg5$qccum.step4
bg5 %<>% mutate(indel = nchar(alleles) > 9)

pops = unique(bg5$popcoh)
i=2
a = bg4 %>% filter(popcoh == pops[i], !pass)
ggqq(a$p_value[a$trustpvalue], reduce = .9)

plots = list()
for(p in pops) {
  for(steps in c(paste0('qccum.step', 1:4), 'qcpass')) {
    subs = bg5 %>% filter(popcoh == p) %>% filter_(steps)
    plots[[length(plots)+1]] = ggqq(subs$p_value, reduce = .9, title=paste(p, steps))
  }
}

png(paste0('plots/',chrom,'/10_evaluate_variant_qc_batch_qq.png'), height=1600, width=2400)
do.call('grid.arrange', c(plots, nrow=5, as.table=FALSE))
dev.off()

breaks = rev(10^(seq(0,log10(min(bg5$variant_qc.AF)),len=10)))
bg5 %<>% mutate(afbin = cut(variant_qc.AF, breaks))
                #qccumstep = ifelse(!qc.step0, 0, ifelse(!qc.step1, 1, ifelse(!qc.step2, 2, ifelse(!qc.step3, 3, ifelse(!qc.step4, 4, 5))))))

a = bind_rows(bg5 %>% filter(qccum.step0) %>% group_by(popcoh, afbin, indel) %>% summarize(n=n(), lambdagc = get.lambda(p_value)) %>% mutate(step=0),
              bg5 %>% filter(qccum.step1) %>% group_by(popcoh, afbin, indel) %>% summarize(n=n(), lambdagc = get.lambda(p_value)) %>% mutate(step=1),
              bg5 %>% filter(qccum.step2) %>% group_by(popcoh, afbin, indel) %>% summarize(n=n(), lambdagc = get.lambda(p_value)) %>% mutate(step=2),
              bg5 %>% filter(qccum.step3) %>% group_by(popcoh, afbin, indel) %>% summarize(n=n(), lambdagc = get.lambda(p_value)) %>% mutate(step=3),
              bg5 %>% filter(qccum.step4) %>% group_by(popcoh, afbin, indel) %>% summarize(n=n(), lambdagc = get.lambda(p_value)) %>% mutate(step=4))

a %>% filter(!indel) %>%
ggplot(aes(afbin, lambdagc, fill=as.factor(step))) + geom_bar(position='dodge', stat='identity') + facet_grid(step ~ popcoh, scales='free') + geom_hline(yintercept=1) + geom_text(aes(afbin, lambdagc, label=n)) + theme(axis.text.x=element_text(angle=90, hjust=0))
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_batch_lambdagc.pdf'), width=18, heigh=12)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Evaluate batch GWAS
# compare AC accross batches
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/', chrom, '/variant_qc_ac_files/variant_qc_batch_ac_*EUR* ~/Downloads/'))

aclist = list()
for(f in list.files(path='~/Downloads/', pattern='variant_qc_batch_ac_controlsonly_AFR.+')) {
  aclist[[f]] = read_table2(paste0('~/Downloads/', f))
}

# aclist = list()
# for(f in list.files(path='~/Downloads/', pattern='variant_qc_batch_ac_controlsonly_AFR.+')) {
#   aclist[[f]] = read_table2(paste0('~/Downloads/', f))
# }

num_samples = sapply(aclist, function(x) max(x$homref + x$het + x$homvar))
use = which(num_samples > 50)
nams = gsub('\\.tsv\\.gz', '', gsub('.+ccdg_', '', names(use)))

n = 100000
n = nrow(aclist[[1]])
thresh = 1e-4

qc_table %<>% mutate(indel = nchar(alleles) > 9)
passvars = qc_table %>% filter(qccum.step4 == 'true', !indel, variant_qc.AF <= 0.05) %$% unique(paste0(locus, alleles))
excall = which(!paste0(aclist[[1]]$locus, aclist[[1]]$alleles) %in% passvars)
# exclude variants where homvar > homref
excall = union(excall, which(aclist[[1]][,'homref'] < aclist[[1]][,'homvar']))

plotlist = list()
pvallist = list()
lambdagc = c()
pvaldf = data.frame()

count = 0
for(j in 1:(length(use)-1)) {
  for(k in (j+1):length(use)) {
  print(count)
  count = count+1
  xy = which(grepl('chrX|chrY', aclist[[1]]$locus))
  sel = setdiff(1:n, xy)
  vars = c('homref','het','homvar')
  #vars = c('homref','het','homvar', 'missing')
  p1mat = as.matrix(aclist[[use[j]]][sel, vars])
  p2mat = as.matrix(aclist[[use[k]]][sel, vars])
  g1mat = cbind(2*p1mat[,1]+p1mat[,2], 2*p1mat[,3]+p1mat[,2])
  g2mat = cbind(2*p2mat[,1]+p2mat[,2], 2*p2mat[,3]+p2mat[,2])
  
  #exc = which((rowSums(p1mat) == 0) | (rowSums(p2mat) == 0) | (rowSums(p1mat > 0) < 2 & rowSums(p2mat > 0) < 2))
  exc = c()
  
  exc = union(exc, excall)
  
  varindices = setdiff(1:nrow(p1mat), exc)
  pvals = rep(NA, nrow(p1mat))
  #pvals[varindices] = sapply(varindices, function(i) tryCatch({chisq.test(cbind(g1mat[i,], g2mat[i,]))$p.value}, warning=function(w) NA))
  pvals[varindices] = sapply(varindices, function(i) tryCatch({poisson.test(c(min(g1mat[i,1], g1mat[i,2]), min(g2mat[i,1], g2mat[i,2])), c(sum(g1mat[i,]), sum(g2mat[i,])))$p.value}, warning=function(w) NA))
  pvals[pvals > 1-1e-10] = NA
  pvallist[[count]] = pvals
  
  if(length(na.omit(pvals)) < 10) {
    plotlist[[count]] = ggplot() + geom_point()
    lambdagc[[count]] = NA
    next
  }
  # title=paste(nams[j], nams[k], sep='\n')
  plotlist[[count]] = ggqq(pvals, title=length(na.omit(pvals)), labs=F)
  lambdagc[[count]] = get.lambda(pvals)
  
  for(l in which(pvals < thresh)) {
    # print(l)
    # print(pvals[l])
    # print(p1mat[l,])
    # print(p2mat[l,])
    pvaldf = rbind(pvaldf, data.frame(n1 = nams[j], n2 = nams[k], locus=aclist[[use[j]]][l, 'locus'], alleles=aclist[[use[j]]][l, 'alleles'], pval=pvals[l],  homref1=p1mat[l, 1], het1=p1mat[l, 2], homvar1=p1mat[l, 3], homref2=p2mat[l, 1], het2=p2mat[l, 2], homvar2=p2mat[l, 3], stringsAsFactors = F))
  }
  
}
}

namplots = list()
for(nam in nams){
  namplots[[length(namplots) + 1]] = ggplot() + geom_text(aes(0, 0), label=nam, size=3) + theme(panel.background = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + xlab('') + ylab('')
}

len = length(use)
a=expand.grid(1:len, 1:len); b = ifelse(apply(a, 1, function(x) x[1] <= x[2]), NA, 1); b2 = b; b2[!is.na(b)] = cumsum(b2[!is.na(b)])
laym = matrix(b2, len)
laym[upper.tri(laym)] = t(laym)[upper.tri(t(laym))] + max(laym, na.rm=T)
diag(laym) = (1:nrow(laym)) + max(laym, na.rm=T)
#laym[upper.tri(laym)] = laym[lower.tri(laym)] = NA

plotlist2 = c(plotlist, plotlist, namplots)

png(paste0('plots/',chrom,'/10_evaluate_variant_qc_batch_qqplots/10_evaluate_variant_qc_batch_qq_ac_controlsonly_AFR_rare_SNVs_norel.png'), width=1600, height=1200)
grid.arrange(grobs=plotlist2, layout_matrix=laym, as.table=TRUE)
dev.off()




j=1; k=6

xy = which(grepl('chrX|chrY', aclist[[1]]$locus))
sel = setdiff(1:n, xy)
vars = c('homref','het','homvar')
#vars = c('homref','het','homvar', 'missing')
p1mat = as.matrix(aclist[[use[j]]][sel, vars])
p2mat = as.matrix(aclist[[use[k]]][sel, vars])

exc = which((rowSums(p1mat) == 0) | (rowSums(p2mat) == 0) | (rowSums(p1mat > 0) < 2 & rowSums(p2mat > 0) < 2))

exc = union(exc, excall)

calcind = setdiff(1:nrow(p1mat), exc)
pvals = rep(NA, nrow(p1mat))
pvals[calcind] = sapply(calcind, function(i) chisq.test(cbind(p1mat[i,], p2mat[i,]))$p.value)

topdiff = aclist[[use[j]]][sel,][head(order(pvals), 10),]
diffvars = left_join(topdiff, bg5, by=c('locus', 'alleles'))





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# break down variant QC
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mincallrate = 0.98
hwep = 0.001
hwep = 1e-9

rt = read_table2('~/Downloads/variant_qc_table_onep_subs.txt')
rt = read_table2('~/Downloads/variant_qc_table_onep.txt.gz')


# statistics for each step (stratified by variant class, not cumulative)
a = rt %>%
  transmute(locus, alleles, ac=variant_qc.AC, af=variant_qc.AF, snv=nchar(alleles) == 9, afclass = ifelse(ac <= 1, 'singleton', ifelse(af <= 0.05, 'rare', 'common')), excesshet_pass = !grepl('ExcessHet', filters), AC0_pass = ac>0, QD_pass=info.QD > 4, VQSR_pass = !grepl('VQSR', filters), cr_AFR_pass = callrate.AFR > mincallrate, cr_AMR_pass = callrate.AMR > mincallrate, cr_EAS_pass = callrate.EAS > mincallrate, cr_EUR_pass = callrate.EUR > mincallrate, cr_SAS_pass = callrate.SAS > mincallrate, cr_FIN_pass = callrate.FIN > mincallrate, cr_PUR_pass = callrate.PUR > mincallrate, hwe_AFR_pass = hwe.AFR > hwep, hwe_AMR_pass = hwe.AMR > hwep, hwe_EAS_pass = hwe.EAS > hwep, hwe_EUR_pass = hwe.EUR > hwep, hwe_SAS_pass = hwe.SAS > hwep, hwe_FIN_pass = hwe.FIN > hwep, hwe_PUR_pass = hwe.PUR > hwep) %>%
  group_by(snv, afclass) %>%
  summarize_at(vars(matches('_pass')), mean) %>%
  gather(filter, pass, -snv, -afclass) %>%
  mutate(variant = ifelse(snv, 'SNV', 'indel'), population=ifelse(!grepl('hwe_|cr_', filter), '', gsub('_[a-z]+$', '', gsub('^[a-z]+_', '', filter))), filtertype=gsub('_.+', '', filter), filtertype2 = recode(filtertype, 'AC0'='0. AC=0', 'excesshet'='1. Excess Het', 'VQSR'='1. VQSR', 'QD'='2. QD < 4', 'cr'='3. Call rate',  'hwe'='4. HWE'))

ggplot(a, aes(x=population, y=1-pass, fill=interaction(variant, afclass))) + geom_bar(position='dodge', stat='identity') + facet_grid(interaction(variant, afclass) ~ filtertype2, scale='free_x', space='free_x') + theme(axis.text.x = element_text(angle=90, hjust=0), strip.text.x = element_text(angle = 90, vjust=0), strip.text.y = element_text(angle = 0), legend.title = element_blank(), legend.position = 'none', panel.background = element_blank()) + xlab('') + ylab('fraction of variants excluded')
ggsave(paste0('plots/',chrom,'/05_variant_qc_exclusion_criteria_noncumulative.pdf'), width=5, height=5)


# cumulative statistics (stratified)
a = rt %>%
  transmute(locus, alleles, ac=variant_qc.AC, af=variant_qc.AF, snv=nchar(alleles) == 9, afclass = ifelse(ac <= 1, 'singleton', ifelse(af <= 0.05, 'rare', 'common')), total_pass=TRUE, AC0_pass = ac>0, excesshet_pass = AC0_pass & !grepl('ExcessHet', filters), VQSR_pass = excesshet_pass & !grepl('VQSR', filters), QD_pass=VQSR_pass & info.QD > 4, cr_pass = QD_pass & lowestcallrate > mincallrate, hwe_pass = cr_pass & lowestphwe > hwep) %>%
  group_by(snv, afclass) %>%
  summarize_at(vars(matches('_pass')), function(x) sum(x)/length(x)) %>%
  gather(filter, pass, -snv, -afclass) %>%
  mutate(variant = ifelse(snv, 'SNV', 'indel'), filtertype=gsub('_.+', '', filter), filtertype2 = recode(filtertype, 'total'=' total', 'AC0'='0. AC=0', 'excesshet'='1. Excess Het', 'VQSR'='1. VQSR', 'QD'='2. QD < 4', 'cr'='3. Call rate',  'hwe'='4. HWE'), filtertype3=factor(filtertype2, levels=rev(levels(factor(filtertype2)))))

ggplot(a, aes(x=filtertype3, y=pass, fill=interaction(variant, afclass))) + geom_bar(position='dodge', stat='identity') + facet_wrap(afclass ~ variant, scale='free_x', nrow=3) + theme(axis.text.x = element_text(angle=90, hjust=0), legend.title = element_blank(), panel.background=element_blank(), legend.position = 'none') + xlab('') + ylab('number of variants remaining') + coord_flip() 
ggsave(paste0('plots/',chrom,'/05_variant_qc_exclusion_criteria_cumulative_freex_percentage.pdf'), width=5, height=6)


write.table(select(ungroup(a), afclass, variant, filtertype, pass, filtertype3), 'out/05_variant_qc_exlcusion_criteria_cumulative.tsv', quote=F, row.names=F, col.names=T, sep='\t')



# cumulative statistics (not stratified, display % of variants remaining and % of loci remaining)
a = rt %>%
  transmute(locus, alleles, ac=variant_qc.AC, af=variant_qc.AF, snv=nchar(alleles) == 9, afclass = ifelse(ac <= 1, '  singleton', ifelse(af <= 0.05, ' rare', 'common')), total_pass=TRUE, AC0_pass = ac>0, excesshet_pass = AC0_pass & !grepl('ExcessHet', filters), VQSR_pass = excesshet_pass & !grepl('VQSR', filters), QD_pass=VQSR_pass & info.QD > 4, cr_pass = QD_pass & lowestcallrate > mincallrate, hwe_pass = cr_pass & lowestphwe > hwep) %>%
  group_by(snv, afclass) %>%
  summarize(total=length(unique(locus[total_pass])), AC0=length(unique(locus[AC0_pass])), excesshet=length(unique(locus[excesshet_pass])), VQSR=length(unique(locus[VQSR_pass])), QD=length(unique(locus[QD_pass])), cr=length(unique(locus[cr_pass])), hwe=length(unique(locus[hwe_pass]))) %>%
  gather(filter, pass, -snv, -afclass) %>%
  mutate(variant = ifelse(snv, ' SNV', 'indel'), filtertype=gsub('_.+', '', filter), filtertype2 = recode(filtertype, 'total'=' total', 'AC0'='0. AC=0', 'excesshet'='1. Excess Het', 'VQSR'='1. VQSR', 'QD'='2. QD < 4', 'cr'='3. Call rate',  'hwe'='4. HWE'), filtertype3=factor(filtertype2, levels=rev(levels(factor(filtertype2)))))

ggplot(a, aes(x=filtertype3, y=pass*100, fill=interaction(variant, afclass))) + geom_bar(position='dodge', stat='identity') + facet_wrap(afclass ~ variant, scale='fixed', nrow=3) + theme(axis.text.x = element_text(angle=90, hjust=0), legend.title = element_blank(), panel.background=element_blank(), legend.position = 'none') + xlab('') + ylab('number of loci remaining') + coord_flip() 
ggsave(paste0('plots/',chrom,'/05_variant_qc_exclusion_criteria_cumulative_fixed_loci.pdf'), width=5, height=6)

ggplot(a, aes(x=filtertype3, y=pass*100, fill=interaction(variant, afclass))) + geom_bar(position='dodge', stat='identity') + facet_wrap(variant ~ afclass, scale='fixed', nrow=2) + theme(axis.text.x = element_blank(), legend.title = element_blank(), panel.background=element_blank(), legend.position = 'none') + xlab('') + ylab('number of variants remaining') + coord_flip() 
ggsave(paste0('plots/',chrom,'/05_variant_qc_exclusion_criteria_cumulative_fixed_loci_poster.pdf'), width=5, height=3.5)


# pairwise center cigar plots
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/onep/variant_qc_ac_files/variant_qc_batch_af_*AFR* ~/Downloads/'))

aflist = list()
for(f in list.files(path='~/Downloads/', pattern='variant_qc_batch_af_AFR.+')) {
  aflist[[f]] = read_table2(paste0('~/Downloads/', f))
}

minaf = 0.001

afs = do.call(rbind, aflist)

afs2 = afs %>% spread(coh, variant_qc.AF)
afs2 %<>% filter(Broad > minaf | Baylor > minaf | NYGC > minaf | WashU > minaf)

afs2 %>% filter(qccum.step4 == 'true') %>% ggplot(aes(Baylor, Broad)) + geom_bin2d(bins=150) + scale_fill_viridis(name = "count", trans = "log10")

cohs = c('Baylor', 'Broad', 'NYGC', 'WashU')
#cohs = c('Baylor', 'Broad')
plotlist = list()
count = 0
for(qc in 0:1) {
  for(j in 1:(length(cohs)-1)) {
    for(k in (j+1):length(cohs)) {
      print(count)
      count = count+1
      if(qc) {
        dat = afs2 %>% filter(qccum.step4 == 'true')
      } else {
        dat = afs2
      }
      plotlist[[count]] = dat %>% ggplot(aes(Baylor, Broad)) + geom_bin2d(bins=150) + scale_fill_viridis(name = "count", trans = "log10") + theme_bw() + theme(legend.position = 'none', axis.text = element_blank()) + xlab('') + ylab('')
      #plotlist[[count]] = dat %>% ggplot(aes(Baylor, Broad)) + geom_point(size=.6) + theme_bw() + theme(legend.position = 'none', axis.text = element_blank()) + xlab('') + ylab('')
  
    }
  }
}

namplots = list()
for(nam in cohs){
  namplots[[length(namplots) + 1]] = ggplot() + geom_text(aes(0, 0), label=nam, size=16) + theme(panel.background = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + xlab('') + ylab('')
}

len = length(cohs)
a=expand.grid(1:len, 1:len); b = ifelse(apply(a, 1, function(x) x[1] <= x[2]), NA, 1); b2 = b; b2[!is.na(b)] = cumsum(b2[!is.na(b)])
laym = matrix(b2, len)
laym[upper.tri(laym)] = t(laym)[upper.tri(t(laym))] + max(laym, na.rm=T)
diag(laym) = (1:nrow(laym)) + max(laym, na.rm=T)
#laym[upper.tri(laym)] = laym[lower.tri(laym)] = NA

plotlist2 = c(plotlist, namplots)

png(paste0('plots/',chrom,'/10_evaluate_variant_qc_cigarplots_center_AFR2.png'), width=1600, height=1200)
grid.arrange(grobs=plotlist2, layout_matrix=laym, as.table=TRUE)
dev.off()

p2 = afs2 %>% filter(qccum.step4 == 'true') %>% ggplot(aes(Baylor, Broad)) + geom_bin2d(bins=150) + scale_fill_viridis(name = "count", trans = "log10") + theme_bw() + theme(legend.position = 'none') + xlab('Baylor') + ylab('Broad')
p1 = afs2 %>% ggplot(aes(Baylor, Broad)) + geom_bin2d(bins=150) + scale_fill_viridis(name = "count", trans = "log10") + theme_bw() + theme(legend.position = 'none') + xlab('Baylor') + ylab('Broad')

ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_cigarplots_center_AFR_BaylorBroad_preQC.png'), p1, width=5, height=4)
ggsave(paste0('plots/',chrom,'/10_evaluate_variant_qc_cigarplots_center_AFR_BaylorBroad_postQC.png'), p2, width=5, height=4)


