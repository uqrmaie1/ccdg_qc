library(tidyverse)
library(viridis)

setwd('~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/')

chrom = 'chr20'
chrom = 'onep'
chrom = 'allchr'

femaleMaxF = 0.4
maleMinF = 0.6

femaleMaxY = 1400
maleMinY = 2000


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# copy files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gsutil = '/Users/robert/google-cloud-sdk/bin/gsutil'
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/out/king/',chrom,'/relatedness_gt0025.kin0.gz ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/', chrom, '/sample_sex_fstat.txt ~/Downloads/'))
#system(paste0(gsutil, ' cp gs://ccdg-qc-multi/out/plink_sexcheck/', chrom, '/ldpruned_common2.sexcheck ~/Downloads/'))
#system(paste0(gsutil, ' cp gs://ccdg-qc-multi/out/plink_sexcheck/', chrom, '/ldpruned_common2.sexcheck ~/Downloads/'))
# from later stage for sex check
#system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/', chrom, '/sample_qc_info_postqc.txt ~/Downloads/'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ibd = as.data.frame(read_table2('~/Downloads/relatedness_gt0025.kin0.gz'))
pheno = read.table('~/Dropbox/postdoc/projects/ccdg_qc/data/phenotypes/combined_phenotypes_broad.txt', stringsAsFactors = F, h=T, comment.char = '', sep='\t')
pheno = read.table('~/Dropbox/postdoc/projects/ccdg_qc/data/phenotypes/combined_phenotypes_feb2018.txt', stringsAsFactors = F, h=T, comment.char = '', sep='\t')
pheno = read.table('~/Dropbox/postdoc/projects/ccdg_qc/data/phenotypes/combined_phenotypes_feb2018_plus_taichi.txt', stringsAsFactors = F, h=T, comment.char = '', sep='\t')
pheno[pheno$sex == 'not recorded' & !is.na(pheno$sex), 'sex'] = NA
sex = read.table(paste0('~/Downloads/sample_sex_fstat.txt'), stringsAsFactors = F, h=T)
#sex2 = read.table(paste0('~/Downloads/ldpruned_common2.sexcheck'), stringsAsFactors = F, h=T)
# from later stage for sex check
# sampleqc = read.table(paste0('~/Downloads/sample_qc_info_postqc.txt'), h=T, stringsAsFactors = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# check sex
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#sex2 = sex %>% mutate(IID=ID, 'F'=sexFstat)
pheno2 = pheno %>%
  right_join(rename(sex, 'FSTAT'='sexFstat'), by=c('subject_id'='ID')) %>%
  #right_join(rename(sex2, 'FSTAT'='F'), by=c('subject_id'='IID')) %>%
  #left_join(select(sampleqc, s, ydp, super), by=c('subject_id'='s')) %>%
  mutate(xfemale = ifelse(FSTAT <= femaleMaxF, T, ifelse(FSTAT >= maleMinF, F, NA)),
                  yfemale = ifelse(ydp <= femaleMaxY, T, ifelse(ydp >= maleMinY, F, NA)),
                  reportedfemale = ifelse(is.na(sex), NA, ifelse(sex=='female', T, F)),
                  predictedfemale = ifelse(is.na(xfemale) | is.na(yfemale), NA, ifelse(xfemale & yfemale, T, ifelse(!xfemale & !yfemale, F, NA))))


ggplot(pheno2, aes(sex, FSTAT, col=xfemale)) + geom_jitter(width=.2) + theme(panel.background = element_blank(), text=element_text(size=16), axis.line = element_line()) + xlab('reported sex') + ylab('chrX inbreeding coefficient') #+ facet_wrap(~ super)
ggsave(paste0('plots/',chrom,'/02_ibd_sex_check.pdf'))

ggplot(pheno2, aes(FSTAT, ydp, col=predictedfemale)) + geom_point() + facet_wrap(~ sex) + theme(panel.background = element_rect(fill=NA, color ='black'), text=element_text(size=16), axis.line = element_line()) + xlab('chrX inbreeding coefficient') + ylab('chrY mean depth')
ggsave(paste0('plots/',chrom,'/02_ibd_sex_check_ydepth.pdf'), height=6, width=8)


sex_mismatch = pheno2 %>% filter (!is.na(reportedfemale) & reportedfemale != predictedfemale | xfemale != yfemale) %>% select(subject_id, xfemale, yfemale, reportedfemale, predictedfemale)
sex_mismatch %>% write.table(paste0('~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/',chrom,'/sex_mismatch.txt'), quote=F, col.names=T, row.names=F)

sexannot = pheno2 %>% transmute(id=subject_id,
                         sexcheck = paste0('x', ifelse(xfemale, 'female', 'male'),
                                           '_y', ifelse(yfemale, 'female', 'male'),
                                           '_reported', ifelse(reportedfemale, 'female', 'male')))
write.table(sexannot, paste0('~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/',chrom,'/sex_annot.txt'), quote=F, col.names=T, row.names=F)
system(paste0(gsutil, ' cp ~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/',chrom,'/sex_annot.txt gs://ccdg-qc-multi/qc_measures/', chrom, '/'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# king ibd
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ggplot(ibd, aes(IBS0, Kinship)) + geom_bin2d(bins=200) + scale_fill_viridis(name = "count", trans = "log10") + geom_hline(yintercept=c(.354, .177, .0884), linetype=1:3) + theme(panel.background = element_blank())
ggsave(paste0('plots/',chrom,'/02_ibd_pihat_ibs0_king_',chrom,'.pdf'), height=5, width=5)

relpairs = ibd %>% dplyr::filter(Kinship > .0884)
exc = unique(relpairs[,2])

write.table(relpairs, paste0('out/ibd_greater_0884_pairs_',chrom,'.txt'), quote=F, row.names=F, col.names=T, sep='\t')
write.table(exc, paste0('out/ibd_greater_0884_',chrom,'.txt'), quote=F, row.names=F, col.names=F)

system(paste0(gsutil, ' cp /Users/robert/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/out/ibd_greater_0884_',chrom,'.txt gs://ccdg-qc-multi/out/king/',chrom,'/'))
system(paste0(gsutil, ' cp /Users/robert/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/out/ibd_greater_0884_pairs_',chrom,'.txt gs://ccdg-qc-multi/out/king/',chrom,'/'))
system(paste0(gsutil, ' cp /Users/robert/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/qc_measures/',chrom,'/sex_mismatch.txt gs://ccdg-qc-multi/qc_full_data/qc_measures/',chrom,'/'))



# ibd2 = relpairs %>% left_join(select(pheno, ID1='ID', pop1='FINAL_POP')) %>% left_join(select(pheno, ID2='ID', pop2='FINAL_POP')) %>% mutate(poppair = paste(pmin(pop1, pop2), pmax(pop1, pop2), sep='_'))
# ggplot(ibd2, aes(IBS0, Kinship, col=poppair)) + geom_point()
