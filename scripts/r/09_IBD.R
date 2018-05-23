library(plyr)
library(tidyverse)
library(magrittr)
library(viridis)
library(openxlsx)
library(data.table)

setwd('~/Dropbox/postdoc/projects/ccdg_qc/qc_full_data/')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# process IBD variant files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cdvars = read.xlsx('../data/mark_IBD/CDmeta_122017.xlsx')
ibdvars = read.xlsx('../data/mark_IBD/IBDmeta_122017.xlsx')

locus = gsub(':[A-Z].+', '', ibdvars[,1])
cat(paste0('chr', locus, '-', gsub('.+:', '', locus)), sep='\n')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in liftover BED coordinates
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd38 = read.table('../data/mark_IBD/cdvars.bed')
ibd38 = read.table('../data/mark_IBD/ibdvars.bed')

write.table(gsub('x-.+', '', rbind(ibd38, cd38)[,1]), '~/Dropbox/postdoc/projects/ccdg_qc/data/mark_IBD/ibdvars.txt', quote=F, row.names=F, col.names=F)
gsutil = '/Users/robert/google-cloud-sdk/bin/gsutil'
system(paste0(gsutil, ' cp ~/Dropbox/postdoc/projects/ccdg_qc/data/mark_IBD/ibdvars.txt gs://ccdg-qc-multi/data/IBD/'))

locusmap = rbind(data.frame(grch37=gsub(':[A-Z].+', '', ibdvars[,1]), grch38=gsub('-.+', '', ibd38[,1]), stringsAsFactors = F),
                 data.frame(grch37=gsub(':[A-Z].+', '', cdvars[,1]), grch38=gsub('-.+', '', cd38[,1]), stringsAsFactors = F))
locusmap = locusmap[!duplicated(locusmap),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compare GRCh37 to GRCh38 alleles
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

system(paste0(gsutil, ' cp gs://ccdg-qc-multi/data/IBD/vartable.tsv ~/Downloads/'))
newalleles = read.table('~/Downloads/vartable.tsv', h=T, stringsAsFactors = F)
newalleles %<>% mutate(alleles = gsub(',', ':', gsub('\"|\\[|\\]', '', alleles)))

oldalleles = rbind(cdvars[,1,drop=F], ibdvars[,1,drop=F]) %>% distinct %>% transmute(locus=gsub(':[A-Z].+', '', VARID), oldalleles=gsub('.+[0-9]:', '', VARID))


allelemap = left_join(locusmap, oldalleles, by=c('grch37'='locus')) %>% left_join(newalleles, by=c('grch38'='locus'))

mean(allelemap[,3] == allelemap[,4])

system(paste0(gsutil, ' cp gs://ccdg-qc-multi/data/IBD/samples_file.tsv ~/Downloads/'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample counts
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

system(paste0(gsutil, ' cp gs://ccdg-qc-multi/data/IBD/samples_file.tsv ~/Downloads/'))
stable = read.table('~/Downloads/samples_file.tsv', stringsAsFactors = F, h=T)
table(stable$Population, stable$diagnosis_specific)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merge CCDG logreg results with previous results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

system(paste0(gsutil, ' cp gs://ccdg-qc-multi/data/IBD/logreg_results.tsv ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/data/IBD/logreg_results_AFR.tsv ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/data/IBD/samples_file.tsv ~/Downloads/'))
ass = read.table('~/Downloads/logreg_results.tsv', stringsAsFactors = F, h=T)
ass = read.table('~/Downloads/logreg_results_AFR.tsv', stringsAsFactors = F, h=T)

# new CD columns
ass %>%
  mutate(alleles=gsub(',', ':', gsub('\"|\\[|\\]', '', alleles))) %>%
  left_join(allelemap, by=c('locus'='grch38', 'alleles'='alleles')) %>%
  mutate(VARID = paste(grch37, oldalleles, sep=':')) %>%
  right_join(cdvars, by='VARID') %>%
  mutate(index=cumsum(!duplicated(locus)), origall=alleles==oldalleles) %>%
  arrange(index, as.numeric(!origall)) %>%
  mutate(locrepeat = ave(locus==locus, locus, FUN=cumsum)) %>% #filter(locus=='chr10:73913343') %>% head
  mutate(caseac=cdhomvar*2+cdhet, casean=(cdhomvar+cdhet+cdhomref)*2, conac=2*conthomvar+conthet, conan=(conthomvar+conthet+conthomref)*2, caseaf=caseac/casean, conaf=conac/conan, missing=(contmissing+cdmissing)/(casean+conan+contmissing+cdmissing)) %>%
  transmute(index, locus, alleles, caseaf, conaf, caseac, casean, conac, conan, missing, beta=logreg_cd.beta, sebeta=logreg_cd.standard_error, z=beta/sebeta, pval=logreg_cd.p_value, or=exp(beta), locrepeat) %>%
  setDT %>%
  data.table::dcast(locus ~ locrepeat, value.var = c('index', 'alleles', 'caseaf', 'conaf', 'caseac', 'casean', 'conac', 'conan', 'missing', 'beta', 'sebeta', 'z', 'pval', 'or')) %>%
  arrange(index_1) %>% dplyr::select(-matches('index')) %>%
  dplyr::select(locus, ends_with('_1'), ends_with('_2'), ends_with('_3')) %>%
  write.xlsx('../data/mark_IBD/out/cd_columns.xlsx')
  
# new IBD columns
ass %>%
  mutate(alleles=gsub(',', ':', gsub('\"|\\[|\\]', '', alleles))) %>%
  left_join(allelemap, by=c('locus'='grch38', 'alleles'='alleles')) %>%
  mutate(VARID = paste(grch37, oldalleles, sep=':')) %>%
  right_join(ibdvars, by='VARID') %>%
  mutate(index=cumsum(!duplicated(locus)), origall=alleles==oldalleles) %>%
  arrange(index, as.numeric(!origall)) %>%
  mutate(locrepeat = ave(locus==locus, locus, FUN=cumsum)) %>%
  mutate(caseac=metahomvar*2+metahet, casean=(metahomvar+metahet+metahomref)*2, conac=2*conthomvar+conthet, conan=(conthomvar+conthet+conthomref)*2, caseaf=caseac/casean, conaf=conac/conan, missing=(contmissing+metamissing)/(casean+conan+contmissing+metamissing)) %>%
  transmute(index, locus, alleles, caseaf, conaf, caseac, casean, conac, conan, missing, beta=logreg_meta.beta, sebeta=logreg_meta.standard_error, z=beta/sebeta, pval=logreg_meta.p_value, or=exp(beta), locrepeat) %>%
  setDT %>%
  data.table::dcast(locus ~ locrepeat, value.var = c('index', 'alleles', 'caseaf', 'conaf', 'caseac', 'casean', 'conac', 'conan', 'missing', 'beta', 'sebeta', 'z', 'pval', 'or')) %>%
  arrange(index_1) %>% dplyr::select(-matches('index')) %>%
  dplyr::select(locus, ends_with('_1'), ends_with('_2'), ends_with('_3')) %>%
  write.xlsx('../data/mark_IBD/out/ibd_columns.xlsx')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AFR case control numbers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

a = read.table('~/Downloads/samples_file.tsv', h=T, stringsAsFactors = F, sep='\t'); a %<>% filter(Population == 'AFR', diagnosis_class == 'ibd' | casecont == 'control')
a %>% group_by(study_nickname, diagnosis_specific) %>% summarize(count=n()) %>% spread(diagnosis_specific, count, fill=0) %>% as.data.frame


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read PLINK files, test PRS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(plink2R)

system(paste0(gsutil, ' cp gs://ccdg-qc-multi/qc_measures/pca/allchr/perpop/pca_scores_AFR.tsv ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/data/IBD/ibdvars* ~/Downloads/'))
pl = read_plink('~/Downloads/ibdvars')
pcaafr = read.table('~/Downloads/pca_scores_AFR.tsv', h=T, stringsAsFactors = F)

oldbetas = ass %>%
  mutate(alleles=gsub(',', ':', gsub('\"|\\[|\\]', '', alleles))) %>%
  left_join(allelemap, by=c('locus'='grch38', 'alleles'='alleles')) %>%
  mutate(VARID = paste(grch37, oldalleles, sep=':')) %>%
  right_join(ibdvars, by='VARID') %>%
  transmute(V2=paste(locus, alleles, sep=':'), BETA)


famannot = left_join(pl$fam, stable, by=c('V1'='s'))
bimannot = left_join(pl$bim, oldbetas, by='V2')

excsnps = which(apply(pl$bed, 2, function(x) mean(is.na(x))) > .05 | is.na(bimannot$BETA))
predmat = scale(pl$bed[,-excsnps], scale=F)
predmat[is.na(predmat)] = 0

scores = predmat %*% bimannot$BETA[-excsnps]
famannot$scores = scores[,1]
famannot %<>% left_join(pcaafr, by=c('V1'='s'))
famannot$scorescorrected = resid(lm(scores ~ sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=famannot, na.action=na.exclude))


tt = famannot %>% filter(Population %in% c('other', 'AFR', 'PUR')) %>% #filter(Population == 'AFR') %>%
  group_by(Population) %>% summarize(pval = t.test(scores[diagnosis_specific == 'control'], scores[diagnosis_specific != 'control'])$p.value)

famannot %>% filter(Population %in% c('other', 'AFR', 'PUR')) %>%
  ggplot(aes(diagnosis_specific == 'control', scores)) + geom_violin() + geom_boxplot(width=.4) + facet_wrap(~ Population) + geom_text(data=tt, aes(1.5, 5, label=round(pval, 3)))
ggsave('~/Dropbox/postdoc/projects/ccdg_qc/data/mark_IBD/PRS_violins.pdf')

famannot %>% filter(Population == 'AFR') %>% summarize(pval = t.test(scorescorrected[diagnosis_specific == 'control'], scorescorrected[diagnosis_specific != 'control'])$p.value) %$% pval
# 0.01652786







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Manuel genes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("http://bioconductor.org/workflows.R")
workflowInstall("liftOver")
library(rtracklayer)
library(gwascat)

system(paste0(gsutil, ' cp gs://ccdg-qc-multi/data/IBD/ibdvars* ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/data/IBD/logreg_results_AFR.tsv ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/data/IBD/manuel_genes.txt ~/Downloads/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/data/IBD/logreg_results__vep38.tsv ~/Downloads/'))


vep = read_table2('~/Downloads/logreg_results__vep38.tsv')
ass = read.table('~/Downloads/logreg_results_AFR.tsv', stringsAsFactors = F, h=T)
manuel38 = read.table('~/Downloads/manuel_genes.txt', stringsAsFactors = F, h=F) %>% transmute(chr = gsub(':.+', '', V1), start=gsub('.+:', '', gsub('-.+', '',V1)), end=gsub('.+-', '', V1), gene=c('NOD2', 'IL23R', 'CARD9'))

chainfile = '~/Downloads/hg38ToHg19.over.chain'
ch = import.chain(chainfile)

ass %<>% mutate(chromosome = gsub(':.+', '', locus), start=gsub('.+:', '', locus), end=start)
cur = ass %>% dplyr::select(chromosome, start, end) %>% makeGRangesFromDataFrame
cur19 = as.data.frame(liftOver(cur, ch))
ov = as.data.frame(findOverlaps(makeGRangesFromDataFrame(manuel38), cur, select='all'))

ass %>% mutate(variant_class=vep$vep.variant_class, locus37 = paste(cur19$seqnames, cur19$start, sep=':')) %>% mutate(gene = manuel38[ov[match(1:n(), ov[,2]), 1], 'gene']) %>% filter(!is.na(gene)) %>%
  mutate(alleles=gsub(',', ':', gsub('\"|\\[|\\]', '', alleles))) %>%
  mutate(ibdac=metahomvar*2+metahet, cdac=cdhomvar*2+cdhet, ibdan=(metahomvar+metahet+metahomref)*2, cdan=(cdhomvar+cdhet+cdhomref)*2, conac=2*conthomvar+conthet, conan=(conthomvar+conthet+conthomref)*2, ibdaf=ibdac/ibdan, cdaf=cdac/cdan, conaf=conac/conan, missing=(contmissing+metamissing)/(ibdan+cdan+conan+contmissing+metamissing)) %>%
  transmute(gene, locus38=locus, locus37, alleles, af_ibd=ibdaf, af_cd=cdaf, af_con=conaf, ac_ibd=ibdac, ac_cd=cdac, ac_con=conac, an_ibd=ibdan, an_cd=cdan, an_con=conan, missing, beta_ibd=logreg_meta.beta, beta_cd=logreg_cd.beta, sebeta_ibd=logreg_meta.standard_error, sebeta_cd=logreg_cd.standard_error, z_ibd=beta_ibd/sebeta_ibd, z_cd=beta_cd/sebeta_cd, pval_ibd=logreg_meta.p_value, pval_cd=logreg_cd.p_value, or_ibd=exp(beta_ibd), or_cd=exp(beta_cd), variant_class, inlist = ifelse(paste(gsub('chr', '', locus37), alleles, sep=':') %in% c(ibdvars[,1], cdvars[,1]), 'x', '')) %>%
  filter((af_ibd > 0 | af_con > 0) & missing < .02) %>% mutate(pval_cd = ifelse(ac_cd + ac_con > 2, pval_cd, NA), pval_ibd = ifelse(ac_ibd + ac_con > 2, pval_ibd, NA) ) %>% arrange(pmin(pval_ibd, pval_cd)) %>%
  write.xlsx('../data/mark_IBD/out/IBD_CCDG_Freeze1_NOD2_IL23R_CARD9.xlsx')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# March 2018 format immune ccdg phenotypes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#phen = read.table('~/Dropbox/postdoc/projects/ccdg_qc/data/phenotypes/combined_phenotypes_feb2018_plus_taichi.txt', sep='\t', stringsAsFactors = F, h=T)

immune = readxl::read_xlsx('~/Dropbox/postdoc/projects/ccdg_qc/data/phenotypes/IMMUNE_CCDG_WGS_MANIFEST_MARCH2018.xlsx', sheet=2, range=readxl::cell_cols("A:S"))

immune %<>% transmute(subject_id = gsub('H_VD-', '', SUBJID), sex=recode(SEX, FEMALE='female', MALE='male', MalLE='male', UNKNOWN='unknown'), race_ethnicity=gsub(' ', '_', `RACE/ETHNICITY`), center=SITE, project=COHORT, study_nickname=COHORT, age_baseline=NA, age_disease=NA, height=NA, target_coverage=NA, casecont=ifelse(grepl('CONTROL', ANALYSIS_CAT), 'control', 'case'), diagnosis_class=ifelse(grepl('T1D', ANALYSIS_CAT), 't1d', 'ibd'), diagnosis_specific = ifelse(casecont == 'control', 'control', recode(DIAGNOSIS, CD='cd', CONTROL='control', 'IBD-U'='ibdu', 'IBDU'='ibdu', UC='uc')))

#write.table(dplyr::select(immune, subject_id, sex, race_ethnicity, casecont, diagnosis_class, diagnosis_specific), '../data/mark_IBD/immune_ccdg_wgs_manifest_march2018_formatted.tsv', quote=F, col.names=T, row.names=F, sep='\t')




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# v2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

system(paste0(gsutil, ' cp gs://ccdg-qc-multi/data/IBD/v2/vartable.tsv ~/Dropbox/postdoc/projects/ccdg_qc/data/mark_IBD/v2/'))
newalleles = read.table('~/Dropbox/postdoc/projects/ccdg_qc/data/mark_IBD/v2/vartable.tsv', h=T, stringsAsFactors = F)
newalleles %<>% mutate(alleles = gsub(',', ':', gsub('\"|\\[|\\]', '', alleles)))

oldalleles = rbind(cdvars[,1,drop=F], ibdvars[,1,drop=F]) %>% distinct %>% transmute(locus=gsub(':[A-Z].+', '', VARID), oldalleles=gsub('.+[0-9]:', '', VARID))


allelemap = left_join(locusmap, oldalleles, by=c('grch37'='locus')) %>% left_join(newalleles, by=c('grch38'='locus'))


#system(paste0(gsutil, ' cp gs://ccdg-qc-multi/data/IBD/v2/logreg_results.tsv ~/Dropbox/postdoc/projects/ccdg_qc/data/mark_IBD/v2/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/data/IBD/v2/logreg_results_AFRICAN_AMERICAN.tsv ~/Dropbox/postdoc/projects/ccdg_qc/data/mark_IBD/v2/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/data/IBD/v2/logreg_results_HISPANIC.tsv ~/Dropbox/postdoc/projects/ccdg_qc/data/mark_IBD/v2/'))
system(paste0(gsutil, ' cp gs://ccdg-qc-multi/data/IBD/v2/samples_file.tsv ~/Downloads/'))
#ass = read.table('~/Dropbox/postdoc/projects/ccdg_qc/data/mark_IBD/v2/logreg_results.tsv', stringsAsFactors = F, h=T)

pop = 'AFRICAN_AMERICAN'
pop = 'HISPANIC'
ass = read.table(paste0('~/Dropbox/postdoc/projects/ccdg_qc/data/mark_IBD/v2/logreg_results_', pop, '.tsv'), stringsAsFactors = F, h=T)

# new CD columns
ass %>%
  mutate(alleles=gsub(',', ':', gsub('\"|\\[|\\]', '', alleles))) %>%
  left_join(allelemap, by=c('locus'='grch38', 'alleles'='alleles')) %>%
  mutate(VARID = paste(grch37, oldalleles, sep=':')) %>%
  right_join(cdvars, by='VARID') %>%
  mutate(index=cumsum(!duplicated(locus)), origall=alleles==oldalleles) %>%
  arrange(index, as.numeric(!origall)) %>%
  mutate(locrepeat = ave(locus==locus, locus, FUN=cumsum)) %>% #filter(locus=='chr10:73913343') %>% head
  mutate(caseac=cdhomvar*2+cdhet, casean=(cdhomvar+cdhet+cdhomref)*2, conac=2*conthomvar+conthet, conan=(conthomvar+conthet+conthomref)*2, caseaf=caseac/casean, conaf=conac/conan, missing=(contmissing+cdmissing)/(casean+conan+contmissing+cdmissing)) %>%
  transmute(index, locus, alleles, caseaf, conaf, caseac, casean, conac, conan, missing, beta=logreg_cd.beta, sebeta=logreg_cd.standard_error, z=beta/sebeta, pval=logreg_cd.p_value, or=exp(beta), locrepeat) %>%
  setDT %>%
  filter(locrepeat == 1 | is.na(locrepeat)) %>% select(-locrepeat) %>%
  #data.table::dcast(locus ~ locrepeat, value.var = c('index', 'alleles', 'caseaf', 'conaf', 'caseac', 'casean', 'conac', 'conan', 'missing', 'beta', 'sebeta', 'z', 'pval', 'or')) %>%
  #arrange(index_1) %>% dplyr::select(-matches('index')) %>%
  #dplyr::select(locus, ends_with('_1'), ends_with('_2'), ends_with('_3')) %>%
  write.xlsx(paste0('../data/mark_IBD/v2/out/cd_columns_', pop, '.xlsx'))

# new IBD columns
ass %>%
  mutate(alleles=gsub(',', ':', gsub('\"|\\[|\\]', '', alleles))) %>%
  left_join(allelemap, by=c('locus'='grch38', 'alleles'='alleles')) %>%
  mutate(VARID = paste(grch37, oldalleles, sep=':')) %>%
  right_join(ibdvars, by='VARID') %>%
  mutate(index=cumsum(!duplicated(locus)), origall=alleles==oldalleles) %>%
  arrange(index, as.numeric(!origall)) %>%
  mutate(locrepeat = ave(locus==locus, locus, FUN=cumsum)) %>%
  mutate(caseac=metahomvar*2+metahet, casean=(metahomvar+metahet+metahomref)*2, conac=2*conthomvar+conthet, conan=(conthomvar+conthet+conthomref)*2, caseaf=caseac/casean, conaf=conac/conan, missing=(contmissing+metamissing)/(casean+conan+contmissing+metamissing)) %>%
  transmute(index, locus, alleles, caseaf, conaf, caseac, casean, conac, conan, missing, beta=logreg_meta.beta, sebeta=logreg_meta.standard_error, z=beta/sebeta, pval=logreg_meta.p_value, or=exp(beta), locrepeat) %>% 
  setDT %>%
  filter(locrepeat == 1 | is.na(locrepeat)) %>% select(-locrepeat) %>%
  #data.table::dcast(locus ~ locrepeat, value.var = c('index', 'alleles', 'caseaf', 'conaf', 'caseac', 'casean', 'conac', 'conan', 'missing', 'beta', 'sebeta', 'z', 'pval', 'or')) %>%
  #arrange(index_1) %>% dplyr::select(-matches('index')) %>%
  #dplyr::select(locus, ends_with('_1'), ends_with('_2'), ends_with('_3')) %>%
  write.xlsx(paste0('../data/mark_IBD/v2/out/ibd_columns_', pop, '.xlsx'))





