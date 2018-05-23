import hail as hl
import sys
import timeit

start = timeit.default_timer()

chrom = str(sys.argv[1])

hl.init(log='/hail.log', min_block_size=2048)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# input
qced_vds_file = 'gs://ccdg-qc-multi/vds/qced/' + chrom + '/post_qc2.vds'
lcr_file = 'gs://ccdg-qc-multi/data/LCR-hs38.bed.gz'
samples_to_keep_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/04_sample_qc_2_keep.txt'
king_pedigree = 'gs://ccdg-qc-multi/out/king/' + chrom + '/ldpruned_common2updateparents.fam'
manifest_file = 'gs://ccdg-qc-multi/data/CCDG_Freeze_1_Subsetting_Manifest_2017-12-06_fixedheader.txt'
predicted_pop_file = 'gs://ccdg-qc-multi/qc_measures/pca/' + chrom + '/predicted_populations.tsv'
combined_pheno_file = 'gs://ccdg-qc-multi/data/combined_phenotypes_feb2018_plus_taichi.txt'
rel_exclusion_file = 'gs://ccdg-qc-multi/out/king/' + chrom + '/ibd_greater_0884_' + chrom + '.txt'

# output
variant_qc_stats_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/variant_qc_numbers.txt'
mendel_all1_file = 'gs://ccdg-qc-multi/qc_measures/mendel/' + chrom + '/all1.txt'
mendel_all2_file = 'gs://ccdg-qc-multi/qc_measures/mendel/' + chrom + '/all2.txt'
mendel_perfam1_file = 'gs://ccdg-qc-multi/qc_measures/mendel/' + chrom + '/perfam1.txt'
mendel_perfam2_file = 'gs://ccdg-qc-multi/qc_measures/mendel/' + chrom + '/perfam2.txt'
mendel_persample1_file = 'gs://ccdg-qc-multi/qc_measures/mendel/' + chrom + '/persample1.txt'
mendel_persample2_file = 'gs://ccdg-qc-multi/qc_measures/mendel/' + chrom + '/persample2.txt'
mendel_pervariant1_file = 'gs://ccdg-qc-multi/qc_measures/mendel/' + chrom + '/pervariant1.txt'
mendel_pervariant2_file = 'gs://ccdg-qc-multi/qc_measures/mendel/' + chrom + '/pervariant2.txt'
tdt_pre_table_file = 'gs://ccdg-qc-multi/qc_measures/tdt/' + chrom + '/tdt_pre_table.txt'
tdt_post_table_file = 'gs://ccdg-qc-multi/qc_measures/tdt/' + chrom + '/tdt_post_table.txt'
variants_table_pre_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/variants_table_pre.tsv.gz'
variants_table_post_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/variants_table_post.tsv.gz'

variant_qc_batch_gwas_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/variant_qc_batch_gwas.txt.gz'
variant_qc_batch_gwas_sample_file_prefix = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/variant_qc_batch_sample_files/variant_qc_batch_gwas_'
variant_qc_ac_files_prefix = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/variant_qc_ac_files/variant_qc_batch_ac_'
variant_qc_af_files_prefix = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/variant_qc_ac_files/variant_qc_batch_af_'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ped = hl.Pedigree.read(king_pedigree)

vds_post = hl.read_matrix_table(qced_vds_file)
#vds_post = hl.sample_qc(vds_post.filter_rows(vds_post.qcpass))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filter samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

samples_to_keep = hl.import_table(samples_to_keep_file, no_header=True).key_by('f0')

vds_post = vds_post.filter_cols(hl.is_defined(samples_to_keep[vds_post.s]), keep=True)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate mendel errors
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# all1, perfam1, persample1, pervariant1 = hl.methods.mendel_errors(vds_pre, ped)
# all2, perfam2, persample2, pervariant2 = hl.methods.mendel_errors(vds_post, ped)

# all1.export(mendel_all1_file)
# all2.export(mendel_all2_file)
# perfam1.export(mendel_perfam1_file)
# perfam2.export(mendel_perfam2_file)
# persample1.export(mendel_persample1_file)
# persample2.export(mendel_persample2_file)
# pervariant1.export(mendel_pervariant1_file)
# pervariant2.export(mendel_pervariant2_file)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate TDT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# tdt_pre = hl.methods.transmission_disequilibrium_test(vds_pre, ped)
# tdt_post = hl.methods.transmission_disequilibrium_test(vds_post, ped)
# tdt_pre.export(tdt_pre_table_file)
# tdt_post.export(tdt_post_table_file)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sex / cohort GWAS to test for batch effects
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# annotate with cohort name

table = hl.import_table(manifest_file, impute=True).key_by('Sample_ID')
vds_post = vds_post.annotate_cols(**table[vds_post.s])

table = hl.import_table(predicted_pop_file, impute=True).key_by('Sample')
vds_post = vds_post.annotate_cols(**table[vds_post.s])

table = hl.import_table(combined_pheno_file, impute=True).key_by('subject_id')
vds_post = vds_post.annotate_cols(casecont = table[vds_post.s].casecont)





# custompops = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'FIN', 'PUR']
# cohorts = set(vds_post.Filename_Prefix.collect())
# print(cohorts)
# valid_combinations = []

#vds_post = vds_post.filter_cols((vds_post.casecont == 'control'))

# vds_post = hl.variant_qc(vds_post)
# vds_post = vds_post.filter_rows(vds_post.variant_qc.AC > 10)
# #vds_post = vds_post.sample_rows(0.1)

# #for p in custompops:
# for p in ['EAS', 'EUR', 'SAS', 'FIN']:
# 	print('GWAS for population ' + p + '...')
# 	for coh in cohorts:
# 		cohsub = vds_post.filter_cols((vds_post.Filename_Prefix == coh) & (vds_post.Population == p))
# 		if cohsub.count_cols() > 10:
# 			valid_combinations = valid_combinations + [(p, coh)]
# 			vds_post = vds_post.annotate_cols(batch = hl.cond(vds_post.Population != p, hl.or_missing(False, 1), hl.cond(vds_post.Filename_Prefix == coh, 1, 0)))
# 			vds_post.cols().select('s', 'Filename_Prefix', 'Population', 'batch').export(variant_qc_batch_gwas_sample_file_prefix + p + '_' + coh + '.tsv')
# 			# write out assignment of samples to batch
# 			vds_post = hl.linear_regression(vds_post, [vds_post.batch], vds_post.GT.n_alt_alleles(), covariates=[vds_post.pc1, vds_post.pc2, vds_post.pc3, vds_post.pc4, vds_post.pc5, vds_post.pc6, vds_post.pc7, vds_post.pc8, vds_post.pc9, vds_post.pc10], root='linreg_' + p + '_' + coh)

# vds_post.rows().flatten().export(variant_qc_batch_gwas_file)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# count genotype calls per variant, population and cohort
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# rel_exclusion = hl.import_table(rel_exclusion_file, no_header=True).key_by('f0')
# vds_post = vds_post.filter_cols(hl.is_defined(rel_exclusion[vds_post.s]), keep=False)

# cohorts = set(vds_post.Filename_Prefix.collect())
# custompops = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'FIN', 'PUR']

# for p in ['AFR']:
# #for p in custompops:
# 	print('count alleles for pop ' + p + '...')
# 	for coh in cohorts:
# 		cohsub = vds_post.filter_cols((vds_post.Filename_Prefix == coh) & (vds_post.Population == p))
# 		if cohsub.count_cols() > 20:
# 			print('count alleles for cohort ' + coh + '...')
# 			rt = cohsub.annotate_rows(homref=hl.agg.count_where(cohsub.GT.is_hom_ref()),
# 				                       het=hl.agg.count_where(cohsub.GT.is_het()),
# 				                       homvar=hl.agg.count_where(cohsub.GT.is_hom_var()),
# 				                       missing=hl.agg.count_where(~hl.is_defined(cohsub.GT))).rows().select('locus', 'alleles', 'homref', 'het', 'homvar', 'missing').annotate(pop=p, coh=coh)
# 			rt.export(variant_qc_ac_files_prefix + p + '_' + coh + '.tsv.gz')
# 			#rt.export(variant_qc_ac_files_prefix + p + '_' + coh + '.tsv.gz')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate AF by center
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cohorts = set(vds_post.Center.collect())
custompops = ['AFR', 'EUR']

#for p in ['AFR']:
for p in custompops:
	print('count alleles for pop ' + p + '...')
	for coh in cohorts:
		cohsub = vds_post.filter_cols((vds_post.Center == coh) & (vds_post.Population == p))
		if cohsub.count_cols() > 20:
			print('count alleles for cohort ' + coh + '...')
			cohsub = hl.variant_qc(cohsub)
			rt = cohsub.rows().flatten().select('locus', 'alleles', 'variant_qc.AF', 'qccum.step4').annotate(pop=p, coh=coh)
			rt.export(variant_qc_af_files_prefix + p + '_' + coh + '.tsv.gz')







# print runtime
stop = timeit.default_timer()
print("runtime: " + str(stop - start) + " seconds")



