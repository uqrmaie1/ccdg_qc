import hail as hl
import sys
import timeit

start = timeit.default_timer()

chrom = str(sys.argv[1])

hl.init(log='/hail.log', min_block_size=2048, default_reference='GRCh38')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# input
vds_ldpruned_common_file = 'gs://ccdg-qc-multi/vds/qced/' + chrom + '/ldpruned_common.vds'
sample_sex_fstat_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/sample_sex_fstat.txt'
mhc_chr8inv_file = 'gs://ccdg-qc-multi/data/MHC_invchr8_longLDreg_liftover_to_GRCh38.txt'
rel_exclusion_file = 'gs://ccdg-qc-multi/out/king/' + chrom + '/ibd_greater_0884_' + chrom + '.txt'
samples_to_keep_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/01_sample_qc_keep.txt'
predicted_pop_file = 'gs://ccdg-qc-multi/qc_measures/pca/' + chrom + '/predicted_populations.tsv'


# output
pca_value_prefix = 'gs://ccdg-qc-multi/qc_measures/pca/' + chrom + '/perpop/pca_values_'
pca_score_prefix = 'gs://ccdg-qc-multi/qc_measures/pca/' + chrom + '/perpop/pca_scores_'
pca_loadings_prefix = 'gs://ccdg-qc-multi/qc_measures/pca/' + chrom + '/perpop/pca_loadings_'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## interval list
mhc_chr8inv = hl.import_locus_intervals(mhc_chr8inv_file)
##
rel_exclusion = hl.import_table(rel_exclusion_file, no_header=True).key_by('f0')

vds = hl.read_matrix_table(vds_ldpruned_common_file)
vds = vds.select_entries(vds.GT)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filter samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

samples_to_keep = hl.import_table(samples_to_keep_file, no_header=True).key_by('f0')
vds = vds.filter_cols(hl.is_defined(samples_to_keep[vds.s]), keep=True)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# annotate samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

table = hl.import_table(predicted_pop_file, impute=True).key_by('Sample')
vds = vds.annotate_cols(**table[vds.s])





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pca
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


vds_pca = vds.filter_cols(hl.is_defined(samples_to_keep[vds.s]), keep=True)
vds_pca = vds_pca.filter_rows(hl.is_defined(mhc_chr8inv[vds_pca.locus]) & (vds_pca.locus.contig != 'chrX') & (vds_pca.locus.contig != 'chrY'), keep=False)



custompops = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'FIN', 'PUR']

for p in custompops:
	vdspop = vds_pca.filter_cols(vds_pca.Population == p, keep=True)
	#vdspop = vds_pca.filter_samples_expr('sa.pcs.Population == "' + p + '"')
	# eigenvalues, scores, loadings = (vdspop.filter_samples_list(vdspop.sample_ids[:num])
	# .pca_of_normalized_genotypes(k=20, compute_loadings=True)
	# )
	eigenvalues, scores, loadings = hl.hwe_normalized_pca(vdspop, k=10)

	# with hadoop_write(pca_value_prefix + p + str(num) + '.tsv') as f:
	# 	for val in eigenvalues:
	# 		print>>f, val
	# scores.flatten().export(pca_score_prefix + p + str(num) + '.tsv')
	# loadings.flatten().export(pca_loadings_prefix + p + str(num) + '.tsv')
	
	with hl.utils.hadoop_open(pca_value_prefix + p + '.tsv', 'w') as f:
		for val in eigenvalues:
			f.write(str(val) + '\n')
	scores.flatten().export(pca_score_prefix + p + '.tsv')







# print runtime
stop = timeit.default_timer()
print("runtime: " + str(stop - start) + " seconds")

