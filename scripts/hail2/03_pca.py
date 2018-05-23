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
vds_1kg_file = 'gs://ccdg-qc-multi/data/1000genomes/vds/hail2_ALL.GRCh38.genotypes.20170504.vds'
mhc_chr8inv_file = 'gs://ccdg-qc-multi/data/MHC_invchr8_longLDreg_liftover_to_GRCh38.txt'
rel_exclusion_file = 'gs://ccdg-qc-multi/out/king/' + chrom + '/ibd_greater_0884_' + chrom + '.txt'
samples_to_keep_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/01_sample_qc_keep.txt'


# output
pca_value_file = 'gs://ccdg-qc-multi/qc_measures/pca/' + chrom + '/pca_values.tsv'
pca_score_file = 'gs://ccdg-qc-multi/qc_measures/pca/' + chrom + '/pca_scores.tsv'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## interval list
#mhc_chr8inv = hl.import_table(mhc_chr8inv_file, no_header=True).key_by('f0')
mhc_chr8inv = hl.import_locus_intervals(mhc_chr8inv_file)

##
rel_exclusion = hl.import_table(rel_exclusion_file, no_header=True).key_by('f0')

vds = hl.read_matrix_table(vds_ldpruned_common_file)
onekg = hl.read_matrix_table(vds_1kg_file)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rename samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print('annotate 1KG...')
onekg = onekg.annotate_cols(s = '1KG_' + onekg.s)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filter samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print('filter samples...')
samples_to_keep = hl.import_table(samples_to_keep_file, no_header=True).key_by('f0')
vds = vds.filter_cols(hl.is_defined(samples_to_keep[vds.s]), keep=True)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# joing with 1000 genomes vds
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print('combining data with 1KG...')
vds = onekg.union_cols(vds)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# downsample to 80,000 variants
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#num_rows = 317938 # from onep vds_ldpruned_common01_file bim file
#num_rows = 211924 # from allchr vds_ldpruned_common01_file bim file


print('number of variants in VDS:')
num_rows = vds.count_rows()
print(num_rows)

prob = min(1, 80000 / num_rows)
vds = vds.sample_rows(prob)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pca
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print('filter VDS...')
vds = vds.filter_cols(hl.is_defined(rel_exclusion[vds.s]), keep=False)
vds = vds.filter_rows(hl.is_defined(mhc_chr8inv[vds.locus]), keep=False)
vds = vds.filter_rows((vds.locus.contig == "chrX") | (vds.locus.contig == "chrY"), keep=False)

print('PCA...')
eigenvalues, scores, loadings = hl.hwe_normalized_pca(vds, k=10)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write output
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

with hl.utils.hadoop_open(pca_value_file, 'w') as f:
	for val in eigenvalues:
		f.write(str(val) + '\n')
scores.flatten().export(pca_score_file)

# print runtime
stop = timeit.default_timer()
print("runtime: " + str(stop - start) + " seconds")

