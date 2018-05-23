import hail as hl
import sys
import timeit

start = timeit.default_timer()


hl.init(log='/hail.log', min_block_size=2048, default_reference='GRCh38')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# input
mhc_chr8inv_file = 'gs://ccdg-qc-multi/data/MHC_invchr8_longLDreg_liftover_to_GRCh38.txt'
vds_1kg_file = 'gs://ccdg-qc-multi/data/1000genomes/vds/hail2_ALL.GRCh38.genotypes.20170504.vds'
onekg_annot_file = 'gs://ccdg-qc-multi/data/1000genomes/onekgpops.tsv'

# output
pca_value_prefix = 'gs://ccdg-qc-multi/1kg/pca_values_'
pca_score_prefix = 'gs://ccdg-qc-multi/1kg/pca_scores_'
pca_loadings_prefix = 'gs://ccdg-qc-multi/1kg/pca_loadings_'

onekg_plink_prefix = 'gs://ccdg-qc-multi/1kg/onekg_ldpruned'
onekg_ldpruned_file = 'gs://ccdg-qc-multi/1kg/onekg_ldpruned.vds'





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## interval list
mhc_chr8inv = hl.import_locus_intervals(mhc_chr8inv_file)

onekg = hl.read_matrix_table(vds_1kg_file)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# annotate samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


table = hl.import_table(onekg_annot_file, impute=True).key_by('Sample')
onekg = onekg.annotate_cols(**table[onekg.s])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filter
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

onekg = onekg.filter_rows(hl.is_defined(mhc_chr8inv[onekg.locus]), keep=False)
onekg = onekg.filter_rows((onekg.locus.contig == "chrX") | (onekg.locus.contig == "chrY"), keep=False)

# filter 5% AF
onekg = hl.variant_qc(onekg)
onekg = onekg.filter_rows(onekg.variant_qc.AF > 0.05, keep=True)

# unphase
onekg2 = onekg.annotate_entries(GT=hl.call(onekg.GT[0], onekg.GT[1], phased=False))
print(onekg2.GT.phased.show())

onekg = onekg2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ld prune
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

onekg = hl.ld_prune(onekg, n_cores=800, r2=0.2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write vds
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

onekg.write(onekg_ldpruned_file, overwrite=True)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write plink
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('export plink')
hl.export_plink(onekg, onekg_plink_prefix, fam_id=onekg.s, id=onekg.s)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pca
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# print('PCA...')
# eigenvalues, scores, loadings = hl.hwe_normalized_pca(onekg, k=10)

# with hl.utils.hadoop_open(pca_value_prefix + '_all.txt', 'w') as f:
# 	for val in eigenvalues:
# 		f.write(str(val) + '\n')
# scores.flatten().export(pca_score_prefix + '_all.txt')


# onekgeur = onekg.filter_cols(onekg.super == 'EUR', keep=True)


# print('PCA EUR...')
# eigenvalues, scores, loadings = hl.hwe_normalized_pca(onekgeur, k=10)

# with hl.utils.hadoop_open(pca_value_prefix + '_eur.txt', 'w') as f:
# 	for val in eigenvalues:
# 		f.write(str(val) + '\n')
# scores.flatten().export(pca_score_prefix + '_eur.txt')




# print runtime
stop = timeit.default_timer()
print("runtime: " + str(stop - start) + " seconds")

