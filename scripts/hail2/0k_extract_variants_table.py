import hail as hl
import sys
import timeit


hl.init(log='/hail.log', min_block_size=2048)



vds = hl.read_matrix_table('gs://ccdg-qc-multi/vds/qced/onep/post_qc2.vds')


print(vds.rows().fields)
vds.rows().select('locus', 'alleles', 'qccum', 'qcpass').flatten().export('gs://ccdg-qc-multi/qc_measures/onep/hail2_variants_table_post.tsv.gz')




