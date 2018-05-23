import hail as hl
import sys
import timeit


hl.init(log='/hail.log', min_block_size=2048, default_reference='GRCh38')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# extract chr20
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# input
vds_splitmulti_file = 'gs://ccdg-qc-multi/vds/splitmulti/allchr/splitmulti.vds'

vds = hl.read_matrix_table(vds_splitmulti_file)

vds = vds.filter_rows(vds.locus.contig == "chr20")

vds.naive_coalesce(5000).write('gs://ccdg-qc-multi/vds/splitmulti/chr20/splitmulti.vds', overwrite=True)

print(vds.count())



