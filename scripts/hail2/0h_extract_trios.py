import hail as hl
import sys
import timeit

start = timeit.default_timer()

hl.init(log='/hail.log', min_block_size=2048)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# input
rawvds_file = 'gs://ccdg-qc-multi/out/full_set.vds'
rawvds_file = 'gs://ccdg-qc-multi/vds/raw/hail2_allchr.vds'

samples_to_exclude_file = 'gs://ccdg-qc/data/Broad_CCDG_Freeze1_WGS_Samples_to_Retract_Nov2017_FINAL.txt'
trioid_file = 'gs://ccdg-qc-multi/out/king/trioids.txt'

# output
vds_trio_file = 'gs://ccdg-qc-multi/vds/trios/splitmulti.vds'



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vds = hl.read_matrix_table(rawvds_file)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filter samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

samples_to_keep = hl.import_table(trioid_file, no_header=True)

vds = vds.filter_cols(hl.is_defined(samples_to_keep.key_by('f0')[vds.s]), keep=True)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filter variants
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#vds = vds.filter_variants_expr('v.nAltAlleles() == 1')
vds = vds.filter_rows(hl.len(vds.alleles) == 2, keep=True)

vds = hl.variant_qc(vds)
vds.filter_rows(vds.variant_qc.AC > 0)

vds.naive_coalesce(1000).write(vds_trio_file, overwrite=True)




# print runtime
stop = timeit.default_timer()
print("runtime: " + str(stop - start) + " seconds")

