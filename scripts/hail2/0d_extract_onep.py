import hail as hl
import sys
import timeit

start = timeit.default_timer()

#chrom = str(sys.argv[1])

hl.init(log='/hail.log', min_block_size=2048)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# input 
vds_file = 'gs://ccdg-qc-multi/vds/raw/hail2_allchr.vds'
onep_file = 'gs://ccdg-qc-multi/out/onep_variants_table.tsv'

# output
vds_onep_file = 'gs://ccdg-qc-multi/vds/raw/hail2_onep.vds'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create subset
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

vds = hl.read_matrix_table(vds_file)

onep = hl.import_table(onep_file, no_header=True).key_by('f0')
onep2 = onep.transmute(**hl.parse_variant(onep.f0, reference_genome=hl.genetics.GenomeReference.GRCh38())).key_by('locus', 'alleles')
vds = vds.filter_rows(hl.is_defined(onep2[vds.locus, vds.alleles]), keep=True)

vds.write(vds_onep_file, overwrite=True)


# print runtime
stop = timeit.default_timer()
print("runtime: " + str(stop - start) + " seconds")

