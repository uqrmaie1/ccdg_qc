import hail as hl
import sys
import timeit

start = timeit.default_timer()


hl.init(log='/hail.log', min_block_size=2048, default_reference='GRCh38')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# input
onekg_b37_file = 'gs://gnomad-public/truth-sets/hail-0.2/1000G_omni2.5.b37.mt'

onekg = hl.read_matrix_table(onekg_b37_file)

print(onekg.count())


# print runtime
stop = timeit.default_timer()
print("runtime: " + str(stop - start) + " seconds")
