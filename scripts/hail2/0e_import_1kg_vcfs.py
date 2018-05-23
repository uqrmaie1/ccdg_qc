import hail as hl
import sys
import timeit

start = timeit.default_timer()

hl.init(log='/hail.log', min_block_size=2048)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# input
# if VCF
onekg_vcf_file = 'gs://ccdg-qc-multi/data/1000genomes/ALL.chr*_GRCh38.genotypes.20170504.vcf.gz'

# output
onekg_vds_file = 'gs://ccdg-qc-multi/data/1000genomes/vds/hail2_ALL.GRCh38.genotypes.20170504.vds'


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#filelist = [rawvcf_file + str(i) + '_preQC.vcf.gz' for i in range(20360+1)]
# if VCF
#vds = hc.import_vcf(rawvcf_file + '*bgz', force_bgz=True)
# if VDS
vds = hl.import_vcf(onekg_vcf_file, force_bgz=True, reference_genome=hl.genetics.GenomeReference.GRCh38(), contig_recoding = dict((x, "chr"+x) for x in list(map(str, range(1,23))) + ["X", "Y"]))
print("Import done!")

vds.write(onekg_vds_file, overwrite=True)


# print runtime
stop = timeit.default_timer()
print("runtime: " + str(stop - start) + " seconds")

