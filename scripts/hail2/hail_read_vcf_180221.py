import hail as hl
import random
import timeit

start = timeit.default_timer()


hl.init(log='/hail.log', min_block_size=1024)
#hc = hail.HailContext(log='/hail.log', default_reference='GRCh38', min_block_size=64)


#input files
filelist = ['gs://fc-c11923ef-b807-4f9f-8386-1e2eebaebb1e/ccdg_preQC/ccdg_' + str(i) + '_preQC.vcf.gz' for i in range(20360+1)]
#random.seed(123)
#filelist = ['gs://fc-c11923ef-b807-4f9f-8386-1e2eebaebb1e/ccdg_preQC/ccdg_' + str(i) + '_preQC.vcf.gz' for i in random.sample(range(20000), 4)]


samples_to_exclude_file = 'gs://ccdg-qc/data/Broad_CCDG_Freeze1_WGS_Samples_to_Retract_Nov2017_FINAL.txt'
rename_file = 'gs://ccdg-qc/data/rename_eight.txt'


lcr_file = 'gs://ccdg-qc-multi/data/LCR-hs38.bed.gz'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# import LCR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lcr = hl.import_interval_list(lcr_file, reference_genome=hl.genetics.GenomeReference.GRCh38())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# import VCF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("Importing VDS...")
#vds = hc.import_vcf(filelist, header_file='gs://ccdg-qc/data/new_header.0b50850d3629e5cb94eec1df404a55ee.vcf', force_bgz=True)
vds = hl.import_vcf(filelist, header_file='gs://ccdg-qc-multi/data/new_header.0b50850d3629e5cb94eec1df404a55ee_GTADDPGQPLPGTPID.vcf', force_bgz=True, reference_genome=hl.genetics.GenomeReference.GRCh38())
print("Import done!")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# exclude samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

samples_to_exclude = hl.import_table(samples_to_exclude_file)
samples_to_exclude = samples_to_exclude.rename({'Collaborator Sample ID' : 'ID'})
samples_to_exclude_list = [row.ID for row in samples_to_exclude.collect()]

#vds = vds.annotate_globals(samples_to_exclude = set(samples_to_exclude_list))
#vds = vds.filter_cols(vds.samples_to_exclude.contains(vds.s)) 

vds = vds.filter_cols(hl.is_defined(samples_to_exclude.key_by('ID')[vds.s]), keep=False)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rename samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mapping_table = hl.import_table(rename_file)
mapping_dict = {row.old_id: row.new_id for row in mapping_table.collect()}
#vds = vds.rename(mapping_dict)

mapping = hl.broadcast(mapping_dict)
vds = vds.annotate_cols(s = hl.bind(mapping.get(vds.s), lambda s: hl.cond(hl.is_missing(s), vds.s, s)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filter LCR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

vds = vds.filter_rows(hl.is_defined(lcr[vds.locus]), keep=False)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write VDS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

vds.write('gs://ccdg-qc-multi/vds/raw/hail2_allchr.vds', overwrite=True)
#vds.write('gs://ccdg-qc-multi/vds/raw/hail2_10GB_subset.vds', overwrite=True)


#hl.export_vcf(vds, 'gs://ccdg-qc-multi/vcf/raw/hail2_10GB_subset.vcf.bgz', parallel='separate_header')

# print runtime

stop = timeit.default_timer()

print("runtime: " + str(stop - start) + " seconds")



# start 15:49
