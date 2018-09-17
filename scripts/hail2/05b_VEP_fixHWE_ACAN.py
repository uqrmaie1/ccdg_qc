import hail as hl
import sys
import timeit
import os

start = timeit.default_timer()

chrom = str(sys.argv[1])

hl.init(log='/hail.log', min_block_size=2048, default_reference='GRCh38')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# input
qced_vds_file = 'gs://ccdg-qc-multi/vds/qced/' + chrom + '/post_qc2.vds'
samples_to_keep_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/04_sample_qc_2_keep.txt'

# output
qced_vds_excl_file = 'gs://ccdg-qc-multi/vds/qced/' + chrom + '/post_qc2_excl.vds'
sample_qc_info_postqc_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/5b_sample_qc_info_postqc.txt'
variant_qc_table_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/5b_variant_qc_table.txt.gz'



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define constants
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mincallrate = 0.98
hwep = 0.000001

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

vds = hl.read_matrix_table(qced_vds_file)

samples_to_keep = hl.import_table(samples_to_keep_file, no_header=True).key_by('f0')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fix HWE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

step0 = vds.variant_qc.AC > 0
step1 = vds.filters.size() == 0
step2 = vds.info.QD > 4
step3 = vds.lowestcallrate >= mincallrate
step4 = vds.lowestphwe >= hwep


qc_struct = hl.Struct(step0 = step0, step1 = step1, step2 = step2, step3=step3, step4=step4)
#vds = vds.annotate_rows(qc = qc_struct)

qccumul = hl.Struct(step0 = step0,
                    step1 = step0 & step1,
                    step2 = step0 & step1 & step2,
                    step3 = step0 & step1 & step2 & step3,
                    step4 = step0 & step1 & step2 & step3 & step4)

vds = vds.annotate_rows(qc = qc_struct, qccum = qccumul, qcpass = qccumul.step4)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# remove samples and variants
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("remove samples...")
vds = vds.filter_cols(hl.is_defined(samples_to_keep[vds.s]), keep=True)


print("remove variants...")
vds = vds.filter_rows(vds.qcpass)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# perform sample and variant qc on remaining variants
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("sample QC...")
vds = hl.sample_qc(vds)

print("variant QC...")
vds = hl.variant_qc(vds)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add VEP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

vds = hl.vep(vds, "gs://ccdg-qc-multi/data/vep85-GRCh38-gcloud.json")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write output VDS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("writing VDS...")
vds.write(qced_vds_excl_file, overwrite=True)


print("writing sample QC...")
vds.cols().flatten().export(sample_qc_info_postqc_file)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write variants table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#vds = hl.read_matrix_table(qced_vds_excl_file)

print("write variants table...")
#vds.rows().select('info', 'callrate', 'lowestcallrate', 'hwe', 'lowestphwe', 'qc', 'qccum', 'qcpass', 'variant_qc', 'filters').flatten().export(variant_qc_table_file)
vds.rows().select('info', 'callrate', 'lowestcallrate', 'hwe', 'lowestphwe', 'qc', 'qccum', 'qcpass', 'variant_qc', 'filters', 'vep').flatten().export(variant_qc_table_file)



# print runtime

stop = timeit.default_timer()
print("runtime: " + str(stop - start) + " seconds")













