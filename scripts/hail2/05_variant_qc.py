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
vds_splitmulti_file = 'gs://ccdg-qc-multi/vds/splitmulti/' + chrom + '/splitmulti.vds'
lcr_file = 'gs://ccdg-qc-multi/data/LCR-hs38.bed.gz'
samples_for_variantqc_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/04_sample_qc_2_keep.txt'
samples_to_keep_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/01_sample_qc_keep.txt'
predicted_pop_file = 'gs://ccdg-qc-multi/qc_measures/pca/' + chrom + '/predicted_populations.tsv'
predicted_sex_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/sample_sex_fstat.txt'
sex_mismatch_file = 'gs://ccdg-qc-multi/qc_full_data/qc_measures/' + chrom + '/sex_mismatch.txt'


# output
qced_vds_file = 'gs://ccdg-qc-multi/vds/qced/' + chrom + '/post_qc2.vds'
qced_vds_stringent_file = 'gs://ccdg-qc-multi/vds/qced/' + chrom + '/post_stringent_qc.vds'
variant_qc_stats_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/variant_qc_numbers.txt'
variant_qc_table_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/variant_qc_table.txt.gz'
sample_qc_info_postqc_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/sample_qc_info_postqc.txt'
sample_qc_info_poststringentqc_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/sample_qc_info_poststringentqc.txt'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define constants
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mincallrate = 0.98
hwep = 0.000000001

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# vds_pre = hl.read_matrix_table(vds_splitmulti_file)

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # filter samples
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# samples_to_keep = hl.import_table(samples_to_keep_file, no_header=True).key_by('f0')
# vds_pre = vds_pre.filter_cols(hl.is_defined(samples_to_keep[vds_pre.s]), keep=True)


# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # annotate samples
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# table = hl.import_table(predicted_pop_file, impute=True).key_by('Sample')
# vds_pre = vds_pre.annotate_cols(**table[vds_pre.s])

# table = hl.import_table(predicted_sex_file, impute=True).key_by('ID')
# vds_pre = vds_pre.annotate_cols(**table[vds_pre.s])



# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #
# # variant qc
# #
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# samples_for_variantqc = hl.import_table(samples_for_variantqc_file, no_header=True).key_by('f0')
# vds = vds_pre.filter_cols(hl.is_defined(samples_for_variantqc[vds_pre.s]), keep=True)
# vds = hl.variant_qc(vds)


# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # step3.1: set missing:
# #    DP > 400
# #    Heterozygous AND [ (AD alt/DP) < 20% OR PL ref < 20 OR (AD ref + AD alt) / DP < 90% ]
# #    Homozygous ref AND [ GQ < 20 OR (AD ref / DP) < 0.9 ]
# #    Homozygous alt AND [ PL ref < 20 OR (AD alt / DP) < 0.9 ]
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# print("starting step 3.1...")

# #vds = vds.annotate_rows(dp10 = hl.agg.count_where(vds.DP < 10), dp400 =  hl.agg.count_where(vds.DP > 400))
# #vds = vds.annotate_cols(ydp = hl.agg.count_where((vds.locus.contig == 'chrY') & (vds.DP > 0)))

# #vds = hl.variant_qc(vds)
# vds = (vds.filter_entries(
#      (vds.DP > 400) |
#      (vds.DP < 10) |
#      (vds.GT.is_hom_ref() & ((vds.AD[0] / vds.DP < 0.9) | (vds.GQ < 20))) |
#      (vds.GT.is_hom_var() & ((vds.AD[1] / vds.DP < 0.9) | (vds.PL[0] < 20) | ((vds.variant_qc.AF > 0.5) & (vds.GQ < 20)))) |
#      (vds.GT.is_het() & ( ((vds.AD[0] + vds.AD[1]) / vds.DP < 0.9) |
#           (vds.AD[1] / vds.DP < 0.20) |
#           (vds.PL[0] < 20) |
#           (~ vds.isFemale & ( vds.locus.in_x_nonpar() |
#                vds.locus.in_y_nonpar() )))) |
#      (vds.locus.in_y_nonpar() & vds.isFemale)
#      , keep=False))
# vds = hl.variant_qc(vds)

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # step3.2: population specific call rate > .98
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# print("starting step 3.2...")

# custompops = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'FIN', 'PUR']

# # pops with at least 500 samples
# vds = vds.annotate_globals(mypops = custompops)
# pops = custompops + ['other']

# vds = vds.annotate_cols(pp2 = hl.cond(vds.mypops.contains(vds.Population), vds.Population, 'other'))

# callratelist = [hl.agg.fraction(hl.agg.filter(vds.pp2 == pop, hl.is_defined(vds.GT))) for pop in custompops]

# vds = vds.annotate_rows(
#   callrate=hl.Struct(**dict(zip(pops[:-1], callratelist))),
#   lowestcallrate = hl.min(hl.array(callratelist))
# )
# #vds = vds.annotate_rows(lowestcallrate = hl.min(vds.callrate))



# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # step4: population specific HWE P-value > 1e-09
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# print("starting step 4...")


# hwelist = [hl.agg.hardy_weinberg(hl.agg.filter(vds.pp2 == pop, vds.GT))['p_hwe'] for pop in custompops]
# vds = vds.annotate_rows(
#   hwe=hl.Struct(**dict(zip(pops[:-1], hwelist))),
#   lowestphwe = hl.min(hl.array(hwelist))
# )
# #vds = vds.annotate_rows(lowestphwe = hl.min(vds.hwe))


# step0 = vds.variant_qc.AC > 0
# step1 = vds.filters.size() == 0
# step2 = vds.info.QD > 4
# step3 = vds.lowestcallrate >= mincallrate
# step4 = vds.lowestphwe >= hwep


# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # annotate vds with filters
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# qc_struct = hl.Struct(step0 = step0, step1 = step1, step2 = step2, step3=step3, step4=step4)
# #vds = vds.annotate_rows(qc = qc_struct)

# qccumul = hl.Struct(step0 = step0,
#                     step1 = step0 & step1,
#                     step2 = step0 & step1 & step2,
#                     step3 = step0 & step1 & step2 & step3,
#                     step4 = step0 & step1 & step2 & step3 & step4)

# vds = vds.annotate_rows(qc = qc_struct, qccum = qccumul, qcpass = qccumul.step4)


# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # annotate vds_pre
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# print("copying annotations...")
# vds_post = vds_pre.annotate_rows(**vds[(vds_pre.locus, vds_pre.alleles), :])

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # downcode
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# print("downcoding...")
# vds_post = vds_post.select_entries(vds_post.GT)



# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # write output VDS
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# print("writing output files...")
# vds_post.write(qced_vds_file, overwrite=True)

vds_post = hl.read_matrix_table(qced_vds_file)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# perform sample qc on remaining variants
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# print("sample QC...")
# vds_sampleqc = hl.sample_qc(vds_post.filter_rows(vds_post.qccum.step4))
# vds_sampleqc.cols().flatten().export(sample_qc_info_postqc_file)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write variants table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("write variants table...")
vds_post.rows().select('locus', 'alleles', 'info', 'callrate', 'lowestcallrate', 'hwe', 'lowestphwe', 'qc', 'qccum', 'qcpass', 'variant_qc', 'filters').flatten().export(variant_qc_table_file)



# print runtime

stop = timeit.default_timer()
print("runtime: " + str(stop - start) + " seconds")



