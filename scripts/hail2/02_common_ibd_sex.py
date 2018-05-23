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
samples_to_keep_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/01_sample_qc_keep.txt'
par_file = 'gs://ccdg-qc-multi/data/par_chrX_hg38_v2.txt'
#raw_vds_file = 'gs://ccdg-qc-multi/vds/raw/hail2_' + chrom + '.vds'

# output
sample_sex_fstat_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/sample_sex_fstat.txt'
ibd_results_file = 'gs://ccdg-qc-multi/qc_measures/' + chrom + '/ibd_results.txt'
vds_common_file = 'gs://ccdg-qc-multi/vds/qced/' + chrom + '/common.vds'
vds_ldpruned_common_file = 'gs://ccdg-qc-multi/vds/qced/' + chrom + '/ldpruned_common.vds'
vds_ldpruned_common_plink = 'gs://ccdg-qc-multi/vds/qced/' + chrom + '/ldpruned_common'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

par = hl.import_locus_intervals(par_file)
vds = hl.read_matrix_table(vds_splitmulti_file)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filter samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

samples_to_keep = hl.import_table(samples_to_keep_file, no_header=True).key_by('f0')
vds = vds.filter_cols(hl.is_defined(samples_to_keep[vds.s]), keep=True)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# variant QC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print("variant QC...")

vds = vds.filter_entries(((vds.locus.contig != "chrX") & (vds.locus.contig != "chrY")) &
	(((vds.AD[0] + vds.AD[1]) / vds.DP < 0.9) |
	(vds.GT.is_hom_ref() & ((vds.AD[0] / vds.DP < 0.9) | (vds.GQ < 20))) |
	(vds.GT.is_het() & ((vds.AD[1] / vds.DP < 0.20) | (vds.PL[0] < 20))) |
	(vds.GT.is_hom_var() & ((vds.AD[1] / vds.DP < 0.9) | (vds.PL[0] < 20))) |
	(vds.DP > 200)), keep=False)

vds = hl.variant_qc(vds)
vds = vds.select_entries(vds.GT)

vds1 = vds.filter_rows((vds.variant_qc.AF > 0.01) &
						(vds.variant_qc.AF < 0.99) &
						(vds.filters.size() == 0) &
						(~vds.was_split) &
						((vds.locus.contig == "chrX") | (vds.locus.contig == "chrY") |
						((vds.info.QD > 4) &
						(vds.variant_qc.call_rate > 0.99) &
						(vds.variant_qc.dp_mean > 20))), keep=True)





# print("writing common VDS...")
# vds1.write(vds_common1_file, overwrite=True)


#vds1 = hl.read_matrix_table(vds_common1_file)

vds5 = vds1.filter_rows((vds1.variant_qc.AF >= 0.05) & (vds1.variant_qc.AF <= 0.95))
vds5.write(vds_common_file, overwrite=True)

#vds5 = hl.read_matrix_table(vds_common_file)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sex imputation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


print("sex imputation...")
vdsnopar = vds5.filter_rows(hl.is_defined(par[vds5.locus]), keep=False)
vdsnopar = vdsnopar.annotate_cols(ydp = hl.agg.count_where((vdsnopar.locus.contig == 'chrY') & (hl.is_defined(vdsnopar.GT))))

vdsx = vdsnopar.filter_rows((vdsnopar.locus.contig == "chrX") & (vdsnopar.variant_qc.AF >= 0.05) & (vdsnopar.variant_qc.AF <= 0.95))
ct = hl.impute_sex(vdsx.GT, female_threshold=0.6, male_threshold=0.7)
vdsct = vdsnopar.cols()
ct = ct.annotate(ydp = vdsct[ct.s].ydp)

(ct.select(ID=ct.s, sexFstat=ct.f_stat, isFemale=ct.is_female, ydp=ct.ydp)
.export(sample_sex_fstat_file))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ld pruning
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


print("LD pruning...")
vds5_ldp = hl.ld_prune(vds5, n_cores=1600, r2=0.1)
#vds5_ldp = hl.ld_prune(vds5, n_cores=60, r2=0.2, window=1000000, memory_per_core=512)

print("writing LD pruned VDS...")
vds5_ldp.write(vds_ldpruned_common_file, overwrite=True)
hl.export_plink(vds5_ldp, vds_ldpruned_common_plink, fam_id=vds5_ldp.s, id=vds5_ldp.s)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IBD analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# use king until pcrelate works

#vds.ibd(min=0.1).flatten().rename({'ibd.Z0': 'Z0', 'ibd.Z1': 'Z1', 'ibd.Z2': 'Z2', 'ibd.PI_HAT': 'PI_HAT'}).export(ibd_results_file)


# print runtime
stop = timeit.default_timer()
print("runtime: " + str(stop - start) + " seconds")



