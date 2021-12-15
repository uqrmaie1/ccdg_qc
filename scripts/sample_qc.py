import argparse
import logging
from typing import Optional, Union
import hail as hl
import pickle

from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.annotations import bi_allelic_expr
from gnomad.sample_qc.filtering import (
    compute_stratified_sample_qc,
    merge_sample_qc_expr,
)
from gnomad.sample_qc.sex import get_ploidy_cutoffs
from gnomad.sample_qc.relatedness import (
    compute_related_samples_to_drop,
    DUPLICATE_OR_TWINS,
    get_relationship_expr,
    PARENT_CHILD,
    SIBLINGS,
    UNRELATED,
)
from gnomad.sample_qc.ancestry import pc_project, assign_population_pcs

from ccdg_qc.resources import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("CCDG_sample_qc")
logger.setLevel(logging.INFO)


def get_qc_vds(
    data_type: str = "genomes",
    split: bool = False,
    autosome_only: bool = False,
    interval_qc: bool = False,
) -> hl.vds.VariantDataset:
    """
    Wrapper function to get ccdg vds with desired filtering
    :param data_type: Whether data is from genomes or exomes, default is genomes
    :param split: Perform split on VDS, default is False
    :param autosome_only: Whether to filter to variants in autosome, default is False
    :param interval_qc: Whether to filter to high quality intervals for exomes QC, default is False
    :return: ccdg vds with chosen filters
    """
    vds = get_ccdg_vds(data_type)
    if data_type == "exomes" and interval_qc:
        logger.info("Filtering CCDG exomes VDS to high quality intervals...")
        int_ht = hl.read_table(
            get_ccdg_results_path(data_type=data_type, result=f"intervals_{INTERVAL_DP}x")
        )
        int_ht = int_ht.filter(int_ht['to_keep'])
        vds = hl.vds.filter_intervals(
            vds, intervals=int_ht.interval.collect(), keep=True
        )

    if split:
        logger.info("Splitting multi-allelic sites in CCDG %s VDS...", data_type)
        vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    if autosome_only:
        logger.info("Filtering CCDG %s VDS to autosomes...", data_type)
        var_ht = filter_to_autosomes(vds.variant_data).rows()
        vds = hl.vds.filter_variants(vds, var_ht, keep=True)

    return vds


def compute_sample_qc(data_type: str = "genomes") -> hl.Table:
    """
    Perform sample QC on the split VDS table using `compute_stratified_sample_qc`.
    :param data_type: Whether data is from genomes or exomes, default is genomes
    :return: Table containing sample QC metrics
    :rtype: hl.Table
    """
    logger.info("Computing sample QC on CCDG %s VDS", data_type)
    vds = get_qc_vds(
        data_type=data_type,
        autosome_only=True,
        split=True,
        interval_qc=True,
    )
    # Use modified compute_stratified_sample_qc with the vds option
    sample_qc_ht = compute_stratified_sample_qc(
        vds,
        strata={
            "bi_allelic": bi_allelic_expr(vds.variant_data),
            "multi_allelic": ~bi_allelic_expr(vds.variant_data),
        },
        tmp_ht_prefix=get_ccdg_results_path(data_type=data_type, result="sample_qc")[
            :-3
        ],
    )

    return sample_qc_ht.repartition(100)

def compute_sex(data_type: str = "genomes", aaf_threshold: float = 0.001, f_stat_cutoff: float = 0.5) -> hl.Table:
    """
    Perform sample QC on the split VDS table using `compute_stratified_sample_qc`.
    :param data_type: Whether data is from genomes or exomes, default is genomes
    :param aaf_threshold: Minimum alternate allele frequency to be used in f-stat calculations.
    :param f_stat_cutoff: f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX samples are below cutoff and XY are above cutoff.
    :return: Table containing imputed sex metrics
    :rtype: hl.Table
    """
    logger.info("Computing sample QC on CCDG %s VDS", data_type)
    vds = get_qc_vds(
        data_type=data_type,
        autosome_only=False,
        split=True,
        interval_qc=True
    )

    calling_intervals = get_calling_intervals_ht(data_type)
    calling_intervals = calling_intervals.filter(calling_intervals.to_keep)

    sex_ht = annotate_sex(
        vds,
        included_intervals=calling_intervals,
        normalization_contig = "chr20",
        reference_genome = "GRCh38",
        gt_expr = "GT",
        f_stat_cutoff = f_stat_cutoff,
        aaf_threshold = aaf_threshold,
    )

    return sex_ht


def annotate_sex(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    is_sparse: bool = True,
    excluded_intervals: Optional[hl.Table] = None,
    included_intervals: Optional[hl.Table] = None,
    normalization_contig: str = "chr20",
    reference_genome: str = "GRCh38",
    sites_ht: Optional[hl.Table] = None,
    aaf_expr: Optional[str] = None,
    gt_expr: str = "GT",
    f_stat_cutoff: float = 0.5,
    aaf_threshold: float = 0.001,
) -> hl.Table:
    """
    Impute sample sex based on X-chromosome heterozygosity and sex chromosome ploidy.

    Return Table with the following fields:
        - s (str): Sample
        - chr20_mean_dp (float32): Sample's mean coverage over chromosome 20.
        - chrX_mean_dp (float32): Sample's mean coverage over chromosome X.
        - chrY_mean_dp (float32): Sample's mean coverage over chromosome Y.
        - chrX_ploidy (float32): Sample's imputed ploidy over chromosome X.
        - chrY_ploidy (float32): Sample's imputed ploidy over chromosome Y.
        - f_stat (float64): Sample f-stat. Calculated using hl.impute_sex.
        - n_called (int64): Number of variants with a genotype call. Calculated using hl.impute_sex.
        - expected_homs (float64): Expected number of homozygotes. Calculated using hl.impute_sex.
        - observed_homs (int64): Expected number of homozygotes. Calculated using hl.impute_sex.
        - X_karyotype (str): Sample's chromosome X karyotype.
        - Y_karyotype (str): Sample's chromosome Y karyotype.
        - sex_karyotype (str): Sample's sex karyotype.

    :param mtds: Input MatrixTable or VariantDataset
    :param bool is_sparse: Whether input MatrixTable is in sparse data format
    :param excluded_intervals: Optional table of intervals to exclude from the computation.
    :param included_intervals: Optional table of intervals to use in the computation. REQUIRED for exomes.
    :param normalization_contig: Which chromosome to use to normalize sex chromosome coverage. Used in determining sex chromosome ploidies.
    :param reference_genome: Reference genome used for constructing interval list. Default: 'GRCh38'
    :param sites_ht: Optional Table to use. If present, filters input MatrixTable to sites in this Table prior to imputing sex,
                    and pulls alternate allele frequency from this Table.
    :param aaf_expr: Optional. Name of field in input MatrixTable with alternate allele frequency.
    :param gt_expr: Name of entry field storing the genotype. Default: 'GT'
    :param f_stat_cutoff: f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX samples are below cutoff and XY are above cutoff.
    :param float aaf_threshold: Minimum alternate allele frequency to be used in f-stat calculations.
    :return: Table of samples and their imputed sex karyotypes.
    """
    logger.info("Imputing sex chromosome ploidies...")

    is_vds = isinstance(mtds, hl.vds.VariantDataset)
    if is_vds:
        if excluded_intervals is not None:
            raise NotImplementedError(
                "excluded_intervals is not used when imputing sex chromosome ploidy for VDS"
            )
        ploidy_ht = hl.vds.impute_sex_chromosome_ploidy(
            mtds,
            calling_intervals=included_intervals,
            normalization_contig=normalization_contig,
        )
        ploidy_ht = ploidy_ht.rename(
            {"x_ploidy": "chrX_ploidy", "y_ploidy": "chrY_ploidy"}
        )
        mt = mtds.variant_data
    else:
        mt = mtds
        if is_sparse:
            ploidy_ht = impute_sex_ploidy(
                mt, excluded_intervals, included_intervals, normalization_contig
            )
        else:
            raise NotImplementedError(
                "Imputing sex ploidy does not exist yet for dense data."
            )

    x_contigs = get_reference_genome(mt.locus).x_contigs
    logger.info("Filtering mt to biallelic SNPs in X contigs: %s", x_contigs)
    if "was_split" in list(mt.row):
        mt = mt.filter_rows((~mt.was_split) & hl.is_snp(mt.alleles[0], mt.alleles[1]))
    else:
        mt = mt.filter_rows(
            (hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1])
        )

    mt = hl.filter_intervals(
        mt,
        [
            hl.parse_locus_interval(contig, reference_genome=reference_genome)
            for contig in x_contigs
        ],
        keep=True,
    )

    if sites_ht is not None:
        if aaf_expr == None:
            logger.warning(
                "sites_ht was provided, but aaf_expr is missing. Assuming name of field with alternate allele frequency is 'AF'."
            )
            aaf_expr = "AF"
        logger.info("Filtering to provided sites")
        mt = mt.annotate_rows(**sites_ht[mt.row_key])
        mt = mt.filter_rows(hl.is_defined(mt[aaf_expr]))

    logger.info("Calculating inbreeding coefficient on chrX")
    sex_ht = hl.impute_sex(
        mt[gt_expr],
        aaf_threshold=aaf_threshold,
        male_threshold=f_stat_cutoff,
        female_threshold=f_stat_cutoff,
        aaf=aaf_expr,
    )

    logger.info("Annotating sex ht with sex chromosome ploidies")
    sex_ht = sex_ht.annotate(**ploidy_ht[sex_ht.key])

    logger.info("Inferring sex karyotypes")
    x_ploidy_cutoffs, y_ploidy_cutoffs = get_ploidy_cutoffs(sex_ht, f_stat_cutoff)
    sex_ht = sex_ht.annotate_globals(
        x_ploidy_cutoffs=hl.struct(
            upper_cutoff_X=x_ploidy_cutoffs[0],
            lower_cutoff_XX=x_ploidy_cutoffs[1][0],
            upper_cutoff_XX=x_ploidy_cutoffs[1][1],
            lower_cutoff_XXX=x_ploidy_cutoffs[2],
        ),
        y_ploidy_cutoffs=hl.struct(
            lower_cutoff_Y=y_ploidy_cutoffs[0][0],
            upper_cutoff_Y=y_ploidy_cutoffs[0][1],
            lower_cutoff_YY=y_ploidy_cutoffs[1],
        ),
        f_stat_cutoff=f_stat_cutoff,
    )
    return sex_ht.annotate(
        **get_sex_expr(
            sex_ht.chrX_ploidy, sex_ht.chrY_ploidy, x_ploidy_cutoffs, y_ploidy_cutoffs
        )
    )


def compute_relatedness(
    data_type: str = "genomes",
    overwrite: bool = False,
) -> hl.Table:
    """
    Perform sample QC on the split VDS table using `compute_stratified_sample_qc`.
    :param data_type: Whether data is from genomes or exomes, default is genomes
    :param overwrite: Whether to overwrite the file
    :return: Table table after running pc_relate
    :rtype: hl.Table
    """
    logger.info("Computing relatedness table on CCDG %s VDS", data_type)
    pca_var_ht = hl.read_table(get_pca_variants_path())
    mt = hl.vds.to_dense_mt(get_qc_vds(data_type))
    mt = mt.filter_rows(hl.is_defined(pca_var_ht[mt.row_key]))
    eig, scores, _ = hl.hwe_normalized_pca(mt.GT, k=10, compute_loadings=False)
    scores = scores.checkpoint(
        get_ccdg_results_path(data_type=data_type, result="pc_scores"),
        overwrite=overwrite,
        _read_if_exists=not overwrite,
    )
    relatedness_ht = hl.pc_relate(
        mt.GT,
        min_individual_maf=0.01,
        scores_expr=scores[mt.col_key].scores,
        block_size=4096,
        min_kinship=0.05,
        statistics="all",
    )
    return relatedness_ht.checkpoint(
        get_ccdg_results_path(data_type=data_type, result="relatedness"),
        overwrite=overwrite,
        _read_if_exists=(not overwrite),
    )


def annotate_relatedness(
    relatedness_ht: hl.Table,
    first_degree_kin_thresholds: float = (0.1767767, 0.4),
    second_degree_kin_cutoff: float = 0.1,
    ibd0_0_max: float = 0.05,
) -> hl.Table:
    relatedness_ht = relatedness_ht.annotate(
        relationship=get_relationship_expr(
            kin_expr=relatedness_ht.kin,
            ibd0_expr=relatedness_ht.ibd0,
            ibd1_expr=relatedness_ht.ibd1,
            ibd2_expr=relatedness_ht.ibd2,
            first_degree_kin_thresholds=tuple(first_degree_kin_thresholds),
            second_degree_min_kin=second_degree_kin_cutoff,
            ibd0_0_max=ibd0_0_max,
        )
    )
    relatedness_ht = relatedness_ht.annotate_globals(
        min_individual_maf=0.01,
        min_emission_kinship=0.05,
        ibd0_0_max=ibd0_0_max,
        second_degree_kin_cutoff=second_degree_kin_cutoff,
        first_degree_kin_thresholds=tuple(first_degree_kin_thresholds),
    )
    return relatedness_ht


def main(args):
    data_type = "exomes" if args.exomes else "genomes"
    hl.init(log=f"/ccdg_sample_qc_{data_type}.log")
    # gcloud compute scp wlu-m:/hard_filter_genomes.log .
    if args.sample_qc:
        compute_sample_qc(data_type).write(
            get_ccdg_results_path(data_type=data_type, result="sample_qc_all"),
            overwrite=args.overwrite,
        )

    if args.impute_sex:
        compute_sex(data_type).write(
            get_ccdg_results_path(data_type=data_type, result="sex"),
            overwrite=args.overwrite,
        )
    # elif args.reannotate_sex:
    #     reannotate_sex(
    #         args.min_cov,
    #         (args.upper_x, (args.lower_xx, args.upper_xx), args.lower_xxx),
    #         ((args.lower_y, args.upper_y), args.lower_yy),
    #     ).write(
    #         get_ccdg_results_path(data_type=data_type, result="sex"),
    #         overwrite=args.overwrite,
    #     )
    ##### Wait for more information
    # if args.compute_hard_filters:
    #     compute_hard_filters(args.min_cov).write(
    #         hard_filtered_samples.path, overwrite=args.overwrite
    #     )

    if args.run_pc_relate or args.reannotate_relatedness:
        if args.run_pc_relate:
            logger.warning(
                "PC-relate requires SSDs and doesn't work with preemptible workers!"
            )
            relatedness_ht = compute_relatedness(
                data_type,
                overwrite=args.overwrite,
            )
        else:
            relatedness_ht = hl.read_table(
                get_ccdg_results_path(data_type=data_type, result="relatedness")
            ).checkpoint(
                "gs://ccdg/tmp/relatedness_ht_checkpoint.ht", overwrite=True
            )  # Copy HT to temp location to overwrite annotation
        relatedness_ht = annotate_relatedness(
            relatedness_ht,
            first_degree_kin_thresholds=tuple(args.first_degree_kin_thresholds),
            second_degree_min_kin=args.second_degree_kin_cutoff,
            ibd0_0_max=args.ibd0_0_max,
        )
        relatedness_ht.write(
            get_ccdg_results_path(data_type=data_type, result="relatedness"),
            overwrite=args.overwrite,
        )

    if args.compute_related_samples_to_drop:
        relatedness_ht = hl.read_table(
            get_ccdg_results_path(data_type=data_type, result="relatedness")
        )
        related_samples_to_remove = hl.maximal_independent_set(
            relartedness_ht.i, pairs.j, False
        ).checkpoint(
            get_ccdg_results_path(data_type=data_type, result="related_samples"),
            overwrite=args.overwrite,
        )

    if args.update_variant_filtered_pca_mt:
        pca_var_ht = hl.read_table(get_pca_variants_path())
        mt = hl.vds.to_dense_mt(get_qc_vds(data_type, split=True))
        mt = mt.filter_rows(hl.is_defined(pca_var_ht[mt.row_key])).checkpoint(
            get_pca_variants_path(ld_pruned=True, data=f"ccdg_{data_type}", mt=True),
            overwrite=args.overwrite,
            _read_if_exists=(not args.overwrite),
        )

    if args.run_pc_project:
        ## TODO: Rank samples and hard filter samples
        mt = hl.read_matrix_table(
            get_pca_variants_path(ld_pruned=True, data=f"ccdg_{data_type}", mt=True)
        )

        pca_loadings = hl.read_table(path_to_gnomad_loadings)

        pca_ht = hl.experimental.pc_project(
            mt.GT,
            pca_loadings.loadings,
            pca_loadings.pca_af,
        )

        pca_ht.checkpoint(
            get_ccdg_results_path(
                data_type=data_type, result="gnomad_pc_project_scores"
            ),
            overwrite=args.overwrite,
        )

        # related_ht = hl.read_table(
        #     get_ccdg_results_path(data_type=data_type, result="related_samples")
        # )
        #
        # related_mt = mt.filter_cols(hl.is_defined(related_mt[mt.col_key]), keep=True)
        # pca_mt = mt.filter_cols(hl.is_defined(related_mt[mt.col_key]), keep=False)

        # pca_ht = hl.experimental.pc_project(
        #     pca_mt.GT, pca_loadings.loadings, pca_loadings.pca_af
        # )
        # pca_mt = pca_mt.annotate_cols(scores=pca_ht[pca_mt.col_key].scores)
        #
        # related_ht = hl.experimental.pc_project(
        #     related_mt.GT, pca_loadings.loadings, pca_loadings.pca_af
        # )
        # related_mt = related_mt.annotate_cols(
        #     scores=related_ht[related_mt.col_key].scores
        # )

    if args.assign_pops:
        with hl.hadoop_open(
            path_to_gnomad_rf,
            "rb",
        ) as f:
            fit = pickle.load(f)

        # Reduce the scores to only those used in the RF model, this was 6 for v2 and 16 for v3.1
        n_pcs = fit.n_features_
        pca_ht = hl.read_table(
            get_ccdg_results_path(
                data_type=data_type, result="gnomad_pc_project_scores"
            )
        )
        pca_ht = pca_ht.annotate(scores=pca_ht.scores[:n_pcs])
        pop_ht, rf_model = assign_population_pcs(
            pca_ht,
            pc_cols=pca_ht.scores,
            fit=fit,
        )

        pop_ht = pop_ht.checkpoint(
            get_ccdg_results_path(data_type=data_type, result="pop_assignment"),
            overwrite=args.overwrite,
            _read_if_exists=not args.overwrite,
        )
        pop_ht.transmute(
            **{f"PC{i + 1}": pop_ht.pca_scores[i] for i in range(n_pcs)}
        ).export(
            get_ccdg_results_path(data_type=data_type, result="pop_assignment")[:-2]
            + "tsv"
        )

        with hl.hadoop_open(
            get_ccdg_results_path(data_type=data_type, result="pop_RF_fit")[:-2]
            + "pickle",
            "wb",
        ) as out:
            pickle.dump(rf_model, out)

    if args.calculate_inbreeding:
        qc_mt = hl.read_matrix_table(
            get_pca_variants_path(ld_pruned=True, data=f"ccdg_{data_type}", mt=True)
        )
        pop_ht = hl.read_table(
            get_ccdg_results_path(data_type=data_type, result="pop_assignment"),
        )
        qc_mt = qc_mt.annotate_cols(pop=pop_ht[qc_mt.col_key].pop)
        qc_mt = qc_mt.annotate_rows(
            call_stats_by_pop=hl.agg.group_by(
                qc_mt.pop, hl.agg.call_stats(qc_mt.GT, qc_mt.alleles)
            )
        )
        inbreeding_ht = (
            qc_mt.annotate_cols(
                inbreeding=hl.agg.inbreeding(
                    qc_mt.GT, qc_mt.call_stats_by_pop[qc_mt.pop].AF[1]
                )
            )
            .cols()
            .select("inbreeding")
        )
        inbreeding_ht.write(
            get_ccdg_results_path(data_type=data_type, result="inbreeding"),
            overwrite=args.overwrite,
        )

    if args.apply_stratified_filters or args.apply_regressed_filters:
        filtering_qc_metrics = args.filtering_qc_metrics.split(",")
        sample_qc_ht = hl.read_table(
            get_ccdg_results_path(data_type=data_type, result="sample_qc_bi_allelic")
        )
        pc_scores = hl.read_table(
            get_ccdg_results_path(data_type=data_type, result="pc_scores")
        )
        sample_qc_ht = sample_qc_ht.select(
            scores=pc_scores[sample_qc_ht.key]["scores"],
        )
        pop_ht = hl.read_table(
            get_ccdg_results_path(data_type=data_type, result="pop_assignment"),
        )

        if "inbreeding" in filtering_qc_metrics:
            inbreeding_ht = hl.read_table(
                get_ccdg_results_path(data_type=data_type, result="inbreeding")
            )[sample_qc_ht.key]
            sample_qc_ht = sample_qc_ht.annotate(
                inbreeding=inbreeding_ht.inbreeding.f_stat
            )

        if args.apply_regressed_filters:
            n_pcs = args.regress_n_pcs
            residuals_ht = compute_qc_metrics_residuals(
                ht=sample_qc_ht,
                pc_scores=sample_qc_ht.scores[:n_pcs],
                qc_metrics={
                    metric: sample_qc_ht[metric] for metric in filtering_qc_metrics
                },
            )
            residuals_ht = residuals_ht.filter(
                hl.is_missing(hard_filtered_samples.ht()[residuals_ht.key])
            )
            stratified_metrics_ht = compute_stratified_metrics_filter(
                ht=residuals_ht,
                qc_metrics=dict(residuals_ht.row_value),
                metric_threshold={
                    "n_singleton_residual": (math.inf, 8.0),
                    "r_het_hom_var_residual": (math.inf, 4.0),
                },
            )

            residuals_ht = residuals_ht.annotate(
                **stratified_metrics_ht[residuals_ht.key]
            )
            residuals_ht = residuals_ht.annotate_globals(
                **stratified_metrics_ht.index_globals(),
                n_pcs=n_pcs,
            )
        else:
            logger.info(
                "Computing stratified QC metrics filters using metrics: "
                + ", ".join(filtering_qc_metrics)
            )
            sample_qc_ht = sample_qc_ht.annotate(qc_pop=pop_ht[sample_qc_ht.key].pop)
            # TODO: compute hard-filtered samples
            sample_qc_ht = sample_qc_ht.filter(
                hl.is_missing(hard_filtered_samples.ht()[sample_qc_ht.key])
            )
            stratified_metrics_ht = compute_stratified_metrics_filter(
                sample_qc_ht,
                qc_metrics={
                    metric: sample_qc_ht[metric] for metric in filtering_qc_metrics
                },
                strata={"qc_pop": sample_qc_ht.qc_pop},
                metric_threshold={"n_singleton": (4.0, 8.0)},
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", help="overwrite", action="store_true")
    parser.add_argument(
        "--exomes",
        help="Sample QC on CCDG exomes. Otherwise, default is CCDG genomes.",
        action="store_true",
    )
    parser.add_argument(
        "--sample_qc", help="Assigns pops from PCA", action="store_true"
    )
    parser.add_argument(
        "--impute_sex",
        help="Runs sex imputation. Also runs sex karyotyping annotation.",
        action="store_true",
    )
    parser.add_argument(
        "--reannotate_sex",
        help="Runs the sex karyotyping annotations again, without re-computing sex imputation metrics.",
        action="store_true",
    )
    parser.add_argument(
        "--compute_hard_filters",
        help="Computes samples to be hard-filtered",
        action="store_true",
    )
    parser.add_argument("--run_pc_relate", help="Run PC-relate", action="store_true")
    parser.add_argument(
        "--reannotate_relatedness",
        help="Runs the relatedness annotation without re-running pc-relate",
        action="store_true",
    )
    parser.add_argument(
        "--first_degree_kin_thresholds",
        help="First degree kinship threshold for filtering a pair of samples with a first degree relationship. \
        Default = (0.1767767, 0.4); \
        Defaults taken from Bycroft et al. (2018)",
        nargs=2,
        default=(0.1767767, 0.4),
        type=float,
    )
    parser.add_argument(
        "--second_degree_kin_cutoff",
        help="Minimum kinship threshold for filtering a pair of samples with a second degree relationship\
        in PC relate and filtering related individuals. (Default = 0.1) \
        Bycroft et al. (2018) calculates 0.08838835 but from evaluation of the distributions v3 has used 0.1",
        default=0.1,
        type=float,
    )
    parser.add_argument(
        "--ibd0_0_max",
        help="IBD0 cutoff to determine parent offspring vs full sibling (Default = 0.05) \
        Default is adjusted from theoretical values; parent-offspring should have an IBD0 = 0. \
        Full siblings should have an IBD0 = 0.25.",
        default=0.05,
    )
    parser.add_argument(
        "--compute_related_samples_to_drop",
        help="Flags related samples to drop",
        action="store_true",
    )
    parser.add_argument(
        "--update_variant_filtered_pca_mt",
        help="Update densified pca variant filtered MatrixTable",
        action="store_true",
    )
    parser.add_argument(
        "--run_pc_project",
        help="Project CCDG data onto gnomAD loadings",
        action="store_true",
    )
    parser.add_argument(
        "--assign_pops", help="Assigns pops from PCA", action="store_true"
    )
    parser.add_argument(
        "--calculate_inbreeding",
        help="Calculate sample level inbreeding",
        action="store_true",
    )
    parser.add_argument(
        "--filtering_qc_metrics",
        help="List of QC metrics for filtering.",
        default=",".join(
            [
                "n_snp",
                "n_singleton",
                "r_ti_tv",
                "r_insertion_deletion",
                "n_insertion",
                "n_deletion",
                "r_het_hom_var",
                "n_transition",
                "n_transversion",
            ]
        ),
    )  # used in v3 'n_het', 'n_hom_var',
    parser.add_argument(
        "--apply_stratified_filters",
        help="Compute per pop filtering.",
        action="store_true",
    )
    parser.add_argument(
        "--apply_regressed_filters",
        help="Computes qc_metrics adjusted for pop.",
        action="store_true",
    )
    parser.add_argument(
        "--regress_n_pcs",
        help="Number of PCs to use for qc metric regressions",
        default=10,
        type=int,
    )
    args = parser.parse_args()
    main(args)
