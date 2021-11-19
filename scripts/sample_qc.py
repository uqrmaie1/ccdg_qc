import argparse
import logging
from typing import Optional, Union
import hail as hl

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
from gnomad.sample_qc.ancestry import pc_project

from ccdg_qc.resources import *


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("CCDG_sample_qc")
logger.setLevel(logging.INFO)


def get_qc_vds(
    data_type: str = "genomes",
    pre_filter_variants: bool = False,
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
    vds = hl.vds.read_vds(get_vds_path(data_type))
    if (data_type == "exomes") & interval_qc:
        logger.info("Filtering CCDG exomes VDS to high quality intervals...")
        int_ht = hl.read_table(get_ccdg_results_path(data_type=data_type, result="high_qual_intervals"))
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
    )
    sample_qc_ht = compute_stratified_sample_qc(
        vds,
        strata={
            "bi_allelic": bi_allelic_expr(vds.variant_data),
            "multi_allelic": ~bi_allelic_expr(vds.variant_data),
        },
        tmp_ht_prefix=get_ccdg_results_path(data_type=data_type, result="sample_qc")[:-3],
    )

    return sample_qc_ht.repartition(100)


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
    if args.update_interval_qc_ht:
        int_ht = exomes_qc_intervals_ht(overwrite=args.overwrite)

    if args.sample_qc:
        compute_sample_qc(data_type).write(
            get_ccdg_results_path(data_type=data_type, result="sample_qc_all"),
            overwrite=args.overwrite,
        )
    # Modify when hl.vds.impute_sex() is available
    # if args.impute_sex:
    #     compute_sex().write(
    #         get_ccdg_results_path(data_type=data_type, result="sex"),
    #         overwrite=args.overwrite,
    #     )
    # elif args.reannotate_sex:
    #     reannotate_sex(
    #         args.min_cov,
    #         (args.upper_x, (args.lower_xx, args.upper_xx), args.lower_xxx),
    #         ((args.lower_y, args.upper_y), args.lower_yy),
    #     ).write(
    #         get_ccdg_results_path(data_type=data_type, result="sex"),
    #         overwrite=args.overwrite,
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
        related_samples_to_remove = hl.maximal_independent_set(relartedness_ht.i, pairs.j, False).checkpoint(
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

    if args.pc_project:
        ## TODO: Rank samples and hard filter samples
        mt = hl.read_matrix_table(get_pca_variants_path(ld_pruned=True, data=f"ccdg_{data_type}", mt=True))
        related_ht = hl.read_table(get_ccdg_results_path(data_type=data_type, result="related_samples"))

        related_mt = mt.filter_cols(hl.is_defined(related_mt[mt.col_key]), keep=True)
        pca_mt = mt.filter_cols(hl.is_defined(related_mt[mt.col_key]), keep=False)
        pca_loadings = hl.read_table("gs://gnomad-public-requester-pays/release/3.1/pca/gnomad.v3.1.pca_loadings.ht")

        pca_ht = hl.experimental.pc_project(pca_mt, pca_loadings.loadings, pca_loadings.pca_af)
        pca_mt = pca_mt.annotate_cols(scores=pca_ht[pca_mt.col_key].scores)

        related_ht = hl.experimental.pc_project(related_mt, pca_loadings.loadings, pca_loadings.pca_af)
        related_mt = related_mt.annotate_cols(scores=related_ht[related_mt.col_key].scores)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", help="overwrite", action="store_true")
    parser.add_argument(
        "--exomes",
        help="Sample QC on CCDG exomes. Otherwise, default is CCDG genomes.",
        action="store_true",
    )
    parser.add_argument(
        "--update_interval_qc_ht",
        help="Update Hail Table of high quality CCDG exomes intervals",
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
    args = parser.parse_args()
    main(args)
