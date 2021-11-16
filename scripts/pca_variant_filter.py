import argparse
import logging
from typing import Optional, Union

import hail as hl

from gnomad.resources.grch38.gnomad import public_release
from gnomad.utils.reference_genome import get_reference_genome

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("pca_variant_determination")
logger.setLevel(logging.INFO)


# TODO: Coordinate with Tim to make a hl.vds.filter_intervals that we can use instead
# https://github.com/broadinstitute/gnomad_methods/blob/3536f87e249f0804f6762facce468597a9c441c6/gnomad/utils/filtering.py#L180
def filter_to_autosomes(
    t: Union[hl.Table, hl.MatrixTable, hl.vds.VariantDataset]
) -> Union[hl.Table, hl.MatrixTable, hl.vds.VariantDataset]:
    """
    Filter Table, MatrixTable or VariantDataset to autosome contigs only.

    This assumes that the input MT or VDS variant_data MT contains a field named `locus` of type Locus

    :param t: Input MatrixTable/Table/VariantDataset
    :return: MatrixTable/Table/VariantDataset subset to autosomes
    """
    if isinstance(t, hl.vds.VariantDataset):
        reference = get_reference_genome(t.variant_data.locus)
    else:
        reference = get_reference_genome(t.locus)
    autosomes = hl.parse_locus_interval(
        f"{reference.contigs[0]}-{reference.contigs[21]}", reference_genome=reference
    )

    if isinstance(t, hl.vds.VariantDataset):
        return hl.vds.filter_intervals(t, [autosomes])
    else:
        return hl.filter_intervals(t, [autosomes])


# TODO: Might need to rethink how we define "high quality" CCDG exomes
def exomes_qc_intervals_ht(
    broad_prop_lower: float = 0.9,
    base_prop_lower: float = 0.5,
    overwrite: bool = False,
    read_if_exist: bool = True,
) -> hl.Table:
    """
    Filter to high quality intervals for exomes QC.

    :param broad_prop_lower: Lower bound for proportion of broad samples defined per interval
    :param base_prop_lower: Lower bound for proportion of bases defined per interval per sample
    :param overwrite: Whether to overwrite CCDG exomes interval QC HT
    :param read_if_exist: Whether to read existing CCDG exomes interval QC HT
    :return: high quality CCDG exomes interval table
    """
    ht = get_sample_manifest_ht("exomes")
    int_mt = hl.read_matrix_table(exomes_interval_mt_path)
    # Compute QC metrics
    int_mt = int_mt.annotate_rows(
        end_pos=hl.if_else(
            int_mt.interval.includes_end,
            int_mt.interval.end.position + 1,
            int_mt.interval.end.position,
        ),
        start_pos=hl.if_else(
            int_mt.interval.includes_start,
            int_mt.interval.start.position,
            int_mt.interval.end.position + 1,
        ),
    )
    int_mt = int_mt.annotate_rows(int_len=int_mt.end_pos - int_mt.start_pos)
    int_mt = int_mt.annotate_cols(center=ht[int_mt.col_id]["cohort"])
    int_mt = int_mt.annotate_entries(prop_base=int_mt.n_bases / int_mt.int_len)
    int_mt = int_mt.annotate_cols(
        washu=hl.if_else(
            (hl.is_defined(int_mt.center) & (int_mt.center == "sccs")), True, False
        )
    )
    n_broad = int_mt.aggregate_cols(hl.agg.count_where(~int_mt.washu))
    int_mt = int_mt.annotate_entries(
        broad_defined=(~(int_mt.washu)) & (int_mt.prop_base > base_prop_lower)
    )
    int_mt = int_mt.annotate_rows(
        broad_defined_prop=hl.agg.count_where(int_mt.broad_defined) / n_broad
    )
    # Filter intervals
    int_mt = int_mt.filter_rows(int_mt.broad_defined_prop > broad_prop_lower)
    int_ht = int_mt.rows().checkpoint(
        f'{get_sample_qc_root(data_type="exomes", mt=False)}/ccdg_exomes_high_qual_intervals.ht',
        overwrite=overwrite,
        _read_if_exists=read_if_exist,
    )

    return int_ht


# TODO: NEED TO HANDLE BI ALLELIC FILTER AFTER AF/CALLRATE calc, and take into account gnomad
# TODO: UKBB (high-qual intervals, high-ish callrate filter on exomes)
# TODO: LD prune gnomaD
# TODO: Do we want to filter out lcrs and segdups?
# TODO: How to handle adj filtering?
# TODO: add in option to only do the ld prune after other MAF/callrate filters
def determine_pca_variants(
    autosomes_only: bool = True,
    snv_only: bool = True,
    gnomad_v3_ac_filter: Optional[int] = None,
    high_qual_ccdg_exome_interval_only: bool = False,
    high_qual_ukbb_exome_interval_only: bool = False,
    af_cutoff: float = 0.001,
    callrate_cutoff: float = 0.99,
    ccdg_exome_callrate_cutoff: float = 0.99,  # TBD
    ukbb_exome_callrate_cutoff: float = 0.99,  # TBD
    ld_pruning: bool = True,
    ld_pruning_dataset: str = "ccdg_genomes",
    ld_r2: float = 0.1,
    read_if_exist: bool = True,
    overwrite: bool = True,
) -> hl.Table:
    """
    Determine a diverse set of variants for relatedness/ancestry PCA using CCDG, gnomAD v3, and UK Biobank.

    :param autosomes_only: Whether to filter to variants in autosomes
    :param snv_only: Whether to filter to SNVs
    :param gnomad_v3_ac_filter: Optional lower bound of AC for variants in gnomAD v3 genomes
    :param high_qual_ccdg_exome_interval_only: Whether to filter to high quality intervals in CCDG exomes
    :param high_qual_ukbb_exome_interval_only: Whether to filter to high quality intervals in UKBB 455K exomes
    :param af_cutoff: lower bound for combined MAF computed from CCDG and gnomAD v3 genomes
    :param callrate_cutoff: lowr bound for combined callrate computed from CCDG and gnomAD v3 genomes
    :param ccdg_exome_callrate_cutoff: lower bound for CCDG exomes callrate
    :param ukbb_exome_callrate_cutoff: lower bound for UKBB exomes callrate
    :param ld_pruning: Whether to conduct LD pruning
    :param ld_pruning_dataset: Which dataset is used for LD pruning, 'ccdg_genomes' or 'gnomAD_genomes'
    :param ld_r2: LD pruning cutoff
    :param read_if_exist: Whether to read existing variant HT
    :param overwrite: Whether to overwrite variant HT
    :return: HT with desired variats for PCA
    """
    logger.info("Loading gnomAD v3 release HT...")
    gnomad_ht = public_release("genomes").ht()
    gnomad_ht = gnomad_ht.select(
        gnomad_AC=gnomad_ht.freq[0].AC, gnomad_AN=gnomad_ht.freq[0].AN,
    )

    logger.info("Getting CCDG genome and exome sample counts...")
    ccdg_genome_count = hl.vds.read_vds(
        get_ccdg_vds_path("genomes")
    ).variant_data.count_cols()
    ccdg_exome_count = hl.vds.read_vds(
        get_ccdg_vds_path("exomes")
    ).variant_data.count_cols()

    def _pre_densify_filter(data_type):
        """
        Wrapper function to get ccdg vds with desired filtering.

        :return: ccdg vds with chosen filters
        """
        logger.info("Loading CCDG %s VDS for initial filtering steps...", data_type)
        vds = hl.vds.read_vds(get_ccdg_vds_path(data_type))
        if autosomes_only:
            logger.info("Filtering CCDG %s VDS to autosomes...", data_type)
            vds = filter_to_autosomes(
                vds
            )  # TODO: switch to hl.vds.filter_to_autosomes if/when available

        variant_filter_expr = True
        ht = vds.variant_data.rows()
        if snv_only:
            logger.info("Filtering CCDG %s VDS to SNVs...", data_type)
            variant_filter_expr &= hl.is_snp(ht.alleles[0], ht.alleles[1])

        if gnomad_v3_ac_filter:
            logger.info(
                "Filtering CCDG %s VDS to gnomAD v3.1.2 variants with adj filtered AC > %d...",
                data_type,
                gnomad_v3_ac_filter,
            )
            variant_filter_expr &= gnomad_ht[ht.key].gnomad_AC > gnomad_v3_ac_filter

        vds = hl.vds.filter_variants(vds, ht.filter(variant_filter_expr), keep=True)

        if high_qual_ccdg_exome_interval_only:
            logger.info(
                "Filtering CCDG %s VDS to high quality CCDG exome intervals...",
                data_type,
            )
            interval_qc_ht = (
                exomes_qc_intervals_ht()
            )  # TODO: Create a checkpointed HT list of "good" intervals
            vds = hl.vds.filter_intervals(
                vds, intervals=interval_qc_ht.interval.collect(), keep=True
            )

        if high_qual_ukbb_exome_interval_only:
            if not autosomes_only:
                raise ValueError(
                    "UK Biobank interval QC filtering is only available for autosomes!"
                )

            logger.info(
                "Filtering CCDG %s VDS to high quality (>85% of samples with 20X coverage) UK Biobank exome intervals...",
                data_type,
            )
            interval_qc_ht = hl.read_table(
                interval_qc_path("broad", 7, "autosomes")
            )  # Note: freeze 7 is all included in gnomAD v4
            interval_qc_ht = interval_qc_ht.filter(
                interval_qc_ht[pct_samples_20x] > 0.85
            )
            vds = hl.vds.filter_intervals(
                vds, intervals=interval_qc_ht.interval.collect(), keep=True
            )

        logger.info("Densifying filtered CCDG %s VDS...", data_type)
        mt = hl.vds.to_dense_mt(vds)
        mt = mt.annotate_rows(
            **{
                f"ccdg_{data_type}_AC": hl.agg.sum(mt.LGT.n_alt_alleles()),
                f"ccdg_{data_type}_AN": hl.agg.count_where(hl.is_defined(mt.LGT)) * 2,
            }
        )
        ht = mt.rows().checkpoint(
            f"{get_sample_qc_root(data_type=data_type, mt=False)}/variant_ccdg_{data_type}_af_callrate.ht",
            overwrite=overwrite,
            _read_if_exists=read_if_exist,
        )

        return ht

    logger.info(
        "Creating Table with joint gnomAD v3.1.2 and CCDG genome allele frequencies and callrate...",
        data_type,
    )
    ccdg_exomes_ht = _pre_densify_filter("exomes")
    ccdg_genomes_ht = _pre_densify_filter("genomes")
    ht = ccdg_exomes_ht.join(ccdg_genomes_ht, how="inner")
    ht = ht.annotate(**gnomad_ht[mt.row_key])
    ht = ht.annotate(
        joint_AC=ht.ccdg_genomes_AC + ht.gnomad_AC,
        joint_AN=ht.ccdg_genomes_AN + ht.gnomad_AN,
    )
    total_genome_an = (gnomad_ht.freq_sample_count + ccdg_genome_count) * 2
    ht = ht.annotate(
        joint_AF=ht.joint_AC / ht.joint_AN, joint_callrate=ht.joint_AN / total_genome_an
    )

    logger.info(
        "Filtering variants to combined gnomAD v3.1.2 and CCDG genome AF of %d and callrate of %d, and CCDG exome callrate of %d...",
        af_cutoff,
        callrate_cutoff,
        ccdg_exome_callrate_cutoff,
    )
    ht = ht.filter(
        (ht.joint_AF > af_cutoff)
        & (ht.joint_callrate > callrate_cutoff)
        & (ht.ccdg_exome_AN / ccdg_exome_count > ccdg_exome_callrate_cutoff)
        #  TODO: UKBB exome callrate filter
    ).checkpoint(
        get_pca_variants_path_ht(ld_pruned=False),
        overwrite=overwrite,
        _read_if_exists=read_if_exist,
    )

    # TODO: Might need to be split if multi-allelic, add in b-alleleic filter?
    if ld_pruning:
        ht = hl.read_table(get_pca_variants_path_ht(ld_pruned=False))
        if ld_pruning_dataset == "ccdg_genomes":
            vds = hl.vds.read_vds(get_ccdg_vds_path("genomes"))
            vds = hl.vds.filter_variants(vds, ht, keep=True)
            mt = hl.vds.to_dense_mt(vds)
            mt = mt.filter_rows(hl.len(mt.alleles) == 2)
            mt = mt.annotate_entries(GT=hl.experimental.lgt_to_gt(mt.LGT, mt.LA))
            ht = hl.ld_prune(mt.GT, r2=ld_r2)
            ht.checkpoint(
                get_pca_variants_path_ht(data=ld_pruning_dataset, ld_pruned=True),
                overwrite=overwrite,
                _read_if_exists=read_if_exist,
            )
        elif ld_pruning_dataset == "gnomAD_genomes":
            mt = get_gnomad_v3_mt(split=~bi_allelic_only)
        # else:
        # Throw error

    ht = ht.annotate_globals(
        af_cutoff=af_cutoff,
        callrate_cutoff=callrate_cutoff,
        ccdg_exome_callrate_cutoff=ccdg_exome_callrate_cutoff,
        ukbb_exome_callrate_cutoff=ukbb_exome_callrate_cutoff,
        ld_r2=ld_r2,
        ld_pruning_dataset=ld_pruning_dataset,
    )

    return ht.select()


def main(args):
    hl.init(log=f"/variant_filter.log")

    filtered_variants = determine_pca_variants(
        autosome_only=args.autosome_only,
        snv_only=args.snv_only,
        gnomad_v3_ac_filter=args.gnomad_v3_ac_filter,
        high_qual_ccdg_interval_only=args.high_qual_ccdg_interval_only,
        high_qual_ukbb_interval_only=args.high_qual_ukbb_interval_only,
        af_cutoff=args.af_cutoff,
        callrate_cutoff=args.callrate_cutoff,
        ccdg_exome_callrate_cutoff=args.ccdg_exome_callrate_cutoff,
        ukbb_exome_callrate_cutoff=args.ukbb_exome_callrate_cutoff,
        ld_pruning=args.ld_pruning,
        ld_pruning_dataset=args.ld_pruning_dataset,
        ld_r2=args.ld_r2,
        read_if_exist=args.read_if_exist,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--autosome_only", help="Filter to autosomes", action="store_true"
    )
    parser.add_argument("--snv_only", help="Filter to SNVs", action="store_true")
    parser.add_argument(
        "--gnomad_v3_ac_filter",
        type=int,
        help="Filter to variants with AC above this value in gnomAD v3",
        default=None,
    )
    parser.add_argument(
        "--high_qual_ccdg_interval_only",
        help="Filter to high quality interval in CCDG exomes",
        action="store_true",
    )
    parser.add_argument(
        "--high_qual_ukbb_interval_only",
        help="Filter to high quality interval in UKBB exomes",
        action="store_true",
    )
    parser.add_argument(
        "--af_cutoff",
        type=float,
        help="Filter to variants with combined MAF above this value",
        default=0.001,
    )
    parser.add_argument(
        "--callrate_cutoff",
        type=float,
        help="Filter to variants with combined callrate above this value",
        default=0.99,
    )
    parser.add_argument(
        "--ccdg_exome_callrate_cutoff",
        type=float,
        help="Filter to variants with callrate above this value in CCDG exomes",
        default=0.99,
    )
    parser.add_argument(
        "--ukbb_exome_callrate_cutoff",
        type=float,
        help="Filter to variants with callrate above this value in UKBB exomes",
        default=0.99,
    )
    parser.add_argument("--ld_pruning", help="Apply LD pruning", action="store_true")
    parser.add_argument(
        "--ld_pruning_dataset",
        type=str,
        help="Dataset to apply LD pruning with",
        default="ccdg_genomes",
    )
    parser.add_argument("--ld_r2", type=float, help="LD pruning cutoff", default=0.1)
    parser.add_argument("--read_if_exist", help="Read if exist", action="store_true")
    parser.add_argument("--overwrite", help="Overwrite", action="store_true")
    args = parser.parse_args()
    main(args)
