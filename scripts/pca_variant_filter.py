import argparse
import logging
from typing import Optional, Union

import hail as hl

from gnomad.resources.grch38.gnomad import public_release as gnomad_public_release
from gnomad.utils.sparse_mt import densify_sites, filter_ref_blocks
from gnomad.utils.reference_genome import get_reference_genome
from gnomad_qc.resources.annotations import (
    last_END_position as gnomad_last_END_position,
)
from gnomad_qc.resources.basics import get_gnomad_v3_mt
from ukbb_qc.resources.basics import release_ht_path as ukbb_release_ht_path
from ukbb_qc.resources.sample_qc import interval_qc_path as ukbb_interval_qc_path
from ukbb_qc.resources.sample_qc import meta_ht_path as ukbb_meta_ht_path

from ccdg_qc.resources import (
    get_ccdg_vds_path,
    get_pca_variants_path_ht,
    get_sample_manifest_ht,
)

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
) -> hl.Table:
    """
    Filter to high quality intervals for exomes QC.

    :param broad_prop_lower: Lower bound for proportion of broad samples defined per interval
    :param base_prop_lower: Lower bound for proportion of bases defined per interval per sample
    :param overwrite: Whether to overwrite CCDG exomes interval QC HT
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
        _read_if_exists=(not overwrite),
    )

    return int_ht


# Note: this was previous QC variant filtering method https://github.com/broadinstitute/gnomad_methods/blob/35066ffc01d63ac2d7b20e069ea6703013ae9da7/gnomad/sample_qc/pipeline.py#L110
# TODO: Do we want to filter out lcrs and segdups?
# TODO: How to handle adj filtering?
# TODO: Variants passing hard thresholds? QD >= 2, FS <= 60 and MQ >= 30
# TODO: inbreeding_coeff_threshold and min_hardy_weinberg_threshold (needs to be done on dense)?
# TODO: Might need to think about adding in a few extra options to reuse already created files so we can add in an option to only do the ld prune after other MAF/callrate filters. Is this properly handled by `read_if_exist` or do we need more?
# TODO: Rethink names of functions, files, and parameters
def determine_pca_variants(
    autosomes_only: bool = True,
    snv_only: bool = True,
    bi_allelic_only: bool = True,
    gnomad_v3_ac_filter: Optional[int] = None,
    high_qual_ccdg_exome_interval_only: bool = False,
    high_qual_ukbb_exome_interval_only: bool = False,
    af_cutoff: float = 0.001,  # TODO: Konrad mentioned that he might want to lower this
    callrate_cutoff: float = 0.99,
    ccdg_exome_callrate_cutoff: float = 0.99,  # TODO: What parameter should this start with?
    ukbb_exome_callrate_cutoff: float = 0.99,  # TODO: What parameter should this start with?
    ld_pruning: bool = True,
    ld_pruning_dataset: str = "ccdg_genomes",
    ld_r2: float = 0.1,
    overwrite: bool = True,
) -> None:
    """
    Determine a diverse set of variants for relatedness/ancestry PCA using CCDG, gnomAD v3, and UK Biobank.

    :param autosomes_only: Whether to filter to variants in autosomes
    :param snv_only: Whether to filter to SNVs
    :param bi_allelic_only: Whether to filter variants that are bi-allelic in both CCDG and gnomAD v3
    :param gnomad_v3_ac_filter: Optional lower bound of AC for variants in gnomAD v3 genomes
    :param high_qual_ccdg_exome_interval_only: Whether to filter to high quality intervals in CCDG exomes
    :param high_qual_ukbb_exome_interval_only: Whether to filter to high quality intervals in UKBB 455K exomes
    :param af_cutoff: Lower bound for combined MAF computed from CCDG and gnomAD v3 genomes
    :param callrate_cutoff: Lower bound for combined callrate computed from CCDG and gnomAD v3 genomes
    :param ccdg_exome_callrate_cutoff: Lower bound for CCDG exomes callrate
    :param ukbb_exome_callrate_cutoff: Lower bound for UKBB exomes callrate
    :param ld_pruning: Whether to conduct LD pruning
    :param ld_pruning_dataset: Which dataset is used for LD pruning, 'ccdg_genomes' or 'gnomAD_genomes'
    :param ld_r2: LD pruning cutoff
    :param overwrite: Whether to overwrite variant HT
    :return: Table with desired variants for PCA
    """
    logger.info("Loading gnomAD v3.1.2 release HT and UK Biobank 455K release HT ...")
    gnomad_ht = gnomad_public_release("genomes").ht()
    gnomad_ht = gnomad_ht.select(
        gnomad_was_split=gnomad_ht.was_split,
        gnomad_AC=gnomad_ht.freq[0].AC,
        gnomad_AN=gnomad_ht.freq[0].AN,
    )
    ukbb_ht = hl.read_table(ukbb_release_ht_path("broad", 7))
    ukbb_ht = ukbb_ht.select(
        ukbb_AC=gnomad_ht.freq[0].AC, ukbb_AN=gnomad_ht.freq[0].AN,
    )
    ukbb_meta_ht = hl.read_table(ukbb_meta_ht_path("broad", 7))

    # Only count samples used in the UK Biobank exome frequency calculations
    ukbb_exome_count = ukbb_meta_ht.filter(
        ukbb_meta_ht.sample_filters.high_quality
        & hl.is_defined(mt.meta.ukbb_meta.batch)
        & ~mt.meta.sample_filters.related
    ).count()

    logger.info("Getting CCDG genome and exome sample counts...")
    ccdg_genome_count = hl.vds.read_vds(
        get_ccdg_vds_path("genomes")
    ).variant_data.count_cols()
    ccdg_exome_count = hl.vds.read_vds(
        get_ccdg_vds_path("exomes")
    ).variant_data.count_cols()

    def _initial_filter(data_type):
        """
        Get Table of CCDG variants passing desired filters.

        Possible filters are:
            - Autosomes only
            - SNVs only
            - gnomAD v3.1.2 AC filter
            - CCDG high quality exome intervals
            - UK Biobank high quality exome intervals

        After densification of the VDS, rows are annotated with:
            - ccdg_{data_type}_was_split
            - ccdg_{data_type}_AC
            - ccdg_{data_type}_AN

        The filtered and annotated rows are returned as a Table and are also checkpointed

        :return: Table of CCDG filtered variants
        """
        logger.info(
            "Loading CCDG %s VDS and splitting multi-allelics for initial filtering steps...",
            data_type,
        )
        vds = hl.vds.read_vds(get_ccdg_vds_path(data_type))
        vds = hl.vds.split_mulit(vds)

        if autosomes_only:
            logger.info("Filtering CCDG %s VDS to autosomes...", data_type)
            vds = filter_to_autosomes(
                vds
            )  # TODO: Switch to hl.vds.filter_to_autosomes if/when available

        ht = vds.variant_data.rows()
        variant_filter_expr = True
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
            )  # TODO: Create a checkpointed HT list of "good" intervals to use here instead
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
                ukbb_interval_qc_path("broad", 7, "autosomes")
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
                f"ccdg_{data_type}_was_split": mt.was_split,
                f"ccdg_{data_type}_AC": hl.agg.sum(mt.LGT.n_alt_alleles()),
                f"ccdg_{data_type}_AN": hl.agg.count_where(hl.is_defined(mt.LGT)) * 2,
            }
        )
        ht = mt.rows().checkpoint(
            f"{get_sample_qc_root(data_type=data_type, mt=False)}/variant_ccdg_{data_type}_af_callrate.ht",
            overwrite=overwrite,
            _read_if_exists=(not overwrite),
        )

        return ht

    logger.info(
        "Creating Table with joint gnomAD v3.1.2 and CCDG genome allele frequencies and callrate...",
        data_type,
    )
    ccdg_exomes_ht = _initial_filter("exomes")
    ccdg_genomes_ht = _initial_filter("genomes")
    ht = ccdg_exomes_ht.join(ccdg_genomes_ht, how="inner")
    ht = ht.annotate(**gnomad_ht[mt.row_key], **ukbb_ht[mt.row_key])
    ht = ht.annotate(
        joint_biallelic=~ht.ccdg_genomes_was_split & ~ht.gnomad_was_split,
        joint_AC=ht.ccdg_genomes_AC + ht.gnomad_AC,
        joint_AN=ht.ccdg_genomes_AN + ht.gnomad_AN,
    )
    total_genome_an = (gnomad_ht.freq_sample_count + ccdg_genome_count) * 2
    ht = ht.annotate(
        joint_AF=ht.joint_AC / ht.joint_AN, joint_callrate=ht.joint_AN / total_genome_an
    )
    ht = ht.checkpoint(
        f"{get_sample_qc_root(data_type='', mt=False)}ancestry_pca_joint.ht",
        overwrite=overwrite,
        _read_if_exists=(not overwrite),
    )

    logger.info(
        "Filtering variants to combined gnomAD v3.1.2 and CCDG genome AF of %d and callrate of %d, CCDG exome callrate "
        "of %d, and UK Biobank exome callrate of %d....",
        af_cutoff,
        callrate_cutoff,
        ccdg_exome_callrate_cutoff,
        ukbb_exome_callrate_cutoff,
    )

    variant_filter_expr = True
    if bi_allelic_only:
        variant_filter_expr = ht.joint_biallelic

    variant_filter_expr &= (
        (ht.joint_AF > af_cutoff)
        & (ht.joint_callrate > callrate_cutoff)
        & (ht.ccdg_exome_AN / (ccdg_exome_count * 2) > ccdg_exome_callrate_cutoff)
        & (ht.ukbb_AN / (ukbb_exome_count * 2) > ukbb_exome_callrate_cutoff)
    )
    ht = ht.annotate_globals(
        autosomes_only=autosomes_only,
        snv_only=snv_only,
        bi_allelic_only=bi_allelic_only,
        gnomad_v3_ac_filter=gnomad_v3_ac_filter,
        high_qual_ccdg_exome_interval_only=high_qual_ccdg_exome_interval_only,
        high_qual_ukbb_exome_interval_only=high_qual_ukbb_exome_interval_only,
        af_cutoff=af_cutoff,
        callrate_cutoff=callrate_cutoff,
        ccdg_exome_callrate_cutoff=ccdg_exome_callrate_cutoff,
        ukbb_exome_callrate_cutoff=ukbb_exome_callrate_cutoff,
    )
    ht.filter(variant_filter_expr).checkpoint(
        get_pca_variants_path_ht(ld_pruned=False),
        overwrite=overwrite,
        _read_if_exists=(not overwrite),
    )

    if ld_pruning:
        logger.info("Creating Table after LD pruning of %s...", ld_pruning_dataset)
        ht = hl.read_table(get_pca_variants_path_ht(ld_pruned=False))
        if ld_pruning_dataset == "ccdg_genomes":
            vds = hl.vds.read_vds(get_ccdg_vds_path("genomes"))
            vds = hl.vds.split_multi(vds)
            vds = hl.vds.filter_variants(vds, ht, keep=True)
            mt = hl.vds.to_dense_mt(vds)
        elif ld_pruning_dataset == "gnomad_genomes":
            mt = get_gnomad_v3_mt(split=True)
            mt = densify_sites(mt, ht, hl.read_table(gnomad_last_END_position.path))
            mt = filter_ref_blocks(mt)
        else:
            ValueError(
                "Only options for LD pruning are 'ccdg_genomes' and `gnomad_genomes`"
            )

        ht = hl.ld_prune(mt.GT, r2=ld_r2)
        ht = ht.annotate_globals(ld_r2=ld_r2, ld_pruning_dataset=ld_pruning_dataset)
        ht.checkpoint(
            get_pca_variants_path_ht(data=ld_pruning_dataset, ld_pruned=True),
            overwrite=overwrite,
            _read_if_exists=(not overwrite),
        )


def main(args):
    hl.init(log=f"/variant_filter.log")

    determine_pca_variants(
        autosomes_only=args.autosomes_only,
        bi_allelic_only=args.bi_allelic_only,
        snv_only=args.snv_only,
        gnomad_v3_ac_filter=args.gnomad_v3_ac_filter,
        high_qual_ccdg_exome_interval_only=args.high_qual_ccdg_interval_only,
        high_qual_ukbb_exome_interval_only=args.high_qual_ukbb_interval_only,
        af_cutoff=args.af_cutoff,
        callrate_cutoff=args.callrate_cutoff,
        ccdg_exome_callrate_cutoff=args.ccdg_exome_callrate_cutoff,
        ukbb_exome_callrate_cutoff=args.ukbb_exome_callrate_cutoff,
        ld_pruning=args.ld_pruning,
        ld_pruning_dataset=args.ld_pruning_dataset,
        ld_r2=args.ld_r2,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--autosomes-only", help="Filter to autosomes", action="store_true"
    )
    parser.add_argument("--snv-only", help="Filter to SNVs", action="store_true")
    parser.add_argument(
        "--bi-allelic-only",
        help="Filter to variants that are bi-allelic in both CCDG and gnomAD v3",
        action="store_true",
    )
    parser.add_argument(
        "--gnomad-v3-ac-filter",
        type=int,
        help="Filter to variants with AC above this value in gnomAD v3",
        default=10,
    )
    parser.add_argument(
        "--high-qual-ccdg-interval-only",
        help="Filter to high quality interval in CCDG exomes",
        action="store_true",
    )
    parser.add_argument(
        "--high-qual-ukbb-interval-only",
        help="Filter to high quality interval in UKBB exomes",
        action="store_true",
    )
    parser.add_argument(
        "--af-cutoff",
        type=float,
        help="Filter to variants with combined MAF above this value",
        default=0.001,
    )
    parser.add_argument(
        "--callrate-cutoff",
        type=float,
        help="Filter to variants with combined callrate above this value",
        default=0.99,
    )
    parser.add_argument(
        "--ccdg-exome-callrate-cutoff",
        type=float,
        help="Filter to variants with callrate above this value in CCDG exomes",
        default=0.99,
    )
    parser.add_argument(
        "--ukbb-exome-callrate-cutoff",
        type=float,
        help="Filter to variants with callrate above this value in UKBB exomes",
        default=0.99,
    )
    parser.add_argument("--ld-pruning", help="Apply LD pruning", action="store_true")
    parser.add_argument(
        "--ld--pruning-dataset",
        type=str,
        help="Dataset to apply LD pruning with",
        default="ccdg_genomes",
        choices=["ccdg_genomes", "gnomad_genomes"]
    )
    parser.add_argument("--ld-r2", type=float, help="LD pruning cutoff", default=0.1)
    parser.add_argument("--overwrite", help="Overwrite", action="store_true")
    args = parser.parse_args()
    main(args)
