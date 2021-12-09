from .results import *

# TODO: Coordinate with Tim to make a hl.vds.filter_to_autosomes that we can use instead
# https://github.com/broadinstitute/gnomad_methods/blob/3536f87e249f0804f6762facce468597a9c441c6/gnomad/utils/filtering.py#L180
def filter_to_autosomes(
    mtds: Union[hl.Table, hl.MatrixTable, hl.vds.VariantDataset]
) -> Union[hl.Table, hl.MatrixTable, hl.vds.VariantDataset]:
    """
    Filter Table, MatrixTable or VariantDataset to autosome contigs only.

    This assumes that the input MT or VDS variant_data MT contains a field named `locus` of type Locus

    :param mtds: Input MatrixTable/Table/VariantDataset
    :return: MatrixTable/Table/VariantDataset subset to autosomes
    """
    if isinstance(mtds, hl.vds.VariantDataset):
        reference = get_reference_genome(mtds.variant_data.locus)
    else:
        reference = get_reference_genome(mtds.locus)
    autosomes = hl.parse_locus_interval(
        f"{reference.contigs[0]}-{reference.contigs[21]}", reference_genome=reference
    )

    if isinstance(mtds, hl.vds.VariantDataset):
        return hl.vds.filter_intervals(mtds, [autosomes])
    else:
        return hl.filter_intervals(mtds, [autosomes])


def annotate_exomes_interval_mt(pct_bases_defined: float = 0.8):
    int_mt = hl.read_matrix_table(exomes_interval_mt_path)
    ht = get_sample_manifest_ht("exomes")
    lst = hl.import_locus_intervals(interval_path, reference_genome="GRCh38")

    int_mt = int_mt.annotate_rows(
        end_pos=hl.if_else(
            int_mt.interval.includes_end,
            int_mt.interval.end.position + 1,
            int_mt.interval.end.position,
        ),
        start_pos=hl.if_else(
            int_mt.interval.includes_start,
            int_mt.interval.start.position,
            int_mt.interval.start.position + 1,
        ),
    )
    int_mt = int_mt.annotate_rows(
        int_len=int_mt.end_pos - int_mt.start_pos,
        target=lst.index(int_mt.interval, all_matches=True).target,
    )
    int_mt = int_mt.annotate_cols(center=ht[int_mt.col_id]["cohort"])
    int_mt = int_mt.annotate_cols(
        washu=hl.if_else(
            (hl.is_defined(int_mt.center) & (int_mt.center == "sccs")), True, False
        )
    )

    int_mt = int_mt.annotate_entries(
        **{
            f"pct_bases_{INTERVAL_DP}x": int_mt[f"n_bases_dp_gte_{INTERVAL_DP}"]
            / int_mt.int_len
        }
    )
    int_mt = int_mt.annotate_cols(
        **{
            f"total_bases_{INTERVAL_DP}x": hl.agg.sum(
                int_mt[f"n_bases_dp_gte_{INTERVAL_DP}"]
            )
        }
    )
    int_mt = int_mt.annotate_rows(
        **{
            f"mean_pct_bases_{INTERVAL_DP}x": hl.agg.mean(
                int_mt[f"pct_bases_{INTERVAL_DP}x"]
            )
        },
        **{
            f"pct_broad_defined_{INTERVAL_DP}x": hl.agg.filter(
                ~(int_mt.washu),
                hl.agg.fraction(
                    int_mt[f"pct_bases_{INTERVAL_DP}x"] > pct_bases_defined
                ),
            )
        },
        **{
            f"pct_washu_defined_{INTERVAL_DP}x": hl.agg.filter(
                int_mt.washu,
                hl.agg.fraction(
                    int_mt[f"pct_bases_{INTERVAL_DP}x"] > pct_bases_defined
                ),
            )
        },
    )
    return int_mt


def interval_target_sum_ht():
    int_ht = hl.read_table(
        get_ccdg_results_path(data_type="exomes", result=f"intervals_{INTERVAL_DP}x")
    )
    int_ht = int_ht.explode(int_ht.target)
    int_ht = int_ht.annotate(target2=int_ht.target.split("\|"))
    int_ht = int_ht.explode(int_ht.target2)
    target_ht = int_ht.group_by("target2").aggregate(
        total_len=hl.agg.sum(int_ht.int_len),
        filtered_len=hl.agg.filter(int_ht.to_keep, hl.agg.sum(int_ht.int_len)),
    )
    target_ht = target_ht.annotate(
        percent_len=target_ht.filtered_len / target_ht.total_len
    )
    target_ht = target_ht.order_by(hl.desc(target_ht.total_len))
    return target_ht


def ccdg_interval_qc_ht(
    pct_samples_lower: float = 0.8,
    overwrite: bool = False,
) -> hl.Table:
    """
    Generate CCDG exomes interval table for exomes QC.

    :param pct_bases_lower: Percent of bases covered over this value in this the interval
    :param overwrite: Whether to overwrite CCDG exomes interval QC HT
    :return: CCDG exomes interval table
    """
    int_ht = annotate_exomes_interval_mt().rows()

    int_ht = int_ht.annotate(
        to_keep=int_ht[f"pct_broad_defined_{INTERVAL_DP}x"] > pct_samples_lower
    )
    int_ht = int_ht.checkpoint(
        get_ccdg_results_path(data_type="exomes", result=f"intervals_{INTERVAL_DP}x"),
        overwrite=overwrite,
        _read_if_exists=(not overwrite),
    )
    return int_ht
