# TODO: Coordinate with Tim to make a hl.vds.filter_intervals that we can use instead
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
        get_ccdg_results_path(data_type="exomes", result="high_qual_intervals"),
        overwrite=overwrite,
        _read_if_exists=(not overwrite),
    )

    return int_ht
