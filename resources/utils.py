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


# TODO: Might need to rethink how we define "high quality" CCDG exome intervals
# TODO: Wait for an updated interval matrix table with DP included (n_bases_over_20x per sample per interval)
def exomes_qc_intervals_ht(
    pct_bases_20x: float = 0.8,
    overwrite: bool = False,
) -> hl.Table:
    """
    Generate CCDG exomes interval table for exomes QC.

    :param pct_bases_20x: Percent of bases with coverage greater than 20x over the interval
    :param overwrite: Whether to overwrite CCDG exomes interval QC HT
    :return: CCDG exomes interval table
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
            int_mt.interval.start.position + 1,
        ),
    )
    int_mt = int_mt.annotate_rows(int_len=int_mt.end_pos - int_mt.start_pos)
    int_mt = int_mt.annotate_cols(center=ht[int_mt.col_id]["cohort"])
    int_mt = int_mt.filter_cols(int_mt.center != "sccs")
    int_mt = int_mt.annotate_entries(
        pct_bases_20x=int_mt.n_bases_over_20x / int_mt.int_len
    ) # presume the field will be called 'n_bases_over_20x'

    n_broad = int_mt.count_cols()
    int_mt = int_mt.annotate_entries(int_defined=(int_mt.pct_bases_20x > pct_bases_20x))
    int_mt = int_mt.annotate_rows(
        pct_broad_defined=hl.agg.count_where(int_mt.int_defined) / n_broad
    ).checkpoint(
        get_ccdg_results_path(data_type="exomes", result="intervals"),
        overwrite=overwrite,
        _read_if_exists=(not overwrite),
    )
    return int_ht
