ccdg_bucket = "gs://ccdg"
exomes_interval_mt_path = "gs://ccdg-30day-temp/dking/iG4Vz/total-final-result.mt"


def get_ccdg_vds_path(data_type: str = "genomes") -> str:
    """
    Return path to CCDG VDS.

    :param data_type: Whether data is from CCDG genomes or exomes, default is genomes
    :return: path to CCDG vds
    """
    data_type = (
        "wgs_136k_recombine" if data_type == "genomes" else "split_200k_ccdg_exome"
    )
    return f"{ccdg_bucket}/vds/{data_type}.vds"


def get_sample_qc_root(data_type: str = "genomes", mt: bool = False) -> str:
    """
    Return path to CCDG sample QC root folder.

    :param data_type: Whether data is from CCDG genomes or exomes, default is 'genomes'
    :param mt: Whether path is for a MatrixTable, default is False
    :return: root of CCDG sample QC path
    """
    return f"{ccdg_bucket}/sample_qc/{'mt' if mt else 'ht'}/{data_type}"


def get_sample_manifest_ht(data_type: str = "exomes") -> str:
    """
    Return path to CCDG sample manifest table.

    :param data_type: Whether data is for CCDG genomes or exomes, default is 'genomes'
    :return: path to CCDG sample manifest table
    """
    data_type = "wgs" if data_type == "genomes" else "wes"
    ht = hl.read_table(f"{ccdg_bucket}/wlu/ccdg_{data_type}_sample.ht")
    return ht


def get_joint_pca_variants_ht_path() -> str:
    """
    Return path to a table of variants with joint gnomAD v3.1.2 and CCDG genome allele frequencies and callrate:

    :return: path to variants with joint gnomAD v3.1.2 and CCDG genome allele frequencies and callrate
    """
    return f"{get_sample_qc_root(data_type='', mt=False)}ancestry_pca_joint.ht"


def get_pca_variants_path(
    ld_pruned: bool = True,
    data: str = "ccdg_genomes",
    mt: bool = False,
) -> str:
    """
    Return path to filtered variants for ancestry PCA.

    :param ld_pruned: Whether the variants have been ld pruned
    :param data: Whether data is from CCDG genomes or gnomAD v3 genomes, default is 'ccdg_genomes'
    :param mt: Whether path is for a MatrixTable, default is False -> hl.Table
    :return: path to filtered MatrixTable or Table
    """
    ld_pruning_flag = f"{data}_ld_pruned" if ld_pruned else "pre_ld_pruning"
    return f"{get_sample_qc_root(data_type='', mt=mt)}{ld_pruning_flag}_combined_variants.{'mt' if mt else 'ht'}"


def get_ccdg_results_path(
    data_type: str = "genomes", mt: bool = False, result: str = None
) -> str:
    """
    Return path to CCDG results.
            Available results are:
            - "pre_filtered_variants" - generated from pca_variant_filter.py -> determine_pca_variants() -> _initial_filter()
            - "high_qual_intervals"
            - "sample_qc_all"
            - "sample_qc_bi_allelic"
            - "sample_qc_multi_allelic"
            - "sex"
            - "pc_scores"
            - "relatedness"

    :param data_type: Whether data is from CCDG genomes or exomes, default is 'genomes'
    :param mt: Whether path is for a MatrixTable, default is False -> hl.Table
    :param result: Type of result
    :return: path to CCDG result MatrixTable or Table
    """
    return f"{get_sample_qc_root(data_type=data_type, mt=mt)}/ccdg_{data_type}_{result}.{'mt' if mt else 'ht'}"
