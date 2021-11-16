ccdg_bucket = "gs://ccdg"
exomes_interval_mt_path = "gs://ccdg-30day-temp/dking/iG4Vz/total-final-result.mt"


def get_ccdg_vds_path(data_type: str = "genomes") -> str:
    """
    Return path to CCDG vds
    :param data_type: Whether data is from genomes or exomes, default is genomes
    :return: path to CCDG vds
    """
    data_type = (
        "wgs_136k_recombine" if data_type == "genomes" else "split_200k_ccdg_exome"
    )
    return f"{ccdg_bucket}/vds/{data_type}.vds"


def get_sample_qc_root(data_type: str = "genomes", mt: bool = False) -> str:
    """
    Return path to CCDG sample QC root folder
    :param data_type: Whether data is from genomes or exomes, default is 'genomes'
    :param mt: Whether path is for a MatrixTable, default is False
    :return: root of CCDG sample QC path
    """
    return f"{ccdg_bucket}/sample_qc/{'mt' if mt else 'ht'}/{data_type}"


def get_sample_manifest_ht(data_type: str = "exomes") -> str:
    """
    Return path to CCDG sample manifest table
    :param data_type: Whether data is for genomes or exomes, default is 'genomes'
    :return: path to CCDG sample manifest table
    """
    data_type = "wgs" if data_type == "genomes" else "wes"
    ht = hl.read_table(f"{ccdg_bucket}/wlu/ccdg_{data_type}_sample.ht")
    return ht


def get_pca_variants_path_ht(ld_pruned: bool = True, data: str = "ccdg_genomes") -> str:
    """
    Return path to filtered variants
    :param ld_pruned: Whether the variatns have been ld pruned
    :return: path to filtered variant ht
    """
    ld_pruning_flag = f"{data}_ld_pruned" if ld_pruned else "pre_ld_pruning"
    return f'{get_sample_qc_root(data_type="", mt=False)}{ld_pruning_flag}_combined_variants.ht'
