import hail as hl

ccdg_bucket = "gs://ccdg"
INTERVAL_DP = 10
exomes_interval_mt_path = (
    f"gs://ccdg/interval_qc/ccdg_exome_interval_qc_dp_gte_{INTERVAL_DP}.mt"
)
exomes_interval_list = "gs://gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.interval_list"
genomes_interval_ht_path = (
    "gs://ccdg/wlu/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list.ht"
)

path_to_gnomad_rf = (
    "gs://gnomad-public-requester-pays/release/3.1/pca/gnomad.v3.1.RF_fit.pkl"
)
path_to_gnomad_loadings = (
    "gs://gnomad-public-requester-pays/release/3.1/pca/gnomad.v3.1.pca_loadings.ht"
)


# def get_ccdg_vds_path(data_type: str = "genomes") -> str:
#     """
#     Return path to CCDG VDS.
#
#     :param data_type: Whether data is from CCDG genomes or exomes, default is 'genomes'
#     :return: path to CCDG VDS
#     """
#     data_type = (
#         "wgs_136k_recombine" if data_type == "genomes" else "split_200k_ccdg_exome"
#     )
#     return f"{ccdg_bucket}/vds/{data_type}.vds"


def get_ccdg_vds(data_type: str = "genomes") -> hl.vds.VariantDataset:
    """
    Return CCDG VDS.

    :param data_type: Whether data is from CCDG genomes or exomes, default is 'genomes'
    :return: CCDG VDS
    """
    vds_path = (
        f"{ccdg_bucket}/vds/wgs_136k_recombine.vds"
        if data_type == "genomes"
        else f"{ccdg_bucket}/vds/split_200k_ccdg_exome.vds"
    )
    if data_type == "genomes":
        vds = hl.vds.read_vds(vds_path)
    else:

        def read_ccdg_vds(path):
            refmt = hl.read_matrix_table(path + "/reference_data")
            refmt = refmt.transmute_rows(ref_allele=refmt.alleles[0][:1])
            return hl.vds.VariantDataset(
                refmt, hl.read_matrix_table(path + "/variant_data")
            )

        vds = read_ccdg_vds(vds_path)

    return vds


def get_calling_intervals_ht(data_type) -> hl.Table:
    """
    Return calling interval Hail Table.

    :param data_type: Whether data is from CCDG genomes or exomes, default is 'genomes'
    :return: calling interval Hail Table
    """
    if data_type == "genomes":
        calling_intervals = hl.read_table(genomes_interval_ht_path)
    else:
        calling_intervals = hl.read_table(
            get_ccdg_results_path(
                data_type="exomes", result=f"intervals_{INTERVAL_DP}x"
            )
        )
    return calling_intervals


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
            - "intervals"
            - "sample_qc_all"
            - "sample_qc_bi_allelic"
            - "sample_qc_multi_allelic"
            - "sex"
            - "pc_scores"
            - "relatedness"
            - "impute_ploidy"

    :param data_type: Whether data is from CCDG genomes or exomes, default is 'genomes'
    :param mt: Whether path is for a MatrixTable, default is False -> hl.Table
    :param result: Type of result
    :return: path to CCDG result MatrixTable or Table
    """
    return f"{get_sample_qc_root(data_type=data_type, mt=mt)}/ccdg_{data_type}_{result}.{'mt' if mt else 'ht'}"
