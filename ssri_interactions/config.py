from pathlib import Path
import os


def _get_base_dir(project_dirname: str) -> Path:
    paths = [
        p
        for p in Path(__file__).absolute().parents
        if p.name.lower() == project_dirname.lower()
    ]
    if len(paths) == 0:
        basepath = Path(os.getcwd())
    else:
        basepath = paths[0]
    return basepath


def get_default_data_dir(project_dirname="SSRI Interactions") -> Path:
    basepath = _get_base_dir(project_dirname)
    return basepath / "data"


def get_default_fig_dir(project_dirname="SSRI Interactions") -> Path:
    basepath = _get_base_dir(project_dirname)
    return basepath / "figs"


def get_default_table_dir(project_dirname="SSRI Interactions") -> Path:
    basepath = _get_base_dir(project_dirname)
    return basepath / "tables"


def get_default_derived_data_dir(project_dirname="SSRI Interactions") -> Path:
    return get_default_data_dir(project_dirname=project_dirname) / "derived"


class Config:
    data_dir = get_default_data_dir("SSRI Interactions")
    derived_data_dir = get_default_derived_data_dir("SSRI Interactions")
    fig_dir = get_default_fig_dir("SSRI Interactions")
    table_dir = get_default_table_dir("SSRI Interactions")

    eeg_states = ("sw", "act")


class ExperimentInfo:
    block_names = (
        "pre",
        "base_shock",
        "post_base_shock",
        "chal",
        "chal_shock",
        "post_chal_shock",
        "way",
        "pre",
        "chal",
        "way",
    )
    group_names = (
        "citalopram_continuation",
        "chronic_citalopram",
        "chronic_saline",
        "chronic_saline_",
        "citalopram_discontinuation",
    )
    cit_sessions = (
        "chronic_08",
        "chronic_01",
        "hamilton_02",
        "hamilton_14",
        "hamilton_23",
        "hamilton_08",
        "hamilton_20",
    )
    sal_sessions = (
        "hamilton_01",
        "hamilton_13",
        "hamilton_07",
        "hamilton_19",
        "chronic_09",
    )
    discontunation_sessions = (
        "hamilton_12",
        "hamilton_05",
        "hamilton_18",
        "hamilton_17",
        "hamilton_25",
        "hamilton_26",
        "hamilton_28",
        "hamilton_29",
        "hamilton_30"
    )


    neuron_type_order = ["sr", "sir", "ff"]

    way_sessions = (
        "hamilton_09",
        "hamilton_31",
        "hamilton_38",
        "hamilton_37",
        "hamilton_36",
        "hamilton_32",
        "acute_15",
        "acute_01",
        "acute_12",
        "acute_11",
    )
