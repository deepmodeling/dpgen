"""Shared dpdata-based first-principles post-processing."""

import glob
import os

import dpdata
import numpy as np

_FP_CONFIG = {
    "pwscf": {"format": "qe/pw/scf", "source": "output", "type_map": True},
    "abacus": {
        "format": "abacus/scf",
        "source": ".",
        "marker": "INPUT",
        "normalize_type_map": True,
    },
    "siesta": {"format": "siesta/output", "source": "output"},
    "gaussian": {
        "format": "gaussian/log",
        "source": "output",
        "type_map": True,
        "multi_param": "use_clusters",
        "atom_pref_param": "use_atom_pref",
    },
    "cp2k": {
        "format": "cp2kdata/e_f",
        "source": "output",
        "type_map": True,
        "multi": True,
    },
    "pwmat": {
        "format": "pwmat/output",
        "source": "OUT.MLMD",
        "type_map": True,
        "one_frame": True,
    },
    "cpx": {
        "format": "qe/cp/traj",
        "source_param": "input_fn",
        "marker": "output",
        "multi": True,
    },
    "custom": {
        "format_param": "output_fmt",
        "source_param": "output_fn",
        "multi": True,
    },
}


def _group_fp_tasks(work_path: str) -> dict[str, list[str]]:
    """Group FP task directories by DP-GEN system index."""
    grouped = {}
    for task in sorted(glob.glob(os.path.join(work_path, "task.*"))):
        system_index = os.path.basename(task).split(".")[1]
        grouped.setdefault(system_index, []).append(task)
    return grouped


def _config_value(config: dict, name: str, jdata: dict) -> str:
    """Resolve a static config value or one stored in ``fp_params``."""
    if name in config:
        return config[name]
    return jdata["fp_params"][config[f"{name}_param"]]


def _source_for_task(config: dict, task: str, jdata: dict) -> str | None:
    """Return the file or directory consumed by dpdata for one FP task."""
    source_name = _config_value(config, "source", jdata)
    marker = config.get("marker", source_name)
    if not os.path.exists(os.path.join(task, marker)):
        return None
    return os.path.normpath(os.path.join(task, source_name))


def post_fp_dpdata(work_path: str, jdata: dict) -> None:
    """Collect a dpdata-backed FP style into DeepMD training data."""
    config = _FP_CONFIG[jdata["fp_style"]]
    output_fmt = _config_value(config, "format", jdata)
    type_map = jdata["type_map"] if config.get("type_map") else None
    atom_pref_param = config.get("atom_pref_param")
    use_multi_systems = config.get("multi", jdata.get(config.get("multi_param"), False))

    for system_index, tasks in _group_fp_tasks(work_path).items():
        all_sys = None
        for task in tasks:
            source = _source_for_task(config, task, jdata)
            if source is None:
                continue

            system = dpdata.LabeledSystem(source, fmt=output_fmt, type_map=type_map)
            nframes = len(system)
            if not nframes or (config.get("one_frame") and nframes != 1):
                continue
            if config.get("normalize_type_map"):
                system.data["atom_types"] = np.asarray(
                    system.data["atom_types"], dtype=int
                )
                system.apply_type_map(jdata["type_map"])
            if atom_pref_param and jdata.get(atom_pref_param, False):
                system.data["atom_pref"] = np.load(os.path.join(task, "atom_pref.npy"))

            if all_sys is None:
                all_sys = (
                    dpdata.MultiSystems(system, type_map=jdata["type_map"])
                    if use_multi_systems
                    else system
                )
            else:
                all_sys.append(system)

        if all_sys is not None:
            sys_data_path = os.path.join(work_path, f"data.{system_index}")
            all_sys.to_deepmd_raw(sys_data_path)
            all_sys.to_deepmd_npy(sys_data_path, set_size=len(tasks))
