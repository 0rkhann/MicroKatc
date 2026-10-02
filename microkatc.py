"""MicroKatc command line.

    python microkatc.py run <study.yaml> [--fresh]    run every analysis the study requests
    python microkatc.py check <study.yaml>            check the study without running anything

Results go to a results/ folder next to the study file. Exit codes: 0 success, 1 invalid study
or stale results, 2 failure during the run.

Run one study per process: COPASI resolves relative output paths against the first folder it
saved a model in, so a second study run in the same Python process would write into the first
study's results/ (measured 2026-10-02).
"""

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

from study import load_study
from study_yaml import StudyError

REPO = Path(__file__).resolve().parent
RUN_INFO = "run_info.json"


class StaleResultsError(Exception):
    """results/ was made from a different version of the study"""


def study_hash(study):
    """SHA-256 of the study file and of every energy input it uses"""
    digest = hashlib.sha256(study.path.read_bytes())
    if study.files is not None:
        for name in sorted(set(study.species()) | set(study.transition_states())):
            digest.update(name.encode())
            digest.update((study.files / f"{name}.out").read_bytes())
    return digest.hexdigest()


def _git_commit(folder):
    result = subprocess.run(
        ["git", "-C", str(folder), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    )
    return result.stdout.strip() or None


def _version(package):
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return None


def _write_run_info(results, info):
    (results / RUN_INFO).write_text(json.dumps(info, indent=2) + "\n")


def prepare_results(study, fresh=False):
    """Creates results/, or checks that the existing one was made from this exact study"""
    results = study.results_dir
    if fresh and results.exists():
        shutil.rmtree(results)
    digest = study_hash(study)
    info_path = results / RUN_INFO
    if info_path.is_file():
        if json.loads(info_path.read_text()).get("study_sha256") != digest:
            raise StaleResultsError(
                f"{results} was made from a different version of {study.path.name}; delete it or run with --fresh"
            )
    elif results.exists() and any(results.iterdir()):
        raise StaleResultsError(
            f"{results} has no {RUN_INFO}; delete it or run with --fresh"
        )
    results.mkdir(parents=True, exist_ok=True)
    thermochange = os.environ.get("thermochange")
    info = {
        "study": str(study.path),
        "study_sha256": digest,
        "microkatc_commit": _git_commit(REPO),
        "thermochange_commit": _git_commit(thermochange)
        if thermochange and study.files
        else None,
        "python": platform.python_version(),
        "python-copasi": _version("python-copasi"),
        "numpy": _version("numpy"),
        "scipy": _version("scipy"),
        "started": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "finished": None,
        "duration_s": None,
    }
    _write_run_info(results, info)
    return results, info


def initial_concentrations(study, initial_M, studied_value):
    """Initial concentrations for one simulation: the analysis' values, product at 0, studied species set"""
    return {
        **initial_M,
        study.product: 0.0,
        study.studied_species: float(studied_value),
    }


def run_analyses(study):
    """Runs every requested analysis in the current folder (results/)"""
    import matplotlib

    matplotlib.use("Agg")
    import numpy as np

    from apparent_activation_energy import ApparentEaAnalysis, DRCAnalysis
    from auxiliary_functions import AuxiliaryFunctions
    from microkinetics_simulation import MicroKinetics
    from plotting_functions import PlotFunctions
    from thermochemistry import write_barrier_tables

    for message in write_barrier_tables(study):
        print(message)
    reactions = study.reactions()
    concentrations = np.array(study.studied_concentrations_M)
    T, studied, grid = study.temperature_K, study.studied_species, PlotFunctions.grid

    ea = study.activation_energy
    if ea is not None:
        analysis = ApparentEaAnalysis(
            T,
            np.array(ea.temperatures_K),
            concentrations,
            studied,
            initial_concentrations(study, ea.initial_M, concentrations[0]),
            ea.simulation_time_s,
            ea.sampling_time_h,
            reactions,
            study.output_step_s,
        )
        df_flux, df_rate = analysis.df_flux, analysis.df_rate
        analysis.plot_ln_ri_vs_1_over_T(
            df_flux, concentrations[0], *grid(len(reactions)), (15, 15)
        )
        forming = [
            r for r, s in zip(reactions, study.steps) if study.product in s.species()
        ]
        analysis.plot_Ea_vs_c0_flux_based(
            df_flux, forming, *grid(len(forming)), (15, 10), log_x=True
        )
        analysis.plot_Ea_vs_c0_rate_based(
            df_rate, [study.product], 1, 1, (15, 10), log_x=True
        )

    drc = study.degree_of_rate_control
    if drc is not None:
        analysis = DRCAnalysis(
            concentrations,
            np.array([T]),
            drc.simulation_time_s,
            initial_concentrations(study, drc.initial_M, concentrations[0]),
            studied,
            reactions,
            drc.e_shift_kcal_mol,
            drc.cores,
            product=study.product,
            time_step=study.output_step_s,
        )
        analysis.plot_drc_vs_c0(
            *grid(len(reactions)), (15, 15), analysis.df_drc, T, log_x=True
        )

    mk = study.microkinetics
    if mk is not None:
        cycles = {label: list(members) for label, members in study.cycles.items()}
        analysis = MicroKinetics(
            T,
            mk.simulation_time_s,
            concentrations,
            studied,
            initial_concentrations(study, mk.initial_M, concentrations[0]),
            cycles,
            [
                study.max_product_M(initial_concentrations(study, mk.initial_M, c))
                for c in concentrations
            ],
            list(mk.plot_species),
            mk.catalyst_snapshot_h,
            mk.conversion,
            study.output_step_s,
            product=study.product,
        )
        analysis.plot_catalyst_concentration_vs_reactant(log_x=True, figsize=(15, 10))
        analysis.plot_concentration_evolution(
            T, True, True, *grid(len(concentrations)), (20, 16)
        )
        analysis.plot_reactant_vs_product_conversion(log_x=True, figsize=(15, 10))
        shares = AuxiliaryFunctions.get_concentrations_of_catalyst(
            analysis.simulations_dfs, mk.catalyst_snapshot_h, cycles
        )
        thresholds = [m * mk.conversion for m in analysis.max_product_M]
        times = AuxiliaryFunctions.compute_time_of_product_conversion_given_reactant_concentration(
            analysis.simulations_dfs, concentrations, studied, thresholds, study.product
        )
        summary = {
            "concentrations_M": concentrations.tolist(),
            "catalyst_M": dict(zip(cycles, shares)),
            # None where the product never reached the threshold in the simulated time
            "t99_h": [{c: t for t, c in times}.get(c) for c in concentrations],
        }
        Path("microkinetics.json").write_text(json.dumps(summary, indent=1) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="microkatc.py", description="Microkinetic analysis of catalytic cycles"
    )
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="run every analysis the study requests")
    run.add_argument("study")
    run.add_argument(
        "--fresh", action="store_true", help="delete the study's results/ first"
    )
    commands.add_parser(
        "check", help="check the study without running anything"
    ).add_argument("study")
    args = parser.parse_args(argv)

    try:
        study = load_study(args.study)
    except StudyError as error:
        print(f"{args.study}: invalid study", file=sys.stderr)
        for message in error.messages:
            print(f"  {message}", file=sys.stderr)
        return 1
    for warning in study.warnings:
        print(warning)
    if args.command == "check":
        if study.files is None:  # typed energies: the barriers cost nothing to check
            from thermochemistry import barrier_messages, barrier_table, gibbs_energies

            gibbs = gibbs_energies(study, study.temperature_K)
            table = barrier_table(study, gibbs)
            for message in barrier_messages(study, gibbs, table, study.temperature_K):
                print(message)
        print(
            f"{args.study}: OK ({len(study.steps)} steps, {len(study.species())} species)"
        )
        return 0

    try:
        results, info = prepare_results(study, args.fresh)
    except StaleResultsError as error:
        print(error, file=sys.stderr)
        return 1
    started, previous = time.time(), os.getcwd()
    os.chdir(results)
    try:
        run_analyses(study)
    except Exception as error:
        traceback.print_exc()
        print(f"{study.path}: run failed: {error}", file=sys.stderr)
        return 2
    finally:
        os.chdir(previous)
    info.update(
        finished=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        duration_s=round(time.time() - started, 1),
    )
    _write_run_info(results, info)
    print(f"Results in {results}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
