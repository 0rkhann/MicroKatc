"""A MicroKatc study: the YAML file read, checked and turned into one Study object (spec sections 3 and 5)"""

import difflib
import math
import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from steps import StepSyntaxError, parse_step
from study_yaml import StudyError, load_yaml

RANGE = {"from": None, "to": None, "points": None}
SCHEMA = {
    "name": None,
    "species": {
        "files": None,
        "energies_kcal_mol": None,
        "product": None,
        "overall_reaction": None,
        "cycles": None,
        "vibrational_correction": None,
    },
    "steps": None,
    "conditions": {
        "temperature_K": None,
        "studied_species": None,
        "studied_range_M": RANGE,
        "output_step_s": None,
    },
    "analyses": {
        "activation_energy": {
            "initial_M": None,
            "temperatures_K": RANGE,
            "simulation_time_s": None,
            "sampling_time_h": None,
        },
        "degree_of_rate_control": {
            "initial_M": None,
            "e_shift_kcal_mol": None,
            "simulation_time_s": None,
            "cores": None,
        },
        "microkinetics": {
            "initial_M": None,
            "simulation_time_s": None,
            "catalyst_snapshot_h": None,
            "conversion": None,
            "plot_species": None,
        },
    },
}
THERMOCHANGE_SCRIPT = os.path.join("formatters", "formatted_energy_outputter.sh")


@dataclass(frozen=True)
class ActivationEnergy:
    initial_M: dict
    temperatures_K: tuple
    simulation_time_s: float
    sampling_time_h: float


@dataclass(frozen=True)
class DegreeOfRateControl:
    initial_M: dict
    e_shift_kcal_mol: float
    simulation_time_s: float
    cores: int


@dataclass(frozen=True)
class Microkinetics:
    initial_M: dict
    simulation_time_s: float
    catalyst_snapshot_h: float
    conversion: float
    plot_species: tuple


@dataclass(frozen=True)
class Study:
    """Everything the analyses need, validated; built only by load_study"""

    path: Path
    name: str
    files: object  # Path to the output-file folder, or None with typed energies
    energies_kcal_mol: object  # {species: {temperature_K: G}}, or None with files
    product: str
    overall: object  # steps.Step
    cycles: dict
    vibrational_correction: str
    steps: tuple
    temperature_K: float
    studied_species: str
    studied_concentrations_M: tuple
    output_step_s: float
    activation_energy: object = None
    degree_of_rate_control: object = None
    microkinetics: object = None
    warnings: tuple = field(default_factory=tuple)

    @property
    def results_dir(self):
        return self.path.parent / "results"

    def species(self):
        """Every species in the steps, first-seen order (transition states excluded)"""
        seen = []
        for step in self.steps:
            seen += [name for name in step.species() if name not in seen]
        return seen

    def transition_states(self):
        return [step.ts for step in self.steps if step.ts is not None]

    def temperatures(self):
        """Every temperature at which Gibbs energies are needed, as floats"""
        needed = {self.temperature_K}
        if self.activation_energy is not None:
            needed.update(self.activation_energy.temperatures_K)
        return sorted(needed)

    def reactions(self):
        """The steps as copasi_helper equations, in file order (COPASI names them r01, r02, ...)"""
        return [step.copasi_equation() for step in self.steps]

    def max_product_M(self, initial_M):
        """Largest product concentration the overall reaction allows from these initial concentrations"""
        limiting = min(
            initial_M.get(name, 0.0) / n for name, n in self.overall.reactants
        )
        return limiting * self.overall.coefficient(self.product, "products")


def _unknown_keys(doc, schema, prefix, errors):
    for key, value in doc.items():
        where = f"{prefix}{key}"
        if key not in schema:
            close = difflib.get_close_matches(str(key), [str(k) for k in schema], n=1)
            hint = f"; did you mean {close[0]}?" if close else ""
            errors.append(f"{where}: unknown key{hint}")
        elif isinstance(schema[key], dict) and isinstance(value, dict):
            _unknown_keys(value, schema[key], f"{where}.", errors)


def _number(value, where, errors, positive=True, integer=False):
    ok_type = isinstance(value, int) if integer else isinstance(value, (int, float))
    if (
        isinstance(value, bool)
        or not ok_type
        or (isinstance(value, float) and not math.isfinite(value))
    ):
        errors.append(
            f"{where}: must be {'an integer' if integer else 'a number'}, got {value!r}"
        )
        return None
    if positive and value <= 0:
        errors.append(f"{where}: must be positive, got {value}")
        return None
    return float(value) if not integer else int(value)


def _range(value, where, errors, minimum_points=1, log=False):
    if not isinstance(value, dict) or set(value) != {"from", "to", "points"}:
        errors.append(f"{where}: must be {{from: ..., to: ..., points: ...}}")
        return None
    start = _number(value["from"], f"{where}.from", errors)
    stop = _number(value["to"], f"{where}.to", errors)
    points = _number(value["points"], f"{where}.points", errors, integer=True)
    if None in (start, stop, points):
        return None
    if points < minimum_points:
        errors.append(f"{where}.points: needs at least {minimum_points}, got {points}")
        return None
    if points > 1 and not start < stop:
        errors.append(f"{where}: from must be smaller than to")
        return None
    if log:
        return tuple(
            float(v) for v in np.logspace(np.log10(start), np.log10(stop), points)
        )
    return tuple(float(v) for v in np.linspace(start, stop, points))


def _concentrations(value, where, errors):
    if not isinstance(value, dict) or not value:
        errors.append(f"{where}: must map species to initial concentrations in M")
        return {}
    out = {}
    for name, c in value.items():
        number = _number(c, f"{where}.{name}", errors, positive=False)
        if number is not None and number < 0:
            errors.append(f"{where}.{name}: must not be negative")
        elif number is not None:
            out[str(name)] = number
    return out


def _multiple(seconds, step_s):
    return abs(seconds / step_s - round(seconds / step_s)) <= 1e-9


def _analysis(analyses, key, errors):
    """An analysis's settings; an empty key gives {} so every required setting is reported"""
    if key not in analyses:
        return None
    raw = analyses[key]
    if raw is None:
        return {}
    if not isinstance(raw, dict):
        errors.append(f"analyses.{key}: must be a mapping of settings")
        return None
    return raw


def _total_on_grid(total_s, step_s, where, errors):
    """COPASI spaces output points by total/round(total/step), so off-grid totals shift every time"""
    if None not in (total_s, step_s) and not _multiple(total_s, step_s):
        errors.append(
            f"{where}.simulation_time_s: {total_s:g} s is not a multiple of output_step_s ({step_s:g} s)"
        )


def _on_grid(hours, step_s, total_s, where, errors):
    seconds = hours * 3600
    if seconds > total_s:
        errors.append(
            f"{where}: {hours:g} h is after the end of the simulation ({total_s:g} s = {total_s / 3600:.2g} h)"
        )
    elif not _multiple(seconds, step_s):
        errors.append(
            f"{where}: {hours:g} h is not a multiple of output_step_s ({step_s:g} s)"
        )


def _energies(
    value, where, errors, needed_temperatures, working_temperature, ea_requested
):
    """{species: {T: G}}; a single number is valid at the working temperature only"""
    if not isinstance(value, dict) or not value:
        errors.append(f"{where}: must map species to Gibbs energies in kcal/mol")
        return {}
    out = {}
    for name, entry in value.items():
        if isinstance(entry, dict):
            table = {}
            for temperature, g in entry.items():
                t = _number(temperature, f"{where}.{name} temperature", errors)
                v = _number(g, f"{where}.{name}.{temperature}", errors, positive=False)
                if t is not None and v is not None:
                    table[t] = v
        else:
            v = _number(entry, f"{where}.{name}", errors, positive=False)
            table = {working_temperature: v} if v is not None else {}
        missing = [
            t
            for t in needed_temperatures
            if not any(abs(t - known) < 1e-6 for known in table)
        ]
        if missing and table:
            purpose = "activation_energy" if ea_requested else "the working temperature"
            listed = ", ".join(f"{t:g}" for t in missing)
            errors.append(f"{where}.{name}: needs values at {listed} K for {purpose}")
        out[str(name)] = table
    return out


def load_study(path, check_energy_sources=True):
    """Reads and checks a study file; raises StudyError listing every problem found.

    check_energy_sources=False skips the output-file and $thermochange checks, for scripts that
    only read the results of a finished run.
    """
    path = Path(path).resolve()
    doc = load_yaml(path)
    errors, warnings = [], []
    _unknown_keys(doc, SCHEMA, "", errors)

    species_block = doc.get("species") if isinstance(doc.get("species"), dict) else {}
    conditions = (
        doc.get("conditions") if isinstance(doc.get("conditions"), dict) else {}
    )
    analyses = doc.get("analyses") if isinstance(doc.get("analyses"), dict) else {}
    for key in ("species", "steps", "conditions", "analyses"):
        if key not in doc:
            errors.append(f"{key}: required")
    for key in ("product", "overall_reaction"):
        if key not in species_block:
            errors.append(f"species.{key}: required")
    if ("files" in species_block) == ("energies_kcal_mol" in species_block):
        errors.append(
            "species: give files or energies_kcal_mol, not both"
            if "files" in species_block
            else "species: give files or energies_kcal_mol"
        )
    if analyses is not None and not any(k in analyses for k in SCHEMA["analyses"]):
        errors.append(
            "analyses: request at least one of activation_energy, degree_of_rate_control, microkinetics"
        )

    # Steps
    steps, raw_steps = [], doc.get("steps", [])
    if not isinstance(raw_steps, list) or not raw_steps:
        errors.append("steps: must be a list of steps such as 'A + B <=> C via TS'")
        raw_steps = []
    seen = {}
    for i, text in enumerate(raw_steps):
        try:
            step = parse_step(text)
        except StepSyntaxError as error:
            errors.append(f"steps[{i}] {error}")
            continue
        left, right = tuple(sorted(step.reactants)), tuple(sorted(step.products))
        if left == right:
            errors.append(f"steps[{i}] {step.text}: both sides are the same")
        key = (frozenset((left, right)), step.ts)  # every step is reversible
        if key in seen:
            errors.append(f"steps[{i}] repeats steps[{seen[key]}] ({step.text})")
        seen.setdefault(key, i)
        steps.append(step)
    names = []
    for step in steps:
        names += [n for n in step.species() if n not in names]
    transition_states = [s.ts for s in steps if s.ts is not None]

    # Conditions
    temperature = (
        _number(conditions.get("temperature_K"), "conditions.temperature_K", errors)
        if "temperature_K" in conditions
        else None
    )
    if "temperature_K" not in conditions:
        errors.append("conditions.temperature_K: required")
    studied = conditions.get("studied_species")
    if studied is not None and not isinstance(studied, str):
        errors.append(
            f"conditions.studied_species: must be one species name, got {studied!r}"
        )
        studied = ""
    if studied is None:
        errors.append("conditions.studied_species: required")
    elif names and studied not in names:
        errors.append(
            f"conditions.studied_species: {studied} does not appear in any step"
        )
    concentrations = None
    if "studied_range_M" in conditions:
        concentrations = _range(
            conditions["studied_range_M"],
            "conditions.studied_range_M",
            errors,
            log=True,
        )
    else:
        errors.append("conditions.studied_range_M: required")
    output_step = _number(
        conditions.get("output_step_s", 1), "conditions.output_step_s", errors
    )

    # Species block
    product = species_block.get("product")
    if product is not None and not isinstance(product, str):
        errors.append(f"species.product: must be one species name, got {product!r}")
        product = None
    if product is not None and product == studied:
        errors.append(
            f"conditions.studied_species: {product} is the product, which always starts at 0"
        )
    if product is not None and names and product not in names:
        errors.append(f"species.product: {product} does not appear in any step")
    overall = None
    if "overall_reaction" in species_block:
        try:
            overall = parse_step(species_block["overall_reaction"], allow_ts=False)
        except StepSyntaxError as error:
            errors.append(f"species.overall_reaction {error}")
    if overall is not None:
        for name in overall.species():
            if names and name not in names:
                errors.append(
                    f"species.overall_reaction: {name} does not appear in any step"
                )
        if product is not None and overall.coefficient(product, "products") == 0:
            errors.append(
                f"species.overall_reaction: the product {product} must be on its right side"
            )
    cycles, raw_cycles = {}, species_block.get("cycles")
    if raw_cycles is not None:
        if not isinstance(raw_cycles, dict) or not all(
            isinstance(v, list) for v in raw_cycles.values()
        ):
            errors.append(
                "species.cycles: must map each cycle label to a list of intermediates"
            )
        else:
            owner = {}
            for label, members in raw_cycles.items():
                for member in members:
                    if names and member not in names:
                        errors.append(
                            f"species.cycles.{label}: {member} does not appear in any step"
                        )
                    if member in owner:
                        errors.append(
                            f"species.cycles.{label}: {member} is also in cycle {owner[member]}"
                        )
                    owner.setdefault(member, label)
                cycles[str(label)] = tuple(members)
    correction = species_block.get("vibrational_correction", "RRHO")
    if correction not in ("RRHO", "Grimme"):
        errors.append(
            f"species.vibrational_correction: must be RRHO or Grimme, got {correction!r}"
        )

    # Analyses
    ea = drc = mk = None
    raw = _analysis(analyses, "activation_energy", errors)
    if raw is not None:
        where = "analyses.activation_energy"
        initial = _concentrations(raw.get("initial_M"), f"{where}.initial_M", errors)
        temps = _range(
            raw.get("temperatures_K"),
            f"{where}.temperatures_K",
            errors,
            minimum_points=3,
        )
        total = _number(
            raw.get("simulation_time_s"), f"{where}.simulation_time_s", errors
        )
        sampling = _number(
            raw.get("sampling_time_h"), f"{where}.sampling_time_h", errors
        )
        _total_on_grid(total, output_step, where, errors)
        if None not in (total, sampling, output_step):
            _on_grid(sampling, output_step, total, f"{where}.sampling_time_h", errors)
        ea = ActivationEnergy(initial, temps or (), total, sampling)
    raw = _analysis(analyses, "degree_of_rate_control", errors)
    if raw is not None:
        where = "analyses.degree_of_rate_control"
        if "initial_M" in raw:
            initial = _concentrations(raw["initial_M"], f"{where}.initial_M", errors)
        elif ea is not None:
            initial = dict(ea.initial_M)
        else:
            initial = {}
            errors.append(
                f"{where}.initial_M: required when activation_energy is not requested"
            )
        cores = _number(
            raw.get("cores", os.cpu_count() or 1),
            f"{where}.cores",
            errors,
            integer=True,
        )
        drc = DegreeOfRateControl(
            initial,
            _number(raw.get("e_shift_kcal_mol"), f"{where}.e_shift_kcal_mol", errors),
            _number(raw.get("simulation_time_s"), f"{where}.simulation_time_s", errors),
            cores,
        )
    raw = _analysis(analyses, "microkinetics", errors)
    if raw is not None:
        where = "analyses.microkinetics"
        initial = _concentrations(raw.get("initial_M"), f"{where}.initial_M", errors)
        total = _number(
            raw.get("simulation_time_s"), f"{where}.simulation_time_s", errors
        )
        snapshot = _number(
            raw.get("catalyst_snapshot_h"), f"{where}.catalyst_snapshot_h", errors
        )
        conversion = _number(raw.get("conversion"), f"{where}.conversion", errors)
        if conversion is not None and not conversion < 1:
            errors.append(
                f"{where}.conversion: must be between 0 and 1, got {conversion}"
            )
        _total_on_grid(total, output_step, where, errors)
        if None not in (total, snapshot, output_step):
            _on_grid(
                snapshot, output_step, total, f"{where}.catalyst_snapshot_h", errors
            )
        plot = raw.get("plot_species", [])
        if not isinstance(plot, list):
            errors.append(f"{where}.plot_species: must be a list of species")
            plot = []
        for name in plot:
            if names and name not in names:
                errors.append(
                    f"{where}.plot_species: {name} does not appear in any step"
                )
        if not cycles:
            errors.append("species.cycles: required for microkinetics")
        mk = Microkinetics(initial, total, snapshot, conversion, tuple(plot))

    # Names used in initial concentrations
    for label, analysis in (
        ("activation_energy", ea),
        ("degree_of_rate_control", drc),
        ("microkinetics", mk),
    ):
        if analysis is None or (
            label == "degree_of_rate_control"
            and "initial_M" not in (analyses[label] or {})
        ):
            continue
        where = f"analyses.{label}.initial_M"
        for name in analysis.initial_M:
            if names and name not in names:
                errors.append(f"{where}.{name}: does not appear in any step")
            if name == studied:
                errors.append(
                    f"{where}.{name}: the studied species is set by conditions.studied_range_M"
                )
            if name == product:
                errors.append(f"{where}.{name}: the product always starts at 0")
        if overall is not None and analysis.initial_M:
            for name, _ in overall.reactants:
                if name != studied and name not in analysis.initial_M:
                    errors.append(
                        f"{where}: give the overall reactant {name}, or no product can form"
                    )
        members = {m for ms in cycles.values() for m in ms}
        if cycles and analysis.initial_M and not members & set(analysis.initial_M):
            errors.append(
                f"{where}: give the initial concentration of at least one catalyst intermediate"
            )

    # Energy sources
    needed = (
        sorted({temperature, *(ea.temperatures_K if ea else ())} - {None})
        if temperature
        else []
    )
    files = energies = None
    everything = names + [ts for ts in transition_states if ts not in names]
    if "files" in species_block and "energies_kcal_mol" not in species_block:
        files = (path.parent / str(species_block["files"])).resolve()
        if check_energy_sources:
            if not files.is_dir():
                errors.append(f"species.files: folder {files} not found")
            else:
                for name in everything:
                    if not (files / f"{name}.out").is_file():
                        errors.append(
                            f"no energy for {name}: {species_block['files']}/{name}.out not found"
                        )
            thermochange = os.environ.get("thermochange")
            if not thermochange:
                errors.append(
                    "$thermochange is not set: export thermochange=/path/to/thermochange"
                )
            elif not os.path.isfile(os.path.join(thermochange, THERMOCHANGE_SCRIPT)):
                errors.append(
                    f"$thermochange={thermochange} does not contain {THERMOCHANGE_SCRIPT}"
                )
    elif "energies_kcal_mol" in species_block and "files" not in species_block:
        energies = _energies(
            species_block["energies_kcal_mol"],
            "species.energies_kcal_mol",
            errors,
            needed,
            temperature,
            ea is not None,
        )
        for name in everything:
            if name not in energies:
                errors.append(
                    f"no energy for {name}: add it to species.energies_kcal_mol"
                )
        for name in energies:
            if everything and name not in everything:
                warnings.append(
                    f"warning: species.energies_kcal_mol.{name} is not used by any step"
                )

    # Warnings
    if cycles:
        members = {m for ms in cycles.values() for m in ms}
        overall_names = set(overall.species()) if overall else set()
        loose = [
            n
            for n in names
            if n not in members and n not in overall_names and n != studied
        ]
        if loose:
            warnings.append(
                f"warning: {', '.join(loose)} belong to no cycle, so the catalyst distribution does not count them"
            )

    if errors:
        raise StudyError(errors)
    return Study(
        path=path,
        name=str(doc.get("name", path.stem)),
        files=files,
        energies_kcal_mol=energies,
        product=product,
        overall=overall,
        cycles=cycles,
        vibrational_correction=correction,
        steps=tuple(steps),
        temperature_K=temperature,
        studied_species=studied,
        studied_concentrations_M=concentrations,
        output_step_s=output_step,
        activation_energy=ea,
        degree_of_rate_control=drc,
        microkinetics=mk,
        warnings=tuple(warnings),
    )
