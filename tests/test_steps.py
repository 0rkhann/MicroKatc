"""Parsing of step and overall-reaction strings (spec section 3.3).

python tests/test_steps.py
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from steps import StepSyntaxError, parse_step


def test_simple_barrierless_step():
    step = parse_step("I1_0L <=> I2_0L + CO")
    assert step.reactants == (("I1_0L", 1),)
    assert step.products == (("I2_0L", 1), ("CO", 1))
    assert step.ts is None
    assert step.copasi_equation() == "I1_0L = I2_0L + CO"


def test_step_with_transition_state():
    step = parse_step("I6_0L + H2 <=> I8_0L  via TS3_0L")
    assert step.ts == "TS3_0L"
    assert step.species() == ["I6_0L", "H2", "I8_0L"]


def test_separators_are_synonyms():
    for text in ("A + B <=> C", "A + B = C", "A + B ⇌ C", "A+B<=>C"):
        assert parse_step(text).copasi_equation() == "A + B = C", text


def test_coefficients_are_expanded():
    for text in ("2 A <=> A2 via TS2", "2*A <=> A2 via TS2", "A + A <=> A2 via TS2"):
        step = parse_step(text)
        assert step.coefficient("A", "reactants") == 2, text
        assert step.copasi_equation() == "A + A = A2", text


def test_errors():
    cases = {
        "I4_0L + CO I5_0L": "no <=> between the two sides",
        "A <=> B <=> C": "more than one <=>",
        "A <=> ": "the right side is empty",
        "1.5 A <=> B": "coefficients must be positive integers",
        "0 A <=> B": "coefficients must be positive integers",
        "A + <=> B": "empty species name",
        "A <=> B via": "via must be followed by one transition-state name",
        "A <=> B via TS1 TS2": "via must be followed by one transition-state name",
        "A b <=> B": "names cannot contain spaces",
        "1A <=> B": "names must start with a letter",
        "C2,x <=> B": "names cannot contain , or /",
        "C2/x <=> B": "names cannot contain , or /",
    }
    for text, message in cases.items():
        try:
            parse_step(text)
        except StepSyntaxError as error:
            assert message in str(error), (text, str(error))
        else:
            raise AssertionError(f"no error for {text!r}")


def test_overall_reaction_rejects_via():
    try:
        parse_step("ete + CO + H2 <=> prod via TS1", allow_ts=False)
    except StepSyntaxError as error:
        assert "via is not allowed here" in str(error)
    else:
        raise AssertionError("expected StepSyntaxError")


def test_species_name_starting_with_via():
    step = parse_step("A + via1 <=> B via TS")
    assert step.reactants == (("A", 1), ("via1", 1)) and step.ts == "TS", step


if __name__ == "__main__":
    for name, test in list(globals().items()):
        if name.startswith("test_"):
            test()
    print("ok")
