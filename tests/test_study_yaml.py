"""Strict YAML loading for study files (spec section 3.6).

python tests/test_study_yaml.py
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from study_yaml import StudyError, load_yaml


def write(text):
    path = os.path.join(tempfile.mkdtemp(), "study.yaml")
    with open(path, "w") as f:
        f.write(text)
    return path


def test_exponent_numbers_without_a_dot_are_numbers():
    doc = load_yaml(write("a: 1e-10\nb: 5e-4\nc: 1.0e-10\nd: 19\ne: 0.1\n"))
    assert doc == {"a": 1e-10, "b": 5e-4, "c": 1e-10, "d": 19, "e": 0.1}
    assert isinstance(doc["d"], int)


def test_species_names_like_no_stay_text():
    doc = load_yaml(
        write("initial_M: {NO: 0.05, ON: 1, Y: 2, yes: 3, off: 4}\nflag: true\n")
    )
    assert doc["initial_M"] == {"NO": 0.05, "ON": 1, "Y": 2, "yes": 3, "off": 4}
    assert doc["flag"] is True


def test_temperature_keys_are_numbers():
    doc = load_yaml(write("TS1: {325: 11.5, 337.5: 11.6}\n"))
    assert doc["TS1"] == {325: 11.5, 337.5: 11.6}


def test_duplicate_key_is_an_error_with_both_lines():
    try:
        load_yaml(
            write("species:\n  energies_kcal_mol:\n    TS1: 1\n    A: 2\n    TS1: 3\n")
        )
    except StudyError as error:
        assert "TS1 appears twice (lines 3 and 5)" in str(error), str(error)
    else:
        raise AssertionError("expected StudyError")


def test_syntax_error_names_the_line():
    try:
        load_yaml(write("a: 1\nb: [1, 2\nc: 3\n"))
    except StudyError as error:
        assert "line" in str(error) and "study.yaml" in str(error), str(error)
    else:
        raise AssertionError("expected StudyError")


def test_top_level_must_be_a_mapping():
    try:
        load_yaml(write("- a\n- b\n"))
    except StudyError as error:
        assert "must be a mapping" in str(error)
    else:
        raise AssertionError("expected StudyError")


if __name__ == "__main__":
    for name, test in list(globals().items()):
        if name.startswith("test_"):
            test()
    print("ok")
