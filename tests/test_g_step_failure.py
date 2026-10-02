"""A failed G step raises with the script's output and leaves no half-written results behind.

Uses a fake get_G_compounds.sh that writes the compound file and then fails, as the real
script does when a species in reactions.csv has no .out file:
    python tests/test_g_step_failure.py   (or: pytest tests)
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

FAKE_SCRIPT = """\
mkdir -p G_values_of_compounds
echo "Compounds,Gibbs Free Energies" > "G_values_of_compounds/G_values_at_$1K_$(printf '%.5e' $2)atm.csv"
echo "KeyError: 'TS3_0L'" >&2
"""


def test_failed_g_step_raises_and_cleans_up():
    os.chdir(tempfile.mkdtemp())
    with open("get_G_compounds.sh", "w") as script:
        script.write(FAKE_SCRIPT)
    from auxiliary_functions import AuxiliaryFunctions

    try:
        AuxiliaryFunctions.calculate_G_values(350.0, 28.7)
    except RuntimeError as error:
        assert "TS3_0L" in str(error)
    else:
        raise AssertionError("expected RuntimeError")
    assert not os.listdir("G_values_of_compounds")  # no stale file for the next run


if __name__ == "__main__":
    test_failed_g_step_raises_and_cleans_up()
    print("ok")
