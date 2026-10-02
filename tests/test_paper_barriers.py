"""Gibbs barriers at 350 K must match Table S3 of the paper's Supporting Information.

Abdullayev et al., ACS Catal. 2025, 15, 4739 (doi:10.1021/acscatal.5c00348). Builds the barrier
table of examples/hydroformylation/study.yaml with the real thermochange, so it needs $thermochange
exported; it is skipped otherwise.
    python tests/test_paper_barriers.py   (or: pytest tests)
"""

import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)

# Rx, Gdir, Ginv in kcal/mol at 350.0 K and 1.0 M (SI Table S3, rounded to 0.1)
TABLE_S3 = [
    ("I1_0L = I2_0L + CO", 11.5, 4.0),
    ("I2_0L + ete = I3_0L", 4.0, 4.1),
    ("I3_0L = I4_0L", 11.5, 18.2),
    ("I4_0L + CO = I5_0L", 4.0, 9.1),
    ("I5_0L = I6_0L", 15.1, 20.3),
    ("I6_0L + CO = I7_0L", 4.0, 7.5),
    ("I6_0L + H2 = I8_0L", 15.8, 6.8),
    ("I8_0L = I9_0L", 11.4, 25.8),
    ("I9_0L = prod + I2_0L", 4.0, 7.5),
    ("I2_0L + PMe3 = I1_1L", 4.0, 19.9),
    ("I1_1L = I2c_1L + CO", 10.8, 4.0),
    ("I2c_1L + ete = I3_1L", 6.8, 4.0),
    ("I3_1L = I4_1L", 13.8, 21.9),
    ("I4_1L + CO = I5_1L", 4.0, 6.6),
    ("I5_1L = I6_1L", 12.7, 21.1),
    ("I6_1L + CO = I7_1L", 4.0, 7.4),
    ("I6_1L + H2 = I8_1L", 13.3, 7.5),
    ("I8_1L = I9t_1L", 10.4, 18.4),
    ("I8_1L = I9c_1L", 13.0, 21.7),
    ("I9t_1L = I2t_1L + prod", 4.0, 8.4),
    ("I9c_1L = prod + I2c_1L", 4.0, 10.6),
    ("I2t_1L + CO = I1_1L", 4.0, 13.7),
]


def test_barriers_match_paper():
    if not os.environ.get("thermochange"):
        print("skipped: export thermochange=/path/to/thermochange to run")
        return
    from study import load_study
    from thermochemistry import barrier_table, gibbs_energies

    study = load_study(os.path.join(ROOT, "examples", "hydroformylation", "study.yaml"))
    df = barrier_table(study, gibbs_energies(study, 350.0))
    assert len(df) == len(TABLE_S3)
    for (rx, gdir, ginv), (_, row) in zip(TABLE_S3, df.iterrows()):
        assert " ".join(row["Rx"].split()) == rx
        assert round(row["Gdir"], 1) == gdir, (rx, row["Gdir"], gdir)
        assert round(row["Ginv"], 1) == ginv, (rx, row["Ginv"], ginv)


if __name__ == "__main__":
    test_barriers_match_paper()
    print("ok")
