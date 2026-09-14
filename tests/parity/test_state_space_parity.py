from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from swim.model import forward_state_space

STATE_SPACE_PATH = (
    Path(__file__).parents[1]
    / "fixtures"
    / "matlab"
    / "port_baseline_v1"
    / "allan_hills"
    / "state_space.mat"
)


@pytest.mark.parity
def test_full_state_space_matches_matlab() -> None:
    expected = loadmat(STATE_SPACE_PATH, simplify_cells=True)
    actual = forward_state_space(
        np.asarray(expected["T_source"]),
        np.asarray(expected["T_site"]),
        evaporation_version="2021",
    )

    comparisons = {
        "d18O_site": actual.delta_18o,
        "dD_site": actual.delta_d,
        "d18Oln_site": actual.delta_18o_log,
        "dDln_site": actual.delta_d_log,
        "dxs_site": actual.deuterium_excess,
        "d17O_xs_site": actual.oxygen_17_excess_per_meg,
        "dlnU_site": actual.logarithmic_deuterium_excess,
        "r_s_site": actual.saturated_mixing_ratio,
        "P_site": actual.pressure_kpa,
    }
    for matlab_name, python_values in comparisons.items():
        absolute_tolerance = 1.0e-7 if matlab_name == "d17O_xs_site" else 1.0e-8
        np.testing.assert_allclose(
            python_values,
            expected[matlab_name],
            rtol=1.0e-9,
            atol=absolute_tolerance,
            equal_nan=True,
            err_msg=matlab_name,
        )
