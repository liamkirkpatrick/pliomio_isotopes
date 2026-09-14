import numpy as np

from swim.reconstruction import (
    seawater_correct_isotopes,
    surface_temperature_from_condensation,
)


def test_surface_condensation_temperature_conversion() -> None:
    condensation = np.array([-40.0, -20.0, 0.0])
    surface = surface_temperature_from_condensation(condensation)

    np.testing.assert_allclose(0.69 * surface - 8.2, condensation)


def test_seawater_correction_is_identity_at_initial_age() -> None:
    delta_18o = np.array([-40.0])
    delta_d = np.array([-310.0])

    corrected = seawater_correct_isotopes(delta_18o, delta_d, [0.0], 0.0)

    np.testing.assert_allclose(corrected.delta_18o_permil, delta_18o)
    np.testing.assert_allclose(corrected.delta_d_permil, delta_d)
