import numpy as np

from swim.isotopes import (
    R17O_VSMOW,
    R18O_VSMOW,
    RD_VSMOW,
    delta_to_ratio,
    linear_deuterium_excess,
    log_delta,
    logarithmic_deuterium_excess,
    oxygen_17_excess,
    ratio_to_delta,
)


def test_vsmow_constants_match_legacy_matlab() -> None:
    assert R18O_VSMOW == 0.00200520
    assert RD_VSMOW == 0.00015576
    assert R17O_VSMOW == 0.0003799


def test_delta_ratio_conversion_round_trip() -> None:
    delta = np.array([-500.0, -50.0, 0.0, 25.0])

    ratio = delta_to_ratio(delta, R18O_VSMOW)

    np.testing.assert_allclose(ratio_to_delta(ratio, R18O_VSMOW), delta)


def test_excess_definitions() -> None:
    delta_18o = np.array([-40.0, -20.0])
    delta_d = np.array([-310.0, -150.0])
    delta_17o = np.array([-21.0, -10.5])

    np.testing.assert_allclose(
        linear_deuterium_excess(delta_d, delta_18o), np.array([10.0, 10.0])
    )
    np.testing.assert_allclose(
        logarithmic_deuterium_excess(delta_d, delta_18o),
        log_delta(delta_d)
        - (-0.0285 * log_delta(delta_18o) ** 2 + 8.47 * log_delta(delta_18o)),
    )
    np.testing.assert_allclose(
        oxygen_17_excess(delta_17o, delta_18o),
        1.0e6
        * (
            np.log1p(delta_17o / 1000.0)
            - 0.528 * np.log1p(delta_18o / 1000.0)
        ),
    )
