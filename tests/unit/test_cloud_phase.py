import numpy as np
import pytest

from swim.cloud_phase import cloud_phase_fractions


def test_active_adjusted_cloud_phase_reference_values() -> None:
    temperature_c = np.array([10.0, 0.0, -10.0, -20.0, -30.0])

    fraction_ice, fraction_liquid = cloud_phase_fractions(
        temperature_c, method="adj"
    )

    np.testing.assert_allclose(
        fraction_liquid,
        np.array(
            [
                0.9999999999990887,
                0.9953674292821418,
                0.9913026370266288,
                0.942675824101131,
                0.36818758226390347,
            ]
        ),
        rtol=1.0e-14,
        atol=1.0e-15,
    )
    np.testing.assert_allclose(fraction_ice + fraction_liquid, 1.0)


def test_mid_remains_the_standalone_legacy_default() -> None:
    default = cloud_phase_fractions(np.array([-20.0, -30.0]))
    explicit = cloud_phase_fractions(np.array([-20.0, -30.0]), method="mid")

    np.testing.assert_array_equal(default[0], explicit[0])
    np.testing.assert_array_equal(default[1], explicit[1])


def test_unknown_cloud_phase_method_is_rejected() -> None:
    with pytest.raises(ValueError, match="Unknown cloud-phase method"):
        cloud_phase_fractions(0.0, method="unknown")  # type: ignore[arg-type]
