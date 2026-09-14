import numpy as np

from swim.fractionation import (
    O17_ICE_EXPONENT_DISTILLATION,
    O17_ICE_EXPONENT_EVAPORATION,
    O17_LIQUID_EXPONENT,
    equilibrium_fractionation_factors,
    kinetic_condensation_factors,
    mixed_phase_effective_fractionation,
    transport_diffusivity_ratios,
)


def test_equilibrium_fractionation_reference_values_at_zero_celsius() -> None:
    factors = equilibrium_fractionation_factors(
        0.0, oxygen_17_ice_exponent=O17_ICE_EXPONENT_DISTILLATION
    )

    np.testing.assert_allclose(factors.deuterium_liquid, 1.1123216522954846)
    np.testing.assert_allclose(factors.deuterium_ice, 1.1335720204365185)
    np.testing.assert_allclose(factors.oxygen_18_liquid, 1.011718982796772)
    np.testing.assert_allclose(factors.oxygen_18_ice, 1.015233348238943)
    np.testing.assert_allclose(factors.oxygen_17_liquid, 1.0061823304352444)
    np.testing.assert_allclose(factors.oxygen_17_ice, 1.0080602260875429)


def test_oxygen_17_ice_exponent_is_explicitly_selectable() -> None:
    distillation = equilibrium_fractionation_factors(
        -30.0, oxygen_17_ice_exponent=O17_ICE_EXPONENT_DISTILLATION
    )
    evaporation = equilibrium_fractionation_factors(
        -30.0, oxygen_17_ice_exponent=O17_ICE_EXPONENT_EVAPORATION
    )

    assert O17_LIQUID_EXPONENT == 0.529
    assert distillation.oxygen_17_liquid == evaporation.oxygen_17_liquid
    assert distillation.oxygen_17_ice > evaporation.oxygen_17_ice


def test_transport_diffusivity_ratios_are_heavy_over_light() -> None:
    ratios = transport_diffusivity_ratios(np.array([10.0, -10.0, -30.0]))

    assert np.all(np.asarray(ratios.hdo_over_h2o) < 1.0)
    assert np.all(np.asarray(ratios.h217o_over_h216o) < 1.0)
    assert np.all(np.asarray(ratios.h218o_over_h216o) < 1.0)


def test_kinetic_factors_equal_one_at_ice_saturation() -> None:
    temperature_c = np.array([0.0, -20.0])
    equilibrium = equilibrium_fractionation_factors(
        temperature_c,
        oxygen_17_ice_exponent=O17_ICE_EXPONENT_DISTILLATION,
    )
    diffusivity = transport_diffusivity_ratios(temperature_c)

    kinetic = kinetic_condensation_factors(
        np.ones(temperature_c.shape), equilibrium, diffusivity
    )

    np.testing.assert_array_equal(kinetic.deuterium_ice, 1.0)
    np.testing.assert_array_equal(kinetic.oxygen_18_ice, 1.0)
    np.testing.assert_array_equal(kinetic.oxygen_17_ice, 1.0)


def test_effective_fractionation_reduces_to_pure_phases() -> None:
    temperature_c = np.array([-20.0, -20.0])
    equilibrium = equilibrium_fractionation_factors(
        temperature_c,
        oxygen_17_ice_exponent=O17_ICE_EXPONENT_DISTILLATION,
    )
    kinetic = kinetic_condensation_factors(
        np.array([1.1, 1.1]),
        equilibrium,
        transport_diffusivity_ratios(temperature_c),
    )
    effective = mixed_phase_effective_fractionation(
        equilibrium,
        kinetic,
        fraction_ice=np.array([0.0, 1.0]),
        fraction_liquid=np.array([1.0, 0.0]),
    )

    np.testing.assert_allclose(
        effective.deuterium,
        np.array(
            [
                np.asarray(equilibrium.deuterium_liquid)[0],
                np.asarray(equilibrium.deuterium_ice)[1]
                * np.asarray(kinetic.deuterium_ice)[1],
            ]
        ),
    )
