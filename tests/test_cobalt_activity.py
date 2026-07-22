"""
Test for Cobalt-60 activity calculation under Y1 irradiation scenario.
Regression test to ensure consistent activity calculations.
"""

import numpy as np
import pytest
import os
import json
from f4epurity.reaction_rate import calculate_reaction_rate
from f4epurity.decay_chain_calc import calculate_total_activity


def test_cobalt_activity_y1_scenario():
    """
    Test that 0.05% Co impurity under Y1 irradiation (1 year) with 1e6s decay
    produces expected Co-60 activity within 1% tolerance.
    
    This test runs the actual F4Epurity calculation pipeline and verifies
    that the logged activity matches the expected value.
    
    Reference conditions:
    - delta_impurity: 0.0005 (0.05%)
    - Flux spectrum: [4987173.5191, 44082393.971, 8885176.6378, 1600249.3564, 96519.789311] n/cm²/s
    - Irradiation: Y1 (1 year)
    - Decay time: 1e6 s
    - Expected Co-60 activity: 233.89 Bq (±1%, i.e., ±2.34 Bq)
    """
    # Load decay data from JSON file
    decay_data_path = os.path.join(os.path.dirname(__file__), "data", "Decay2020.json")
    with open(decay_data_path, "r", encoding="utf-8") as f:
        decay_data = json.load(f)
    
    # Test parameters
    delta_impurity = 0.0005  # 0.05%
    flux_spectrum = np.array([[4987173.5191, 44082393.971, 8885176.6378, 
                               1600249.3564, 96519.789311]])
    irrad_scenario = "Y1"
    decay_time = 1e6  # seconds
    
    # Cross sections for Co-59 -> Co-60 from log
    sigma_eff_co60 = 2.326465  # barns
    sigma_eff_co60m = 3.873241  # barns
    
    # Calculate reaction rates using F4Epurity function
    reaction_rate_co60 = calculate_reaction_rate(delta_impurity, sigma_eff_co60, flux_spectrum)
    reaction_rate_co60m = calculate_reaction_rate(delta_impurity, sigma_eff_co60m, flux_spectrum)
    
    # Calculate number of atoms for Co-59
    # For Co-59: natural abundance ~100%, atomic mass ~59
    # atoms = (delta_impurity * N_A) / atomic_mass = (0.0005 * 6.022e23) / 59
    N_A = 6.02214076e23  # Avogadro's number
    atomic_mass_co59 = 58.933194  # u
    number_of_atoms = (delta_impurity * N_A) / atomic_mass_co59
    
    # Wrap reaction rates in arrays - decay_chain_calc expects array values for mesh cells
    if np.isscalar(reaction_rate_co60):
        rr_co60_array = [float(reaction_rate_co60)]
        rr_co60m_array = [float(reaction_rate_co60m)]
    else:
        rr_co60_array = reaction_rate_co60.tolist()
        rr_co60m_array = reaction_rate_co60m.tolist()

    # Build nuclide dictionary structure matching F4Epurity format
    nuclide_dict = {
        'co59': {
            'atoms': number_of_atoms,
            'reactions': {
                'co060': rr_co60_array,
                'co060m': rr_co60m_array
            }
        }
    }
    
    # Calculate activities using the actual F4Epurity function
    # This will also produce the logged output that the user wants to verify
    activities = calculate_total_activity(nuclide_dict, irrad_scenario, decay_time, decay_data)
    
    # Extract Co-60 activity from results
    co60_activity = activities.get('Co060', [0])[0]  # Note: returned keys are capitalized
    if isinstance(co60_activity, np.ndarray):
        co60_activity = float(co60_activity.flat[0])
    else:
        co60_activity = float(co60_activity)
    
    # Expected value from reference calculation
    expected_activity = 233.89  # Bq
    tolerance = 0.01  # 1%
    
    # Print result for user visibility
    print(f"\n✓ F4Epurity calculation complete:")
    print(f"  Co-60 activity: {co60_activity:.2f} Bq")
    print(f"  Expected: {expected_activity:.2f} ± {expected_activity*tolerance:.2f} Bq")
    
    # Assert within tolerance
    assert co60_activity == pytest.approx(expected_activity, rel=tolerance), \
        f"Co-60 activity {co60_activity:.2f} Bq is not within 1% of expected {expected_activity:.2f} Bq"
    
    print(f"  ✓ Test PASSED (within {tolerance*100}% tolerance)")

