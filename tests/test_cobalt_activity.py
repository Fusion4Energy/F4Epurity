"""
Test for Cobalt-60 activity calculation under Y1 irradiation scenario.
Regression test to ensure consistent activity calculations.
"""

import numpy as np
import pytest
import os
import json
from math import pi
import pandas as pd
from importlib.resources import as_file, files
from f4epurity.reaction_rate import calculate_reaction_rate
from f4epurity.decay_chain_calc import calculate_total_activity
from f4epurity.dose import extract_dose_factors


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
    reaction_rate_co60 = calculate_reaction_rate(sigma_eff_co60, flux_spectrum)
    reaction_rate_co60m = calculate_reaction_rate(sigma_eff_co60m, flux_spectrum)
    
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
    
    
    # Assert within tolerance
    assert co60_activity == pytest.approx(expected_activity, rel=tolerance), \
        f"Co-60 activity {co60_activity:.2f} Bq is not within 1% of expected {expected_activity:.2f} Bq"
    
    print(f"  ✓ Test PASSED (within {tolerance*100}% tolerance)")
    
    # ==================== Dose Calculation Tests ====================
    # Load dose conversion factors
    dose_matrix_file_path = files("f4epurity.resources").joinpath("F4E_dosematrix.xlsx")
    with as_file(dose_matrix_file_path) as fp:
        dose_factors_df = pd.read_excel(fp)
    
    # Extract dose conversion factor for Co-60
    dose_factor = extract_dose_factors('Co060', dose_factors_df)
    
    # Check if dose factor was found and convert to float
    if isinstance(dose_factor, str):
        # Try alternative naming conventions
        for name_variant in ['Co-60', 'co060', 'CO060', 'Co60', 'CO60']:
            dose_factor = extract_dose_factors(name_variant, dose_factors_df)
            if not isinstance(dose_factor, str):
                print(f"  Found dose factor using variant: {name_variant}")
                break
        else:
            pytest.fail(f"Could not find dose conversion factor for Co-60 in any naming format")
    
    # Ensure dose_factor is numeric
    dose_factor = float(dose_factor)
    
    # Calculate dose at source (point source model)
    # dose = dose_factor * activity * 1e6  (convert Sv to μSv)
    dose_at_source = dose_factor * co60_activity * 1e6  # μSv/h/g
    
    # Calculate dose at specific distances using 1/r^2 law for point source
    # dose_at_r = dose_at_source / (4 * pi * r^2)
    distance_1cm = 1.0  # cm
    distance_100cm = 100.0  # cm
    
    dose_at_1cm = dose_at_source / (4 * pi * distance_1cm**2)
    dose_at_100cm = dose_at_source / (4 * pi * distance_100cm**2)
    
    # Expected values computed by hand from the activity value (233.89 Bq) used as input for RadPro
    # These reference values validate the dose calculation chain from activity to dose rate
    expected_dose_1cm = 0.719009897593613  # μSv/h/g at 1 cm
    expected_dose_100cm = 7.16568455521852E-05  # μSv/h/g at 100 cm
    dose_tolerance = 0.05  # 5% tolerance
    
    
    # Assert dose values match expected within 5% tolerance
    assert dose_at_1cm == pytest.approx(expected_dose_1cm, rel=dose_tolerance), \
        f"Dose at 1 cm ({dose_at_1cm:.6e} μSv/h/g) is not within 5% of expected ({expected_dose_1cm:.6e} μSv/h/g)"
    
    assert dose_at_100cm == pytest.approx(expected_dose_100cm, rel=dose_tolerance), \
        f"Dose at 100 cm ({dose_at_100cm:.6e} μSv/h/g) is not within 5% of expected ({expected_dose_100cm:.6e} μSv/h/g)"
    
    print(f"  ✓ Dose tests PASSED (within {dose_tolerance*100}% tolerance)")

