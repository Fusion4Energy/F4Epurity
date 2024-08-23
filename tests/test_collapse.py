import numpy as np
import os
import pytest
import pyvista as pv
from f4epurity.collapse import extract_flux_values, get_flux_from_vtk, extract_xs, perform_collapse, collapse_flux

# Test extract_flux_values function
def test_extract_flux_values():
    # Create a dummy grid for testing
    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])
    cells = np.array([[8, 0, 1, 2, 3, 4, 5, 6, 7], [8, 0, 1, 2, 3, 4, 5, 6, 7], [8, 0, 1, 2, 3, 4, 5, 6, 7]])
    cell_types = np.array([12, 12, 12])  
    grid = pv.UnstructuredGrid(cells, cell_types, points)
    
    # Assigning dummy cell data
    grid.cell_data['ValueBin-000'] = np.array([1, 2, 3])
    grid.cell_data['ValueBin-001'] = np.array([4, 5, 6])
    grid.cell_data['ValueBin-002'] = np.array([7, 8, 9])
    grid.cell_data['ValueBin-003'] = np.array([7, 8, 9])
    grid.cell_data['ValueBin-004'] = np.array([7, 8, 9])
    
    cell_ids = [0, 1, 2]
    expected_output = np.array([[1, 4, 7, 7, 7], [2, 5, 8, 8, 8], [3, 6, 9, 9, 9]])
    assert np.array_equal(extract_flux_values(grid, cell_ids), expected_output)

# Test get_flux_from_vtk 
def test_get_flux_from_vtk():
    filepath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'result_NeutronFlux_Plasma-N_n.s-1.cm-2.vtr')
    x1, y1, z1 = -245, 3451, 0
    x2, y2, z2 = 500, 4200, 50  
    assert isinstance(get_flux_from_vtk(filepath, x1, y1, z1), np.ndarray)
    assert isinstance(get_flux_from_vtk(filepath, x1, y1, z1, x2, y2, z2), np.ndarray)

# Test extract_xs 
def test_extract_xs():
    parent = "co59"
    product = "co060"
    element = "Co"
    # From TENDL17 5-g data
    expected_output = [8.6890E+00, 2.1642E+00, 4.3089E-03, 1.3517E-03, 1.3366E-04]
    assert extract_xs(parent, product, element) == expected_output

# Test perform_collapse
def test_perform_collapse():
    xs_values = [2, 4, 9] 
    flux = np.array([3, 5, 2])  
    expected_output = 4.4
    assert perform_collapse(xs_values, flux) == expected_output

# Test collapse flux
def test_collapse_flux():
    # Dummy xs data
    xs_values = [1.0, 2.0, 3.0, 4.0, 5.0]
    filepath_flux = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'result_NeutronFlux_Plasma-N_n.s-1.cm-2.vtr') 
    x1, y1, z1 = -245, 3451, 0
    # Span two voxels with the line source
    x2, y2, z2 = -245, 3551, 0  

    # Expected result calculated with visit and dummy xs values
    expected_result_point = np.array([2.73654])
    expected_result_line = np.array([2.73654, 2.70801])

    # Extract sigma_eff from the returned tuple
    sigma_eff_1, _ = collapse_flux(xs_values, filepath_flux, x1, y1, z1)
    sigma_eff_2, _ = collapse_flux(xs_values, filepath_flux, x1, y1, z1, x2, y2, z2)

    # Assert that the results are almost equal to the expected results to the 5th decimal place
    np.testing.assert_array_almost_equal(sigma_eff_1, expected_result_point, decimal=5)
    np.testing.assert_array_almost_equal(sigma_eff_2, expected_result_line, decimal=5)