import numpy as np
import os
import pytest
from unittest.mock import patch
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

# Test extract_xs 
def test_extract_xs():
    parent = "co59"
    product = "co060"
    element = "Co"
    # From TENDL17 5-g data
    expected_output = [8.6890E+00, 2.1642E+00, 4.3089E-03, 1.3517E-03, 1.3366E-04]
    assert extract_xs(parent, product, element) == expected_output

def test_perform_collapse():
    # Dummy data
    xs_values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    flux = np.array([1.0, 2.0, 3.0, 4.0, 5.0])

    # Expected result
    total_flux = np.sum(flux)
    expected_result = np.sum(xs_values * flux) / total_flux

    # Call the function
    result = perform_collapse(xs_values, flux)

    np.testing.assert_almost_equal(result, expected_result, decimal=5)

def mock_get_flux_from_vtk(*args, **kwargs):
    # Return a dummy flux spectrum
    return np.array([[1.0, 2.0, 3.0, 4.0, 5.0], [2.0, 3.0, 4.0, 5.0, 6.0]])

@patch('f4epurity.collapse.get_flux_from_vtk', new=mock_get_flux_from_vtk)
def test_collapse_flux():
    # Dummy data
    xs_values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])

    # Expected result
    flux_spectrum = mock_get_flux_from_vtk()
    expected_result = np.array([perform_collapse(xs_values, flux) for flux in flux_spectrum])

    # Call the function with the mock_get_flux_from_vtk function
    result, _ = collapse_flux(xs_values, 'dummy_path', 0, 0, 0)

    np.testing.assert_almost_equal(result, expected_result, decimal=5)