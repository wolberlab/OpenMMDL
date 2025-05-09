import pytest
import pandas as pd
import numpy as np
import re
import os
import xml.etree.ElementTree as ET
from unittest.mock import patch, mock_open
from io import StringIO

from openmmdl.openmmdl_analysis.visualization.pharmacophore import PharmacophoreGenerator


@pytest.fixture
def sample_df():
    """Create a sample dataframe for testing."""
    data = {
        "INTERACTION": ["hbond", "hbond", "hydrophobic", "saltbridge", "saltbridge", "pistacking", "metal"],
        "LIGCOO": ["(1.234, 2.345, 3.456)", "(4.567, 5.678, 6.789)", "(7.890, 8.901, 9.012)", 
                   "(10.123, 11.234, 12.345)", "(13.456, 14.567, 15.678)", "(16.789, 17.890, 18.901)", "0"],
        "PROTCOO": ["(21.234, 22.345, 23.456)", "(24.567, 25.678, 26.789)", "(27.890, 28.901, 29.012)", 
                    "(30.123, 31.234, 32.345)", "(33.456, 34.567, 35.678)", "(36.789, 37.890, 38.901)", "0"],
        "TARGETCOO": ["0", "0", "0", "0", "0", "0", "(40.123, 41.234, 42.345)"],
        "PROTISDON": ["True", "False", "False", "False", "False", "False", "False"],
        "PROTISPOS": ["False", "False", "False", "True", "False", "False", "False"],
        "Acceptor_hbond": [0, 1, 0, 0, 0, 0, 0],
        "Donor_hbond": [1, 0, 0, 0, 0, 0, 0],
        "hydrophobic": [0, 0, 1, 0, 0, 0, 0],
        "PI_saltbridge": [0, 0, 0, 1, 0, 0, 0],
        "NI_saltbridge": [0, 0, 0, 0, 1, 0, 0],
        "pistacking": [0, 0, 0, 0, 0, 1, 0],
        "metal": [0, 0, 0, 0, 0, 0, 1]
    }
    return pd.DataFrame(data)


def test_init(sample_df):
    """Test the initialization of the PharmacophoreGenerator class."""
    generator = PharmacophoreGenerator(sample_df, "test_ligand")
    
    assert generator.df_all.equals(sample_df)
    assert generator.ligand_name == "test_ligand"
    assert generator.comlex_name == "test_ligand_complex"
    assert isinstance(generator.coord_pattern, re.Pattern)
    assert isinstance(generator.clouds, dict)


def test_generate_clouds(sample_df):
    """Test the _generate_clouds method."""
    generator = PharmacophoreGenerator(sample_df, "test_ligand")
    clouds = generator.clouds
    
    # Check that expected interaction types exist in the clouds
    assert "hydrophobic" in clouds
    assert "acceptor" in clouds
    assert "donor" in clouds
    assert "metal" in clouds
    
    # Check that each cloud has the expected properties
    for interaction, cloud in clouds.items():
        assert "coordinates" in cloud
        assert "color" in cloud
        assert "radius" in cloud
        
        # Check color is a list of 3 values
        assert isinstance(cloud["color"], list)
        assert len(cloud["color"]) == 3
        
        # Check radius is a float
        assert isinstance(cloud["radius"], float)
        
        # Check coordinates is a list
        assert isinstance(cloud["coordinates"], list)


def test_to_dict(sample_df):
    """Test the to_dict method."""
    generator = PharmacophoreGenerator(sample_df, "test_ligand")
    clouds_dict = generator.to_dict()
    
    # Should return the same as generator.clouds
    assert clouds_dict == generator.clouds


def test_generate_pharmacophore_centers(sample_df):
    """Test the _generate_pharmacophore_centers method."""
    generator = PharmacophoreGenerator(sample_df, "test_ligand")
    
    # Test with hydrophobic interactions
    interactions = ["hydrophobic"]
    centers = generator._generate_pharmacophore_centers(interactions)
    
    assert "hydrophobic" in centers
    assert len(centers["hydrophobic"]) == 3  # x, y, z coordinates
    assert all(isinstance(coord, float) for coord in centers["hydrophobic"])


def test_generate_pharmacophore_vectors(sample_df):
    """Test the _generate_pharmacophore_vectors method."""
    generator = PharmacophoreGenerator(sample_df, "test_ligand")
    
    # Test with donor interactions
    interactions = ["Donor_hbond"]
    vectors = generator._generate_pharmacophore_vectors(interactions)
    
    assert "Donor_hbond" in vectors
    assert len(vectors["Donor_hbond"]) == 2  # ligand and protein coords
    assert len(vectors["Donor_hbond"][0]) == 3  # x, y, z ligand
    assert len(vectors["Donor_hbond"][1]) == 3  # x, y, z protein
    assert all(isinstance(coord, float) for coord in vectors["Donor_hbond"][0])
    assert all(isinstance(coord, float) for coord in vectors["Donor_hbond"][1])


@patch("builtins.open", new_callable=mock_open)
@patch("xml.etree.ElementTree.ElementTree.write")
def test_generate_md_pharmacophore_cloudcenters(mock_write, mock_open, sample_df):
    """Test the generate_md_pharmacophore_cloudcenters method."""
    generator = PharmacophoreGenerator(sample_df, "test_ligand")
    
    generator.generate_md_pharmacophore_cloudcenters("test_output")
    
    # Check that ElementTree.write was called with the correct filename
    mock_write.assert_called_once_with("test_output.pml", encoding="UTF-8", xml_declaration=True)



@patch("builtins.open", new_callable=mock_open)
@patch("xml.etree.ElementTree.ElementTree.write")
def test_generate_bindingmode_pharmacophore(mock_write, mock_open, sample_df):
    """Test the generate_bindingmode_pharmacophore method."""
    generator = PharmacophoreGenerator(sample_df, "test_ligand")
    
    # Mock binding mode dictionary
    binding_mode = {
        "Acceptor_hbond": {
            "LIGCOO": [[4.567, 5.678, 6.789]],
            "PROTCOO": [[24.567, 25.678, 26.789]]
        },
        "hydrophobic": {
            "LIGCOO": [[7.890, 8.901, 9.012]],
            "PROTCOO": [[27.890, 28.901, 29.012]]
        }
    }
    
    generator.generate_bindingmode_pharmacophore(binding_mode, "test_binding_mode")
    
    # Check that the file was opened for writing
    expected_path = "./Binding_Modes_Markov_States/test_binding_mode.pml"
    mock_write.assert_called_once()


def test_generate_pharmacophore_centers_all_points(sample_df):
    """Test the _generate_pharmacophore_centers_all_points method."""
    generator = PharmacophoreGenerator(sample_df, "test_ligand")
    
    # Test with hydrophobic interactions
    interactions = ["hydrophobic"]
    points = generator._generate_pharmacophore_centers_all_points(interactions)
    
    assert "hydrophobic" in points
    assert len(points["hydrophobic"]) == 1  # One point for the single hydrophobic interaction
    assert len(points["hydrophobic"][0]) == 3  # x, y, z coordinates
    assert all(isinstance(coord, float) for coord in points["hydrophobic"][0])


@patch("builtins.open", new_callable=mock_open)
@patch("xml.etree.ElementTree.ElementTree.write")
def test_generate_point_cloud_pml(mock_write, mock_open, sample_df):
    """Test the generate_point_cloud_pml method."""
    generator = PharmacophoreGenerator(sample_df, "test_ligand")
    
    generator.generate_point_cloud_pml("test_cloud")
    
    # Check that the file was opened for writing
    mock_write.assert_called_once_with("test_cloud.pml", encoding="UTF-8", xml_declaration=True)


def test_format_clouds():
    """Test the _format_clouds method."""
    generator = PharmacophoreGenerator(pd.DataFrame(), "test_ligand")
    
    interaction_coords = {
        "hydrophobic": [[1.0, 2.0, 3.0]],
        "acceptor": [[4.0, 5.0, 6.0]],
    }
    
    formatted_clouds = generator._format_clouds(interaction_coords)
    
    assert "hydrophobic" in formatted_clouds
    assert "acceptor" in formatted_clouds
    
    assert formatted_clouds["hydrophobic"]["coordinates"] == [[1.0, 2.0, 3.0]]
    assert formatted_clouds["acceptor"]["coordinates"] == [[4.0, 5.0, 6.0]]
    
    assert formatted_clouds["hydrophobic"]["color"] == [1.0, 1.0, 0.0]
    assert formatted_clouds["acceptor"]["color"] == [1.0, 0.0, 0.0]
    
    assert formatted_clouds["hydrophobic"]["radius"] == 0.1
    assert formatted_clouds["acceptor"]["radius"] == 0.1


def test_coordinate_parsing():
    """Test the coordinate pattern regex."""
    generator = PharmacophoreGenerator(pd.DataFrame(), "test_ligand")
    
    coord_str = "(1.234, 2.345, 3.456)"
    match = generator.coord_pattern.match(coord_str)
    
    assert match is not None
    x, y, z = map(float, match.groups())
    assert x == 1.234
    assert y == 2.345
    assert z == 3.456


def test_empty_dataframe():
    """Test initialization with an empty dataframe."""
    generator = PharmacophoreGenerator(pd.DataFrame(), "test_ligand")
    
    # Should have empty clouds but not fail
    clouds = generator.clouds
    for interaction in clouds:
        assert len(clouds[interaction]["coordinates"]) == 0
