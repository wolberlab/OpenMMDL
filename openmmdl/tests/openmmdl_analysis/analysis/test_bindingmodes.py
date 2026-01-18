import pandas as pd
import pytest

from openmmdl.openmmdl_analysis.analysis.bindingmodes import BindingModeProcesser


@pytest.fixture
def sample_interaction_list():
    """Minimal interaction_list DataFrame for BindingModeProcesser."""
    return pd.DataFrame(
        [
            # Two hbonds that share the same ligand atom index (ACCEPTORIDX=5)
            {
                "FRAME": 1,
                "INTERACTION": "hbond",
                "Prot_partner": "12GLYA",
                "PROTISDON": True,
                "ACCEPTORIDX": 5,
                "DONORIDX": 0,
            },
            {
                "FRAME": 2,
                "INTERACTION": "hbond",
                "Prot_partner": "13SERA",
                "PROTISDON": True,
                "ACCEPTORIDX": 5,
                "DONORIDX": 0,
            },
            # One hydrophobic contact on ligand carbon 7
            {
                "FRAME": 3,
                "INTERACTION": "hydrophobic",
                "Prot_partner": "14LEUA",
                "LIGCARBONIDX": 7,
            },
        ]
    )


def test_schema_ligand_aggregates_by_ligand_atom(sample_interaction_list):
    """schema='ligand' covers all residues together."""
    bm_lig = BindingModeProcesser(
        pdb_md=None,
        ligand="LIG",
        peptide=None,
        special=None,
        ligand_rings=None,
        interaction_list=sample_interaction_list.copy(deep=True),
        threshold=0,  # include all interactions for this test
        total_frames=10,
        schema="ligand",
    )

    # Two hbonds share ACCEPTORIDX=5 => ONE ligand-centric column
    assert "ligand_5_Acceptor_hbond" in bm_lig.unique_data
    hbond_rows = bm_lig.interaction_list[bm_lig.interaction_list["INTERACTION"] == "hbond"]
    assert (hbond_rows["ligand_5_Acceptor_hbond"] == 1).all()

    # Hydrophobic should also be ligand-centric
    assert "ligand_7_hydrophobic" in bm_lig.unique_data
    hydro_rows = bm_lig.interaction_list[bm_lig.interaction_list["INTERACTION"] == "hydrophobic"]
    assert (hydro_rows["ligand_7_hydrophobic"] == 1).all()


def test_schema_residue_keeps_residue_specificity(sample_interaction_list):
    """schema='residue' covers each residues individually."""
    bm_res = BindingModeProcesser(
        pdb_md=None,
        ligand="LIG",
        peptide=None,
        special=None,
        ligand_rings=None,
        interaction_list=sample_interaction_list.copy(deep=True),
        threshold=0,
        total_frames=10,
        schema="residue",
    )

    # Residue schema produces two distinct columns for the two different residues
    assert "12GLYA_5_Acceptor_hbond" in bm_res.unique_data
    assert "13SERA_5_Acceptor_hbond" in bm_res.unique_data

    # Must NOT collapse to ligand_* naming
    assert "ligand_5_Acceptor_hbond" not in bm_res.unique_data


def test_interactions_all_includes_rare_contacts_for_ligand_schema(sample_interaction_list):
    """interactions_all/unique_data_all must include all observed interactions regardless of threshold."""
    bm_lig_high_thresh = BindingModeProcesser(
        pdb_md=None,
        ligand="LIG",
        peptide=None,
        special=None,
        ligand_rings=None,
        interaction_list=sample_interaction_list.copy(deep=True),
        threshold=50,  # 50% of total_frames=10 => min count 5 (none of these pass)
        total_frames=10,
        schema="ligand",
    )

    # High threshold filters these out of the "main" set
    assert "ligand_5_Acceptor_hbond" not in bm_lig_high_thresh.unique_data
    assert "ligand_7_hydrophobic" not in bm_lig_high_thresh.unique_data

    # But "all interactions" must still contain them
    assert "ligand_5_Acceptor_hbond" in bm_lig_high_thresh.unique_data_all
    assert "ligand_7_hydrophobic" in bm_lig_high_thresh.unique_data_all
