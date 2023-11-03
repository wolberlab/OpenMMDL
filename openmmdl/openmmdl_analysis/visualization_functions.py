import json
import re
import mdtraj as md
import MDAnalysis as mda
import pickle
import nglview as nv
import subprocess
import os

from openmmdl.openmmdl_analysis.barcode_generation import waterids_barcode_generator


def interacting_water_ids(df_all, waterbridge_interactions):
    """Generates a list of all water ids that form water bridge interactions.

    Args:
        df_all (pandas dataframe): dataframe containing all interactions from plip analysis (typicaly df_all)
        waterbridge_interactions (list): list of strings containing the names of all water bridge interactions

    Returns:
        list: list of all unique water ids that form water bridge interactions
    """
    interacting_waters = []
    for waterbridge_interaction in waterbridge_interactions:
        waterid_barcode = waterids_barcode_generator(df_all, waterbridge_interaction)
        for waterid in waterid_barcode:
            if waterid != 0:
                interacting_waters.append(waterid)
    return list(set(interacting_waters))


def save_interacting_waters_trajectory(pdb_file_path, dcd_file_path, interacting_waters, ligname, outputpath='./'):
    """Saves .pdb and .dcd files of the trajectory containing ligand, receptor and all interacting waters.

    Args:
        pdb_file_path (str): path to original pdb file
        dcd_file_path (str): path to original dcd file
        interacting_waters (list): list of all interacting water ids
        ligname (str): name of the ligand in the pdb file
        outputpath (str, optional): filepath to output new pdb and dcd files. Defaults to './'.
    """
    u = mda.Universe(pdb_file_path, dcd_file_path)
    water_atoms = u.select_atoms(f"protein or nucleic or resname {ligname}")

    for water in interacting_waters:
        add_water_atoms = u.select_atoms(f'resname HOH and resid {water}')
        water_atoms = water_atoms + add_water_atoms

    water_atoms.write(f'{outputpath}interacting_waters.pdb')

    with mda.Writer(f'{outputpath}interacting_waters.dcd', water_atoms.n_atoms) as W:
        for ts in u.trajectory:
            W.write(water_atoms)


def cloud_json_generation(df_all):
    """generates dict for visualization of interaction clouds. Later saved as .json file.

    Args:
        df_all (pandas dataframe): dataframe containing all interactions from plip analysis (typicaly df_all)

    Returns:
        dict: dict containing all interaction clouds
    """     
    coord_pattern = re.compile(r'\(([\d.-]+), ([\d.-]+), ([\d.-]+)\)')
    hydrophobe_coords = []
    acceptor_ccords = []
    donor_coords = []
    waterbridge_coords = []
    negative_ionizable_coords = [] 
    positive_ionizable_coords = []
    pistacking_coords = []
    pication_coords = []
    halogen_coords = []
    metal_coords = []

    for index, row in df_all.iterrows():
        coord_match = coord_pattern.match(row['LIGCOO'])
        if coord_match:
            x, y, z = map(float, coord_match.groups())
            x, y, z = round(x, 3), round(y, 3), round(z, 3)
            interaction = row['INTERACTION']
            if interaction == 'hbond':
                if row['PROTISDON'] == 'True':
                    interaction = 'donor'
                else:
                    interaction = 'acceptor'
            if interaction == 'saltbridge':
                if row['PROTISPOS'] == 'True':
                    interaction = 'negative_ionizable'
                else:
                    interaction = 'positive_ionizable'
            if interaction == 'hydrophobic':
                hydrophobe_coords.append([x, y, z])
            if interaction == 'acceptor':
                acceptor_ccords.append([x, y, z])
            if interaction == 'donor': 
                donor_coords.append([x, y, z])
            if interaction == 'waterbridge':
                waterbridge_coords.append([x, y, z])
            if interaction == 'negative_ionizable':
                negative_ionizable_coords.append([x, y, z])
            if interaction == 'positive_ionizable':
                positive_ionizable_coords.append([x, y, z])
            if interaction == 'pistacking':
                pistacking_coords.append([x, y, z])
            if interaction == 'pication':
                pication_coords.append([x, y, z])
            if interaction == 'halogen':
                halogen_coords.append([x, y, z])
            if interaction == 'metal':
                metal_coords.append([x, y, z])


    clouds = {}
    clouds['hydrophobic'] = {"coordinates": hydrophobe_coords, "color": [1.0, 1.0, 0.0], "radius": 0.1}
    clouds['acceptor'] = {"coordinates": acceptor_ccords, "color": [1.0, 0.0, 0.0], "radius": 0.1}
    clouds['donor'] = {"coordinates": donor_coords, "color": [0.0, 1.0, 0.0], "radius": 0.1}
    clouds['waterbridge'] = {"coordinates": waterbridge_coords, "color": [0.0, 1.0, 0.9], "radius": 0.1}
    clouds['negative_ionizable'] = {"coordinates": negative_ionizable_coords, "color": [0.0, 0.0, 1.0], "radius": 0.1}
    clouds['positive_ionizable'] = {"coordinates": positive_ionizable_coords, "color": [1.0, 0.0, 0.0], "radius": 0.1}
    clouds['pistacking'] = {"coordinates": pistacking_coords, "color": [0.0, 0.0, 1.0], "radius": 0.1}
    clouds['pication'] = {"coordinates": pication_coords, "color": [0.0, 0.0, 1.0], "radius": 0.1}
    clouds['halogen'] = {"coordinates": halogen_coords, "color": [1.0, 0.0, 0.9], "radius": 0.1}
    clouds['metal'] = {"coordinates": metal_coords, "color": [1.0, 0.6, 0.0], "radius": 0.1}

    return clouds


def visualization(filepath, ligname, receptor_type='protein', height='1000px', width='1000px'):
    """Generates visualization of the trajectory with the interacting waters and interaction clouds.
    
    Args:
        json_file_path (str): path to .json file containing the interaction clouds
        pdb_file_path (str): path to pdb file (use interacting_waters.pdb for better visualization)
        dcd_file_path (str): path to dcd file (use interacting_waters.dcd for better visualization)
        interacting_waters_file_path (str): path to .pkl file containing the interacting water ids
        receptor_type (str, optional): type of receptor. Defaults to 'protein'.
        height (str, optional): height of the visualization. Defaults to '1200px'.
        width (str, optional): width of the visualization. Defaults to '1200px'.
    
    Returns:
        nglview widget: returns the nglview widget containing the visualization
    """
    with open(f'{filepath}clouds.json') as f:
        data = json.load(f)

    sphere_buffers = []
    for name, cloud in data.items():
        sphere_buffer = {"position": [], "color": [], "radius": []}
        for point in cloud["coordinates"]:
            sphere_buffer["position"] += point
            sphere_buffer["color"] += cloud["color"]
            sphere_buffer["radius"] += [cloud["radius"]]
        sphere_buffers.append(sphere_buffer)
    
    pdb_structure = md.load(f'{filepath}interacting_waters.pdb')
    dcd_trajectory = md.load(f'{filepath}interacting_waters.dcd', top=pdb_structure)
    with open(f'{filepath}interacting_waters.pkl', 'rb') as f:
        interacting_watersids = pickle.load(f)

    view = nv.show_mdtraj(dcd_trajectory)
    view.clear_representations()
    view.add_cartoon(selection=receptor_type)

    for water in interacting_watersids:
        view.add_licorice(selection=f"water and {water}")
    view.add_licorice(selection=ligname)

    for sphere_buffer, name in zip(sphere_buffers, ['hydrophobic', 'acceptor', 'donor', 'waterbridge', 'negative_ionizable', 'positive_ionizable', 'pistacking', 'pication', 'halogen', 'metal']):
        js = (
        f"""
        var params = {sphere_buffer};
        var shape = new NGL.Shape('{name}');
        var buffer = new NGL.SphereBuffer(params);
        shape.addBuffer(buffer);
        var shapeComp = this.stage.addComponentFromObject(shape);
        shapeComp.addRepresentation("buffer");
        """
        )
        view._js(js)
    view.layout.width = width
    view.layout.height = height
    return view


def run_visualization():
    package_dir = os.path.dirname(__file__)
    notebook_path = os.path.join(package_dir, 'visualization.ipynb')
    subprocess.run(['jupyter', 'notebook', notebook_path])
