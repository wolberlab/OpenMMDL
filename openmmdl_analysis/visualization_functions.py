import json
import re
import mdtraj as md
import MDAnalysis as mda
import pickle
import nglview as nv

from openmmdl_analysis_function import waterids_barcode_generator


def interacting_water_ids(df_all, waterbridge_interactions):
    interacting_waters = []
    for waterbridge_interaction in waterbridge_interactions:
        waterid_barcode = waterids_barcode_generator(df_all, waterbridge_interaction)
        for waterid in waterid_barcode:
            if waterid != 0:
                interacting_waters.append(waterid)
    return list(set(interacting_waters))


def save_interacting_waters_trajectory(pdb_file_path, dcd_file_path, interacting_waters, outputpath='./'):
    u = mda.Universe(pdb_file_path, dcd_file_path)
    water_atoms = u.select_atoms("protein or nucleic or resname UNK")

    for water in interacting_waters:
        add_water_atoms = u.select_atoms(f'resname HOH and resid {water}')
        water_atoms = water_atoms + add_water_atoms

    water_atoms.write(f'{outputpath}interacting_waters.pdb')

    with mda.Writer(f'{outputpath}interacting_waters.dcd', water_atoms.n_atoms) as W:
        for ts in u.trajectory:
            W.write(water_atoms)


def cloud_json_generation(df_all):
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
    clouds['hydrophobic'] = {"coordinates": hydrophobe_coords, "color": [1.0, 1.0, 0.0], "radius": 0.05}
    clouds['acceptor'] = {"coordinates": acceptor_ccords, "color": [1.0, 0.0, 0.0], "radius": 0.05}
    clouds['donor'] = {"coordinates": donor_coords, "color": [0.0, 1.0, 0.0], "radius": 0.05}
    clouds['waterbridge'] = {"coordinates": waterbridge_coords, "color": [0.0, 1.0, 0.9], "radius": 0.05}
    clouds['negative_ionizable'] = {"coordinates": negative_ionizable_coords, "color": [0.0, 0.0, 1.0], "radius": 0.05}
    clouds['positive_ionizable'] = {"coordinates": positive_ionizable_coords, "color": [1.0, 0.0, 0.0], "radius": 0.05}
    clouds['pistacking'] = {"coordinates": pistacking_coords, "color": [0.0, 0.0, 1.0], "radius": 0.05}
    clouds['pication'] = {"coordinates": pication_coords, "color": [0.0, 0.0, 1.0], "radius": 0.05}
    clouds['halogen'] = {"coordinates": halogen_coords, "color": [1.0, 0.0, 0.9], "radius": 0.05}
    clouds['metal'] = {"coordinates": metal_coords, "color": [1.0, 0.6, 0.0], "radius": 0.05}

    return clouds


def visualization(json_file_path,pdb_file_path, dcd_file_path, interacting_waters_file_path, receptor_type='protein', height='1200px', width='1200px'):
        
    with open(json_file_path) as f:
        data = json.load(f)

    sphere_buffers = []
    for name, cloud in data.items():
        sphere_buffer = {"position": [], "color": [], "radius": []}
        for point in cloud["coordinates"]:
            sphere_buffer["position"] += point
            sphere_buffer["color"] += cloud["color"]
            sphere_buffer["radius"] += [cloud["radius"]]
        sphere_buffers.append(sphere_buffer)
    
    pdb_structure = md.load(pdb_file_path)
    dcd_trajectory = md.load(dcd_file_path, top=pdb_structure)
    with open(interacting_waters_file_path, 'rb') as f:
        interacting_watersids = pickle.load(f)

    view = nv.show_mdtraj(dcd_trajectory)
    view.clear_representations()
    view.add_cartoon(selection=receptor_type)

    for water in interacting_watersids:
        view.add_licorice(selection=f"water and {water}")
    view.add_licorice(selection="UNK")

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
    view.display(gui=True, style="ngl")

        

