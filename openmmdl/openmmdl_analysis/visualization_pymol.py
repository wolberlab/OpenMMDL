import json
import time
from pymol import cmd, cgo


def visualize_pdb(pdb_path):
    # Load the PDB file
    cmd.load(pdb_path, "molecule")

    # Select and show only proteins and ligands
    cmd.show("cartoon", "polymer")
    cmd.show("sticks", "organic")
    cmd.show("spheres", "solvent")
    cmd.hide("sticks", "hydrogen")


def visualize_pointcloud(json_path, buffer_size=1000):
    with open(json_path, "r") as json_file:
        clusters = json.load(json_file)

    for type, prop in clusters.items():
        r = prop["radius"]
        color = [prop["color"][0], prop["color"][1], prop["color"][2]]
        points_buffer = []

        for i, point in enumerate(prop["coordinates"]):
            points_buffer.extend([cgo.COLOR, *color, cgo.SPHERE, *point, r])

            # Load in chunks of buffer_size
            if (i + 1) % buffer_size == 0:
                cmd.load_cgo(points_buffer, type)
                points_buffer = []  # Reset the buffer
                cmd.refresh()
                time.sleep(0.1)  # Small delay to ensure proper loading

        # Load any remaining points in the buffer
        if points_buffer:
            cmd.load_cgo(points_buffer, type)
            cmd.refresh()


# The function to be called by PyMOL with provided arguments
def main(pdb_path, json_path):
    visualize_pdb(pdb_path)
    visualize_pointcloud(json_path)


cmd.extend("openmmdl_visualization", main)
