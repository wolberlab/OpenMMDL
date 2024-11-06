import json
import re
import MDAnalysis as mda
import pickle
import nglview as nv
import subprocess
import os
import shutil
from typing import List, Dict, Any, Optional, Union
from openmmdl.openmmdl_analysis.barcode_generation import BarcodeGenerator


class TrajectorySaver:
    def __init__(
        self, pdb_md: mda.Universe, ligname: str, special: str, nucleic: bool
    ) -> None:
        """Initializes the TrajectorySaver with an mda.Universe object, ligand name, special residue name and receptor type.

        Args:
            pdb_md (mda.Universe): MDAnalysis Universe object containing the trajectory.
            ligname (str): Name of the ligand in the pdb file.
            special (str): Name of the special residue/ligand in the pdb file (e.g., HEM).
            nucleic (bool): True if the receptor is nucleic, False otherwise.
        """
        self.pdb_md = pdb_md
        self.ligname = ligname
        self.special = special
        self.nucleic = nucleic

    def save_interacting_waters_trajectory(
        self, interacting_waters: List[int], outputpath: str = "./Visualization/"
    ) -> None:
        """Saves .pdb and .dcd files of the trajectory containing ligand, receptor and all interacting waters.

        Args:
            interacting_waters (List[int]): List of all interacting water IDs.
            outputpath (str, optional): Filepath to output new pdb and dcd files. Defaults to './Visualization/'.
        """
        water_atoms = self.pdb_md.select_atoms(
            f"protein or nucleic or resname {self.ligname} or resname {self.special}"
        )

        for water in interacting_waters:
            add_water_atoms = self.pdb_md.select_atoms(f"resname HOH and resid {water}")
            water_atoms = water_atoms + add_water_atoms

        water_atoms.write(f"{outputpath}interacting_waters.pdb")

        with mda.Writer(
            f"{outputpath}interacting_waters.dcd", water_atoms.n_atoms
        ) as W:
            for ts in self.pdb_md.trajectory:
                W.write(water_atoms)

    def save_frame(
        self, frame: int, outpath: str, selection: Optional[str] = None
    ) -> None:
        """Saves a single frame of the trajectory.

        Args:
            frame (int): Number of the frame to save.
            outpath (str): Path to save the frame to.
            selection (Optional[str], optional): MDAnalysis selection string. Defaults to None.
        """
        self.pdb_md.trajectory[frame]
        if selection:
            frame_atomgroup = self.pdb_md.atoms[selection]
        else:
            frame_atomgroup = self.pdb_md.atoms
        frame_atomgroup.write(outpath)


class Visualizer:
    def __init__(
        self,
        md: mda.Universe,
        cloud_path: str,
        ligname: str,
        special: Optional[str] = None,
    ) -> None:
        self.md = md
        self.cloud = self.load_cloud(cloud_path)
        self.ligname = ligname
        self.special = special

    def load_cloud(
        self, cloud_path: str
    ) -> Dict[str, Dict[str, Union[List[float], List[int]]]]:
        """Loads interaction cloud data from a JSON file.

        Args:
            cloud_path (str): Path to the cloud data JSON file.

        Returns:
            Dict[str, Dict[str, Union[List[float], List[int]]]]: The loaded cloud data.
        """
        with open(cloud_path, "r") as f:
            data = json.load(f)
        return data

    def visualize(
        self,
        receptor_type: str = "protein or nucleic",
        height: str = "1000px",
        width: str = "1000px",
    ):
        """Generates visualization of the trajectory with the interacting waters and interaction clouds.

        Args:
            receptor_type (str, optional): Type of receptor. Defaults to 'protein or nucleic'.
            height (str, optional): Height of the visualization. Defaults to '1000px'.
            width (str, optional): Width of the visualization. Defaults to '1000px'.

        Returns:
            nglview widget: Returns an nglview.widget object containing the visualization.
        """

        sphere_buffers = []
        for name, cloud in self.cloud.items():
            sphere_buffer = {"position": [], "color": [], "radius": []}
            for point in cloud["coordinates"]:
                sphere_buffer["position"] += point
                sphere_buffer["color"] += cloud["color"]
                sphere_buffer["radius"] += [cloud["radius"]]
            sphere_buffers.append(sphere_buffer)

        with open(f"interacting_waters.pkl", "rb") as f:
            interacting_watersids = pickle.load(f)

        view = nv.show_mdanalysis(self.md)
        view.clear_representations()
        view.add_cartoon(selection=receptor_type)

        for water in interacting_watersids:
            view.add_licorice(selection=f"water and {water}")
        view.add_licorice(selection=self.ligname)
        if self.special:
            view.add_licorice(selection=self.special)

        for sphere_buffer, name in zip(
            sphere_buffers,
            [
                "hydrophobic",
                "acceptor",
                "donor",
                "waterbridge",
                "negative_ionizable",
                "positive_ionizable",
                "pistacking",
                "pication",
                "halogen",
                "metal",
            ],
        ):
            js = f"""
            var params = {sphere_buffer};
            var shape = new NGL.Shape('{name}');
            var buffer = new NGL.SphereBuffer(params);
            shape.addBuffer(buffer);
            var shapeComp = this.stage.addComponentFromObject(shape);
            shapeComp.addRepresentation("buffer");
            """
            view._js(js)
        view.layout.width = width
        view.layout.height = height
        return view


def run_visualization() -> None:
    """Runs the visualization notebook in the current directory. The visualization notebook is copied from the package directory to the current directory and automatically started."""
    package_dir = os.path.dirname(__file__)
    notebook_path = os.path.join(package_dir, "visualization.ipynb")
    current_dir = os.getcwd()
    shutil.copyfile(notebook_path, f"{current_dir}/visualization.ipynb")
    subprocess.run(["jupyter", "notebook", "visualization.ipynb"])
