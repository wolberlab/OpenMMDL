import json
import re
import pickle
import nglview as nv
import subprocess
import os
import shutil


class Visualizer:
    """
    Visualizes the OpenMMDL Analysis output of the trajectories with annotated interaction clouds and water molecule interactions.

    The class provides tools to render 3D molecular scenes using NGLview, including trajectory
    visualization, ligand highlighting, interacting water molecules, and various interaction cloud types
    (e.g., hydrogen bonds, hydrophobic regions, etc.).

    Attributes
    ----------
    md : MDAnalysis.Universe
        The loaded trajectory and topology data used for visualization.
    cloud : dict
        Dictionary parsed from the JSON file representing interaction clouds (e.g., hydrophobic, donor, etc.).
    ligname : str
        Ligand selection string used in NGLview for highlighting.
    special : str
        Additional selection string for visual emphasis (can be None).
    """
    def __init__(self, md, cloud_path, ligname, special):
        self.md = md
        self.cloud = self._load_cloud(cloud_path)
        self.ligname = ligname
        self.special = special

    def _load_cloud(self, cloud_path):
        """
        Loads interaction clouds from a JSON file.
    
        Parameters
        ----------
        cloud_path : str
            Path to the JSON file containing the interaction clouds.
    
        Returns
        -------
        dict
            Dictionary representing different types of interaction clouds and their properties.
        """
        with open(cloud_path, "r") as f:
            data = json.load(f)

        return data

    def visualize(
        self, receptor_type="protein or nucleic", height="1000px", width="1000px"
    ):
        """
        Generates visualization of the trajectory with the interacting waters and interaction clouds.

        Parameters
        ----------
        receptor_type : str, optional
            Type of receptor. Defaults to 'protein or nucleic'.
        height : str, optional
            Height of the visualization. Defaults to '1000px'.
        width : str, optional
            Width of the visualization. Defaults to '1000px'.

        Returns
        -------
        nglview.NGLWidget
            Returns an nglview.widget object containing the visualization
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


def run_visualization():
    """
    Runs the visualization notebook in the current directory. 
    The visualization notebook is copied from the package directory to the current directory and automaticaly started.
    """
    package_dir = os.path.dirname(__file__)
    notebook_path = os.path.join(package_dir, "visualization.ipynb")
    current_dir = os.getcwd()
    shutil.copyfile(notebook_path, f"{current_dir}/visualization.ipynb")
    subprocess.run(["jupyter", "notebook", "visualization.ipynb"])
