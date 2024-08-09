

class LigandExtractor:
    @staticmethod
    def extract_ligand_name(lig_file_name):
        """
        Extract the ligand name based on the file type.

        Params:
        -------
        lig_file_name: str
            The name of the ligand file (including the extension).

        Returns:
        --------
        lig_name: str
            The extracted ligand name.
        """
        if lig_file_name.endswith(".sdf"):
            lig_name = "UNL"
        elif lig_file_name.endswith(".pdb"):
            lig_name = lig_file_name[:-4]
        else:
            raise ValueError("Unsupported file format. Only .sdf and .pdb are supported.")

        return lig_name
