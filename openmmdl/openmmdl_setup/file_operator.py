import tempfile
import shutil
from werkzeug.utils import secure_filename


class LigandExtractor:
    @staticmethod
    def extract_ligand_name(lig_file_name):
        """
        Extracts the ligand name from the file name based on its extension.

        This method determines the ligand name by stripping the file extension from the provided ligand
        file name. It handles files with `.pdb` extensions by removing the extension, while it assigns a
        default name `"UNL"` to files with `.sdf` extensions.

        Params:
        -------
        lig_file_name: str
            The name of the ligand file, including its extension (e.g., `ligand.pdb` or `ligand.sdf`).

        Returns:
        --------
        lig_name: str
            The extracted ligand name. For `.pdb` files, this is the file name without the `.pdb` extension.
            For `.sdf` files, it returns a default name `"UNL"`.
        """
        if lig_file_name.endswith(".sdf"):
            lig_name = "UNL"
        elif lig_file_name.endswith(".pdb"):
            lig_name = lig_file_name[:-4]
        else:
            raise ValueError(
                "Unsupported file format. Only .sdf and .pdb are supported."
            )

        return lig_name


class FileUploader:
    """
    Handles uploading and temporary storage of files.
    """

    @staticmethod
    def save_uploaded_files(uploadedFiles, request):
        """
        Saves files from the request into temporary storage and updates the given dictionary.

        Parameters
        ----------
        uploadedFiles : dict
            Dictionary to store uploaded files. Keys are field names; values are lists of (temp file, filename) tuples.

        request : object
            The request object containing uploaded files.

        Returns
        -------
        None
        """
        uploadedFiles.clear()
        for key in request.files:
            filelist = []
            for file in request.files.getlist(key):
                temp = tempfile.TemporaryFile()
                shutil.copyfileobj(file, temp)
                filelist.append((temp, secure_filename(file.filename)))
            uploadedFiles[key] = filelist
