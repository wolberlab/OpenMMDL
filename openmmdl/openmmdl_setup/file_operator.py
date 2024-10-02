import tempfile
import shutil
from werkzeug.utils import secure_filename
from typing import Dict, List, Tuple, Any
import tempfile
import shutil
from werkzeug.utils import secure_filename
from werkzeug.datastructures import FileStorage

from typing import Dict, List, Tuple
import os


class LigandExtractor:
    @staticmethod
    def extract_ligand_name(lig_file_name: str) -> str:
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
        if not isinstance(lig_file_name, str):
            raise TypeError("lig_file_name must be a string")

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
    def save_uploaded_files(
        uploadedFiles: Dict[str, List[Tuple[Any, str]]], request: Any
    ) -> None:
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
        if not isinstance(uploadedFiles, dict):
            raise TypeError("uploadedFiles must be a dictionary")

        if not hasattr(request, "files"):
            raise TypeError("request object must have a 'files' attribute")

        uploadedFiles.clear()
        for key in request.files:
            filelist: List[Tuple[Any, str]] = []
            for file in request.files.getlist(key):
                if not isinstance(file, FileStorage):
                    raise TypeError("file must be a FileStorage instance")

                temp = tempfile.TemporaryFile()
                shutil.copyfileobj(file, temp)
                filelist.append((temp, secure_filename(file.filename)))
            uploadedFiles[key] = filelist
