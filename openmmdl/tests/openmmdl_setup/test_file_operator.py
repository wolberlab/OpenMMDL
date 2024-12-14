import pytest
from werkzeug.datastructures import FileStorage
from io import BytesIO
from openmmdl.openmmdl_setup.file_operator import LigandExtractor, FileUploader  # Replace with the actual module name
import tempfile

# Test cases for LigandExtractor
class TestLigandExtractor:
    def test_extract_ligand_name_pdb(self):
        lig_file_name = "ligand.pdb"
        expected_output = "ligand"
        assert LigandExtractor.extract_ligand_name(lig_file_name) == expected_output

    def test_extract_ligand_name_sdf(self):
        lig_file_name = "ligand.sdf"
        expected_output = "UNL"
        assert LigandExtractor.extract_ligand_name(lig_file_name) == expected_output

    def test_extract_ligand_name_invalid_extension(self):
        lig_file_name = "ligand.txt"
        with pytest.raises(ValueError, match="Unsupported file format. Only .sdf and .pdb are supported."):
            LigandExtractor.extract_ligand_name(lig_file_name)

    def test_extract_ligand_name_non_string(self):
        lig_file_name = 12345
        with pytest.raises(TypeError, match="lig_file_name must be a string"):
            LigandExtractor.extract_ligand_name(lig_file_name)



# Custom class to mock the behavior of request.files in Flask
class MockFileMultiDict:
    """
    This class mimics the behavior of Flask's request.files MultiDict for testing purposes.
    """
    def __init__(self, file_dict):
        self.file_dict = file_dict

    def getlist(self, key):
        return self.file_dict.get(key, [])

    def __iter__(self):
        return iter(self.file_dict)

    def items(self):
        return self.file_dict.items()

# Test cases for FileUploader
class TestFileUploader:
    @pytest.fixture
    def fake_request(self):
        """
        Simulates a fake request with files using FileStorage and provides the getlist() method.
        """
        file1 = FileStorage(
            stream=BytesIO(b"dummy content 1"), filename="file1.txt", content_type="text/plain"
        )
        file2 = FileStorage(
            stream=BytesIO(b"dummy content 2"), filename="file2.txt", content_type="text/plain"
        )

        class FakeRequest:
            def __init__(self):
                # Mimic request.files using the MockFileMultiDict
                self.files = MockFileMultiDict({
                    "file_field_1": [file1],
                    "file_field_2": [file2],
                })

        return FakeRequest()

    def test_save_uploaded_files_success(self, fake_request):
        uploadedFiles = {}
        FileUploader.save_uploaded_files(uploadedFiles, fake_request)

        # Verify that files were saved correctly
        assert "file_field_1" in uploadedFiles
        assert "file_field_2" in uploadedFiles
        assert len(uploadedFiles["file_field_1"]) == 1
        assert uploadedFiles["file_field_1"][0][1] == "file1.txt"
        assert len(uploadedFiles["file_field_2"]) == 1
        assert uploadedFiles["file_field_2"][0][1] == "file2.txt"

    def test_save_uploaded_files_invalid_dict(self, fake_request):
        uploadedFiles = []
        with pytest.raises(TypeError, match="uploadedFiles must be a dictionary"):
            FileUploader.save_uploaded_files(uploadedFiles, fake_request)

    def test_save_uploaded_files_invalid_request(self):
        uploadedFiles = {}
        invalid_request = object()  # Request without 'files' attribute
        with pytest.raises(TypeError, match="request object must have a 'files' attribute"):
            FileUploader.save_uploaded_files(uploadedFiles, invalid_request)

    def test_save_uploaded_files_non_filestorage(self, fake_request):
        # Modify fake_request to include an invalid file type (non-FileStorage instance)
        fake_request.files.file_dict["file_field_1"] = [object()]  # Invalid file type
        uploadedFiles = {}
        with pytest.raises(TypeError, match="file must be a FileStorage instance"):
            FileUploader.save_uploaded_files(uploadedFiles, fake_request)
