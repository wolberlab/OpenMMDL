"""Tests for the guided tutorial inside OpenMMDL Setup."""

import io
import re
import zipfile
from html.parser import HTMLParser

import pytest
from flask import render_template, session

import openmmdl.openmmdl_setup.openmmdlsetup as M
from openmmdl.openmmdl_setup.setup_tutorials import (
    PDB_TUTORIAL_FILES,
    PDB_TUTORIAL_ID,
    PDB_TUTORIAL_PAGES,
)

# Template and dummy data for every page that carries tutorial steps.
PAGE_TEMPLATES = {
    "select_file_type": ("selectFileType.html", {}),
    "configure_pdb_file": ("configurePdbFile.html", {}),
    "select_chains": (
        "selectChains.html",
        {"chains": [("A", 756, "Protein"), ("B", 753, "Protein"), ("C", 4, "BMA, NAG, MAN")]},
    ),
    "add_residues": ("addResidues.html", {"spans": [("A", 101, 111, "SER, GLY, ASN")]}),
    "convert_residues": (
        "convertResidues.html",
        {"residues": [("A", "MSE", "12", ["MET", "SEC"], "MET")]},
    ),
    "add_heavy_atoms": ("addHeavyAtoms.html", {"residues": [("A", "LYS", "33", "CD, CE, NZ")]}),
    "add_hydrogens": ("addHydrogens.html", {"unitCell": None, "boundingBox": (8.0, 7.5, 9.0)}),
    "simulation_options": (
        "simulationOptions.html",
        {
            "display_save_script": True,
            "display_processed_pdb": True,
            "display_save_all_files": True,
        },
    ),
}


class _Elements(HTMLParser):
    """Collects (tag, attrs) for every start tag of a document."""

    def __init__(self):
        super().__init__()
        self.elements = []

    def handle_starttag(self, tag, attrs):
        self.elements.append((tag, dict(attrs)))

    handle_startendtag = handle_starttag


_SELECTOR = re.compile(
    r"(#[\w-]+)|(\.[\w-]+)|\[([\w-]+)(\*?=)?['\"]?([^'\"\]]*)['\"]?\]|([a-z][\w-]*)"
)


def _matches(selector, elements):
    """Minimal CSS matcher for the selector forms used by the tutorial config:
    tag, #id, .class and [attr=value] / [attr*=value] filters, in any combination."""
    conditions = []
    for id_, cls, attr, op, value, tag in _SELECTOR.findall(selector):
        if id_:
            conditions.append(lambda a, v=id_[1:]: a.get("id") == v)
        elif cls:
            conditions.append(lambda a, v=cls[1:]: v in (a.get("class") or "").split())
        elif attr:
            if op == "*=":
                conditions.append(lambda a, k=attr, v=value: v in (a.get(k) or ""))
            elif op == "=":
                conditions.append(lambda a, k=attr, v=value: a.get(k) == v)
            else:
                conditions.append(lambda a, k=attr: k in a)
        elif tag:
            conditions.append(lambda a, t=tag, _tag=None: a.get("__tag__") == t)
    for tag, attrs in elements:
        attrs = dict(attrs, __tag__=tag)
        if all(c(attrs) for c in conditions):
            return True
    return False


def _render(page_key):
    template, kwargs = PAGE_TEMPLATES[page_key]
    with M.app.test_request_context():
        session["setupTutorial"] = PDB_TUTORIAL_ID
        session["fileType"] = "pdb"
        session["dataFields"] = []
        return render_template(template, tutorial_page=page_key, **kwargs)


def test_every_tutorial_page_has_a_template():
    assert set(PDB_TUTORIAL_PAGES) == set(PAGE_TEMPLATES)


@pytest.mark.parametrize("page_key", sorted(PDB_TUTORIAL_PAGES))
def test_tutorial_step_targets_exist_on_page(page_key):
    """Every step target must match an element of the rendered page, so that template
    edits cannot silently detach a tutorial step from its control."""
    html = _render(page_key)
    parser = _Elements()
    parser.feed(html)
    for step in PDB_TUTORIAL_PAGES[page_key]["steps"]:
        for selector in filter(None, [step.get("target"), step.get("done", {}).get("selector")]):
            if _matches(selector, parser.elements):
                continue
            # Rows added by page JavaScript (the additional-molecule row) exist only as a
            # template string in the page source; accept the selector's tokens there.
            tokens = [t for t in re.findall(r"[\w-]{3,}", selector) if t not in ("input", "name")]
            assert all(t in html for t in tokens), f"{page_key}: no element for {selector!r}"


@pytest.mark.parametrize("page_key", sorted(PDB_TUTORIAL_PAGES))
def test_tutorial_page_injects_panel_config(page_key):
    html = _render(page_key)
    assert "window.OpenMMDLTutorial" in html
    assert "openmmdl_tutorial.js" in html
    assert 'class="om-tutorial-active"' in html


def test_context_provides_absolute_docs_and_stop_urls():
    with M.app.test_request_context():
        session["setupTutorial"] = PDB_TUTORIAL_ID
        page = M._setup_tutorial_context("configure_pdb_file")
    assert page["docsUrl"].startswith("https://openmmdl.readthedocs.io/en/")
    assert page["docsUrl"].endswith("pdb_path.html#selecting-input-files")
    assert page["stopUrl"] == "/tutorial/stop"
    assert page["filesUrl"] == "/tutorial/pdb-small-molecule/files"


def test_start_and_stop_toggle_the_session_flag():
    client = M.app.test_client()
    assert client.get("/").status_code == 200
    with client.session_transaction() as sess:
        assert "setupTutorial" not in sess

    resp = client.get("/tutorial/pdb-small-molecule")
    assert resp.status_code == 302 and resp.headers["Location"].endswith("/")
    with client.session_transaction() as sess:
        assert sess["setupTutorial"] == PDB_TUTORIAL_ID
    # The flag survives plain page loads (Back / Start Over targets).
    assert b"openmmdl_tutorial.js" in client.get("/").data
    assert b"Exit Tutorial" in client.get("/").data

    assert client.get("/tutorial/stop").status_code == 302
    with client.session_transaction() as sess:
        assert "setupTutorial" not in sess
    assert b"openmmdl_tutorial.js" not in client.get("/").data


def test_tutorial_files_zip_contains_bundle_and_readme():
    client = M.app.test_client()
    resp = client.get("/tutorial/pdb-small-molecule/files")
    assert resp.status_code == 200
    archive = zipfile.ZipFile(io.BytesIO(resp.data))
    assert archive.namelist() == ["README.txt", *PDB_TUTORIAL_FILES]
    readme = archive.read("README.txt").decode()
    for name in PDB_TUTORIAL_FILES:
        assert name in readme
        assert archive.getinfo(name).file_size > 0
