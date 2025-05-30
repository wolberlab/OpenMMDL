**Installation**
===============

This page details how to install **OpenMMDL** on your local device. 

Installation via conda-forge
------------------------------

**OpenMMDL** is implemented in conda-forge and can be installed for linux-based system with the following command line:

.. code-block:: text

    conda install -c conda-forge openmmdl

Installation via repository
------------------------------
**OpenMMDL** can be installed through the cloning of the repository and the installation of the required packages listed in *environment.yml* or *requirements.txt*

The first step of the installation of **OpenMMDL** on your device consists in the cloning of **OpenMMDL** to your home directory.
This can be achieved by using the following command lines:

.. code-block:: text

    cd ~
    git clone https://github.com/wolberlab/OpenMMDL

**OpenMMDL** is written in Python 3.10 and uses several packages, which can be easily installed on a separate environment using conda (we recommend using miniconda):

.. code-block:: text

    cd OpenMMDL
    conda create -n openmmdl -c conda-forge --file requirements.txt

After installation, activate the conda environment:

.. code-block:: text

    conda activate openmmdl

Install the **OpenMMDL** package with the required entry points via pip:

.. code-block:: text

    pip install .

This should lead to the environment having now all the required packages and entry points for **OpenMMDL**.

Installation via Docker
------------------------------

**OpenMMDL** is mainly supported for Linux distribution systems, thus for Windows and MacOS the installation with docker may be preferred, due to docker creating an image with *ubuntu:22.04*.

For this first clone this repository:

.. code-block:: text

    cd ~
    git clone https://github.com/wolberlab/OpenMMDL

change the path to the **OpenMMDL** folder and build the docker image from there:

.. code-block:: text

    cd OpenMMDL
    docker build -t openmmdl_env .

This will build the **OpenMMDL** image with docker. Now that it is build you can access it through an interactive terminal:

.. code-block:: text

    docker run -it --name openmmdl_test openmmdl_env
    
From there you can access all the **OpenMMDL** entry points. Currently due to **OpenMMDL Setup** using flask it can be difficult to access it through the docker image.
