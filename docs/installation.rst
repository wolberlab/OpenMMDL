**Installation**
===============

This page details how to install **OpenMMDL** on your local device. 


The first step of the installation of **OpenMMDL** on your device consists in the cloning of **OpenMMDL** to your home directory.
This can be achieved by using the following command lines:

.. code-block:: text

    cd ~
    git clone https://github.com/wolberlab/OpenMMDL


The next step consists of installing all required dependencies in a new separate environment, which will be used to run the **OpenMMDL** package.
**OpenMMDL** is written in Python 3.10 and uses several packages, which can be easily installed on a separate environment using conda (we recommend using miniconda):

.. code-block:: text

    cd OpenMMDL
    conda create -n openmmdl c conda-forge --file requirements.txt
    
After installation, activate the conda environment and install the **OpenMMDL** package:
(Be sure, you are still within the OpenMMDL folder)

.. code-block:: text

    conda activate openmmdl
    pip install .
    
Now with **OpenMMDL** being installed on the local device, everything is set up to work with **OpenMMDL**.
