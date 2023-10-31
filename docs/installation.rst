**Installation**
===============

This page details how to install openmmdl on your local device. 


The first step of the installation of OpenMMDL on your device, consists in the cloning of OpenMMDL to your home directory.
This can be accoplished by using the following command lines:

.. code-block:: text

    cd ~
    git clone https://github.com/wolberlab/OpenMMDL


The next step consists in installing all required dependencies in a new separate environment, which will be used to run the OpenMMDL package.
OpenMMDL is written in python 3.10 and uses several packages, which can be easily installed on a separate environment using conda (we recommend using miniconda):

.. code-block:: text

    cd OpenMMDL
    conda create -n openmmdl --file requirements.txt
    
After installation, activate the conda environment and install the OpenMMDL package:
(Be sure, you are still within the OpenMMDL folder)

.. code-block:: text

    conda activate openmmdl
    pip install .
    
Now with OpenMMDL being installed on the local device everything is setup to work with OpenMMDL.
