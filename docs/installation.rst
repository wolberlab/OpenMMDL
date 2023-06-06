**Installation**
===============

This page details how to install openmmdl on your local device. 


The first step of the installation of OpenMMDL on your device, consists in the cloning of OpenMMDL to your home directory.
This can be accoplished by using the following command lines:

.. code-block:: text

    cd ~
    git clone https://github.com/wolberlab/OpenMMDL


The next step consists in installing all required dependencies in a new separate environment, which will be used to run the openmmdl-setup and openmmdl-simulation.
OpenMMDL is written in python 3.10 and uses several packages, which can be easily installed on a separate environment using conda (we recommend using miniconda):

.. code-block:: text

    cd OpenMMDL
    conda create -n openmmdl --file requirements.txt
    
After installation, activate the conda environment:


.. code-block:: text

    conda activate openmmdl
    
The final step of the installation of OpenMMDL consists in the creation of alias for openmmdl-simulation and openmmdl-setup:


.. code-block:: text

    echo 'alias openmmdl-setup="python3 ~/OpenMMDL/openmmdl_setup/openmmdlsetup.py"' >> ~/.bashrc
    echo 'alias openmmdl-simulation="python3 ~/OpenMMDL/openmmdl_simulation/openmmdlsimulation.py"' >> ~/.bashrc
    
    
Now with OpenMMDL being installed on the local device and the alias being assigned everything is setup to work with OpenMMDL.
