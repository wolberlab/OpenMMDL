Installation
===============

This page details how to install openmmdl on your local device. 

**Clone the OpenMMDL Github repository**
The first step of the installation of OpenMMDL on your device, consists in the cloning of OpenMMDL to your home directory.
This can be accoplished by using the following command lines:


.. prompt:: bash $

    cd ~
    git clone https://github.com/wolberlab/OpenMMDL

**Install all required dependencies in a separate environment**
The next step consists in installing all required dependencies in a new separate environment, which will be used to run the openmmdl-setup and openmmdl-simulation.
OpenMMDL is written in python 3.10 and uses several packages, which can be easily installed on a separate environment using conda (we recommend using miniconda):



.. prompt:: bash $

    cd OpenMMDL
    conda create -n openmmdl --file requirements.txt
    
