FROM ubuntu:22.04

# Install Miniconda
RUN apt-get update && apt-get install -y wget && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Set environment variables
ENV PATH="/opt/conda/bin:$PATH"

# Create and activate a conda environment
RUN conda create --name openmmdl_env python=3.12 -y && \
    echo "conda activate openmmdl_env" >> ~/.bashrc

# Install OpenMM and other dependencies
RUN conda install -n openmmdl_env -c conda-forge openmmdl -y

# Set environment variables for OpenMM
ENV PATH /opt/conda/envs/openmmdl_env/bin:$PATH
ENV CONDA_DEFAULT_ENV openmmdl_env

# Set working directory
WORKDIR /app

# Copy the OpenMMDL repository that you already cloned locally
COPY . /app

# Install OpenMMDL
RUN pip install .

# Optional: install Jupyter if you want to run notebooks
# RUN pip install jupyterlab

# Default command to run bash
CMD ["/bin/bash"]