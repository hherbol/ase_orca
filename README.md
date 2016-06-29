# ASE Orca

This is currently a rudimentary patch to the Atomic Simulation Enviroment (ASE) to allow for the use of Orca simulations.  Simply copy these files into their respective folders and you should be all set to go.  The examples folder shows example code, and the xyz folder has example xyz files needed for these codes.

## Installation

To install everything, first install pip, then the ASE library, then copy this git repo to the ASE python code.

### Installing Pip

To install pip, follow the guidelines [here](https://pip.pypa.io/en/stable/installing/).  Simply put, you'll download "get-pip.py" and then run it using python.

### Installing ASE

After installing pip, you'll be able to install ASE easily:

    pip install --upgrade --user ase

In linux, you'll find the python files in a folder such as:

    ~/.local/lib/python2.7/site-packages/ase

### Patching ASE for Orca

    cd ~/.local/lib/python2.7/site-packages/ase/
    git init
    git remote add origin git@github.com:hherbol/ase_orca.git
    git fetch
    git checkout -t origin/master

