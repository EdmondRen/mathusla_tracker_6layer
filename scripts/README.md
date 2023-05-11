This directory contains helper functions and scripts to run simulation/tracker.

# Setup the environment:

With `initsim` alias defined in .bashrc, you can do

\# initsim

function initsim() {
        export PYTHIA8=/home/tomren/home/mathusla/pythia8308
        export PYTHIA8DATA=$PYTHIA8/share/Pythia8/xmldoc
        module load qt/5.15.2 gcc/9.3.0 StdEnv/2020   root/6.26.06  eigen/3.3.7
        source ~/GEANT4/install/bin/geant4.sh
}


# Submit jobs

submit.py contains a piece of code to submit a given script to slurm batch system

You can write another script to pass multiple jobs to slurm, for example `submit_jobs_trackerdebug.py`, then run `python submit_jobs_trackerdebug.py`

Existing scripts:

* submit_jobs_singletrack.py: Run simulation
   * Type: muon/pion/electron gun
   * Events: 40000
   * Momentum: 0.1, 0.2, 0.5, 1, 3, 10, 50, 100 GeV/c
   
# Run reconstruction

Use the following command to run tracker
    
    tracker SIMULATION.root OUTPUT_DIR
    
There are shell scripts/ python scripts that can run multiple tracking jobs. For example, `run_recon_singletrack.sh`

# Simulation with ReadFile generator

The `read_file` generator takes pre-defined particles as the input of simulation, effectively acting as a "vertex gun". Helper functions to use the `read_file` generator are in `sim_filereader_helper.py`, and the usage examples can be found in `vertex_gun.ipynb`. The functions are described as following:

### frame_transform()