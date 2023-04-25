# Setup the environment:

With `initsim` alias defined in .bashrc, you can do

\# initsim

function initsim() {
        export PYTHIA8=/home/tomren/home/mathusla/pythia8308
        export PYTHIA8DATA=$PYTHIA8/share/Pythia8/xmldoc
        module load qt/5.15.2 gcc/9.3.0 StdEnv/2020   root/6.26.06  eigen/3.3.7
        source ~/GEANT4/install/bin/geant4.sh
}

# Run simulation

# Submit jobs

submit.py contains a piece of code to submit a given script to slurm batch system

You can write another script to pass multiple jobs to slurm, for example `submit_jobs_trackerdebug.py`, then run `python submit_jobs_trackerdebug.py`

Existing scripts:

* submit_jobs_singletrack.py: Run simulation
   * Type: muon/pion gun
   * Events: 40000
   * Momentum: 0.5, 1, 3, 10, 50, 1000 GeV/c
   
# Run reconstruction

Use the following command to run tracker
    
    tracker SIMULATION.root OUTPUT_DIR
    
There are shell scripts/ python scripts that can run multiple tracking jobs.
