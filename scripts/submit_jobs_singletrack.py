#!/usr/bin/env python3

from submit import *

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--series', nargs='+', help='Series to process')
    parser.add_argument('--debug', action='store_true', help='Only show the commands without actually submitting jobs')
    parser.add_argument('--verbose', action='store_true', help='Show the slurm script content')
    parser.add_argument('--hours', type=int, default=6, help='Requested time in hours per job')
    parser.add_argument('--system', type=str, default='slurm', help='batch system, {slurm, lsf}')
    parser.add_argument('--slurm_account', type=str, default='rrg-mdiamond', help='Slurm system account')
    parser.add_argument('--slurm_partition', type=str, default='', help='Slurm system partition')
    
    args = parser.parse_args()
    Series = args.series #["25201106_180710"]#25201106_174906
    DEBUG     = args.debug
    verbose   =args.verbose
    hours = args.hours
    system = args.system
    slurm_account = args.slurm_account
    slurm_partition = args.slurm_partition    
    
    
    
    Log_Dir = '/project/def-mdiamond/tomren/mathusla/data/fit_study_6layer/log/'
    DataDir = '/project/def-mdiamond/tomren/mathusla/data/fit_study_6layer/'

    simulation='/project/def-mdiamond/tomren/mathusla/Mu-Simulation/simulation '
    tracker='/project/def-mdiamond/tomren/mathusla/MATHUSLA-Kalman-Algorithm/tracker/build/tracker '
    tracker_noscatter='/project/def-mdiamond/tomren/mathusla/MATHUSLA-Kalman-Algorithm_noscatter/tracker/build/tracker '

    EnergyList=[[0.1, 0.2, 0.5, 1, 3, 10, 30, 100], [0.1, 0.2, 0.5, 1, 3, 10, 30, 100], [0.1, 0.2, 0.5, 1, 3, 10, 30, 100]]
    EventCount=40000
    TrackerRuns=1
    Scripts= ['muon_gun_tom_range.mac','pion_gun_tom_range.mac','electron_gun_tom_range.mac']
    Names = ["muon", "pion", "electron"]
    CORES = 1
    

    # make directory for log
    os.system(f"mkdir -p {Log_Dir}")

    for i, sim_script in enumerate(Scripts):
        for energy in EnergyList[i]:
            job_script=f"""mkdir -p {DataDir}/{Names[i]}_{energy}_GeV 
# {simulation} -j1 -q  -o {DataDir}/{Names[i]}_{energy}_GeV  -s {sim_script} energy {energy} count {EventCount}  
for f in {DataDir}/{Names[i]}_{energy}_GeV/*/*/run*.root; do 
    {tracker} $f `dirname $f` 
    mv `dirname $f`/stat0.root `dirname $f`/stat_seedmod.root -f 
    
    {tracker_noscatter} $f `dirname $f` 
    mv `dirname $f`/stat0.root `dirname $f`/stat_noscatter.root -f     
done 
"""
#             if i==0:
#                 continue
                
            # Run longer time for pions at 10 GeV and above.
            if (i in [1,2]) and energy >3:
                hours_mod = 4
            else:
                # continue
                hours_mod = hours
                
            # This is the core function to submit job script
            script_prefix=sim_script.split("_")[0]
            submit_script(job_script, f"{script_prefix}_{energy}", blockidx=0, hours=hours_mod, cores=CORES, log_dir=Log_Dir, job_name="reco", system=system, slurm_account=slurm_account, slurm_partition='', debug=DEBUG, verbose=verbose)
            
if __name__ == "__main__":
    main()            
