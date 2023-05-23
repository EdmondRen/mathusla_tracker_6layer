#!/usr/bin/env python3

from submit import *

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--series', nargs='+', help='Series to process')
    parser.add_argument('--debug', action='store_true', help='Only show the commands without actually submitting jobs')
    parser.add_argument('--verbose', action='store_true', help='Show the slurm script content')
    parser.add_argument('--hours', type=int, default=6, help='Requested time in hours per job')
    parser.add_argument('--system', type=str, default='slurm', help='batch system, {slurm, lsf}')
    parser.add_argument('--slurm_account', type=str, default='def-mdiamond', help='Slurm system account')
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

    EnergyList=[[100]]
    EventCount=160000
    TrackerRuns=1
    Scripts= ['filereader_muon_M4000_P40000_N10000.mac','filereader_mumu_range_CenterModule.mac','filereader_mumu_range_CenterModule_1m.mac']
    Names = ["muon"]
    CORES = 1
    

    # make directory for log
    os.system(f"mkdir -p {Log_Dir}")
    energy = EnergyList[0][0]
    hours_mod = hours
    
    for ijob in range(len(Scripts)):
        sim_script = Scripts[ijob]
        # if ijob==0:
        #     data_dir_sub = f"{DataDir}/XtoMuMu_M4GeV_P40GeV/"
        # elif ijob==1:
        #     data_dir_sub = f"{DataDir}/XtoMuMu_P10GeV_manual/"
        # elif ijob==2:
        if ijob==2:
            data_dir_sub = f"{DataDir}/XtoMuMu_P10GeV_manual_1m/"            
        else:
            continue



        job_script=f"""mkdir -p {data_dir_sub} 
{simulation} -j1 -q  -o {data_dir_sub}  -s {sim_script}  
for f in {data_dir_sub}/*/*/run*.root; do 
    {tracker} $f `dirname $f` 
    mv `dirname $f`/stat0.root `dirname $f`/stat_seedmod.root -f 
done 
"""
        # This is the core function to submit job script
        script_prefix=sim_script.split("_")[0]
        submit_script(job_script, f"{script_prefix}_{energy}_2", blockidx=0, hours=hours_mod, cores=CORES, log_dir=Log_Dir, job_name="reco", system=system, slurm_account=slurm_account, slurm_partition='', debug=DEBUG, verbose=verbose)
    
    
            
if __name__ == "__main__":
    main()            
