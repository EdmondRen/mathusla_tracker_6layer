EnergyList=(1000)
# EnergyList=(1 2 5 10 50 100 1000)

EventCount=200
TrackerRuns=200

DataDir='/project/def-mdiamond/tomren/mathusla/data/fit_study'
Script='muon_gun_tom_range.mac'

simulation='/project/def-mdiamond/tomren/mathusla/Mu-Simulation/simulation '
tracker='/project/def-mdiamond/tomren/mathusla/MATHUSLA-Kalman-Algorithm/tracker/build/tracker'

for energy in ${EnergyList[@]}
do
#   echo "mkdir -p ${DataDir}/muon_${energy}_GeV"
  echo "simulation -j2 -o ${DataDir}/muon_${energy}_GeV -q -s $Script energy $energy count $EventCount"
#   for f in ${DataDir}/muon_${energy}_GeV/*/*/run*.root; do
#     # Run tracker for $TrackerRuns times
#     for ((irun=1; irun<=TrackerRuns;irun++)); do
#         echo "tracker $f `dirname ${f}`" 
#         #Rename the output for a unique index
#         echo "mv `dirname ${f}`/stat0.root `dirname ${f}`/stat${irun}.root"
#     done
#   done
#   mkdir -p ${DataDir}/muon_${energy}_GeV
#   ${simulation} -j2 -o ${DataDir}/muon_${energy}_GeV -q -s $Script energy $energy count $EventCount
  for f in ${DataDir}/muon_${energy}_GeV_manual/*/*/run*.root; do
    # Run tracker for $TrackerRuns times
    for ((irun=1; irun<=TrackerRuns;irun++)); do
        ${tracker} $f `dirname ${f}`
        #Rename the output for a unique index
        mv `dirname ${f}`/stat0.root `dirname ${f}`/stat${irun}.root
    done
  done
done
