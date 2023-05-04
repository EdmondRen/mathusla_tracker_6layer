# EnergyList=(1000)
EnergyList=(0.1 0.2 0.5 1 3 10 30 100)
EnergyList2=(0.5 1 3 10 100)
#EnergyList=(100)

EventCount=200
TrackerRuns=200

DataDir='/project/def-mdiamond/tomren/mathusla/data/fit_study_6layer/'
Script='muon_gun_tom_range.mac'

simulation='/project/def-mdiamond/tomren/mathusla/Mu-Simulation/simulation '
tracker='/project/def-mdiamond/tomren/mathusla/MATHUSLA-Kalman-Algorithm/tracker/build/tracker '
tracker_original='/project/def-mdiamond/tomren/mathusla/MATHUSLA-Kalman-Algorithm_original/tracker/build/tracker '
tracker_noscatter='/project/def-mdiamond/tomren/mathusla/MATHUSLA-Kalman-Algorithm_noscatter/tracker/build/tracker '


for energy in ${EnergyList[@]}
do
#   echo "mkdir -p ${DataDir}/muon_${energy}_GeV"
  # echo "simulation -j2 -o ${DataDir}/muon_${energy}_GeV -q -s $Script energy $energy count $EventCount"
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


#   for f in ${DataDir}/muon_${energy}_GeV/*/*/run*.root; do
#     ${tracker} $f `dirname ${f}`
#     #Rename the output for a unique index
#     mv `dirname ${f}`/stat0.root `dirname ${f}`/stat_seedmod.root -f
    
#     if [[ "${energy}" == "1000" ]]
#     then
#         ${tracker_original} $f `dirname ${f}`
#         mv `dirname ${f}`/stat0.root `dirname ${f}`/stat_original.root -f    
#         ${tracker_noscatter} $f `dirname ${f}`
#         mv `dirname ${f}`/stat0.root `dirname ${f}`/stat_noscatter.root -f                
#     fi
#   done
  
  for f in ${DataDir}/muon_${energy}_GeV/*/*/run*.root; do
    ${tracker} $f `dirname ${f}`
    #Rename the output for a unique index
    mv `dirname ${f}`/stat0.root `dirname ${f}`/stat_seedmod.root -f
  done  
  
done

# for energy in ${EnergyList2[@]}
# do
#   for f in ${DataDir}/pion_${energy}_GeV/*/*/run*.root; do
#     ${tracker} $f `dirname ${f}`
#     #Rename the output for a unique index
#     mv `dirname ${f}`/stat0.root `dirname ${f}`/stat_seedmod.root -f
#   done

# done
