source /cvmfs/larsoft.opensciencegrid.org/products/setup
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh
setup root v6_22_06a -q e19:p383b:prof
setup jobsub_client v1_3_3
export JOBSUB_GROUP=minerva
source /cvmfs/minerva.opensciencegrid.org/minerva/hep_hpc_products/setups
setup cmake v3_9_5
