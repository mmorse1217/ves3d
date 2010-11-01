#PBS -l nodes=1
#PBS -l walltime=600:00:00
#PBS -N GenerateMats
#PBS -m ae
#PBS -M rahimian@gatech.edu
#PBS -r n

matlab -nodesktop -nodisplay -nosplash -r "cd /nethome/arahimian3/Ves3D/matlab/; Generate;" > /dev/null &

wait