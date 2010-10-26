#PBS -l nodes=1
#PBS -l walltime=600:00:00
#PBS -N p45
#PBS -m ae
#PBS -M rahimian@gatech.edu
#PBS -r n

nohup nice matlab -nodesktop -nodisplay -nosplash -r "cd /nethome/arahimian3/Ves3D/matlab/; p = 48; precision = 'single'; Generate;" > /dev/null &

nohup nice matlab -nodesktop -nodisplay -nosplash -r "cd /nethome/arahimian3/Ves3D/matlab/; p = 96; precision = 'single'; Generate;" > /dev/null &

wait