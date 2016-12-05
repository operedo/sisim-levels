### Directivas para el gestor de colas (opcionales)
# Cambiar el nombre del trabajo
#$ -N testsisim-extrae
# Especificar un shell
#$ -S /bin/bash



#export OMP_NUM_THREADS=4
#export EXTRAE_CONFIG_FILE=/scratch/nas/2/$USER/sisimLevels/extrae.xml
#export EXTRAE_HOME=/Soft/extrae/3.0.2
#export LD_PRELOAD=${EXTRAE_HOME}/lib/libomptrace.so
#
#export OMP_SCHEDULE="static"
#WORKING_DIR=/scratch/nas/2/$USER/sisimLevels 
#cd ${WORKING_DIR}
#/scratch/nas/2/$USER/sisimLevels/sisimFortranLevelsPar.exe sisim-large.par > salida_levels_par_16.txt 2>&1 
#exit 0

#SIZE=huge
SIZE=large
WORKING_DIR=/home/$USER/dev/sisim-levels 
cd ${WORKING_DIR}

for intent in 1
do
export OMP_SCHEDULE="static"
for i in 16 8 4 2 1
#for i in 16
do
	export OMP_NUM_THREADS=${i}
	echo "THREADS=${i}"
	echo "INI_TIME=$(date +%s)"
	#/usr/bin/time ${WORKING_DIR}/sisimFortranLevelsPar.exe sisim-large.par > salida_levels_par_${i}.txt 2>&1
	/usr/bin/time ${WORKING_DIR}/sisimFortranLevelsPar.exe sisim-${SIZE}.par > salida_levels_par_${i}.txt 2>&1
	mv sisim.out sisim.out_par_${i}
	echo "FIN_TIME=$(date +%s)"

done

/usr/bin/time ${WORKING_DIR}/sisimFortranSeq.exe sisim-${SIZE}.par > salida_levels_gslib.txt 2>&1
mv sisim.out sisim.out_gslib
/usr/bin/time ${WORKING_DIR}/sisimFortranLevelsSeq.exe sisim-${SIZE}.par > salida_levels_seq.txt 2>&1
mv sisim.out sisim.out_seq



done

exit 0
