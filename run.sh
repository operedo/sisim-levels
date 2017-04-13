#for domain in 800x800x160 400x800x160
for domain in 420x600x400-10
do

for neig in 16 32 
do

for proc in 1 2 4 8 16
do
	export OMP_NUM_THREADS=$proc
	/usr/bin/time ./sisimFortranLevelsPar.exe  sisim-extreme-${neig}-${domain}.par > salida_par_${proc}.txt 2>&1 
	grep 'TIME\|elapsed' salida_par_${proc}.txt > salida_par_${proc}_${neig}_${domain}.txt
done

#/usr/bin/time ./sgsimFortranLevelsSeq.exe  sgsim-extreme-${neig}-600x600x160.par > salida_seq.txt 2>&1 
#grep 'TIME\|elapsed' salida_seq.txt > salida_seq_${neig}.txt

/usr/bin/time ./sisimFortranSeq.exe  sisim-extreme-${neig}-${domain}.par > salida_gslib.txt 2>&1 
grep 'TIME\|elapsed' salida_gslib.txt > salida_gslib_${neig}_${domain}.txt

done

#for neig in 32 16
#do
#
#for proc in 16 8 4 2 1
#do
#	export OMP_NUM_THREADS=$proc
#	/usr/bin/time ./sgsimFortranLevelsPar.exe  sgsim-extreme-${neig}-${domain}.par > salida_par_${proc}.txt 2>&1 
#	grep 'TIME\|elapsed' salida_par_${proc}.txt > salida_par_${proc}_${neig}_${domain}.txt
#done
#
##/usr/bin/time ./sgsimFortranLevelsSeq.exe  sgsim-extreme-${neig}-600x600x160.par > salida_seq.txt 2>&1 
##grep 'TIME\|elapsed' salida_seq.txt > salida_seq_${neig}.txt
#
#/usr/bin/time ./sgsimFortranSeq.exe  sgsim-extreme-${neig}-${domain}.par > salida_gslib.txt 2>&1 
#grep 'TIME\|elapsed' salida_gslib.txt > salida_gslib_${neig}_${domain}.txt
#
#done


done


for domain in 210x600x400-10
do

for neig in 16 32 64 128 
do

for proc in 1 2 4 8 16
do
	export OMP_NUM_THREADS=$proc
	/usr/bin/time ./sisimFortranLevelsPar.exe  sisim-extreme-${neig}-${domain}.par > salida_par_${proc}.txt 2>&1 
	grep 'TIME\|elapsed' salida_par_${proc}.txt > salida_par_${proc}_${neig}_${domain}.txt
done

#/usr/bin/time ./sgsimFortranLevelsSeq.exe  sgsim-extreme-${neig}-600x600x160.par > salida_seq.txt 2>&1 
#grep 'TIME\|elapsed' salida_seq.txt > salida_seq_${neig}.txt

/usr/bin/time ./sisimFortranSeq.exe  sisim-extreme-${neig}-${domain}.par > salida_gslib.txt 2>&1 
grep 'TIME\|elapsed' salida_gslib.txt > salida_gslib_${neig}_${domain}.txt

done

done

