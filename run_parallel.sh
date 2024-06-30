for i in 0 1 2 3 4 5
do
	export TEMPERATURE_FILE=temp_int_scaled7.dat
	export ALPHA = i*0.2+0.1
	export OUTPUT_FILE=output_0.1_$i
	echo "TEMPERATURE_FILE is set to $TEMPERATURE_FILE"
	mpirun -np 4 -oversubscribe ./vampire-parallel
done
