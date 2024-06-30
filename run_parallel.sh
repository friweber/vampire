for i in 0 1 2 3 4 5 6 7 8 9 10 11
do
	export TEMPERATURE_FILE=temp_int_scaled${i}.dat
	export OUTPUT_FILE=output_0.1_$i
	echo "TEMPERATURE_FILE is set to $TEMPERATURE_FILE"
	mpirun -np 10 ./vampire-parallel
done
