for ((i=64; i<=2080; i+=64))
do
	echo $i
	./matrix_vector $i $i
done
