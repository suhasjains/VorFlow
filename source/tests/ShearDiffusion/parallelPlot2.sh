#bash

echo -n "How many files to plot? "
read nfile

for((i=0;i<=$nfile;i++))
	do
		python pickledPlot2.py $i
	done

