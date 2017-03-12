#bash

echo -n "How many files to plot?"
read nfile

for((i=0;i<=$nfile;i++))
	do
		python $VORFLOW_DIR/source/tests/pickledPlot.py $i &
	done

