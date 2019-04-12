liste=""
for rep in $(find . -maxdepth 1 -type d)
do
    if [ $rep != "." ] && [ $rep != "./python_solution_analytique" ] && [ $rep != "./.tmp" ] && [ $rep != "./variability" ] 
    then
	liste=$liste" "$rep
    fi 
done

python ToolBox.py AddNewTestCasesToPrm $liste

