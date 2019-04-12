nb_args=$#
nb_variables=$1

variables_name=""
variables_value=""


# le nom du fichier est le dernier argument
# $@ permet de recuperer la liste complete des arguments
idx=1
for arg in $@
do

    idx_max1=$((1+$nb_variables))
    if [ $idx -le $idx_max1 ]
    then
	if [ $idx -gt 1 ]
	then
	    variables_name=$variables_name" "$arg
	fi
    fi

    if [ $idx -gt $idx_max1 ]
    then
	if [ $idx -lt $nb_args ]
	then
	    variables_value=$variables_value" "$arg
	fi
    fi

    if [ $idx = $nb_args ]
    then 
	filename=$arg
    fi

    idx=$((idx+1))

done

here=$PWD
cd ..
if [ -f $filename ]
then
    echo "file "$filename" already exist"
else
    echo "# "$variables_name" absolute_error  relative_error" > $filename
fi
cd $here
echo "python ToolBox.py add_error_in_function_of "$variables_value" "$filename > post_run



#for i in `seq 1 $nb_val` ;
