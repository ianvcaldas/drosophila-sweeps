
for slimscript in $(ls ./*.slim)
do
    for n in {1..300}
    do
        echo "${slimscript} ${n}"
        ../../bin/slim3.7 ${slimscript} &> /dev/null
    done
done
