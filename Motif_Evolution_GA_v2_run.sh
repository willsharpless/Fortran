#! /bin/bash
#
for i in 678 799 1056 1189 1278
do 
    echo "Beginning Support Community Evolution v2"
    echo $(date)

    gfortran -g -c -fbacktrace -fcheck=all -w Motif_Evolution_GA_v2_${i}.f90 
    if [ $? -ne 0 ]; then
    echo "Compile error."
    exit
    fi
    #
    gfortran Motif_Evolution_GA_v2_${i}.o M_scramble.o rkf45_WASedit.o -framework Accelerate
    if [ $? -ne 0 ]; then
    echo "Load error."
    exit
    fi
    rm Motif_Evolution_GA_v2_${i}.o
    #
    mv a.out Motif_Evolution_GA_v2_${i}
    # ./Motif_Evolution_GA_v2_${i}
    # exit

    for ii in 1 2 3
    do
        echo "Beginning run " ${ii} " with support community " ${i} " at " $(date)
        ./Motif_Evolution_GA_v2_${i} > Evolution_runs/Motif_Evolution_v2_${i}_300its_run${ii}_011021.txt
        echo "Completed run " ${ii} " with support community " ${i} " at " $(date)
    done
done
