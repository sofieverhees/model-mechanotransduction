#!/bin/bash
# run different versions of the coupled model 

T = 10
dt = 0.05

for mesh_name in 'cell_substrate' 'lamellipodium' #'cell_substrate_empty-nucleus' 
do
    for coupling in 1 2 3 4 #1 2 3 4
    do
        for C1 in 1 #1 
        do 
            for E in 0.001 0.01 0.1 0.3 1 5.7 10 50 100 7000000 #0.001 0.01 0.1 0.3 1 5.7 10 50 100 7000000
            do
                for twoDstim in 'yes' 'no' #'yes' 'no'
                do
                    for partfixed in 'yes' 'no' #'yes' 'no'
                    do 
                        python main_coupled.py -T $T -dt $dt -C1 $C1 -bigE $E -mesh_name $mesh_name -coupling $coupling -twoDstim $twoDstim -partfixed $partfixed
                    done
                done
            done
        done
    done
    for coupling in 2 4 #2 4
    do
        for C1 in 0.5 2.0 #0.5 2.0 
        do 
            for E in 0.001 0.01 0.1 0.3 1 5.7 10 50 100 7000000 #0.001 0.01 0.1 0.3 1 5.7 10 50 100 7000000
            do
                for twoDstim in 'yes' 'no' #'yes' 'no'
                do
                    for partfixed in 'yes' 'no' #'yes' 'no'
                    do 
                        python main_coupled.py -T $T -dt $dt -C1 $C1 -bigE $E -mesh_name $mesh_name -coupling $coupling -twoDstim $twoDstim -partfixed $partfixed
                    done
                done
            done
        done
    done
done
