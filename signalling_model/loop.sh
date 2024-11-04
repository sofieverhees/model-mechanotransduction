#!/bin/bash
# run different versions of the coupled model

for mesh_name in 'cell_substrate'
do
    for T in 10
    do
        for dt in 0.05
        do
            for D1 in 40 100 #40 100
            do
                for E in 0.001 0.01 0.1 0.3 1 5.7 10 50 100 7000000 #0.001 0.01 0.1 0.3 1 5.7 10 50 100 7000000
                do
                    for twoDstim in '2D' '2xD' '3D' #'2D' '2xD' '3D' 
                    do 
                        python rho_model_reduced_wFAK.py -T $T -dt $dt -E $E -D1 $D1 -mesh_name $mesh_name -twoDstim $twoDstim
                    done
                done
            done
        done
    done
done
