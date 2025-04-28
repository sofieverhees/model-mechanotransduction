#!/bin/bash
# run different versions of the coupled model 

T=10
dt=0.05

for mesh_name in 'cell_substrate' #'cell_substrate' 'lamellipodium'
do
    for coupling in 1 2 3 4 #1 2 3 4
    do
        for C1 in 1 #1 
        do 
            for E in 0.1 5.7 7000000 #0.001 0.01 0.1 0.3 1 5.7 10 50 100 7000000
            do
                for twoDstim in 'no' #'yes' 'no'
                do
                    for partfixed in 'no' #'yes' 'no'
                    do 
                        for th_lamb in 0.1 1 10 #0.1 1 10 
                        do
                            for th_mu in 0.1 1 10 #0.1 1 10 
                            do 
                                python main_viscoelastic.py -T $T -dt $dt -C1 $C1 -bigE $E -mesh_name $mesh_name -coupling $coupling -twoDstim $twoDstim -partfixed $partfixed -th_lamb $th_lamb -th_mu $th_mu
                            done
                        done
                    done
                done
            done
        done
    done
done
