#!/bin/bash
# run different versions of the coupled model (next: no, yes/no)

for mesh_name in 'cell_substrate_empty-nucleus' #'cell_substrate' 'lamellipodium' 'cell_substrate_nucleus'
do
    for coupling in 1 2 3 4 #1 2 3 4
    do
        for T in 10
        do
            for dt in 0.05
            do
                for C1 in 1 #0 1 0.5 2.0 
                do 
                    for E in 0.001 0.01 0.1 0.3 1 5.7 10 50 100 7000000 #0.001 0.01 0.1 0.3 1 5.7 10 50 100 7000000
                    do
                        for twoDstim in 'no' #'yes' 'no'
                        do
                            for partfixed in 'no' #'yes' 'no'
                            do 
                                for omega in 0.1 1 10 #0.1 1 10
                                do 
                                    python coupled_elas_rho_model_reduced_wFAK_empty-nucleus.py -T $T -dt $dt -C1 $C1 -bigE $E -mesh_name $mesh_name -coupling $coupling -twoDstim $twoDstim -partfixed $partfixed -omega $omega
                                done
                            done
                        done
                    done
                done
            done
        done
    done
done
