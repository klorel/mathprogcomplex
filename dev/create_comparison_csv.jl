

function create_csv_scaling(folder)
    csv1 = readdlm(joinpath("..","knitro_runs","scaling0_nb_iter_$(folder).csv"), ';')
    csv2 = readdlm(joinpath("..","knitro_runs","scaling1_nb_iter_$(folder).csv"), ';')
    csv3 = readdlm(joinpath("..","knitro_runs","scaling2_nb_iter_$(folder).csv"), ';')
    csv4 = readdlm(joinpath("..","knitro_runs","scaling3_nb_iter_$(folder).csv"), ';')

    f = open(joinpath("..", "knitro_runs","comparaisons_scaling_$folder.csv"),"w")
    write(f, "Step 1 scaling 0;Step 2 scaling 0;Step 3 scaling 0;Step 1 scaling 1;Step 2 scaling 1;Step 3 scaling 1;Step 1 scaling 2;Step 2 scaling 2;Step 3 scaling 2;Step 1 scaling 3;Step 2 scaling 3;Step 3 scaling 3; Step 1 best; Step 2 best;Step 3 best\n ")

    for i in 2:size(csv1,1)
        min1 = min(csv1[i,1], csv2[i,1], csv3[i,1], csv4[i,1])
        min2 = min(csv1[i,2], csv2[i,2], csv3[i,2], csv4[i,2])
        min3 = min(csv1[i,3], csv2[i,3], csv3[i,3], csv4[i,3])
        write(f, "$(csv1[i,1]);$(csv1[i,2]);$(csv1[i,3]);$(csv2[i,1]);$(csv2[i,2]);$(csv2[i,3]);$(csv3[i,1]);$(csv3[i,2]);$(csv3[i,3]);$(csv4[i,1]);$(csv4[i,2]);$(csv4[i,3]);$min1;$min2;$min3 \n")
    end

end

create_csv_scaling("Phase_0_IEEE14")
