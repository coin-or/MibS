BEGIN {
   solver = "MIBS";
   experiment_num = 0;
   counter = -1;
   num_instances = 0;
}


/LL Data File:/ {
   instance_name = $4
      instance[num_instances] = instance_name;
      counter = num_instances;
      num_instances++;
}

/nodes processed:/ {
   nodesPro[experiment_num, counter] = $6;
}

/fully processed:/ {
   nodesPro[experiment_num, counter] = $7;
}

/CPU time:/{
    timeCpu[experiment_num, counter] = $5;
}

/wall-clock/ {
    timeWall[experiment_num, counter] = $5;
}

/gap/ {
   gap[experiment_num, counter] = $6;
   if (gap[experiment_num, counter] > .0001){
      nodesPro[experiment_num, counter] = -nodesPro[experiment_num, counter];
      timeCpu[experiment_num, counter] = -timeCpu[experiment_num, counter];
   }else{
      good[counter] = 1;
   }
}

/feasibility/ {
   feas_time[experiment_num, counter] = $5
}


/Cost/ {
    cost[experiment_num, counter] = $3
}

/VF) =/ {
    timeCpuVF[experiment_num, counter] = $7;
}

/UB) =/ {
    timeCpuUB[experiment_num, counter] = $7;
}
    

/ALPS did not find a solution/{
    nodesPro[experiment_num, counter] = "-0.01";
    timeCpu[experiment_num, counter] = "-0.01";
    cost[experiment_num, counter] = "-0.01";
}

END {
    file = sprintf("%s.summary", "UBTime");
    for (j = 0; j < num_instances; j++){

	              #printf("%-30s %7d   %7.2f   %7.2f\n", instance[j], nodes[0, j], cost[0, j], time[0,j]) >> file;
	printf("%7.3f\n", timeCpuUB[0,j]) >> file;
    }
}
