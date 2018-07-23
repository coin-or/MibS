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

/Number of integer UL Variables:/ {
   numULIntVar[experiment_num, counter] = $6;
}

/Number of integer LL Variables:/ {
    numLLIntVar[experiment_num, counter] = $6;
}

END {
    file1 = sprintf("%s.summary", "r1LeqR2List");
    file2 = sprintf("%s.summary", "r1GeR2List");
    for (j = 0; j < num_instances; j++){
	if (numULIntVar[0,j] <= numLLIntVar[0,j]) {
	    printf("%-30s\n", instance[j]) >> file1;
       }
	if (numULIntVar[0,j] > numLLIntVar[0,j]) {
	    printf("%-30s\n", instance[j]) >> file2;
	}
    }
}
