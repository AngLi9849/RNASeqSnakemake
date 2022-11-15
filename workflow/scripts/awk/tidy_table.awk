{for (i=1;i<=NF;i++) {if(i==1) {printf "%-100s\t",$i} else if(i<NF) {printf"%-20s\t",$i} else {printf "%-20s\n",$i}}}
