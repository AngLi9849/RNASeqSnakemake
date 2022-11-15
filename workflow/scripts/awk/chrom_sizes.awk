$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($1,2,100)"\t"; } $0 !~ ">" {c+=length($0);} END { print c; }

