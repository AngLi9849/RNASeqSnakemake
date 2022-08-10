cat resources/annotations/ensembl_homo_sapiens.GRCh38.104_genome.gtf | awk -F'\t' 'match($0,/transcript_name "ACTB-[^"]*"/,a) && $3=="transcript" && match($0,/transcript_biotype "([^"]*)"/,b) && match($0,/transcript_name "([^"]*)"/,n) && match($0,/tag "([^"]*)"/,tag) && match($0,/transcript_id "([^"]*)"/,id) {if ($0 ~ "transcript_support_level" ) {match($0,/transcript_support_level "([^"]).*"/,t)} else {t[1]=6} ; print n[1],tag[1],t[1], b[1], id[1]}' | sort | uniq