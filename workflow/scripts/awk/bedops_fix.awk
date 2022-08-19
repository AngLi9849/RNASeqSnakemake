FNR > 5 { if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }
