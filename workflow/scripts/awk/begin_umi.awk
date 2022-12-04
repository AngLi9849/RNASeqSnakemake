$0~/^@/ {
  split($1,a,":") ;
  if (gensub("@","","g",a[1]) ~ /[^ACTGN]/) {
    print
  } else {
    id=$1 ;
    match(a[1],/@([ACTGN]*)/,b) ;
    $1=gensub(b[1]":","",1,id)"_"b[1] ;
    print
  }
}

$0 !~ /^@/ {
  print
}

