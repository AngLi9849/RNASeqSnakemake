($7 == "gene" || $7 == "trscrpt") { print }

$6=="-" && ($7 != "gene" && $7 != "trscrpt") {
  if ($7 ~ "long_") {
    print
  } 
  else if (gene!=$4) {
    for (i in save) print save[i] ;
    delete save ;
    gene=$4 ; a=$3 ; b=$2 ;  n=1 ; v=1 ; k=$12 ; $13=k ; $12=n ; $14=v ; $8=$8""n ; $9=$9" "n ; save[1]=$0 ; 
#    if (k>1 && ($7=="exon" || $7=="intron")) { 
#      $7="first_"$7 ; print ;
#    } 
  }
  else if (gene==$4) {
    a=$3 ;  
    if (a<b) {
      for (i in save) print save[i] ;
      delete save ;
      n+=1 ; b=$2 ; v=1 ; k=$12 ; $13=k ; $12=n ; $14=v ; $8=$8""n ; $9=$9" "n ; save[1]=$0
    }
    else if (a>=b) {
      n+=0 ; b=b  ; k=k  ;
      if (k>$12 && ($7=="exon" || $7=="intron") ) {
        $7="long_"$7 ; $8=$4":"$7 ; $9=$8" "$7 ; $12=1 ; $13=1 ; $14=1 ; print
      }
      else if (k==$12) {
        v+=1 ;  $13=k ; $12=n ; $14=v ; $8=$8""n ; $9=$9" "n ; save[length(save)+1]=$0 ;
      }
      else if (k<$12) {
        load=$0 ;
        for (i in save) {
          $0=save[i] ; $7="long_"$7 ; $8=$4":"$7 ; $9=$8" "$7 ; $12=1 ; $13=1 ; $14=1 ; print ;
        } ;
        delete save ;
        $0 = load ;
        n+=0 ; b=$2 ; v=1 ; k=$12 ; $13=k ; $12=n ; $14=v ; $8=$8""n ; $9=$9" "n ; save[1]=$0 ;
      };
    } ;
#    if (k==1 && n>1 && ($7=="exon" || $7=="intron")) {
#      $7="last_"$7 ; print
#    } ;
  }
}

$6=="+" && ($7 != "gene" && $7 != "trscrpt") {
  if ($7 ~ "long_") {
    print
  }
  else if (gene!=$4) {
    for (i in save) print save[i] ;
    delete save ;
    gene=$4 ; c=$3 ; d=$2 ; w=$12 ; m=1 ; $13=m ; $8=$8""w ; $9=$9" "w ; save[1]=$0
 #   if (w>1 && ($7=="exon" || $7=="intron")) {
 #     $7="last_"$7 ; print ;
 #   }
  }
  else if (gene==$4) {
    c=$3 ;
   if (c<d) {
      for (i in save) print save[i] ;
      delete save ;
      m+=1 ;  w=$12 ; d=$2 ; $13=m ; $8=$8""w ; $9=$9" "w ; save[1]=$0
    }
    else if (c>=d) {
      m+=0 ; w=w ;
      if (w>$12 && ($7=="exon" || $7=="intron") ) {
        $7="long_"$7 ; $8=$4":"$7 ; $9=$8" "$7 ; $12=1 ; $13=1 ; $14=1 ; print
      }
      else if (w==$12) {
        w=$12 ; $13=m ; $8=$8""w ; $9=$9" "w ; save[length(save)+1]=$0
      }
      else if (w < $12) {
        load=$0 ; 
        for (i in save) {
          $0=save[i] ; $7="long_"$7 ; $8=$4":"$7 ; $9=$8" "$7 ; $12=1 ; $13=1 ; $14=1 ; print ;
        } ;
        delete save ;
        $0 = load ;
        m+=0 ; w=$12 ; d=$2 ; $13=m ; $8=$8""w ; $9=$9" "w ; save[1]=$0 ;
      } ;
    } ; 
#    if (w==1 && m>1 && ($7=="exon" || $7=="intron")) {
#      $7="first_"$7 ; print
#    } ;
  }
}

END {
  for (i in save) print save[i] ;
}
