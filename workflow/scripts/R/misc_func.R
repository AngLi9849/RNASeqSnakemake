# common_prefix: Finds the common prefixing string elements and connect them by space from a vector/list of strings (str_ls) with elements separated by specified delimiter (split)
# Example: Returns "TEST SUBJECT" with common_prefix( str_ls = c("TEST_SUBJECT_1","TEST_SUBJECT_2","TEST_SUBJECT_3"), split = "_" )

common_prefix <- function(str_ls,split) {
  str_ls_split <- sapply(
    str_ls,function(x) { strsplit(x,split)}
  )
  str_ls_min <- min(sapply(str_ls_split,length))
  min_str_ls <- data.frame(
    sapply(
      str_ls_split,function(x){
        x[1:str_ls_min]
      }))
  prefix_n <- min(sapply(c(1:str_ls_min), function(x) {
    ifelse(length(unique(c(min_str_ls[x,])))!=1,x,Inf)
  }))-1
  prefix <- ifelse(
    prefix_n<=0,
    "",
    paste(str_ls_split[[1]][1:prefix_n],collapse=" ")
  )
  return(prefix)
}


# unique_suffix: Return a list of space-separated unique suffixing string elements from a vector/list of strings (str_ls) with elements separated by specified delimiter (split)
# Example: Returns c("1 11","2 22","3 33") with unique_suffix( str_ls = c("TEST_SUBJECT_1_11","TEST_SUBJECT_2_22","TEST_SUBJECT_3_33"), split = "_" )

unique_suffix <- function(str_ls,split) {
  str_ls_split <- sapply(
    str_ls,function(x) { strsplit(x,split)}
  )
  str_ls_min <- min(sapply(str_ls_split,length))
  min_str_ls <- data.frame(
    sapply(
      str_ls_split,function(x){
        x[1:str_ls_min]
      }))
  suffix_n <- min(sapply(c(1:str_ls_min), function(x) {
    ifelse(length(unique(c(min_str_ls[x,])))!=1,x,Inf)
  }))
  suffix_ls <- sapply(str_ls_split, function(x){
    ifelse(
      length(x) <= suffix_n-1, 
      "", 
      paste(x[suffix_n:length(x)],collapse=" ")
    )
  }) 
  return(suffix_ls) 
}

