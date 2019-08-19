# helper functions for MS peack picking pipeline

# dilution filter function
# ds - MS feature intensity matrix (rows - samples, columns - features)
# dilFs - dilution factors of included dilution serie
dil_fct=function(ds, dilFs=c(1,2,4,8,16)){
  out=apply(ds, 2, function(x, dil=dilFs){
    cor(x, dil)
  })
  return(as.numeric(out))
}


# relative standard deviation function
# x - MS feature
rsd_fct=function(x){ sd(x, na.rm=T) / mean(x, na.rm=T) }




