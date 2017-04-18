data=read.table("drutes.conf/REdual/hinim.in")
bc=-100
for(i in 1:(length(data[,2])-1)){
  if(data[i,2]<=bc && data[i+1,2]>bc){
    grad = (data[i+1,3] - data[i,3])/(data[i+1,2] - data[i,2])
    value =  data[i,3] + grad*(bc - data[i,2])
    break
  }
}

write(value,"hbc.in")
