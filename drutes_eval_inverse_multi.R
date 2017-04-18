# read reference when source is called
ref1=read.table("ref1_2D.in", quote="\"")
timeref=ref1[,1]
ref1=ref1[,2]
ref2=read.table("ref2_2D.in", quote="\"")
ref2=ref2[,2]
ref3=read.table("ref3_2D.in", quote="\"")
ref3=ref3[,2]


eval_fun= function(ln_id,pop,gens,weight,output=T){
  system2("./drut_opti_c.sh",wait=T) # change according to variant
  waitallfiles=F
  while(waitallfiles){
    boolwaitvec=c()
    for(i in 1:ln_id){
      file=paste(i,'/out/Re_dual_totH_f_theta_f-14.out',sep='')
      boolwaitvec[i]=file.exists(file)
      if(!boolwaitvec[i]){
        file=paste(i,'/out/Re_dual_totH_f_theta_f-14.sci',sep='')
        boolwaitvec[i]=file.exists(file)
      }
      if(!boolwaitvec[i]){
        file=paste(i,'/out/Re_dual_totH_f_theta_f-14.msh',sep='')
        boolwaitvec[i]=file.exists(file)
      }
#      print('in bool loop')
    }
    waitallfiles=any(!boolwaitvec)
  }
  
  result=matrix(ncol=13,nrow=ln_id)
  print('done with drutes')
  for(i in 1:ln_id){
    t=0
    told=1e30
    # TDR 1 with 6 obs point around TDR 1
    while(t!=told){
      told=t
      obs1=read.table(paste(i,"/out/obspt_Re_dual_totH_f-1.out",sep=''), quote="\"",comment.char="#", sep="")
      timeobs1=obs1[,1]
      t=length(timeobs1)
    }
    obs1=obs1[,3]
    obs1m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-1.out",sep=''), quote="\"",comment.char="#", sep="")
    obs1m=obs1m[,3]
    obs1=obs1*weight[i]+obs1m*(1-weight[i])
    #
    obs2=read.table(paste(i,"/out/obspt_Re_dual_totH_f-2.out",sep=''), quote="\"",comment.char="#", sep="")
    obs2=obs2[,3]
    obs2m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-2.out",sep=''), quote="\"",comment.char="#", sep="")
    obs2m=obs2m[,3]
    obs2=obs2*weight[i]+obs2m*(1-weight[i])
    #
    obs3=read.table(paste(i,"/out/obspt_Re_dual_totH_f-3.out",sep=''), quote="\"",comment.char="#", sep="")
    obs3=obs3[,3]
    obs3m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-3.out",sep=''), quote="\"",comment.char="#", sep="")
    obs3m=obs3m[,3]
    obs3=obs3*weight[i]+obs3m*(1-weight[i])
    #
    obs4=read.table(paste(i,"/out/obspt_Re_dual_totH_f-4.out",sep=''), quote="\"",comment.char="#", sep="")
    obs4=obs4[,3]
    obs4m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-4.out",sep=''), quote="\"",comment.char="#", sep="")
    obs4m=obs4m[,3]
    obs4=obs4*weight[i]+obs4m*(1-weight[i])
    #
    obs5=read.table(paste(i,"/out/obspt_Re_dual_totH_f-5.out",sep=''), quote="\"",comment.char="#", sep="")
    obs5=obs5[,3]
    obs5m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-5.out",sep=''), quote="\"",comment.char="#", sep="")
    obs5m=obs5m[,3]
    obs5=obs5*weight[i]+obs5m*(1-weight[i])
    #
    obs6=read.table(paste(i,"/out/obspt_Re_dual_totH_f-6.out",sep=''), quote="\"",comment.char="#", sep="")
    obs6=obs6[,3]
    obs6m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-6.out",sep=''), quote="\"",comment.char="#", sep="")
    obs6m=obs6m[,3]
    obs6=obs6*weight[i]+obs6m*(1-weight[i])
    #
    print('read all obs file for TDR1')
    TDR1sim=colMeans(rbind(obs1,obs2,obs3,obs4,obs5,obs6))
    # TDR 2 with 6 obs point around TDR 2
    obs7=read.table(paste(i,"/out/obspt_Re_dual_totH_f-7.out",sep=''), quote="\"",comment.char="#", sep="")
    obs7=obs7[,3]
    obs7m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-7.out",sep=''), quote="\"",comment.char="#", sep="")
    obs7m=obs7m[,3]
    obs7=obs7*weight[i]+obs7m*(1-weight[i])
    #
    obs8=read.table(paste(i,"/out/obspt_Re_dual_totH_f-8.out",sep=''), quote="\"",comment.char="#", sep="")
    obs8=obs8[,3]
    obs8m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-8.out",sep=''), quote="\"",comment.char="#", sep="")
    obs8m=obs8m[,3]
    obs8=obs8*weight[i]+obs8m*(1-weight[i])
    #
    obs9=read.table(paste(i,"/out/obspt_Re_dual_totH_f-9.out",sep=''), quote="\"",comment.char="#", sep="")
    obs9=obs9[,3]
    obs9m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-9.out",sep=''), quote="\"",comment.char="#", sep="")
    obs9m=obs9m[,3]
    obs9=obs9*weight[i]+obs9m*(1-weight[i])
    #
    obs10=read.table(paste(i,"/out/obspt_Re_dual_totH_f-10.out",sep=''), quote="\"",comment.char="#", sep="")
    obs10=obs10[,3]
    obs10m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-10.out",sep=''), quote="\"",comment.char="#", sep="")
    obs10m=obs10m[,3]
    obs10=obs10*weight[i]+obs10m*(1-weight[i])
    #
    obs11=read.table(paste(i,"/out/obspt_Re_dual_totH_f-11.out",sep=''), quote="\"",comment.char="#", sep="")
    obs11=obs11[,3]
    obs11m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-11.out",sep=''), quote="\"",comment.char="#", sep="")
    obs11m=obs11m[,3]
    obs11=obs11*weight[i]+obs11m*(1-weight[i])
    #
    obs12=read.table(paste(i,"/out/obspt_Re_dual_totH_f-12.out",sep=''), quote="\"",comment.char="#", sep="")
    obs12=obs12[,3]
    obs12m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-12.out",sep=''), quote="\"",comment.char="#", sep="")
    obs12m=obs12m[,3]
    obs12=obs12*weight[i]+obs12m*(1-weight[i])
    #
    TDR2sim=colMeans(rbind(obs7,obs8,obs9,obs10,obs11,obs12))
    # TDR 3 with 6 obs point around TDR 3
    #overwriting above to limit memory use
    obs13=read.table(paste(i,"/out/obspt_Re_dual_totH_f-13.out",sep=''), quote="\"",comment.char="#", sep="")
    obs13=obs13[,3]
    obs13m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-13.out",sep=''), quote="\"",comment.char="#", sep="")
    obs13m=obs13m[,3]
    obs13=obs13*weight[i]+obs13m*(1-weight[i])
    #
    obs14=read.table(paste(i,"/out/obspt_Re_dual_totH_f-14.out",sep=''), quote="\"",comment.char="#", sep="")
    obs14=obs14[,3]
    obs14m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-14.out",sep=''), quote="\"",comment.char="#", sep="")
    obs14m=obs14m[,3]
    obs14=obs14*weight[i]+obs14m*(1-weight[i])
    #
    obs15=read.table(paste(i,"/out/obspt_Re_dual_totH_f-15.out",sep=''), quote="\"",comment.char="#", sep="")
    obs15=obs15[,3]
    obs15m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-15.out",sep=''), quote="\"",comment.char="#", sep="")
    obs15m=obs15m[,3]
    obs15=obs15*weight[i]+obs15m*(1-weight[i])
    #
    obs16=read.table(paste(i,"/out/obspt_Re_dual_totH_f-16.out",sep=''), quote="\"",comment.char="#", sep="")
    obs16=obs16[,3]
    obs16m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-16.out",sep=''), quote="\"",comment.char="#", sep="")
    obs16m=obs16m[,3]
    obs16=obs16*weight[i]+obs16m*(1-weight[i])
    #
    obs17=read.table(paste(i,"/out/obspt_Re_dual_totH_f-17.out",sep=''), quote="\"",comment.char="#", sep="")
    obs17=obs17[,3]
    obs17m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-17.out",sep=''), quote="\"",comment.char="#", sep="")
    obs17m=obs17m[,3]
    obs17=obs17*weight[i]+obs17m*(1-weight[i])
    #
    obs18=read.table(paste(i,"/out/obspt_Re_dual_totH_f-18.out",sep=''), quote="\"",comment.char="#", sep="")
    obs18=obs18[,3]
    obs18m=read.table(paste(i,"/out/obspt_Re_dual_totH_m-18.out",sep=''), quote="\"",comment.char="#", sep="")
    obs18m=obs18m[,3]
    obs18=obs18*weight[i]+obs18m*(1-weight[i])
    #
    TDR3sim=colMeans(rbind(obs13,obs14,obs15,obs16,obs17,obs18))
    print('result time')
    any(is.na(c(ref1,ref2,ref3)))
    any(is.na(timeref))
    any(is.na(c(TDR1sim,TDR2sim,TDR3sim)))
    any(is.na(timeobs1))
    res1=log10(RMSE(ref1,timeref,TDR1sim,timeobs1))
    res2=log10(RMSE(ref2,timeref,TDR2sim,timeobs1))
    res3=log10(RMSE(ref3,timeref,TDR3sim,timeobs1))
    any(is.na(c(res1,res2,res3)))
    result[i,1]=mean(c(res1,res2,res3))
    result[i,2]=res1
    result[i,3]=res2
    result[i,4]=res3
    result[i,5]=cor_p(ref1,timeref,TDR1sim,timeobs1)
    result[i,6]=cor_p(ref2,timeref,TDR2sim,timeobs1)
    result[i,7]=cor_p(ref3,timeref,TDR3sim,timeobs1)
    result[i,8]=means(ref1,timeref,TDR1sim,timeobs1)/mean(ref1)
    result[i,9]=means(ref2,timeref,TDR2sim,timeobs1)/mean(ref2)
    result[i,10]=means(ref3,timeref,TDR3sim,timeobs1)/mean(ref3)
    result[i,11]=vars(ref1,timeref,TDR1sim,timeobs1)/var(ref1)
    result[i,12]=vars(ref2,timeref,TDR2sim,timeobs1)/var(ref2)
    result[i,13]=vars(ref3,timeref,TDR3sim,timeobs1)/var(ref3)

    if(output){
      ths1=approxTDR(ref1,timeref,TDR1sim,timeobs1)
      ths2=approxTDR(ref2,timeref,TDR2sim,timeobs1)
      ths3=approxTDR(ref3,timeref,TDR3sim,timeobs1)
      write.csv(ths1,paste(pop[i],"_gen",gens,"TDR1sim.csv",sep=''))
      write.csv(ths2,paste(pop[i],"_gen",gens,"TDR2sim.csv",sep=''))
      write.csv(ths3,paste(pop[i],"_gen",gens,"TDR3sim.csv",sep=''))
    }
  }
  print('done eval')
  return(result)#
}

RMSE=function(real,timeref,sim, simtime){
  ln_time=length(timeref)-1
  newsim=approx(y=sim,x=simtime,xout=timeref[2:length(timeref)])
  result=sqrt(sum((c(sim[1],newsim$y)-real)^2)/ln_time)
  return(result)
}

cor_p=function(real,timeref,sim, simtime){
  ln_time=length(timeref)-1
  newsim=approx(y=sim,x=simtime,xout=timeref[2:length(timeref)])
  result=cor(c(sim[1],newsim$y),real,method="pearson")
  return(result)
}

means=function(real,timeref,sim, simtime){
  ln_time=length(timeref)-1
  newsim=approx(y=sim,x=simtime,xout=timeref[2:length(timeref)])
  result=mean(c(sim[1],newsim$y))
  return(result)
}

vars=function(real,timeref,sim, simtime){
  ln_time=length(timeref)-1
  newsim=approx(y=sim,x=simtime,xout=timeref[2:length(timeref)])
  result=var(c(sim[1],newsim$y))
  return(result)
}

approxTDR=function(real,timeref,sim, simtime){
  ln_time=length(timeref)-1
  newsim=approx(y=sim,x=simtime,xout=timeref[2:length(timeref)])
  result=c(sim[1],newsim$y)
  return(result)
}