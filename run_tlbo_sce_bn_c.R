source("tlbo_full_vec_sce_bn_c_multi.R")
source('drutes_eval_inverse_multi.R')


seed=sample.int(1e+8,size=1)
set.seed(seed)
write(seed,"rand_seed.out")

maxeval=1# number of parallel oeprataions 
system(paste('./opti_setup.sh', maxeval))
pop=1 #population

#parameters
alphamat1min=0.005
alphamat1max=0.05
# nmat1=$3
nmat1min=1.3
nmat1max=2.3
# #thetamat1 based on thetas=0.457
# #anomaly will also be assigned to fracture
# alphamat2=$4
alphamat2min=0.005
alphamat2max=0.07
# nmat2=$5
nmat2min=1.2
nmat2max=3
# thsmat2=$6
thsmat2min=0.68 #0.7
thsmat2max=0.72 # 0.9  
# # lower soil will also be assigned tof racture
# alphamat3=$7
alphamat3min=0.005
alphamat3max=0.05
# nmat3=$8
nmat3min=1.3
nmat3max=2.3
# conductivity Ksmat2 and Ksmat3 will also be assigned to fracture medium
# Ksmat1=${9}
Ksmat1min=0.1 #cm/h
Ksmat1max=5 #2
# Ksmat2=${10}
Ksmat2min=0.1 #cm/h
Ksmat2max=5
# Ksmat3=${11}
Ksmat3min=0.1 #cm/h
Ksmat3max=10
# # dual fracture at the top
# alphafrac1=${12}
alphafrac1min=0.05
alphafrac1max=0.15
# nfrac1=${13}
nfrac1min=2
nfrac1max=3
# thsfrac1=${14}
thsfrac1min=0.3
thsfrac1max=0.5
# Ksfrac1=${15}
Ksfrac1min=15 #cm/h
Ksfrac1max=30
# # exchange values alpha. Beta assumed to be 0.4 and gamma 15.0.
# aex=${16}
aexmin=0.5
aexmax=1.5
# wfrac=${17}
wfracmin=0.1
wfracmax=0.3
#
Ksexmin=1e-7
Ksexmax=1e-4
#{19}
winfmin=0.5
winfmax=0.95

#{20}
alphamat4min=0.005
alphamat4max=0.07
#{21}
nmat4min=1.2
nmat4max=3
#{22}
thsmat4min=0.45 #0.7
thsmat4max=0.52 # 0.9
#{23}
Ksmat4min=0.1 #cm/h
Ksmat4max=5

xmin=c(alphamat1min,nmat1min,alphamat2min,nmat2min,thsmat2min,alphamat3min,nmat3min,
       Ksmat1min,Ksmat2min,Ksmat3min,alphafrac1min,nfrac1min,thsfrac1min,Ksfrac1min,
       aexmin,wfracmin,Ksexmin,winfmin,alphamat4min,nmat4min,thsmat4min,Ksmat4min)

xmax=c(alphamat1max,nmat1max,alphamat2max,nmat2max,thsmat2max,alphamat3max,nmat3max,
       Ksmat1max,Ksmat2max,Ksmat3max,alphafrac1max,nfrac1max,thsfrac1max,Ksfrac1max,
       aexmax,wfracmax,Ksexmax,winfmax,alphamat4max,nmat4max,thsmat4max,Ksmat4max)
print(xmax)
TLBO=TLBO_sce_bn_fv_R(pop,2,length(xmin),xmin=xmin,xmax=xmax,1000,F,maxeval,0.8)
write.csv(TLBO,paste('inverse_tlbo_sce_bn.csv',sep=''))
