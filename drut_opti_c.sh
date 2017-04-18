#!/bin/bash

function run_drutes {

  cd $1
  #upper layer
  alphamat1=$2
  nmat1=$3
#thetamat1 based on thetas=0.457
  #anomaly will also be assigned to fracture
  alphamat2=$4
  nmat2=$5
  thsmat2=$6
  # lower soil will also be assigned tof racture
  alphamat3=$7
  nmat3=$8
  # conductivity Ksmat2 and Ksmat3 will also be assigned to fracture medium
  Ksmat1=${9}
  Ksmat2=${10}
  Ksmat3=${11}
  # dual fracture at the top
  alphafrac1=${12}
  nfrac1=${13}
  thsfrac1=${14}
  Ksfrac1=${15}
  # exchange values alpha. Beta assumed to be 0.4 and gamma 15.0.
  aex=${16}
  wfrac=${17}
  Ksex=${18}
  winf=${19}
  echo $winf 
  alphamat4=${20}
  nmat4=${21}
  thsmat4=${22}
  Ksmat4=${23}
  # we assume thetatotal=0.457, therefore we can derive the matrix water content as follows
  thsmat1=`echo "scale=12; (0.457-$wfrac*$thsfrac1)/(1-$wfrac)" | bc`
#   substitution of parameters into input files for drutes 
  sed -e 's/!alphamat1/'$alphamat1'/g' -e 's/!nmat1/'$nmat1'/g' -e 's/!thsmat1/'$thsmat1'/g' -e 's/!alphamat3/'$alphamat3'/g' -e 's/!nmat3/'$nmat3'/g' -e 's/!Ksmat1/'$Ksmat1'/g' -e 's/!Ksmat3/'$Ksmat3'/g' -e 's/!alphafrac/'$alphafrac1'/g' -e 's/!nfrac/'$nfrac1'/g' -e 's/!thsfrac/'$thsfrac1'/g' -e 's/!Ksfrac/'$Ksfrac1'/g' -e 's/!aex/'$aex'/g' -e 's/!wfrac/'$wfrac'/g' -e 's/!Ksex/'$Ksex'/g' drutes.conf/REdual/dualtempconst1D.conf > drutes.conf/REdual/dual.conf
  cp drutes.conf/global1D_2w.conf drutes.conf/global.conf
  cp drutes.conf/REdual/pet1D_2w.in drutes.conf/REdual/102.bc
  sed -e 's/!winf/'$winf'/g' drutes.conf/REdual/dual_bc1D.conf > drutes.conf/REdual/dual_bc.conf
  nice -n19 bin/drutes > /dev/null
  cp out/Re_dual_totH_m_total_head_m-10.dat drutes.conf/REdual/hinim.in
  cp out/Re_dual_totH_f_total_head_f-10.dat drutes.conf/REdual/hinif.in
  nice -n19 Rscript hini.R
  read hini < hbc.in
  sed -e 's/!winf/'$winf'/g' -e 's/!hini/'$hini'/g' drutes.conf/REdual/dual_bc2D.conf > drutes.conf/REdual/dual_bc.conf
  cp drutes.conf/global2D.conf drutes.conf/global.conf
  sed -e 's/!alphamat1/'$alphamat1'/g' -e 's/!nmat1/'$nmat1'/g' -e 's/!thsmat1/'$thsmat1'/g' -e 's/!alphamat2/'$alphamat2'/g' -e 's/!nmat2/'$nmat2'/g' -e 's/!thsmat2/'$thsmat2'/g' -e 's/!alphamat3/'$alphamat3'/g' -e 's/!nmat3/'$nmat3'/g' -e 's/!Ksmat1/'$Ksmat1'/g' -e 's/!Ksmat2/'$Ksmat2'/g' -e 's/!Ksmat3/'$Ksmat3'/g' -e 's/!Ksmat4/'$Ksmat4'/g' -e 's/!alphafrac/'$alphafrac1'/g' -e 's/!nfrac/'$nfrac1'/g' -e 's/!thsfrac/'$thsfrac1'/g' -e 's/!Ksfrac/'$Ksfrac1'/g' -e 's/!aex/'$aex'/g' -e 's/!wfrac/'$wfrac'/g' -e 's/!Ksex/'$Ksex'/g' -e 's/!alphamat4/'$alphamat4'/g' -e 's/!nmat4/'$nmat4'/g' -e 's/!thsmat4/'$thsmat4'/g' drutes.conf/REdual/dualtempconst2D.conf > drutes.conf/REdual/dual.conf
  cp drutes.conf/REdual/pet2D.in drutes.conf/REdual/101.bc
  nice -n19 bin/drutes #> /dev/null
  echo 'process': $1 'done' > objfnc.val
  cd ..
}


rm -f drutes.vals

#count the number of processes       
let nproc=0
while read l a b c d e f g h i j k m n o p q r s t u v w; do
    if [[  $l == "p"  ]]; then
      let nproc=nproc+1
    fi
  done < pars.in
  
#execute drutes function in parallel
let z=0
while read l a b c d e f g h i j k m n o p q r s t u v w
  do
    if [[  $l == "p"  ]]; then
      let z=z+1
      if [[ $z -lt $nproc ]] ; then
        run_drutes $z $a $b $c $d $e $f $g $h $i $j $k $m $n $o $p $q $r $s $t $u $v $w&
      else
        run_drutes $z $a $b $c $d $e $f $g $h $i $j $k $m $n $o $p $q $r $s $t $u $v $w
      fi
    fi
  done < pars.in
 

#at the end of drutes function file obj.val is created, if all files for each process exist then we have finished 
let z=0    
while [[ $z -lt $nproc ]]
  do
  let z=0
    for i in `seq 1 $nproc`
      do
        FILE="$i/objfnc.val"
        if [ -f $FILE ];
          then
             let z=z+1
        fi
    done    
done

for i in `seq 1 $nproc` ; do
  val=`cat $i/objfnc.val`
  echo $i $val >> drutes.vals
done
