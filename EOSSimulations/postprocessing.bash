charged=$1
numberSimulations=$2
if [ $charged -eq 0 ]; then
    EqDIR="Equilibrium_unchargedSurface-longrun"
elif [ $charged -eq 1 ]; then
    EqDIR="Equilibrium_chargedSurface-longrun"
fi

script=$EqDIR"/runEquilibrium.bash"
r_ion=$(cat $script | awk -F "=" '/r_ion=/{print $2}')
eta=$(cat $script | awk -F "=" '/eta=/{print $2}')
elemCharge=$(cat $script | awk -F "=" '/elemCharge=/{print $2}')
kT=$(cat $script | awk -F "=" '/kT=/{print $2}')
D=$(echo | awk '{printf "%.14g\n", '$kT'/(6*atan2(0,-1)*'$eta'*'$r_ion')}')
boxHeight=$(cat $script | awk -F "=" '/boxHeight=/{print $2}')
tau=$(echo | awk '{printf "%.14g\n", '$r_ion'^2/'$D'}')
tauEq=$(echo | awk '{printf "%.14g\n", '$boxHeight'^2/'$D'}')
source scalings.bash
tauEqDL=$(echo | awk '{printf "%.14g\n", '$tauEq'/'$time'}')
dc=$(cat $script | awk -F "=" '/dc=/{print $2}')
vacuumPermitivity=$(cat $script | awk -F "=" '/vacuumPermitivity=/{print $2}')
conc=$(cat $script | awk -F "=" '/conc=/{print $2}')
NA=$(cat $script | awk -F "=" '/NA=/{print $2}')

debyeLength=$(echo | awk '{print sqrt('$dc'*'$vacuumPermitivity'*'$kT'/(2*('$elemCharge')^2*'$conc'*'$NA'))}')

Emax=$(echo | awk '{print '$kT$'/('$elemCharge'*'$debyeLength')*'$charge'*'$length'/'$energy'}')

currentSIcoeff=$(echo | awk '{print '$charge'/('$length'^2*'$time')}')
ESIcoeff=$(echo | awk '{print '$energy'/('$charge'*'$length')}')

rm Current-vs-E-SI.dat 2> /dev/null || true

listOfEmaxCoeff="0.01 0.05 0.1 0.2 0.5 1"
for EmaxCoeff in $listOfEmaxCoeff
do

if [ $charged -eq 0 ]; then
    DIR="EOS_unchargedSurface_EmaxCoeff=$EmaxCoeff"
elif [ $charged -eq 1 ]; then
    DIR="EOS_chargedSurface_EmaxCoeff=$EmaxCoeff"
fi
cd $DIR
externalFieldx=$(echo | awk '{print '$EmaxCoeff'*'$Emax'}')
externalFieldxSI=$(echo | awk '{print '$externalFieldx'*'$ESIcoeff'}')

for i in $(seq 1 $numberSimulations)
do

bash ../tools/qvmeasure.sh pos.$i.dat data.main $i

done
cat qvsum.*.dat | datamash -W mean 1 sstdev 1 mean 2 sstdev 2 mean 3 sstdev 3 | awk '{print $1, $2/sqrt('$numberSimulations'), $3, $4/sqrt('$numberSimulations'), $5, $6/sqrt('$numberSimulations')}' > qvsumAverage.dat
cat qvsumAverage.dat | awk '{print '$currentSIcoeff'*$1,'$currentSIcoeff'*$2,'$currentSIcoeff'*$3,'$currentSIcoeff'*$4,'$currentSIcoeff'*$5,'$currentSIcoeff'*$6}' > currentSI.dat

cat currentSI.dat | awk '{print '$externalFieldxSI', $1, $2}' >> ../Current-vs-E-SI.dat

cd ..
done

if [ $charged -eq 0 ]; then
    rm Current-vs-E-SI-unchargedSurface.dat 2> /dev/null || true
    mv Current-vs-E-SI.dat ./Current-vs-E-SI-unchargedSurface.dat
elif [ $charged -eq 1 ]; then
    rm Current-vs-E-SI-chargedSurface.dat 2> /dev/null || true
    mv Current-vs-E-SI.dat ./Current-vs-E-SI-chargedSurface.dat
fi

# bash tools/computeIonDensityProfilexyz.sh $DIR/displacement.$i.dat $DIR/data.main $dir | awk '{print $1, $2 > "displacementHistogramplus'$dir'.dat"; print $1, $3 > "displacementHistogramminus'$dir'.dat"}'
