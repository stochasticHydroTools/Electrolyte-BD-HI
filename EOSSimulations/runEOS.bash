charged=$1
numberSimulations=$2
if [ $charged -eq 0 ]; then
    EqDIR="Equilibrium_unchargedSurface-longrun"
elif [ $charged -eq 1 ]; then
    EqDIR="Equilibrium_chargedSurface-longrun"
fi

wetFraction=0.5
script="runEquilibrium.bash"
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
dtCoeff=0.05
dt=$(echo | awk '{print '$dtCoeff'*'$tau'/'$time'}')

dc=$(cat $script | awk -F "=" '/dc=/{print $2}')
vacuumPermitivity=$(cat $script | awk -F "=" '/vacuumPermitivity=/{print $2}')
conc=$(cat $script | awk -F "=" '/conc=/{print $2}')
NA=$(cat $script | awk -F "=" '/NA=/{print $2}')

debyeLength=$(echo | awk '{print sqrt('$dc'*'$vacuumPermitivity'*'$kT'/(2*('$elemCharge')^2*'$conc'*'$NA'))}')

# disable fold and add an external field
Emax=$(echo | awk '{print '$kT$'/('$elemCharge'*'$debyeLength')*'$charge'*'$length'/'$energy'}')
fold=0

Tot_Sim_Time=$(echo | awk '{print 1*'$tauEqDL'}')
Snap_Time=$(echo $H | awk '{print 0.01*'$tauEqDL'}')
Relax_Time=$(echo | awk '{print 20*'$tau'/'$time'}')

numberSteps=$(echo $Tot_Sim_Time $dt | awk '{print int($1/$2)}')
printSteps=$(echo $Snap_Time $dt | awk '{print int($1/$2)}')
relaxSteps=$(echo $Relax_Time $dt | awk '{print int($1/$2)}')

cd $EqDIR
bash ../tools/equilibrium-configuration-extract.bash pos.dat $numberSimulations $tauEqDL
cd ..

listOfEmaxCoeff="0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1"
for EmaxCoeff in $listOfEmaxCoeff
do
    
externalFieldx=$(echo | awk '{print '$EmaxCoeff'*'$Emax'}')

cp $EqDIR/data.main ./
# make chnages to the data.main
awk '/fold/{gsub(/1/, '$fold')};{print}' data.main > data.main-update
rm data.main && mv data.main-update ./data.main
awk '/externalField/{$0="externalField '$externalFieldx' 0 0"};{print}' data.main >> data.main-update
rm data.main && mv data.main-update ./data.main
awk '/numberSteps/{$0="numberSteps '$numberSteps'"};{print}' data.main >> data.main-update
rm data.main && mv data.main-update ./data.main
awk '/printSteps/{$0="printSteps '$printSteps'"};{print}' data.main >> data.main-update
rm data.main && mv data.main-update ./data.main
awk '/relaxSteps/{$0="relaxSteps '$relaxSteps'"};{print}' data.main >> data.main-update
rm data.main && mv data.main-update ./data.main
awk '/dt/{$0="dt '$dt'"};{print}' data.main >> data.main-update
rm data.main && mv data.main-update ./data.main
awk '/wetFraction/{$0="wetFraction '$wetFraction'"};{print}' data.main >> data.main-update
rm data.main && mv data.main-update ./data.main

if [ $charged -eq 0 ]; then
    DIR="EOS_unchargedSurface_EmaxCoeff=$EmaxCoeff"
elif [ $charged -eq 1 ]; then
    DIR="EOS_chargedSurface_EmaxCoeff=$EmaxCoeff"
fi
mkdir -p $DIR

if ! test -f ../build/slab_$HOSTNAME; then echo "ERROR: This script needs slab_$HOSTNAME to be compiled and available at ../buil/ " >/dev/stderr; exit 1; fi

for i in $(seq 1 $numberSimulations)
do
    {
    	cat $EqDIR/eqpos$i.dat | ../DPPoissonTests/poisson_$HOSTNAME data.main > $DIR/pos.$i.dat 2> log;
    } || {
    	echo "ERROR: UAMMD Failed, here is the log" > /dev/stderr
    	cat log > /dev/stderr
    	exit 1
    }
done

mv data.main $DIR
mv log $DIR
done
