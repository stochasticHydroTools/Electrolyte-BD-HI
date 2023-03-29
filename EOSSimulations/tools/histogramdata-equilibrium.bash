# this script cuts the pos.dat file from (i-\tau)H^2/D to iH^2/D
# execute while in the Equilibrium long run directory: bash ../tools/histogramdata-equlibrium.bash tau

period=$1
numberSimulations=$2
tauEqDL=$3
for i in $(seq 1 $numberSimulations)
do
    t1DL=$(echo | awk '{print '$i'-'$period'}')
    t2DL=$(echo | awk '{print '$i'}')
    bash ../tools/poscut.bash pos.dat data.main $t1DL $t2DL $tauEqDL
    bash ../tools/computeIonDensityProfile.sh poscut$t1DL-$t2DL.dat data.main | awk '{print $1, $2 > "densityplus'$i'.dat"; print $1, $3 > "densityminus'$i'.dat"}'
    rm poscut$t1DL-$t2DL.dat
done
