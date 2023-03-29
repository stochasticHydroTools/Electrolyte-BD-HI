positions=$1
numberSimulations=$2
tauEqDL=$3
datamain="data.main"
N=$(grep -Eo '^numberParticles[[:space:]].*' $datamain | awk '{print $2}')
dt=$(grep -Eo '^dt[[:space:]].*' $datamain | awk '{print $2}')
printSteps=$(grep -Eo '^printSteps[[:space:]].*' $datamain | awk '{print $2}')
Nprint=$N

STEPS=$(cat $positions | grep -c "#")

for i in $(seq 1 $numberSimulations)
do
    step=$(echo 1 | awk '{print int('$i'*'$tauEqDL'/('$dt'*'$printSteps'))}')
    awk '/#/{x=NR+'$Nprint';next}(x>=NR){if (NR>('$step'-1)*('$N'+1)+1 && NR<='$step'*('$N'+1)) print $0}' pos.dat > eqpos$i.dat
done
