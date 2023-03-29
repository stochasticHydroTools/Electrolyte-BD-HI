positions=$1
datamain=$2
t1DL=$3 # in terms of H^2/D
t2DL=$4 # in terms of H^2/D
tauEqDL=$5

N=$(grep -Eo '^numberParticles[[:space:]].*' $datamain | awk '{print $2}')
dt=$(grep -Eo '^dt[[:space:]].*' $datamain | awk '{print $2}')
printSteps=$(grep -Eo '^printSteps[[:space:]].*' $datamain | awk '{print $2}')
Nprint=$N

t1=$(echo | awk '{print '$t1DL'*'$tauEqDL'}')
t2=$(echo | awk '{print '$t2DL'*'$tauEqDL'}')
step1=$(echo | awk '{print int('$t1'/('$dt'*'$printSteps'))}')
step2=$(echo | awk '{print int('$t2'/('$dt'*'$printSteps'))}')
echo $step1
echo $step2

awk '/#/{x=NR+'$Nprint'}(x>=NR){if (NR>'$step1'*('$N'+1) && NR<='$step2'*('$N'+1)) print $0 NR}' pos.dat > poscut$t1DL-$t2DL.dat
