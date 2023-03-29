# This script calculates the \left(\sum_{i=1}^{N}\Delta x_i q_i/\tau\right)/V
positions=$1
datamain=$2
index=$3
N=$(grep -Eo '^numberParticles[[:space:]].*' $datamain | awk '{print $2}')
H=$(grep -Eo '^H[[:space:]].*' $datamain | awk '{print $2}')
Lxy=$(grep -Eo '^Lxy[[:space:]].*' $datamain | awk '{print $2}')
V=$(echo $H $Lxy | awk '{print $1*$2^2}')

printedSTEPS=$(cat $positions | grep -c "#")
numberSteps=$(grep -Eo '^numberSteps[[:space:]].*' $datamain | awk '{print $2}')
dt=$(grep -Eo '^dt[[:space:]].*' $datamain | awk '{print $2}')
t=$(echo | awk '{print '$numberSteps'*'$dt'}')

awk '{if (NR>1 && NR<=1+'$N') print $0}' $positions > xini.dat
awk '{if (NR>1+('$printedSTEPS'-1)*('$N'+1) && NR<=1+'$printedSTEPS'*('$N'+1)) print $0}' $positions > xend.dat

# awk '(NR==FNR){xold[NR]=$1;yold[NR]=$2;zold[NR]=$3;next}{print ($1-xold[FNR]),($2-yold[FNR]),($3-zold[FNR]),$4}' xini.dat xend.dat > displacement.$index.dat

awk '(NR==FNR){xold[NR]=$1;yold[NR]=$2;zold[NR]=$3;next}{print $4*($1-xold[FNR])/('$t'*'$V'),$4*($2-yold[FNR])/('$t'*'$V'),$4*($3-zold[FNR])/('$t'*'$V')}' xini.dat xend.dat | awk '{sumx += $1;sumy += $2;sumz += $3} END {print sumx, sumy, sumz}' > qvsum.$index.dat

rm xini.dat
rm xend.dat
