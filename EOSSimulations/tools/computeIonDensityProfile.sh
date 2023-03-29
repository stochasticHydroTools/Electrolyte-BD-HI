positions=$1
datamain=$2
nbins=200
nsteps=$(cat $positions | grep -c "#");
H=$(grep -Eo '^H[[:space:]].*' $datamain | awk '{print $2}')
Hhalf=$(echo $H | awk '{print $1*0.5}')
N=$(grep -Eo '^numberParticles[[:space:]].*' $datamain | awk '{print $2}')

function histogram {
nbins=$1
upper=$2
lower=$3

awk 'BEGIN{min='$lower'; max='$upper'; nbins='$nbins';}{b=int(($1-min)/(max-min)*nbins)}$3>0{hplus[b] += $2;}$3<0{hminus[b]+=$2}
     END{for(i=0;i<nbins;i++){
	   z=(i+0.5)/nbins*(max-min)+min;
	   print z/'$Hhalf', hplus[i]*1.0, hminus[i]*1.0;
	 }
	}'
}

cat $positions |
    grep -v "#" |
    awk '{printf "%.14g %.14g %.14g\n", $3,('$nbins'/('$N'*'$nsteps')),$4}' |
    histogram $nbins $Hhalf -$Hhalf
