positions=$1
datamain=$2
dir=$3 # 1:x, 2:y, 3:z
nbins=200
nsteps=1 #$(cat $positions | grep -c "#");
N=$(grep -Eo '^numberParticles[[:space:]].*' $datamain | awk '{print $2}')

if [ $dir -eq 1 ]; then
    halfsize=$(grep -Eo '^Lxy[[:space:]].*' $datamain | awk '{print 0.5*$2}')
elif [ $dir -eq 2 ]; then
    halfsize=$(grep -Eo '^Lxy[[:space:]].*' $datamain | awk '{print 0.5*$2}')
elif [ $dir -eq 3 ]; then
    halfsize=$(grep -Eo '^H[[:space:]].*' $datamain | awk '{print 0.5*$2}')
else
    echo "ERROR x,y,z?"
    exit
fi

function histogram {
nbins=$1
upper=$2
lower=$3

awk 'BEGIN{min='$lower'; max='$upper'; nbins='$nbins';}{b=int(($1-min)/(max-min)*nbins)}$3>0{hplus[b] += $2;}$3<0{hminus[b]+=$2}
     END{for(i=0;i<nbins;i++){
	   z=(i+0.5)/nbins*(max-min)+min;
	   print z/'$halfsize', hplus[i]*1.0, hminus[i]*1.0;
	 }
	}'
}


if [ $dir -eq 1 ]; then
    cat $positions |
	grep -v "#" |
	awk '{printf "%.14g %.14g %.14g\n", $1,('$nbins'/('$N'*'$nsteps')),$4}' |
	histogram $nbins $halfsize -$halfsize
elif [ $dir -eq 2 ]; then
    cat $positions |
	grep -v "#" |
	awk '{printf "%.14g %.14g %.14g\n", $2,('$nbins'/('$N'*'$nsteps')),$4}' |
	histogram $nbins $halfsize -$halfsize
elif [ $dir -eq 3 ]; then
    cat $positions |
	grep -v "#" |
	awk '{printf "%.14g %.14g %.14g\n", $3,('$nbins'/('$N'*'$nsteps')),$4}' |
	histogram $nbins $halfsize -$halfsize
fi
