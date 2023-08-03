charged=$1
numberSimulations=$2 # simulation will be run for $numberSimulation*tau_Eq, where tau_Eq is the approximate time required for the equilibrium
wetFraction=$3

# Dimensional parameters
boxHeight=10e-9
boxWidth=20e-9
r_ion=2.5e-10
eta=8.9e-4
dc=78.3
dcTop=2.55 # polystyrene
dcBottom=2.55
#mol/m^3 (mM)
conc=1000
if [ $charged -eq 0 ]; then
    zetaPot=0
elif [ $charged -eq 1 ]; then
    zetaPot=-100e-3
fi
kT=4.11e-21
elemCharge=1.6e-19
D=$(echo | awk '{printf "%.14g\n", '$kT'/(6*atan2(0,-1)*'$eta'*'$r_ion')}')
tau=$(echo | awk '{printf "%.14g\n", '$r_ion'^2/'$D'}')
NA=6.022e23
vacuumPermitivity=8.8541878e-12
debyeLength=$(echo | awk '{print sqrt('$dc'*'$vacuumPermitivity'*'$kT'/(2*('$elemCharge')^2*'$conc'*'$NA'))}')

echo $debyeLength $boxHeight | awk '{print "H/lambda = " $2/$1}'
surfaceChargeDimensional=$(echo $dc $vacuumPermitivity $debyeLength $zetaPot | awk '{print $1*$2*$4/$3}')

# scalings
source scalings.bash
###### dimensionless parameters for data.main ######
####################################################
hydrodynamicRadius=$(echo | awk '{print '$r_ion'/'$length'}')
hxy_stokes=-1 #grid size for DPStokes will be selected automatically
temperature=$(echo | awk '{print '$kT'/'$energy'}')
viscosity=$(echo | awk '{print '$eta'*'$length'^3/('$energy'*'$time')}')
dtCoeff=0.2
dt=$(echo | awk '{print '$dtCoeff'*'$tau'/'$time'}')
gwCoeff=0.25
gw=$(echo | awk '{print '$gwCoeff'*'$r_ion'/'$length'}')
H=$(echo | awk '{print '$boxHeight'/'$length'}')
Lxy=$(echo | awk '{print '$boxWidth'/'$length'}')
permitivity=$(echo | awk '{printf "%.14g\n", '$dc'*'$vacuumPermitivity'*'$energy'*'$length'/('$charge'^2)}')
permitivityTop=$(echo | awk '{printf "%.14g\n", '$dcTop'*'$vacuumPermitivity'*'$energy'*'$length'/('$charge'^2)}')
permitivityBottom=$(echo | awk '{printf "%.14g\n", '$dcBottom'*'$vacuumPermitivity'*'$energy'*'$length'/('$charge'^2)}')

#Repulsive potential parameters
U0=$(echo | awk '{print 0.6*'$kT'/'$energy'}')
r_m=$(echo | awk '{print 1.7*'$r_ion'/'$length'}')
p=2
sigma=$(echo | awk '{print 2*'$r_ion'/'$length'}')

#Wall params
wall_U0=$(echo | awk '{print 0.25*'$kT'/'$energy'}')
wall_r_m=$(echo | awk '{print 1.7*'$r_ion'/'$length'}')
wall_p=3
wall_sigma=$(echo | awk '{print 2*'$r_ion'/'$length'}')

#finding the equilibrium conditions in the absence of electric field  
fold=1
externalFieldx=0
externalFieldy=0
externalFieldz=0

####################################################
BrownianUpdateRule=Leimkuhler #EulerMaruyama

echo "Dimensional surface charge (C/m^2) = " $surfaceChargeDimensional
N_wall=$(echo | awk '{print int(sqrt(('$surfaceChargeDimensional')^2)*('$boxWidth')^2/'$elemCharge')}')
cond=$(echo $zetaPot | awk '{if ($1 > 0.0) {print "True"} else {print "False"} }')
if [ $cond = True ]; then
    q=-1
else
    q=1
fi

numberParticles=$(echo | awk '{print 2*int('$boxHeight'*'$boxWidth'^2*'$conc'*'$NA')+'$N_wall'}')
echo "Number of ions to screen the surface charge = " $N_wall
echo "Total number of ions = " $numberParticles

function createDataMain {

cat<<EOF > $1

numberParticles $numberParticles
gw $gw

Lxy $Lxy
H $H

permitivity $permitivity
permitivityTop $permitivityTop
permitivityBottom $permitivityBottom

externalField $externalFieldx $externalFieldy $externalFieldz

fold $fold

temperature $temperature
viscosity $viscosity
hydrodynamicRadius $hydrodynamicRadius
wetFraction $wetFraction
hxy_stokes $hxy_stokes

numberSteps $numberSteps
printSteps $printSteps
relaxSteps $relaxSteps
dt $dt

outfile /dev/stdout
readFile /dev/stdin

#Repulsive potential parameters
U0 $U0
r_m $r_m
p $p
sigma $sigma

wall_U0 $wall_U0
wall_r_m $wall_r_m
wall_p $wall_p
wall_sigma $wall_sigma

#DPPoisson splitting, number of grid cells. You may switch this parameter to "split" instead. 
Nxy 72

imageDistanceMultiplier 1
BrownianUpdateRule $BrownianUpdateRule

velocityfile fluidVelocity

EOF
}

# initial relaxation parameters
dt_s=$dt
dt=$(echo | awk '{print 0.001*'$tau'/'$time'}')
numberSteps=1
printSteps=1
relaxSteps=500

createDataMain data.main.relax

# actual simulation parameters
dt=$dt_s
tauEq=$(echo | awk '{printf "%.14g\n", '$boxHeight'^2/'$D'}') # time required for the equilibrium to achieve
tauEqDL=$(echo | awk '{printf "%.14g\n", '$tauEq'/'$time'}')
Tot_Sim_Time=$(echo | awk '{print '$numberSimulations'*'$tauEqDL'}') # total simulation time
Snap_Time=$(echo $H | awk '{print 0.01*'$tauEqDL'}') # how often data will be printed
Relax_Time=$(echo | awk '{print 20*'$tau'/'$time'}') # further initial relaxation

numberSteps=$(echo $Tot_Sim_Time $dt | awk '{print int($1/$2)}')
printSteps=$(echo $Snap_Time $dt | awk '{print int($1/$2)}')
relaxSteps=$(echo $Relax_Time $dt | awk '{print int($1/$2)}')

createDataMain data.main

if ! test -f ../build/slab; then echo "ERROR: This script needs slab to be compiled and available at ../build/ " >/dev/stderr; exit 1; fi

if [ $charged -eq 0 ]; then
    DIR="Equilibrium_unchargedSurface-longrun"
elif [ $charged -eq 1 ]; then
    DIR="Equilibrium_chargedSurface-longrun"
fi

mkdir -p $DIR

bash tools/init.sh $N_wall $q > $DIR/initpos.dat

{
    cat $DIR/initpos.dat | ../build/slab data.main.relax > $DIR/initpos-relaxed.dat 2> log.relax;
    echo "Done with relaxation"
    mv fluidVelocity fluidVelocity.relax;  
    cat $DIR/initpos-relaxed.dat | grep -v '#' | ../build/slab data.main > $DIR/pos.dat 2> log;
} || {
    echo "ERROR: UAMMD Failed, here is the log" > /dev/stderr
    cat log > /dev/stderr
    exit 1
}

mv log* $DIR/
mv data.main* $DIR/
mv fluidVelocity* $DIR/
cp runEquilibrium.bash $DIR/
cp tools/init.sh $DIR/

period=0.25 # the extend of time (fraction of tauEq) that will be used to cmpute the histrograms
cd $DIR
bash ../tools/histogramdata-equilibrium.bash $period $numberSimulations $tauEqDL
cd ..
