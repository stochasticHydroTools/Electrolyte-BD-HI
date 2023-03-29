N_wall=$1
q=$2 
datamain="data.main"
H=$(grep -Eo "^H[[:space:]].*" $datamain | awk '{print $2}')
Lxy=$(grep -Eo "^Lxy[[:space:]].*" $datamain | awk '{print $2}')
numberParticles=$(grep -Eo '^numberParticles[[:space:]].*' $datamain | awk '{print $2}')
N_rest=$(echo $numberParticles $N_wall | awk '{print $1-$2}')

seed=$(head -100 /dev/urandom | cksum | awk '{print $1}')
seq $N_wall |
    awk 'BEGIN{srand('$seed')}
        {acc=0;
         while(acc==0){
                H='$H';
                x=(rand()-0.5)*'$Lxy';
                y=(rand()-0.5)*'$Lxy';
                z=(rand()-0.5)*(H-2.0);
                p=1.0;
                Z=rand();
                if(Z<p){
                    acc=1;
                    print x,y,z,'$q';
                }
          }
        }'

seed=$(head -100 /dev/urandom | cksum | awk '{print $1}')
seq $N_rest |
    awk 'BEGIN{srand('$seed')}
        {acc=0;
         while(acc==0){
                H='$H';
                x=(rand()-0.5)*'$Lxy';
                y=(rand()-0.5)*'$Lxy';
                z=(rand()-0.5)*(H-2.0);
                p=1.0;
                Z=rand();
                if(Z<p){
                    acc=1;
                    print x,y,z,(NR%2-0.5)*2;
                }
          }
        }'

