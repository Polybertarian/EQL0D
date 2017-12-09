#!/bin/bash/
# Hard-coded constants
M_NEUTRON=1.0086649670000E+00;
N_AVOGADRO=6.0220434469282E-01;
MEV=1.6021765314E-13;

# Prepare MATLAB file
LIBRARY_NAME=$(basename $(dirname $1))
echo "clearvars;" > "$LIBRARY_NAME".m

# Grab nuclide data from .dec file
grep "8457    1\|8457    2\|8457    3" $1 | awk 'BEGIN{FIELDWIDTHS = "11 11 11 11 11 11 14"}{print $1" "$2" "$3" "$4" "$5" "$6" "$7}' |  sed 's/E+/+/g' |  sed 's/E-/-/g' | sed 's/+/E+/g' | sed 's/-/E-/g' | column -t > data.tmp

# Grab other data from .xsdata file
grep "300" $2 | awk '{print $4$5" "$6}' | column -t | uniq > data2.tmp

# Grab list of nuclides ZAI
echo "ZAI1=[" >> "$LIBRARY_NAME".m
grep "8457  1" data.tmp | awk '{print $1"*10+"$4}' >> "$LIBRARY_NAME".m
echo "];" >> "$LIBRARY_NAME".m

#Grab atomic weight ratio
echo "AWR=[" >> "$LIBRARY_NAME".m
grep "8457  1" data.tmp | awk '{print $2}' >> "$LIBRARY_NAME".m
echo "];" >> "$LIBRARY_NAME".m
# Convert AWR to Atomic Weight 
echo "AW=AWR*$M_NEUTRON;" >> $LIBRARY_NAME.m

# Grab half-life
echo "T12=[" >> "$LIBRARY_NAME".m
grep "8457  2" data.tmp | awk '{print $1}' >> "$LIBRARY_NAME".m
echo "];" >> "$LIBRARY_NAME".m

# Grab decay energy
echo "EV=[" >> "$LIBRARY_NAME".m
grep "8457  3" data.tmp | awk '{print "("$1"+"$3"+"$5")"}' >> "$LIBRARY_NAME".m
echo "];" >> "$LIBRARY_NAME".m
# Convert to Joules
echo "Q=EV*$MEV;" >> $LIBRARY_NAME.m

# Grab other list of nuclides ZAI
echo "ZAI2=[" >> "$LIBRARY_NAME".m
awk '{print $1}' data2.tmp >> "$LIBRARY_NAME".m
echo "];" >> "$LIBRARY_NAME".m

# Grab other list of nuclides ZAI
echo "AW2=[" >> "$LIBRARY_NAME".m
awk '{print $2}' data2.tmp >> "$LIBRARY_NAME".m
echo "];" >> "$LIBRARY_NAME".m

echo "AW(ismember(ZAI1,ZAI2))=AW2(ismember(ZAI2,ZAI1));" >> $LIBRARY_NAME.m
echo "AW2(ismember(ZAI2,ZAI1))=[];" >> $LIBRARY_NAME.m
echo "ZAI2(ismember(ZAI2,ZAI1))=[];" >> $LIBRARY_NAME.m
echo "ZAI1=[ZAI1;ZAI2];[I,J]=sort(ZAI1);ZAI1=ZAI1(J);" >> $LIBRARY_NAME.m
echo "AW=[AW;AW2];AW=AW(J);" >> $LIBRARY_NAME.m
echo "Q=[Q;zeros(size(ZAI2))];Q=Q(J);" >> $LIBRARY_NAME.m
echo "T12=[T12;zeros(size(ZAI2))];T12=T12(J);" >> $LIBRARY_NAME.m

echo "DAT.T12=T12;" >> $LIBRARY_NAME.m
echo "DAT.Q=Q;" >> $LIBRARY_NAME.m
echo "DAT.ZAI0=ZAI1;" >> $LIBRARY_NAME.m

# Convert to g/atom 
echo "DAT.AMASS=AW/$N_AVOGADRO/1.0E+24;" >> $LIBRARY_NAME.m

# Remove single neutron from list 
echo "DAT.AMASS(1)=[];DAT.Q(1)=[];DAT.T12(1)=[];DAT.ZAI0(1)=[];" >> $LIBRARY_NAME.m
echo "clearvars -except DAT" >> $LIBRARY_NAME.m
echo "save $LIBRARY_NAME.mat" >> "$LIBRARY_NAME".m

#matlab -nodesktop -nosplash -nojvm -r "run('$LIBRARY_NAME.m'),exit"

#rm -rf data.tmp data2.tmp "$LIBRARY_NAME".m

exit 0


