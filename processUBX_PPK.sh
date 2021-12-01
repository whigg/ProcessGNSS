#!/bin/bash

yes="y"
no="n"
ubx="UBX"
o=".obs"
n=".nav"

echo Enter filename, without extension
read filename

BASEDIR=$(dirname "$0")
CONVBINpath=$BASEDIR"/RTKLIB-demo5/app/consapp/convbin/gcc"
RNX2RTKPpath=$BASEDIR"/RTKLIB-demo5/app/consapp/rnx2rtkp/gcc"

echo Do you want to process position? y/n
read response

if [ $response == $yes ]; then
    echo Enter base station filename [.obs], without extension and with wildcards iff necessary
    read correction
    CORRECTIONFILE="$correction$o"
    OBFILE="$filename$o"
    NAVFILE="$filename$n"
    OB=$BASEDIR"/tempData/"$OBFILE
    NAV=$BASEDIR"/tempData/"$NAVFILE
    COR=$BASEDIR"/tempData/"$CORRECTIONFILE
    echo $OBFILE
    echo $NAVFILE
    echo $baseECEF
    cd ~
    cd $RNX2RTKPpath
    ./rnx2rtkp -k ppkStatic.conf -u -o $filename".pos" $OB $COR $NAV ##############
    mv $filename".pos" $BASEDIR"/tempData/"
    mv $filename"_events.pos" $BASEDIR"/tempData/"
    
fi

echo Do you want to process statistics? y/n
read response

if [ $response == $yes ]; then
    echo Python will automatically process the .pos file \in /tempData
    cd $BASEDIR
    source env/bin/activate
    python3 dataAnalysis.py 
fi