#!/bin/bash

yes="y"
no="n"
ubx="UBX"
o=".obs"
n=".nav"

echo Enter filename, without extension
read filename

echo Do you want to process RAW to RINEX? y/n
read response

BASEDIR=$(dirname "$0")
CONVBINpath=$BASEDIR"/RTKLIB-demo5/app/consapp/convbin/gcc"
RNX2RTKPpath=$BASEDIR"/RTKLIB-demo5/app/consapp/rnx2rtkp/gcc"

if [ $response == $yes ]; then

    ubxFILE="$filename.$ubx"
    echo RINEX now processing: $ubxFILE
    cd $CONVBINpath
    
    RAWPATH=$BASEDIR"/tempData/"$ubxFILE
    echo $RAWPATH

    ./convbin -r ubx $RAWPATH -v 2.11 -os
fi


echo Do you want to process position? y/n
read response

if [ $response == $yes ]; then
    OBFILE="$filename$o"
    NAVFILE="$filename$n"
    
    OB=$BASEDIR"/tempData/"$OBFILE
    NAV=$BASEDIR"/tempData/"$NAVFILE
    echo $OBFILE
    echo $NAVFILE
    cd ~
    cd $RNX2RTKPpath
    ./rnx2rtkp -p 0 -m 15 -t -e -o $filename".pos" $OB $NAV
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
