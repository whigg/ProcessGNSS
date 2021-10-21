#!/bin/bash

yes="y"
no="n"
dat="dat"
o=".obs"
n=".gps"

echo Enter filename, without .dat extension
read filename

echo Do you want to process RAW to RINEX? y/n
read response

BASEDIR=$(dirname "$0")

if [ $response == $yes ]; then

    datFILE="$filename.$dat"
    echo RINEX now processing: $datFILE notice: currently configured only for GPS nav messages
    cd $BASEDIR

    ./teqc -O.obs L1+l2+ca+p2+s1+s2 +nav $filename".gps" $filename".dat" > $filename".obs" 
fi


echo Do you want to process position? y/n
read response

if [ $response == $yes ]; then
    OBFILE="$filename$o"
    NAVFILE="$filename$n"
    
    OB=$BASEDIR"/"$OBFILE
    NAV=$BASEDIR"/"$NAVFILE
    echo $OBFILE
    echo $NAVFILE
    cd $BASEDIR
    cd ..
    cd RTKLIB-demo5/app/consapp/rnx2rtkp/gcc
    ./rnx2rtkp -p 0 -m 15 -e -u -o $filename".pos" $OB $NAV
    mv $filename".pos" $BASEDIR
    mv $filename"_events.pos" $BASEDIR
    
 fi

echo Do you want to process statistics? y/n
read response

if [ $response == $yes ]; then
    echo Python will process the .pos file moved to /tempData
    cd $BASEDIR
    cd ..
    source env/bin/activate
    python3 dataAnalysis.py
 fi
