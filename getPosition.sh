#!/bin/bash
# Ask the user for their name
echo Input file name ie "filename.ubx"
read filename
echo RINEX now processing: $filename

RTKpath='Documents/GitHub/GNSSProcess/RTKLIB-demo5/app/consapp/convbin/gcc'

cd $RTKpath
base="/Users/derekpickell/Documents/GitHub/GNSSProcess/tempData/"
filepath=$base$filename
echo $filepath

./convbin -r ubx $filepath
# change settings to output rinex2

echo Do you want to process position? y/n
read response

if [$response == "y"]; then
#	cd ../..
#	cd rnx2rtkp/gcc
#	./rnx2rtkp 


