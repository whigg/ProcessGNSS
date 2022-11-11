#!/bin/zsh
#converts whatever file is in directory to proper name/format etc. 

#CORE DIRECTORIES - CHANGE ACCORDINGLY
GNSSDIR="/Users/derek/Documents/GNSSIR/"
BASEDIR=$(pwd)
export EXE=$GNSSDIR"/exe"
export REFL_CODE=$GNSSDIR"/io"
export ORBITS=$GNSSDIR"/io"

echo Ensure file is in current directory. Enter file name:
read filename

#GET STATION DETAILS AND TIMES, MOVE FILE & RENAME
station=${filename:0:4}
year=${filename:5:4}
yearShort=${filename:7:2}
month=${filename:9:2}
day=${filename:11:2}
echo station: $station year: $year month: $month day: $day

RINEXDIR=$GNSSDIR"io/rinex/"$station"/"$year
echo moving file to GNSSIR folder directory: $RINEXDIR
mkdir -p $RINEXDIR
cp $filename $RINEXDIR

cd $GNSSDIR"/gnssrefl/gnssrefl"
source env/bin/activate
doy=$(ymd $year $month $day)
rename=$station$doy"0."$yearShort"o"
echo doy: $doy, renaming file to: $rename
mv $RINEXDIR"/"$filename $RINEXDIR"/"$rename

#CREATE SNR FILE AND VIEW GNSS-IR
# echo creating snr file...
# rinex2snr $station $year $doy -nolook True -orb gnss -overwrite True

# echo visualizing reflections...
# quickLook $station $year $doy -fr 1

line=11
xyz=$(sed -n "$line p" $RINEXDIR"/"$rename)
xyz_parsed=${xyz:1:39}
echo xyz_parsed
echo getting station approximate location [$xyz_parsed] to make JSON
# make_json_input bamb 1501145.6369 -1191418.2057  6066411.5157 -xyz True -allfreq True
# gnssir $station $year $doy # -plt True 