#!/bin/sh

#printf "\n\nDownloading and making databases available for CoMW\n\n"
#Downloading the databases and unzipping them
#wget -A zip -r -l 1 -nd -P ./databases/ https://mzacomw.au.dk
#for i in ./databases/*.zip
#do 
#	unzip -o $i -d ./databases/
#	rm -rf $i
#done


printf "\n\nDownloading and installing sword for CoMW\n\n"
git clone --recursive https://github.com/rvaser/sword.git sword
cd sword/
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make

printf "\n\nCoMW installed, ready to use."