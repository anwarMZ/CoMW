
#!/bin/sh

printf "\n\nCloning and making Sword\n\n"
#Cloning the sword into sword directory
git clone --recursive https://github.com/rvaser/sword.git sword
cd sword/
make 

printf "\n\nDownloading and making databases available for CoMW\n\n"
#Downloading the databases and unzipping them
wget -A zip -r -l 1 -nd -P ./databases/ https://mzacomw.au.dk
for i in ./databases/*.zip
do 
	unzip -o $i -d ./databases/
	rm -rf $i
done
printf "\n\nCoMW installed, ready to use."

