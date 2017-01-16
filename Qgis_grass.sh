#!/bin/sh
# add to the repository the following ines
# first the Backup!
echo "A backup of your source list is made \n"
sudo cp /etc/apt/sources.list /etc/apt/sources.list.backup
echo "\n Backup Done \n"
# installing GRASS
echo "\n STARTING:::Installation of GRASS version 7\n"
sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
sudo add-apt-repository ppa:grass/grass-stable
sudo apt-get update
sudo apt-get install grass7

echo "\n DONE:::Installation of GRASS version 7\n"

echo "\n STARTING:::Installation of QGIS and plugins\n"
sudo add-apt-repository "deb http://qgis.org/debian trusty main"
sudo add-apt-repository "deb-src http://qgis.org/debian trusty main"

# add public keys
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-key 3FF5FFCAD71472C4

gpg --export --armor 3FF5FFCAD71472C4 | sudo apt-key add -


sudo apt-get update
sudo apt-get install qgis
sudo apt-get install python-qgis
sudo apt-get install qgis-plugin-grass

echo "\n DONE:::Installation QGIS and plugins\n"
