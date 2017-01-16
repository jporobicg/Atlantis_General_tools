#!/bin/sh
# add to the repository the following ines
# first the Backup!
print "A backup of your source list is made \n"
sudo cp /etc/apt/sources.list /etc/apt/sources.list.backup

sudo add-apt-repository "deb http://qgis.org/debian trusty main"
sudo add-apt-repository "deb-src http://qgis.org/debian trusty main"

# add public keys
gpg --keyserver keyserver.ubuntu.com --recv DD45F6C3
gpg --export --armor DD45F6C3 | sudo apt-key add -

sudo apt-get update
sudo apt-get install qgis
sudo apt-get python-qgis
sudo apt-get qgis-plugin-grass
