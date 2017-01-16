#!/bin/sh
# know if something is installed the basic staff for Atlantis()
dpkg -l | grep build-essentials
dpkg -l | grep autoconf
dpkg -l | grep subversion


# if I do not have somo of this library
sudo apt-get install build-essential   # ok
sudo apt-get install autoconf          # ok
sudo apt-get install subversion        # ok


#install proj
/home/demiurgo/proj-4.8.0
wget http://download.osgeo.org/proj/proj-4.9.1.tar.gz
sudo ./configure
sudo make
sudo make install


#  install libxml2-dev
sudo apt-get install libxml2-dev

cd  /home/demiurgo/Documents/2015/Atlantis_Model

# checkout Atlantis
svn co  https://svnserv.csiro.au/svn/ext/atlantis/Atlantis/trunk/

# Checkout SETas
svn co  https://svnserv.csiro.au/svn/ext/atlantis/runFiles/trunk/SETas_model_New

# update
svn update [filename - foldername]

#  install netcdf library
sudo apt-get install libnetcdf-dev
# install gawkl
sudo apt-get install gawk

# build atlantis
cd  /home/demiurgo/Documents/2015/Atlantis_Model/trunk/atlantis



aclocal
autoconf
automake -a
sudo chmod +x configure
./configure

make
sudo make install
