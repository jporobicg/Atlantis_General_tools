#!/bin/sh
# Instalar netcdf -Javier Porobic
#------------------
# install zlib y HDF5
# Variables
home=$HOME
ini=pwd

# Create instalation  folder
sudo mkdir /usr/programs/
sudo mkdir /usr/programs/netcdf

# Download programs from this ftps
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/hdf5-1.8.13.tar.gz
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/zlib-1.2.8.tar.gz
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.3.tar.gz
wget ftp://cirrus.ucsd.edu/pub/ncview/ncview-2.1.4.tar.gz
wget ftp://ftp.unidata.ucar.edu/pub/udunits/udunits-2.2.19.tar.gz


# Copy the files to the folder
sudo mv  *.tar.gz /usr/programs/netcdf/

# Unzip the compress files
cd /usr/programs/netcdf
#sudo su
# you need to be superuser ,  loop to uncompress the files
for file in `ls *.gz`
do
    sudo tar -xzf $file
done
# remove the compressed files
sudo rm  *.tar.gz

#  Zlbi
cd zlib-1.2.8
sudo ./configure --prefix=$home/local
sudo make check
sudo make install
cd ..

# HDF5
cd hdf5-1.8.13
sudo ./configure --with-zlib=$home/local --prefix=$home/local
sudo  make
sudo make check
sudo make install
cd ..

# NETCDF 4  -  You need to put the variables to indicate the position of hdf5 and
# zlib
cd netcdf-4.3.3
sudo CPPFLAGS=-I$home/local/include LDFLAGS=-L$home/local/lib ./configure --prefix=$home/local
sudo make check
sudo make install
export PATH=/usr/programs/netcdf/netcdf-4.3.3/ncdump/:$PATH

#cd $ini

# # Instalacion NCview
# #-------------------

#  Instalation Of Udunits
sudo apt-get install udunits-bin


# # antes de instalar NCVIEW debes instalar netcdf+4 y HDF-5

cd ncview-2.1.4
sudo ./configure --with-nc-config=$home/local/bin/nc-config
sudo make
sudo make install

# by default Ubuntu doesnt have some lybraries, you can check if you want, but
# is more easy only install it
sudo apt-get install xorg-dev
sudo apt-get install apt-file
sudo apt-file update   #update the cache
sudo apt-get install libnetpbm10-dev
sudo apt-get install libudunits2-0 libudunits2-dev udunits-bin

# Instalation of netcdf4 in R

wget http://cran.r-project.org/src/contrib/ncdf4_1.13.tar.gz
# it is nescesary to deffine the library
sudo R CMD INSTALL --configure-args="--with-nc-config=/home/demiurgo/local/bin/nc-config" ncdf4_1.13.tar.gz
