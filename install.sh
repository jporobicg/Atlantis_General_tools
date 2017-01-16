#!/bin/sh
#--------------#
#   PASO 1     #
#--------------#

#  instalar desde repositorio
#   sudo apt-get install gmt

#  o
        sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
        sudo apt-get update
        sudo apt-get install gmt
        sudo apt-get install gmt-tutorial
        sudo apt-get install gmt-doc
        sudo apt-get install gmt-example


#--------------#
#   PASO 2     #
#--------------#


#  Path for netcdf & GMT

#  escribir en la terminal
 emacs ~/.bashrc
# Pegar esto al final del documento

        export NETCDFHOME=/usr/lib
        NETCDF_PREFIX=$NETCDFHOME
        GMTHOME=/usr/lib/gmt
        export GMTHOME
        PATH="$PATH:$GMTHOME/bin"
        export PATH
        MANPATH="$MANPATH:$GMTHOME/man"
        export MANPATH

#  Download coastal lines

wget ftp://ftp.iag.usp.br/pub/gmt/gshhg-gmt-2.3.4.tar.gz

#  Uncompress coastal lines

tar -zxvf gshhg-gmt-2.3.1.tar.gz

#  Move to 'usr' folder
mv  -v /gshhg-gmt-2.3.1/binned_*.nc /usr/share/gmt/coast
