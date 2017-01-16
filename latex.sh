#!/bin/sh
## descargar Latex intstaler

wget http://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz
mkdir /home/demiurgo/Documents/.latex

#  extraer archivos
# instalar LAtex
cd  /home/demiurgo/Documents/.latex/install-tl-20150209/
 sudo ./install-tl  #cambiar las opciones
PATH=/usr/local/texlive/2014/bin/x86_64-linux:$PATH

PATH=/usr/local/texlive/2014/bin/i386-linux:$PATH; export PATH
# install ghostcript

wget
http://downloads.ghostscript.com/public/binaries/ghostscript-9.15-linux-x86_64.tgz

# problems with Swaeve
