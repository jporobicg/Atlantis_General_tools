#!/bin/sh
sudo mkdir /usr/programs/ncBrowser/
cd /usr/programs/ncBrowser/
wget ftp://ftp.epic.noaa.gov/java/ncBrowse/ver140jars/ncBrowse.jar
wget ftp://ftp.epic.noaa.gov/java/ncBrowse/ver140jars/lashandler.jar
wget ftp://ftp.epic.noaa.gov/java/ncBrowse/ver140jars/noaa_pmel.jar
wget ftp://ftp.epic.noaa.gov/java/ncBrowse/ver140jars/netcdfAll.jar
wget ftp://ftp.epic.noaa.gov/java/ncBrowse/ver140jars/visad.jar
wget ftp://ftp.epic.noaa.gov/java/ncBrowse/ver140jars/png.jar

#crear un archivo con estas caracteristicas en > /usr/local/bin
#!/bin/sh
#
# JAVA_HOME is the installation directory
# for Java 2 version 1.4
#
PATH=$JAVA_HOME/bin:$PATH
#
# classpath for Unix cannot have spaces between components
#
classpath="/usr/programs/ncBrowser/ncBrowse.jar:\
/usr/programs/ncBrowser/visad.jar:/usr/programs/ncBrowser/noaa_pmel.jar:\
/usr/programs/ncBrowser/netcdfAll.jar:/usr/programs/ncBrowser/dods.jar:\
/usr/programs/ncBrowser/png.jar:/usr/programs/ncBrowser/lashandler.jar"

java -cp $classpath ncBrowse.Browser $1

chmod +x /usr/local/bin/ncBrowse
