sudo apt-get install perl python3 python3-pip cpanminus wget autoconf automake libtool ghostscript
cpanm Graph
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.17.tar.gz
tar -zxvf ViennaRNA-2.4.17.tar.gz
cd ViennaRNA-2.4.17
./configure
make
sudo make install