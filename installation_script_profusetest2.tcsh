git clone https://github.com/tedholzman/profillicSimulation.git
cd profillicSimulation
ln -s ../HMMoC-BFloat-Algebra
ln -s ../prolific
ln -s ../profillic
ln -s ../profuse
ln -s ../seqan-trunk
ln -s ../hmmer
ln -s /usr/local/lib boost-lib
ln -s /usr/local/include boost-include
make profusetest2
profusetest2