# How to use it

c++ -std=c++11 -I./include -O3 src/grid.cpp -o a.out <br>
./a.out structure/KAXQIL_clean_14.cif forcefield/UFF.def 298.0 12.0 Xe 0.12 40

Widom value of enthalpy for KAXQIL_clean : -44.6297 kJ/mol

# Benchmark
c++ -std=c++11 -I./include -O3 src/ads_grid.cpp -o ads_grid.out

./ads_grid.out structure/KAXQIL_clean_14.cif forcefield/UFF.def 298.0 12.0 Xe 0.12 20 0.8 <br>
KAXQIL_clean_14,-44.626,0.0300725,0.190666 <br>
./ads_grid.out structure/IFY.cif forcefield/zeolite_ff.def 298.0 12.0 Xe 0.12 20 0.8 <br>
IFY,-26.1962,6.81103e-05,0.205935 <br>
./ads_grid.out structure/AEI.cif forcefield/zeolite_ff.def 298.0 12.0 Xe 0.12 20 0.8 <br>
AEI,-22.2866,5.06784e-05,0.227736 <br>
./ads_grid.out structure/IFY.cif forcefield/zeolite_ff.def 298.0 12.0 Xe 0.12 20 0.8 <br>
IFY,-26.1962,6.81103e-05,0.216043 <br>

# Acknowledgement
