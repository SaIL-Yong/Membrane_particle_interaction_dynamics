###Membrane_Dynamcs_Code_Development
To run the code in Unix/mac system, run the following commands:
mkdir build

cd build

cmake ..

cmake --build .

to compile a debug version

mkdir debug

cd debug

cmake -DCMAKE_BUILD_TYPE=Debug ..

cmake --build .

###Two existing branch : 1) For membrane dynamics/energy minimization only.
###                      2) For Particle Interaction.
