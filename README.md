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

For generating video from image sequence, using ffmpeg

example syntax

ffmpeg -r 10 -i img%05d.jpg -c:v libx264 -pix_fmt yuvj420p -crf 18 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" video.mp4
-crf (0-51) controls vidoe quality (0 is lossless)
-r is the frame rate (per second)
