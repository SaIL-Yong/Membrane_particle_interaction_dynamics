### **Step-1: Copying/Downloading the source code**



Clone the repository in you local machine/cluster using the following command line:

git clone [https://github.com/SaIL-Yong/Membrane\_particle\_interaction\_dynamics.git](https://github.com/SaIL-Yong/Membrane_particle_interaction_dynamics.git)



### **Step-2: Compiling and Building Executable**



1. cd "path\_to\_git\_repo/code\_base"
2. mkdir build
3. cd build
4. cmake ..
5. make -j N ### N=1.4.8, ...



This process will generate an executable **MemDynamics**



### **Step-3: Running Simulation**



1. Copy MemDynamics to the folder you want to run the simulation
2. Save **run\_parameters.txt**  file in the same folder (For detailed parameter choice and selection, see our manuscript)
3. Save your particle mesh file and vesicle mesh file in **.off** format

4\. For running on the local machine: type in the command line : ./MemDynamics

