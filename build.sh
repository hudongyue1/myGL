# build
project_Name=$1
echo $projec_tName
c++ -o ./run/${project_Name:-miniRayTracer} ./code/${project_Name:-miniRayTracer}.cpp -O3 -std=c++11

# run
./run/${project_Name:-miniRayTracer}