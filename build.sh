projec_tName=$1
echo $projec_tName
c++ -o ./run/${projec_tName:-miniRayTracer} ./code/${projec_tName:-miniRayTracer}.cpp -O3 -std=c++11