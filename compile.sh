if [[ -f "./main" ]]; then
  rm ./main
fi
# g++  -I/home/adam/OneDrive/PhD/Programming/C++/Numerical\ Integrator/Headers -I/usr/include/python3.11/ -lpython3.11 ./main.cpp -O3 -o main -lboost_iostreams -lboost_system -lboost_filesystem
g++ -I./Headers -I./Examples ./main.cpp -O3 -o main -lboost_iostreams -lboost_system -lboost_filesystem
if [[ -f "./main" ]]; then
  ./main
fi
