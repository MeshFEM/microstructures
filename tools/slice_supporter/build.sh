g++ -O3 -I$MeshFEM -std=c++11 RemoveUnsupported.cc -I/opt/local/include -L/opt/local/lib -lX11 -o RemoveUnsupported
g++ -O3 -std=c++11 SliceSupporter.cc -I/opt/local/include -L/opt/local/lib -lX11 -o SliceSupporter
