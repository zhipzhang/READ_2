#INSTALL
First tar -zxvf hessioxxx.tar.gz && cd hessioxxx && make

Then set the HESS in CMakeLists to the ./hessioxxx/lib/libhessio.so

mkdir build && cmake ..

cmake --build .

Finally mv libclass.pcm to the compiled/lib

excutable will be in compiled/bin
