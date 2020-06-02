# [WIP] - A C lang implementation of Fresnel transform for disk shaped 2D binary signal, with WASM compilation support

```bash
emconfigure ./configure --disable-fortran --libdir=/home/nic/CLionProjects/untitled1/emsfftw --prefix=/home/nic/CLionProjects/untitled1/emsfftw
emmake make
emmake make install

emcc wasmfres.c emsfftw/libfftw3.a -o wasmfresnel.js -O2 -s WASM=1 -s "EXTRA_EXPORTED_RUNTIME_METHODS=['ccall']" -s ASSERTIONS=1 -s ALLOW_MEMORY_GROWTH=1 -s EXPORTED_FUNCTIONS="['_fresnelCircle']"

```
