# geant_xray_detector

## To build:
```
source <your geant4 build path>/bin/geant4.sh
cd build
cmake ../
make -j8
cd ../
```

## To run:
```
build/main photon_sweep.mac
```
## To analyze output, or to change sweep parameters:
look at photon_angular_response_sweep.ipynb
