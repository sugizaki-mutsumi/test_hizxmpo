# Example of HiZ WXM MPO x-ray tracing
-  X-ray reflection physics codes are based on XRTG4.
 https://xraytracer.com
- Source codes in source/src are modified for G4 11.3 and HiZ WXM MPO.

##  Requirements
- Geant4 ver 11.3
- ROOT (v6.32.02)
- (python numpy, astropy, matplotlib, deveoloped in anaconda enviromnent) 

## How to build
```
 $ mkdir build
 $ cd build
 $ cmake ../source
 $ make
 $ make install
```

## How to run
### Move to TestBench directory
```
 $ cd TestBench
```
### Edit parameters in macro file *test_hizxmpo.mac* 
 - gamma-ray source shape, dimension, energy
 - output filename
 - number of beam photons

### Run in Batch mode
 ```
 $ ../bin/test_hizxmpo test_hizxmpo.mac
 (takes a few minutes to simulate 100000 events)
```

### Plot PSF
```
 $ python -i plt_psfsct.py
 or
 $ jupyter-lab plt_psfsct.ipynb
  
```

