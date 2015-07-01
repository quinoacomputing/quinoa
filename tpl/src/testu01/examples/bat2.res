========= Summary results of Rabbit =========

 Version:          TestU01 1.2.3
 File:             vax.bin
 Number of bits:   1048576
 Number of statistics:  38
 Total CPU time:   00:00:01.75
 The following tests gave p-values outside [0.001, 0.9990]:
 (eps  means a value < 1.0e-300):
 (eps1 means a value < 1.0e-15):

       Test                          p-value
 ----------------------------------------------
  4  AppearanceSpacings              1.1e-4
  7  Fourier1                         eps  
  8  Fourier3                      3.2e-213
 13  HammingCorr, L = 64            1 - eps1
 16  HammingIndep, L = 32             eps  
 17  HammingIndep, L = 64             eps  
 24  RandomWalk1 M                    eps  
 24  RandomWalk1 J                    eps  
 ----------------------------------------------
 All other tests were passed
