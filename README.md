# MB-PIPNet model
1. The source codes to perform energy and gradient calculations on MB-PIPNet model are provided in "src-getpot" and "src-getpot-pbc"

2. The PIP 1b and 2b used in the example are provided in "PIP_bases" including "bemsa1b.f90" and "bemsa2b.f90".
   These bases can be efficiently and systematically generated using our MSA software: "https://github.com/szquchen/MSA-2.0" and "https://github.com/PaulLHouston/PESPIP"

3. The NN structures are specified in "NN_struct.dat" and trained parameters are in "WB01.txt"

4. Here are the instructions to run the MB-PIPNet model:

A: for gas phase system
   A.1 Go to "src-getpot" folder
   A.2 Check Makefile. We use ifort compiler with mkl libarary. User are free to switch to gfortran.
   A.3 Compile the code through "make getpot.x"
   A.4 In the example, we are calling the MB-PIPNet model for energy and gradient calculation of water trimer
       run the code through "./getpot.x"
   A.5 Example output includes:
 Potential energy (a.u.)
    -0.02487896
 Gradient: Analytical V.S. Numerical, (a.u.)
    -0.00008382    -0.00008377
    
    -0.00043106    -0.00043108
    
     0.00000294     0.00000295
     
     0.00073955     0.00073956
     
     0.00007022     0.00007016
     
    -0.00010496    -0.00010496
    
    -0.00021726    -0.00021728
    
    -0.00055136    -0.00055138
     0.00057574     0.00057580
    -0.00005914    -0.00005914
     0.00024315     0.00024315
    -0.00008335    -0.00008335
     0.00043638     0.00043639
    -0.00004389    -0.00004389
    -0.00086811    -0.00086813
    -0.00060949    -0.00060948
     0.00043136     0.00043134
     0.00055309     0.00055311
    -0.00011503    -0.00011503
     0.00076031     0.00076033
    -0.00002945    -0.00002952
    -0.00000032    -0.00000032
     0.00002193     0.00002193
    -0.00023955    -0.00023956
    -0.00032440    -0.00032442
     0.00041327     0.00041334
    -0.00048676    -0.00048675

B: for condensed phase system
   B.1 Go to "src-getpot-pbc" folder
   B.2 Check Makefile. We use ifort compiler with mkl libarary. User are free to switch to gfortran.
   B.3 Compile the code through "make getpot.x"
   B.4 In the example, we are calling the MB-PIPNet model for energy and gradient calculation of water trimer
       run the code through "./getpot.x 256mer.xyz"
   B.5 Example output includes:
 Potential energy (a.u.)
    -3.93541185
 Gradient: Analytical V.S. Numerical, (a.u.)
     0.00393497     0.00393496
    -0.00463793    -0.00463794
    -0.00219358    -0.00219358
     0.02122628     0.02122625
    -0.02020203    -0.02020197
    -0.00165442    -0.00165442
    -0.01637745    -0.01637742
     0.00145604     0.00145604
     0.01932483     0.01932480
    -0.00085289    -0.00085287
    -0.00064634    -0.00064635
     0.01506490     0.01506496
     0.00314290     0.00314290
    -0.00022201    -0.00022194
     0.01677918     0.01677918
    -0.01533535    -0.01533535
     0.00675339     0.00675338
    -0.00477681    -0.00477683
     0.00025506     0.00025506
    -0.01163454    -0.01163454
    -0.01575561    -0.01575570
    -0.02721358    -0.02721357
    -0.02553368    -0.02553362
     0.00371682     0.00371683
     0.00589336     0.00589335
    -0.00641001    -0.00641001
    -0.00767874    -0.00767874
    -0.01688862    -0.01688858
     0.00079908     0.00079909
     0.00906781     0.00906777
    -0.01311159    -0.01311159
     ......

