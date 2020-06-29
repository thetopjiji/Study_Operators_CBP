# Study_Operators_CBP

The article is already accepted by the conference Artificial Evolution 2019 (Mulhouse, France).

This work is dedicated to investigating the influences of different operators for CBP and we try to find the most suitable operator in memetic algorithm framework for this problem.


The source code is shared.

Compile: 
g++ BCP_mem.cpp BCP_main.cpp BCP_X.cpp -O3 -lm -Wall -o BCP_1

Run:
BCP_1 -i V-nos6.mtx.rnd --seed 0 -rep 0 -alb 3 -pop 20
