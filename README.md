# Study_Operators_CBP

The article is already accepted by the conference Artificial Evolution 2019 (Mulhouse, France).

This work is dedicated to investigating the influences of different operators for CBP and we try to find the most suitable operator in memetic algorithm framework for this problem.


The source code is shared.

Compile: 
g++ BCP_mem.cpp BCP_main.cpp BCP_X.cpp -O3 -lm -Wall -o BCP_1

Run:
BCP_1 -i V-nos6.mtx.rnd --seed 0 -rep 0 -alb 3 -pop 20

--seed 0: the seed is set to 0, you could change it as you like. It is just for reproducing the results.
-rep 0  : not a parameter for the algorithm, you could set it as 0
-alb 3  : the lower bound of the graph
-pop 20 : the number of the population
