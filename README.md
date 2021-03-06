# HW3_Fixed-outline-Floorplan-Design_B-star-tree

This is a floorplaner with floorplan represented by **B-star tree** structure (http://cc.ee.ntu.edu.tw/~ywchang/Papers/boundary-btree.pdf) using **simulated annealing algorithm**.<br> 

A minor difference is that the cost function is modified in order to implement the fix-outline constraint. The cost function includes the W and H penalty when violating FO. And also the aspect ratio of the floorplan is added into the cost function.<br>

## How to compile <br>
In directory ```./src```, enter the following command, <br>
```
$ make
```
It will generate the executable file ```./main``` in ```./HW3/bin/```. <br>
  
To remove it, please enter the following command, <br>
```
$ make clean
```
## How to Run
In directory ```./src```, enter the following command, <br>
```
$ ../bin/[exe] [hardblocks file] [nets file] [pl file] [dead space_ratio]
```
E.g.
```
$ ../bin/hw3 ../testcase/n100.hardblocks ../testcase/n100.nets ../testcase/n100.pl 0.1 
```
