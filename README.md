# HW3_Fixed-outline-Floorplan-Design_B-star-tree

B star tree algorithm with simulated annealing.<br> 

A minor difference is that I modified the cost function in order to implement the fix-outline constraint. My cost function includes the W and H penalty when violating FO. And also the aspect ratio of the floorplan is added into the cost function.

How to compile-
```
$ make
```
```
$ ../bin/hw3 ../testcase/n100.hardblocks ../testcase/n100.nets ../testcase/n100.pl ../output/n100.floorplan 0.1
```
