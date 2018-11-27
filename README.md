# parallel-matrix-multiplication
To compute |X^k| in parallel with Cannon’s method.
——————————
1. This lab assignment is a matrix multiplication (MM) problem using Cannon’s method. The matrix multiplication is fundamental to many high-performance applications. 
2. We are getting familiar with MPI through this assignment. We also find so many useful MPI calls which are need to be used in usual MPI program. We used MPI_Cart_create—makes a new communicator to which Cartesian topology information has been attached in this lab.
3. Cannon algorithm is an excellent parallel algorithm. On multi-CPU processors, Cannon algorithm can improve computational efficiency and utilize computer resources efficiently.
4. Through this lab, we also learned the differences between two-dimensional array and the pointer, including the creation of dynamic two-dimensional arrays, initialization, and passing as formal parameters. The formal parameters are a pointer to an array and the length of the array.
5. We learned to solve some compile problems. For example, when using the head file #include <math.h>, we had to add -lm to compile the .c file.
