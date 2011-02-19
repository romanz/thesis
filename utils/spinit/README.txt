SPINIT - sparse matrix faster initialization utility
====================================================

INTRODUCTION
------------
This utility can be useful when you need to update a sparse matrix values 
without changing its sparsity pattern.

Since MATLAB uses CSC (Compressed Sparse Column) format for sparse matries, 
the non-zeroes are stored in an array, which can be safely overwritten in
order to update the matrix, while preserving the sparsity pattern.


INSTALLATION
------------
Run make.m to compile and build MEX extension.

USAGE
-----
See spinit.m for usage example.

AUTHOR
------
Roman Zeyde (roman.zeyde@gmail.com)
