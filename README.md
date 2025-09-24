# Parallel_Merge_Sort_in_Pthread

Implementation of the Parallel Merge Sort algorithm in Pthreads.


Compile: 
```console
gcc -O2 pth_msort_test.c pth_msort.c -lpthread -lm
```
Execute: 
```console
./a.out M
```

Note that $N=2^M$ and $M$ can be between $24$ and $30$.
