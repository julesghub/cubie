# cubie
C code for generating a cubed sphere. 
It outputs the sphere as a globe.vtu file, viewable in paraview.

This can optionally use the METIS library and development header to create a partition on the mesh. To use this install METIS, set the appropriate METIS path in the makefile and tweak the c code.

To compile:
    'make'

To run:
    './go -l 11 -d 6 -r 10 -i 5'

where:
l is the length resolution 
d is the depth resolution
r is the outer radius
i is the inner radius
