# cubie

Cubie
-----

C code for generating a cubed sphere, see article (http://dx.doi.org/10.1006/jcph.1996.0047)

    C. Ronchia, R. Iaconoa, P.S. Paoluccib, The “Cubed Sphere”: A New Method for the Solution of Partial Differential Equations in Spherical Geometry, Journal of Computational Physics, Volume 124, Issue 1, 1 March 1996, Pages 93–114.

It outputs a globe.vtu file that defines the sphere with:
 - vertex coordinate information
 - vertex connectivity 
 - boundary vertices
 - partition information (optional if METIS is used to compile code)

The globe.vtu file is viewable in paraview. 

Optional
--------
The code optionally use the METIS library and development header to create a partition on the mesh. To use first install METIS development libraries, set the appropriate METIS path in the makefile and tweak the c code.

Compile
-------
To compile:
    'make'

To run:
    './go -l 11 -d 6 -r 10 -i 5'

where:
l - number of elements across a sixth 
d - number of elements in depth
r - the outer radius
i - the inner radius

To clean the output:
    'make clean'
