# BeadPack

Particle tracking in beadpacks in two and three dimensions, including advection, diffusion, and reaction processes. 

## Dependencies

A compiler compatible with C++17 is required (for example, an up-to-date version of the gnu compiler g++ will do).
Required external libraries:
- boost : C++ library implementing useful features and extensions (boost.org)
- CGAL (v5 or later) : C++ interpolation library (cgal.org)
- Eigen (v3 or later) : C++ matrix algebra library (eigen.tuxfamily.org)
- GMP : C arbitrary-precision arithmetic library (gmplib.org)
- MPFR : C multiple-precision arithmetic with correct rounding library (mpfr.org)
- nanoflann : C++ KDTree searching library (github.com/jlblancoc/nanoflann)

## Documentation

See the manual file Manual.pdf in the man directory for further details and help or contact Tomás Aquino.

## Citations

Please link to this repository (https://github.com/tcAquino/BeadPack). If relevant, please cite the original publication: 
- Aquino, Tomás, and Tanguy Le Borgne. "The chemical continuous time random walk framework for upscaling transport limitations in fluid–solid reactions." Advances in Water Resources 154 (2021): 103981.

## License

License information can be found in the file LICENSE.txt.

## About the author

This repository was developed by Tomás Aquino and first uploaded as part of the ChemicalWalks Marie Sk􏰀lodowska Curie Individual Fellowship project, funded by the European Commission under the H2020 programme (project ID: 838426, project website: https://chemicalwalks.wordpress.com/).

The author thanks Joris Heyman, who was closely involved in the development and testing.

Tomás Aquino can be contacted at tomas.c.aquino@gmail.com.
