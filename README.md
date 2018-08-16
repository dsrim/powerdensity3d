## Reconstruction of isotropic and anisotropic conductivities from power densities in three dimensions

This repository contains MATLAB code illustrating the reconstruction algorithms that recovers the isotropic
and anisotropic conductivity from power density measurements, on a cube-shaped domain. This code was
developed and tested in MATLAB version R2016b.

``run/`` directory contains two runfiles 
* ``runIso`` implements the isotropic reconstruction
* ``runAniso`` implements the anisotropic reconstruction

Both of the runfiles
* computes the solutions and power densities using FDM (5-pt stencil) or FEM (quadratic, PDE Toolbox)
* computes the reconstruction using methods outlined in the referenced paper

``plot/`` directory contains plotting tools of the output. For example, all figures in the referenced paper 
can be reproduced by running ``plotExp1.m``, ``plotExp2.m``, ``plotExp3.m``.

### Required software

* MATLAB PDE Toolbox for solving the anisotropic forward problem.
* Manifold-valued Image Restoration Toolbox (MVIRT) for plotting anisotropic coefficient.
 A modified version of ``helpers/plotSPD.m`` is used.
 [[Github repo]](https://github.com/kellertuer/MVIRT)

### References

* [[*Imaging of isotropic and anisotropic conductivities from power densities in three dimensions*]](http://iopscience.iop.org/article/10.1088/1361-6420/aabe5a/meta)
<br> F. Monard and D. Rim <br> *Inverse Probl.* (2018) **34** (7), 075005 [[arXiv]](http://arxiv.org/abs/1711.03137)
* *MVIRT, A toolbox for manifold-valued image registration.*
<br>Bergmann, R (2017).  <br>
IEEE International Conference on Image Processing, IEEE ICIP 2017, Beijing, China, September 17–20, 2017
 [[Github repo]](https://github.com/kellertuer/MVIRT)
* *A second order non-smooth variational model for restoring manifold-valued images* <br>
M. Bačák, R. Bergmann, G. Steidl, A. Weinmann (2016). <br>
SIAM Journal on Scientific Computing. 38, (1), A567–A597. [[doi]](http://dx.doi.org/10.1137/15M101988X) [[www]](http://arxiv.org/pdf/1506.02409v2.pdf)


[![DOI](https://zenodo.org/badge/115078808.svg)](https://zenodo.org/badge/latestdoi/115078808)
