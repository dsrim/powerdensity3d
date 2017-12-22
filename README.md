## Reconstruction of isotropic and anisotropic conductivities from power densities in three dimensions

This repository contains MATLAB code illustrating the reconstruction algorithms that recovers the isotropic
and anisotropic conductivity from power density measurements, on a cube-shaped domain.

``run/`` directory contains two runfiles 
* ``runIso`` implements the isotropic reconstruction
* ``runAniso`` implements the anisotropic reconstruction

Both of the runfiles
* computes the solutions and power densities using FDM (5-pt stencil) or FEM (quadratic, PDE Toolbox)
* computes the reconstruction using methods outlined in the referenced paper

``plot/`` directory contains plotting tools of the output. For example, all figures in the referenced paper 
can be reproduced by running ``plotExp1.m``, ``plotExp2.m``, ``plotExp3.m``.

### References

*Imaging of isotropic and anisotropic conductivities from power densities in three dimensions*
<br> F. Monard and D. Rim <br> Preprint (2017) [arXiv](http://arxiv.org/abs/1711.03137)
