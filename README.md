# DiracBilinears.jl 

[![Build Status](https://github.com/TatsuyaMiki/DiracBilinears.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/TatsuyaMiki/DiracBilinears.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://TatsuyaMiki.github.io/DiracBilinears.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://TatsuyaMiki.github.io/DiracBilinears.jl/dev/)



This package provides tools for computing microscopic physical quantities in relativistic quantum theory called Dirac bilinears, which are the fundamental 16 independent quantities derived from the Dirac field.
By connecting to the external first-principles calculation package Quantum ESPRESSO, Wannier90, and wan2respack, it can evaluate 
- spatial distributions of the bilinears
- Wannier matrix elements of the bilinears

## Installation

> [!CAUTION]
> This part will be replaced with instructions for the official Julia package manager.

You can download DiracBilinears.jl from this repository:
```sh
git clone https://github.com/TatsuyaMiki/DiracBilinears.jl.git
```
To install the package, run the following commands:
```sh
cd DiracBilinears.jl
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```
<!-- 
[TO BE REPLACED]
DiracBilinears.jl can be installed with the Julia package manager as
```
julia -e 'import Pkg; Pkg.add("DiracBilinears")'
```
-->

## Usage

### Spatial distribution

The following files are required:
- seedname.scf.in
- seedname.nscf.in

Before evaluation of the bilinears, you need to perform DFT calculations by [Quantum ESPRESSO] to obtain the Bloch functions.
```sh
QE/bin/pw.x < seedname.scf.in > seedname.scf.out
QE/bin/pw.x < seedname.nscf.in > seedname.nscf.out
```
[Quantum ESPRESSO]: https://www.quantum-espresso.org


To compute the spatial distribution of electron chirality $\tau^Z(\boldsymbol{r})$ on a uniform $32\times 32\times 32$ spatial mesh ``nrmseh=(32,32,32)``, use the following commands:
```Julia
using DiracBilinears
τz = calc_density(calc="τz", nrmesh=(32,32,32), qedir=QEDIR)
```
The physical quantities to be calculated can be selected by specifying the option ``calc``, whose list is provided in the following table.
|``calc`` | Physical quantities |
|:---:|:---:|
|"ρ" (or "rho")|| electron density |
|"ms"| magnetization |
|"j"| current |
|"∇ms" (or "nabla_ms")|pseudoscalar|
|"τz" (or "tau_z")|electron chirality|
|"∇ρ" (or "nabla_rho")| dradient part of polarization|
|"ps"|spin-derived electric polarization|  

You can generate the xsf file after the ``calc_density()`` calculation:
```Julia
write_density(τz; qedir=QEDIR, savefile=FILENAME)
```
The file extension of FILENAME must be set to ".xsf".




### Wannier matrix element

The following files are required:
- seedname.scf.in
- seedname.pw2wan.in
- seedname.win.ref
- conf.toml

For Wannier matrix elements, you need to perform DFT calculations by [Quantum ESPRESSO] to obtain the Bloch functions.
```sh
QE/bin/pw.x < seedname.scf.in > seedname.scf.out
```
Then, you need to generate the Wannier functions using [Wannier90] and [wan2respack] ([RESPACK]) after the scf and the nscf calculations.
The following command is for these two packages (for more information, see the GitHub page of [wan2respack]):
```sh
python Wan2respack/bin/wan2respack.py -pp conf.toml
QE/bin/pw.x < seedname.nscf_wannier.in > seedname.nscf_wannier.out
Wanier90/wannier90.x -pp seedname
QE/bin/pw2wannier90.x < seedname.pw2wan.in > seedname.pw2wan.out
Wanier90/wannier90.x seedname
python Wan2respack/bin/wan2respack.py conf.toml
```
[Quantum ESPRESSO]: https://www.quantum-espresso.org
[wan2respack]: https://github.com/respack-dev/wan2respack/tree/main
[Wannier90]: https://wannier.org
[RESPACK]: https://sites.google.com/view/kazuma7k6r


To compute the matrix elements, the path of two directories dir-wfn-full (``WFNDIR``) and dir-wan (``WANDIR``) need to be specified.
These directories are generated by [wan2respack] calculation.
The following command is used to calculate the matrix elements of electron chirality for $5\times 5 \times 5$ unit cells:
```Julia
using DiracBilinears
rs, degen = calc_rgrid(mpmesh=(5,5,5))
τz, rs = calc_wannier_matrix(calc="τz", rgrid=rs, wfndir=WFNDIR, wandir=WANDIR)
```
You can select the option ``calc`` from the following table.
|``calc`` | Physical quantities |
|:---:|:---:|
|"ρ" (or "rho")| electron density |
|"ms"| magnetization |
|"j"| current |
|"τz" (or "tau_z", "chirality")|electron chirality|
|"ps"|spin-derived electric polarization| 

The Wannier matrix elements can be saved to a text file using the following command:
```Julia
write_wannier_matrix(τz, rs; savefile=FILENAME)
```
