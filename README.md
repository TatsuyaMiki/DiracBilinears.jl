# DiracBilinears.jl 

[![Build Status](https://github.com/TatsuyaMiki/DiracBilinears.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/TatsuyaMiki/DiracBilinears.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://TatsuyaMiki.github.io/DiracBilinears.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://TatsuyaMiki.github.io/DiracBilinears.jl/dev/)



This package provides tools for computing microscopic physical quantities in relativistic quantum theory called Dirac bilinears, which are the fundamental 16 independent quantities derived from the Dirac field.
By connecting to the external first-principles calculation package Quantum ESPRESSO, Wannier90, and wan2respack, it can evaluate 
- spatial distributions of the bilinears
- Wannier matrix elements of the bilinears

## Folder structure

```
  |--LICENSE         
  |--README.md       
  |
  |--docs
  |
  |--examples
  |    |
  |    |--README.md
  |    |--BaTiO3
  |    |--Te
  |
  |--src
  |
  |--test
```

The src directory includes the following files:
- DiracBilinear.jl (The main code.)
- density.jl (The code for spatial distribution.)
- wannier_matrix.jl (The code for Wannier matrix elements.)
- for_wan2respack.jl, read_files.jl (The code for connecting to Quantum ESPRESSO and wan2respack.)
- module.jl (Additional modules)

## Requirements
You will need the following packages.

- [Quantum ESPRESSO](https://www.quantum-espresso.org)
- [Wannier90](http://www.wannier.org)
- [wan2respack](https://github.com/respack-dev/wan2respack/tree/spinor) ("spinor" branch on GitHub)

## Installation

DiracBilinears.jl can be installed with the Julia package manager:
```
julia -e 'import Pkg; Pkg.add("DiracBilinears")'
```

## Usage

### Spatial distribution

The following files are required:
- seedname.scf.in (The input file for scf calculation.)
- seedname.nscf.in (The input file for nscf calculation.)

Before evaluation of the bilinears, you need to perform DFT calculations by [Quantum ESPRESSO] to obtain the Bloch functions.
```sh
QE/bin/pw.x < seedname.scf.in > seedname.scf.out
QE/bin/pw.x < seedname.nscf.in > seedname.nscf.out
```
[Quantum ESPRESSO]: https://www.quantum-espresso.org


To compute the spatial distribution of electron chirality, use the following commands:
```Julia
using DiracBilinears
τz = calc_density(calc="τz", qedir="$outdir/seedname.save")
```
Here, ``qedir`` should be set to the output directory from the Quantum ESPRESSO calculation, which contains the Bloch wave functions.
The variable ``$outdir`` corresponds to the directory specified in the nscf input file.
You can select the physical quantity to calculate by specifying the ``calc`` option. 
The available options are listed in the table below.
|``calc`` | Physical quantities |
|:---:|:---:|
|"ρ" (or "rho")|| electron density |
|"ms"| magnetization |
|"j"| current |
|"∇ms" (or "nabla_ms")| pseudoscalar|
|"τz" (or "tau_z")| electron chirality|
|"∇ρ" (or "nabla_rho")| gradient part of polarization|
|"ps"|spin-derived electric polarization|  

You can generate the xsf file for [XCrySDen] plot after the ``calc_density()`` calculation:
```Julia
write_density(τz; qedir=QEDIR, savefile=FILENAME)
```
The file extension of FILENAME must be set to ".xsf".

[XCrySDen]: http://www.xcrysden.org/



### Wannier matrix element

The following files are required:
- seedname.scf.in (The input file for scf calculation.)
- seedname.pw2wan.in (The input file for pw2wannier90.x.)
- seedname.win.ref (The reference file for generating the input file of Wannier90.)
- conf.toml (The configuration file for wan2respack.)

For the Wannier matrix elements, you first need to perform DFT calculations using [Quantum ESPRESSO] to obtain the Bloch functions. 
You can do this by running:
```sh
QE/bin/pw.x < seedname.scf.in > seedname.scf.out
```
After completing the scf calculation, you can generate the Wannier functions using [Wannier90] and [wan2respack] ([RESPACK])　(for more information, see the GitHub page of [wan2respack]). 
First, run:
```sh
python Wan2respack/bin/wan2respack.py -pp conf.toml
```
This command will generate an input file ``seedname.nscf_wannier.in``. 
In this file, please replace ``calc="scf"`` with ``calc="nscf"``.
Next, calculate the Wannier functions by running:
```
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


To compute the matrix elements, you need to specify the paths of two directories, "dir-wfn-full" (referred to as ``WFNDIR``) and "dir-wan" (referred to as ``WANDIR``). 
These directories are generated by the [wan2respack] calculation. 
To calculate the matrix elements of electron chirality, use the following commands:
```Julia
using DiracBilinears
rs, degen = calc_rgrid(rfile=FILENAME)
τz = calc_wannier_matrix(calc="τz", rgrid=rs, wfndir=WFNDIR, wandir=WANDIR)
```
Please set ``rfile`` to the input file of Wannier90 or the input file of scf calculation. 
You can choose the type of ``calc`` from the table shown below.
|``calc`` | Physical quantities |
|:---:|:---:|
|"ρ" (or "rho")| electron density |
|"ms"| magnetization |
|"j"| current |
|"τz" (or "tau_z", "chirality")| electron chirality|
|"ps"| spin-derived electric polarization| 

The Wannier matrix elements can be saved to a text file using the following command:
```Julia
write_wannier_matrix(τz, rs, degen; savefile=FILENAME)
```


## License and citation

This package is released under the GNU General Public License version 3 (GPL v3). See ``LICENSE`` for details.

If you find this package useful in your research, please consider citing the following papers:
- T. Miki, H. Ikeda, M.-T. Suzuki, and S. Hohisno, [Phys. Rev. Lett. **134**, 226401 (2025)](https://doi.org/10.1103/PhysRevLett.134.226401).
- T. Miki, H.-Y. Chen, T. Koretsune, Y. Nomura, [Comput. Phys. Commun. **317**, 109857 (2025)](https://doi.org/10.1016/j.cpc.2025.109857).

