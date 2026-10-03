# DDOF: Differential Diffusion for OpenFOAM

## Introduction
OpenFOAM supports differential (species-specific) diffusion through the `FickianFourier` thermophysical transport model, but it requires binary diffusion coefficients for every species pair, which are tedious to prepare. This repository provides Python tools (based on Cantera) that generate these coefficients directly from a reaction mechanism, together with example cases. In OpenFOAM-10, an additional library is supplied to extend support to LES and tabulated coefficients. In OpenFOAM-14, the native implementation is used and the workflow is much simpler.

## Versions
- [OpenFOAM-10](OpenFOAM-10/readme.md)
- [OpenFOAM-14](OpenFOAM-14/readme.md)
