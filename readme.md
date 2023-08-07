# Using the builtin differential diffusion in OpenFOAM-10
## Introduction
The differential diffusion transportation is availble since OpenFOAM-9. However, it is never easy to use that. A LONG and LARGE table for diffusion coefficients should be generaetd beforehand, which can be time-costing. Besides, there is no tutorial (until 2023/08/07) on the internet showing how to make this works. This repository will offer you a easy way to use that differential diffusion functionality.

Generally, we offer two methods. 
1. Generate that crazy table for you with a piece of Python code, using the diffusion values from Cantera. <span style="color:red;">DNS and RANS is supported natively in OpenFOAM, our library can add support for LES.</span>

2. Generate a much simplified table containing polynomial coefficients fitting the binary diffusion values of all species pair, and let the OpenFOAM compute the diffusion coefficients during runtime. <span style="color:red;"> Our library add supports for DNS, LES and RANS</span>

Note, the method 2 is definitely more accurate, but might be slower since recalculating the coefficients is more expensive than looking up the table.

## Compile the code
You should, of course, install OpenFOAM-10 and enable the environment.
Go to the repository directory and run
```bash
wmake fluidReactionThermo
```

Then everything is ready.

## Usage
### Generate differential diffusion table.
- Our Python class require the Cantera for Python. Please install [Cantera](https://cantera.org/install/index.html) beforehand. Install through conda is the recommended way.
- Convert the reaction mechanism into the cantera format. You might need to change      the file name.
```bash
ck2yaml --input=chem.inp --thermo=therm.dat --transport=tran.dat
```
- Put the generated mechanism file(usually 'mech.yaml') into the 'ddtransport'          folder. If you install with conda, please make sure you enable the environment with Cantera. Then run 
    ```Python
    python getDiffTable.py
    ```
    to generate the complex table (for method 1).
    
    Or run
    ```Python
    python getDiffPoly.py
    ```
to generate the simplified table (for method 2).

### Put the table in the constant folder of your case

### Configure the thermophysicalTransport file
See the following example. You can find more examples from 'Cases' folder
```cpp
laminar
{
    // FickianFourier, TabulatedFickianFourier for DNS
    // FickianEddyDiffusivity, TabulatedFickianEddyDiffusivity for RANS and LES
    // Add "Tabulated" to use a simplified table
    model          FickianFourier;

    mixtureDiffusionCoefficients no;
    
    Prt             0.7;
    Sct             0.7;

    // D           for method 1, complex table
    // DPolyCoeffs for method 2, simplified table
    DPolyCoeffs // [m^2/s] 
    {
        // #include "transportProperties"         // Use method 1, complex table
        #include "transportPropertiesPolynomials" // Use method 2, simplified table
    }
}
```
### Don't forget to add the library to your 'controlDict' file
```cpp
// For simplified or complex
libs ("libtabulatedFluidReactionThermophysicalTransportModels.so");
```

## Examples
We offer a 1D hydrogen flame case using DNS, LES and RANS model, with both simplified and complex tables.
- 'Cases/DNS-builtin'  (DNS, complex table)
- 'Cases/DNS'          (DNS, simplified table)
- 'Cases/LES-builtin'  (LES, complex table)
- 'Cases/LES'          (LES, simplified table)
- 'Cases/RANS-builtin' (RANS, complex table)
- 'Cases/RANS'         (RANS, simplified table)

The complex table is too large to be uploaded to Github. Please generate by yourself.

In LES/RANS cases, the turbulent kinetic energy is set to be zero, so there is no turbulent diffusion. The cases are configured for testing purposes.

The result is accurate.


## Another OpenFOAM solver supporting differential diffusion.
See [reacitngDNS](https://github.com/ZSHtju/reactingDNS_OpenFOAM). reactingDNS (in OpenFOAM-8) is accurate and well-tested in many cases.

In this repository, the differential diffusion is achieved in a library, instead of a solver. Consequently, you can enable it in all reacting solvers, and integrate it with other functionalities (e.g. turbulence, dynamic mesh, spray clouds...)


## Reference
<span style="color:red;"> to be added in the future</span>