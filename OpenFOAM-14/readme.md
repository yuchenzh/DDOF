# Differential diffusion in OpenFOAM-14

OpenFOAM-14 supports differential diffusion natively, so no extra library needs to be compiled. The only thing you need is the binary diffusion coefficient file.

## Generate the coefficient file
Install [Cantera](https://cantera.org/install/index.html) for Python (conda is recommended), put your mechanism (e.g. `mech.yaml`, converted with `ck2yaml` if needed) into `ddtransport`, adjust the settings in `getDiffPoly.py`, and run
```bash
cd ddtransport
python getDiffPoly.py
```
This generates `transportPropertiesPolynomials`. Copy it into the `constant` folder of your case (Remember to adjust the mechanism file).

## Example case
See `Cases/ch4` for details. Key files:
- `constant/transportPropertiesPolynomials`: for transport polynomial coefficients
- `constant/thermophysicalTransport` to set the differential diffusion model


