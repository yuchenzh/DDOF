import cantera as ct 
import numpy as np



# ===== Setting ====== #
mech = "./h2.yaml"
Tlow = 300
Thigh = 3000
grids = 100


#---------------------#

def transportTemplate(species_name,muPoly,kappaPoly,oneOverKappaPoly):
    muString    = "muCoeffs<8>       ("
    for i in range(len(muPoly) + 1):
        muString += f"{muPoly[i]:.5e}" + " "
    muString += ");\n"

    kappaString = "kappaCoeffs<8>    ("
    for i in range(len(kappaPoly) + 1)  :
        kappaString += f"{kappaPoly[i]:.5e}" + " "
    kappaString += ");\n"
    
    oneOverKappaString = "oneOverKappaCoeffs<8>    ("
    for i in range(len(kappaPoly) + 1)  :
        oneOverKappaString += f"{oneOverKappaPoly[i]:.5e}" + " "
    oneOverKappaString += ");\n"
    
    string = ""
    string += species_name + "\n"
    string += "{\n"
    string += "    transport\n"
    string += "    {\n"
    string += "        " + muString
    string += "        " + kappaString
    string += "        " + oneOverKappaString
    string += "    }\n"
    string += "}\n"
    return string

with open("thermoPolynomial","w+") as f:
    gas = ct.Solution("./h2.yaml")
    for species in gas.species_names:
        reactants_species = "{}:1".format(species)
        TList = np.linspace(Tlow, Thigh, grids)
        muList = []
        kappaList = []
        oneOverKappaList = []
        for T in TList:
            gas.TPY = T, 101325,reactants_species
            gas.transport_model = "mixture-averaged"
            mu = gas.viscosity
            kappa = gas.thermal_conductivity
            muList.append(mu)
            kappaList.append(kappa)
            oneOverKappaList.append(1./kappa)
            
        pmu = np.polyfit(TList, muList, 7)
        pmu = np.poly1d(pmu)
        pkappa = np.polyfit(TList, kappaList, 7)
        pkappa= np.poly1d(pkappa)
        poneOverKappa = np.polyfit(TList, oneOverKappaList, 7)
        poneOverKappa = np.poly1d(poneOverKappa)
        
        #print(transportTemplate(species,pmu,pkappa,poneOverKappa))

        f.write(transportTemplate(species,pmu,pkappa,poneOverKappa))
    