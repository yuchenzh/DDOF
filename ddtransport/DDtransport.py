import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import os

# generate logscale array from low to high
def logspace(low, high, stepRatio):
    pArray = []
    pArray.append(low)
    low *= (1+stepRatio)
    while (low < high):
        pArray.append(low)
        low *= (1+stepRatio)
    pArray.append(high)
    
    return np.array(pArray)

def writeArrayWithColons(array):
    toWrite = ""
    for ele in array:
        toWrite += str(ele) + "\t"
    # remove the last tabulation
    toWrite = toWrite[:-1]
    
    toWrite = "(" + toWrite + ")"
    
    return toWrite

class DDTransport:
    def __init__(self, mechenism, prange, Trange, accuracy = 0.05, upper = False, debug = False):
        self.mechanism = mechenism
        self.plow  = prange[0]
        self.phigh = prange[1]
        self.Tlow  = Trange[0]
        self.Thigh = Trange[1]
        self.debug = debug
        self.accuracy = accuracy
        self.upper = upper
        
        self.pArray = logspace(self.plow, self.phigh, 0.5*accuracy)
    
        # initialize gas. The gas is used to retrieve fitting coeffs
        self.gas = ct.Solution(self.mechanism)
        self.gas.TPX = 300, 101325, "O2:1"
        self.transport_model = "Mix"

    def generateTableForASpeciesNamePair(self, speciesA, speciesB, maxRefinementLevel = 5):
        baseD, TArray = self.generateTableForASpeciesNamePairAtPlow(speciesA, speciesB, maxRefinementLevel)
        
        # note D = T^1.5/p * sum(coeffs*lnT). We only have tocalculate the D at the lowest pressure
        
        # generate the D matrix [Np*NT]
        Dmatrix = np.zeros((len(self.pArray), len(baseD)))
        px, py = np.shape(Dmatrix)
        
        for i in range(px):
            Dmatrix[i,:] = baseD * self.plow/self.pArray[i]
        return Dmatrix, TArray
            
            
        
        
    def constructTArrayWithMinMaxAndStep(self, step):
        TArray = np.arange(self.Tlow, self.Thigh, step)
        if (TArray[-1] < self.Thigh):
            TArray = np.append(TArray, self.Thigh)
        return TArray
        
        
    def generateTableForASpeciesNamePairAtPlow(self, speciesA, speciesB, maxRefinementLevel = 5):
        # basic settings of interval and accuracy
        TInterval = 100
        finalAccuracy = self.accuracy
        TArray = self.constructTArrayWithMinMaxAndStep(TInterval)
        
        # generate the D array
        error = 0.5
        Darray = []
        
        refinementLevel = 0
        while((error > finalAccuracy) and (refinementLevel <= maxRefinementLevel)):
            Darray = []
            for T in TArray:
                Darray.append(self.D(self.plow, T, speciesA, speciesB))
            error = self.checkError(Darray)
            
            
            if (error > finalAccuracy):
                TInterval *= 0.5
                TArray = self.constructTArrayWithMinMaxAndStep(TInterval)
            
                
            if (self.debug):
                print("SpeciesPair = " , speciesA  + "-" + speciesB , "error = ", error, "refinementLevel = ", refinementLevel)
            
            refinementLevel += 1
        return np.array(Darray), np.array(TArray)
                
    def constructEntryForSpeciesPair(self, speciesA, speciesB, maxRefinementLevel = 5):
        # retrieve data
        Dmatrix, TArray = self.generateTableForASpeciesNamePair(speciesA, speciesB, maxRefinementLevel)
        pArray = self.pArray
        
        # generate other relevant informations
        if (self.upper):
            nameij = speciesA.upper() + "-" + speciesB.upper()+ "\n"
        else:
            nameij = speciesA + "-" + speciesB+ "\n"

        startBracket = "{" + "\n"
        endBracket   = "}" + "\n"
        typeEntry = "type   uniformTable;\n"
        lowEntry  = "low   ({}   {});\n".format(pArray[0], TArray[0])
        highEntry = "high  ({}   {});\n".format(pArray[-1], TArray[-1])
        
        valueEntry = "\n" + "values" + "\n" + str(np.size(pArray)) + "\t" +  str(np.size(TArray)) + "\n"
        tempValueEntry = ""
        for i in range(len(pArray)):
            # write data part
            tempValueEntry += (writeArrayWithColons(Dmatrix[i,:]) + "\n")
            
        #add colons
        valuEntryWithColons = "(\n" + tempValueEntry + ");\n"
            
            
        #write all
        entry = nameij + startBracket + typeEntry + lowEntry + highEntry + valueEntry + valuEntryWithColons + endBracket
        return entry
        

        
    
    
        
    
    def checkError(self, array):
        error = np.abs(np.diff(array))
        # append the last element as the last element of error to make the length of error the same as array
        error = np.append(error, error[-1])
        
        return np.max(error/array)
         


    
    def D(self, p, T, species_i, species_j):
        index_i = self.gas.species_index(species_i)
        index_j = self.gas.species_index(species_j)
        
        coeffs = self.gas.get_binary_diff_coeffs_polynomial(index_i, index_j)
        sqrtT = np.sqrt(T)
        lnT   = np.log(T) ** np.arange(0,5)
        
        value = np.sum(coeffs*lnT)*sqrtT*T/p
        
        return value
    
    def write(self, maxRefinementLevel = 5):
        # retrieve species list
        speciesList = self.gas.species_names
        
        # create the target file, if exists, delete it
        if (os.path.exists("transportProperties")):
            os.remove("transportProperties")
            
        f = open("transportProperties", "w")
        
        count = 0
        totalNumber = len(speciesList) * (len(speciesList) - 1) / 2
        
        # generate the entry for each species pair
        for i in range(len(speciesList)):
            for j in range(i, len(speciesList)):
                entry = self.constructEntryForSpeciesPair(speciesList[i], speciesList[j], maxRefinementLevel)
                
                # write the current entry to file, using the append mode
                f = open("transportProperties", "a")
                f.write(entry)
                f.close()
                
                count += 1
                print("Progress {}/{}".format(count, totalNumber), end = "\r")
        
        print("transportProperties file generated successfully")
        
        
    def writeDiffPolyNomials(self):
        # retrieve species list
        speciesList = self.gas.species_names
        
        # create the target file, it exists, delete it
        if (os.path.exists("transportPropertiesPolynomials")):
            os.remove("transportPropertiesPolynomials")
        
        f = open("transportPropertiesPolynomials", "w")
        
        count = 0
        totalNumber = len(speciesList) * (len(speciesList) - 1) / 2
        
        # print out each species pair
        for i in range(len(speciesList)):
            for j in range(i, len(speciesList)):
                nameij = ""
                if (self.upper):
                    nameij = speciesList[i].upper() + "-" + speciesList[j].upper()
                else:
                    nameij = speciesList[i] + "-" + speciesList[j]
                contents = writeArrayWithColons(self.gas.get_binary_diff_coeffs_polynomial(i, j))
                entry = "{}\t\t{{polyCoeffs\t{};}}\n".format(nameij, contents)
                
                
                # write the current entry to file, using the append mode
                f = open("transportPropertiesPolynomials", "a")
                f.write(entry)
                f.close()
                
                count +=1
                print("Progress {}/{}".format(count, totalNumber), end = "\r")
        print("transportPropertiesPolynomials file generated successfully")
                
                
        
        
