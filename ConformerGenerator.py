from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import AllChem
import sys

uLimPrune = 0.05
numConf = 10000
fileName = "AmineSmiles.txt" # assumes a text file has rows of   SMILESID    NAME  // seperated by space, not commas

name = []
smiles = []

with open(fileName,'r') as file:
    for line in file:
        lin = line.split('#')[0]
        if len(lin.split()) > 0:
            smiles.append(lin.split()[0])
            name.append(lin.split()[1])

def Conformers(mol,name, ffType = 'MMFF94s'):

    print("Generating conformers for {}".format(name))
    
    allConformers  = AllChem.EmbedMultipleConfs(mol, numConfs=numConf,pruneRmsThresh=uLimPrune, numThreads=0)
    
    print("Number of Conformers is {} out of an initial requirement of {}".format(len(allConformers),numConf))
    print("Beginning MMFF94 optimization of conformers")

    for i, conformer in enumerate(allConformers):
        
        if "mmff" in ffType.lower():
            mp = AllChem.MMFFGetMoleculeProperties(mol,mmffVariant=ffType)
            ff = AllChem.MMFFGetMoleculeForceField(mol,mp, confId=conformer)
        elif "uff" in ffType.lower():
            ff = AllChem.UFFGetMoleculeForceField( conformer, confId=conformer)
        else:
            print("Damn, no FF type has been given.")

        ff.Minimize()
        thisEnergy = ff.CalcEnergy()
        
        if i == 0:
            lowestEnergy = thisEnergy
            lowID = conformer
            index = i
            lowFF = ff
        elif thisEnergy < lowestEnergy:
            lowestEnergy = thisEnergy
            lowID = conformer
            index = i
            lowFF = ff
 

    print(index, lowFF.CalcEnergy())
    print(Chem.MolToPDBBlock(mol,confId=lowID),file=open(name + '.pdb', 'w+'))
    print(Chem.MolToMolBlock(mol,confId=lowID),file=open(name + '.mol', 'w+'))
    
####################################################
#
#        Start Script Here
#
####################################################

for i, smile in enumerate(smiles):
    mol = Chem.MolFromSmiles(smile)
    mol = Chem.AddHs(mol)
    # generate conformers and write out the lowest energy conformer to pdb/mol
    Conformers(mol,name[i])
