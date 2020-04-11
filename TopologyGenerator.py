import os
import sys
import subprocess as sub
import re
import fileinput
import numpy as np
import argparse
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('-sm',   type = str,                             help='name of special molecule')
parser.add_argument('-t',    type = str,    default = "B3LYP",          help='QM theory name, default is B3LYP')
parser.add_argument('-b',    type = str,    default = "cc-pVTZ",      help='QM basis name, default is cc-pVTZ')
parser.add_argument('-p',    type = str,    default = "esp",         help='abbreviation of population analysis call, default is esp')
parser.add_argument('-ff',   type = str,    default = "gaff2",        help='force field, default is gaff2')
parser.add_argument('-dat',  type = str,    default = 'ff99SB',      help='files used by GAFF for generating parameters, default is ff99SB')
parser.add_argument('-mult', type = int,    default = 1,             help='QM multiplicity, default is 1')
parser.add_argument('-gauss',type = str,    default = "gaussInput",  help='name gaussian input file, default is gaussInput')
parser.add_argument('-log',  type = str,    default = "gjf",         help='postfix of gaussian input file, default is gjf')
parser.add_argument('-i',    type = bool,   default = False,         help = 'Include background charge?, default is False and for water')
parser.add_argument('-cm',   type = str,                             help = 'charge method, i.e., resp, bcc, HI, CI, mbis etc...')
parser.add_argument('-g',    type = str,    default = "g16",         help = 'name of QM program i.e., g09 or g16. default is g16')
parser.add_argument('-chge', type = int,    default = 0,           help = 'charge of molecule, default is 0')
parser.add_argument('-qot',  type = str,    default = "MP2",          help = 'geometry optimization QM level of theory. Default is MP2')
parser.add_argument('-qob',  type = str,    default = "AUG-cc-pVTZ",      help = 'geometry optimization QM basis set. Default is aug-cc-pVTZ')
parser.add_argument('-res',  type = str,    default = "MOL",         help = 'default residue name of molecule, default is MOL')
parser.add_argument('-GOpt', type = bool,   default = False,         help = 'whether or not do run a geometry optimization, default is False')
parser.add_argument('-mo',   type = int,    default = 3000,          help = 'memory to use for Gaussian Optimization, default is 3000 mb')
parser.add_argument('-npp',  type = int,    default = 12,             help = 'number of processors to use for Gaussian, default is 12')
parser.add_argument('-dial', type = str,    default = "Water",        help = 'dielectric for SCRF calculation, default is water. Enter numerical value')

# for info on parsing see https://docs.python.org/3/library/argparse.html#dest

def _typeFilePresent(molName):
    files = os.listdir(os.getcwd())
    type = []
    for file in files:
        if molName in file.split(".")[0]:
            if file.endswith("mol2"):
                type.append("mol2")
            elif file.endswith("pdb"):
                type.append("pdb")
            elif file.endswith("gro"):
                type.append("gro")
            elif file.endswith("g96"):
                type.append("g96")
    if len(type) < 1:
        print('No structure file was found for molecule {}'.format(molName) )
    return type

def _obabel(fileArray):
    for ele in fileArray:
        if ele.endswith("mol2"): 
            cmd = 'module load openbabel; obabel -imol2 {0:s} -opdb {1:s}.pdb'.format(ele,ele.split('.')[0]) 
            p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0]
            print( p.decode() ) # disable once this code works
            name = ele.split('.')[0] + '.pdb'
            break
            
    return  name    

def _obabelGaussChkToPDB(name):
    #cmd1 = 'module load gaussian; formchk {0:s}.chk {1:s}.fchk; module unload gaussian'.format(name+"Opt",name)
    #cmd2 = 'module load openbabel; babel -ifchk {0:s}.fchk -opdb {1:s}.pdb'.format(name,name)
    cmd2 = 'module load gaussian; newzmat -ichk -opdb {0:s} {1:s}'.format(name+"Opt",name)
    #print("creating fchk file")
    #p1 = sub.Popen(cmd1, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0]
    #print( p1.decode() )
    print("creating pdb file")
    p2 = sub.Popen(cmd2, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0]
    print( p2.decode() )
    print("Finished creating pdb file")
        
def _gaussianInputFile(molName, QMtheory, QMbasis, pop,charge,mult):
    f = open(molName+".dat",'w')
    f.write( '%mem={}MB\n'.format(mo) )
    f.write( '%nproc={} \n'.format(npp) )
    f.write('%chk=molecule\n')
    f.write(' \n')
    if chargeMethod.lower() == "resp":
        if implicitBackground:
            f.write('#T {0}{1} pop={2} TEST density=current SCRF=({3}) SCF=tight iop(6/33=2) iop(6/42=6) iop(6/50=1)\n'.format(QMtheory, QMbasis,pop, dialKeyword))
        else:
            f.write('#T {0}{1} pop={2} TEST density=current SCF=tight iop(6/33=2) iop(6/42=6) iop(6/50=1)\n'.format(QMtheory, QMbasis,pop))
    elif QMtheory == "AM1" or QMtheory == "PM6":
        f.write('#T {0} \n'.format(QMtheory) )
    else:
        if implicitBackground:
            f.write('#T {0}{1} pop={2} density=current SCRF=({3})\n'.format(QMtheory, QMbasis,pop, dialKeyword) )
        else:
            f.write('#T {0}{1} pop={2} density=current\n'.format(QMtheory, QMbasis,pop) )

    # Density=Current Opt # other options
    f.write(' \n')
    f.write('Remark Section')
    f.write(' \n \n'.format())
    f.write('{0} {1} ! Molecule specification\n'.format(charge,mult))

    atype = []

    getLength = open(molName + ".pdb")
    length = 0
    for line in getLength:
        if len(line.split()) > 0:
            if ("HETATM" in line) or ("ATOM" in line.split()[0] ):
                length += 1
    getLength.close()

    coords = np.zeros(( length,3))       
    coords,atype = _getPDBMolCoords(molName,length)

    for i in range(length):
        f.write('{0} {1:7.7} {2:7.7} {3:7.7} \n'.format(atype[i], coords[i][0], coords[i][1], coords[i][2]))
    f.write(' \n')
    if chargeMethod == "resp":
        if implicitBackground and input.dial.lower() != "water":
            f.write( 'eps={} \n\n'.format(input.dial) )
        f.write(molName + ".gesp")
        f.write(' \n \n'.format())
        f.write(molName + ".gesp")
        f.write(' \n \n'.format())
    elif implicitBackground and input.dial.lower() != "water":
        #f.write( 'Alpha=0.5 \n' )
        f.write( 'eps={} \n\n'.format(input.dial) )
        #f.write( 'epsinf=1.77 \n\n' )
    f.close()

def _gaussianOpt(molName, QMOptTheory, QMOptBasis,charge,mult):
    f = open(molName + "Opt" +".dat",'w')
    f.write( '%mem={}MB\n'.format(mo) )
    f.write( '%nproc={} \n'.format(npp) )
    f.write('%chk={0:s}Opt\n'.format(molName) )
    f.write(' \n')
    if implicitBackground:
        f.write('#T {0}/{1} Opt Freq Density=Current SCRF=(Solvent=Water)\n'.format(QMOptTheory, QMOptBasis)) # use implicit solvent
    else:
        f.write('#T {0}{1} Opt Freq Density=Current \n'.format(QMOptTheory, QMOptBasis,pop))
    # Density=Current Opt # other options
    f.write(' \n')
    f.write('Remark Section')
    f.write(' \n \n'.format())
    f.write('{0} {1} ! Molecule specification\n'.format(charge,mult))

    atype = []

    getLength = open(molName + ".pdb")
    length = 0
    for line in getLength:
        if len(line.split()) > 0:
            if ("HETATM" in line.rstrip().split()[0] ) or ("ATOM" in line.rstrip().split()[0] ):
                length += 1

                    
    getLength.close()

    coords = np.zeros(( length,3))       
    coords,atype = _getPDBMolCoords(molName,length)

    for i in range(length):
        f.write('{0} {1:7.7} {2:7.7} {3:7.7} \n'.format(atype[i], coords[i][0], coords[i][1], coords[i][2]))
    f.write(' \n')

    f.close()


def _getPDBMolCoords(molName,length):
    pdb = open(molName+".pdb") 
    coords = np.zeros(( length,3))
    atype = [] 
    iter = 0
    for line in pdb:
        if len(line.split()) > 0:
            if ("HETATM" in line.rstrip().split()[0] ) or ("ATOM" in line.rstrip().split()[0] ): 
                if "mol" in line.lower():  # some pdb's skip the "mol" column, and instead have atom types at end. 
                    coords[iter] =  line.split()[5:8] 
                else:
                    coords[iter] =  line.split()[4:7] 
                atype.append( line.split()[2][0].upper() ) # take first letter of string in column 3, and capitalize it
                iter += 1

    pdb.close()
    return coords, atype 

def _getGaussianCharges(gauss_output):

    filename = open(gauss_output)            # this is the log file output by Gaussian
    loop = True
    count = 0
    while loop and count < 10: # sometimes ESP is output more than once, we need the last one

        count += 1
        if pop == 'ESP':
            for line in filename:
                if 'ESP charges:' in line: 
                    gout = []
                    break
            
            for line in filename:  # This keeps reading the file
                if 'Sum of ESP charges' in line: 
                    break
                if ("" == line):
                    loop = False # reached EOF, exit

                gout.append(line)
    
        elif 'MBS' in pop:
            for line in filename:
                if 'Mulliken charges:' in line:
                    gout = []
                    break

            for line in filename:
                if 'Sum of Mulliken Charges' in line:
                    break
                if ("" == line):
                    loop = False # reached EOF, exit

                gout.append(line)

        elif 'CM5' == chargeMethod or 'Hirshfeld' == chargeMethod:
            for line in filename:
                if ('Hirshfeld charges' in line) or ('CM5 charges' in line):
                    gout = []
                    break

            for line in filename:
                if 'Tot' in line:
                    break
                if ("" == line):
                    loop = False # reached EOF, exit

                if "Q-H" not in line:
                    gout.append(line)
        elif 'hi' in chargeMethod.lower() or 'ci' in chargeMethod.lower():
            for line in filename:
                if ('Iterated Hirshfeld charges, spin densities' in line) or ('CM5 charges' in line):
                    gout = []
                    break

            for line in filename:
                if 'Tot' in line:
                    break
                if ("" == line):
                    loop = False # reached EOF, exit

                if "Q-H" not in line:
                    gout.append(line) 
    goutsplit=[]
    print(gout)
    for  x in gout:              # split stats into chunks
        y = x.split()
        if len(y) > 2:      # first few lines read in are not needed and are only single elements i.e. [i]
            goutsplit.append(y)
    newCharges = np.zeros(len(gout))
    for i in range( len(gout) ):
        if chargeMethod.lower() == "ci": #take the Iterated CM5 charge version
            newCharges[i] = float(goutsplit[i][7])
        else:
            newCharges[i] = float(goutsplit[i][2])
            

    return newCharges

def _updateMol2Charges(filename,newCharges):

    file = open(filename,'r')
    fileout = open("tempMol2.txt",'w')
    for line in file:
        fileout.write(line)
        if "@<TRIPOS>ATOM" in line:
            break
    for i, line in enumerate(file):
        if "@<TRIPOS>BOND" in line:
            fileout.write(line)
            break
        newline = line[0:72] + str(newCharges[i])
        fileout.write(newline + '\n')
    for line in file:
        fileout.write(line)

    file.close()
    fileout.close()
    """ Now remove original mol2, and rename temp to that name"""
    os.remove(filename)
    os.rename("tempMol2.txt",filename)
 
"Generates the Horton input file so that Horton can be called for partial charges. " 
def GenerateHortonFile(name):

    f = open(name + ".py",'w')

    f.write( '#!/usr/bin/env python \n\n' )
    f.write( 'import numpy as np \n\n' )
    f.write( 'from horton import * \n\n' )
    f.write( "fn_fchk = 'molecule.fchk' \n" )
    f.write( 'mol = IOData.from_file(fn_fchk) \n\n')
    f.write( '# Partition the density with the Becke scheme \n')
    f.write( "grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, mode='only', agspec='ultrafine') \n" )
    f.write( 'moldens = mol.obasis.compute_grid_density_dm(mol.get_dm_full(), grid.points) \n' )
    f.write( 'wpart = MBISWPart(mol.coordinates, mol.numbers, mol.pseudo_numbers, grid, moldens) \n')
    f.write( 'wpart.do_charges() \n')
    f.write( '# Write the result to a file \n')
    f.write( "np.savetxt('hortonCharges.txt', wpart['charges']) \n")

    f.close()
def GetDDEC6ChargesFromFile():

    f = open("DDEC6_even_tempered_net_atomic_charges.xyz",'r')
    charges = []
    for line in f:
        if len( line.strip().split()) > 1:
            break
    iter = -1
    for line in f:
        iter += 1
        if iter == ( len( atype )  ):
            break
        charges.append( float( line.split()[4] ) )
    f.close()
    return charges

def GenerateDDEC6File():

    f = open("job_control.txt",'w')
    f.write( '<atomic densities directory complete path> \n' )
    f.write( '/home/bkelly08/Programs/chargemol_09_26_2017/atomic_densities/ \n' )
    f.write( '</atomic densities directory complete path> \n\n')
    f.write( '<input filename>\n' )
    f.write( '{0}.wfx \n'.format(molName) )
    f.write( '</input filename> \n\n' )
    f.write( '<charge type> \nDDEC6 \n</charge type> \n\n' )
    f.write( '<compute BOs> \n.false. \n</compute BOs>' )
    
    f.close()
    
def GenerateWFXFile():
    # Gaussian will generate the wfx file when it is called
    filename = "coords"
    
    GenerateZMatrixFile(filename)
 
    cmd = 'module load gaussian; g16 zMatrixInput.dat; module unload gaussian'
    p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0] 
    
def GenerateZMatrixFile(filename):
  
    PrintXYZ(atype, coords, filename + ".xyz")    # creates xyz file with atom names and coordinates
    ConvertXYZ2ZMat(filename)                     # creates Zmat file that is by default ".com"
    
    zMat = CopyZMatrix(filename + ".com")

    f = open("zMatrixInput.dat",'w')
    f.write( ' %chk={0}.chk \n '.format(molName) )
    f.write( '%mem={}MB\n'.format(mo) )
    f.write( '%nproc={} \n'.format(npp) )
    f.write( '#P {0}/{1} geom=connectivity guess=mix CHARGE DENSITY=CURRENT \n'.format(QMtheory,QMbasis) )    # PW91PW91
    f.write( '# scf=(fermi,conver=8,maxcycle=400) density=current output=wfx \n\n' )
    f.write( 'Truly witty remark \n\n' )
    f.write( '{0} {1} \n'.format(charge, mult) )
    for line in zMat:
        f.write( '{0} \n'.format( str(line) ) )
        
    if background == True: # This is taken from QMMM_Interface, and so there shouldn't be solCoords...
        for i, item in enumerate(solCoords) :
            f.write('{0:7.7} {1:7.7} {2:7.7} {3:7.7} \n'.format(solCoords[i][0], solCoords[i][1], solCoords[i][2], chargeGroup[i] ) )
    f.write(' \n')
    f.write( '{0}.wfx \n\n'.format( molName ) )
    f.close()

def ConvertXYZ2ZMat(filename):
    cmd = 'module load gaussian; newzmat -ixyz {0:s}; module unload gaussian'.format(filename)
    p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0]
    
def CopyZMatrix(file):
    f = open(file,'r')

    zMatrixStuff = []
    
    for line in f:
        ln = line.split(",")
        if len(ln) > 1 :#ln[0] == sumQQ and ln[1] == mult:
            break
    for line in f:
        zMatrixStuff.append(line.strip())
        
    f.close()
    return zMatrixStuff
    
def PrintXYZ(atomNames, coords, filename):

    f = open(filename,'w')

    for i in range(len(atomNames)):
        f.write('{0} {1:7.7} {2:7.7} {3:7.7} \n'.format(atomNames[i][:1], coords[i][0], coords[i][1], coords[i][2]))
    f.write(' \n')
    
    f.close()

##############################################################################
#
#
#
#                               START THE CODE
#
#
##############################################################################      
"""
Need to setup some user defined variables
   - name of molecule
   - charge on molecule
   
"""

input = parser.parse_args()

molName      = input.sm    # solute molecule
pop          = input.p     # population analysis method i.e., ESP or MBS (mulliken)
QMtheory     = input.t     # "HF"  #HF/6-31G*  # B3LYP # HF
if input.t == "AM1" or input.t == "PM6":
    QMbasis = ""
else:
    QMbasis      = "/" + input.b     # "6-31G(d)" #"AUG-cc-pVTZ", "6-31G*"
QMprogram    = input.g     # g09/g16
mult         = input.mult  # multiplicity of molecule
ffType       = input.ff    # gaff
#goutput      = input.o
inputname    = input.gauss
end          = input.log
chargeMethod = input.cm    # Needed to distinguish Iterative Hirshfeld, Iterative CM5, ddec6, MBIS etc...

charge       = input.chge
QMOptTheory  = input.qot
QMOptBasis   = "/" + input.qob
residueName  = input.res
GaussOpt     = input.GOpt  # do a geometry optimization first?
mo           = input.mo    # amount of memory to use in Gaussian
npp          = input.npp   # number of processors to use in Gaussian
ffDat        = input.dat
implicitBackground = input.i   # use implicit SCRF background solvent in QM

if input.dial.lower() == "water":
    dialKeyword = "Solvent=Water"
else:
    dialKeyword = "Solvent=Water,read"

print("Geometry Optimization? ", GaussOpt)

background = False
if chargeMethod.lower() == "hi" or chargeMethod.lower() == "ci":
    pop = "iterhirsh"
####################################################################

###### Check if geometry optimization is required and do if so #####

if GaussOpt:
    """ need to optimize structure """
    print("Request to optimize structure acknowledged")
    OptName = molName + "Opt"

    print("Creating input for optimizing structure in gaussian")
    _gaussianOpt(molName, QMOptTheory, QMOptBasis,charge,mult)

    print("Running gaussian to optimize structure")
    cmd = 'module load gaussian; g16 {0:s}.dat; module unload gaussian'.format(OptName)
    p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0]
    print( p.decode() )
    
    print("Converting gaussian checkpoint to pdb file")
    shutil.copyfile(molName+".pdb",molName+"_original.pdb")
    _obabelGaussChkToPDB(molName)  # create pdb from gaussian checkpoint file
    #molName = OptName
    shutil.copyfile(molName+"Opt.log","optimizedLogFile.log")

    cmd = 'python ~/scripts/python/modify_pdb.py -sm {}'.format(molName)  # this modifies the pdb to a nicer format... atoms get numbered rather than a column of C C H C O C N, you have C1 C2 H1 C3 O1 C4 N1
    p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0]
    print( p.decode() )
####################################################################
coordFileType = []
    
coordFileType = _typeFilePresent(molName)
if "pdb" in coordFileType:
    """ skip ahead right to antechamber, but get that name"""
    for ele in coordFileType:
        if "pdb" in ele:
            coordFile = ele 

else: 
    """ Convert file to pdb """
    coordFile = _obabel(coordFileType) 
print("coordFile is: ",coordFile)	
coordFile = molName + "." + coordFile        
### We would like our molecule to be called MOL rather than LIG or UNL etc """

replaceWords = ["LIG","RES","UNL","UNK"]
for i in range( len( replaceWords ) ):
    with fileinput.FileInput(coordFile, inplace=True) as file:
        for line in file:
            print(line.replace(replaceWords[i], residueName ), end='')  # residueName is MOL by default

# get elements for each atom
atype = []

getLength = open(molName + ".pdb")
length = 0
for line in getLength:
    if len(line.split()) > 0:
        if ("HETATM" in line) or ("ATOM" in line.split()[0] ):
            length += 1
getLength.close()

coords = np.zeros(( length,3))       
coords,atype = _getPDBMolCoords(molName,length)   
########################################
#
#    Proceed to Antechamber or Horton or Other
#
########################################    

""" create gaussian input file"""
if chargeMethod.lower() == "resp":
    print("Step 1: Create Gaussian input file for generating gesp file")
    _gaussianInputFile(molName, QMtheory, QMbasis, pop,charge,mult)

    """ Now run gaussian to get .gesp file """
    print("Step 1b: Running Gaussian...")
    cmd = 'module load gaussian; g16 {0:s}.dat; module unload gaussian'.format( molName )
    p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0] 
    print(p.decode())
        
    """ Now run Antechamber again to get resp charges (they are in mol2 output file) """   
    print("Step 2: Running Antechamber to get resp charges...")    
    cmd = 'module load gcc/5.4.0; module load amber/18; antechamber -i {0:s}.gesp -fi gesp -o {1:s}.mol2 -fo mol2 -c {2:s} -at {3:s} -nc {4} -pf y'.format(molName,molName, chargeMethod, ffType, charge)
    p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0] 
    print(p.decode())    
elif ("bcc") in chargeMethod.lower():
    print("Running Antechamber to get AM1BCC charges...")
    cmd = 'module load gcc/5.4.0; module load amber/18;antechamber -i {0:s}.pdb -fi pdb -o {1:s}.mol2 -fo mol2 -c {2:s} -at {3:s} -nc {4} -pf y '.format(molName,molName, chargeMethod, ffType, charge)
    p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0]
    print( p.decode() )
elif chargeMethod.lower() == "hi" or chargeMethod.lower() == "ci":

    """ this is a workaround... first calculate mol2 with bcc charges, then call gaussian, get Actual charges and replace mol2 with these charges """
    print("Running hirshfeld charge calculation...")
    cmd = 'module load gcc/5.4.0; module load amber/18;antechamber -i {0:s}.pdb -fi pdb -o {1:s}.mol2 -fo mol2 -c bcc -at {3:s} -nc {4} -pf y'.format(molName,molName, chargeMethod, ffType, charge)
    p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0]
    print( p.decode() )

    _gaussianInputFile(molName, QMtheory, QMbasis, pop,charge,mult)

    """ Now run gaussian to get charge file """
    cmd = 'module load gaussian; g16 {0:s}.dat; module unload gaussian'.format( molName )
    p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0] 
    print(p.decode())
    newCharges = _getGaussianCharges(molName +".log")     # gets the new charges from the gaussian log file

    _updateMol2Charges(molName + ".mol2",newCharges)     # updates the charges in the mol2 file with the new charges

elif chargeMethod.lower() == "mbis":
    # Minimum Basis Iterative Stockholder charges
    # step 1) call antechamber for a dummy force-field with dummy bcc charges
    # step 2) call gaussian to generate a chk file.
    # step 3) convert chk file to fchk file
    # step 4) call Horton, give fchk file, get back Horton charges (likely MBIS)
    # step 5) modify mol2 file created by Antechamber, put in the Horton charges
    # step 6) done later, but proceed with Antechamber FF generation, it doesn't know we changed the charges ;)
    print( "Minimum Basis Iterative Stockholder charges" )
    ### this is a workaround... first calculate mol2 with bcc charges, then call gaussian, get Actual charges and replace mol2 with these charges 
    
    # Step 1)
    print( "Step 1) call antechamber for dummy mol2 file" )
    cmd = 'module load gcc/5.4.0; module load amber/18;antechamber -i {0:s}.pdb -fi pdb -o {1:s}.mol2 -fo mol2 -c bcc -at {3:s} -nc {4} -pf y'.format(molName,molName, chargeMethod, ffType,charge)
    p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0]
    print( p.decode() )
    
    # Step 2)
    print( "Step 2) make Gaussian input File...")
    _gaussianInputFile(molName, QMtheory, QMbasis, pop,charge,mult)

    ### Now run gaussian to get charge file """
    print( "Step 2b) Running Gaussian...")
    cmd = 'module load gaussian; g16 {0:s}.dat; module unload gaussian'.format( molName )
    p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0] 
    #print(p.decode())
    
    # Step 3)
    print( "Step 3) Converting Gaussian checkpoint file to Formatted checkpoint...")
    # convert the (Gaussian?) checkpoint to a formatted checkpoint file
    cmd = 'module load gaussian; formchk {0:s}.chk {1:s}.fchk; module unload gaussian'.format("molecule","molecule")
    p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0]
    #print( p.decode() )

    # Step 4)
    print( "Step 4) Generating Horton Inputfile.")
    GenerateHortonFile("mbis") # generates the... Horton ... file... it is kind of obvious..., takes as input the name of the python script to be produced
    # Now run Horton
    print( "Step 4b) Running Horton...")
    cmd = 'source ~/ENV/bin/activate; export PYTHONPATH=~/.local/lib/python2.7/site-packages/:$PYTHONPATH;python mbis.py'
    p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0]  
    print( p.decode() )
   
    # A text file should now be present in the folder with the charges. It should also be called hortonCharges.txt
    file = open("hortonCharges.txt",'r')
    newCharges = []
    for iter, charge in enumerate(file):
        newCharges.append( format(float(charge.split()[0] ), ' 20.20f' ) ) # can customize with , ' 10.8f') etc...
    file.close()
    # Step 5)
    print( "Updating mol2 file with MBIS charges.")
    _updateMol2Charges(molName + ".mol2",newCharges)     # updates the charges in the mol2 file with the new charges
    
elif chargeMethod.lower() == "ddec6":
    ### Same thing here, need to call Antechamber, make a mol2 file, then call Chargemol - get ddec6 charges, and then replace the mol2 charges with ddec6 charges
    
    # Step 1)
    print( "Step 1) call antechamber for dummy mol2 file" )
    cmd = 'module load gcc/5.4.0; module load amber/18;antechamber -i {0:s}.pdb -fi pdb -o {1:s}.mol2 -fo mol2 -c bcc -at {3:s} -nc {4:d} -pf y'.format(molName,molName, chargeMethod, ffType, charge)
    p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0]
    print( p.decode() )
    
    print("Step 2a: Calling Gaussian for WFX file")
    GenerateWFXFile()
    GenerateDDEC6File()
    print("Step 2b: Calling Chargemol for DDEC6 charges")
    cmd = 'module load gcc/6.4.0 openmpi/2.1.1; chmod u+x Chargemol_09_02_2017_linux_parallel; ./Chargemol_09_02_2017_linux_parallel'
    p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0] 
    print( p.decode() )
    GetDDEC6ChargesFromFile()
    print("Finished calling DDEC6 Charges")
    
    print("Step 3: putting ddec6 charges in mol2 file")
    file = open("DDEC6_even_tempered_net_atomic_charges.xyz",'r')
    newCharges = []
    for line in file:
        if "Nonperiodic system" in line:
            break
    
    for line in file:
        if "Chargemol" in line:
            break
        if len(line.split()) > 3:
            newCharges.append( float(line.split()[4] ) ) # can customize with , ' 10.8f') etc...
    file.close()
    print("newcharges: ", newCharges)
    # Step 4)
    print( "Step 4: Updating mol2 file with ddec6 charges.")
    _updateMol2Charges(molName + ".mol2",newCharges)     # updates the charges in the mol2 file with the new charges
    
else:
    print("No chargeMethod listed. Please select either 'resp' or 'bcc' or some new fangled method of your own desire")

### Now use parmchk """       
cmd = 'module load gcc/5.4.0; module load amber/18;parmchk2 -i {0:s}.mol2 -f mol2 -o {1:s}.frcmod'.format(molName,molName)
p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0] 
print(p.decode())      
        
########################################
#
#   Time to use leap
#
########################################        

### edit tleap.in i.ee, replace all NAMES with actual molecule name"""
f = open("tleap.in",'w')
#f.write('source oldff/leaprc.{0:s} \n'.format(ffDat) )
f.write('source leaprc.{0:s} \n'.format(ffType) )
f.write('loadamberparams {0:s}.frcmod\n'.format(molName) )
f.write('MOL = loadmol2 {0:s}.mol2 \n'.format(molName) )
f.write('check MOL \n')
f.write('saveoff MOL {0:s}.lib \n'.format(molName) )
f.write('saveamberparm MOL {0:s}.prmtop {1:s}.inpcrd \n'.format(molName,molName) )
f.write('quit\n')
f.close()
#cmd = 'module load gcc/5.4.0; module load amber/18; tleap -s -f tleap.in'
cmd = 'module load gcc/5.4.0; module load amber/18;tleap -s -f tleap.in '
p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0] 
print(p.decode())


########################################
#
#   Time to use acpype to get GMX.top
#
########################################
cmd = 'module load nixpkgs/16.09 intel/2016.4 openmpi/2.1.1 python/2.7.14 scipy-stack/2019a gaussian; python ~/acpype.tgz/acpype/acpype.py -p {0:s}.prmtop -x {1:s}.inpcrd'.format(molName,molName)
p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0]   
print(p.decode())  

########################################
#
#   Change name of MOL files
#
########################################
if "bcc" in QMtheory.lower() or "pm" in QMtheory.lower():  # if semiempirical, don't write out basis... there isn't one...
    cmd = 'cp MOL_GMX.top {}_{}_{}_{}_GMX.top'.format(molName, chargeMethod.lower(), QMtheory,ffType )       
    p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0]
else:
    cmd = 'cp MOL_GMX.top {}_{}_{}_{}_{}_{}_{}_GMX.top'.format(molName, chargeMethod.lower()
    , QMtheory, str(input.b),ffType,implicitBackground,input.dial )       
    p = sub.Popen(cmd, shell=True, stderr = sub.STDOUT, stdout = sub.PIPE).communicate()[0]

