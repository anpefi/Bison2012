#! /usr/bin/env python
'''
GAS (Genetic Age-structured Simulations)  v.0.4.3 (09/10/09)
Coded by Andres Perez-Figueroa and Tiago Antao
Fitted simuPOP 0.9.3
'''


#Classes
class OptionError(Exception):
    """
    Class for exception handling in options reading
    """
    def __init__(self, param, message):
        self.param=param
        self.message=message

#IMPORTS.  
from sys import exit
import os
import bz2
try:
    from simuOpt import *
    setOptions(quiet=True)
    from simuPOP import *  #Check version??? Only important if source code will be distributed
    from simuUtil import *
except:
    print "simuPOP (or one of its components) import failed. Please check your simuPOP installation."
    exit(2)
try:
    from math import log, exp, sqrt
    from random import *
    from statlib import stats
except:
    print "Error importing math/random/statlib python modules. Please check them"
    exit(2)
try:
    import numpy
except:
    print "numpy import failed. Python's numpy package is required to run this program. Please check if numpy is installed for your python"
    exit(2)
try:
    import matplotlib
    matplotlib.use('Agg') #Avoids crash when running on systems without a X display
    from pylab import *
except:
    usePylab = False
    print "matplotlib import failed. Output graphics won't be generated"
else:
    usePylab = True


    
#Globals!
obsN0 = 0
obsLambda = 1
propAlphas = 0.0
Vkmale = []
fatherslist = [] #experimental

#Options dictionary. It will allow to set options in command line and get the parameters from a file
options = [
    {'longarg':'name=',
     'default':'Bison',
     'allowedTypes': [types.StringType],
     'label':'Name of the simulation',
     'description': 'Base name for configuration (.cfg) and prefix for resulting files and figures',
    },
    {'longarg':'reps=',
     'default':2,
     'label':'Replicates',
     'description':'Number of replicates to run the simulation. max: 10000',
     'allowedTypes': [types.IntType],
     'validate': valueBetween(1,10000),
    },
    {'longarg':'time=',
     'default':10,
     'label':'Years',
     'description':'Time units (years) of evolution. max: 1000',
     'allowedTypes': [types.IntType],
     'validate': valueBetween(1,1000),
    },
    {'longarg':'burnin=',
     'default':10,
     'label':'Burn-in',
     'description':'Lenght of burn-in stage',
     'allowedTypes': [types.IntType],
     'validate': valueBetween(0,999),
    },
    {'longarg':'nMSats=',
     'default':10,
     'label':'number of Microsats',
     'description':'Number of microsatellite-type loci (multi-allelic). max:100',
     'allowedTypes': [types.IntType],
     'validate': valueBetween(0,100),
    },
    {'longarg':'maxAllele=',
     'default':20,
     'label':'maximum of alleles',
     'description':'Maximum number of alleles in microsatellite loci. max:99',
     'allowedTypes': [types.IntType],
     'validate': valueBetween(2,99),
    },
    {'longarg':'initFreq=',
     'default':'DIRICHLET',
     'label':'Frequencies initalization',
     'description':'Choose model for intialization of allelic frequencies in microsat loci.',
     'chooseOneOf':['EQUAL', 'DIRICHLET','FILE'],
     'validate': valueOneOf(['EQUAL', 'DIRICHLET','FILE']),
    },
    {'longarg':'fileFreq=',
     'default':'',
     'label':'Frequencies datafile',
     'description':'File name for user-defined frequencies. Used only if initFreq = FILE. ',
     'allowedTypes': [types.StringType],
     'validate': valueOr(valueValidFile(),valueOneOf([''])),
     # TODO: filename input needs validation in order to avoid exceptions
    },
    {'longarg':'nSNPs=',
     'default':0,
     'label':'number of SNPs',
     'description':'Number of SNP-type loci (bi-allelic). max:100',
     'allowedTypes': [types.IntType],
     'validate': valueBetween(0,100),
    },
    {'longarg':'initialN=',
     'default':430,
     'label':'initial N',
     'description':'Initial population size. max:10000',
     'allowedTypes': [types.IntType],
     'validate': valueBetween(1,100000),
    },
    {'longarg':'constSize',
     'default':True,
     'label':'Constant size',
     'allowedTypes': [types.BooleanType],
     'description':'Force population with constant size over time',
    },
    {'longarg':'cohortN=',
     'default':100,
     'label':'Cohort size N0',
     'description':'Cohort size. Number of newborns generated if constSize is selected and during burnin',
     'allowedTypes': [types.IntType],
     'validate': valueBetween(1,100000),
    },
    {'longarg':'saveGens=',
     'default':'0,1,2,3,6',
     'label':'Years to save',
     'description':'List of years to save population in a genepop-format file. Use with caution, very Huge files!',
     'allowedTypes': [types.StringType],
     #'validate': valueIsList(),
    },
    {'longarg':'maxAge=',
     'default':20,
     'label':'Max age',
     'description':'maximun age of individuals. Note that 0 is the age of newborns (max: 50)',
     'allowedTypes': [types.IntType],
     'validate': valueBetween(0,50),
    },
    {'longarg':'adultAge=',
     'default':2,
     'label':'Adult age',
     'description':'age of individuals when they become adults',
     'allowedTypes': [types.IntType],
     'validate': valueBetween(0,50),
    },
    {'longarg':'ASinit=',
     'default':'25.445,13.679,10.162,8.839,8.208,7.142,6.003,5.339,4.543,3.374,2.712,2.042,1.058,0.600,0.322,0.201,0.112,0.082,0.060,0.044,0.033,0.000',
     'label':'Initial Age Structure',
     'description':'List of probability of age of individuals. They are relative, there is no need to sum up to 1. Separate number by commas. Should be as many as ages',
     'allowedTypes': [types.StringType],
     #'validate': valueIsList(),
    },
    {'longarg':'matingMode=',
     'default':'MONOGFEMS1',
     'label':'Mating mode',
     'description':'Mode of mating: RANDOM: random mating with replacement; MONOGFEMS1: monogamic females with one offspring per female',
     'chooseOneOf':['MONOGFEMS1', 'RANDOM'],
     'validate': valueOneOf(['MONOGFEMS1', 'RANDOM']),
    },
    {'longarg':'sMales=',
     #'default':'0.580,0.640,0.819,0.829,0.835,0.837,0.836,0.831,0.822,0.809,0.791,0.767,0.736,0.698,0.650,0.593,0.527,0.453,0.375,0.298,0.226',
     'default':'0.550,0.760,0.89,0.95,0.89,0.86,0.91,0.87,0.76,0.82,0.77,0.54,0.57,0.55,0.64,0.57,0.75,0.75,0.75,0.75,0.75',
     'label':'survivorship (s) on males',
     'description':'List of age-specific survivorship in males. Separate number by commas. Should be as many as ages',
     'allowedTypes': [types.StringType],
     #'validate': valueIsList(),
    },
    {'longarg':'sFemales=',
     #'default':'0.580,0.640,0.819,0.829,0.835,0.837,0.836,0.831,0.822,0.809,0.791,0.767,0.736,0.698,0.650,0.593,0.527,0.453,0.375,0.298,0.226',
     'default':'0.550,0.760,0.89,0.95,0.89,0.86,0.91,0.87,0.76,0.82,0.77,0.54,0.57,0.55,0.64,0.57,0.75,0.75,0.75,0.75,0.75',
     'label':'survivorship (s) on females',
     'description':'List of age-specific survivorship in females. Separate number by commas. Should be as many as ages',
     'allowedTypes': [types.StringType],
     #'validate': valueIsList(),
    },
     {'longarg':'bMales=',
     'default':'0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1',
     'label':'Age-dependent probability of mate males',
     'description':'List of age-specific probability of mates in males. They are relative, there is no need to sum up to 1. Separate number by commas. Should be as many as ages',
     'allowedTypes': [types.StringType],
     #'validate': valueIsList(),
    },
    {'longarg':'bFemales=',
     #'default':'0,0,0.275279324,0.398912121,0.51899086,0.620106432,0.696777653,0.750634658,0.785834983,0.806277304,0.814572581,0.811838379,0.797703091,0.770298949,0.726313808,0.661503159,0.572485695,0.460581969,0.336261303,0.218573287,0.125647857',
     'default':'0,0,0.7,0.7,0.95,0.93,0.93,0.93,0.93,0.92,0.92,0.92,0.92,0.92,0.8,0.8,0.572485695,0.460581969,0.336261303,0.218573287,0.125647857',
     'label':'Fecundity (b) on females',
     'description':'List of age-specific fecundity in females. If the population is growing they mean the average number of offspring per indivdual or each age. Separate number by commas. Should be as many as ages',
     'allowedTypes': [types.StringType],
     #'validate': valueIsList(),
    },
    {'longarg':'culling',
     'default':False,
     'label':'Culling',
     'allowedTypes': [types.BooleanType],
     'description':'Allow culling',
    },
    {'longarg':'cullObjective=',
     'default':0,
     'label':'Culling objective',
     'description':'Population size after culling',
     'allowedTypes': [types.IntType],
     'validate': valueBetween(0,10000),
    },
    {'longarg':'cullTarget=',
     'default':'ALL',
     'label':'Culling target',
     'description':'target for culling: ALL: random culling; YOUNG: Random culling of non-adult individuals; ADULTS: random culling of adult individuals',
     'chooseOneOf':['ALL', 'YOUNG','ADULTS','ADULTS10'],
     'validate': valueOneOf(['ALL', 'YOUNG','ADULTS','ADULTS10']),
    },
    {'longarg':'cullMode=',
     'default':'INTERVAL',
     'label':'Culling mode',
     'description':'Mode for culling: INTERVAL: culling occurs every x years; THRESHOLD: culling occurs when population size is larger than a threshold',
     'chooseOneOf':['INTERVAL', 'THRESHOLD'],
     'validate': valueOneOf(['INTERVAL', 'THRESHOLD']),
    },
    {'longarg':'cullLevel=',
     'default':0,
     'label':'Culling threshold',
     'description':'Population size threshold to do culling if mode = THRESHOLD',
     'allowedTypes': [types.IntType],
     'validate': valueBetween(0,10000),
    },
    {'longarg':'cullInterval=',
     'default':0,
     'label':'Culling Interval',
     'description':'Interval of time between culling events if mode = INTERVAL',
     'allowedTypes': [types.IntType],
     'validate': valueBetween(0,10000),
    },
    {'longarg':'alphaMating',
     'default':False,
     'label':'Alpha males',
     'allowedTypes': [types.BooleanType],
     'description':'Enable the alpha males mating system',
    },
    {'longarg':'alphaOff=',
     'default':0.0,
     'label':'Prop. offspring by alphas',
     'description':'Proportion of offspring sired by alphas',
     'allowedTypes': [types.FloatType],
     'validate': valueBetween(0.0,1.0),
    },
    {'longarg':'alphaMinAge=',
     'default':0,
     'label':'Minimum age of alphas',
     'description':'Minimum age of an individual who can become an alpha',
     'allowedTypes': [types.IntType],
     'validate': valueBetween(0,50),
    },
    {'longarg':'alphaMaxAge=',
     'default':0,
     'label':'Maximum age of alphas',
     'description':'Maximum age of an individual who can become an alpha',
     'allowedTypes': [types.IntType],
     'validate': valueBetween(0,50),
    },
    {'longarg':'alphaDuration=',
     'default':0,
     'label':'duration of alpha',
     'description':'Duration of alpha trait',
     'allowedTypes': [types.IntType],
     'validate': valueBetween(0,50),
    },
    {'longarg':'alphaProp=',
     'default':0.0,
     'label':'Prop. of alphas',
     'description':'Proportion of males being alphas',
     'allowedTypes': [types.FloatType],
     'validate': valueBetween(0.0,1.0),
    },
    {'longarg':'saveDecays',
     'default':False,
     'allowedTypes': [types.BooleanType],
     'description':'Save AD, He and Ho in a file for posterior analysis of distribution',
    },
    {'longarg':'cullAge=',
     'default':3,
     'description':'Threshold age for culling adults (ADULTS10)',
     'allowedTypes': [types.IntType],
     'validate': valueBetween(0,30),
    }
    
]



#**********************************************************************************************************************#
#Functions:
#**********************************************************************************************************************#
def getOptions(details=__doc__):
    pars = simuOpt(options,__doc__)
    pars.getParam(nCol=3)
    
    #Automatically save configurations
    name = pars.name
    if not os.path.isdir(name):
        os.makedirs(name)
    pars.saveConfig( os.path.join(name, name+'.cfg'))
    #return the rest of the parameters
    return True

def createSinglePop(popSize, genomeLoci):
    '''Set the options for create a simuPOP's single population '''
    preOps = [pyExec('cprog=repli*gens')]
    inOps = []
    pop = population(size=popSize, ploidy=2,
        loci=genomeLoci,
        infoFields=['age', 'fitness', 'alphaField', 'father_idx', 'mother_idx','father_ref', 'mother_ref','ref','ageFather', 'ageMother','cohort'],
        #maxAllele = maxAlleleN  maxAllele is not working on simuPOP 0.9
    )
    oExpr = '"%s/sims/%%d/pop-%d.gen" %% (gen-burnin)' %  (filePrefix,repli)
 
    return pop, preOps, inOps, oExpr


def createGenome(numMSats, numSNPs):
    '''
    Initialize the genome
    '''
    loci = (numMSats+numSNPs)*[1]
    preOps = []
    
    if numMSats > 0:
        if iniFreqOption=='EQUAL':
            #Same frequency for every allele
            preOps.append( initByFreq( alleleFreq = [(1.0/maxAlleleN)] * maxAlleleN , loci=xrange(0, numMSats)        ))
        if iniFreqOption=='DIRICHLET':
            #Dirichlet distrubution
            diri = numpy.random.mtrand.dirichlet([1.0]*maxAlleleN,numMSats)
            preOps= [initByFreq(alleleFreq = diri[x].tolist(), loci=[x]) for x in xrange(numMSats)]
        if iniFreqOption=='FILE':
            #Read frecuencies data from a file
            loadData(filename)
            
            preOps= [initByFreq(alleleFreq = frecs[x][:].tolist(), loci=[x]) for x in xrange(numMSats)]

    if numSNPs > 0:    #TODO: Multiple options to initializa SNPs frequencies
        preOps.append( initByFreq( alleleFreq = [0.5, 0.5] + [0.0] * (maxAlleleN-2), loci=xrange(numMSats, numMSats+numSNPs)        ))

    preOps.append(initSex( sex=[Male,Female])) #There is no need to initialize sex outside of Freqs, but i prefer that
    inOps = []
    return loci, preOps, inOps


def createSim(pop, reps, mateOp):
    '''
    Establish basic parameters for simulator class (simuPOP)
    '''
    sim = simulator (
        pop,
        mateOp,
        rep = reps
    )
    return sim

def demo(gen, oldSize=[]):
    '''
    Demographic function to obtain the number of offspring generation.
    THIS IS A DURING MATING OPERATOR!
    '''
    global obsLambda
    global obsN0
    newborns = 0.0
    survivals = 0.0
    prop = []
    lambda_model=0.0
    if myTracer:
        print "Entering demo()"
        print "Age Struct in demo() (males, females): ", pop.vars()['ageStructM'], pop.vars()['ageStructF']
    #number of expected Deaths:
    pop.vars()['natDeaths'] = 0.0
    for x in xrange(maxAge+1):
        pop.vars()['natDeaths'] += ((1.0 - sFemale[x]) * float(pop.vars()['ageStructF'][x])) + ((1.0 - sMale[x]) * float(pop.vars()['ageStructM'][x]))
    pop.vars()['natDeaths']=round(pop.vars()['natDeaths'])

    newborns = cohortSize
    #If offspring per female is limited to 1, then there is a need of check that newborns can't be higher than mature females
    if matingMode == 'MONOGFEMS1':
        matureF = 0
        for x in xrange(minMatingAge, maxAge+1):
            matureF += pop.vars()['ageStructF'][x]
        if matureF<cohortSize:
            newborns=matureF
            if myTracer:
                print"More newborns than mature females, readjusting N0"

    if myTracer:
        print"expected Deaths: ",int(pop.vars()['natDeaths'])
    for x in xrange(maxAge):
        survivals += pop.vars()['ageStruct'][x]
        # Note this not include individual of the last Age (The ages are increased by one!) who will die
    #Check
    #print gen, oldSize[0],survivals
    if gen == 0:
        return numIndivs   #This cause of t0 = t1 (Not important)
    elif constSize==True or (constSize==False and gen<burnin):
        if obsN0>0:
            obsLambda=newborns/obsN0
        obsN0=newborns  #This is equal to cohortSize or matureF
        #print newborns, survivals+ newborns
        return survivals+ newborns
    newborns=0.0
    for x in xrange(maxAge+1):
        newborns += (bFemale[x] * float(pop.vars()['ageStructF'][x]))   # of offspring determined by females
    if matingMode == 'MONOGFEMS1' and newborns>matureF:
        newborns=matureF
    if obsN0>0:
        obsLambda=newborns/obsN0
    obsN0=newborns

    #if pop.vars()['ageStruct'][0]>0:
     #   lambda_model= newborns / float(pop.vars()['ageStruct'][0])
    #print "lambda = ",lambda_model
    
    if myTracer:
        print "survivals:", survivals
        print "newborns: ",newborns
        print "newpopsize: ", int(newborns+survivals)
        #print "Age Struct in demo() (males, females): ", pop.vars()['ageStructM'], pop.vars()['ageStructF']
    return [int(newborns+survivals)]

def checkDemo(pop):
    '''
    Computes, every time unit, the observed bx and lx in order to obtain the right T and calculate Ne by Felsenstein
    '''
    if myTracer:
        print "Entering checkDemo()"
    
    # computing observed sx (overall)
    sx = [0.0] * (maxAge+1)
    for x in xrange(maxAge+1):
        #print pop.vars()['ageStruct'][x], pop.vars()['removed'][x]
        if pop.vars()['ageStruct'][x] > 0:
            sx[x]= float(pop.vars()['ageStruct'][x] - pop.vars()['removed'][x]) / float(pop.vars()['ageStruct'][x])
        else:
            sx[x]=0.0001 
    
    # Computing lx
    for x in xrange(maxAge+1):
        if x == 0:
            pop.vars()['lx'][0]=1.0
            #pop.vars()['lxM'][0]=1.0 #For males   = lx
            pop.vars()['lxF'][0]=1.0 #For females
        else:
            pop.vars()['lx'][x]=pop.vars()['lx'][x-1]*sx[x-1]
            #pop.vars()['lxM'][x]=pop.vars()['lxM'][x-1]*pop.vars()['ageStructM'][x-1]
            pop.vars()['lxF'][x]=pop.vars()['lxF'][x-1]*pop.vars()['ageStructF'][x-1]
        #Here is the problem when culling of youngs is very strong (i.e. all individuals of a early age class are removed
        #In that case we are removing a whole cohort. The base probem is arisen because the lx calculed is assuming equilibrium in age structure.
    if myTracer:
        print""
        print "lx: ", pop.vars()['lx']
        print "Observed Offspring per age (all): ", pop.vars()['offspring'] #testing

    #Computing the observed T
    
#    s1=0.0
#    s2=0.0
#    for x in xrange(maxAge+1):
#        pop.vars()['removed'][x] = 0
#        s1 += x * pop.vars()['lx'][x] * pop.vars()['offspringF'][x]
#        s2 += pop.vars()['lx'][x] * pop.vars()['offspringF'][x]
#    if s2>0:
#        pop.vars()['T'] = s1/s2
#    else:
#        pop.vars()['T'] = T

    #Option 2. Computing observed T as average age of parents
    agesParents = []
    for x in xrange(pop.popSize()):
        if pop.individual(x).info('age')==0:
            agesParents.append(pop.individual(x).intInfo('ageFather'))
            agesParents.append(pop.individual(x).intInfo('ageMother'))
    #print agesParents
    Tobs= (float(sum(map(int,agesParents)))/float(len(agesParents)))
    #print Tobs
    pop.vars()['T'] = Tobs
    return True

def reportAge(pop):
    '''
    Store the age structure in the population
    '''
    if myTracer:
        print "Entering reportAge()"
    pop.vars()['ASprop'] = [0.0]*(maxAge+1)  #for saving in the pop namespace

    
    lst = [0] * (maxAge +1)
    lstM = [0] * (maxAge +1)
    lstF = [0] * (maxAge +1)
    #propAges = [0.0] * (maxAge +1)
    tot = 0.0
    totM = 0.0
    totF = 0.0
    for i in pop.individuals():
        lst[i.intInfo('age')] +=1
        tot += 1.0
        if i.sex()==1:
            lstM[i.intInfo('age')] +=1
            totM += 1.0
        else:
            lstF[i.intInfo('age')] +=1
            totF += 1.0
    
    for a in xrange(maxAge+1):
        pop.vars()['ageStruct'][a] = lst[a]
        pop.vars()['ageStructM'][a] = lstM[a]
        pop.vars()['ageStructF'][a] = lstF[a]
        pop.vars()['ASprop'][a]=float(lst[a])/tot
    
    if myTracer:
        print "Age structure: ", pop.vars()['ASprop'] #Reporting only overall age structure
        #print "Age structure (before deaths): ", ageStruct, sum(map(int,ageStruct)), pop.popSize()
    
    return True

def reportParentAge(pop):
    '''
    This function is just to check that bx rules!
    '''
    
    if myTracer:
        print "Entering reportParentAge()"
    pop.vars()['offspringM'] = [0.0]*(maxAge+1)
    pop.vars()['offspringF'] = [0.0]*(maxAge+1)
    pop.vars()['offspring'] = [0.0]*(maxAge+1)
    '''
    c = [0.0]*(maxAge+1)
    for i in pop.individuals():
        c[i.intInfo('age')] += 1.0
    '''
    for i in pop.individuals():
        a = i.intInfo('age')
        if a == 0:
            x = i.intInfo('ageMother')
            if x >= 0:
                if pop.vars()['ageStructF'][x]>0:
                    pop.vars()['offspringF'][x] += 1.0/(pop.vars()['ageStructF'][x]*2.0)
                    pop.vars()['offspring'][x] += 1.0/(pop.vars()['ageStruct'][x])
                    #print pa
                #print "mother ",x," of age ", pa
                x = i.intInfo('ageFather')
                if x >=0:
                    if pop.vars()['ageStructM'][x]>0:
                        pop.vars()['offspringM'][x] += 1.0/(pop.vars()['ageStructM'][x]*2.0)
                        pop.vars()['offspring'][x] += 1.0/(pop.vars()['ageStruct'][x])
                        #print pa
                        #print "father ",x," of age ", pa

    if myTracer:
        print "Offspring per age (female, male): ", pop.vars()['offspringF'], pop.vars()['offspringM']
    
    return True

def culling(pop):
    '''
    Simulation of culling procedure: the random kill of individuals when population size reaches some
    give threshold. A percetage of individuals is removed. Several possible regimes: ALL, JUVENILES, ADULTS
    regarding of which individuals could be removed
    '''
    #Culling is a Postmating operator, so ages are increased by 1 ????  BE CAREFUL WITH THIS IF CULLING IS MOVED
    total = pop.popSize()
    #Check if culling is allowed
    if myTracer:
        print "Entering culling()"
    if cullOpt == False or total < cullObjective:  #There is not sense to do culling if the population is lower than desired after culling
        pop.vars()['Culled'] = 0
        return True
    #Check if population size is above threshold if THRESHOLD mode.
    if cullMode=='THRESHOLD' and total < cullLevel:
        pop.vars()['Culled'] = 0
        return True
    #Check if it is time to do culling when INTERVAL mode is selected
    if cullMode=='INTERVAL' and (sim.gen() - burnin)%cullInterval != 0:
        pop.vars()['Culled'] = 0
        return True
    #Cull individuals by removing then from population. Store idx in remove
    remove = [] #Since revision 19, it estores all removable indivs and then disorder them and remove x
    
    for rn in xrange(pop.popSize()):
        if cullTarget=='ALL':  #Newborns are included
            remove.append(rn)
            
        elif cullTarget=='YOUNG' and (pop.individual(rn).intInfo('age') < cullAge ):
            remove.append(rn)
            
        elif cullTarget=='ADULTS' and (pop.individual(rn).intInfo('age') > cullAge ):
            remove.append(rn)

        elif cullTarget=='ADULTS10' and (pop.individual(rn).intInfo('age') > cullAge ): 
            remove.append(rn)
            
    if len(remove) > (total-cullObjective):
        shuffle(remove)
        remove=remove[:(total-cullObjective)]
        pass

    #Set removed by age
    for x in remove:
        pop.vars()['removed'][pop.individual(x).intInfo('age')] += 1
        if pop.individual(x).sex()==1:
             Vkmale.append(fatherslist.count(pop.individual(x).intInfo('ref')))
    pop.removeIndividuals(remove)
    pop.vars()['Culled'] = len(remove)
    '''After culling a new age structure should be calculated in order to obtain the right
    number of offspring given the new situation of the population'''
    #reportAge(pop) #Not needed if culling is postMating

    if myTracer:
        pass
    
    return True

def initAlphas(pop):
    ''' Set initial distribution of alpha males '''

    #determine number of alphable males:
    males = [x for x in xrange(pop.popSize()) if pop.individual(x).sex() == 1
        and pop.individual(x).info('age') >= alphaMinAge and pop.individual(x).info('age') <= alphaMaxAge]
    expected = int(numIndivs*alphaProp/2)
    #print males
    if expected <len(males):
        males = males[:(expected-1)]

    for x in males:
        pop.individual(x).setInfo(1,'alphaField')
    #print "alphas ", males
    return True

def checkAlphas(pop):  #TODO: Only if alpha mating is allowed
    if alphaMating==False:
        return True
    global propAlphas
    if myTracer:
        print "Entering checkAlphas()"
    expected = int(alphaProp*pop.dvars().numOfMale)
    
    # 'old' alphas
    count = 0
    for x in xrange(pop.popSize()):
        a = pop.individual(x).info('alphaField')
        if pop.individual(x).sex() == 1 and a > 0:
            if a < alphaDuration:
                pop.individual(x).setInfo((a+1),'alphaField')
                count += 1
            else:
                pop.individual(x).setInfo(-1,'alphaField') #remove from alpha list forever
    if myTracer:
        print " ", count, " alphas remaining of ", expected, " expected"
    propAlphas = float(count)
    # New alphas
    if count<expected:
        alphables = [x for x in xrange(pop.popSize()) if pop.individual(x).sex() == 1
        and pop.individual(x).info('age') >= alphaMinAge and pop.individual(x).info('age') <= alphaMaxAge
        and pop.individual(x).info('alphaField') == 0]

        if len(alphables)>0:
            new = expected - count
            if new > len(alphables):
                new = len(alphables)
            count = 0
            for x in xrange(new):
                i = randint(0,len(alphables))
                pop.individual(alphables[i]).setInfo(1, 'alphaField')
                count +=1
                del alphables[i]
                if len(alphables)==0:
                    break
            if myTracer:
                print " ", count, "new alphas"
    propAlphas += float(count)
    propAlphas = propAlphas/float(pop.dvars().numOfMale)
    
    return True

def saver(pop, expr):  #To be removed after extraction of saveGenepopAge call to a "safe" location
    #SaveFstat(pop,outputExpr='"out-%d-%d.txt" % (rep, gen)')

    saveGenepopAge(pop,outputExpr=expr)
    return True


def setFitness(pop):
    '''
    Set the called Fitness for simuPOP random mating. It is in fact the relative fecundity
    of each individual
    '''
    for i in pop.individuals():
        a = i.intInfo('age')  #Age were already increased so the real age is a-1

        if i.sex() == 1:
            i.setInfo(bMale[a],'fitness')
        else:
            i.setInfo(bFemale[a],'fitness')
   
    return True
    

def zip(gen):
    '''zip the genepop files to save disk space'''
    gen=int(str(gen)) #Normalize time number
    gen= gen -burnin
    files =  filter(lambda x : x.endswith('.gen'),
            os.listdir(filePrefix + os.sep + 'sims' + os.sep + str(gen)))
    for file in files:
        bz=bz2.BZ2File(filePrefix + os.sep + 'sims' + os.sep + str(gen) + os.sep + file + '.bz2', 'w')
        f=open(filePrefix + os.sep + 'sims' + os.sep + str(gen) + os.sep + file )
        l = f.readline()
        while l<>'':
            bz.write(l)
            l = f.readline()
        f.close()
        bz.close()
        os.remove(filePrefix + os.sep + 'sims' + os.sep + str(gen) + os.sep + file)

def zipper(pop):
    if sim.gen()>1:
        try:
            zip(sim.gen() - 2)
        except:
            print "zipper error in ", sim.gen()-2
    return True

def evolveSim(sim, gens, numMarkers,
        genPreOps, genInOps, popPreOps, popInOps, reportOps,
        oExpr, saveGens):
    sim.evolve (
        preOps = genPreOps + popPreOps,
        ops = genInOps + popInOps + reportOps+  traceOps + [parentsTagger()], 
#            [recombinator(rate=0.1),],
        gen = gens,

    )



def naturalDeath(pop):
    '''
    This function tries to simulate the demographic survival from every age class & sex.
    '''
    if myTracer:
         print "Entering naturalDeath()"
    dead = []
    props = [0.0]*(maxAge+1) #Proportion of survivals per age
 
    #Calculating number of deaths (the same as in demo())
    pop.vars()['natDeaths'] = 0.0
    for x in xrange(maxAge+1):
        pop.vars()['natDeaths'] += ((1.0 - sFemale[x]) * float(pop.vars()['ageStructF'][x])) + ((1.0 - sMale[x]) * float(pop.vars()['ageStructM'][x]))
    pop.vars()['natDeaths']=round(pop.vars()['natDeaths'])

    sum = [0.0]*(maxAge+1)
    #PostMating operator, so ages are increased by one
    if sim.gen()>0:
        d=0
        i=0
        
        while d<int(pop.vars()['natDeaths']):

            p = pop.individual(i).intInfo('age') - 1
            if p>=0:
                sum[p] +=1.0

                if pop.individual(i).sex()==1 and uniform(0,1) > sMale[p]:
                    if dead.count(i)==0:
                        dead.append(i)
                        d += 1
                elif pop.individual(i).sex()==2 and uniform(0,1) > sFemale[p]:
                    if dead.count(i)==0:
                        dead.append(i)
                        d += 1
                else:
                    props[p] +=1.0
            i+=1
            if i==pop.popSize():
                i=0
        for i in xrange(maxAge+1):
            if sum[i] > 0:
                props[i] = props[i]/sum[i]
            else:
                props[i]=0
        #Set removed by age
        for x in dead:
            pop.vars()['removed'][pop.individual(x).intInfo('age')-1] += 1
            if pop.individual(x).sex()==1 and sim.gen()>burnin:
                Vkmale.append(fatherslist.count(pop.individual(x).intInfo('ref')))
                #Vkmale.append(pop.individual(x).intInfo('numOffs'))
        pop.removeIndividuals(dead)

        if myTracer:
            #natD = len(dead)
            print "Natural deaths: ", len(dead)
            print "survivorship per age: ", props
    return True

def incAge(pop):
     '''
        Temporary function to increase age and get the age of parents.
        Used with simuPOP 0.9.2 due to a bug avoiding the right use of InfoExec rountine
     '''
     if myTracer:
         print "Entering incAge()"
     # Increase age of individuals and save current age in ageParents infoFields
     ageIdx = pop.infoIdx('age')
     for ind in pop.individuals():
        ind.setInfo(ind.info(ageIdx), 'ageFather')
        ind.setInfo(ind.info(ageIdx), 'ageMother')
        ind.setInfo(ind.info(ageIdx) + 1, ageIdx)
     return True

def doAge(pop, maxAge, minMatingAge):
    """
    I need to explain what happens here
    """
    #NOTE that ages were increased before reach this.
    
    #Population divided into groups according to their age
    # Two groups: (0,0) age <=maxAge; (0,1) age > maxAge, individuals that should die after mating
    pop.setVirtualSplitter(infoSplitter('age', cutoff=[maxAge + 0.1]))
    # Clone the whole population, even the last age individuals, who should be removed
    clonePop= [cloneMating(subPop=(0, 0), weight=-1, ops=[inheritTagger(infoFields=['age']),inheritTagger(infoFields=['fitness']), inheritTagger(infoFields=['ref']), inheritTagger(infoFields=['alphaField'])])]
    if matingMode == 'RANDOM':
        # mate the mature individuals. All individuals are in the pool to be parents, but as juveniles have fitness=0 they not be choosen
        matePop = [randomMating(selectionField='fitness')]
    if matingMode == 'MONOGFEMS1':
        matePop = [homoMating(pyParentsChooser(generatorMonogFem), mendelianOffspringGenerator(numOffspring=1), weight=1)]

    #Combine both mating schemes to produce next generation
    mateList = ( clonePop +  matePop)
  
    return heteroMating(mateList, subPopSize=demo, shuffleOffspring=False)

def generatorMonogFem(pop,sp):
    '''
    Chooserparents for species with (or without) alpha males and where females are restricted to a unique mating
    Like Bison
    '''
    if myTracer:
        print "generator() initialization"
    males = [x for x in xrange(pop.popSize()) if pop.individual(x).sex() == 1
        and pop.individual(x).info('age') >= minMatingAge and pop.individual(x).info('alphaField')==0]
    alphaMales = [x for x in xrange(pop.popSize()) if pop.individual(x).sex() == 1
        and pop.individual(x).info('age') >= minMatingAge and pop.individual(x).info('alphaField')>0]
    females = [x for x in xrange(pop.popSize()) if pop.individual(x).sex() == 2
        and pop.individual(x).info('age') >= minMatingAge]
    ages=map(int, pop.indInfo('age'))
    ages = [ages[x]-1 for x in xrange(len(ages))] #Yeah, because ages are increased by 1 during mating
    if myTracer:
        print "males ", males
        print "alphamales ", alphaMales
        print "females ", females
        print ages
    
    choosen = []
    #Set probablity for each alpha male.
    if len(alphaMales) > 0:
        probA = [ bMale[ ages[x] ] for x in alphaMales  ]
        suma = sum(map(float,probA))
        probA = [ probA[x]/suma for x in xrange(len(alphaMales)) ]
        for x in xrange(1,len(probA)):
            probA[x] = probA[x]+probA[x-1] #Cumulative prob
    #Set probablity for normal male.
    if len(males) > 0:
        probM = [ bMale[ ages[x] ] for x in males  ]
        suma = sum(map(float,probM))
        probM = [ probM[x]/suma for x in xrange(len(probM)) ]
        for x in xrange(1,len(probM)):
            probM[x] = probM[x]+probM[x-1]   #Cumulative prob
   
    while True:
        #Set probablity for each  remaining female.
        
        probF = [ bFemale[ ages[x] ] for x in females  ]
        suma = sum(map(float,probF))
        probF = [ probF[x]/suma for x in xrange(len(probF)) ]
        
        if len(females)>1:
            for x in xrange(1,len(probF)):
                probF[x] = probF[x]+probF[x-1]
        if len(females)<1: #no more females available
            yield (-99, -99) #error codes
        #random choose of a female
        r = uniform(0,1)
        for x in xrange(len(females)):
            if probF[x]>r:
                female = females[x]
                del females[x] #removing that female for next rounds  PROBLEM!!! need check for avoiding empty females
                break
        #Choose the set of males to be used
        #if uniform(0,1)< alphaOff:
        if alphaMating and r < alphaOff:
            #Choose a random father between alpha males
            r = uniform(0,1)
            for x in xrange(len(alphaMales)):
                if probA[x]>r:
                    male = alphaMales[x]
                    break
        else:
             #Choose a random father between non-alpha males
         
            r = uniform(0,1)
            for x in xrange(len(males)):
                if probM[x]>r:
                    male = males[x]
                    break
        if myTracer:
            pass
            print female, male
       
        fatherslist.append(pop.individual(male).intInfo('ref'))
#        print pop.individual(male).intInfo('ref'), pop.individual(male).intInfo('father_ref'), pop.individual(male).intInfo('numOffs')
        yield(male,female) #The two individuals passed to the pyMating
        
        


def setRefs(pop, refer):
    if myTracer:
        print "Entering setRefs()"
    for i in xrange(pop.popSize()):
            if pop.individual(i).intInfo('ref') == 0:
                pop.individual(i).setInfo(refer[0],'ref')
                refer[0] += 1
    return True






def saveGenepopAge(pop, output='', outputExpr='', maxAllele=0, loci=[], shift=1,
    combine=None):
    '''
    Save population in Genepop format, including info on individual age added to the label.
    Adapted from saveFstat.
    '''
    if output != '':
        file = output
    elif outputExpr != '':
        file = eval(outputExpr, globals(), pop.vars() )
    else:
        raise exceptions.ValueError, "Please specify output or outputExpr"
    # open file
   
    try:
        f = open(file, "w")
    except exceptions.IOError:
        raise exceptions.IOError, "Can not open file " + file + " to write."
    #
    # file is opened.

    np = pop.numSubPop()
    
    if loci == []:
        loci = xrange(pop.totNumLoci())
    nl = len(loci)
    
    if maxAllele != 0:
        nu = maxAllele
    else:
        nu = maxAlleleN #pop.maxAllele()
    
    nd = 2
    # write the first line
    f.write( 'Sim %d %d %d %d\n' % (np, nl, nu, nd) )
    # following lines with loci name.
    for loc in loci:
        f.write( pop.locusName(loc) +"\n");
    gs = pop.totNumLoci()
    for sp in xrange(0, pop.numSubPop()):
        # genotype of subpopulation sp, individuals are
        # rearxranged in perfect order
        # here is where the age classes infoField should be
        f.write('Pop\n')
        gt = pop.genotype(sp)
        for ind in xrange(0, pop.subPopSize(sp)):
            i = pop.individual(ind, sp)
            f.write("%d(%d), " % (sp+1,int(i.info('age'))))
            p1 = 2*gs*ind          # begining of first hemo copy
            p2 = 2*gs*ind + gs     # second
            for al in loci: # allele
                # change from 0 based allele to 1 based allele
                if combine is None:
                    ale1 = gt[p1+al] + shift
                    ale2 = gt[p2+al] + shift
                    f.write('%%0%dd%%0%dd ' % (nd, nd) % (ale1, ale2))
                else:
                    f.write('%%0%dd' % nd % combine([gt[p1+al], gt[p2+al]]))
            f.write( "\n")
    f.close()





def obtainAgeParents(fields):
    '''
    get the infoFields ageFather from both parents. Despite of the name of the infoField, it represents,
    at mating time, the individual's age. The returning values are stored in ageFather and ageMother of the new
    individual.
    '''
    if len(fields) > 2:  #There is a mother (Individuals cloned don't have :)
        agesParents = [fields[0]]+[fields[2]]
    else:
        agesParents = [fields[0]]+[-1]
    return agesParents

def obtainRefParents(fields):
    #print fields
    if len(fields) > 2:  #There is a mother
        refParents = [fields[0]]+[fields[2]]
    else:
        refParents = [fields[0]]+[-1]
    return refParents

def statHe(pop):
    He = [pop.dvars().expHetero[x] for x in xrange(numSNPs + numMSats) ]
    aveHe = sum(map(float, He)) / float(len(He))
    return aveHe

def statHo(pop):
    Ho = pop.dvars().HeteroFreq
    aveHo = sum(map(float, Ho)) / float(len(Ho))
    return aveHo

def statA(pop):
    A = [pop.dvars().numOfAlleles[x] for x in xrange(numSNPs + numMSats) ]
    aveA = sum(map(float, A)) / float(len(A))
    return aveA

def statNeFel(pop):
    ''' 
    Calculate Ne following Felsenstein 1971, equation 24.
    see Waples and Yokota 2007 for details.
    '''
    '''
    Ne = sum(V(k)) / (1+sum(Fels(k))

    V(k)= v(k)*N(k)
    v(k) = q(k)/l(k)
    q(k) = sum[k->]( l(k)*b(k) ) ( sum( l(k)*b(k))

    Fels(k) = l(k)*s(k)*(1-s(k))*v(k+1)*v(k+1)
    '''
    V = [0.0]*(maxAge+1)
    N = [0.0]*(maxAge+1)
    v = [0.0]*(maxAge+1)
    q = [0.0]*(maxAge+1)
    b = [0.0]*(maxAge+1)
    s = [0.0]*(maxAge+1)  # Obtener a apartir del lx observado
    lxb = [0.0]*(maxAge+1)
    Fels = [0.0]*(maxAge+1)
    for k in xrange(maxAge+1):
        b[k]= pop.vars()['offspringF'][k]
        lxb[k] = pop.vars()['lx'][k] * b[k]
        N[k] = pop.vars()['ageStruct'][0] * pop.vars()['lx'][k]
        if pop.vars()['lx'][k] > 0 and k<maxAge:
            s[k] = pop.vars()['lx'][k+1]/pop.vars()['lx'][k]
        else:
            s[k]=0.0
    for k in xrange(maxAge+1):
        try:
            q[k] = sum(map(float, lxb[k:])) / sum(map(float, lxb))
        except ZeroDivisionError:
            q[k]=0.0
        if pop.vars()['lx'][k] > 0:
            v[k] = q[k] / pop.vars()['lx'][k]
        else:
            v[k]=0.0
        V[k] = v[k]*N[k]
        
    for k in xrange(maxAge):
        Fels[k] = pop.vars()['lx'][k]*s[k]*(1-s[k])*v[k+1]*v[k+1]

    #Calculating the current Generation lenght (T)
#    pop.vars()['T'] = sum(map(float,q))
    Ne = sum(map(float, V)) / (1.0 + sum(map(float, Fels)) )
    
    return Ne

def acumStats(pop):
    global obsLambda, propAlphas, obsN0
    t = sim.gen()
    if myTracer:
        print "acumStats()"
    
    genN[t][repli]=pop.popSize()
    genN0[t][repli]=obsN0
    genMales[t][repli]=pop.dvars().propOfMale
    genLambda[t][repli] = obsLambda
    genAlphas[t][repli] = propAlphas
    genA[t][repli]=statA(pop)
    #Allelic diversity loss
    if t >= burnin:
        genAD[t][repli] = (genA[t][repli] - 1.0) / (genA[burnin][repli] - 1.0)
    genHe[t][repli]=statHe(pop)
    genHo[t][repli]=statHo(pop)
    if t >= burnin:
        genHeL[t][repli] = genHe[t][repli] / genHe[burnin][repli]
        genHoL[t][repli] = genHo[t][repli] / genHo[burnin][repli]
    if t>0:
        genNeFel[t][repli]=statNeFel(pop)
    else:
        genNeFel[t][repli]=-99.0
    for a in xrange(maxAge+1):
        genAS[t][a] += pop.vars()['ASprop'][a]

    if cullOpt:
        genCulled[t][repli] = pop.vars()['Culled']
    #Observed T (following Carey 1993)
    genT[t][repli]= pop.vars()['T']
    T = round(pop.vars()['T'])
    # NeH: Ne obtained from He loss (from burnin)
    if t>=(T):
        try:
            genNeH[t][repli]= genHe[t-T][repli] / (2*(genHe[t-T][repli] - genHe[t][repli]))
        except:
            print "Error! #0002"
            print t, T, repli
    if t==(gens-1):
        try:
            globalNeH[repli] +=  1.0/( -2.0 * exp( ( log(genHe[t][repli]) - log(genHe[burnin][repli])) / ((t-burnin)/T)  ) + 2.0)
        except:
            print "Error! #0004"
        if repli == reps-1:
            #globalNeH[0] = globalNeH[0]/float(reps)
            try:
                print 'Overall NeH:  %.2f  (%.2f) ' % (stats.mean(globalNeH), stats.sem(globalNeH) )
            except:
                print "Error!  #0001"
        
    return True

def outputStats():
    '''
    Print a table on screen.
    Saves the main table of results in a text file.
    
    '''
    aveN=genN.mean(1)
    aveN0=genN0.mean(1)
    aveLambda=genLambda.mean(1)
    aveMales=genMales.mean(1)
    aveAlphas=genAlphas.mean(1)
    aveA=genA.mean(1)
    aveAD=genAD.mean(1)
    aveHe=genHe.mean(1)
    aveHo=genHo.mean(1)
    aveHeL=genHeL.mean(1)
    aveHoL=genHoL.mean(1)
    aveNeH=genNeH.mean(1)
    aveNeFel=genNeFel.mean(1)
    aveT=genT.mean(1)
    if cullOpt and myTracer:
        masked_Culled = ma.masked_values(genCulled, 0)
        aveCulled = masked_Culled.mean(1)
        print "averaged culled individuals per year (only if culling happens)"
        print aveCulled

    #Standard errors of the mean
    stdN=genN.std(1)/sqrt(reps)
    stdN0=genN0.std(1)/sqrt(reps)
    stdLambda=genLambda.std(1)/sqrt(reps)
    stdMales=genMales.std(1)/sqrt(reps)
    stdAlphas=genAlphas.std(1)/sqrt(reps)
    stdA=genA.std(1)/sqrt(reps)
    stdAD=genAD.std(1)/sqrt(reps)
    stdHe=genHe.std(1)/sqrt(reps)
    stdHo=genHo.std(1)/sqrt(reps)
    stdHeL=genHeL.std(1)/sqrt(reps)
    stdHoL=genHoL.std(1)/sqrt(reps)
    stdNeH=genNeH.std(1)/sqrt(reps)
    stdNeFel=genNeFel.std(1)/sqrt(reps)
    stdT=genT.std(1)/sqrt(reps)

    aveAS = [0.0]*(gens)

    for t in xrange(gens):
       aveAS[t] = genAS[t]/reps
       

    

    f = open(filePrefix+os.sep+filePrefix+"-RESULTS.csv","w")
    f.write('Simulation Name:;;  %s \n' % (filePrefix) )
    f.write('Loci (msats);%d;Loci (SNPs);%d;Initialization;%s;' % (numMSats, numSNPs, iniFreqOption) )
    if iniFreqOption == 'FILE':
        f.write('%s\n' % (filename))
    else:
        f.write('max Alleles per locus; %d\n' % (maxAlleleN))
    f.write('replicates;%d;;time units;%d;;burning;%d\n' % (reps, gens-burnin-1, burnin))
    f.write('\nDemographic parameters\n')
    f.write('initial N;%d;;N0;%d;' % (numIndivs, cohortSize))
    if constSize:
        f.write('Constant for all evolution (Stable population size)\n')
    else:
        f.write('used only during burn-in\n')
    f.write('Ages:;;')
    for k in xrange(maxAge+1):
        f.write('%d;' % (k))
    f.write('\n')
    f.write('Survival (males):;;')
    for k in xrange(maxAge+1):
        f.write('%.3f;' % (sMale[k]))
    f.write('\n')
    f.write('Survival (females):;;')
    for k in xrange(maxAge+1):
        f.write('%.3f;' % (sFemale[k]))
    f.write('\n')
    f.write('Birth rate (females):;;')
    for k in xrange(maxAge+1):
        f.write('%.3f;' % (bFemale[k]))
    f.write('\n')
    f.write('Mating prob (males):;;')
    for k in xrange(maxAge+1):
        f.write('%.3f;' % (bMale[k]))
    f.write('\n')
    if alphaMating:
        f.write('Mating mode:;;%s;%s\n' % (matingMode,'with alpha males'))
    else:
        f.write('Mating mode:;;%s;%s\n' % (matingMode,'without alpha males'))
    #Tables of alpha and culling parameters
    if alphaMating:
        f.write('\nAlpha Mating parameters\n')
        f.write('perc. offspring:;;%.1f;;' % (alphaOff))
        f.write('perc. of alphas:;;%.1f\n' % (alphaProp))
        f.write('Min age:;%d;;' % (alphaMinAge))
        f.write('Max age:;%d;;' % (alphaMaxAge))
        f.write('duration:;%d;;' % (alphaDuration))
    if cullOpt:
        f.write('\nCulling parameters\n')
        f.write('objective:;;%d;;' % (cullObjective))
        f.write('target:;;%s\n' % (cullTarget))
        f.write('Mode:;%s;;' % (cullMode))
        if cullMode=='LEVEL':
            f.write('level:;%.2f\n' % (cullLevel))
        if cullMode=='INTERVAL':
            f.write('interval:;%d\n' % (cullInterval))
    try:
        f.write('\nOverall NeH:;%.2f;stderr;%.2f \n' % (stats.mean(globalNeH), stats.sem(globalNeH) ))
    except:
        print "Error!  #0001b"
    try:
        f.write('Overall k(males):;%.2f;stderr;%.2f \n' % (stats.mean(globalK), stats.sem(globalK) ))
        f.write('Overall Vk(males):;%.2f;stderr;%.2f \n' % (stats.mean(globalVk), stats.sem(globalVk) ))
        f.write('Overall prop of males siring:;%.4f;stderr;%.4f \n' % (stats.mean(globalDads), stats.sem(globalDads) ))
        f.write('\nSummary TABLE for saved generations\n')
        f.write("t;N;N(se);PropMales;PropMales(se);AD;AD(se);He-Loss;He-Loss(se);Ho-Loss;Ho-Loss(se);T;T(se);Ne[Fel];Ne[Fel](se);Ne[LDNe]; Ne[LDNe](-CI);Ne[LDNe](+CI)\n")
    except:
         print "Error!  #0002"
    for t in saveGens:
        f.write('%d;%.1f;%.1f;%.4f;%.4f;%.2f;%.2f;%.4f;%.4f;%.4f;%.4f;%.1f;%.1f;%.1f;%.1f\n' % ((t-burnin),
            aveN[t], stdN[t], aveMales[t], stdMales[t], aveAD[t], stdAD[t],aveHeL[t], stdHeL[t],
            aveHoL[t], stdHoL[t],  aveT[t], stdT[t], aveNeFel[t], stdNeFel[t]) )
    f.write('\n\nFull results\n\n')
    f.write("t;N;N(se);N0;N0(se);lambda;lambda(se);PropMales;PropMales(se);alphas;alphas(se);A;A(se);AD;AD(se);He;He(se);HeD;HeD(se);Ho;Ho(se);HoD;HoD(se);T;T(se);Ne[Fel];Ne[Fel](se);Ne[H]; Ne[H](se)\n")
    for t in xrange(burnin,gens):
         f.write('%d;%.1f;%.1f;%.1f;%.1f;%.3f;%.3f;%.4f;%.4f;%.2f;%.2f;%.2f;%.2f;%.4f;%.4f;%.4f;%.4f;%.4f;%.4f;%.4f;%.4f;%.4f;%.4f;%.1f;%.1f;%.1f;%.1f;%.1f;%.1f\n' % ((t-burnin),
            aveN[t], stdN[t],aveN0[t], stdN0[t],aveLambda[t], stdLambda[t], aveMales[t], stdMales[t], aveAlphas[t], stdAlphas[t],
            aveA[t], stdA[t], aveAD[t], stdAD[t],aveHe[t], stdHe[t],aveHeL[t], stdHeL[t],
            aveHo[t], stdHo[t], aveHoL[t], stdHoL[t], aveT[t], stdT[t], aveNeFel[t], stdNeFel[t], aveNeH[t], stdNeH[t]) )
    f.close()

#    print globalK
    print 'Overall k(males):  %.2f (%.2f) ' % (stats.mean(globalK), stats.sem(globalK) )
    print 'Overall Vk(males):   %.2f  (%.2f)' % (stats.mean(globalVk), stats.sem(globalVk) )
#    print globalDads
    print 'Overall prop of males siring:   %.3f  (%.3f)' % (stats.mean(globalDads), stats.sem(globalDads) )
#

    if saveDecays:
        #Save data for all replicates in saved generations for AD and HeL
        #From revision 27 also for Ne[H] and N
        f = open(filePrefix+os.sep+filePrefix+"-AData.dat","w")
        for t in saveGens:
            for r in xrange(reps):
                f.write('%d %d %f\n' % (t-burnin,r,genAD[t][r]))
        f.close()
        f = open(filePrefix+os.sep+filePrefix+"-HeData.dat","w")
        for t in saveGens:
            for r in xrange(reps):
                f.write('%d %d %f\n' % (t-burnin,r,genHeL[t][r]))
        f.close()
        f = open(filePrefix+os.sep+filePrefix+"-NData.dat","w")
        for t in saveGens:
            for r in xrange(reps):
                f.write('%d %d %f\n' % (t-burnin,r,genN[t][r]))
        f.close()
        f = open(filePrefix+os.sep+filePrefix+"-NeData.dat","w")
        for t in saveGens:
            for r in xrange(reps):
                f.write('%d %d %f\n' % (t-burnin,r,genNeH[t][r]))
        f.close()

        f = open(filePrefix+os.sep+filePrefix+"-globalNeData.dat","w")
        for r in xrange(reps):
            f.write('%d %d %f\n' % (t-burnin,r,globalNeH[r]))
        f.close()

    #Draw charts, if Pylab is loaded
    if usePylab and myTracer==False:
        print "Generating figures..."
        maxN= ( int( round(max(map(float,aveN[burnin:]))) / 100 ) + 1 ) * 100
        plotStats(aveN[burnin:], "N", "t","Population size",filePrefix,"-PopSize",0,maxN)
        plotStats(aveAD[burnin:], "AD", "t","Allelic Diversity",filePrefix,"-ADecay")
        maxNeFel= ( int( round(max(map(float,aveNeFel[burnin:]))) / 100 ) + 1 ) * 100
        maxNeH = ( int( round(max(map(float,aveNeH[burnin:]))) / 100 ) + 1 ) * 100
        plotStats(aveNeH[burnin:], "Ne", "t","Ne (He loss)",filePrefix,"-NeH", 0, maxNeH)
        plotStats(aveNeFel[burnin:], "Ne", "t","Ne (Felsenstein's)",filePrefix,"-NeFel", 0, maxNeFel)
        plotStats2(aveHeL[burnin:], aveHoL[burnin:], "H(t)/H(0)", "t","He & Ho diversity","He","Ho",filePrefix,"-HDecay",0, 1.2)
        maxNe = max(maxNeFel,maxNeH)
        plotStats2(aveNeH[burnin:], aveNeFel[burnin:], "Ne", "t","Effective size","Ne(He)","Ne(Fel)",filePrefix,"-Ne",0, maxNe)
        if len(saveGens)>1:
            plotBoxplot(genAD, "AD", "t","Allelic diversity",filePrefix,"-AD-boxplot")
            plotSuperBoxplot(genAD,aveAD[burnin:], "AD", "t","Allelic diversity",filePrefix,"-AD-SBP")
            plotBoxplot(genHeL, "H(t)/H(0)", "t","He diversity",filePrefix,"-He-boxplot")
            plotBoxplot(genHoL, "H(t)/H(0)", "t","Ho diversity",filePrefix,"-Ho-boxplot")
            plotBoxplot(genN, "N", "t","Population Size",filePrefix,"-PopSize-boxplot")
        

#        for x in xrange(gens):
#            relHe.append(aveHe[x]/aveHe[burnin])
#            relHo.append(aveHo[x]/aveHo[burnin])
#            relA.append((aveA[x]-1)/(aveA[burnin] - 1))
#        plotStats2(relHe[burnin:], relHo[burnin:], "H", "t","He & Ho",filePrefix,"-HDecay")
#        plotStats(relA[burnin:], "AD", "t","Allele Diversity Decay",filePrefix,"-ADecay")
        plotAS(aveAS[burnin:], filePrefix, "-Ages")
        print "...finished"

    return True

def plotStats(stat, yLabel, xLabel, title, prefix, short, lmin=0, lmax=1 ):
    '''
    plot any List on its length
    '''
    timescale = [x for x in xrange(len(stat))]
    fig = figure(1)
    fig.suptitle(title)
    xlabel(xLabel)
    ylabel(yLabel)
    plot(timescale, stat)
    ylim(lmin,lmax)
    savefig(prefix+os.sep+prefix+short+".png")
    #show()
    close(1)
    return True

def plotBoxplot(stat, yLabel, xLabel, title, prefix, short, lmin=0, lmax=1 ):
    '''

    '''
    data = []
    gener = []
    fig = figure(1)
    for g in saveGens:
        if (g-burnin)>0:
            data.append(stat[g,:])
            gener.append(g-burnin)
            fig.suptitle(title)
    xlabel(xLabel)
    ylabel(yLabel)
    boxplot(data,0,positions=gener)
    #ylim(lmin,lmax)
    savefig(prefix+os.sep+prefix+short+".png")
    close(1)
    return True

def plotSuperBoxplot(stat1, stat2, yLabel, xLabel, title, prefix, short, lmin=0, lmax=1 ):
    '''

    '''
    data = []
    gener = []
    timescale = [x for x in xrange(len(stat2))]
    fig = figure(1)
    for g in saveGens:
        data.append(stat1[g,:])
        gener.append(g-burnin)
        fig.suptitle(title)
    xlabel(xLabel)
    ylabel(yLabel)
    boxplot(data,0,positions=gener)
    plot(timescale, stat2)
    #ylim(lmin,lmax)
    savefig(prefix+os.sep+prefix+short+".png")
    close(1)
    return True

def plotAS(stat, prefix, short ):
    timescale = xrange(gens-burnin)
    fig = figure(1)
    fig.suptitle('Age Structure')
    xlabel('t')
    #axis.set_ylim(0.0,1.0) #trying to set the axis scale. Check matplotlib documentation
    for a in xrange(maxAge+1):
        line=[]
        for t in xrange(gens-burnin):
            line.append(stat[t][a])
        
        plot(timescale, line, label=str(a))
    #legend(loc=0)
    savefig(prefix+os.sep+prefix+short+".png")
    #show()
    close(1)
    return True



def plotStats2(stat1, stat2, yLabel, xLabel, title, label1, label2, prefix, short, lmin=0.0, lmax=1.0 ):
    timescale = [x for x in xrange(len(stat1))]
    fig = figure(2)
    fig.suptitle(title)
    xlabel(xLabel)
    ylabel(yLabel)
    
    plot(timescale, stat1, 'b', label=label1)
    plot(timescale, stat2, 'g', label=label2)
    ylim(lmin,lmax)
    legend()
    savefig(prefix+os.sep+prefix+short+".png")
    #show()
    close(2)
    return True

def obtainNe(pop):
    He = pop.dvars().expHetero
    meanExpHe = (sum(map(float,He))/len(He))
    print "He = ", meanExpHe, "  He0= ", He0[0]
    g= sim.gen()
    if g==0:
        He0[0]=meanExpHe
    meanExpHe = meanExpHe / param[0]
    print "He = ", meanExpHe
    t=g
    if t > 0:
        try:
            Ne =  (1.0 / (2.0*(1.0 - exp(log(meanExpHe)/float(t))) ) )
        except:
            Ne = -99.0
    else:
        Ne = 0.0
    print "Ne(1) = ", Ne
    return True

def allowSel(geno, gen, fields):
    return fields[0]

def loadData(filename):
    """
    Reads frequency data from a file
    """
    try:
        f = open(filename)
    except exceptions.IOError:
        print "Can not open file " + filename + " to read. Exiting"
        sys.exit(2)

    l=0
    s = f.readlines()
    for a in xrange(len(s)):
        ss = s[a].split()
        ss = ss[1:] #Ignore the locus name
        for i in xrange(len(ss)):
            frecs[l][i]=float(ss[i])
            #print float(ss[i])
        l += 1

    
    f.close()
    return True

def checkOptions(opt, par1, par2=None, par3=None):
    """
    ON PROGRESS
    """

    if opt == 'filame':
        if par2=='FILE': #Then the file must exist
            return True
    return True


if __name__ == '__main__':
    

    myTracer = False  #Make this variable True for print some tracers. Do not use with big numbers of individuals or replicates
    traceOps = []
    if myTracer == True:
        TurnOnDebug(DBG_MATING)
        #traceOps = [dumper(structure=False, stage=PostMating, infoFields=[ 'AlphaField'])]

   
  
#    (filePrefix, reps, gens, burnin, numMSats, maxAlleleN, iniFreqOption, numSNPs, numIndivs, constSize,
#        saveGens, maxAge, minMatingAge, ASinit, matingMode,
#        sMale, sFemale, bMale, bFemale,
#        cullOpt, cullObjective, cullTarget, cullMode, cullLevel, cullInterval,
#        alphaMating, alphaOff, alphaMinAge, alphaMaxAge, alphaDuration, alphaProp) = getOptions()
    #Get parameters
    pars = simuOpt(options,'GAS: Genetic Age-structured Simulations   v.0.4.2 (10/09/2009)')
    if not pars.getParam(nCol=3):
        sys.exit(1)
    #Automatically save configurations
    name = pars.name
    if not os.path.isdir(name):
        os.makedirs(name)
    pars.saveConfig( os.path.join(name, name+'.cfg'))

    #assign parameters (temp)
    filePrefix = pars.name
    reps = pars.reps
    gens = pars.time
    burnin = pars.burnin
    numMSats = pars.nMSats
    maxAlleleN = pars.maxAllele
    iniFreqOption = pars.initFreq
    if iniFreqOption == 'FILE':
        filename = pars.fileFreq
    
    numSNPs = pars.nSNPs
    numIndivs = pars.initialN
    cohortSize = pars.cohortN
    constSize = pars.constSize
    saveGensTemp = pars.saveGens
    maxAge= pars.maxAge
    minMatingAge= pars.adultAge
#    ASinit= pars.ASinit
    ASinit=eval(pars.ASinit)
    matingMode= pars.matingMode
    sMale= eval(pars.sMales)   #TODO: evaluate this to avoid wrong inputs
    sFemale= eval(pars.sFemales)
    bMale= eval(pars.bMales)
    bFemale= eval(pars.bFemales)
    cullOpt= pars.culling
    cullObjective= pars.cullObjective
    cullTarget= pars.cullTarget
    cullMode= pars.cullMode
    cullLevel= pars.cullLevel
    cullInterval= pars.cullInterval
    alphaMating= pars.alphaMating
    alphaOff= pars.alphaOff
    alphaMinAge= pars.alphaMinAge
    alphaMaxAge= pars.alphaMaxAge
    alphaDuration= pars.alphaDuration
    alphaProp= pars.alphaProp
    saveDecays=pars.saveDecays
    cullAge=pars.cullAge
    
    # A checking of inputs should be done here! (TODO)
    # Setting generations to ve saved
    saveGens=saveGensTemp.split(',')
    try:
        for x in xrange(len(saveGens)):
            if int(saveGens[x])>gens:
                raise OptionError, ('years to save','You\'ve tried to save a year beyond the limit')
            saveGens[x] = int(saveGens[x])+burnin
    except OptionError, e:
        print "Error parsing parameters in ",e.param,": ",e.message
        exit(0)
    except:
        saveGens = [0]
    gens = gens + burnin +1
    frecs = numpy.zeros((numMSats, maxAlleleN))
    

    He0 = 0.0 #To store the He(0) to standarize the He's in every replicate(Global)
    sumHe = 0.0

    
   
    # Output headers. 
    
    print 'Simulation Name:  %s ' % (filePrefix)
    print '******************************************'
    print 'Loci (msats): %d\tLoci (SNPs): %d\tInitialization:  %s' % (numMSats, numSNPs, iniFreqOption),
    if iniFreqOption == 'FILE':
        print ' -> %s' % (filename)
    else:
        print ' (max Alleles per msat locus: %d)' % (maxAlleleN)
    print 'replicates: %d\ttime units: %d\tburning: %d' % (reps, gens-burnin-1, burnin)
    print '\nDemographic parameters'
    print '----------------------'
    print 'initial N: %d\tN0: %d  ' % (numIndivs, cohortSize),
    if constSize:
        print 'Constant for all evolution (Stable population size)'
    else:
        print 'used only during burn-in'
    print '\nAges:\t',
    for k in xrange(maxAge+1):
        print '%d  ' % (k),
    print ' '
    print 'Survival (males):\t',
    for k in xrange(maxAge+1):
        print '%.3f  ' % (sMale[k]),
    print ' '
    print 'Survival (females):\t',
    for k in xrange(maxAge+1):
        print '%.3f  ' % (sFemale[k]),
    print ' '
    print 'Birth rate (females):\t',
    for k in xrange(maxAge+1):
        print '%.3f  ' % (bFemale[k]),
    print ' '
    print 'Mating prob (males):\t',
    for k in xrange(maxAge+1):
        print '%.3f  ' % (bMale[k]),
    print ' '
    if alphaMating:
        print '\nMating mode: %s %s' % (matingMode,'with alpha males')
    else:
        print '\nMating mode: %s %s' % (matingMode,'without alpha males')
    #Tables of alpha and culling parameters
    if alphaMating:
        print '\nAlpha Mating parameters'
        print '-----------------------'
        print 'prop. offspring: %.1f  ' % (alphaOff),
        print 'prop. of alphas: %.1f ' % (alphaProp)
        print 'Min age: %d   ' % (alphaMinAge),
        print 'Max age: %d  ' % (alphaMaxAge),
        print 'duration: %d' % (alphaDuration)
    if cullOpt:
        print '\nCulling parameters'
        print '------------------'
        print 'objective: %d\t\t' % (cullObjective)
        print 'target: %s INDIVIDUALS' % (cullTarget)
        print 'Mode: %s' % (cullMode)
        if cullMode=='THRESHOLD':
            print 'level: %.2f' % (cullLevel)
        if cullMode=='INTERVAL':
            print 'interval: %d ' % (cullInterval)
   
    # Globals   
    tempSize = [numIndivs]
   
    
    

    
    # Calculate the equilibrium age structure for the leslie matrix
    initAgeStructM = [0.0]*(maxAge+1)
    initAgeStructF = [0.0]*(maxAge+1)
    for x in xrange(maxAge+1):
#        if x == 0:
#            initAgeStructM[0]=1.0 #For males   = lx
#            initAgeStructF[0]=1.0 #For females
#        else:
        initAgeStructM[x]=ASinit[x]
        initAgeStructF[x]=ASinit[x]

    if myTracer:
        print "lx (m): ", initAgeStructM
        print "lx (f): ", initAgeStructF
    # Obtain T (generation interval) and Lambda (population growth rate) from the Leslie parameters for FEMALES!
    births = [0.0]*(maxAge+1)
    bxl = [0.0]*(maxAge+1)
    acum = 0.0
    acumbxl = 0.0
    for x in xrange(maxAge+1):
        births[x] = bFemale[x]*initAgeStructF[x]*(float(numIndivs)/2.0)
        acum += births[x]

        bxl[x] = bFemale[x]*initAgeStructF[x]
        acumbxl += bxl[x]
    # Remember: In this stage initAgeStruct is equal to lx
    expLambda = acum / ( initAgeStructF[0]*float(numIndivs))
    q = 0.0
    qx = 0.0
    V = [0.0]*(maxAge+1)
    for x in xrange(maxAge+1):
        qx = sum(float(bxl[i]) for i in xrange(x, maxAge+1)) / expLambda
        #print qx
        q += qx
        if initAgeStructF[x] > 0:
            V[x] = qx / initAgeStructF[x]
    #print V
    
    #print "Expected Lambda = {0:.3f}".format(expLambda)
    #print 'Expected Lambda = %.3f' % (expLambda)
    T = q
    #print "T (generation interval) = {0:.3f}".format(T)
    #print 'T (generation interval) = %.3f' % (T)

    
    


    #for x in xrange(maxAge+1):


    sumTot = 0.0
    for x in xrange(len(initAgeStructF)):
        sumTot += initAgeStructF[x]
    for x in xrange(len(initAgeStructF)):
        initAgeStructF[x] = initAgeStructF[x]/sumTot
    if myTracer: print "Initial age structure (Females): ", initAgeStructF
    acum = 0.0
    for x in xrange(len(initAgeStructF)):
        acum += initAgeStructF[x]
        initAgeStructF[x] = acum

    sumTot = 0.0
    for x in xrange(len(initAgeStructM)):
        sumTot += initAgeStructM[x]
    for x in xrange(len(initAgeStructM)):
        initAgeStructM[x] = initAgeStructM[x]/sumTot
    if myTracer: print "Initial age structure (Males): ", initAgeStructM
    acum = 0.0
    for x in xrange(len(initAgeStructM)):
        acum += initAgeStructM[x]
        initAgeStructM[x] = acum
        
    

    refer=[1] #counter for global index, dev purpouses

    # Global variables used to obtain averages across replicates and generations
    
    genN = numpy.zeros((gens, reps))
    genN0 = numpy.zeros((gens, reps))
    genLambda = numpy.zeros((gens, reps))
    genMales = numpy.zeros((gens, reps))
    genAlphas = numpy.zeros((gens, reps))
    genA = numpy.zeros((gens, reps))
    genAD = numpy.zeros((gens, reps))
    genHe = numpy.zeros((gens, reps))
    genHeL = numpy.zeros((gens, reps))
    genHo = numpy.zeros((gens, reps))
    genHoL = numpy.zeros((gens, reps))
    genNeFel = numpy.zeros((gens, reps))
    genNeH = numpy.zeros((gens, reps))
    genT = numpy.zeros((gens, reps))
    if cullOpt:
        genCulled = numpy.zeros((gens, reps))
    ASp = [0.0] * (maxAge+1)
    genAS = numpy.zeros((gens, (maxAge+1)))
    #genAS is a matrix to accumulate age props OVER REPLICATES
   
    globalNeH = [0.0]*(reps)
    globalK = [0.0]*(reps)
    globalVk = [0.0]*(reps)
    globalDads = [0.0]*(reps)



    try:
        os.mkdir(filePrefix + os.sep + 'sims')
    except OSError:
        pass #OK
    for gen in saveGens:
        try:
            #print gen
            os.mkdir(filePrefix + os.sep + 'sims' + os.sep + str(gen-burnin))
        except OSError:
            pass #OK

   
    progress = simuProgress("Running", reps*gens)
    cprog = 0
    for repli in xrange(reps):

        Vkmale = []
        fatherslist = []
        obsN0=cohortSize
        obsLambda=1.0
        refer=[1] #Now, reference number is for every replicate
        

        (loci, genPreOps, genInOps) = createGenome(numMSats, numSNPs)

        (pop, popPreOps, popInOps, oExpr) = createSinglePop(numIndivs, loci)

        pop.vars()['ageStruct'] = [0] * (maxAge + 1)
        pop.vars()['ageStructM'] = [0] * (maxAge + 1)
        pop.vars()['ageStructF'] = [0] * (maxAge + 1)
        pop.vars()['offspringF'] = [0] * (maxAge + 1)
        for x in xrange(maxAge+1):
            pop.vars()['offspringF'][x]=bFemale[x]   #Intialization to theoretical values
        pop.vars()['lx'] = [0] * (maxAge + 1)
        pop.vars()['lxM'] = [0] * (maxAge + 1)
        pop.vars()['lxF'] = [0] * (maxAge + 1)
        pop.vars()['removed'] = [0] * (maxAge + 1) #To store the number of removed individuals (deaths+culling) per age
        pop.vars()['T'] = T
        if cullOpt:
            pop.vars()['Culled'] = 0


        agess = []
        for x in xrange(numIndivs):
            rnd = uniform(0.0,1.0)
            for a in xrange(maxAge+1):
                if  rnd < initAgeStructM[a] and pop.individual(x).sex()==1:
                    agess.append(a)
                    break
                if  rnd < initAgeStructF[a] and pop.individual(x).sex()==2:
                    agess.append(a)
                    break

        pop.setIndInfo(agess, 'age')

        #initAlphas(pop)
        setRefs(pop, refer)
        mateOp = doAge(pop, maxAge, minMatingAge)


        #Reporting options and mix functions. They define the running sequences
        reportOps = [
            # INITIALIZACION
            pyOperator(func=initAlphas, at=[0], stage=PreMating),
            #pyOperator(fun=setRefs, param=refer, at=[0], stage=PreMating),
            # PREMATING
            stat(Fst=xrange(numSNPs + numMSats),
                alleleFreq=xrange(numSNPs + numMSats), heteroFreq=[0],
                expHetero=xrange(numSNPs + numMSats), numOfAlleles=xrange(numSNPs + numMSats),
                popSize=True, numOfMale=True, stage=PreMating),

            pyOperator(func=reportAge, step=1, stage=PreMating),
            pyOperator(func=acumStats, step=1, stage=PreMating),
            pyOperator(func=checkAlphas, step=1, stage=PreMating),
            infoExec('if age==0: cohort = sim.gen()', stage=PreMating), #recording when individuals born
            pyOperator(func=setFitness, step=1, stage=PreMating, infoFields=['fitness']),
            infoExec('ageFather = age', stage=PreMating),
            infoExec('ageMother = age', stage=PreMating),
            #pyOperator(func=incAge, stage=PreMating, begin=1), #Deprecated. Used when were a bug in infoExec
            infoExec('if sim.gen()>0: age += 1', stage=PreMating), #Everybody's birthday is just after mating

            #DURING MATING
            pySelector(loci=[0],func=allowSel, infoFields=['fitness','fitness']), #Obsolete
            pyTagger(func=obtainAgeParents,infoFields=['ageFather','ageMother']),
            pyTagger(func=obtainRefParents,infoFields=['father_ref','mother_ref']), #Re-do. To allow calculation of Vk

            #POSTMATING
            
            pyOperator(func=naturalDeath, step=1, stage=PostMating),
            
            pyOperator(func=reportParentAge, step=1, stage=PostMating),  #TRACkING
            pyOperator(func=checkDemo, step=1, stage=PostMating),
            
            pyOperator(func=setRefs, param=refer, step=1, stage=PostMating),

            stat(popSize=True, stage=PostMating), #Obtain the new pop size after births/deaths to evaluate if culling is needed.
            pyOperator(func=culling, begin=burnin, step=1, stage=PostMating), #Here??

           
            pyOperator(func=saver, param=oExpr, at=saveGens, stage=PostMating),
            pyExec('cprog += 1'),
            pyExec('progress.update(cprog)'),
            #pyOperator(func=zipper, at=saveGens, rep=0),


        ]
        sim = createSim(pop, 1, mateOp)
        evolveSim(sim, gens, numSNPs + numMSats,
            genPreOps, genInOps, popPreOps, popInOps, reportOps,
            oExpr, saveGens)
#        print Vkmale
        print fatherslist
        globalK[repli]= mean(Vkmale)
        globalVk[repli]= var(Vkmale)
        counter=0
#        fvk = open("Check-Vk.csv","a")
#        print "kmean: ", globalK[repli]
        for x in xrange(len(Vkmale)):
#            fvk.write('%d;;  %d \n' % (repli, Vkmale[x]) )
            if Vkmale[x]>0:
                counter += 1
        globalDads[repli]= float(counter)/float(len(Vkmale))
      # del pop
        del sim
        
     #End replicates
    for g in saveGens:
        zip(g)

    #zip(gen - 1)
    progress.done()
    outputStats()

    
   
    #raw_input("Press Enter to close this program")
    exit(0)
    

    # Program ends here
