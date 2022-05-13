import casefoam

# Functions for updated variables in generated cases.
def update_nEquivalentParticles(nEqPar):
    nEqPar = 'nEquivalentParticles            ' + nEqPar
    return {
        'constant/dsmcProperties': {'#!stringManipulation':
                                {'nEquivalentParticles            2e16': '%s' %nEqPar}}
    }

def update_maxRadialWeightingFactor(maxRadWeight):
    maxRadWeight = 'maxRadialWeightingFactor            ' + maxRadWeight
    return {
        'constant/dsmcProperties': {'#!stringManipulation':
                                {'maxRadialWeightingFactor    20': '%s' %maxRadWeight}}
    }

# Define case structure
baseCase = '2d-wedge'
nEqPars = ['nEquivalentParticles_01', 'nEquivalentParticles_02',
                  'nEquivalentParticles_03', 'nEquivalentParticles_04',
                  'nEquivalentParticles_05', 'nEquivalentParticles_06']

maxRadWeights = ['maxRadialWeightingFactor_01', 'maxRadialWeightingFactor_02',
                 'maxRadialWeightingFactor_03', 'maxRadialWeightingFactor_04']

caseStructure = [nEqPars, maxRadWeights]

caseData = {
    'nEquivalentParticles_01': update_nEquivalentParticles('2e16'),
    'nEquivalentParticles_02': update_nEquivalentParticles('1e16'),
    'nEquivalentParticles_03': update_nEquivalentParticles('5e15'),
    'nEquivalentParticles_04': update_nEquivalentParticles('2.5e15'),
    'nEquivalentParticles_05': update_nEquivalentParticles('1.2e15'),
    'nEquivalentParticles_06': update_nEquivalentParticles('6e14'),
    'maxRadialWeightingFactor_01':  update_maxRadialWeightingFactor('125'),
    'maxRadialWeightingFactor_02':  update_maxRadialWeightingFactor('250'),
    'maxRadialWeightingFactor_03':  update_maxRadialWeightingFactor('500'),
    'maxRadialWeightingFactor_04':  update_maxRadialWeightingFactor('1000')
}

# Generate cases
casefoam.mkCases(baseCase, caseStructure, caseData, hierarchy='tree',writeDir='Cases')