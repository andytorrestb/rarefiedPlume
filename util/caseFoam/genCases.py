import casefoam

baseCase = '1d'
caseStructure = [['deltaT_02', 'deltaT_03', 'deltaT_04'],
                ['alpha_02', 'alpha_03', 'alpha_04']]

def update_timestep(deltaT):
    deltaT = 'deltaT          ' + deltaT
    return {
        'system/controlDict': {'#!stringManipulation':
                                {'deltaT          1.0E-05': '%s' %deltaT}}
    }

def update_maxCo(maxCo):
    maxCo = 'maxCo           ' + maxCo
    return {
        'system/controlDict': {'#!stringManipulation':
                                {'maxCo           0.5': '%s' %maxCo}}
    }

def update_alpha(alpha):
    maxCo = 'alpha                                     ' + alpha
    return {
        'constant/dsmcProperties': {'#!stringManipulation':
                                {'alpha                                     1.36': '%s' %maxCo}}
    }

caseData = {
    'deltaT_02': update_timestep('6.4e-4'),
    'deltaT_03': update_timestep('6.4e-5'),
    'deltaT_04': update_timestep('1.0e-5'),
    'alpha_02':  update_alpha('1.36'),
    'alpha_03':  update_alpha('1.26'),
    'alpha_04':  update_alpha('1.16')
}

# generate cases+
casefoam.mkCases(baseCase, caseStructure, caseData, hierarchy='tree',writeDir='Cases')