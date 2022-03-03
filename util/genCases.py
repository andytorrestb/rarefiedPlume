import casefoam

baseCase = '1d'
caseStructure = [['deltaT_02', 'deltaT_03', 'deltaT_04'],
                ['endTime_02', 'endTime_03', 'endTime_04']]

def update_timestep(deltaT):
    deltaT = 'deltaT          ' + deltaT
    return {
        'system/controlDict': {'#!stringManipulation':
                                {'deltaT          1.0E-05': '%s' %deltaT}}
    }


def update_endTime(endTime):
    maxCo = 'endTime         ' + endTime
    return {
        'system/controlDict': {'#!stringManipulation':
                                {'endTime         2.0E-01': '%s' %maxCo}}
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
    'endTime_02':  update_endTime('1.0E-01'),
    'endTime_03':  update_endTime('2.0E-01'),
    'endTime_04':  update_endTime('3.0E-01')
}

# generate cases+
casefoam.mkCases(baseCase, caseStructure, caseData, hierarchy='tree',writeDir='Cases')