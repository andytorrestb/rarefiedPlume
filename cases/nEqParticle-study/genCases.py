import casefoam

baseCase = '2d-wedge'

caseStructure = [['nEquivalentParticles_01', 'nEquivalentParticles_02', 'nEquivalentParticles_03'],
                ['vertical', 'diagonol', 'horizontal']]

def update_nEquivalentParticles(nEqPar):
    nEqPar = 'nEquivalentParticles            ' + nEqPar
    return {
        'system/controlDict': {'#!stringManipulation':
                                {'nEquivalentParticles            2e16': '%s' %nEqPar}}
    }

def update_sampleDict(direction):
    if (direction == 'horizontal'):
        print(direction)
        start = 'start ( 1.05 0.0 0.0 );'
        end = 'end ( 20.0 0.0 0.0 );'
        return {
            'system/sampleDict': {'#!stringManipulation':
                                    {'start ( 0.0 1.05 0.0 );': '%s' %start,
                                     'end ( 0.0 20.0 0.0 );': '%s' %end}}
        }
    elif (direction == 'diagonol'):
        print(direction)
        start = 'start ( 0.7425 0.7425 0.0 );'
        end = 'end ( 14.1421 14.1421 0.0 );'
        return {
            'system/sampleDict': {'#!stringManipulation':
                                    {'start ( 0.0 1.05 0.0 );': '%s' %start,
                                     'end ( 0.0 20.0 0.0 );': '%s' %end}}
        }
    elif (direction == 'vertical'):
        print(direction)
        start = 'start ( 0.0 1.05 0.0 );'
        end = 'end ( 0.0 20.0 0.0 );'
        return {
            'system/sampleDict': {'#!stringManipulation':
                                    {'start ( 0.0 1.05 0.0 );': '%s' %start,
                                     'end ( 0.0 20.0 0.0 );': '%s' %end}}
        }
    else:
        print("No Valid input")
        return

caseData = {
    'nEquivalentParticles_01': update_nEquivalentParticles('2e16'),
    'nEquivalentParticles_02': update_nEquivalentParticles('1e16'),
    'nEquivalentParticles_03': update_nEquivalentParticles('5e15'),
    'vertical':  update_sampleDict('vertical'),
    'diagonol':  update_sampleDict('diagonol'),
    'horizontal':  update_sampleDict('horizontal')
}

# generate cases+
casefoam.mkCases(baseCase, caseStructure, caseData, hierarchy='tree',writeDir='Cases')