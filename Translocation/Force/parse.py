import subprocess
import numpy as np

def readParameters(filename='parameters'):

    parameterFile = open(filename, 'r')

    # Remove lines starting with '#'
    lines = [x.split(' ') for x in parameterFile.read().split('\n') if len(x)>0 and x[0] != '#']

    # Remove '=' and spaces
    lines = [[y.strip() for y in x if y!='=' and y!=''] for x in lines]

    parameters = {}

    for line in lines:
        if len(line) < 3:
            parameters[line[0]] = line[1]
        else:
            parameters[line[0]] = line[1:]
   
    parameterFile.close()

    return parameters
