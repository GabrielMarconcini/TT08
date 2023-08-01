#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 17:38:18 2021

@author: Bruno Ramoa
@affiliation: Institute for Polymers and Composites, University of Minho, Portugal
"""

# Script for MMS in OpenFOAM
import sympy as sym
from sympy import sin, cos, exp, pi, sqrt
import pyMMSFoam as mms
from pyMMSFoam import x,y,z,t
import sys

# Manufactured solutions
u = 1 + 0.5*sin(pi*x)*cos(pi*y)
v = 1 - 0.5*cos(pi*x)*sin(pi*y)
w = 0
p = 1 + cos(pi*x)
Y1 = sin(pi*(x + y))**2 # + sin(pi*y)**2
Y2 = cos(pi*(x + y))**2 # + cos(pi*y)**2

U = sym.Matrix([u,v,w])

# Momentum balance equation
nu = 0.01
R = nu*(mms.grad(U) + mms.grad(U).T)
S_U  = mms.div(U*U.T) - mms.div(R) + mms.grad(p)


# Species transport equation
D  = 0.01
S_Y1 = mms.div(Y1*U) - mms.laplacian(Y1, D)
S_Y2 = mms.div(Y2*U) - mms.laplacian(Y2, D)



# # Uncomment what is needed

# # Generate fvOptions
orig_stdout = sys.stdout
file        = open('system/fvOptions', 'w')
sys.stdout  = file
mms.generateFvOptions(S_U, "momentumSource", "U")
sys.stdout  = orig_stdout
file.close()

orig_stdout = sys.stdout
file        = open('system/firstSpecieSource', 'w')
sys.stdout  = file
mms.generateFvOptions(S_Y1, "firstSpecieSource", "Y1")
sys.stdout  = orig_stdout
file.close()

orig_stdout = sys.stdout
file        = open('system/secondSpecieSource', 'w')
sys.stdout  = file
mms.generateFvOptions(S_Y2, "secondSpecieSource", "Y2")
sys.stdout  = orig_stdout
file.close()

# # Generate boundary conditions
# # Velocity
orig_stdout = sys.stdout
file        = open('0.orig/include/U_BC', 'w')
sys.stdout  = file
mms.generateDirichletBoundaries(U, "U")
mms.generateNeumannBoundaries(U, "U")
sys.stdout  = orig_stdout
file.close()

# # Pressure
orig_stdout = sys.stdout
file        = open('0.orig/include/p_BC', 'w')
sys.stdout  = file
mms.generateDirichletBoundaries(p, "p")
mms.generateNeumannBoundaries(p, "p")
sys.stdout  = orig_stdout
file.close()

# # Species
orig_stdout = sys.stdout
file        = open('0.orig/include/Y1_BC', 'w')
sys.stdout  = file
mms.generateDirichletBoundaries(Y1, "Y1")
mms.generateNeumannBoundaries(Y1, "Y1")
sys.stdout  = orig_stdout
file.close()

orig_stdout = sys.stdout
file        = open('0.orig/include/Y2_BC', 'w')
sys.stdout  = file
mms.generateDirichletBoundaries(Y2, "Y2")
mms.generateNeumannBoundaries(Y2, "Y2")
sys.stdout  = orig_stdout
file.close()


# # Generate functionObjects
orig_stdout = sys.stdout
file        = open('system/functionObj', 'w')
sys.stdout  = file
mms.generateFunctionObject(U, "U")
mms.generateFunctionObject(p, "p")
mms.generateFunctionObject(Y1, "Y1")
mms.generateFunctionObject(Y2, "Y2")
sys.stdout  = orig_stdout
file.close()

