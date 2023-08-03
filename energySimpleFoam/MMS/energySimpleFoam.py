#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 17:38:18 2021

@author: Bruno Ramoa
@affiliation: Institute for Polymers and Composites, University of Minho, Portugal
"""

# Script for MMS in OpenFOAM
import sympy as sym
from sympy import sin, cos, exp, pi, sqrt, tanh
import pyMMSFoam as mms
from pyMMSFoam import x,y,z,t
import sys

u = 1 + 0.5*y*tanh(2*pi*(x*x + y*y))
v = 1 - 0.5*x*tanh(2*pi*(x*x + y*y))
w = 0
p = 0.5 + tanh(2*pi*x*x) 
T = 1 + 0.2*x*(tanh(2*pi*x*y) + 1.5)

U = sym.Matrix([u,v,w])

# Momentum balance equation
nu = 0.5

R = nu*(mms.grad(U) + mms.grad(U).T)

S = mms.div(U*U.T) - mms.div(R) + mms.grad(p)

# Energy balance equation
alpha = 1

S_T = mms.div(U*T) - mms.laplacian(T, alpha) 

## Test
#print(S_T)
#print("\n")
#print(mms.laplacian(T, alphaEff))

# # Uncomment what is needed

# # Generate fvOptions
orig_stdout = sys.stdout
file = open('system/fvOptions', 'w')
sys.stdout = file
mms.generateFvOptions(S,"momentumSource", "U")
mms.generateFvOptions(S_T, "sourceTerm", "T")

sys.stdout  = orig_stdout
file.close()

# # # Generate boundary conditions
# # Velocity
orig_stdout = sys.stdout
file = open('0.orig/include/U_bc', 'w')
sys.stdout = file
mms.generateDirichletBoundaries(U, "U")
mms.generateNeumannBoundaries(U, "U")

sys.stdout  = orig_stdout
file.close()

# # Pressure
orig_stdout = sys.stdout
file = open('0.orig/include/p_bc', 'w')
sys.stdout = file
mms.generateDirichletBoundaries(p, "p")
mms.generateNeumannBoundaries(p, "p")

sys.stdout  = orig_stdout
file.close()

# # Temperature
orig_stdout = sys.stdout
file = open('0.orig/include/T_bc', 'w')
sys.stdout  = file
mms.generateDirichletBoundaries(T, "T")
mms.generateNeumannBoundaries(T, "T")
sys.stdout  = orig_stdout
file.close()

# # Generate functionObjects
orig_stdout = sys.stdout
file = open('system/funcobj', 'w')
sys.stdout = file
mms.generateFunctionObject(U, "U")
mms.generateFunctionObject(p, "p")
mms.generateFunctionObject(T, "T")

sys.stdout  = orig_stdout
file.close()
