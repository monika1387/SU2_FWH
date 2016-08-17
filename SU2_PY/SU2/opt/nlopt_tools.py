#!/usr/bin/env python 

## \file scipy_tools.py
#  \brief tools for interfacing with scipy
#  \author T. Lukaczyk, F. Palacios
#  \version 4.2.0 "Cardinal"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2016 SU2, the open-source CFD code.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------
#  Imports
# -------------------------------------------------------------------

import os, sys, shutil, copy

from .. import eval as su2eval
from numpy import array, zeros
from numpy.linalg import norm
sys.path.append('%s/externals/nlopt/swig'%os.environ['SU2_HOME'])
import nlopt

# -------------------------------------------------------------------
#  Scipy SLSQP
# -------------------------------------------------------------------

def scipy_slsqp(project,optmethod='SLSQP',x0=None,xu=None, xl=None,its=100,accu=1e-10,grads=True):
    """ result = scipy_slsqp(project,x0=[],xb=[],its=100,accu=1e-10)
    
        Runs the Scipy implementation of SLSQP with 
        an SU2 project
        
        Inputs:
            project - an SU2 project
            x0      - optional, initial guess
            xb      - optional, design variable bounds
            its     - max outer iterations, default 100
            accu    - accuracy, default 1e-10
        
        Outputs:
           result - the outputs from scipy.fmin_slsqp
    """

    # import scipy optimizer
    from scipy.optimize import fmin_slsqp

    # handle input cases
    if x0 is None: x0 = []
    if xl is None: xl = [-inf]*n_dv
    if xu is None: xu = [inf]*n_dv

    # number of design variables
    dv_size = project.config['DEFINITION_DV']['SIZE']
    n_dv = sum( dv_size)
    project.n_dv = n_dv
    
    # Initial guess
    if not x0: x0 = [0.0]*n_dv
    
    # prescale x0
    dv_scales = project.config['DEFINITION_DV']['SCALE']
    k = 0
    for i, dv_scl in enumerate(dv_scales):
        for j in range(dv_size[i]):
            x0[k] =x0[k]/dv_scl;
            k = k + 1

    # scale accuracy
    obj = project.config['OPT_OBJECTIVE']
    obj_scale = []
    for this_obj in obj.keys():
        obj_scale = obj_scale + [obj[this_obj]['SCALE']]
    
    accu = accu*obj_scale[0]

    # scale accuracy
    eps = 1.0e-04

    n_eqcons  = len( project.config['OPT_CONSTRAINT']['EQUALITY'])
    n_ieqcons = len( project.config['OPT_CONSTRAINT']['INEQUALITY'])

    if optmethod == 'SLSQP':
        opt = nlopt.opt(nlopt.LD_SLSQP, n_dv)
    elif optmethod == 'MMA':
        opt = nlopt.opt(nlopt.LD_MMA, n_dv)
    elif optmethod == 'AUGLAG_MMA':
        opt = nlopt.opt(nlopt.AUGLAG_EQ, n_dv)
        opt_local = nlopt.opt(nlopt.LD_MMA, n_dv)
        opt_local.set_xtol_rel(1e-01)
        opt.set_local_optimizer(opt_local)
    else:
        sys.exit('Optimizer %s not found'%optmethod)

    opt.set_min_objective(lambda x, grad: obj_f(x, grad, project))
    opt.set_xtol_rel(1e-04)
    opt.set_maxeval(its)
    if n_ieqcons > 0:
        opt.add_inequality_mconstraint(lambda cons, x, grad: con_cieq(cons, x, grad, project), [1e-08]*n_ieqcons)
    if n_eqcons > 0:
        opt.add_equality_mconstraint(lambda cons, x, grad: con_ceq(cons, x, grad, project), [1e-08]*n_eqcons)

    outputs = opt.optimize(x0)

    # Done
    result = opt.last_optimize_result()
    print opt.get_ftol_rel()
    return outputs
    
    
def obj_f(x, grad, project):
    """ obj = obj_f(x,project)
        
        Objective Function
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         minimize f(x), list[nobj]
        scipy_slsqp: minimize f(x), float
    """

    obj_list = project.obj_f(x)
    obj = 0
    for this_obj in obj_list:
        obj = obj+this_obj
    if grad.size > 0:
        grad_out = obj_df(x, project)
        for i in range(0, grad.size):
            grad[i] = grad_out[i]

    return obj

def obj_df(x,project):
    """ dobj = obj_df(x,project)
        
        Objective Function Gradients
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         df(x), list[nobj x dim]
        scipy_slsqp: df(x), ndarray[dim]
    """    
    
    dobj_list = project.obj_df(x)
    dobj=[0.0]*len(dobj_list[0])
    
    for this_dobj in dobj_list:
        idv=0
        for this_dv_dobj in this_dobj:
            dobj[idv] = dobj[idv]+this_dv_dobj;
            idv+=1
    dobj = array( dobj )
    return dobj

def con_ceq(cons, x,grad, project):
    """ cons = con_ceq(x,project)
        
        Equality Constraint Functions
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         ceq(x) = 0.0, list[nceq]
        scipy_slsqp: ceq(x) = 0.0, ndarray[nceq]
    """
    
    cons_out = project.con_ceq(x)
    n_eqcons  = len(project.config['OPT_CONSTRAINT']['EQUALITY'])

    for i in range(0, neqcons):
        cons[i] = cons_out[i]

    if grad.size > 0:
        grad_out = con_dceq(x, project)
        for i in range(0,n_eqcons):
            for j in range(0, grad[i].size):
                grad[i][j] = grad_out[i][j]

def con_dceq(x, project):
    """ dcons = con_dceq(x,project)
        
        Equality Constraint Gradients
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         dceq(x), list[nceq x dim]
        scipy_slsqp: dceq(x), ndarray[nceq x dim]
    """
    
    n_eqcons  = len(project.config['OPT_CONSTRAINT']['EQUALITY'])
    dcon = zeros([n_eqcons, len(x)])
    deqcons = project.con_dceq(x)

    for i_eqcons, eqcons in enumerate(deqcons):
        dcon[i_eqcons] = array(eqcons)

    return dcon

def con_cieq(cons, x,grad, project):
    """ cons = con_cieq(x,project)
        
        Inequality Constraints
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         cieq(x) < 0.0, list[ncieq]
        scipy_slsqp: cieq(x) > 0.0, ndarray[ncieq]
    """
    
    cons_out = project.con_cieq(x)
    n_ieqcons  = len(project.config['OPT_CONSTRAINT']['INEQUALITY'])

    for i in range(0, n_ieqcons):
        cons[i] = cons_out[i]

    if grad.size > 0:
        grad_out = con_dcieq(x, project)
        for i in range(0,n_ieqcons):
            for j in range(0, grad[i].size):
                grad[i][j] = grad_out[i][j]

def con_dcieq(x,project):
    """ dcons = con_dcieq(x,project)
        
        Inequality Constraint Gradients
        SU2 Project interface to scipy.fmin_slsqp
        
        su2:         dcieq(x), list[ncieq x dim]
        scipy_slsqp: dcieq(x), ndarray[ncieq x dim]
    """
    
    n_ieqcons  = len(project.config['OPT_CONSTRAINT']['INEQUALITY'])
    dcon = zeros([n_ieqcons, len(x)])
    dieqcons = project.con_dcieq(x)

    for i_ieqcons, ieqcons in enumerate(dieqcons):
        dcon[i_ieqcons] = array(ieqcons)

    return dcon
