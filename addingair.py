#!/usr/bin/env python
# encoding: utf-8
r"""
Two-dimensional variable-coefficient acoustics
==============================================

Solve the variable-coefficient acoustics equations in 2D:

.. math:: 
    p_t + K(x,y) (u_x + v_y) & = 0 \\ 
    u_t + p_x / \rho(x,y) & = 0 \\
    v_t + p_y / \rho(x,y) & = 0.

Here p is the pressure, (u,v) is the velocity, :math:`K(x,y)` is the bulk modulus,
and :math:`\rho(x,y)` is the density.

This example shows how to solve a problem with variable coefficients.
The left and right halves of the domain consist of different materials.
"""
import matplotlib.pyplot as plt
import numpy as np

def setup(kernel_language='Fortran', use_petsc=False, outdir='./_output', 
          solver_type='classic', time_integrator='SSP104', lim_type=2, 
          disable_output=False, num_cells=(1000,52)):
    """
    Example python script for solving the 2d acoustics equations.
    """
    from clawpack import riemann
    
    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='classic':
        solver=pyclaw.ClawSolver2D(riemann.vc_acoustics_2D)
        solver.dimensional_split=False
        solver.limiters = pyclaw.limiters.tvd.MC
    elif solver_type=='sharpclaw':
        solver=pyclaw.SharpClawSolver2D(riemann.vc_acoustics_2D)
        solver.time_integrator=time_integrator
        if time_integrator=='SSPLMMk2':
            solver.lmm_steps = 3
            solver.cfl_max = 0.25
            solver.cfl_desired = 0.24

    solver.bc_lower[0]=pyclaw.BC.extrap
    solver.bc_upper[0]=pyclaw.BC.extrap
    solver.bc_lower[1]=pyclaw.BC.wall
    solver.bc_upper[1]=pyclaw.BC.wall
    solver.aux_bc_lower[0]=pyclaw.BC.extrap
    solver.aux_bc_upper[0]=pyclaw.BC.extrap
    solver.aux_bc_lower[1]=pyclaw.BC.extrap
    solver.aux_bc_upper[1]=pyclaw.BC.extrap

    x = pyclaw.Dimension(0.,1000000.,num_cells[0],name='x')   #Scale of Graph    0.,10000.
    y = pyclaw.Dimension(-5200.,0.0,num_cells[1],name='y')  #-7000.,100.
    domain = pyclaw.Domain([x,y])

    num_eqn = 3
    num_aux = 6 # density, sound speed
    state = pyclaw.State(domain,num_eqn,num_aux)

    grid = state.grid
    X, Y = grid.p_centers

    rho_air = 1.225 #Density of Air  
    rho_1 = 1025.5 # Density in layer 1 from 0-200 m
    rho_2 = 1026.5 # Density in layer 2 from 2-1000 m
    rho_3 = 1028.  # Density in layer 3 from 1000-000 m
    rho_4 = 1028.  # Density in layer 4 from 4000-5300 m

    import numpy as np
    import csv

    # f = open('tempsalinity.csv')
    # csv_f = csv.reader(f)
 
    reader = csv.reader(open("tempsalinity.csv","rb"),delimiter=',')
    result = np.array(list(reader)).astype('float')
    depth= result[0]
    salinity = result [1]
    temperature = result [2]
    speed_arr=np.ndarray((51,), float)
  
    for counter in range (0,51):
        T = temperature[counter]
        S = salinity[counter]
        D = depth[counter]
        speed = 1448.96 + 4.591*T-5.304e-2*(T**2)+2.374e-4*(T**3)+1.340*(S-35)+1.630e-2*D+1.675e-7*(D**2)-1.025e-2*T*(S-35)-7.139e-13*T*(D**3)
        speed_arr[counter]=speed
    y_values = y.centers[::-1]
    yinterp = np.flipud(np.interp(-y_values, depth, speed_arr))
    y_values = np.flipud(y_values)
    import math
    for j in xrange(num_cells[1]):
        for i in xrange(num_cells[0]):
            if y_values[j] > 0:
                state.aux[1, i, j] = 343.
            else:
                state.aux[1, i, j] = yinterp[j]
            if y_values[j] > 0:
                 state.aux[0, i, j] = rho_air
            elif y_values[j] < 0 and y_values[j] >= -200:
                 state.aux[0, i, j] = rho_1
            elif y_values[j] < -200 and y_values[j] >= -1000:
                state.aux[0, i, j] = rho_2
            elif y_values[j] < -1000 and y_values[j] >= -4000:
                state.aux[0, i, j] = rho_3
            elif y_values[j] < -4000 and y_values[j] >= -5400:
                state.aux[0, i, j] = rho_4

    # Set initial condition
    x0 =5000.0 ; y0 = -1000;
    r = np.sqrt((X-x0)**2 + (Y-y0)**2)
    width = 500.0; rad = 500
    amplitude = 1.0
    #ADJUSTING THE PRESSURE
    state.q[0, :, :] = amplitude * np.exp(-(X - x0)**2 / width**2)
    #ADJUSTING THE PRESSURE
    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.
    claw = pyclaw.Controller()
    claw.keep_copy = True
    if disable_output:
        claw.output_format = None 
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.tfinal = 1000
    claw.num_output_times = 100
    claw.write_aux_init = True
    claw.setplot = "./setplot.py"
    if use_petsc:
        claw.output_options = {'format':'binary'}

    return claw


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup, "./setplot.py")