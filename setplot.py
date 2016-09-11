#!/usr/bin/env python

from __future__ import print_function

import numpy
import matplotlib.pyplot as plt

import cmocean

def mark_interface(current_data):
    import matplotlib.pyplot as plt

    plt.plot((0.,100000.),(0.,0.),'-.k',linewidth=2)
    # plt.plot((0.,100000.),(-200.,-200.),'-k',linewidth=2)
    # plt.plot((0.,100000.),(-1000.,-1000.),'-k',linewidth=2)
    # plt.plot((0.,100000.),(-4000.,-4000.),'-k',linewidth=2)
    # plt.plot((0.,100000.),(-6000.,-6000.),'-k',linewidth=2)

def density(cd):
    return cd.aux[0, :, :] - 1000.0

def sound_speed(cd):
    return cd.aux[1, :, :]

def log_pressure(cd):
    return numpy.log(cd.q[0, :, :] + 1.0)

def setplot(plotdata):
    """ 
    Plot solution using VisClaw.

    This example shows how to mark an internal boundary on a 2D plot.
    """ 

    y_axes = [-5200., 0.0]
    x_axes_init = [0.0, 1000.0]
    rho_bounds = [25, 28]
    c_bounds = [1480, 1525]

    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    
    # Figure for pressure
    plotfigure = plotdata.new_plotfigure(name='Pressure', figno=0)
    plotfigure.kwargs = {'figsize':(16, 6)}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Pressure'
    plotaxes.scaled = False      # so aspect ratio is 1
    # plotaxes.xlimits = [0.0, 10000]
    plotaxes.ylimits = y_axes
    plotaxes.afteraxes = mark_interface

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = cmocean.cm.haline
    plotitem.add_colorbar = True
    plotitem.pcolor_cmin = 0
    plotitem.pcolor_cmax = 0.1

    # Figure for sound speed and density
    plotfigure = plotdata.new_plotfigure(name="Density and Sound Speed", 
                                         figno=2)
    plotfigure.kwargs = {'figsize':(16, 6)}

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = "subplot(1, 2, 1)"
    plotaxes.title = "Density ($kg/m^3$)"
    plotaxes.scaled = False
    plotaxes.xlimits = x_axes_init
    plotaxes.ylimits = y_axes
    plotaxes.afteraxes = mark_interface

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = density
    plotitem.pcolor_cmap = cmocean.cm.dense
    plotitem.add_colorbar = True
    plotitem.pcolor_cmin = rho_bounds[0]
    plotitem.pcolor_cmax = rho_bounds[1]
    plotitem.colorbar_label = r"$\rho - 1000~(kg/cm^3)$"

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = "subplot(1, 2, 2)"
    plotaxes.title = "Speed of Sound ($m/s$)"
    plotaxes.scaled = False
    plotaxes.xlimits = x_axes_init
    plotaxes.ylimits = y_axes
    plotaxes.afteraxes = mark_interface

    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.colorbar_label = "$c~(m/s)$"
    plotitem.plot_var = sound_speed
    plotitem.pcolor_cmap = cmocean.cm.haline
    plotitem.add_colorbar = True
    plotitem.pcolor_cmin = c_bounds[0]
    plotitem.pcolor_cmax = c_bounds[1]

    
    return plotdata