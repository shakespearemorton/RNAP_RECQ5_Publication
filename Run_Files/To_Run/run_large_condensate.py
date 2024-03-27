import numpy as np
import pandas as pd
import sys
import os

import simtk.unit as unit
import simtk.openmm as mm
import simtk.openmm.app as app

temperature = 300
box_a = 500
box_b = 500
box_c = 500
output_dcd = 'slab.dcd'
output_interval = 100000
steps = 100000000

system_xml = 'system.xml'
with open(system_xml, 'r') as f:
    system = mm.XmlSerializer.deserialize(f.read())
box_vec_a = np.array([box_a, 0, 0])*unit.nanometer
box_vec_b = np.array([0, box_b, 0])*unit.nanometer
box_vec_c = np.array([0, 0, box_c])*unit.nanometer
system.setDefaultPeriodicBoxVectors(box_vec_a, box_vec_b, box_vec_c)
print('System box vectors are:')
print(system.getDefaultPeriodicBoxVectors())

top = app.PDBFile('start.pdb').getTopology()

npt_final_state_xml = 'NPT_final_state.xml'
with open(npt_final_state_xml, 'r') as f:
    npt_final_state = mm.XmlSerializer.deserialize(f.read())
init_coord = np.array(npt_final_state.getPositions().value_in_unit(unit.nanometer))
# move the geometric center of atom coordinates to box center
init_coord -= np.mean(init_coord, axis=0)
init_coord += 0.5*np.array([box_a, box_b, box_c])

start_temperature = 150
collision = 1/unit.picosecond
timestep = 10*unit.femtosecond
integrator = mm.NoseHooverIntegrator(start_temperature*unit.kelvin, collision, timestep)
platform_name = 'CUDA'
platform = mm.Platform.getPlatformByName(platform_name)
properties = {'Precision': 'mixed'}
simulation = app.Simulation(top, system, integrator, platform, properties)
simulation.context.setPositions(init_coord)
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(start_temperature*unit.kelvin)
n_iterations = 100
for temperature_i in np.linspace(start_temperature, temperature, n_iterations):
    integrator.setTemperature(temperature_i*unit.kelvin)
    simulation.step(1000)

dcd_reporter = app.DCDReporter(output_dcd, output_interval, enforcePeriodicBox=True)
state_data_reporter = app.StateDataReporter(sys.stdout, output_interval, step=True, time=True, potentialEnergy=True,
                                            kineticEnergy=True, totalEnergy=True, temperature=True, speed=True)
simulation.reporters.append(dcd_reporter)
simulation.reporters.append(state_data_reporter)
simulation.step(steps)
with open('equilibrated_config.pdb','w') as f:
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    top.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
    app.PDBFile.writeFile(top, state.getPositions(), f)
