import numpy as np
import sys
import os

import openmm.unit as unit
import openmm as mm
import openmm.app as app
#export OPENMM_PLUGIN_DIR=~/anaconda3/envs/omm/lib/plugins
temperature = 300
box_a = 100
box_b = 100
box_c = 100
output_dcd = 'rg2.dcd'
output_interval = 10000
steps = 1000000000

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
print(top)

init_coord = np.array(app.PDBFile('start.pdb').getPositions().value_in_unit(unit.nanometer))
init_coord -= np.mean(init_coord, axis=0)
init_coord += 0.5*np.array([box_a, box_b, box_c])

start_temperature = 150
collision = 1/unit.picosecond
timestep = 10*unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(start_temperature*unit.kelvin, collision, timestep)
platform_name = 'CUDA'
platform = mm.Platform.getPlatformByName(platform_name)
simulation = app.Simulation(top, system, integrator, platform)
simulation.context.setPositions(init_coord)
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(start_temperature*unit.kelvin)
n_iterations = 10
dcd_reporter = app.DCDReporter(output_dcd, output_interval, enforcePeriodicBox=True)
state_data_reporter = app.StateDataReporter(sys.stdout, output_interval, step=True, time=True, potentialEnergy=True,
                                            kineticEnergy=True, totalEnergy=True, temperature=True, speed=True)
for temperature_i in np.linspace(start_temperature, temperature, n_iterations):
	integrator.setTemperature(temperature_i*unit.kelvin)
	simulation.step(1000)


simulation.reporters.append(dcd_reporter)
simulation.reporters.append(state_data_reporter)
simulation.step(steps)
