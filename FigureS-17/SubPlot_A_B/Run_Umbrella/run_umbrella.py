import numpy as np
import pandas as pd
import sys
import os
import openmm.unit as unit
import openmm as mm
import openmm.app as app

class COMSeparationReporter():
    def __init__(self, file, reportInterval, g1_indices, g2_indices):
        self._out = open(file, 'w')
        self._reportInterval = reportInterval
        self.g1_indices = g1_indices
        self.g2_indices = g2_indices
        self.saved_data = None


    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        if self.saved_data is None:
            steps = self._reportInterval - simulation.currentStep % self._reportInterval - 1
        else:
            steps = 1
        return (steps, True, False, False, False)

    def report(self, simulation, state):
        positions = state.getPositions(asNumpy=True)

        # Calculate the masses and centers of mass for each group
        g1_masses = self.getMasses(system, self.g1_indices)
        g2_masses = self.getMasses(system, self.g2_indices)
        g1_com = self.getCenterOfMass(positions[self.g1_indices], g1_masses)
        g2_com = self.getCenterOfMass(positions[self.g2_indices], g2_masses)

        # Calculate the separation distance
        separation = np.linalg.norm(g1_com - g2_com)
        r = simulation.context.getParameter('r0')
        self._out.write(f'{simulation.currentStep}\t{r}\t{separation}\n')
        self._out.flush()

    def getMasses(self, system, indices):
        masses = [system.getParticleMass(i).value_in_unit(unit.dalton) for i in indices]
        return np.array(masses)

    def getCenterOfMass(self, positions, masses):
        total_mass = np.sum(masses)
        com = np.sum(positions * masses[:, None], axis=0) / total_mass
        return com

#export OPENMM_PLUGIN_DIR=~/anaconda3/envs/omm/lib/plugins
temperature = 300
box_a = 300
box_b = 300
box_c = 300
output_dcd = 'box.dcd'
output_interval = 100
steps = 10000

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
timestep = 1*unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(start_temperature*unit.kelvin, collision, timestep)
platform_name = 'CPU'
platform = mm.Platform.getPlatformByName(platform_name)
simulation = app.Simulation(top, system, integrator, platform)
simulation.context.setPositions(init_coord)
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(start_temperature*unit.kelvin)
n_iterations = 10
state_data_reporter = app.StateDataReporter('pot_eng_nneg.txt', output_interval, step=True, potentialEnergy=True)
comrepo = COMSeparationReporter('local_sep_nneg.txt', output_interval,[0],[1])  # Replace 100 with desired reporting interval

for temperature_i in np.linspace(start_temperature, temperature, n_iterations):
	integrator.setTemperature(temperature_i*unit.kelvin)
	simulation.step(1000)
simulation.step(steps/10)
simulation.reporters.append(comrepo)
simulation.reporters.append(state_data_reporter)
r_range = np.linspace(4,0.3,10)
for i in r_range:
    simulation.context.setParameter('r0',i)
    simulation.step(steps)

