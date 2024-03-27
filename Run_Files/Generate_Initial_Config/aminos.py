import openmm as mm
import openmm.app as app
import openmm.unit as unit

amino_acids = {
    1: 'MET',
    2: 'GLY',
    3: 'LYS',
    4: 'THR',
    5: 'ARG',
    6: 'ALA',
    7: 'ASP',
    8: 'GLU',
    9: 'TYR',
    10: 'VAL',
    11: 'LEU',
    12: 'GLN',
    13: 'TRP',
    14: 'PHE',
    15: 'SER',
    16: 'HIS',
    17: 'ASN',
    18: 'PRO',
    19: 'CYS',
    20: 'ILE',
    21: 'SEP',
    22: 'SRI'
}

mass = {
    'MET': 131.199997,
    'GLY': 57.049999,
    'LYS': 128.199997,
    'THR': 101.099998,
    'ARG': 156.199997,
    'ALA': 71.080002,
    'ASP': 115.099998,
    'GLU': 129.100006,
    'TYR': 163.199997,
    'VAL': 99.07,
    'LEU': 113.199997,
    'GLN': 128.100006,
    'TRP': 186.199997,
    'PHE': 147.199997,
    'SER': 87.080002,
    'HIS': 137.100006,
    'ASN': 114.099998,
    'PRO': 97.120003,
    'CYS': 103.099998,
    'ILE': 113.199997,
    'SEP': 165.03,
    'SRI': 1
}

size = {
    'MET': 6.18,
    'GLY': 4.5,
    'LYS': 6.36,
    'THR': 5.62,
    'ARG': 6.56,
    'ALA': 5.04,
    'ASP': 5.58,
    'GLU': 5.92,
    'TYR': 6.46,
    'VAL': 5.86,
    'LEU': 6.18,
    'GLN': 6.02,
    'TRP': 6.78,
    'PHE': 6.36,
    'SER': 5.18,
    'HIS': 6.08,
    'ASN': 5.68,
    'PRO': 5.56,
    'CYS': 5.48,
    'ILE': 6.18,
    'SEP': 6.36,
    'SRI': 5.8
}

hydropathy = {
    'MET': 0.676471,
    'GLY': 0.57353,
    'LYS': 0.382354,
    'THR': 0.588236,
    'ARG': 0.558824,
    'ALA': 0.602942,
    'ASP': 0.294119,
    'GLU': 0.558824,
    'TYR': 0.897059,
    'VAL': 0.664707,
    'LEU': 0.720589,
    'GLN': 0,
    'TRP': 1,
    'PHE': 0.82353,
    'SER': 0.588236,
    'HIS': 0.764707,
    'ASN': 0.588236,
    'PRO': 0.758824,
    'CYS': 0.64706,
    'ILE': 0.705883,
    'SEP': 0.162,
    'SRI': 0.602
}

charge = {
    'MET': 0,
    'GLY': 0,
    'LYS': 1,
    'THR': 0,
    'ARG': 1,
    'ALA': 0,
    'ASP': -1,
    'GLU': -1,
    'TYR': 0,
    'VAL': 0,
    'LEU': 0,
    'GLN': 0,
    'TRP': 0,
    'PHE': 0,
    'SER': 0,
    'HIS': 0,
    'ASN': 0,
    'PRO': 0,
    'CYS': 0,
    'ILE': 0,
    'SEP': -2,
    'SRI': 4
}


reverse = {v: k for k, v in amino_acids.items()}
_kcal_to_kj = 4.184
epsilon=0.2*_kcal_to_kj
mu=1,
delta=0.08
use_pbc = True
ldby=1*unit.nanometer
dielectric_water=80.0
NA = unit.AVOGADRO_CONSTANT_NA # Avogadro constant
kB = unit.BOLTZMANN_CONSTANT_kB  # Boltzmann constant
EC = 1.602176634e-19*unit.coulomb # elementary charge
VEP = 8.8541878128e-12*unit.farad/unit.meter # vacuum electric permittivity
GAS_CONST = 1.0*unit.BOLTZMANN_CONSTANT_kB*unit.AVOGADRO_CONSTANT_NA # gas constant
