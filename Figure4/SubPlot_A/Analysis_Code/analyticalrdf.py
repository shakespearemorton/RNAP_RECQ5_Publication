from numpy import sqrt, pi, zeros, ones, arctan
from random import random, seed
from time import time

seed(time())

def SphereCutVol(Rs, A, B):
    Root = sqrt(Rs**2 - A**2 - B**2)
    Vcut = (1/6)*Rs**3*(pi - 2*arctan(A*B/(Rs*Root)))
    Vcut += (1/2) * (arctan(A/Root) - pi/2)*(Rs**2*B - 1/3*B**3)
    Vcut += (1/2) * (arctan(B/Root) - pi/2)*(Rs**2*A - 1/3*A**3)
    Vcut += (1/3) * A * B * Root
    return Vcut

def OctVolume(Rs, xb, yb, zb):
    if xb**2 + yb**2 + zb**2 < Rs**2:
        return xb * yb * zb
    VOctant = (1/8) * (4/3) * pi * Rs**3
    for B in [xb, yb, zb]:
        if B < Rs:
            VOctant -= (pi/4)*(2/3*Rs**3 - B*Rs**2 + 1/3*B**3)
    for (a, b) in [(xb, yb), (xb, zb), (yb, zb)]:
        if a**2 + b**2 < Rs**2:
            VOctant += SphereCutVol(Rs, a, b)
    return VOctant

def SphereVolume(Rs, BoxBounds):
    [Xmin, Xmax, Ymin, Ymax, Zmin, Zmax] = BoxBounds
    VSphere = 0
    for xb in [Xmin, Xmax]:
        for yb in [Ymin, Ymax]:
            for zb in [Zmin, Zmax]:
                VSphere += OctVolume(Rs, abs(xb), abs(yb), abs(zb))
    return VSphere

def ShellVolume(Rmin, Rmax, BoxBounds):
    Rmin = max(Rmin, 0)
    InnerShell = SphereVolume(Rmin, BoxBounds)
    OuterShell = SphereVolume(Rmax, BoxBounds)
    Volume = OuterShell - InnerShell
    return Volume

def RDF_AnalyticNorm(Particles, r, dr):
    Global_Gr = zeros(len(r))
    NonEmptyShells = zeros(len(r))
    MaxDist = r[-1] + dr/2
    XList = [p[0] for p in Particles]
    YList = [p[1] for p in Particles]
    ZList = [p[2] for p in Particles]
    BoxBounds = [min(XList), max(XList), min(YList), max(YList), min(ZList), max(ZList)]
    Lx = BoxBounds[1] - BoxBounds[0]
    Ly = BoxBounds[3] - BoxBounds[2]
    Lz = BoxBounds[5] - BoxBounds[4]
    MeanDensity = len(Particles) / (Lx * Ly * Lz)

    for CentralP in range(len(Particles)):
        Local_Gr = zeros(len(r))
        for Neighbour in range(len(Particles)):
            if CentralP != Neighbour:
                dx = Particles[CentralP][0] - Particles[Neighbour][0]
                dy = Particles[CentralP][1] - Particles[Neighbour][1]
                dz = Particles[CentralP][2] - Particles[Neighbour][2]
                d = sqrt(dx**2 + dy**2 + dz**2)
                IdxList = [k for k in range(len(r)) if abs(r[k] - d) <= dr/2]
                for Pos in IdxList:
                    Local_Gr[Pos] += 1
        LocalBox = [BoxBounds[0] - Particles[CentralP][0],
                    BoxBounds[1] - Particles[CentralP][0],
                    BoxBounds[2] - Particles[CentralP][1],
                    BoxBounds[3] - Particles[CentralP][1],
                    BoxBounds[4] - Particles[CentralP][2],
                    BoxBounds[5] - Particles[CentralP][2]]
        for RIdx in range(len(r)):
            SVolume = ShellVolume(r[RIdx] - dr/2, r[RIdx] + dr/2, LocalBox)
            if SVolume > 0.0:
                Local_Gr[RIdx] /= SVolume
                NonEmptyShells[RIdx] += 1
        Local_Gr /= MeanDensity
        Global_Gr += Local_Gr
    for k in range(len(Global_Gr)):
        if NonEmptyShells[k] != 0:
            Global_Gr[k] /= NonEmptyShells[k]
    return Global_Gr
