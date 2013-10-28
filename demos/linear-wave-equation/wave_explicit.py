# NOTE This is a demo, not a regression test
import sys
from firedrake import *

output = True

mesh = Mesh("wave_tank.msh")
# Plumb the space filling curve into UnitSquareMesh after the call to
# gmsh. Doru knows how to do this.

T = 10

dt = 0.001
t = 0
fs = FunctionSpace(mesh, 'Lagrange', 2)
p = Function(fs)
phi = Function(fs)

u = TrialFunction(fs)
v = TestFunction(fs)

outfs = FunctionSpace(mesh, 'Lagrange', 1)

forcing = Function(fs)

bc = [DirichletBC(fs, forcing, 1)]

#p.interpolate(Expression("exp(-40*((x[0]-.5)*(x[0]-.5)+(x[1]-.5)*(x[1]-.5)))"))

if output:
    outfile = File("out.pvd")
    phifile = File("phi.pvd")

    outfile << project(p, outfs)
    #phifile << project(phi, outfs)

# Mass matrix
m = u * v * dx

lump_mass = False

step = 0
while t <= T:
    step += 1

    phi -= dt / 2 * p

    forcing.vector()[:] = sin(2*pi*5*t)

    solve(u * v * dx == v * p * dx + dt * inner(nabla_grad(v), nabla_grad(phi)) * dx,
          p, bcs=bc)

    phi -= dt / 2 * p

    t += dt

    sys.stdout.write("\r"+str(t)+"   ")
    sys.stdout.flush()
    if output and step % 10 == 0:
        outfile << project(p, outfs)
        #phifile << project(phi, outfs)
