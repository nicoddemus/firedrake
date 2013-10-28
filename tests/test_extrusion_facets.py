import pytest

from firedrake import *

power = 0
layers = 11


@pytest.fixture
def f():
    mesh = UnitSquareMesh(2 ** power, 2 ** power)
    m = ExtrudedMesh(mesh, layers, layer_height=1.0 / (layers - 1))
    fs = FunctionSpace(m, "CG", 1)
    f = Function(fs)
    f.interpolate(Expression("x[0]"))
    return f


def test_extrusion_external_integral_horizontal(f):
    assert abs(assemble(f * ds_h) - 2.0) < 1.0e-14

def NONtest_extrusion_external_integral_vertical(f):
    assert abs(assemble(f * ds_v) - 2.0) < 1.0e-14

def NONtest_extrusion_external_integral(f):
    assert abs(assemble(f * ds) - 2.0) < 1.0e-14

def NONtest_extrusion_internal_integral_horizontal(f):
    assert abs(assemble(f('+') * dS_h) - 1.0 / (2.0 ** 0.5)) < 1.0e-14

def NONtest_extrusion_internal_integral_vertical(f):
    assert abs(assemble(f('+') * dS_v) - 1.0 / (2.0 ** 0.5)) < 1.0e-14

def NONtest_extrusion_internal_integral(f):
    assert abs(assemble(f('+') * dS) - 1.0 / (2.0 ** 0.5)) < 1.0e-14

def NONtest_extrusion_internal_integral_unit_tri():
    t_mesh = UnitTriangleMesh()
    t = ExtrudedMesh(t_mesh, layers, layer_height=1.0 / (layers - 1))
    V = FunctionSpace(t, 'CG', 1)
    u = Function(V)
    u.interpolate(Expression("x[0]"))
    assert abs(assemble(u('+') * dS) < 1.0e-14)

if __name__ == '__main__':
    import os
    ff = f()
    from IPython import embed
    embed()
    test_extrusion_external_integral_horizontal(ff)
    #pytest.main(os.path.abspath(__file__))
