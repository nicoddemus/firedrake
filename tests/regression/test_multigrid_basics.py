from firedrake import *
import pytest
import numpy as np


@pytest.fixture(scope='module')
def dg0():
    mesh = UnitSquareMesh(1, 1)
    mh = MeshHierarchy(mesh, 3)
    fs = FunctionSpaceHierarchy(mh, 'DG', 0)

    return FunctionHierarchy(fs)


@pytest.fixture(scope='module')
def dg0_noreorder():
    mesh = UnitSquareMesh(1, 1, reorder=False)
    mh = MeshHierarchy(mesh, 3, reorder=False)
    fs = FunctionSpaceHierarchy(mh, 'DG', 0)

    return FunctionHierarchy(fs)


@pytest.fixture(scope='module')
def cg1_noreorder():
    mesh = UnitSquareMesh(1, 1, reorder=False)
    mh = MeshHierarchy(mesh, 3, reorder=False)
    fs = FunctionSpaceHierarchy(mh, 'CG', 1)

    return FunctionHierarchy(fs)


@pytest.fixture(scope='module')
def cg1():
    mesh = UnitSquareMesh(1, 1)
    mh = MeshHierarchy(mesh, 3)
    fs = FunctionSpaceHierarchy(mh, 'CG', 1)

    return FunctionHierarchy(fs)


def test_prolongation_dg0(dg0):
    for f in dg0:
        f.assign(0)
    dg0[0].assign(1)

    for f in dg0[1:]:
        assert np.allclose(f.dat.data, 0)

    for i in range(len(dg0) - 1):
        dg0.prolong(i)

    for f in dg0:
        assert np.allclose(f.dat.data, 1)


def test_restriction_dg0(dg0):
    for f in dg0:
        f.assign(0)
    dg0[-1].assign(1)

    for f in dg0[:-1]:
        assert np.allclose(f.dat.data, 0)

    for i in range(len(dg0) - 1, 0, -1):
        dg0.restrict(i)

    for f in dg0:
        assert np.allclose(f.dat.data, 1)


@pytest.mark.xfail(reason="Inverse cell map not correct")
def test_restriction_dg0_coords(dg0):
    for f in dg0:
        f.assign(0)
    dg0[-1].interpolate(Expression("x[0]"))

    for i in range(len(dg0) - 1, 0, -1):
        dg0.restrict(i)

    for f in dg0:
        tmp = Function(f)
        tmp.interpolate(Expression("x[0]"))
        assert np.allclose(tmp.dat.data, f.dat.data)


def test_restriction_dg0_coords_no_reorder(dg0_noreorder):
    for f in dg0_noreorder:
        f.assign(0)
    dg0_noreorder[-1].interpolate(Expression("x[0]"))

    for i in range(len(dg0_noreorder) - 1, 0, -1):
        dg0_noreorder.restrict(i)

    for f in dg0_noreorder:
        tmp = Function(f)
        tmp.interpolate(Expression("x[0]"))
        assert np.allclose(tmp.dat.data, f.dat.data)


@pytest.mark.xfail(reason="Inverse cell map not correct")
def test_prolongation_cg1_coords(cg1):
    for f in cg1:
        f.assign(0)

    cg1[0].interpolate(Expression("x[0]"))

    for i in range(len(cg1) - 1):
        cg1.prolong(i)

    for f in cg1:
        tmp = Function(f)
        tmp.interpolate(Expression("x[0]"))
        assert np.allclose(tmp.dat.data, f.dat.data)


def test_prolongation_cg1_coords_no_reorder(cg1_noreorder):
    for f in cg1_noreorder:
        f.assign(0)

    cg1_noreorder[0].interpolate(Expression("x[0]"))

    for i in range(len(cg1_noreorder) - 1):
        cg1_noreorder.prolong(i)

    for f in cg1_noreorder:
        tmp = Function(f)
        tmp.interpolate(Expression("x[0]"))
        assert np.allclose(tmp.dat.data, f.dat.data)


if __name__ == '__main__':
    import os
    pytest.main(os.path.abspath(__file__))
