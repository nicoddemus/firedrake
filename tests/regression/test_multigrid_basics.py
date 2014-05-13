from firedrake import *
import pytest
import numpy as np


@pytest.fixture(scope='module')
def g():
    mesh = UnitSquareMesh(1, 1)
    mh = MeshHierarchy(mesh, 3)
    fs = FunctionSpaceHierarchy(mh, 'DG', 0)

    return FunctionHierarchy(fs)


def test_prolongation_dg0(g):
    for f in g:
        f.assign(0)
    g[0].assign(1)

    for f in g[1:]:
        assert np.allclose(f.dat.data, 0)

    for i in range(len(g) - 1):
        g.prolong(i)

    for f in g:
        assert np.allclose(f.dat.data, 1)


def test_restriction_dg0(g):
    for f in g:
        f.assign(0)
    g[-1].assign(1)

    for f in g[:-1]:
        assert np.allclose(f.dat.data, 0)

    for i in range(len(g) - 1, 0, -1):
        g.restrict(i)

    for f in g:
        assert np.allclose(f.dat.data, 1)


@pytest.mark.xfail(reason="Inverse cell map not correct")
def test_restriction_dg0_x(g):
    for f in g:
        f.assign(0)
    g[-1].interpolate(Expression("x[0]"))

    for i in range(len(g) - 1, 0, -1):
        g.restrict(i)

    for f in g:
        tmp = Function(f)
        tmp.interpolate(Expression("x[0]"))
        assert np.allclose(tmp.dat.data, f.dat.data)


if __name__ == '__main__':
    import os
    pytest.main(os.path.abspath(__file__))
