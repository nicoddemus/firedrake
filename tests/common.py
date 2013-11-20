import pytest

from firedrake import *


@pytest.fixture(scope='module', params=['extruded','non-extruded'])
def mesh(request):
    return {'extruded': ExtrudedMesh(UnitSquareMesh(1,1), 5, layer_height=0.25),
            'non-extruded': UnitSquareMesh(3, 3)}[request.param]


def extmesh(nx, ny, nz):
    return ExtrudedMesh(UnitSquareMesh(nx, ny), nz+1, layer_height=1.0/nz)


@pytest.fixture(scope='module')
def cg1(mesh):
    return FunctionSpace(mesh, "CG", 1)


@pytest.fixture(scope='module')
def cg2(mesh):
    return FunctionSpace(mesh, "CG", 2)


@pytest.fixture(scope='module')
def dg0(mesh):
    return FunctionSpace(mesh, "DG", 0)


@pytest.fixture(scope='module')
def dg1(mesh):
    return FunctionSpace(mesh, "DG", 1)


@pytest.fixture(scope='module')
def vcg1(mesh):
    return VectorFunctionSpace(mesh, "CG", 1)


@pytest.fixture(scope='module')
def cg1cg1(cg1):
    return cg1 * cg1


@pytest.fixture(scope='module')
def cg1dg0(cg1, dg0):
    return cg1 * dg0


@pytest.fixture(scope='module')
def cg2dg1(cg2, dg1):
    return cg2 * dg1


@pytest.fixture(scope='module')
def cg1vcg1(mesh, cg1, vcg1):
    return cg1 * vcg1
