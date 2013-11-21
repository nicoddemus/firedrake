import pytest

from firedrake import *


@pytest.fixture(scope='session')
def ext_mesh():
    return ExtrudedMesh(UnitSquareMesh(1, 1), 11, layer_height=0.1)


@pytest.fixture(scope='session')
def mesh(request):
    return UnitSquareMesh(3, 3)


def extmesh(nx, ny, nz):
    return ExtrudedMesh(UnitSquareMesh(nx, ny), nz + 1, layer_height=1.0 / nz)


@pytest.fixture(scope='session')
def cg1(mesh):
    return FunctionSpace(mesh, "CG", 1)


@pytest.fixture(scope='session')
def cg2(mesh):
    return FunctionSpace(mesh, "CG", 2)


@pytest.fixture(scope='session')
def dg0(mesh):
    return FunctionSpace(mesh, "DG", 0)


@pytest.fixture(scope='session')
def dg1(mesh):
    return FunctionSpace(mesh, "DG", 1)


@pytest.fixture(scope='session')
def vcg1(mesh):
    return VectorFunctionSpace(mesh, "CG", 1)


@pytest.fixture(scope='session')
def cg1cg1(cg1):
    return cg1 * cg1


@pytest.fixture(scope='session')
def cg1dg0(cg1, dg0):
    return cg1 * dg0


@pytest.fixture(scope='session')
def cg2dg1(cg2, dg1):
    return cg2 * dg1


@pytest.fixture(scope='session')
def cg1vcg1(mesh, cg1, vcg1):
    return cg1 * vcg1


@pytest.fixture(scope='session')
def ecg1(ext_mesh):
    return FunctionSpace(ext_mesh, "CG", 1)


@pytest.fixture(scope='session')
def ecg2(ext_mesh):
    return FunctionSpace(ext_mesh, "CG", 2)


@pytest.fixture(scope='session')
def edg0(ext_mesh):
    return FunctionSpace(ext_mesh, "DG", 0)


@pytest.fixture(scope='session')
def edg1(ext_mesh):
    return FunctionSpace(ext_mesh, "DG", 1)


@pytest.fixture(scope='session')
def evcg1(ext_mesh):
    return VectorFunctionSpace(ext_mesh, "CG", 1)


@pytest.fixture(scope='session')
def ecg1cg1(ecg1):
    return ecg1 * ecg1


@pytest.fixture(scope='session')
def ecg1dg0(ecg1, edg0):
    return ecg1 * edg0


@pytest.fixture(scope='session')
def ecg2dg1(ecg2, edg1):
    return ecg2 * edg1


@pytest.fixture(scope='session')
def ecg1vcg1(ext_mesh, ecg1, evcg1):
    return ecg1 * evcg1
