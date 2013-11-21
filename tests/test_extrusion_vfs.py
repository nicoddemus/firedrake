import pytest
from firedrake import *


def vfs_1():
    m = UnitSquareMesh(1, 1)
    mesh = ExtrudedMesh(m, 11, layer_height=0.1)

    expr1 = Expression("x[2]")
    exactfs1 = FunctionSpace(mesh, "CG", 1)
    exact1 = Function(exactfs1)
    exact1.interpolate(expr1)
    return abs(1.0 / 3.0 - assemble(exact1 * exact1 * dx))  # 1/3


def vfs_2():
    m = UnitSquareMesh(1, 1)
    mesh = ExtrudedMesh(m, 11, layer_height=0.1)
    expr2 = Expression(("0.0", "0.0", "x[2]"))
    exactfs2 = VectorFunctionSpace(mesh, "CG", 1)
    exact2 = Function(exactfs2)
    exact2.interpolate(expr2)
    return abs(1.0 / 3.0 - assemble(dot(exact2, exact2) * dx))  # 1/3


def vfs_3():
    m = UnitSquareMesh(1, 1)
    mesh = ExtrudedMesh(m, 11, layer_height=0.1)

    V = VectorFunctionSpace(mesh, "CG", 1)
    v = Function(V)
    v.assign(1.0)
    p = TestFunction(V)
    F = inner(p, v) * dx

    out = assemble(F)
    return abs(3.0 - sum(sum(out.dat.data)))


def vfs_4():
    m = UnitSquareMesh(1, 1)
    mesh = ExtrudedMesh(m, 11, layer_height=0.1)
    V = VectorFunctionSpace(mesh, "CG", 1)
    expr = Expression(("x[0]", "x[1]", "x[2]"))
    ff = project(expr, V)  # BOOM
    return abs(ff.dat.data - mesh._coordinates)


def vfs_5():
    m = UnitSquareMesh(1, 1)
    mesh = ExtrudedMesh(m, 11, layer_height=0.1)

    V1 = FunctionSpace(mesh, "CG", 1)
    f = Function(V1)
    f.interpolate(Expression("x[0]"))
    W = V1 * V1
    u, p = TrialFunctions(W)
    v, q = TestFunctions(W)
    a = u * v * dx + p * q * dx
    L = f * v * dx + f * q * dx

    out = Function(W)
    solve(a == L, out, solver_parameters={
          'pc_type': 'jacobi', 'ksp_type': 'cg'})
    return np.array([max(abs(out.dat.data[0] - f.dat.data)), max(abs(out.dat.data[1] - f.dat.data))])


def vfs_6():
    m = UnitSquareMesh(1, 1)
    mesh = ExtrudedMesh(m, 11, layer_height=0.1)

    V1 = FunctionSpace(mesh, "CG", 1)
    V2 = FunctionSpace(mesh, "CG", 2)

    f = Function(V1)
    g = Function(V2)
    f.interpolate(Expression("x[0]"))
    g.interpolate(Expression("x[0]"))

    W = V1 * V2
    u, p = TrialFunctions(W)
    v, q = TestFunctions(W)

    a = u * v * dx + p * q * dx
    L = f * v * dx + g * q * dx

    out = Function(W)
    solve(a == L, out, solver_parameters={
          'pc_type': 'jacobi', 'ksp_type': 'cg'})
    return np.array([max(abs(out.dat.data[0] - f.dat.data)), max(abs(out.dat.data[1] - g.dat.data))])


def not_vfs_7():
    m = UnitSquareMesh(1, 1)
    mesh = ExtrudedMesh(m, 11, layer_height=0.1)

    V1 = FunctionSpace(mesh, "CG", 1)
    V2 = VectorFunctionSpace(mesh, "CG", 1)

    f = Function(V1)
    g = Function(V2)
    f.interpolate(Expression("x[0]"))
    g.interpolate(Expression(("-x[1]", "x[0]")))

    W = V1 * V2
    u, p = TrialFunctions(W)
    v, q = TestFunctions(W)
    a = u * v * dx + dot(p, q) * dx
    L = f * v * dx + dot(g, q) * dx

    out = Function(W)
    solve(a == L, out, solver_parameters={
          'pc_type': 'jacobi', 'ksp_type': 'cg'})
    print out.dat.data


def test_vector_function_spaces():
    assert vfs_1() < 1.0e-14
    assert vfs_2() < 1.0e-14
    assert vfs_3() < 1.0e-13
    assert (vfs_4() < 1.0e-14).all()
    assert (vfs_5() < 1.0e-14).all()
    assert (vfs_6() < 1.0e-14).all()


if __name__ == '__main__':
    import os
    pytest.main(os.path.abspath(__file__))
