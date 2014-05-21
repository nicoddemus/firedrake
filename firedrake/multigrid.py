import numpy as np

from pyop2 import op2
import pyop2.coffee.ast_base as ast

import dmplex
import function
import functionspace
import mesh


__all__ = ['MeshHierarchy', 'FunctionSpaceHierarchy', 'FunctionHierarchy']


class MeshHierarchy(mesh.Mesh):
    """Build a hierarchy of meshes by uniformly refining a coarse mesh"""
    def __init__(self, m, refinement_levels, reorder=True):
        """
        :arg m: the coarse mesh to refine
        :arg refinement_levels: the number of levels of refinement
        :arg reorder: whether to reorder the refined meshes
        """
        m._plex.setRefinementUniform(True)
        dm_hierarchy = m._plex.refineHierarchy(refinement_levels)
        for dm in dm_hierarchy:
            dm.removeLabel("boundary_faces")
            dm.markBoundaryFaces("boundary_faces")
            dm.removeLabel("exterior_facets")
            dm.removeLabel("interior_facets")
            dm.removeLabel("boundary_ids")
            dm.removeLabel("op2_core")
            dm.removeLabel("op2_non_core")
            dm.removeLabel("op2_exec_halo")

        self._hierarchy = [m] + [mesh.Mesh(None, name="%s_refined_%d" % (m.name, i + 1),
                                           plex=dm, distribute=False, reorder=reorder)
                                 for i, dm in enumerate(dm_hierarchy)]

        self._ufl_cell = m.ufl_cell()
        # Simplex only
        factor = 2 ** self.ufl_cell().topological_dimension()
        self._c2f_cells = []

        for mc, mf in zip(self._hierarchy[:-1], self._hierarchy[1:]):
            if not hasattr(mc, '_cell_new_to_old'):
                o, n = dmplex.get_entity_renumbering(mc._plex, mc._cell_numbering, "cell")
                mc._cell_old_to_new = o
                mc._cell_new_to_old = n
            if not hasattr(mf, '_cell_old_to_new'):
                o, n = dmplex.get_entity_renumbering(mf._plex, mf._cell_numbering, "cell")
                mf._cell_old_to_new = o
                mf._cell_new_to_old = n

            cback = mc._cell_new_to_old
            fforward = mf._cell_old_to_new
            ofcells = np.dstack([(cback * factor) + i for i in range(factor)]).flatten()
            fcells = fforward[ofcells]
            self._c2f_cells.append(fcells.reshape(-1, factor))

    def __iter__(self):
        for m in self._hierarchy:
            yield m

    def __len__(self):
        return len(self._hierarchy)

    def __getitem__(self, idx):
        return self._hierarchy[idx]


class FunctionSpaceHierarchy(object):
    """Build a hierarchy of function spaces.

    Given a hierarchy of meshes, this constructs a hierarchy of
    function spaces, with the property that every coarse space is a
    subspace of the fine spaces that are a refinement of it.
    """
    def __init__(self, mesh_hierarchy, family, degree):
        """
        :arg mesh_hierarchy: a :class:`.MeshHierarchy` to build the
             function spaces on.
        :arg family: the function space family
        :arg degree: the degree of the function space
        """
        self._mesh_hierarchy = mesh_hierarchy
        self._hierarchy = [functionspace.FunctionSpace(m, family, degree)
                           for m in self._mesh_hierarchy]

        self._map_cache = {}
        self._ufl_element = self[0].ufl_element()

    def __len__(self):
        return len(self._hierarchy)

    def __iter__(self):
        for fs in self._hierarchy:
            yield fs

    def __getitem__(self, idx):
        return self._hierarchy[idx]

    def ufl_element(self):
        return self._ufl_element

    def cell_node_map(self, level, bcs=None):
        """A :class:`pyop2.Map` from cells on a coarse mesh to the
        corresponding degrees of freedom on a the fine mesh below it.

        :arg level: the coarse level the map should be from.
        :arg bcs: optional iterable of :class:`.DirichletBC`\s
             (currently ignored).
        """
        if not 0 <= level < len(self) - 1:
            raise RuntimeError("Requested coarse level %d outside permissible range [0, %d)" %
                               (level, len(self) - 1))
        try:
            return self._map_cache[level]
        except KeyError:
            pass
        Vc = self._hierarchy[level]
        Vf = self._hierarchy[level + 1]

        family = self.ufl_element().family()
        degree = self.ufl_element().degree()

        c2f = self._mesh_hierarchy._c2f_cells[level]
        if family == "Discontinuous Lagrange":
            if degree != 0:
                raise RuntimeError
            arity = Vf.cell_node_map().arity * c2f.shape[1]
            map_vals = Vf.cell_node_map().values_with_halo[c2f].flatten()

            map = op2.Map(Vc.mesh().cell_set,
                          Vf.node_set,
                          arity,
                          map_vals)

            self._map_cache[level] = map
            return map

        if family == "Lagrange":
            if degree != 1:
                raise RuntimeError
            map_vals = dmplex.p1_coarse_fine_map(Vc, Vf, c2f)

            arity = map_vals.shape[1]
            map = op2.Map(Vc.mesh().cell_set, Vf.node_set, arity, map_vals)

            self._map_cache[level] = map
            return map


class FunctionHierarchy(object):
    """Build a hierarchy of :class:`~.Function`\s"""
    def __init__(self, fs_hierarchy):
        """
        :arg fs_hierarchy: the :class:`~.FunctionSpaceHierarchy` to build on.

        `fs_hierarchy` may also be an existing
        :class:`FunctionHierarchy`, in which case a copy of the
        hierarchy is returned.
        """
        if isinstance(fs_hierarchy, FunctionHierarchy):
            self._function_space = fs_hierarchy.function_space()
        else:
            self._function_space = fs_hierarchy

        self._hierarchy = [function.Function(f) for f in fs_hierarchy]

    def __iter__(self):
        for f in self._hierarchy:
            yield f

    def __len__(self):
        return len(self._hierarchy)

    def __getitem__(self, idx):
        return self._hierarchy[idx]

    def function_space(self):
        return self._function_space

    def cell_node_map(self, i):
        return self._function_space.cell_node_map(i)

    def prolong(self, level):
        """Prolong from a coarse to the next finest hierarchy level.

        :arg level: The coarse level to prolong from"""

        if not 0 <= level < len(self) - 1:
            raise RuntimeError("Requested coarse level %d outside permissible range [0, %d)" %
                               (level, len(self) - 1))
        fs = self[level].function_space()
        family = fs.ufl_element().family()
        degree = fs.ufl_element().degree()

        if family == "Discontinuous Lagrange":
            if degree != 0:
                raise RuntimeError("Can only prolong P0 fields, not P%dDG" % degree)
            self._prolong_dg0(level)
        elif family == "Lagrange":
            if degree != 1:
                raise RuntimeError("Can only prolong P1 fields, not P%d" % degree)
            self._prolong_cg1(level)
        else:
            raise RuntimeError("Prolongation only implemented for P0DG and P1")

    def restrict(self, level):
        """Restrict from a fine to the next coarsest hierarchy level.

        :arg level: The fine level to restrict from"""

        if not 0 < level < len(self):
            raise RuntimeError("Requested fine level %d outside permissible range [1, %d)" %
                               (level, len(self)))

        fs = self[level].function_space()
        family = fs.ufl_element().family()
        degree = fs.ufl_element().degree()
        if family == "Discontinuous Lagrange":
            if degree == 0:
                self._restrict_dg0(level)
            else:
                raise RuntimeError("Can only restrict P0 fields, not P%dDG" % degree)
        elif family == "Lagrange":
            if degree != 1:
                raise RuntimeError("Can only restrict P1 fields, not P%d" % degree)
            self._restrict_cg1(level)
        else:
            raise RuntimeError("Restriction only implemented for P0DG and P1")

    def _prolong_dg0(self, level):
        c2f_map = self.cell_node_map(level)
        coarse = self[level]
        fine = self[level + 1]
        if not hasattr(self, '_prolong_kernel'):
            k = ast.FunDecl("void", "prolong_dg0",
                            [ast.Decl(coarse.dat.ctype, "**coarse"),
                             ast.Decl(fine.dat.ctype, "**fine")],
                            body=ast.c_for("fdof", c2f_map.arity,
                                           ast.Assign(ast.Symbol("fine", ("fdof", 0)),
                                                      ast.Symbol("coarse", (0, 0))),
                                           pragma=None),
                            pred=["static", "inline"])
            self._prolong_kernel = op2.Kernel(k, "prolong_dg0")
        op2.par_loop(self._prolong_kernel, coarse.cell_set,
                     coarse.dat(op2.READ, coarse.cell_node_map()),
                     fine.dat(op2.WRITE, c2f_map))

    def _restrict_dg0(self, level):
        c2f_map = self.cell_node_map(level - 1)
        coarse = self[level - 1]
        fine = self[level]
        if not hasattr(self, '_restrict_kernel'):
            k = ast.FunDecl("void", "restrict_dg0",
                            [ast.Decl(coarse.dat.ctype, "**coarse"),
                             ast.Decl(fine.dat.ctype, "**fine")],
                            body=ast.Block([ast.Decl(coarse.dat.ctype, "tmp", init=0.0),
                                            ast.c_for("fdof", c2f_map.arity,
                                                      ast.Incr(ast.Symbol("tmp"),
                                                               ast.Symbol("fine", ("fdof", 0))),
                                                      pragma=None),
                                            ast.Assign(ast.Symbol("coarse", (0, 0)),
                                                       ast.Div(ast.Symbol("tmp"),
                                                               c2f_map.arity))]),
                            pred=["static", "inline"])
            self._restrict_kernel = op2.Kernel(k, "restrict_dg0")

        op2.par_loop(self._restrict_kernel, coarse.cell_set,
                     coarse.dat(op2.WRITE, coarse.cell_node_map()),
                     fine.dat(op2.READ, c2f_map))

    def _prolong_cg1(self, level):
        c2f_map = self.cell_node_map(level)
        coarse = self[level]
        fine = self[level + 1]
        if not hasattr(self, '_prolong_kernel'):
            # Only 2D for now.
            # Due to smart map, fine field is:
            # u_f = A u_c
            # Where A is [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]
            k = """void prolong_cg1(double **coarse, double **fine)
            {
                static const double coeffs[6][3] = {{1, 0, 0},
                                                    {0, 1, 0},
                                                    {0, 0, 1},
                                                    {0, 0.5, 0.5},
                                                    {0.5, 0, 0.5},
                                                    {0.5, 0.5, 0}};
                for ( int i = 0; i < 6; i++ ) {
                    fine[i][0] = 0;
                    for ( int j = 0; j < 3; j++ ) {
                        fine[i][0] += coeffs[i][j] * coarse[j][0];
                    }
                }
            }"""
            self._prolong_kernel = op2.Kernel(k, "prolong_cg1")
        op2.par_loop(self._prolong_kernel, coarse.cell_set,
                     coarse.dat(op2.READ, coarse.cell_node_map()),
                     fine.dat(op2.WRITE, c2f_map))

    def _restrict_cg1(self, level):
        c2f_map = self.cell_node_map(level - 1)
        coarse = self[level - 1]
        fine = self[level]
        if not hasattr(self, '_restrict_kernel'):
            element = coarse.function_space().fiat_element
            import FIAT
            quadrature = FIAT.make_quadrature(element.ref_el, 2)
            weights = quadrature.get_weights()
            points = quadrature.get_points()

            fine_basis = element.tabulate(0, points).values()[0]
            coarse_basis = None

            k = """
            #include "firedrake_geometry.h"
            static inline void restrict_cg1(double **coarse, double **fine, double **coordinates)
            {
            const double fine_basis[4][3] = %(fine_basis)s;
            const double weight[4] = %(weight)s;
            double coarse_[3] = {0};

            double J[4];
            compute_jacobian_triangle_2d(J, coordinates);
            double K[4];
            double detJ;
            compute_jacobian_inverse_triangle_2d(K, detJ, J);
            const double det = fabs(detJ);
            for ( int fcell = 0; fcell < 4; fcell++ ) {
                for ( int ip = 0; ip < 4; ip++ ) {
                    double fine_coeff = 0;
                    for ( int i = 0; i < 3; i++ ) {
                        fine_coeff += fine[fcell*3 + i][0] * fine_basis[ip][i];
                    }
                    for ( int i = 0; i < 3; i++ ) {
                        coarse_[i] += (weight[ip] * fine_coeff * fine_basis[ip][i]) * det/4.0;
                    }
                }
            }
            for ( int i = 0; i < 3; i++ ) {
                coarse[i][0] = coarse_[i];
            }
            }
            """ % {"fine_basis": "{{" + "},\n{".join([", ".join(map(str, x)) for x in fine_basis.T])+"}}",
                   "weight": "{"+", ".join(["%s" % w for w in weights]) +"}"}

            k = op2.Kernel(k, 'restrict_cg1', include_dirs=["/data/lmitche1/src/firedrake/firedrake"])
            c2f = self._function_space._mesh_hierarchy._c2f_cells[level - 1]

            arity = fine.cell_node_map().arity * c2f.shape[1]
            map_vals = fine.cell_node_map().values_with_halo[c2f].flatten()

            map_ = op2.Map(coarse.cell_set,
                          fine.node_set,
                          arity,
                          map_vals)


            ccoords = coarse.function_space().mesh().coordinates
            op2.par_loop(k, coarse.cell_set,
                         coarse.dat(op2.WRITE, coarse.cell_node_map()),
                         fine.dat(op2.READ, map_),
                         ccoords.dat(op2.READ, ccoords.cell_node_map(), flatten=True))
