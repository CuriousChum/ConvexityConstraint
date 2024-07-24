import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from math import comb

from collections import defaultdict
from scipy.spatial import ConvexHull
import re
import ast


def _generate_sum_exactly_k(n_vars: int, k: int):
    if n_vars < 1:
        raise Exception("cannot sum less than 0 elements")
    elif n_vars == 1:
        yield (k,)
    else:
        for i in range(k + 1):
            for alpha in _generate_sum_exactly_k(n_vars - 1, k - i):
                yield alpha + (i,)


def _generate_sum_upto_k(n_vars: int, k: int):
    if n_vars < 1:
        raise Exception("cannot sum less than 0 elements")
    for cur_sum in range(k + 1):
        for alpha in _generate_sum_exactly_k(n_vars, cur_sum):
            yield alpha


def create_adj_list_and_neighbor_dct(faces):
    adj_list = defaultdict(list)
    nbr_dct = defaultdict(set)
    for u, v, w in faces:
        if (v, w) not in adj_list[u] and (w, v) not in adj_list[u]:
            adj_list[u].append((v, w))
        nbr_dct[u].add(v)
        nbr_dct[u].add(w)
        if (w, u) not in adj_list[v] and (u, w) not in adj_list[v]:
            adj_list[v].append((w, u))
        nbr_dct[v].add(w)
        nbr_dct[v].add(u)
        if (u, v) not in adj_list[w] and (v, u) not in adj_list[w]:
            adj_list[w].append((u, v))
        nbr_dct[w].add(u)
        nbr_dct[w].add(v)
    return adj_list, nbr_dct


def min_points_for_poly(dim: int, degree: int):
    return comb(dim + degree, degree)


def polynomial_vct(degree: int, coord: np.array,
                   grad=False):
    # NOTE: Check if same as np.polynomial vandermonde func
    dim = coord.shape[0]
    if grad:
        polynomial_terms = np.ones(
            (min_points_for_poly(dim, degree), dim), dtype=np.double)
        for j, alpha in enumerate(_generate_sum_upto_k(dim, degree)):
            for i, alpha_i in enumerate(alpha):
                partial_powers = np.array(alpha)
                partial_powers[i] = max(partial_powers[i] - 1, 0)

                partial_factors = np.ones_like(partial_powers)
                partial_factors[i] = alpha_i

                polynomial_terms[j, i] *= np.prod(
                    partial_factors
                    * (coord ** partial_powers),
                )
    else:
        polynomial_terms = np.fromiter(
            (np.prod(coord ** alpha)
             for alpha in _generate_sum_upto_k(dim, degree)),
            dtype=np.double
        )

    return polynomial_terms


def polynomial_eval(degree: int, coeff: np.array, coord: np.array,
                    grad=False):
    return polynomial_vct(degree, coord, grad).T @ coeff


def __experimental_find_containing_entity_from_coords(coords, nodes,
                                                      adj_list):
    dists = np.linalg.norm(coords[:, np.newaxis] - nodes, axis=2)
    print(dists)
    closest_indices = np.argmin(dists, axis=1)
    print(closest_indices)

    closest_nodes = nodes[closest_indices]
    raise NotImplementedError
    return closest_nodes


def find_containing_entity(coord: np.array, nodes: np.array,
                           adj_list: dict):
    dist = np.linalg.norm(nodes - coord, axis=1)
    closest_idx = np.argmin(dist)
    if dist[closest_idx] == 0:
        # the entity containing coord is a vertex/node
        return (closest_idx,)

    closest_node = nodes[closest_idx]
    for i, (nbr1, nbr2) in enumerate(adj_list[closest_idx]):
        nbr1_vct = nodes[nbr1] - closest_node
        nbr2_vct = nodes[nbr2] - closest_node
        coefs = np.linalg.inv(
            np.column_stack([nbr1_vct, nbr2_vct])
        ) @ (coord - closest_node)

        if np.all(coefs >= 0):
            return closest_idx, nbr1, nbr2
    raise Exception("Containing entity not found")


class LagrangeBasis:
    def __init__(self, nodes, adj_list, degree: int = 1):
        self.nodes = nodes
        self.degree = degree
        self.adj_list = adj_list

    def eval_at(self, coords):
        res = np.zeros((coords.shape[0], self.nodes.shape[0]),
                       dtype=np.double)
        for i, coord in enumerate(coords):
            vertices = find_containing_entity(coord, self.nodes,
                                              self.adj_list)
            cur_basis_vals = res[i]
            if len(vertices) == 1:
                cur_basis_vals[vertices[0]] = 1
            elif len(vertices) == 2:
                coord_vct = coord - self.nodes[vertices[0]]
                edge_vct = self.nodes[vertices[1]]
                - self.nodes[vertices[0]]
                scale_factor = coord_vct / edge_vct
                assert (np.allclose(
                    scale_factor,
                    np.ones_like(scale_factor)
                        * np.average(scale_factor)
                        ))
                scale_factor = np.average(scale_factor)
                cur_basis_vals[vertices[0]] = 1 - scale_factor
                cur_basis_vals[vertices[1]] = scale_factor
            else:
                for i, vert_i in enumerate(vertices):
                    coord_vct = coord - self.nodes[vert_i]
                    edge_1_vct = self.nodes[
                        vertices[(i + 1) % len(vertices)]]
                    - self.nodes[vert_i]
                    edge_2_vct = self.nodes[
                        vertices[(i + 2) % len(vertices)]]
                    - self.nodes[vert_i]
                    # find a, b such that a*u + b*v = c, here u,v are edges
                    # c is coord_vct
                    coefs = np.linalg.inv(
                        np.column_stack([edge_1_vct, edge_2_vct])) \
                        @ coord_vct

                    cur_basis_vals[vert_i] = 1 - np.sum(coefs)

        return res


class FunctionalEstimateOnTriangulation:
    def __init__(self, points, func_vals, faces,
                 degree=1,
                 grad_est_at_nodes=None,
                 hess_est_at_nodes=None,
                 ):
        self.points = points
        self.faces = faces
        self.func_vals = func_vals
        self.dim = points.shape[1]
        self.degree = degree
        self.adj_list, self.neighbors = create_adj_list_and_neighbor_dct(faces)

        self.grad_est_at_nodes = grad_est_at_nodes
        self.hess_est_at_nodes = hess_est_at_nodes
        self.basis = LagrangeBasis(self.points, self.adj_list, degree=degree)

    def in_hull(self, coord):
        raise NotImplementedError

    def get_neighbors_by_levels(self, idx, level=1):
        cur_level_idx, tmp = {idx}, set()
        res = set()
        while level > 0:
            while cur_level_idx:
                cur_idx = cur_level_idx.pop()
                for nbr in self.neighbors[cur_idx]:
                    res.add(nbr)
                    tmp.add(nbr)
            level -= 1
            cur_level_idx = tmp.copy()
        return res

    def calc_pointwise_grad_est_from_func(self, idx, f_vals):
        min_req_pts = min_points_for_poly(self.dim, self.degree + 1)
        level = 1
        nbr_indices = self.get_neighbors_by_levels(idx, level)
        while len(nbr_indices) < 16:
            level += 1
            nbr_indices = self.get_neighbors_by_levels(idx, level)
        A = np.empty((min_req_pts, len(nbr_indices))).T
        b = np.empty(len(nbr_indices))
        for i, nbr_i in enumerate(nbr_indices):
            nbr_node = self.points[nbr_i]
            A[i] = polynomial_vct(self.degree + 1, nbr_node)
            b[i] = f_vals[nbr_i]
        coefs = np.linalg.pinv(A) @ b
        return polynomial_eval(self.degree + 1, coefs, self.points[idx],
                               grad=True)

    def calculate_pointwise_grad_est(self, idx):
        return self.calc_pointwise_grad_est_from_func(idx, self.func_vals)

    def calc_entire_grad_est_from_func(self, f_vals):
        grads = np.empty((len(self.points), self.dim),
                         dtype=np.double
                         )
        for i in range(len(self.points)):
            grads[i] = self.calc_pointwise_grad_est_from_func(
                i, self.func_vals)
        return grads

    def calculate_entire_grad_est(self, return_val=False):
        self.grad_est_at_nodes = self.calc_entire_grad_est_from_func(
            self.func_vals)
        return self.grad_est_at_nodes.copy() if return_val else None

    def calc_ppr(self, grad_est, coords):
        basis_vals = self.basis.eval_at(coords)
        return basis_vals @ grad_est

    def calculate_grad_ppr(self, coords):
        if self.grad_est_at_nodes is None:
            self.calculate_entire_grad_est()
        return self.calc_ppr(self.grad_est_at_nodes, coords)

    def calculate_entire_hess_est(self, return_val=False):
        if self.grad_est_at_nodes is None:
            self.calculate_entire_grad_est()
        grad_x, grad_y = self.grad_est_at_nodes.T
        hess_x = np.empty((len(self.points), self.dim),
                          dtype=np.double
                          )
        hess_y = np.empty((len(self.points), self.dim),
                          dtype=np.double
                          )

        for i in range(len(self.points)):
            hess_x[i] = self.calc_pointwise_grad_est_from_func(i, grad_x)
            hess_y[i] = self.calc_pointwise_grad_est_from_func(i, grad_y)

        self.hess_est_at_nodes = np.stack([hess_x, hess_y], axis=1)
        return self.hess_est_at_nodes.copy() if return_val else None


def get_points_and_ppr_on_square(fe: FunctionalEstimateOnTriangulation,
                                 num: int):
    points = np.linspace(1, 2, num=num)
    points = np.linspace(
        np.column_stack((points, np.ones_like(points))),
        np.column_stack((points, np.zeros_like(points) + 2)),
        num=num
    ).reshape(-1, 2)
    hess_x = fe.hess_est_at_nodes[:, :, 0]
    hess_y = fe.hess_est_at_nodes[:, :, 1]

    ppr_hess_x = fe.calc_ppr(hess_x, points)
    ppr_hess_y = fe.calc_ppr(hess_y, points)

    print(ppr_hess_x.shape, ppr_hess_y.shape)
    return points, np.stack([ppr_hess_x, ppr_hess_y], axis=1)


def plot_points_and_faces(points, faces):
    points = np.asanyarray(points)
    edges = np.array(
        [((a, b), (b, c), (c, a)) for a, b, c in faces]
    ).reshape(-1, 2)
    plt.plot(points[edges.T][0], points[edges.T][1])


def parse_wolfram_format(filename: str):
    with open(filename, "r") as f:
        data = f.read()
        parts = re.findall(r"\"(\w+)\": {([0-9\.{}, \-+e]*)}[},]", data)
        dct = {}
        for name, val in parts:
            dct[name] = val.replace("{", "[").replace("}", "]")
            dct[name] = ast.literal_eval(dct[name])


def parabola(coords):
    return np.linalg.norm(np.asanyarray(coords)) ** 2


if __name__ == "__main__":
    with open("cvx.txt", "r") as f:
        data = f.read()
        parts = re.findall(r"\"(\w+)\" ?-> {([0-9\.{}, \-+e]*)}[},]", data)
        dct = {}
        for name, val in parts:
            dct[name] = val.replace("{", "[").replace("}", "]")
            dct[name] = ast.literal_eval(dct[name])

    for i in range(len(dct["points"])):
        ori_pts = np.array(
            [(x, y, val) for _, val, [x, y] in dct["points"]],
            dtype=np.double)
        dct["points"][i][1] = parabola(
            dct["points"][i][2]) - dct["points"][i][1]

    cvx_hull = ConvexHull(ori_pts)
    to_plot = np.array([(x, y, val) for _, val, [x, y] in dct["points"]])

    fig, ax = plt.subplots(1, 1,
                           # subplot_kw={"projection": "3d"}
                           )

    x, y, z = np.array(sorted((x, y, z) for x, y, z in to_plot)).T

    # ax.plot_trisurf(x, y, z, vmin=z.min())
    edges = np.array([((a, b), (a, c), (b, c))
                      for (a, b, c), _, _ in dct["faces"]]) \
        .reshape(-1, 2)

    faces = np.array(
        [(a, b, c) for (a, b, c), _, _ in dct["faces"]],
        dtype=np.int64)

    fig2, ax2 = plt.subplots(1, 1, subplot_kw={"projection": "3d"})

    points = np.column_stack([x, y])
    fn_est = FunctionalEstimateOnTriangulation(
        points,
        z,
        faces,
        # get_delaunay_triangulation_faces(points),
    )

    grad_est = fn_est.calculate_entire_grad_est(return_val=True)
    hess_est = fn_est.calculate_entire_hess_est(return_val=True)
    dethess_est = [np.linalg.det(h) for h in hess_est]
    hessrank_est = [np.linalg.matrix_rank(h, tol=2e-1) for h in hess_est]

    # ax2.plot_trisurf(x, y, z, vmin=z.min(), color="green")
    # ax2.plot_trisurf(x, y, dethess_est, color="green")
    # ax2.scatter(x, y, hessrank_est)

    points, hess_ppr_est = get_points_and_ppr_on_square(fn_est, 100)
    dethess_ppr_est = np.array([np.linalg.det(h) for h in hess_ppr_est])
    hessrank_ppr_est = np.array([np.linalg.matrix_rank(h, tol=5e-1)
                                 for h in hess_ppr_est])

    fig3, ax3 = plt.subplots(1, 1, subplot_kw={"projection": "3d"})
    fig4, ax4 = plt.subplots(1, 1, subplot_kw={"projection": "3d"})

    ax2.scatter(points[:, 0], points[:, 1], hessrank_ppr_est)
    ax3.scatter(points[:, 0], points[:, 1], dethess_ppr_est)

    # ax2.scatter(x, y, dethess_est, color="green")
    ax2.scatter(x, y, hessrank_est, color="green")
    ax4.scatter(x, y, dethess_est, color="green")

    # u_from = np.column_stack([x[edges.T[0]], y[edges.T[0]],
    #                           np.zeros_like(edges.T[0])])
    # u_to = np.column_stack([x[edges.T[1]], y[edges.T[1]],
    #                         np.zeros_like(edges.T[1])])
    # lc = np.stack([u_from, u_to], axis=1)
    # print(lc[:, :, :2])
    # ax2.plot(lc[:, :, 0], lc[:, :, 1], lc[:, :, 2])

    # fig3, ax3 = plt.subplots()
    # ax3.scatter(grad_est[:, 0], grad_est[:, 1])
    # ax3.set(aspect='equal')

    fig4, ax4 = plt.subplots()
    ax4.plot(x[edges.T], y[edges.T], color="blue")
    ax4.set(aspect='equal')

    plt.show()
