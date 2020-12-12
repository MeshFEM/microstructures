////////////////////////////////////////////////////////////////////////////////
#include "jacobians.h"
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace micro_quadfoam {

// Let Q = [v0, v1, v2, v3] be a quad.
// We compute the best Jacobian J such that J maps the reference square [-1,1]²
// onto Q. To remove translation ambiguities, let us assume that J maps the centroid
// of the unit square onto the centroid of Q. Now, we have an overdetermined
// system of equations that we can solve via linear least-squares.
//
//               ┌       ┐
//               │ ex0 … │
//               │ ey0 … │
//             x └       ┘
//     ┌       ┐             ┌       ┐
// J = │ x1 x2 │           = │ ux0 … │
//     │ y1 y2 │             │ uy0 … │
//     └       ┘             └       ┘
//
// We rewrite this system as a matrix-vector product where J is written a column-vector.
//
//                     ┌    ┐
//                     │ x1 │
//                     │ x2 │
//                     │ y1 │
//                     │ y2 │
//                   x └    ┘
// ┌                 ┐           ┌     ┐
// │ ex0 ey0         │         = │ ux0 │
// │         ex0 ey0 │           │ uy0 │
// │ ex1 ey1         │           │ ux1 │
// │        …        │           │  …  │
// └                 ┘           └     ┘
//

void jacobians(const Eigen::MatrixXd &V, const Eigen::MatrixXi &Q, Eigen::MatrixXd &J, Eigen::MatrixXd &P)
{
	assert(V.cols() <= 3);
	Eigen::Matrix<double, 8, 4> A;
	A << -1, -1, 0, 0,
		0, 0, -1, -1,
		1, -1, 0, 0,
		0, 0, 1, -1,
		1, 1, 0, 0,
		0, 0, 1, 1,
		-1, 1, 0, 0,
		0, 0, -1, 1;
	Eigen::Matrix<double, 8, 1> b;
	Eigen::Matrix<double, 4, 1> x;
	J.resize(Q.rows(), 4);
	P.resize(Q.rows() * 4, 3);
	P.setZero();
	for (int q = 0; q < Q.rows(); ++q) {
		// Barycenter
		Eigen::RowVector2d c(0, 0);
		for (int lv = 0; lv < 4; ++lv) { c += V.row(Q(q, lv)).head<2>(); }
		c /= (double) Q.cols();

		// Normalize quad
		Eigen::Matrix<double, 4, 2> v;
		v <<
			V.row(Q(q, 0)).head<2>() - c,
			V.row(Q(q, 1)).head<2>() - c,
			V.row(Q(q, 2)).head<2>() - c,
			V.row(Q(q, 3)).head<2>() - c;

		// Compute Jacobian
		b.segment<2>(0) = v.row(0);
		b.segment<2>(2) = v.row(1);
		b.segment<2>(4) = v.row(2);
		b.segment<2>(6) = v.row(3);
		x = A.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(b);
		// The determinant should only become negative if we have inverted quads to begin with?
		if (x[0]*x[3] - x[1]*x[2] < 0) {
			std::cerr << "Jacobian with negative determinant." << std::endl;
			assert(false);
		}
		J.row(q) = x.transpose();
		b = A*x;
		P.row(4*q+0).head<2>() = b.segment<2>(0).transpose() + c;
		P.row(4*q+1).head<2>() = b.segment<2>(2).transpose() + c;
		P.row(4*q+2).head<2>() = b.segment<2>(4).transpose() + c;
		P.row(4*q+3).head<2>() = b.segment<2>(6).transpose() + c;
	}
}

////////////////////////////////////////////////////////////////////////////////

} // namespace micro_quadfoam
