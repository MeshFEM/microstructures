#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace micro_quadfoam {

//
// @brief         { Return the closest 2D Jacobian associated to each input quads, potentially
//                shifting quad indices to avoid flips in the Jacobian. }
//
// @param[in]     V     { #V x 3 input points }
// @param[in]     Q     { #Q x 4 input quads }
// @param[out]    J     { #Q x 4 output Jacobian x1 x2 y1 y2 }
// @param[out]    P     { 4*#Q x 3 output parallelogram corners }
//
//     ┌       ┐
// J = │ x1 x2 │
//     │ y1 y2 │
//     └       ┘
//
void jacobians(const Eigen::MatrixXd &V, const Eigen::MatrixXi &Q,
	Eigen::MatrixXd &J, Eigen::MatrixXd &P);

} // namespace micro_quadfoam
