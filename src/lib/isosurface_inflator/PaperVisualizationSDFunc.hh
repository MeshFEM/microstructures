#ifndef PAPERVISUALIZATIONSDFUNC_HH
#define PAPERVISUALIZATIONSDFUNC_HH

#include "SignedDistance.hh"
#include "Joint.hh"
#include <stdexcept>

namespace SD = SignedDistance;

struct PaperVisualizationSDFunc : public SignedDistanceRegion<3> {
    using Real = double;
    enum Mode { EDGE1, EDGE2, EDGE3, SHORT_EDGE1, BLEND_FULL, BLEND_HULL, HULL, JOINT_1, JOINT_2, JOINT_UNION, REDUCED_SMOOTH, REDUCED_SMOOTH_BULGE };
    Mode mode = EDGE1;
    SD::Primitives::InflatedEdge<double> edge1, edge2, edge3,
                                         edge4, edge5, edge6,
                                         shortEdge1;
    std::unique_ptr<Joint<double>> joint1, joint2,
                                   standalone_joint;
    double smoothParam;

    std::vector<Point3D> centers, centers2;
    std::vector<double> radii;

    PaperVisualizationSDFunc(double r1, double r2, double r3, double r4, double s)
    {
        centers.push_back(Point3D(0.0, 0.0, 0.0));
        centers.push_back(Point3D(0.4, 0.0, 0.0));
        centers.push_back(Point3D(0.0, 0.3, 0.0));
        centers.push_back(Point3D(0.0, 0.0, 0.3));

        centers2.push_back(Point3D(0.4, 0.0, 0.0));
        centers2.push_back(Point3D(0.0, 0.0, 0.0));
        centers2.push_back(Point3D(0.4, 0.3, 0.0));
        centers2.push_back(Point3D(0.4, 0.0, 0.3));

        setParameters(r1, r2, r3, r4, s);
    }

    void setParameters(double r1, double r2, double r3, double r4, double s) {
        radii = { r1, r2, r3, r4 };
        edge1.set(centers[0], centers[1], radii[0], radii[0]);
        edge2.set(centers[0], centers[2], radii[0], radii[1]);
        edge3.set(centers[0], centers[3], radii[0], radii[2]);
        shortEdge1.set(centers[0], 0.75 * centers[1], radii[0], radii[2]);
        smoothParam = s;

        edge4.set(centers2[0], centers2[1], radii[0], radii[0]);
        edge5.set(centers2[0], centers2[2], radii[0], radii[1]);
        edge6.set(centers2[0], centers2[3], radii[0], radii[2]);


        joint1 = Future::make_unique<Joint<double>>(centers,  radii, s, JointBlendMode::HULL);
        joint2 = Future::make_unique<Joint<double>>(centers2, radii, s, JointBlendMode::HULL);

        auto standaloneHullCenters = centers;
        standaloneHullCenters[1][0] *= 0.75;
        standalone_joint = Future::make_unique<Joint<double>>(standaloneHullCenters, radii, s, JointBlendMode::HULL);
    }

    virtual void boundingSphere(Point3D &c, double &r) const override {
        c.setZero();
        r = 3.0;
    }

    virtual const BBox<PointNd<3>> & boundingBox() const override {
    	static auto box = BBox<Point3D>(Point3D(-1, -1, -1), Point3D(2, 2, 2));
    	return box;
    }

    virtual Real signedDistance(const Point3D &p) const override {
        const double maxOverlapSmoothingAmt = 0.02;

        if (mode == Mode::EDGE1) return edge1.signedDistance(p);
        if (mode == Mode::EDGE2) return edge2.signedDistance(p);
        if (mode == Mode::EDGE3) return edge3.signedDistance(p);
        if (mode == Mode::SHORT_EDGE1) return shortEdge1.signedDistance(p);

        // Blended version of the stand-alone joint
        if ((mode == Mode::BLEND_FULL) || (mode == Mode::BLEND_HULL)) {
            double smoothAmt = (mode == Mode::BLEND_FULL) ? smoothParam
                                                          : standalone_joint->smoothingAmt(p);
            std::vector<double> dists = { shortEdge1.signedDistance(p), edge2.signedDistance(p), edge3.signedDistance(p) };
            return SD::exp_smin_reparam_accurate(dists, smoothAmt);
        }
        if (mode == Mode::HULL)       return standalone_joint->blendingHull().signedDistance(p);

        if (mode == Mode::JOINT_1) {
            double smoothAmt = (mode == Mode::BLEND_FULL) ? smoothParam
                                                          : joint1->smoothingAmt(p);
            std::vector<double> dists = { edge1.signedDistance(p),
                                          edge2.signedDistance(p),
                                          edge3.signedDistance(p) };
            return SD::exp_smin_reparam_accurate(dists, smoothAmt);
        }
        if (mode == Mode::JOINT_2) {
            double smoothAmt = (mode == Mode::BLEND_FULL) ? smoothParam
                                                          : joint2->smoothingAmt(p);
            std::vector<double> dists = { edge4.signedDistance(p),
                                          edge5.signedDistance(p),
                                          edge6.signedDistance(p) };
            return SD::exp_smin_reparam_accurate(dists, smoothAmt);
        }
        if (mode == Mode::JOINT_UNION) {
            std::vector<double> dists1 = { edge1.signedDistance(p), edge2.signedDistance(p), edge3.signedDistance(p) };
            std::vector<double> dists2 = { edge4.signedDistance(p), edge5.signedDistance(p), edge6.signedDistance(p) };

            double hardDist1 = *std::min_element(dists1.begin(), dists1.end()),
                   hardDist2 = *std::min_element(dists2.begin(), dists2.end());

            double smoothDist1 = SD::exp_smin_reparam_accurate(dists1, joint1->smoothingAmt(p)),
                   smoothDist2 = SD::exp_smin_reparam_accurate(dists2, joint2->smoothingAmt(p));

            double smoothEffect1 = hardDist1 - smoothDist1,
                   smoothEffect2 = hardDist2 - smoothDist2;
            assert(smoothEffect1 >= -1e-9);
            assert(smoothEffect2 >= -1e-9);
            smoothEffect1 = std::max(smoothEffect1, 0.0);
            smoothEffect2 = std::max(smoothEffect2, 0.0);

            double meanSmoothEffectSq = smoothEffect1 * smoothEffect2;

            double overlapSmoothAmt = maxOverlapSmoothingAmt * tanh(1000.0 * meanSmoothEffectSq);
            return SD::exp_smin_reparam_accurate(smoothDist1, smoothDist2, overlapSmoothAmt);
        }

        if (mode == Mode::REDUCED_SMOOTH) {
            std::vector<double> dists1 = { edge1.signedDistance(p), edge2.signedDistance(p), edge3.signedDistance(p) };
            std::vector<double> dists2 = { edge4.signedDistance(p), edge5.signedDistance(p), edge6.signedDistance(p) };

            double hardDist1 = *std::min_element(dists1.begin(), dists1.end()),
                   hardDist2 = *std::min_element(dists2.begin(), dists2.end());

            double smoothDist1 = SD::exp_smin_reparam_accurate(dists1, joint1->smoothingAmt(p) / 3),
                   smoothDist2 = SD::exp_smin_reparam_accurate(dists2, joint2->smoothingAmt(p) / 3);

            double smoothEffect1 = hardDist1 - smoothDist1,
                   smoothEffect2 = hardDist2 - smoothDist2;
            assert(smoothEffect1 >= -1e-9);
            assert(smoothEffect2 >= -1e-9);
            smoothEffect1 = std::max(smoothEffect1, 0.0);
            smoothEffect2 = std::max(smoothEffect2, 0.0);

            double meanSmoothEffectSq = smoothEffect1 * smoothEffect2;

            double overlapSmoothAmt = maxOverlapSmoothingAmt * tanh(1000.0 * meanSmoothEffectSq);
            return SD::exp_smin_reparam_accurate(smoothDist1, smoothDist2, overlapSmoothAmt);
        }

        if (mode == Mode::REDUCED_SMOOTH_BULGE) {
            std::vector<double> dists1 = { edge1.signedDistance(p), edge2.signedDistance(p), edge3.signedDistance(p) };
            std::vector<double> dists2 = { edge4.signedDistance(p), edge5.signedDistance(p), edge6.signedDistance(p) };
            double smoothDist1 = SD::exp_smin_reparam_accurate(dists1, joint1->smoothingAmt(p) / 3),
                   smoothDist2 = SD::exp_smin_reparam_accurate(dists2, joint2->smoothingAmt(p) / 3);
            return SD::exp_smin_reparam_accurate(smoothDist1, smoothDist2, maxOverlapSmoothingAmt);
        }

        throw std::runtime_error("Unimplemented");
    }

    virtual bool isInside(const Point3D &p) const override {
        return signedDistance(p) <= 0;
    }
};

#endif /* end of include guard: PAPERVISUALIZATIONSDFUNC_HH */
