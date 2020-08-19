////////////////////////////////////////////////////////////////////////////////
// Isometries.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Represents an isometry, intended to be used as an element of a symmetry
//      group. These can also be used to map points.
//      These isometries are not necessarily orientation-preserving.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  06/26/2015 18:01:53
////////////////////////////////////////////////////////////////////////////////
#ifndef ISOMETRIES_HH
#define ISOMETRIES_HH
#include "InflatorTypes.hh"
#include <MeshFEM/Future.hh>
#include <memory>
#include <iostream>
#include <cassert>
#include <string>

namespace Symmetry {
    enum class Axis : unsigned int { X = 0, Y = 1, Z = 2, ANY = 255 };
    inline std::string axisName(Axis a) {
        switch(a) {
            case Axis::X: return "X";
            case Axis::Y: return "Y";
            case Axis::Z: return "Z";
            case Axis::ANY: return "ANY";
        }
        return "Unknown";
    }
}

// Elements of a symmetry group
struct Isometry {
    struct Operation;

    // Default constructor: identity isometry
    Isometry() { }
    Isometry(std::unique_ptr<const Operation> &&op) { operations.emplace_back(std::move(op)); }

    Isometry(Isometry &&b) : operations(std::move(b.operations)) { }
    Isometry(const Isometry &b) {
        operations.reserve(b.operations.size());
        for (auto &op : b.operations)
            operations.emplace_back(op->clone());
    }

    Isometry &operator=(Isometry &&b) { operations = std::move(b.operations); return *this; }
    Isometry &operator=(const Isometry &b) {
        operations.clear();
        operations.reserve(b.operations.size());
        for (auto &op : b.operations)
            operations.emplace_back(op->clone());
        return *this;
    }

    bool isIdentity() const { return operations.empty(); }
    bool hasTranslation() const {
        for (const auto &op : operations)
            if (op->isTranslation()) return true;
        return false;
    }
    bool hasPermutation() const {
        for (const auto &op : operations)
            if (op->isPermutation()) return true;
        return false;
    }
    bool hasReflection() const {
        for (const auto &op : operations)
            if (op->isReflection()) return true;
        return false;
    }
    bool hasReflection(const Symmetry::Axis axis = Symmetry::Axis::ANY) const {
        for (const auto &op : operations)
            if (op->isReflection(axis)) return true;
        return false;
    }
    bool affectsAxis(const Symmetry::Axis axis) const {
        for (const auto &op : operations)
            if (op->affectsAxis(axis)) return true;
        return false;
    }

    // Compose "op" on the right (op will be performed before operations[]).
    Isometry compose(std::unique_ptr<const Operation> &&op) const {
        Isometry result(*this);
        result.operations.emplace_back(std::move(op));
        return result;
    }

    Isometry compose(const std::unique_ptr<const Operation> &op) const { return compose(op->clone()); }
    Isometry compose(const Operation *op) const { return compose(op->clone()); }

    // Performed in right-to-left: operations[n - 1], ..., operations[0]
    template<typename Real>
    Point3<Real> apply(Point3<Real> p) const {
        for (size_t i = operations.size(); i > 0; --i) {
            operations[i - 1]->apply(p);
        }
        return p;
    }

    // Compose the isometry with the map from pattern parameters to a point.
    // posMap is assumed to have nParams + 1 columns, where the last column
    // holds a constant translation.
    // @return  transformed copy of posMap
    Eigen::Matrix3Xd xformMap(Eigen::Matrix3Xd posMap) const {
        for (size_t i = operations.size(); i > 0; --i)
            operations[i - 1]->xformMap(posMap);
        return posMap;
    }

    void print(std::ostream &os) const {
        bool first = true;
        if (operations.size() == 0) os << "identity";
        for (auto &op : operations) {
            if (!first) os << " * ";
            op->print(os);
            first = false;
        }
    }

    struct Operation {
        template<typename Real>
        void apply(Point3<Real> &p) const { return applyOperation(this, p); }
        virtual void print(std::ostream &os) const = 0;
        virtual std::unique_ptr<Operation> clone() const = 0;
        virtual bool isTranslation() const { return false; }
        virtual bool isReflection(const Symmetry::Axis /*axis*/ = Symmetry::Axis::ANY) const { return false; }
        virtual bool isPermutation() const { return false; }
        virtual bool affectsAxis(const Symmetry::Axis /*axis*/) const = 0;
        virtual void xformMap(Eigen::Matrix3Xd &posMap) const = 0;
        virtual ~Operation() { }
    };

    struct Translation : public Operation {
        Translation(double x, double y, double z) : t(x, y, z) { }
        virtual ~Translation() { }
        virtual void print(std::ostream &os) const override { os << "t(" << t[0] << ", " << t[1] << ", " << t[2] << ")"; }
        virtual std::unique_ptr<Operation> clone() const override { return Future::make_unique<Translation>(*this); }
        virtual bool isTranslation() const override { return true; }
        virtual bool affectsAxis(const Symmetry::Axis axis) const override { return t[size_t(axis)] != 0; }
        virtual void xformMap(Eigen::Matrix3Xd &posMap) const override {
            assert(posMap.cols() > 1);
            posMap.col(posMap.cols() - 1) += t;
        }
        Vector3<double> t;
    };

    struct Reflection : public Operation {
        Reflection(Symmetry::Axis a) : a(a) { }
        virtual ~Reflection() { }
        virtual void print(std::ostream &os) const override { os << "r(" << axisName(a) << ")"; }
        virtual std::unique_ptr<Operation> clone() const override { return Future::make_unique<Reflection>(*this); }
        virtual bool isReflection(const Symmetry::Axis axis = Symmetry::Axis::ANY) const override {
            return (axis == Symmetry::Axis::ANY) || (axis == a);
        }
        virtual bool affectsAxis(const Symmetry::Axis axis) const override { return axis == a; }
        virtual void xformMap(Eigen::Matrix3Xd &posMap) const override {
            posMap.row(static_cast<unsigned int>(a)) *= -1.0;
        }
        Symmetry::Axis a;
    };

    struct Permutation : public Operation {
        Permutation(Symmetry::Axis a1, Symmetry::Axis a2) : a1(a1), a2(a2) { }
        virtual ~Permutation() { }
        virtual void print(std::ostream &os) const override { os << "p(" << axisName(a1) << ", " << axisName(a2) << ")"; }
        virtual std::unique_ptr<Operation> clone() const override { return Future::make_unique<Permutation>(*this); }
        virtual bool isPermutation() const override { return true; }
        virtual bool affectsAxis(const Symmetry::Axis axis) const override { return (axis == a1) || (axis == a2); }
        virtual void xformMap(Eigen::Matrix3Xd &posMap) const override {
            posMap.row(static_cast<unsigned int>(a1)).swap(posMap.row(static_cast<unsigned int>(a2)));
        }
        Symmetry::Axis a1, a2;
    };

    // Unfortunately virtual template methods are not allowed, so we're forced
    // to implement polymorphism with a dynamic cast in applyOperation()
    template<typename Real>
    static void applyOperation(const Operation *op, Point3<Real> &p) {
        if (const Translation *trans = dynamic_cast<const Translation *>(op)) {
            for (size_t i = 0; i < 3; ++i) p[i] += trans->t[i];
        }
        else if (const Reflection *refl = dynamic_cast<const Reflection *>(op)) {
            p[static_cast<unsigned int>(refl->a)] *= -1;
        }
        else if (const Permutation *perm = dynamic_cast<const Permutation *>(op)) {
            std::swap(p[static_cast<unsigned int>(perm->a1)],
                      p[static_cast<unsigned int>(perm->a2)]);
        }
        else { assert("Unrecognized operation!"); } // Impossible
    }

    static std::unique_ptr<Translation> translation(double x, double y, double z)         { return Future::make_unique<Translation>(x, y, z); }
    static std::unique_ptr< Reflection>  reflection(Symmetry::Axis a)                     { return Future::make_unique< Reflection>(a); }
    static std::unique_ptr<Permutation> permutation(Symmetry::Axis a1, Symmetry::Axis a2) { return Future::make_unique<Permutation>(a1, a2); }

    std::vector<std::unique_ptr<const Operation>> operations;
};

inline std::ostream &operator<<(std::ostream &os, const Isometry &i) {
    i.print(os);
    return os;
}

#endif /* end of include guard: ISOMETRIES_HH */
