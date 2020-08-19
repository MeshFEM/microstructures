#ifndef PRINTABILITY_HH
#define PRINTABILITY_HH

namespace PatternOptimization {
namespace Constraints {

template<size_t N>
struct Printability : Constraint<N> {
    using Base = Constraint<N>;

    using   SField = ScalarField<Real>;
    using Jacobian = typename Base::Jacobian;

    template<class _Iterate>
    Printability(const _Iterate &it)
        : Base(ConstraintType::INEQUALITY)
    {
        if (!it.isParametric())
            throw std::runtime_error("Printability constraints only support parametric optimization.");
        // **Note**: nlopt enforces constraints of the form c_i <= 0, while we
        // defined our printability constraints in the form c_i >= 0.
        // We convert by simply negating the constraint matrix:

        // Note: this means the violation reported should be negated.
        m_C = -it.selfSupportingConstraints();

        const auto &params = it.params();
        assert(size_t(m_C.cols()) == params.size() + 1);

        // Use position map to place points; we need to construct homogeneous
        // param vector.
        m_paramVec = Eigen::VectorXd::Zero(params.size() + 1);
        for (size_t i = 0; i < params.size(); ++i) m_paramVec[i] = params[i];
        m_paramVec[params.size()] = 1.0;

        if (m_C.rows() == 0) throw std::runtime_error("No self-supporting constraints found.");
    }

    virtual size_t dimension() const override { return m_C.rows(); }
    virtual SField  evaluate() const override { return SField(m_C * m_paramVec); }

    // Jacobian of m_c [p] wrt p is just m_C with last column removed.
    //                 [1]
    virtual Jacobian jacobian(const std::vector<VectorField<Real, N>> &/* bdrySVels */) const override {
        assert(m_C.cols() > 1);
        return Jacobian(m_C.block(0, 0, m_C.rows(), m_C.cols() - 1));
    }

    virtual void writeDescription(std::ostream &os, const std::string &name) const override {
        // Positive values mean violation
        os << "max printability violation: " << (m_C * m_paramVec).maxCoeff() << std::endl;
        Base::writeDescription(os, name);
    }

    virtual ~Printability() { }

private:
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> m_C;
    Eigen::VectorXd                                     m_paramVec;
};

// Configuration to be applied by iterate factory
template<class _Sim>
struct IFConfigPrintability : public IFConfig {
    static constexpr size_t N = _Sim::N;

    template<class _Iterate>
    void configIterate(const std::unique_ptr<_Iterate> &it, ObjectiveTermNormalizations &/* normalizations */) const {
        static_assert(_Iterate::_N == N, "Mismatch in problem dimensions.");
        auto pc = Future::make_unique<Printability<N>>(*it);
        it->addConstraint("Printability", std::move(pc));
    }
};

}} // end namespace PatternOptimization::Constraints


#endif /* end of include guard: PRINTABILITY_HH */
