#include "gradient_descent.hh"
#include <MeshFEM/EdgeFields.hh>

void optimize_gd(ScalarField<Real> &params,
        const PatternOptimization::BoundConstraints &bds,
        PatternOptimization::IterateManagerBase &im,
        const PatternOptimization::OptimizerConfig &oconfig,
        const std::string &outPath)
{
    im.setOutPath(outPath);

    for (size_t i = 0; i < oconfig.niters; ++i) {
        // To std::vector<double>
        std::vector<Real> paramsCopy;
        params.getFlattened(paramsCopy);

        auto &it = im.get(params.size(), params.data());
        it.writeDescription(std::cout);
        std::cout << std::endl;

        // Report
        im.updateAndReport(paramsCopy);

        params += it.steepestDescent() * oconfig.gd_step;

        // Apply bound constraints
        for (size_t p = 0; p < im.numParameters(); ++p) {
            if (bds.hasUpperBound.at(p)) params[p] = std::min(params[p], bds.upperBound.at(p));
            if (bds.hasLowerBound.at(p)) params[p] = std::max(params[p], bds.lowerBound.at(p));
        }

        if (outPath != "") it.writeMeshAndFields(outPath + "_" + std::to_string(i));
    }
}

void validatePerturbation(ScalarField<Real> &perturbation) {
    double maxParameterChange = 0.01;
    for (size_t i = 0; i < perturbation.domainSize(); i++) {
        if (perturbation(i) > maxParameterChange)
            perturbation(i) = maxParameterChange;
        else if (perturbation(i) < -maxParameterChange)
            perturbation(i) = -maxParameterChange;
    }
}

// gradient descent with custom line search
void optimize_gd_smartstep(ScalarField<Real> &params,
                           const PatternOptimization::BoundConstraints &bds,
                           PatternOptimization::IterateManagerBase &im,
                           const PatternOptimization::OptimizerConfig &oconfig,
                           const std::string &outPath)
{
    double step_size = oconfig.gd_step;
    int success_iterations = 0;
    int failure_iterations = 0;
    std::vector<Real> oldParamsCopy;
    std::vector<Real> steepestDescentCopy;
    double val = std::numeric_limits<Real>::max();
    double bestVal = std::numeric_limits<Real>::max();
    double previousVal = std::numeric_limits<Real>::max();
    double minimumStep = 1e-10;
    size_t maxTries = 100;
    size_t tries = 0;

    im.setOutPath(outPath);

    for (size_t i = 0; i < oconfig.niters; ++i) {
        //std::cout << "[GD] New gradient descent iteration" << std::endl;
        if (step_size < minimumStep) {
            if (tries < maxTries) {
                step_size = oconfig.gd_step;
                tries++;
            }
            else {
                return;
            }
        }

        try {
            // To std::vector<double>
            std::vector<Real> paramsCopy;
            params.getFlattened(paramsCopy);

            auto &it = im.get(params.size(), params.data());

            val = it.evaluate();
            if (val <= (previousVal + 1e-10)) {
                //std::cout << "[GD] Previous val: " << previousVal << std::endl;
                previousVal = val;
                it.steepestDescent().getFlattened(steepestDescentCopy);
                //std::cout << "[GD] New val: " << previousVal << std::endl;
                params.getFlattened(oldParamsCopy);
                if (val < bestVal) {
                    bestVal = val;
                    //std::cout << "[GD] Best iteration: " << i << std::endl;
                    //std::cout << "[GD] Best value: " << val << std::endl;
                }
                success_iterations++;
                failure_iterations=0;
            }
            else {
                success_iterations=0;
                failure_iterations++;
            }

            // Report
            im.updateAndReport(paramsCopy);

            ScalarField<Real> perturbation = it.steepestDescent() * step_size;
            params += perturbation;

            // Apply bound constraints
            for (size_t p = 0; p < im.numParameters(); ++p) {
                if (bds.hasUpperBound.at(p)) params[p] = std::min(params[p], bds.upperBound.at(p));
                if (bds.hasLowerBound.at(p)) params[p] = std::max(params[p], bds.lowerBound.at(p));
            }
        }
        catch (...) {
            //std::cout << "[GD] Exploded!" << std::endl;
            step_size = step_size / 2;
            success_iterations = 0;
            failure_iterations = 0;

            // Recover!
            ScalarField<Real> oldParams(oldParamsCopy);
            ScalarField<Real> steepestDescent(steepestDescentCopy);
            ScalarField<Real> newPerturbation = steepestDescent * step_size;
            params = oldParams + newPerturbation;

            // Apply bound constraints
            for (size_t p = 0; p < im.numParameters(); ++p) {
                if (bds.hasUpperBound.at(p)) params[p] = std::min(params[p], bds.upperBound.at(p));
                if (bds.hasLowerBound.at(p)) params[p] = std::max(params[p], bds.lowerBound.at(p));
            }
        }

        if (success_iterations >= 10) {
            step_size *= 2;
            success_iterations = 0;
            failure_iterations = 0;
        }
        else if (failure_iterations >= 1) {
            //std::cout << "[GD] Failure" << std::endl;
            step_size = step_size / 2;
            success_iterations = 0;
            failure_iterations = 0;

            // Recover!
            ScalarField<Real> oldParams(oldParamsCopy);
            ScalarField<Real> steepestDescent(steepestDescentCopy);
            ScalarField<Real> newPerturbation = steepestDescent * step_size;
            params = oldParams + newPerturbation;

            // Apply bound constraints
            for (size_t p = 0; p < im.numParameters(); ++p) {
                if (bds.hasUpperBound.at(p)) params[p] = std::min(params[p], bds.upperBound.at(p));
                if (bds.hasLowerBound.at(p)) params[p] = std::max(params[p], bds.lowerBound.at(p));
            }
        }

        //std::cout << "[GD] Step size is now: " << step_size << std::endl;
        //std::cout << "[GD] End of iteration val: " << previousVal << std::endl;
    }
}