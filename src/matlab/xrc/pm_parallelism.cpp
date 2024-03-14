////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   This file is part of ParaMonte: Parallel Monte Carlo and Machine Learning library.
//
//   Copyright (C) 2012-present, The Computational Data Science Lab
//
//   https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "mex.hpp"
#include "mexAdapter.hpp"
extern "C" {
    int32_t getImageCountMPI(void);
}
class MexFunction : public matlab::mex::Function {
    matlab::data::ArrayFactory factory;
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        // Function implementation
        checkArguments(outputs, inputs);
        //matlab::data::TypedArray<double> imageCount = getImageCountMPI();
        double imageCount = (double) getImageCountMPI();
        outputs[0] = factory.createScalar(imageCount);
    }
    void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        matlab::data::ArrayFactory factory;
        if (inputs.size() != 1) {
            matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({factory.createScalar("One inputs required")}));
        }
        //if (inputs[0].getNumberOfElements() != 1) {
        //    matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({factory.createScalar("Input function name must be a scalar")}));
        //}
        if (inputs[0].getType() != matlab::data::ArrayType::CHAR) {
            matlabPtr->feval(u"error", 0, std::vector<matlab::data::Array>({factory.createScalar("Input multiplier must be a noncomplex scalar double")}));
        }
    }
};