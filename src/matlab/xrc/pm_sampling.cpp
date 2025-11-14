// pm_sampling_mex.cpp  — C++ MEX modernization of pm_sampling.c
//
// Matches semantics of the original C MEX file while using the C++ MEX API.
// - Input args (same order/requirements):
//   0) method   : char row vector (expects "ParaDRAM")
//   1) fun      : function_handle
//   2) ndim     : scalar integer
//   3) input    : char row vector (passed through to runParaDRAM)
// - No outputs; errors are thrown like the original.
//
// Build:  mex -v -R2018a CXXFLAGS="$CXXFLAGS -std=c++17" pm_sampling_mex.cpp

#include <cstdint>
#include <string>
#include <vector>
#include <stdexcept>

#include "mex.hpp"
#include "mexAdapter.hpp"

using matlab::mex::ArgumentList;
using matlab::data::Array;
using matlab::data::ArrayFactory;
using matlab::data::ArrayType;
using matlab::data::CharArray;
using matlab::data::TypedArray;
using matlab::engine::MATLABEngine;

///////////////////////////////////////////////////////////////////////////////////////////////////
// External sampler symbol (name alias kept to mirror the original C file)
///////////////////////////////////////////////////////////////////////////////////////////////////
#define runParaDRAM runParaDRAMD

#if OMP_ENABLED
extern "C" int32_t runParaDRAM( double(*getLogFunc)( double logFuncState[]
                                                   , int32_t ndim
                                                   , int32_t njob
                                                   , double *avgTimePerFunCallComp
                                                   , double *avgTimePerFunCallComm)
                              , const int32_t ndim
                              , const char* input );
#else
extern "C" int32_t runParaDRAM( double(*getLogFunc)( double state[]
                                                   , int32_t ndim )
                              , const int32_t ndim
                              , const char* input );
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
// MexFunction class
///////////////////////////////////////////////////////////////////////////////////////////////////
class MexFunction : public matlab::mex::Function {
public:
    // Store engine and function handle so static callbacks can use them
    static MexFunction* self;

    MexFunction() { self = this; }
    ~MexFunction() override = default;

    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        engine_ = getEngine();
        ArrayFactory f;

        // === Validate counts ===
        if (!outputs.empty()) {
            err(u"Internal ParaMonte MATLAB library error occurred: Too many output arguments.");
        }
        if (inputs.size() != 4) {
            err(u"Internal ParaMonte MATLAB library error occurred: input variable mismatch.");
        }

        // === Parse input #0: method (row char) ===
        if (inputs[0].getType() != ArrayType::CHAR || inputs[0].getDimensions().size() != 2 ||
            inputs[0].getDimensions()[0] != 1) {
            err(u"Input must be a row vector of characters.");
        }
        CharArray methodChar = inputs[0];
        const std::u16string method16 = methodChar.toUTF16();
        // We only accept "ParaDRAM" like the original
        if (method16.size() < 8 || method16.substr(0,8) != u"ParaDRAM") {
            err(u"Internal ParaMonte MATLAB library error occurred: Invalid input sampling method.");
        }

        // === Parse input #1: MATLAB function handle ===
        if (inputs[1].getType() != ArrayType::HANDLE_OBJECT_REF) {
            err(u"The second input argument must be a function handle.");
        }
        funcHandle_ = inputs[1]; // keep a reference-managed handle

        // === Parse input #2: ndim (scalar) ===
        if (inputs[2].getType() != ArrayType::DOUBLE || inputs[2].getNumberOfElements() != 1) {
            err(u"Internal ParaMonte MATLAB library error occurred: Input #2 (ndim) must be a scalar.");
        }
        const int32_t ndim = static_cast<int32_t>( static_cast<TypedArray<double>>(inputs[2])[0] );

        // === Parse input #3: input string (row char) ===
        if (inputs[3].getType() != ArrayType::CHAR || inputs[3].getDimensions().size() != 2 ||
            inputs[3].getDimensions()[0] != 1) {
            err(u"Internal ParaMonte MATLAB library error occurred: Input #3 must be a row char vector.");
        }
        CharArray inChar = inputs[3];
        inputStr_ = inChar.toAscii();        // ASCII/UTF-8; mirrors mxArrayToString behavior
        const char* c_input = inputStr_.c_str();

        // === Call external sampler ===
        int32_t stat = 0;
#if OMP_ENABLED
        stat = runParaDRAM(&MexFunction::getLogFuncOMP, ndim, c_input);
#else
        stat = runParaDRAM(&MexFunction::getLogFuncScalar, ndim, c_input);
#endif
        if (stat != 0) {
            err(u"Mex:ParaMonte: Runtime Error Occurred.");
        }
    }

private:
    // Convenience error helper
    [[noreturn]] void err(std::u16string msg) {
        ArrayFactory f;
        engine_->feval(u"error", 0, std::vector<Array>{ f.createScalar(msg) });
        throw std::runtime_error("unreachable");
    }

    // === Static callbacks that match runParaDRAM expectations ===
#if OMP_ENABLED
    // getLogFunc for OMP-enabled multi-job path
    static double getLogFuncOMP(double logFuncState[], int32_t ndim, int32_t njob,
                                double* avgTimePerFunCallComp, double* avgTimePerFunCallComm)
    {
        // logFuncState layout: for each job j, slot 0 is for logf output, slots 1..ndim are the state
        // We need to build an ndim-by-njob matrix of states (without the first column per job).
        ArrayFactory f;
        auto eng = self->engine_;

        // Construct state matrix (ndim x njob)
        matlab::data::TypedArray<double> state = f.createArray<double>({ static_cast<size_t>(ndim),
                                                                         static_cast<size_t>(njob) });
        {
            // Fill column-wise like original code
            // original loops: for (ijob) for (idim = ndimp1*ijob+1 .. ndimp1*(ijob+1)-1)
            const int32_t ndimp1 = ndim + 1;
            size_t idx = 0;
            for (int32_t j = 0; j < njob; ++j) {
                const int32_t base = ndimp1 * j;
                for (int32_t d = 1; d < ndimp1; ++d) {
                    // MATLAB column-major order: state(d, j+1)
                    // We can set by linear iterator in createArray order (column-major)
                    state.begin()[idx++] = logFuncState[base + d];
                }
            }
        }

        // Call feval(funHandle, state) expecting 3 outputs:
        //   1) logf (1 x njob double)
        //   2) avgTimePerFunCallComp (scalar)
        //   3) avgTimePerFunCallComm (scalar)
        std::vector<Array> in{ self->funcHandle_, state };
        auto out = eng->feval(u"feval", 3, in);

        // Parse outputs
        {
            // out[0] : vector of length njob (row or column); read linearly
            TypedArray<double> logf = out[0];
            // Write back to logFuncState: position 0 in each job’s block
            for (int32_t j = 0; j < njob; ++j) {
                logFuncState[j * (ndim + 1)] = logf.begin()[static_cast<size_t>(j)];
            }
        }

        *avgTimePerFunCallComp = static_cast<TypedArray<double>>(out[1])[0];
        *avgTimePerFunCallComm = static_cast<TypedArray<double>>(out[2])[0];

        // The original returned a dummy double (mold = 0). Preserve that.
        return 0.0;
    }
#else
    // getLogFunc for scalar (non-OMP) path
    static double getLogFuncScalar(double stateBuf[], int32_t ndim)
    {
        ArrayFactory f;
        auto eng = self->engine_;

        // Build ndim x 1 column vector
        TypedArray<double> state = f.createArray<double>({ static_cast<size_t>(ndim),
                                                           static_cast<size_t>(1) });
        // Copy data (column vector)
        {
            auto it = state.begin();
            for (int32_t i = 0; i < ndim; ++i) *it++ = stateBuf[i];
        }

        // feval(funHandle, state) → 1 output (scalar logf)
        std::vector<Array> in{ self->funcHandle_, state };
        auto out = eng->feval(u"feval", 1, in);

        double logf = static_cast<TypedArray<double>>(out[0])[0];
        return logf;
    }
#endif

private:
    std::shared_ptr<MATLABEngine> engine_;
    Array funcHandle_;      // MATLAB function handle (kept alive for callbacks)
    std::string inputStr_;  // backing storage to keep c_str() valid during sampler call
};

// Define the static self pointer
MexFunction* MexFunction::self = nullptr;
