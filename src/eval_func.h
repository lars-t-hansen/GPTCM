/* header file for evaluation functions (i.e. log densities) */

#ifndef EVAL_FUNC_H
#define EVAL_FUNC_H

#include <stdio.h>
#include <RcppArmadillo.h>


// typedef std::vector<double> stdvec;

// In the following, most of the pointers are pointers to armadillo storage.  They should instead be
// references to armadillo data structures.  Too much fuss to fix immediately.

// Initializer for per-computation constants

struct initS {
    unsigned int L{0};
    unsigned int N{0};
    unsigned int p{0};
    double kappaA{0.0};
    double kappaB{0.0};
    double kappa{0.0};
    bool invGamma{false};
    double v0Sq{0.0};
    double vSq{0.0};
    double tau0Sq{0.0};
    double w0Sq{0.0};
    bool dirichlet{false};
    const unsigned int *gammaIndicator{nullptr};
    const double *datTheta;
    const unsigned int * datEvent{nullptr};
    const double * datTime{nullptr};
    const double *datProportion{nullptr};
    const double *datProportionConst{nullptr};
};

struct dataS {
    dataS(const initS& init)
        : L(init.L)
        , N(init.N)
        , p(init.p)
        , kappaA(init.kappaA)
        , kappaB(init.kappaB)
        , kappa(init.kappa)
        , invGamma(init.invGamma)
        , vSq(init.vSq)
        , v0Sq(init.v0Sq)
        , tau0Sq(init.tau0Sq)
        , w0Sq(init.w0Sq)
        , dirichlet(init.dirichlet)
        , gammaIndicator(init.gammaIndicator)
        , datTheta(init.datTheta)
        , datEvent(init.datEvent)
        , datTime(init.datTime)
        , datProportion(init.datProportion)
        , datProportionConst(init.datProportionConst)
          // Per-iteration inputs, loop-dependent
        , old_par(NAN)
        , jj(UINT_MAX)
        , l(UINT_MAX)
        , tauSq(NAN)
        , wSq(NAN)
        , datX{nullptr}
        , logMu_l{nullptr}
        , logAlpha_l{nullptr}
        , alphaRowsum{nullptr}
          // Computation state - updated in outer function only
        , datMu{nullptr}
        , alphas{nullptr}
        , weibullS{nullptr}
        , weibullLambda{nullptr}
          // Computation state - updated in inner function too
        , currentPars{nullptr}
        {}

    // constant parameters
    const unsigned int L;
    const unsigned int N;
    const unsigned int p;
    const double kappaA;
    const double kappaB;
    const double kappa;
    const bool invGamma;
    const double vSq;
    const double v0Sq;
    const double tau0Sq;
    const double w0Sq;
    const bool dirichlet;
    const unsigned int * const gammaIndicator;
    const double * const datTheta;
    const unsigned int * const datEvent;
    const double * const datTime;
    const double * const datProportion;
    const double * const datProportionConst;

    // Per-iteration inputs

    // old_par: input to log_dens_betasFull and log_dens_zetasFull, old value of
    // current parameter.
    double old_par;

    // jj: input to various functions, the current outermost iteration index along the p dimension.
    unsigned int jj;

    // l: input to various functions, the current outermost iteration index along the L dimension.
    unsigned int l;

    // tauSq: input to various functions.
    double tauSq;

    // wSq: input to various functions.
    double wSq;

    // datX: input parameters for ... something
    const double *datX;

    // logMu_l: Computed from current parameters, but per-iteration
    const double *logMu_l;

    // logAlpha_l: Computed from current parameters, but per-iteration
    const double *logAlpha_l;

    // alphaRowsum: Computed from alpha, but per-iteration
    const double *alphaRowsum;

    // Computation state - input/output variables.  Data dependencies flow through these and
    // the values escape to the driver code.

    // These are updated only in the outer function, not inside the evaluator!  All updates are
    // visible in arms_gibbs.cpp.  eval_func.cpp will either make a const data structure mapping the
    // memory or will make a copy if it needs to be writable.
    double *datMu;
    double *alphas;
    double *weibullS;
    double *weibullLambda;
    double *currentPars;
};

class EvalFunction
{
public:

    static double log_dens_xis
    (
        double par,
        void *abc_data
    );

    static double log_dens_betas
    (
        double par,
        void *abc_data
    );

    static double log_dens_betasFull
    (
        double par,
        void *abc_data
    );

    static double log_dens_zetas
    (
        double par,
        void *abc_data
    );

    static double log_dens_zetasFull
    (
        double par,
        void *abc_data
    );
    /*
    static double log_dens_phi
    (
        double par,
        void *abc_data
    );
    */

    static double log_dens_kappa
    (
        double par,
        void *abc_data
    );

};

#endif
