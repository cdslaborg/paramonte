!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This module contains classes and procedures for computing the properties related to the covariance matrices of a random sample.<br>
!>
!>  \details
!>
!>  Variance
!>  ========
!>
!>  Variance is the squared deviation from the mean of a random variable.<br>
!>  The variance is also often defined as the square of the standard deviation.<br>
!>  Variance is a measure of dispersion, meaning it is a measure of how far a set of numbers is spread out from their average value.<br>
!>  It is the second central moment of a distribution, and the covariance of the random variable with itself.<br>
!>  It is frequently represented by \f$\Sigma\f$, \f$\sigma^2\f$, \f$s^2\f$, \f$\up{Var}(X)\f$, or \f$\mathbb{V}(X)\f$.<br>
!>
!>  Variance as a measure of dispersion
!>  -----------------------------------
!>
!>  An advantage of variance as a measure of dispersion is that it is more amenable to algebraic manipulation
!>  than other measures of dispersion such as the expected absolute deviation.<br>
!>  For example, the variance of a sum of uncorrelated random variables is equal to the sum of their variances.<br>
!>  A disadvantage of the variance for practical applications is that, unlike the standard deviation,
!>  its units differ from the random variable, which is why the standard deviation is more commonly
!>  reported as a measure of dispersion once the calculation is finished.<br>
!>
!>  Population vs. sample variance
!>  ------------------------------
!>
!>  There are two distinct concepts that are both called *variance*.<br>
!>  One, as discussed above, is part of a theoretical probability distribution and is **defined by an equation**.<br>
!>  The other variance is a **characteristic of a set of observations**.<br>
!>  When variance is calculated from observations, those observations are typically measured from a real world system.<br>
!>  If all possible observations of the system are present then the calculated variance is called the **population variance**.<br>
!>  Normally, however, only a subset is available, and the variance calculated from this is called the **sample variance**.<br>
!>  The variance calculated from a sample is considered an estimate of the full population variance.<br>
!>  There are multiple ways to calculate an **estimate** of the population variance, as discussed in the section below.<br>
!>  The two kinds of variance are closely related.<br>
!>  To see how, consider that a theoretical probability distribution can be used as a generator of hypothetical observations.<br>
!>  If an infinite number of observations are generated using a distribution, then the sample variance calculated
!>  from that infinite set will match the value calculated using the distribution's equation for variance.<br>
!>
!>  History
!>  -------
!>
!>  The term **variance** was apparently first introduced by Ronald Fisher in his 1918 paper
!>  *The Correlation Between Relatives on the Supposition of Mendelian Inheritance*.<br>
!>
!>  Definition
!>  ----------
!>
!>  The variance of a random variable \f$X\f$ is the expected value of the squared deviation from the mean of \f$X\f$, \f$\mu= \up{E}[X]\f$:<br>
!>  \f{equation}{
!>      \up{Var}(X) = \up{E}\left[(X - \mu)^2\right] ~.
!>  \f}
!>  This definition encompasses random variables that are generated by processes that are discrete, continuous, neither, or mixed.<br>
!>  The variance can also be thought of as the covariance of a random variable with itself:<br>
!>  \f{equation}{
!>      \up{Var}(X) = \up{Cov}(X, X) ~.
!>  \f}
!>  The variance is also equivalent to the second **cumulant** of the probability distribution that generates \f$X\f$.<br>
!>  The variance is typically designated as \f$\up{Var}(X)\f$, or sometimes as \f$V(X)\f$ or \f$\mathbb{V}(X)\f$, or symbolically as \f$\sigma_X^2\f$ or simply \f$\sigma^2\f$.<br>
!>  The expression for the variance can be expanded as follows:<br>
!>  \f{equation}{
!>      \begin{aligned}
!>          \up{Var}(X)
!>          &= \up{E} \left[ (X - \up{E}[X])^2\right] \\
!>          &= \up{E} \left[ X^2 - 2X\up{E}[X] + \up{E} [X]^{2} \right] \\
!>          &= \up{E} \left[ X^2 \right] - 2\up{E}[X]\up{E}[X] + \up{E}[X]^{2} \\
!>          &= \up{E} \left[ X^2 \right] - \up{E}[X]^{2}
!>      \end{aligned}
!>  \f}
!>  In other words, the variance of \f$X\f$ is equal to the mean of the square of \f$X\f$ minus the square of the mean of \f$X\f$.<br<
!>  However, this formulation is never used for computations using floating point arithmetic, because it suffers from
!>  [catastrophic cancellation](https://en.wikipedia.org/wiki/Catastrophic_cancellation) if the two components of the equation are similar in magnitude.<br>
!>
!>  Population variance and sample variance
!>  ---------------------------------------
!>
!>  Real-world observations such as the measurements of yesterday rain throughout the day typically cannot be complete sets of all possible observations that could be made.<br>
!>  As such, the variance calculated from the finite set will in general not match the variance that would have been calculated from the full population of possible observations.<br>
!>  This means that one estimates the mean and variance from a limited set of observations by using an estimator equation.<br>
!>  The estimator is a function of the sample of \f$n\f$ observations drawn without observational bias from the whole population of potential observations.<br>
!>  The simplest estimators for population mean and population variance are simply the mean and variance of the sample, the sample mean and (uncorrected) sample variance.<br>
!>  These are consistent estimators as they converge to the correct value as the number of samples increases), but can be improved.<br>
!>  Estimating the population variance by taking the sample variance is close to optimal in general, but can be improved in two ways.<br>
!>  The sample variance is computed as an average of squared deviations about the (sample) mean, by dividing by \f$n\f$.<br>
!>  However, using values other than \f$n\f$ improves the estimator in various ways.<br>
!>  Four common values for the denominator are \f$n\f$, \f$n − 1\f$, \f$n + 1\f$, and \f$n − 1.5\f$:<br>
!>  <ol>
!>      <li>    \f$n\f$ is the simplest (population variance of the sample),
!>      <li>    \f$n − 1\f$ eliminates bias,
!>      <li>    \f$n + 1\f$ minimizes mean squared error for the normal distribution,
!>      <li>    \f$n − 1.5\f$ mostly eliminates bias in unbiased estimation of standard deviation for the normal distribution.<br>
!>  </ol>
!>  Firstly, if the true population mean is unknown, then the sample variance (which uses the sample mean in place of the true mean) is a biased estimator:<br>
!>  It underestimates the variance by a factor of \f$(n − 1) / n\f$.<br>
!>  Correcting by this factor (dividing by \f$n − 1\f$ instead of \f$n\f$) is called the [Bessel correction](https://en.wikipedia.org/wiki/Bessel%27s_correction).<br>
!>  The resulting estimator is unbiased, and is called the **corrected sample variance** or **unbiased sample variance**.<br>
!>  For example, when \f$n = 1\f$ the variance of a single observation about the sample mean (itself) is obviously zero regardless of the population variance.<br>
!>  If the mean is determined in some other way than from the same samples used to estimate the variance then this bias does not arise
!>  and the variance can safely be estimated as that of the samples about the (independently known) mean.<br>
!>
!>  Biased sample variance
!>  ----------------------
!>
!>  In many practical situations, the true variance of a population is not known a priori and must be computed.<br>
!>  When dealing with extremely large populations, it is not possible to count every object in the population, so the computation must be performed on a sample of the population.<br>
!>  This is generally referred to as sample variance or empirical variance.<br>
!>  Sample variance can also be applied to the estimation of the variance of a continuous distribution from a sample of that distribution.<br>
!>  We take a sample with replacement of \f$n\f$ values \f$X_1, \ldots, X_n\f$ from the population and estimate the variance on the basis of this sample.<br>
!>  Directly taking the variance of the sample data gives the average of the squared deviations:<br>
!>  \f{equation}{
!>      \tilde{\sigma}_X^2 = \frac{1}{n} \sum_{i = 1}^{n} \left( X_i - \overline{X} \right)^2 = \left( \frac{1}{n} \sum_{i = 1}^{n} X_i^2 \right) - \overline{X}^2 = \frac{1}{n^2} \sum_{i, j ~:~ i < j} \left( X_i - X_j \right)^2 ~.<br>
!>  \f}
!>  Here, \f$\overline{X}\f$ denotes the sample mean:<br>
!>  \f{equation}{
!>      \overline{X} = \frac{1}{n} \sum_{i = 1}^n X_i ~.
!>  \f}
!>  Since the \f$X_i\f$ are selected randomly, both \f$\overline{X}\f$ and \f$\tilde{\sigma}_X^2\f$ are random variables.<br>
!>  Their expected values can be evaluated by averaging over the ensemble of all possible samples \f$X_i\f$ of size \f$n\f$ from the population.<br>
!>  For \f$\tilde{\sigma}_X^2\f$ this gives:<br>
!>  \f{equation}{
!>      \begin{aligned}
!>          \up{E}[\tilde{\sigma}_X^2]
!>          &= \up{E} \left[ \frac{1}{n} \sum_{i = 1}^n \left( X_i - \frac{1}{n} \sum_{j = 1}^n X_{j} \right)^2 \right] \\
!>          &= \frac{1}{n} \sum_{i = 1}^n \up{E} \left[ X_i^2 - \frac{2}{n} X_i \sum_{j = 1}^n X_j + \frac{1}{n^2} \sum_{j = 1}^n X_j \sum_{k = 1}^n X_k \right] \\
!>          &= \frac{1}{n} \sum_{i = 1}^n \left( \frac{n - 2}{n} \up{E} \left[ X_i^2 \right] - \frac{2}{n} \sum_{j \neq i} \up{E} \left[ X_i X_j \right] + \frac{1}{n^2} \sum_{j = 1}^n \sum_{k \neq j}^n \up{E} \left[ X_j X_k \right] + \frac{1}{n^2} \sum_{j = 1}^n \up{E}\left[ X_j^2 \right] \right) \\
!>          &= \frac{1}{n} \sum_{i = 1}^n \left[ \frac{n - 2}{n} \left( \sigma^2 + \mu^2 \right) - \frac{2}{n}(n - 1)\mu^2 + \frac{1}{n^2} n (n - 1) \mu^2 + \frac{1}{n} \left(\sigma^2 + \mu^2 \right)\right] \\
!>          &= \frac{n - 1}{n} \sigma^2 ~.
!>      \end{aligned}
!>  \f}
!>  Hence \f$\tilde{\sigma}_X^2\f$ gives an estimate of the population variance that is biased by a factor of \f$\frac{n - 1}{n}\f$.<br>
!>  For this reason, \f$\tilde{\sigma}_X^2\f$ is referred to as the **biased sample variance**.<br>
!>  The bias-correction factor in this case is \f$\xi = \frac{n}{n - 1}\f$.<br>
!>
!>  Unbiased sample variance
!>  ------------------------
!>
!>  Correcting for this bias yields the unbiased sample variance, denoted \f$\sigma^2\f$:<br>
!>  \f{equation}{
!>      \sigma^{2} = \frac{n}{n-1} \tilde{\sigma}_X^2 = \frac{n}{n - 1} \left[ \frac{1}{n} \sum_{i = 1}^n \left(X_i - \overline{X} \right)^2\right] = \frac{1}{n-1}\sum_{i = 1}^n\left(X_i - \overline{X} \right)^2 ~.
!>  \f}
!>  Either estimator may be simply referred to as the sample variance when the version can be determined by context.<br>
!>  The same proof is also applicable for samples taken from a continuous probability distribution.<br>
!>  The use of the term \f$n − 1\f$ is called **the Bessel correction**, and it is also used in sample covariance and the sample standard deviation (the square root of variance).<br>
!>  The square root is a concave function and thus introduces negative bias (by the Jensen inequality), which depends on the distribution, and thus the corrected sample standard deviation is biased.<br>
!>  The unbiased estimation of standard deviation is a technically involved problem, though for the normal distribution using the term \f$n − 1.5\f$ yields an almost unbiased estimator.<br>
!>
!>  Biased weighted sample variance
!>  -------------------------------
!>
!>  Given a **weighted sample** of \f$n\f$ observations \f$X_{1:n}\f$ with a [weighted sample mean](@ref pm_sampleMean) \f$\hat\mu_w\f$,
!>  the biased variance of the weighted sample is computed as,
!>  \f{equation}{
!>      \hat\sigma_w^2 = \frac{1}{\sum_{i=1}^{n} w_i} \sum_{i = 1}^n w_i (X_i - \hat\mu_w)^2 ~,
!>  \f}
!>  where `n = nsam` is the number of observations in the sample, \f$w_i\f$ are the weights of individual data points, and \f$\hat\mu_w\f$ is the **weighted mean** of the sample.<br>
!>  When the sample size is small, the above equation yields a biased estimate of the variance.<br>
!>
!>  Unbiased weighted sample variance
!>  ---------------------------------
!>
!>  There is no unique generic equation for the **unbiased variance** of a **weighted sample**.<br>
!>  However, depending on the types of the weights involved, a few popular definitions exist.<br>
!>
!>  <ol>
!>      <li>    The **unbiased variance** of a sample with **frequency**, **count**, or **repeat** weights can be computed via the following equation,
!>              \f{equation}{
!>                  \hat\sigma_w^2
!>                  = \frac{\xi}{\sum_{i=1}^{n} w_i} \sum_{i=1}^{n} w_i (X_i - \hat\mu_w)^2
!>                  = \frac{1}{\left( \sum_{i=1}^{n} w_i \right) - 1} \sum_{i=1}^{n} w_i (X_i - \hat\mu_w)^2 ~,
!>              \f}
!>              where the bias correction factor \f$\xi\f$ is,
!>              \f{equation}{
!>                  \xi = \frac{\sum_{i=1}^{n} w_i}{\left( \sum_{i=1}^{n} w_i \right) - 1} ~.
!>              \f}
!>              [Frequency weights](@ref pm_sampleWeight::fweight_type) represent the number of duplications of each observation in the sample whose population variance is to be estimated.<br>
!>              Therefore, the **frequency weights** are expected to be **integers** or **whole numbers**.<br>
!>      <li>    The **unbiased variance** of a sample with **reliability weights**, also sometimes confusingly known as **probability weights** or **importance weights**,
!>              can be computed by the following equation,
!>              \f{equation}{
!>                  \hat\sigma_w^2
!>                  = \frac{\xi}{\sum_{i=1}^{n}w_i} \sum_{i=1}^{n} w_i \left(X_i - \hat\mu_w\right)^2
!>                  = \frac{\sum_{i=1}^{n} w_i}{\left(\sum_{i=1}^{n}w_i\right)^2 - \left(\sum_{i=1}^{n}w_i^2\right)} \sum_{i=1}^{n} w_i \left(X_i - \hat\mu_w\right)^2
!>                  ~,
!>              \f}
!>              where the bias correction factor \f$\xi\f$ is,
!>              \f{equation}{
!>                  \xi = \frac{\left(\sum_{i=1}^{n} w_i\right)^2}{\left(\sum_{i=1}^{n}w_i\right)^2 - \left(\sum_{i=1}^{n}w_i^2\right)} ~.
!>              \f}
!>              <ol>
!>                  <li>    [Reliability weights](@ref pm_sampleWeight::rweight_type) weights, also known as **reliability weights** or **sampling weights**
!>                          represent the probability of a case (or subject) being selected into the sample from a population.<br>
!>                  <li>    Application of the term *unbiased* to the above equation is controversial as some believe that bias
!>                          cannot be correct without the knowledge of the sample size, which is lost in normalized weights.<br>
!>                  <li>    Reliability weights are frequently (but not necessarily) normalized, meaning that \f$\sum^{i = 1}_{n} w_i = 1\f$.<br>
!>              </ol>
!>  </ol>
!>
!>  Variance updating
!>  -----------------
!>
!>  Consider a sample \f$X_{1..n}\f$ of size \f$n\f$ with weights \f$w_{1..n}\f$.<br>
!>  The weighted mean of this sample can be expressed as,
!>  \f{equation}{
!>      \large
!>      \hat\mu_{1:n} = \frac{1}{w_{1:n}} \sum_{i = 1}^{n} w_i X_i ~,
!>  \f}
!>  where \f$w_{1:n} = \sum w_{1..n}\f$ and the weighted **biased** variance of the sample can be expressed as,
!>  \f{eqnarray}{
!>      \large
!>      \tilde{\sigma}_{1:n}^2 = \frac{\hat\sigma_{1:n}^2}{\xi_{1:n}}
!>      &=& \frac{1}{w_{1:n}} \sum_{i = 1}^{n} w_i (X_i - \hat\mu_{1:n})^2 ~,
!>      &=& \frac{1}{w_{1:n}} \sum_{i = 1}^{n} w_i X_i^2 - \xi_{1:n} \hat\mu_{1:n}^2 ~,
!>  \f}
!>  where \f$\xi_{1:n}\f$ is the appropriate bias-correction factor (which can be one).<br>
!>  This implies,
!>  \f{equation}{
!>      \large
!>      \sum_{i = 1}^{n} w_i X_i^2 = w_{1:n} \left(\frac{\hat\sigma_{1:n}^2}{\xi_{1:n}} + \hat\mu_{1:n}\right) ~.
!>  \f}
!>  The above can be used to express the variance of the sample \f$X_{1..n+m}\f$ resulting from the
!>  merger of two separate samples \f$X_{1..n}\f$ and \f$X_{n+1..n+m}\f$ as,
!>  \f{equation}{
!>      \large
!>      \frac{\hat\sigma_{1:n+m}^2}{\xi_{1:n+m}} = \frac{1}{w_{1:n+m}}
!>      \left(
!>          w_{1:n} \left( \frac{\hat\sigma_{1:n}^2}{\xi_{1:n}} + \hat\mu_{1:n}^2 \right) +
!>          w_{n+1:n+m} \left(\frac{\hat\sigma_{n+1:n+m}^2}{\xi_{n+1:n+m}} + \hat\mu_{n+1:n+m}^2 \right) -
!>          w_{1:n+m} \hat\mu_{1:n+m}^2
!>      \right) ~.
!>  \f}
!>
!>  \note
!>  Note the effects of bias-correction in computing the variance become
!>  noticeable only for sample sample sizes (i.e., when `nsam` is small).<br>
!>
!>  \note
!>  For a two or higher-dimensional `sample`, if the variance is to be computed for the entire `sample` (as opposed to computing it along
!>  a particular dimension), simply pass `reshape(sample, shape = size(sample))` to the appropriate [getVar](@ref pm_sampleVar::getVar) interface.<br>
!>  Alternatively, a 1D pointer of the same size as the multidimensional sample can be passed to the procedure.<br>
!>
!>  \note
!>  While it is tempting to extend this generic interface to `weight` arguments of type `integer` or `real` of various kinds,
!>  such extensions do not appear to add any tangible benefits beyond making the interface more flexible for the end user.<br>
!>  But such extensions would certainly make the maintenance and future extensions of this interface difficult and complex.<br>
!>  In the case of `integer` (frequency) weights,
!>  <ol>
!>      <li>    The summation `sum(weight)` involved in the computation may lead to an integer-overflow if individual weights are too large.<br>
!>      <li>    Avoiding overflow would then require coercing the weights to `real` before summation, which add an extra layer of unnecessary type coercion.<br>
!>      <li>    Furthermore, according to the coercion rules of the Fortran standard, if an `integer` is multiplied with a `real`,
!>              the `integer` value must be first converted to `real` of the same kind as the real value, then multiplied.<br>
!>      <li>    The type coercion to `real` will have to happen a second time when the weights are multiplied with the data values.<br>
!>      <li>    Each integer-real type coercion costs about a `real` multiplication on modern hardware (See, e.g., [this thread](https://stackoverflow.com/a/28668349/2088694)).<br>
!>  </ol>
!>  By contrast,
!>  <ol>
!>      <li>    Real-valued weights, even if the weights are counts, do not require type coercion if `real` values in the computation are of the same kind as is here.<br>
!>      <li>    The floating-point multiplication tends to be faster than integer multiplication on most modern architecture.<br>
!>      <li>    However, real-valued weight summation is 4-8 times more expensive then `integer` addition, but less than `real` multiplication.<br>
!>  </ol>
!>  Considering all factors in the above, there does not seem to exist any performance benefits with providing dedicated interfaces for `weight` arguments of different type and kind.<br>
!>  The following list compares the cost and latencies of some of the basic operations involving `integer` and `real` numbers.<br>
!>  <ol>
!>      <li>    Central Processing Unit (CPU):
!>              <ol>
!>                  <li>    Integer add: 1 cycle
!>                  <li>    32-bit integer multiply: 10 cycles
!>                  <li>    64-bit integer multiply: 20 cycles
!>                  <li>    32-bit integer divide: 69 cycles
!>                  <li>    64-bit integer divide: 133 cycles
!>              </ol>
!>      <li>    On-chip Floating Point Unit (FPU):
!>              <ol>
!>                  <li>    Floating point add: 4 cycles
!>                  <li>    Floating point multiply: 7 cycles
!>                  <li>    Double precision multiply: 8 cycles
!>                  <li>    Floating point divide: 23 cycles
!>                  <li>    Double precision divide: 36 cycles
!>              </ol>
!>  </ol>
!>
!>  Generalization of variance
!>  ==========================
!>
!>  Variance of complex variables
!>  -----------------------------
!>
!>  If \f$x\f$ is a scalar complex-valued random variable, with values in \f$\mathbb{C}\f$, then its variance is \f$\up{E} \left[(x-\mu )(x-\mu )^{*}\right]\f$, where \f$x^{*}\f$ is the complex conjugate of \f$x\f$.<br>
!>  This variance is a real scalar.<br>
!>
!>  Matrix variance of vector-valued variables
!>  ------------------------------------------
!>
!>  If \f$X\f$ is a vector-valued random variable, with values in \f$\mathbb{R}^{n}\f$, and thought of as a column vector, then a natural generalization of variance is \f$\up{E}\left[(X-\mu)(X-\mu)^{\up{T}}\right]\f$,
!>  where \f$\mu = \up{E}(X)\f$ and \f$X^{\up{T}}\f$ is the transpose of \f$X\f$, and so is a row vector.<br>
!>  The result is a positive semi-definite square matrix, commonly referred to as the **variance-covariance matrix** (or simply as the [covariance matrix](@ref pm_sampleCov)).<br>
!>
!>  Matrix variance of complex vector-valued variables
!>  --------------------------------------------------
!>
!>  If \f$X\f$ is a vector- and complex-valued random variable, with values in \f$\mathbb {C}^{n}\f$, then the covariance matrix is \f$\up{E} \left[(X-\mu)(X-\mu)^{\dagger}\right]\f$,
!>  where \f$X^{\dagger}\f$ is the conjugate transpose of \f$X\f$.<br>
!>  This matrix is also positive semi-definite and square.<br>
!>
!>  Extension of scalar variance to higher dimensions
!>  -------------------------------------------------
!>
!>  Another generalization of variance for vector-valued random variables \f$X\f$, which results in a scalar value rather than in a matrix, is the generalized variance \f$\det(C)\f$,
!>  the determinant of the covariance matrix.<br>
!>  The generalized variance can be shown to be related to the multidimensional scatter of points around their mean.<br>
!>
!>  A different generalization is obtained by considering the variance of the Euclidean distance between the random variable and its mean.<br>
!>  This results in \f$\up{E} \left[(X-\mu)^{\up {T}}(X-\mu)\right] = \up{tr}(C)\f$, which is the **trace of the covariance matrix**.<br>
!>
!>  Summary
!>  -------
!>
!>  <ol>
!>      <li>    Variance is a measure of dispersion in data.<br>
!>      <li>    If the population mean is known a priori independent of the current sample based upon which the variance is to be computed,
!>              then the computed sample variance is an **unbiased estimate** of the true population variance.<br>
!>      <li>    If the population mean is unknown a priori and must be computed via the same sample based upon which the variance is to be computed,
!>              then the computed sample variance is an **biased estimate** of the true population variance.<br>
!>      <li>    This variance bias can be corrected by applying appropriate correction factors to the computed variances.<br>
!>      <li>    The most famous such correction is called the [Bessel correction](https://en.wikipedia.org/wiki/Bessel%27s_correction).<br>
!>  </ol>
!>
!>  \see
!>  [pm_sampling](@ref pm_sampling)<br>
!>  [pm_sampleACT](@ref pm_sampleACT)<br>
!>  [pm_sampleCCF](@ref pm_sampleCCF)<br>
!>  [pm_sampleCor](@ref pm_sampleCor)<br>
!>  [pm_sampleCov](@ref pm_sampleCov)<br>
!>  [pm_sampleConv](@ref pm_sampleConv)<br>
!>  [pm_sampleECDF](@ref pm_sampleECDF)<br>
!>  [pm_sampleMean](@ref pm_sampleMean)<br>
!>  [pm_sampleNorm](@ref pm_sampleNorm)<br>
!>  [pm_sampleQuan](@ref pm_sampleQuan)<br>
!>  [pm_sampleScale](@ref pm_sampleScale)<br>
!>  [pm_sampleShift](@ref pm_sampleShift)<br>
!>  [pm_sampleWeight](@ref pm_sampleWeight)<br>
!>  [pm_sampleAffinity](@ref pm_sampleAffinity)<br>
!>  [pm_sampleVar](@ref pm_sampleVar)<br>
!>  Box and Tiao, 1973, *Bayesian Inference in Statistical Analysis*, Page 421.<br>
!>  *Updating mean and variance estimates: an improved method*, D.H.D. West, 1979.<br>
!>  Geisser and Cornfield, 1963, *Posterior distributions for multivariate normal parameters*.<br>
!>
!>  \test
!>  [test_pm_sampleCov](@ref test_pm_sampleCov)
!>
!>  \bug
!>  \status See \unresolved, See [this page](https://groups.google.com/g/comp.lang.fortran/c/NDE6JKTFbNU) for more information.<br>
!>  \source \gfortran
!>  \desc
!>  Ideally, there should be only one generic interface in this module for computing the biased/corrected/weighted variance.<br>
!>  This requires ability to resolve the different weight types, which requires custom derived types for weights.<br>
!>  Fortran PDTs are ideal for such use cases. However, the implementation of PDTs is far from complete in \gfortran.<br>
!>  \remedy
!>  Given that the importance of \gfortran support, separate generic interfaces were instead developed for different sample weight types.<br>
!>  Once the \gfortran PDT bugs are resolved, the [getVar](@ref pm_sampleVar::getVar) generic interface can be extended to serve as a
!>  high-level wrapper for the weight-specific generic interfaces in this module.<br>
!>
!>  \final
!>
!>  \todo
!>  \pmed
!>  The inclusion of bias correction in the calculation of variance is a
!>  frequentist abomination and shenanigan that must be eliminated in the future.<br>
!>  The correction factor should be computed separately from the actual variance calculation.<br>
!>
!>  \author
!>  \AmirShahmoradi, Nov 24, 2020, 4:19 AM, Dallas, TX<br>
!>  \FatemehBagheri, Thursday 12:45 AM, August 20, 2021, Dallas, TX<br>
!>  \AmirShahmoradi, Monday March 6, 2017, 2:48 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_sampleVar

    use pm_kind, only: SK, IK, LK
    use pm_sampleWeight, only: weight_type
    use pm_sampleWeight, only: fweight_type, fweight
    use pm_sampleWeight, only: rweight_type, rweight

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_sampleVar"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the bias correction factor for the computation of the variance of a (weighted) sample.<br>
    !>
    !>  \details
    !>  When the population mean is replaced with sample mean, multiplying the output of this generic interface with the output (co)variance
    !>  from [setVar](@ref pm_sampleVar::setVar) or [setCov](@ref pm_sampleCov::setCov) yields an unbiased estimate of the sample variance.<br>
    !>  See the documentation of the parent module [pm_sampleVar](@ref pm_sampleVar) for algorithmic details and variance definitions.<br>
    !>
    !>  \param[in]  weisum      :   The input scalar of,<br>
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing either,<br>
    !>                              <ol>
    !>                                  <li>    the size of an unweighted sample, **only if** the input argument `correction` is missing or is not [rweight_type](@ref pm_sampleWeight::rweight_type),
    !>                                  <li>    the quantity `sum(weight)` of a weighted sample, **only if** the input argument `correction` is missing or is not [rweight_type](@ref pm_sampleWeight::rweight_type),
    !>
    !>                              </ol>
    !>                              based upon which the bias correction is computed.<br>
    !>                              The vector `weight(1:nsam)` refers to the weights of a sample of size `nsam`.<br>
    !>  \param[in]  weisqs      :   The input scalar of the same type and kind as `weisum`, containing the sum of squared sample weight: `sum(weight**2)`.<br>
    !>                              (**optional**. If missing, the output corresponds to the Bessel correction (i.e., [frequency weight](@ref pm_sampleWeight::fweight_type) is assumed).<br>
    !>                              If present, the correction corresponding to [reliability weights](@ref pm_sampleWeight::rweight_type) is returned.)
    !>
    !>  \return
    !>  `normfac`               :   The output scalar of the same type and kind as `weisum` containing the bias correction factor of variance equation.<br>
    !>                              Multiplying `normfac` with the sum of the sample weight (or its size, if unweighted) yields in the Bias Correction factor.<br>
    !>
    !>  \interface{getVarCorrection}
    !>  \code{.F90}
    !>
    !>      use pm_sampleVar, only: getVarCorrection
    !>
    !>      normfac = getVarCorrection(weisum) ! unweighted or frequency-weighted sample variance correction
    !>      normfac = getVarCorrection(weisum, weisqs) ! reliability-weighted sample variance correction
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < weisum` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < weisqs` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  Note the effects of bias-correction in computing the variance become noticeable only for very sample sample sizes (i.e., when `nsam` is small).<br>
    !>
    !>  \see
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>
    !>  \example{getVarCorrection}
    !>  \include{lineno} example/pm_sampleVar/getVarCorrection/main.F90
    !>  \compilef{getVarCorrection}
    !>  \output{getVarCorrection}
    !>  \include{lineno} example/pm_sampleVar/getVarCorrection/main.out.F90
    !>  \postproc{getVarCorrection}
    !>  \include{lineno} example/pm_sampleVar/getVarCorrection/main.py
    !>  \vis{getVarCorrection}
    !>  \image html pm_sampleVar/getVarCorrection/getVarCorrection.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_sampleVar](@ref test_pm_sampleVar)
    !>
    !>  \final{getVarCorrection}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX<br>

    interface getVarCorrection

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getVarCorrectionFreq_RK5(weisum) result(correction)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarCorrectionFreq_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)                :: weisum
        real(TKG)                                           :: correction
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getVarCorrectionFreq_RK4(weisum) result(correction)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarCorrectionFreq_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)                :: weisum
        real(TKG)                                           :: correction
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getVarCorrectionFreq_RK3(weisum) result(correction)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarCorrectionFreq_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)                :: weisum
        real(TKG)                                           :: correction
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getVarCorrectionFreq_RK2(weisum) result(correction)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarCorrectionFreq_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)                :: weisum
        real(TKG)                                           :: correction
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getVarCorrectionFreq_RK1(weisum) result(correction)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarCorrectionFreq_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)                :: weisum
        real(TKG)                                           :: correction
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getVarCorrectionReli_RK5(weisum, weisqs) result(correction)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarCorrectionReli_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)                :: weisum, weisqs
        real(TKG)                                           :: correction
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getVarCorrectionReli_RK4(weisum, weisqs) result(correction)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarCorrectionReli_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)                :: weisum, weisqs
        real(TKG)                                           :: correction
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getVarCorrectionReli_RK3(weisum, weisqs) result(correction)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarCorrectionReli_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)                :: weisum, weisqs
        real(TKG)                                           :: correction
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getVarCorrectionReli_RK2(weisum, weisqs) result(correction)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarCorrectionReli_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)                :: weisum, weisqs
        real(TKG)                                           :: correction
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getVarCorrectionReli_RK1(weisum, weisqs) result(correction)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarCorrectionReli_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)                :: weisum, weisqs
        real(TKG)                                           :: correction
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the variance of the input sample of type `complex` or `real` of shape `(nsam)` or `(ndim, nsam)` or `(nsam, ndim)`
    !>  where `ndim` is the number of data dimensions (the number of data attributes) and `nsam` is the number of data points.<br>
    !>
    !>  \details
    !>  See the documentation of the parent module [pm_sampleVar](@ref pm_sampleVar) for algorithmic details and variance definitions.<br>
    !>
    !>  \param[in]  sample      :   The input `contiguous` array of shape `(nsam)`, `(ndim, nsam)`, or `(nsam, ndim)` of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the sample comprised of `nsam` observations each with `ndim` attributes.<br>
    !>                              If `sample` is a matrix, then the direction along which variance is computed is dictated by the input argument `dim`.<br>
    !>  \param[in]  dim         :   The input scalar `integer` of default kind \IK representing the dimension (`1` or `2`) of the input `sample` along which the mean must be computed.<br>
    !>                              <ol>
    !>                                  <li>    If `dim = 1`, the input `sample` of rank `2` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                                  <li>    If `dim = 2`, the input `sample` of rank `2` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                              </ol>
    !>                              The input `dim` must always be `1` or missing for an input `sample` of rank `1`.<br>
    !>                              (**optional**. If missing, the variance of the whole input `sample` is computed.)
    !>  \param[in]  weight      :   The `contiguous` vector of length `nsam` of,
    !>                              <ol>
    !>                                  <li>    type `integer` of default kind \IK, or
    !>                                  <li>    type `real` of the same kind as the input `sample`,
    !>                              </ol>
    !>                              containing the corresponding weights of individual `nsam` observations in `sample`.<br>
    !>                              (**optional**. default = `getFilled(1, size(sample, dim))` if `dim` is present, or `getFilled(1, size(sample))` if `dim` is missing.)
    !>  \param[in]  correction  :   The input scalar object that can be the following:<br>
    !>                              <ol>
    !>                                  <li>    The constant [fweight](@ref pm_sampleWeight::fweight) or an object of type [fweight](@ref pm_sampleWeight::fweight)
    !>                                          implying a bias correction based on the assumption of **frequency weights** for the sample observations, even if the `weight` argument is missing.<br>
    !>                                          This is the most popular type of correction, also known as the [Bessel correction](https://en.wikipedia.org/wiki/Bessel%27s_correction).<br>
    !>                                  <li>    The constant [rweight](@ref pm_sampleWeight::rweight) or an object of type [rweight_type](@ref pm_sampleWeight::rweight_type)
    !>                                          implying a bias correction based on the assumption of **reliability weights** for the sample observations.<br>
    !>                              </ol>
    !>                              (**optional**. If missing, no bias-correction will be applied to the output `var`.)
    !>
    !>  \return
    !>  `var`                   :   The output object of type `real` of the same kind as the input `sample` of rank `rank(sample) - 1`, representing its variance:
    !>                              <ol>
    !>                                  <li>    If `sample` is a vector, the output `var` is a scalar.<br>
    !>                                  <li>    If `sample` is a matrix, the output `var` is an `allocatable` vector of size `ndim = size(sample, dim)`.<br>
    !>                              </ol>
    !>
    !>  \interface{getVar}
    !>  \code{.F90}
    !>
    !>      use pm_sampleVar, only: getVar
    !>
    !>      var = getVar(sample(:), weight = weight(1 : size(sample)), correction = correction)
    !>      var = getVar(sample(:,:), weight = weight(1 : size(sample)), correction = correction)
    !>
    !>      var = getVar(sample(:), weight = weight(1 : size(sample, dim)), correction = correction)
    !>      var(:) = getVar(sample(:,:), dim, weight = weight(1 : size(sample, dim)), correction = correction)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  This generic interface is merely a flexible wrapper around the generic `subroutine` interface [setVar](@ref pm_sampleVar::setVar).<br>
    !>  As such, all conditions and warnings associated with [setVar](@ref pm_sampleVar::setVar) equally hold for this generic interface.<br>
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  Note the effects of bias-correction in computing the variance become noticeable only for very sample sample sizes (i.e., when `nsam` is small).<br>
    !>
    !>  \note
    !>  For a one-dimensional `sample`, one can also use the concise Fortran syntax to achieve the same goal as with
    !>  the interface `var = getVar(sample(:), mean, correction)` with integer `weight`:<br>
    !>  \code{.F90}
    !>
    !>       mean = sum(sample) / size(sample)
    !>       var = sum((sample-mean)**2) / (size(sample) - 1)
    !>
    !>  \endcode
    !>  But the above concise version will be likely slightly slower as it potentially involves more loops.<br>
    !>
    !>  \note
    !>  For a two or higher-dimensional `sample`, if the variance is to be computed for the entire `sample` (as opposed to computing it along
    !>  a particular dimension), simply pass `reshape(sample, shape = size(sample))` to the appropriate [getVar](@ref pm_sampleVar::getVar) interface.<br>
    !>  Alternatively, a 1D pointer of the same size as the multidimensional sample can be passed to the procedure.<br>
    !>
    !>  \see
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getShifted](@ref pm_sampleShift::getShifted)<br>
    !>  [setShifted](@ref pm_sampleShift::setShifted)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>
    !>  \example{getVar}
    !>  \include{lineno} example/pm_sampleVar/getVar/main.F90
    !>  \compilef{getVar}
    !>  \output{getVar}
    !>  \include{lineno} example/pm_sampleVar/getVar/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleVar](@ref test_pm_sampleVar)
    !>
    !>  \todo
    !>  \pvlow
    !>  The functionality of this interface can be expanded in the future to include the computation of
    !>  the variance of higher dimensional input `sample` and whole `sample` input arrays of arbitrary shape.<br>
    !>
    !>  \final{getVar}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX<br>

    ! ALL D1

    interface getVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getVarALL_WNO_D1_CK5(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK4_ENABLED
    PURE module function getVarALL_WNO_D1_CK4(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK3_ENABLED
    PURE module function getVarALL_WNO_D1_CK3(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK2_ENABLED
    PURE module function getVarALL_WNO_D1_CK2(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK1_ENABLED
    PURE module function getVarALL_WNO_D1_CK1(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getVarALL_WNO_D1_RK5(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK4_ENABLED
    PURE module function getVarALL_WNO_D1_RK4(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK3_ENABLED
    PURE module function getVarALL_WNO_D1_RK3(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK2_ENABLED
    PURE module function getVarALL_WNO_D1_RK2(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK1_ENABLED
    PURE module function getVarALL_WNO_D1_RK1(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getVarALL_WTI_D1_CK5(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK4_ENABLED
    PURE module function getVarALL_WTI_D1_CK4(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK3_ENABLED
    PURE module function getVarALL_WTI_D1_CK3(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK2_ENABLED
    PURE module function getVarALL_WTI_D1_CK2(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK1_ENABLED
    PURE module function getVarALL_WTI_D1_CK1(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getVarALL_WTI_D1_RK5(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK4_ENABLED
    PURE module function getVarALL_WTI_D1_RK4(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK3_ENABLED
    PURE module function getVarALL_WTI_D1_RK3(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK2_ENABLED
    PURE module function getVarALL_WTI_D1_RK2(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK1_ENABLED
    PURE module function getVarALL_WTI_D1_RK1(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getVarALL_WTR_D1_CK5(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK4_ENABLED
    PURE module function getVarALL_WTR_D1_CK4(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK3_ENABLED
    PURE module function getVarALL_WTR_D1_CK3(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK2_ENABLED
    PURE module function getVarALL_WTR_D1_CK2(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK1_ENABLED
    PURE module function getVarALL_WTR_D1_CK1(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getVarALL_WTR_D1_RK5(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK4_ENABLED
    PURE module function getVarALL_WTR_D1_RK4(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK3_ENABLED
    PURE module function getVarALL_WTR_D1_RK3(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK2_ENABLED
    PURE module function getVarALL_WTR_D1_RK2(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK1_ENABLED
    PURE module function getVarALL_WTR_D1_RK1(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getVar

    ! ALL D2

    interface getVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getVarALL_WNO_D2_CK5(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK4_ENABLED
    PURE module function getVarALL_WNO_D2_CK4(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK3_ENABLED
    PURE module function getVarALL_WNO_D2_CK3(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK2_ENABLED
    PURE module function getVarALL_WNO_D2_CK2(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK1_ENABLED
    PURE module function getVarALL_WNO_D2_CK1(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getVarALL_WNO_D2_RK5(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK4_ENABLED
    PURE module function getVarALL_WNO_D2_RK4(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK3_ENABLED
    PURE module function getVarALL_WNO_D2_RK3(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK2_ENABLED
    PURE module function getVarALL_WNO_D2_RK2(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK1_ENABLED
    PURE module function getVarALL_WNO_D2_RK1(sample, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WNO_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getVarALL_WTI_D2_CK5(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK4_ENABLED
    PURE module function getVarALL_WTI_D2_CK4(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK3_ENABLED
    PURE module function getVarALL_WTI_D2_CK3(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK2_ENABLED
    PURE module function getVarALL_WTI_D2_CK2(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK1_ENABLED
    PURE module function getVarALL_WTI_D2_CK1(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getVarALL_WTI_D2_RK5(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK4_ENABLED
    PURE module function getVarALL_WTI_D2_RK4(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK3_ENABLED
    PURE module function getVarALL_WTI_D2_RK3(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK2_ENABLED
    PURE module function getVarALL_WTI_D2_RK2(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK1_ENABLED
    PURE module function getVarALL_WTI_D2_RK1(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTI_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getVarALL_WTR_D2_CK5(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK4_ENABLED
    PURE module function getVarALL_WTR_D2_CK4(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK3_ENABLED
    PURE module function getVarALL_WTR_D2_CK3(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK2_ENABLED
    PURE module function getVarALL_WTR_D2_CK2(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK1_ENABLED
    PURE module function getVarALL_WTR_D2_CK1(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getVarALL_WTR_D2_RK5(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK4_ENABLED
    PURE module function getVarALL_WTR_D2_RK4(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK3_ENABLED
    PURE module function getVarALL_WTR_D2_RK3(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK2_ENABLED
    PURE module function getVarALL_WTR_D2_RK2(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK1_ENABLED
    PURE module function getVarALL_WTR_D2_RK1(sample, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarALL_WTR_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getVar

    ! DIM D1

    interface getVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getVarDIM_WNO_D1_CK5(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK4_ENABLED
    PURE module function getVarDIM_WNO_D1_CK4(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK3_ENABLED
    PURE module function getVarDIM_WNO_D1_CK3(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK2_ENABLED
    PURE module function getVarDIM_WNO_D1_CK2(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK1_ENABLED
    PURE module function getVarDIM_WNO_D1_CK1(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getVarDIM_WNO_D1_RK5(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK4_ENABLED
    PURE module function getVarDIM_WNO_D1_RK4(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK3_ENABLED
    PURE module function getVarDIM_WNO_D1_RK3(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK2_ENABLED
    PURE module function getVarDIM_WNO_D1_RK2(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK1_ENABLED
    PURE module function getVarDIM_WNO_D1_RK1(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getVarDIM_WTI_D1_CK5(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK4_ENABLED
    PURE module function getVarDIM_WTI_D1_CK4(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK3_ENABLED
    PURE module function getVarDIM_WTI_D1_CK3(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK2_ENABLED
    PURE module function getVarDIM_WTI_D1_CK2(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK1_ENABLED
    PURE module function getVarDIM_WTI_D1_CK1(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getVarDIM_WTI_D1_RK5(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK4_ENABLED
    PURE module function getVarDIM_WTI_D1_RK4(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK3_ENABLED
    PURE module function getVarDIM_WTI_D1_RK3(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK2_ENABLED
    PURE module function getVarDIM_WTI_D1_RK2(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK1_ENABLED
    PURE module function getVarDIM_WTI_D1_RK1(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getVarDIM_WTR_D1_CK5(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK4_ENABLED
    PURE module function getVarDIM_WTR_D1_CK4(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK3_ENABLED
    PURE module function getVarDIM_WTR_D1_CK3(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK2_ENABLED
    PURE module function getVarDIM_WTR_D1_CK2(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if CK1_ENABLED
    PURE module function getVarDIM_WTR_D1_CK1(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getVarDIM_WTR_D1_RK5(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK4_ENABLED
    PURE module function getVarDIM_WTR_D1_RK4(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK3_ENABLED
    PURE module function getVarDIM_WTR_D1_RK3(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK2_ENABLED
    PURE module function getVarDIM_WTR_D1_RK2(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

#if RK1_ENABLED
    PURE module function getVarDIM_WTR_D1_RK1(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:)
        real(TKG)                                                               :: var
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getVar

    ! DIM D2

    interface getVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getVarDIM_WNO_D2_CK5(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getVarDIM_WNO_D2_CK4(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getVarDIM_WNO_D2_CK3(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getVarDIM_WNO_D2_CK2(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getVarDIM_WNO_D2_CK1(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getVarDIM_WNO_D2_RK5(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getVarDIM_WNO_D2_RK4(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getVarDIM_WNO_D2_RK3(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getVarDIM_WNO_D2_RK2(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getVarDIM_WNO_D2_RK1(sample, dim, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WNO_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getVarDIM_WTI_D2_CK5(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getVarDIM_WTI_D2_CK4(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getVarDIM_WTI_D2_CK3(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getVarDIM_WTI_D2_CK2(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getVarDIM_WTI_D2_CK1(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getVarDIM_WTI_D2_RK5(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getVarDIM_WTI_D2_RK4(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getVarDIM_WTI_D2_RK3(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getVarDIM_WTI_D2_RK2(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getVarDIM_WTI_D2_RK1(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTI_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        integer(IK)             , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getVarDIM_WTR_D2_CK5(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getVarDIM_WTR_D2_CK4(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getVarDIM_WTR_D2_CK3(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getVarDIM_WTR_D2_CK2(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getVarDIM_WTR_D2_CK1(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getVarDIM_WTR_D2_RK5(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getVarDIM_WTR_D2_RK4(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getVarDIM_WTR_D2_RK3(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getVarDIM_WTR_D2_RK2(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getVarDIM_WTR_D2_RK1(sample, dim, weight, correction) result(var)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarDIM_WTR_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                    :: dim
        class(weight_type)      , intent(in)                    , optional      :: correction
        real(TKG)               , intent(in)    , contiguous                    :: weight(:)
        real(TKG)               , intent(in)    , contiguous                    :: sample(:,:)
        real(TKG)                                                               :: var(size(sample, 3 - dim, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the variance of an input (weighted) sample of type `complex` or `real` of shape `(nsam)` or `(ndim, nsam)` or `(nsam, ndim)`
    !>  where `ndim` is the number of data dimensions (the number of data attributes) and `nsam` is the number of data points.<br>
    !>
    !>  \details
    !>  See the documentation of the parent module [pm_sampleVar](@ref pm_sampleVar) for algorithmic details and variance definitions.<br>
    !>
    !>  \param[out] var         :   The output object of rank `rank(sample) - 1` of type `real` of the same kind as the input `sample` representing its variance:<br>
    !>                              <ol>
    !>                                  <li>    If `sample` is a vector, the output `var` must be a scalar.<br>
    !>                                  <li>    If `sample` is a matrix, the output `var` must be a vector of size `ndim = size(sample, 3 - dim)`.<br>
    !>                              </ol>
    !>  \param[in]  mean        :   The input scalar or `contiguous` vector of shape `(ndim)` of the same type and kind as the input `sample` containing the `sample` mean.<br>
    !>                              <ol>
    !>                                  <li>    If the input `sample` is a 1D array, then `mean` must be a scalar.<br>
    !>                                  <li>    If the input `sample` is a 2D array, then `mean` must be a vector of size `ndim = size(sample, 3 - dim)`
    !>                                          (i.e., computed along the specified input dimension `dim`).<br>
    !>                              </ol>
    !>                              (**optional**. default = `0.` or [getFilled(0., ndim)](@ref pm_arrayFill::getFilled) depending on the rank of the input `sample` and value of `dim`.)
    !>  \param[in]  sample      :   The input `contiguous` array of shape `(nsam)`, `(ndim, nsam)`, or `(nsam, ndim)` of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the sample comprised of `nsam` observations each with `ndim` attributes.<br>
    !>                              If `sample` is a matrix, then the direction along which variance is computed is dictated by the input argument `dim`.<br>
    !>  \param[in]  dim         :   The input scalar `integer` of default kind \IK representing the dimension (`1` or `2`) of the input `sample` along which the mean must be computed.<br>
    !>                              <ol>
    !>                                  <li>    If `dim = 1`, the input `sample` of rank `2` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                                  <li>    If `dim = 2`, the input `sample` of rank `2` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                              </ol>
    !>                              The input `dim` must always be `1` or missing for an input `sample` of rank `1`.<br>
    !>                              (**optional**. If missing, the variance of the whole input `sample` is computed.)
    !>  \param[in]  weight      :   The input `contiguous` vector of length `nsam` of
    !>                              <ol>
    !>                                  <li>    type `integer` of default kind \IK, or
    !>                                  <li>    type `real` of the same kind as the input `sample`,
    !>                              </ol>
    !>                              containing the corresponding weights of individual `nsam` observations in `sample`.<br>
    !>                              (**optional**. default = `getFilled(1, size(sample, dim))` if `dim` is present, or `getFilled(1, size(sample))` if `dim` is missing.)
    !>  \param[in]  weisum      :   The input scalar of the same type and kind as the input `weight` containing `sum(weight)`.<br>
    !>                              This quantity is also a byproduct of computing the mean of a sample and is automatically returned by [setMean()](@ref pm_sampleMean::setMean).<br>
    !>                              (**optional**. It **must** be present **if and only if** the input argument `weight` is also present.)
    !>  \param[in]  correction  :   The input scalar object that can be the following:<br>
    !>                              <ol>
    !>                                  <li>    The constant [fweight](@ref pm_sampleWeight::fweight) or an object of type [fweight](@ref pm_sampleWeight::fweight)
    !>                                          implying a bias correction based on the assumption of **frequency weights** for the sample observations, even if the `weight` argument is missing.<br>
    !>                                          This is the most popular type of correction, also known as the [Bessel correction](https://en.wikipedia.org/wiki/Bessel%27s_correction).<br>
    !>                                  <li>    The constant [rweight](@ref pm_sampleWeight::rweight) or an object of type [rweight_type](@ref pm_sampleWeight::rweight_type)
    !>                                          implying a bias correction based on the assumption of **reliability weights** for the sample observations.<br>
    !>                              </ol>
    !>                              (**optional**. If missing, no bias-correction will be applied to the output `var`.)
    !>
    !>
    !>  \interface{setVar}
    !>  \code{.F90}
    !>
    !>      use pm_sampleVar, only: setVar
    !>
    !>      ! 1D sample unweighted
    !>
    !>      call setVar(var, sample(1 : nsam))
    !>      call setVar(var, mean, sample(1 : nsam))
    !>
    !>      call setVar(var, sample(1 : nsam), dim)
    !>      call setVar(var, mean, sample(1 : nsam), dim)
    !>
    !>      ! 1D sample weighted
    !>
    !>      call setVar(var, sample(1 : nsam), weight(1 : nsam))
    !>      call setVar(var, mean, sample(1 : nsam), weight(1 : nsam), weisum)
    !>
    !>      call setVar(var, sample(1 : nsam), dim, weight(1 : nsam))
    !>      call setVar(var, mean, sample(1 : nsam), dim, weight(1 : nsam), weisum)
    !>
    !>      ! 2D sample unweighted
    !>
    !>      call setVar(var(1 : size(sample, 3 - dim)), sample(:,:), dim)
    !>      call setVar(var(1 : size(sample, 3 - dim)), mean(1 : size(sample, 3 - dim)), sample(:,:), dim)
    !>
    !>      ! 2D sample weighted
    !>
    !>      call setVar(var(1 : size(sample, 3 - dim)), sample(:,:), dim, weight(1 : size(sample, dim)))
    !>      call setVar(var(1 : size(sample, 3 - dim)), mean(1 : size(sample, 3 - dim)), sample(:,:), dim, weight(1 : size(sample, dim)), weisum)
    !>
    !>      ! 2D whole sample unweighted
    !>
    !>      call setVar(var, sample(:,:))
    !>      call setVar(var, mean, sample(:,:))
    !>
    !>      ! 2D whole sample weighted
    !>
    !>      call setVar(var, sample(:,:), weight(1 : size(sample)))
    !>      call setVar(var, mean, sample(:,:), dim, weight(1 : size(sample)), weisum)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < sum(weight)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0. <= weight)` must hold for the corresponding input arguments.<br>
    !>  The condition `1 < size(sample, dim)` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < size(sample, 3 - dim)` must hold for the corresponding input arguments.<br>
    !>  The condition `1 <= dim .and. dim <= rank(sample)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample, 3 - dim) == size(var)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample, 3 - dim) == size(mean)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample, dim) == size(weight)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  Note the effects of bias-correction in computing the variance become noticeable only for very sample sample sizes (i.e., when `nsam` is small).<br>
    !>
    !>  \note
    !>  For a one-dimensional `sample`, one can also use the concise Fortran syntax to achieve the same goal as with
    !>  the interface `var = setVar(sample(:), mean, correction)` with integer `weight`:<br>
    !>  \code{.F90}
    !>
    !>       mean = sum(sample) / size(sample)
    !>       var = sum( (sample-mean)**2 ) / (size(sample) - 1)
    !>
    !>  \endcode
    !>  But the above concise version will be likely slightly slower as it potentially involves more loops.<br>
    !>
    !>  \note
    !>  For a two or higher-dimensional `sample`, if the variance is to be computed for the entire `sample` (as opposed to computing it along
    !>  a particular dimension), simply pass `reshape(sample, shape = size(sample))` to the appropriate [setVar](@ref pm_sampleVar::setVar) interface.<br>
    !>  Alternatively, a 1D pointer of the same size as the multidimensional sample can be passed to the procedure.<br>
    !>
    !>  \see
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [setVarMean](@ref pm_sampleVar::setVarMean)<br>
    !>  [getShifted](@ref pm_sampleShift::getShifted)<br>
    !>  [setShifted](@ref pm_sampleShift::setShifted)<br>
    !>
    !>  \example{setVar}
    !>  \include{lineno} example/pm_sampleVar/setVar/main.F90
    !>  \compilef{setVar}
    !>  \output{setVar}
    !>  \include{lineno} example/pm_sampleVar/setVar/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleVar](@ref test_pm_sampleVar)
    !>
    !>  \todo
    !>  \pmed
    !>  The examples of this generic interface need to be extended to weighted samples.<br>
    !>
    !>  \todo
    !>  \pvlow
    !>  The functionality of this interface can be expanded in the future to include the computation of
    !>  the variance of higher dimensional input `sample` and whole `sample` input arrays of arbitrary shape.<br>
    !>
    !>  \final{setVar}
    !>
    !>  \author
    !>  \FatemehBagheri, Monday 02:15 AM, September 27, 2021, Dallas, TX<br>

    ! AvgDIM_WNO

    interface setVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D1_CK5(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D1_CK4(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D1_CK3(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D1_CK2(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D1_CK1(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D1_RK5(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D1_RK4(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D1_RK3(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D1_RK2(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D1_RK1(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D2_CK5(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D2_CK4(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D2_CK3(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D2_CK2(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D2_CK1(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D2_RK5(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D2_RK4(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D2_RK3(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D2_RK2(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarAvgDIM_WNO_D2_RK1(var, mean, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WNO_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setVar

    ! AvgDIM_WTI

    interface setVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D1_CK5(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D1_CK4(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D1_CK3(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D1_CK2(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D1_CK1(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D1_RK5(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D1_RK4(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D1_RK3(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D1_RK2(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D1_RK1(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D2_CK5(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D2_CK4(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D2_CK3(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D2_CK2(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D2_CK1(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D2_RK5(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D2_RK4(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D2_RK3(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D2_RK2(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarAvgDIM_WTI_D2_RK1(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTI_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setVar

    ! AvgDIM_WTR

    interface setVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D1_CK5(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D1_CK4(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D1_CK3(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D1_CK2(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D1_CK1(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D1_RK5(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D1_RK4(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D1_RK3(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D1_RK2(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D1_RK1(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D2_CK5(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D2_CK4(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D2_CK3(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D2_CK2(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D2_CK1(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D2_RK5(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D2_RK4(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D2_RK3(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D2_RK2(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarAvgDIM_WTR_D2_RK1(var, mean, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgDIM_WTR_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)    , contiguous                :: mean(:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setVar

    ! OrgDIM_WNO

    interface setVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D1_CK5(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D1_CK4(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D1_CK3(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D1_CK2(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D1_CK1(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D1_RK5(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D1_RK4(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D1_RK3(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D1_RK2(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D1_RK1(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D2_CK5(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D2_CK4(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D2_CK3(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D2_CK2(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D2_CK1(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                :: dim
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D2_RK5(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D2_RK4(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D2_RK3(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D2_RK2(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarOrgDIM_WNO_D2_RK1(var, sample, dim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WNO_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setVar

    ! OrgDIM_WTI

    interface setVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D1_CK5(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D1_CK4(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D1_CK3(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D1_CK2(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D1_CK1(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D1_RK5(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D1_RK4(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D1_RK3(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D1_RK2(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D1_RK1(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D2_CK5(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D2_CK4(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D2_CK3(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D2_CK2(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D2_CK1(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D2_RK5(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D2_RK4(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D2_RK3(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D2_RK2(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarOrgDIM_WTI_D2_RK1(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTI_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                :: dim
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setVar

    ! OrgDIM_WTR

    interface setVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D1_CK5(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D1_CK4(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D1_CK3(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D1_CK2(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D1_CK1(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D1_RK5(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D1_RK4(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D1_RK3(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D1_RK2(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D1_RK1(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D2_CK5(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D2_CK4(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D2_CK3(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D2_CK2(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D2_CK1(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D2_RK5(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D2_RK4(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D2_RK3(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D2_RK2(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarOrgDIM_WTR_D2_RK1(var, sample, dim, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgDIM_WTR_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                :: dim
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! AvgALL_WNO

    interface setVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D1_CK5(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D1_CK4(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D1_CK3(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D1_CK2(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D1_CK1(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D1_RK5(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D1_RK4(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D1_RK3(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D1_RK2(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D1_RK1(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D2_CK5(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D2_CK4(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D2_CK3(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D2_CK2(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D2_CK1(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D2_RK5(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D2_RK4(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D2_RK3(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D2_RK2(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarAvgALL_WNO_D2_RK1(var, mean, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WNO_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setVar

    ! AvgALL_WTI

    interface setVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D1_CK5(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D1_CK4(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D1_CK3(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D1_CK2(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D1_CK1(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D1_RK5(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D1_RK4(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D1_RK3(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D1_RK2(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D1_RK1(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D2_CK5(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D2_CK4(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D2_CK3(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D2_CK2(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D2_CK1(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D2_RK5(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D2_RK4(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D2_RK3(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D2_RK2(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarAvgALL_WTI_D2_RK1(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTI_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setVar

    ! AvgALL_WTR

    interface setVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D1_CK5(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D1_CK4(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D1_CK3(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D1_CK2(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D1_CK1(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D1_RK5(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D1_RK4(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D1_RK3(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D1_RK2(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D1_RK1(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D2_CK5(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D2_CK4(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D2_CK3(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D2_CK2(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D2_CK1(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)            , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D2_RK5(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D2_RK4(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D2_RK3(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D2_RK2(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarAvgALL_WTR_D2_RK1(var, mean, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarAvgALL_WTR_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(in)                                :: mean
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setVar

    ! OrgALL_WNO

    interface setVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D1_CK5(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D1_CK4(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D1_CK3(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D1_CK2(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D1_CK1(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D1_RK5(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D1_RK4(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D1_RK3(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D1_RK2(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D1_RK1(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D2_CK5(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D2_CK4(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D2_CK3(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D2_CK2(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D2_CK1(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D2_RK5(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D2_RK4(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D2_RK3(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D2_RK2(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarOrgALL_WNO_D2_RK1(var, sample)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WNO_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setVar

    ! OrgALL_WTI

    interface setVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D1_CK5(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D1_CK4(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D1_CK3(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D1_CK2(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D1_CK1(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D1_RK5(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D1_RK4(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D1_RK3(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D1_RK2(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D1_RK1(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D2_CK5(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D2_CK4(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D2_CK3(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D2_CK2(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D2_CK1(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D2_RK5(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D2_RK4(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D2_RK3(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D2_RK2(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarOrgALL_WTI_D2_RK1(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTI_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                                :: weisum
        integer(IK)             , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setVar

    ! OrgALL_WTR

    interface setVar

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D1_CK5(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D1_CK4(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D1_CK3(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D1_CK2(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D1_CK1(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D1_RK5(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D1_RK4(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D1_RK3(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D1_RK2(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D1_RK1(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D2_CK5(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D2_CK4(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D2_CK3(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D2_CK2(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D2_CK1(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        complex(TKG)            , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D2_RK5(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D2_RK4(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D2_RK3(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D2_RK2(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarOrgALL_WTR_D2_RK1(var, sample, weight, weisum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarOrgALL_WTR_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)               , intent(in)                                :: weisum
        real(TKG)               , intent(in)    , contiguous                :: weight(:)
        real(TKG)               , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)               , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the (weighted) sample variance and mean of an input time series of `nsam` observations,
    !>  or of an input `sample` of `nsam` observations with `ndim` attributes optionally weighted by the input `weight`,
    !>  optionally also `sum(weight)`.<br>
    !>
    !>  \details
    !>  This generic interface is developed to specifically improve the performance of other procedures
    !>  where `sum(weight)` are also needed in addition to the mean of the weighted sample.<br>
    !>  It uses a **one-pass** approach to simultaneously compute the sample mean and variance.<br>
    !>
    !>  \param[out] var         :   The output object of rank `rank(sample) - 1` of type `real` of the same kind as the input `sample` representing its variance:<br>
    !>                                  <ol>
    !>                                      <li>    If `sample` is a vector, the output `var` must be a scalar.<br>
    !>                                      <li>    If `sample` is a matrix, the output `var` must be a vector of size `ndim = size(sample, 3 - dim)`.<br>
    !>                                  </ol>
    !>  \param[out] mean        :   The output object of the same type, kind, and rank as the output `var` containing the `sample` mean.<br>
    !>  \param[in]  sample      :   The input `contiguous` array of shape `(nsam)`, `(ndim, nsam)` or `(nsam, ndim)` of the same type and kind as the output `var`,
    !>                              containing the sample whose `mean` is to be computed.<br>
    !>  \param[in]  dim         :   The input scalar `integer` of default kind \IK representing the dimension (`1` or `2`) of the input `sample` along which the mean must be computed.<br>
    !>                              <ol>
    !>                                  <li>    If `dim = 1`, the input `sample` of rank `2` is assumed to have the shape `(nsam, ndim)`.<br>
    !>                                  <li>    If `dim = 2`, the input `sample` of rank `2` is assumed to have the shape `(ndim, nsam)`.<br>
    !>                              </ol>
    !>                              The input `dim` must always be `1` or missing for an input `sample` of rank `1`.<br>
    !>                              (**optional**. If missing, the variance of the whole input `sample` is computed.)
    !>  \param[in]  weight      :   The input `contiguous` vector of length `nsam` of,
    !>                              <ol>
    !>                                  <li>    type `integer` of default kind \IK, or
    !>                                  <li>    type `real` of the same kind as the input `sample`,
    !>                              </ol>
    !>                              containing the corresponding weights of the data points in the input `sample`.<br>
    !>  \param[out] weisum      :   The output scalar of the same type and kind as the input `weight`, containing `sum(weight)`.<br>
    !>                              (**optional**. It must be present **if and only if** the input argument `weight` is also present.)
    !>  \param[in] meang        :   The input object of the same type, kind, and rank as the output `var` containing the best guess for the `sample` mean.<br>
    !>                              If no good guess is known a priori, `meang` can be set to any observation in `sample`.<br>
    !>                              See the example usage below.<br>
    !>
    !>  \interface{setVarMean}
    !>  \code{.F90}
    !>
    !>      use pm_sampleVar, only: setVarMean
    !>
    !>      ! 1D sample
    !>
    !>      call setVarMean(var, mean, sample(1 : nsam), meang)
    !>      call setVarMean(var, mean, sample(1 : nsam), weight(1 : nsam), weisum, meang)
    !>
    !>      call setVarMean(var, mean, sample(1 : nsam), dim, meang)
    !>      call setVarMean(var, mean, sample(1 : nsam), dim, weight(1 : nsam), weisum, meang)
    !>
    !>      ! 2D sample
    !>
    !>      call setVarMean(var, mean(1 : ndim), sample(:,:), meang(1 : ndim))
    !>      call setVarMean(var, mean(1 : ndim), sample(:,:), weight(1 : size(sample)), weisum, meang(1 : ndim))
    !>
    !>      call setVarMean(var, mean(1 : ndim), sample(:,:), dim, meang(1 : ndim))
    !>      call setVarMean(var, mean(1 : ndim), sample(:,:), dim, weight(1 : size(sample, dim)), weisum, meang(1 : ndim))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `all(0. <= weight)` must hold for the corresponding input arguments.<br>
    !>  The condition `1 <= dim .and. dim <= rank(sample)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample, dim) == size(weight, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample, 3 - dim) == size(var, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample, 3 - dim) == size(mean, 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sample, 3 - dim) == size(meang, 1)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  If the input sample is to be an array of type `integer`, simply convert the sample to
    !>  an array of type `real` of the desired kind for the output `real` mean of the sample.<br>
    !>  There is no point in accepting an input sample of type `integer` since it will be inevitably
    !>  converted to an array of type `real` within the procedure to avoid potential integer overflow.<br>
    !>  Furthermore, an `sample` of type `integer` creates ambiguity about the `kind` of the `real`-valued returned mean by the procedure.<br>
    !>  See the notes in the description of the [pm_sampleMean](@ref pm_sampleMean).<br>
    !>
    !>  \note
    !>  Note that the mean of any one or two-dimensional sample can be simply computed via the Fortran intrinsic routine `sum()`:
    !>  \code{.F90}
    !>      integer(IK)                 :: i
    !>      integer(IK) , parameter     :: NDIM = 3_IK
    !>      integer(IK) , parameter     :: NSAM = 1000_IK
    !>      real(TKG)   , parameter     :: sample(NDIM,NSAM) = reshape([( real(i,RK), i = 1, NSAM )], shape = shape(sample))
    !>      real(TKG)   , allocatable   :: mean(:)
    !>      mean = sum(sample, dim = 1) / size(transpose(sample), dim = 1)  ! assuming the first dimension represents the observations
    !>      mean = sum(sample, dim = 2) / size(sample, dim = 2)             ! assuming the second dimension represents the observations
    !>  \endcode
    !>
    !>  \note
    !>  The mean of a whole multidimensional array can be obtained by either,
    !>  <ol>
    !>      <li>    reshaping the array to a vector form and passing it to this procedure, or
    !>      <li>    mapping the array to a 1-dimensional pointer array of the same size as the `ndim` dimensional array.<br>
    !>  </ol>
    !>  See the examples below.<br>
    !>
    !>  \see
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [setCovMean](@ref pm_sampleCov::setCovMean)<br>
    !>
    !>  \example{setVarMean}
    !>  \include{lineno} example/pm_sampleVar/setVarMean/main.F90
    !>  \compilef{setVarMean}
    !>  \output{setVarMean}
    !>  \include{lineno} example/pm_sampleVar/setVarMean/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleMean](@ref test_pm_sampleMean)
    !>
    !>  \final{setVarMean}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! ALL D1

    interface setVarMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D1_CK5(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D1_CK4(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D1_CK3(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D1_CK2(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D1_CK1(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D1_RK5(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D1_RK4(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D1_RK3(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D1_RK2(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D1_RK1(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D1_CK5(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D1_CK4(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D1_CK3(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D1_CK2(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D1_CK1(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D1_RK5(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D1_RK4(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D1_RK3(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D1_RK2(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D1_RK1(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D1_CK5(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D1_CK4(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D1_CK3(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D1_CK2(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D1_CK1(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D1_RK5(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D1_RK4(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D1_RK3(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D1_RK2(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D1_RK1(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! ALL D2

    interface setVarMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D2_CK5(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D2_CK4(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D2_CK3(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D2_CK2(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D2_CK1(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D2_RK5(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D2_RK4(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D2_RK3(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D2_RK2(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMeanALL_WNODD_D2_RK1(var, mean, sample, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WNODD_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D2_CK5(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D2_CK4(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D2_CK3(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D2_CK2(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D2_CK1(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D2_RK5(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D2_RK4(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D2_RK3(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D2_RK2(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMeanALL_WTISD_D2_RK1(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTISD_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D2_CK5(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D2_CK4(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D2_CK3(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D2_CK2(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D2_CK1(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D2_RK5(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D2_RK4(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D2_RK3(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D2_RK2(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMeanALL_WTRSD_D2_RK1(var, mean, sample, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanALL_WTRSD_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! DIM D1

    interface setVarMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D1_CK5(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D1_CK4(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D1_CK3(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D1_CK2(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D1_CK1(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D1_RK5(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D1_RK4(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D1_RK3(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D1_RK2(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D1_RK1(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D1_CK5(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D1_CK4(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D1_CK3(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D1_CK2(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D1_CK1(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D1_RK5(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D1_RK4(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D1_RK3(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D1_RK2(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D1_RK1(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D1_CK5(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D1_CK4(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D1_CK3(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D1_CK2(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D1_CK1(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:)
        complex(TKG)    , intent(in)                                :: meang
        complex(TKG)    , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D1_RK5(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D1_RK4(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D1_RK3(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D1_RK2(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D1_RK1(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:)
        real(TKG)       , intent(in)                                :: meang
        real(TKG)       , intent(out)                               :: mean
        real(TKG)       , intent(out)                               :: var
        integer(IK)     , intent(in)                                :: dim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! DIM D2

    interface setVarMean

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D2_CK5(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)    , contiguous                :: meang(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D2_CK4(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)    , contiguous                :: meang(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D2_CK3(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)    , contiguous                :: meang(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D2_CK2(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)    , contiguous                :: meang(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D2_CK1(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(in)                                :: dim
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)    , contiguous                :: meang(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D2_RK5(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: meang(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D2_RK4(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: meang(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D2_RK3(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: meang(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D2_RK2(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: meang(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMeanDIM_WNODD_D2_RK1(var, mean, sample, dim, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WNODD_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: meang(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D2_CK5(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(in)                                :: dim
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)    , contiguous                :: meang(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D2_CK4(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(in)                                :: dim
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)    , contiguous                :: meang(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D2_CK3(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(in)                                :: dim
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)    , contiguous                :: meang(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D2_CK2(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(in)                                :: dim
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)    , contiguous                :: meang(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D2_CK1(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(in)                                :: dim
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)    , contiguous                :: meang(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D2_RK5(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(in)                                :: dim
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: meang(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D2_RK4(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(in)                                :: dim
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: meang(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D2_RK3(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(in)                                :: dim
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: meang(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D2_RK2(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(in)                                :: dim
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: meang(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMeanDIM_WTISD_D2_RK1(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTISD_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(in)                                :: dim
        integer(IK)     , intent(out)                               :: weisum
        integer(IK)     , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: meang(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D2_CK5(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D2_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)    , contiguous                :: meang(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D2_CK4(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D2_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)    , contiguous                :: meang(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D2_CK3(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D2_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)    , contiguous                :: meang(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D2_CK2(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D2_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)    , contiguous                :: meang(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D2_CK1(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D2_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        complex(TKG)    , intent(in)    , contiguous                :: sample(:,:)
        complex(TKG)    , intent(in)    , contiguous                :: meang(:)
        complex(TKG)    , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D2_RK5(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D2_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: meang(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D2_RK4(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D2_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: meang(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D2_RK3(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D2_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: meang(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D2_RK2(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D2_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: meang(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMeanDIM_WTRSD_D2_RK1(var, mean, sample, dim, weight, weisum, meang)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanDIM_WTRSD_D2_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)     , intent(in)                                :: dim
        real(TKG)       , intent(out)                               :: weisum
        real(TKG)       , intent(in)    , contiguous                :: weight(:)
        real(TKG)       , intent(in)    , contiguous                :: sample(:,:)
        real(TKG)       , intent(in)    , contiguous                :: meang(:)
        real(TKG)       , intent(out)   , contiguous                :: mean(:)
        real(TKG)       , intent(out)   , contiguous                :: var(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the (weighted) merged variance of a `complex` or `real` sample resulting from the merger of two separate (weighted) samples \f$A\f$ and \f$B\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_sampleVar](@ref pm_sampleVar) for more information and definition online updating of sample variance.<br>
    !>
    !>  \note
    !>  The input and output variances of this generic interface are all **biased** variances.<br>
    !>  A **biased** variance can be readily unbiased by multiplying it with the appropriate bias-correction factors
    !>  detailed in the documentation of [pm_sampleVar](@ref pm_sampleVar).<br>
    !>
    !>  \param[in]      varB        :   The input object of the same type and kind and rank and shape as the input argument `varA`,
    !>                                  containing the **biased** variance of the second sample that must be merged with the first sample.<br>
    !>  \param[in]      varA        :   The input scalar or vector of shape `(1:ndim)` of type `real` of the same kind as that of the input argument `meanDiff`,
    !>                                  containing the **biased** variance of the first sample that must be merged with the second sample.<br>
    !>  \param[in]      meanDiff    :   The input object of,
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the difference of mean of the two samples `meanDiff = meanA - meanB`.<br>
    !>                                  The subtraction order is immaterial.<br>
    !>  \param[in]      fracA       :   The input scalar of type `real` of the same kind as kind of `varA`,
    !>                                  containing the sum of the weights of all points in sample \f$A\f$ divided by the sum of the weights of all points in the merged sample.<br>
    !>                                  If the sample is unweighted, then `fracA` is simply `size(sampleA) / (size(sampleA) + size(sampleB))`.<br>
    !>
    !>  \return
    !>  `varMerged`                 :   The output object of the same type and kind and rank and shape as `meanDiff`,
    !>                                  containing the **biased** variance of the sample resulting form the merger of the two samples.<br>
    !>
    !>  \interface{getVarMerged}
    !>  \code{.F90}
    !>
    !>      use pm_sampleVar, only: getVarMerged
    !>
    !>      ! univariate sample
    !>
    !>      varMerged = getVarMerged(varB, varA, meanDiff, fracA)
    !>
    !>      ! ndim-dimensional sample
    !>
    !>      varMerged = getVarMerged(varB(1:ndim), varA(1:ndim), meanDiff(1:ndim), fracA)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < fracA .and. fracA < 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(varB) == size(varA)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(varA) == size(meanDiff)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getCovMerged](@ref pm_sampleCov::getCovMerged)<br>
    !>  [setCovMerged](@ref pm_sampleCov::setCovMerged)<br>
    !>  [getVarMerged](@ref pm_sampleVar::getVarMerged)<br>
    !>  [setVarMerged](@ref pm_sampleVar::setVarMerged)<br>
    !>  [getMeanMerged](@ref pm_sampleMean::getMeanMerged)<br>
    !>  [setMeanMerged](@ref pm_sampleMean::setMeanMerged)<br>
    !>  [setVarMeanMerged](@ref pm_sampleVar::setVarMeanMerged)<br>
    !>
    !>  \example{getVarMerged}
    !>  \include{lineno} example/pm_sampleVar/getVarMerged/main.F90
    !>  \compilef{getVarMerged}
    !>  \output{getVarMerged}
    !>  \include{lineno} example/pm_sampleVar/getVarMerged/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleVar](@ref test_pm_sampleVar)
    !>
    !>  \final{getVarMerged}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! D0

    interface getVarMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getVarMergedNew_D0_CK5(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D0_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)                                                   :: varMerged
    end function
#endif

#if CK4_ENABLED
    PURE module function getVarMergedNew_D0_CK4(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D0_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)                                                   :: varMerged
    end function
#endif

#if CK3_ENABLED
    PURE module function getVarMergedNew_D0_CK3(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D0_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)                                                   :: varMerged
    end function
#endif

#if CK2_ENABLED
    PURE module function getVarMergedNew_D0_CK2(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D0_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)                                                   :: varMerged
    end function
#endif

#if CK1_ENABLED
    PURE module function getVarMergedNew_D0_CK1(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D0_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)                                                   :: varMerged
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getVarMergedNew_D0_RK5(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D0_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)                                                   :: varMerged
    end function
#endif

#if RK4_ENABLED
    PURE module function getVarMergedNew_D0_RK4(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D0_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)                                                   :: varMerged
    end function
#endif

#if RK3_ENABLED
    PURE module function getVarMergedNew_D0_RK3(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D0_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)                                                   :: varMerged
    end function
#endif

#if RK2_ENABLED
    PURE module function getVarMergedNew_D0_RK2(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D0_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)                                                   :: varMerged
    end function
#endif

#if RK1_ENABLED
    PURE module function getVarMergedNew_D0_RK1(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D0_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)                                                   :: varMerged
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! D1

    interface getVarMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getVarMergedNew_D1_CK5(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)                                                   :: varMerged(size(meanDiff, 1, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getVarMergedNew_D1_CK4(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)                                                   :: varMerged(size(meanDiff, 1, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getVarMergedNew_D1_CK3(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)                                                   :: varMerged(size(meanDiff, 1, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getVarMergedNew_D1_CK2(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)                                                   :: varMerged(size(meanDiff, 1, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getVarMergedNew_D1_CK1(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)                                                   :: varMerged(size(meanDiff, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getVarMergedNew_D1_RK5(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)                                                   :: varMerged(size(meanDiff, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getVarMergedNew_D1_RK4(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)                                                   :: varMerged(size(meanDiff, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getVarMergedNew_D1_RK3(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)                                                   :: varMerged(size(meanDiff, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getVarMergedNew_D1_RK2(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)                                                   :: varMerged(size(meanDiff, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getVarMergedNew_D1_RK1(varB, varA, meanDiff, fracA) result(varMerged)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVarMergedNew_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)                                                   :: varMerged(size(meanDiff, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the (weighted) merged variance of a `complex` or `real` sample resulting from the merger of two separate (weighted) samples \f$A\f$ and \f$B\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_sampleVar](@ref pm_sampleVar) for more information and definition online updating of sample variance.<br>
    !>
    !>  \note
    !>  The input and output variances of this generic interface are all **biased** variances.<br>
    !>  A **biased** variance can be readily unbiased by multiplying it with the appropriate bias-correction factors
    !>  detailed in the documentation of [pm_sampleVar](@ref pm_sampleVar).<br>
    !>
    !>  \param[out]     varMerged   :   The output object of the same type and kind and rank and shape as `varA`,
    !>                                  containing the **biased** variance of the sample resulting form the merger of the two samples.<br>
    !>                                  (**optional**. If missing, the resulting merged variance will be written to the argument `varB`.)
    !>  \param[inout]   varB        :   The input or input/output object of the same type and kind and rank and shape as the input argument `varA`,
    !>                                  containing the **biased** variance of the second sample that must be merged with the first sample.<br>
    !>                                  If the input argument `varMerged` is missing, then `varB` contains the updated **biased** variance of the merged sample on return.<br>
    !>                                  Otherwise, the contents of `varB` remain intact upon return.<br>
    !>  \param[in]      varA        :   The input scalar or `contiguous` vector of shape `(1:ndim)` of type `real` of the same kind as the kind of the input argument `meanDiff`,
    !>                                  containing the **biased** variance of the first sample that must be merged with the second sample.<br>
    !>  \param[in]      meanDiff    :   The input object of
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the difference of mean of the two samples `meanDiff = meanA - meanB`.<br>
    !>                                  The subtraction order is immaterial.<br>
    !>  \param[in]      fracA       :   The input scalar of type `real` of the same kind as kind of `varA`,
    !>                                  containing the ratio of the sum of the weights of all points in sample \f$A\f$ to sum of weights of all points in the merged sample.<br>
    !>                                  If the sample is unweighted, then `fracA` is simply `size(sampleA) / (size(sampleA) + size(sampleB))`.<br>
    !>
    !>  \interface{setVarMerged}
    !>  \code{.F90}
    !>
    !>      use pm_sampleVar, only: setVarMerged
    !>
    !>      ! univariate sample
    !>
    !>      call setVarMerged(           varB, varA, meanDiff, fracA)
    !>      call setVarMerged(varMerged, varB, varA, meanDiff, fracA)
    !>
    !>      ! ndim-dimensional sample
    !>
    !>      call setVarMerged(                   varB(1:ndim), varA(1:ndim), meanDiff(1:ndim), fracA)
    !>      call setVarMerged(varMerged(1:ndim), varB(1:ndim), varA(1:ndim), meanDiff(1:ndim), fracA)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < fracA .and. fracA < 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(varB) == size(varA)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(varA) == size(meanDiff)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(varA) == size(varMearged)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(varA) == size(meanMearged)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getCovMerged](@ref pm_sampleCov::getCovMerged)<br>
    !>  [setCovMerged](@ref pm_sampleCov::setCovMerged)<br>
    !>  [getVarMerged](@ref pm_sampleVar::getVarMerged)<br>
    !>  [setVarMerged](@ref pm_sampleVar::setVarMerged)<br>
    !>  [getMeanMerged](@ref pm_sampleMean::getMeanMerged)<br>
    !>  [setMeanMerged](@ref pm_sampleMean::setMeanMerged)<br>
    !>  [setVarMeanMerged](@ref pm_sampleVar::setVarMeanMerged)<br>
    !>
    !>  \example{setVarMerged}
    !>  \include{lineno} example/pm_sampleVar/setVarMerged/main.F90
    !>  \compilef{setVarMerged}
    !>  \output{setVarMerged}
    !>  \include{lineno} example/pm_sampleVar/setVarMerged/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleVar](@ref test_pm_sampleVar)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! New_D0

    interface setVarMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMergedNew_D0_CK5(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D0_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)       , intent(out)                               :: varMerged
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMergedNew_D0_CK4(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D0_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)       , intent(out)                               :: varMerged
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMergedNew_D0_CK3(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D0_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)       , intent(out)                               :: varMerged
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMergedNew_D0_CK2(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D0_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)       , intent(out)                               :: varMerged
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMergedNew_D0_CK1(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D0_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)       , intent(out)                               :: varMerged
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMergedNew_D0_RK5(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D0_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)       , intent(out)                               :: varMerged
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMergedNew_D0_RK4(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D0_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)       , intent(out)                               :: varMerged
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMergedNew_D0_RK3(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D0_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)       , intent(out)                               :: varMerged
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMergedNew_D0_RK2(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D0_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)       , intent(out)                               :: varMerged
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMergedNew_D0_RK1(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D0_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: varB, varA
        real(TKG)       , intent(out)                               :: varMerged
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! New_D1

    interface setVarMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMergedNew_D1_CK5(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMergedNew_D1_CK4(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMergedNew_D1_CK3(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMergedNew_D1_CK2(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMergedNew_D1_CK1(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(in)                                :: fracA
        complex(TKG)    , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMergedNew_D1_RK5(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMergedNew_D1_RK4(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMergedNew_D1_RK3(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMergedNew_D1_RK2(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMergedNew_D1_RK1(varMerged, varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedNew_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:), varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Old_D0

    interface setVarMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMergedOld_D0_CK5(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D0_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(inout)                             :: varB
        real(TKG)       , intent(in)                                :: varA
        complex(TKG)    , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMergedOld_D0_CK4(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D0_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(inout)                             :: varB
        real(TKG)       , intent(in)                                :: varA
        complex(TKG)    , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMergedOld_D0_CK3(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D0_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(inout)                             :: varB
        real(TKG)       , intent(in)                                :: varA
        complex(TKG)    , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMergedOld_D0_CK2(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D0_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(inout)                             :: varB
        real(TKG)       , intent(in)                                :: varA
        complex(TKG)    , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMergedOld_D0_CK1(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D0_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(inout)                             :: varB
        real(TKG)       , intent(in)                                :: varA
        complex(TKG)    , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMergedOld_D0_RK5(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D0_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(inout)                             :: varB
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMergedOld_D0_RK4(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D0_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(inout)                             :: varB
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMergedOld_D0_RK3(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D0_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(inout)                             :: varB
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMergedOld_D0_RK2(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D0_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(inout)                             :: varB
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMergedOld_D0_RK1(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D0_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(inout)                             :: varB
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(in)                                :: meanDiff
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Old_D1

    interface setVarMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMergedOld_D1_CK5(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(inout) , contiguous                :: varB(:)
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMergedOld_D1_CK4(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(inout) , contiguous                :: varB(:)
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMergedOld_D1_CK3(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(inout) , contiguous                :: varB(:)
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMergedOld_D1_CK2(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(inout) , contiguous                :: varB(:)
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMergedOld_D1_CK1(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(inout) , contiguous                :: varB(:)
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMergedOld_D1_RK5(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(inout) , contiguous                :: varB(:)
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMergedOld_D1_RK4(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(inout) , contiguous                :: varB(:)
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMergedOld_D1_RK3(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(inout) , contiguous                :: varB(:)
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMergedOld_D1_RK2(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(inout) , contiguous                :: varB(:)
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMergedOld_D1_RK1(varB, varA, meanDiff, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMergedOld_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(inout) , contiguous                :: varB(:)
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(in)    , contiguous                :: meanDiff(:)
        real(TKG)       , intent(in)                                :: fracA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the (weighted) merged variance and mean of a `complex` or `real` sample
    !>  resulting from the merger of two separate (weighted) samples \f$A\f$ and \f$B\f$.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_sampleVar](@ref pm_sampleVar) and [pm_sampleMean](@ref pm_sampleMean)
    !>  for more information and definition online updating of sample variance and mean.<br>
    !>
    !>  \note
    !>  The input and output variances of this generic interface are all **biased** variances.<br>
    !>  A **biased** variance can be readily unbiased by multiplying it with the appropriate bias-correction factors
    !>  detailed in the documentation of [pm_sampleVar](@ref pm_sampleVar).<br>
    !>
    !>  \param[out]     varMerged   :   The output object of the same type and kind and rank and shape as `meanA`,
    !>                                  containing the **biased** variance of the sample resulting form the merger of the two samples.<br>
    !>                                  (**optional**. If missing, the resulting merged variance will be written to the argument `varB`.)
    !>  \param[out]     meanMerged  :   The output object of the same type and kind and rank and shape as `meanA`,
    !>                                  containing the mean of the sample resulting form the merger of the two samples.<br>
    !>                                  (**optional**. If missing, the resulting merged mean will be written to the argument `meanB`.)
    !>  \param[inout]   varB        :   The input or input/output object of type `real` of the same kind and rank and shape as the input argument `meanA`,
    !>                                  containing the **biased** variance of the second sample that must be merged with the first sample.<br>
    !>                                  If the input argument `varMerged` is missing, then `varB` contains the updated **biased** variance of the merged sample on return.<br>
    !>                                  Otherwise, the contents of `varB` remain intact upon return.<br>
    !>  \param[inout]   meanB       :   The input or input/output object of the same type and kind and rank and shape as the input argument `meanA`,
    !>                                  containing the mean of the second sample that must be merged with the first sample.<br>
    !>                                  If the input argument `meanMerged` is missing, then `meanB` contains the updated mean of the merged sample on return.<br>
    !>                                  Otherwise, the contents of `meanB` remain intact upon return.<br>
    !>  \param[in]      varA        :   The input scalar or `contiguous` vector of shape `(1:ndim)` of type `real` of the same kind as the kind of the input argument `meanA`,
    !>                                  containing the **biased** variance of the first sample that must be merged with the second sample.<br>
    !>  \param[in]      meanA       :   The input scalar or `contiguous` vector of shape `(1:ndim)` of,
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the mean of the first sample that must be merged with the second sample.<br>
    !>  \param[in]      fracA       :   The input scalar of type `real` of the same kind as kind of `varA`,
    !>                                  containing the ratio of the sum of the weights of all points in sample \f$A\f$ to sum of weights of all points in the merged sample.<br>
    !>                                  If the sample is unweighted, then `fracA` is simply `size(sampleA) / (size(sampleA) + size(sampleB))`.<br>
    !>
    !>  \interface{setVarMeanMerged}
    !>  \code{.F90}
    !>
    !>      use pm_sampleVar, only: setVarMeanMerged
    !>
    !>      ! univariate sample
    !>
    !>      call setVarMeanMerged(varB, meanB, varA, meanA, fracA)
    !>      call setVarMeanMerged(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
    !>
    !>      ! ndim-dimensional sample
    !>
    !>      call setVarMeanMerged(varB(1:ndim), varA(1:ndim), meanDiff(1:ndim), fracA)
    !>      call setVarMeanMerged(varMerged(1:ndim), varB(1:ndim), varA(1:ndim), meanDiff(1:ndim), fracA)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < fracA .and. fracA < 1` must hold for the corresponding input arguments.<br>
    !>  The condition `size(varB) == size(varA)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(varA) == size(meanA)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(varA) == size(meanB)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(varA) == size(varMearged)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(varA) == size(meanMearged)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getCor](@ref pm_sampleCor::getCor)<br>
    !>  [setCor](@ref pm_sampleCor::setCor)<br>
    !>  [getCov](@ref pm_sampleCov::getCov)<br>
    !>  [setCov](@ref pm_sampleCov::setCov)<br>
    !>  [getVar](@ref pm_sampleVar::getVar)<br>
    !>  [setVar](@ref pm_sampleVar::setVar)<br>
    !>  [getMean](@ref pm_sampleMean::getMean)<br>
    !>  [setMean](@ref pm_sampleMean::setMean)<br>
    !>  [getCovMerged](@ref pm_sampleCov::getCovMerged)<br>
    !>  [setCovMerged](@ref pm_sampleCov::setCovMerged)<br>
    !>  [getVarMerged](@ref pm_sampleVar::getVarMerged)<br>
    !>  [setVarMerged](@ref pm_sampleVar::setVarMerged)<br>
    !>  [getMeanMerged](@ref pm_sampleMean::getMeanMerged)<br>
    !>  [setMeanMerged](@ref pm_sampleMean::setMeanMerged)<br>
    !>  [setVarMeanMerged](@ref pm_sampleVar::setVarMeanMerged)<br>
    !>  [setVarMeanMerged](@ref pm_sampleVar::setVarMeanMerged)<br>
    !>
    !>  \example{setVarMeanMerged}
    !>  \include{lineno} example/pm_sampleVar/setVarMeanMerged/main.F90
    !>  \compilef{setVarMeanMerged}
    !>  \output{setVarMeanMerged}
    !>  \include{lineno} example/pm_sampleVar/setVarMeanMerged/main.out.F90
    !>
    !>  \test
    !>  [test_pm_sampleVar](@ref test_pm_sampleVar)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! New_D0

    interface setVarMeanMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMeanMergedNew_D0_CK5(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D0_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(in)                                :: varB
        real(TKG)       , intent(out)                               :: varMerged
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(in)                                :: meanB
        complex(TKG)    , intent(out)                               :: meanMerged
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMeanMergedNew_D0_CK4(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D0_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(in)                                :: varB
        real(TKG)       , intent(out)                               :: varMerged
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(in)                                :: meanB
        complex(TKG)    , intent(out)                               :: meanMerged
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMeanMergedNew_D0_CK3(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D0_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(in)                                :: varB
        real(TKG)       , intent(out)                               :: varMerged
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(in)                                :: meanB
        complex(TKG)    , intent(out)                               :: meanMerged
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMeanMergedNew_D0_CK2(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D0_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(in)                                :: varB
        real(TKG)       , intent(out)                               :: varMerged
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(in)                                :: meanB
        complex(TKG)    , intent(out)                               :: meanMerged
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMeanMergedNew_D0_CK1(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D0_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(in)                                :: varB
        real(TKG)       , intent(out)                               :: varMerged
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(in)                                :: meanB
        complex(TKG)    , intent(out)                               :: meanMerged
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMeanMergedNew_D0_RK5(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D0_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(in)                                :: varB
        real(TKG)       , intent(out)                               :: varMerged
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(in)                                :: meanB
        real(TKG)       , intent(out)                               :: meanMerged
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMeanMergedNew_D0_RK4(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D0_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(in)                                :: varB
        real(TKG)       , intent(out)                               :: varMerged
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(in)                                :: meanB
        real(TKG)       , intent(out)                               :: meanMerged
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMeanMergedNew_D0_RK3(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D0_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(in)                                :: varB
        real(TKG)       , intent(out)                               :: varMerged
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(in)                                :: meanB
        real(TKG)       , intent(out)                               :: meanMerged
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMeanMergedNew_D0_RK2(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D0_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(in)                                :: varB
        real(TKG)       , intent(out)                               :: varMerged
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(in)                                :: meanB
        real(TKG)       , intent(out)                               :: meanMerged
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMeanMergedNew_D0_RK1(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D0_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(in)                                :: varB
        real(TKG)       , intent(out)                               :: varMerged
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(in)                                :: meanB
        real(TKG)       , intent(out)                               :: meanMerged
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! New_D1

    interface setVarMeanMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMeanMergedNew_D1_CK5(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanB(:)
        complex(TKG)    , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMeanMergedNew_D1_CK4(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanB(:)
        complex(TKG)    , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMeanMergedNew_D1_CK3(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanB(:)
        complex(TKG)    , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMeanMergedNew_D1_CK2(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanB(:)
        complex(TKG)    , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMeanMergedNew_D1_CK1(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanB(:)
        complex(TKG)    , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMeanMergedNew_D1_RK5(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
        real(TKG)       , intent(in)    , contiguous                :: meanA(:)
        real(TKG)       , intent(in)    , contiguous                :: meanB(:)
        real(TKG)       , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMeanMergedNew_D1_RK4(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
        real(TKG)       , intent(in)    , contiguous                :: meanA(:)
        real(TKG)       , intent(in)    , contiguous                :: meanB(:)
        real(TKG)       , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMeanMergedNew_D1_RK3(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
        real(TKG)       , intent(in)    , contiguous                :: meanA(:)
        real(TKG)       , intent(in)    , contiguous                :: meanB(:)
        real(TKG)       , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMeanMergedNew_D1_RK2(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
        real(TKG)       , intent(in)    , contiguous                :: meanA(:)
        real(TKG)       , intent(in)    , contiguous                :: meanB(:)
        real(TKG)       , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMeanMergedNew_D1_RK1(varMerged, meanMerged, varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedNew_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(in)    , contiguous                :: varB(:)
        real(TKG)       , intent(out)   , contiguous                :: varMerged(:)
        real(TKG)       , intent(in)    , contiguous                :: meanA(:)
        real(TKG)       , intent(in)    , contiguous                :: meanB(:)
        real(TKG)       , intent(out)   , contiguous                :: meanMerged(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Old_D0

    interface setVarMeanMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMeanMergedOld_D0_CK5(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D0_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(out)                               :: varB
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(out)                               :: meanB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMeanMergedOld_D0_CK4(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D0_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(out)                               :: varB
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(out)                               :: meanB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMeanMergedOld_D0_CK3(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D0_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(out)                               :: varB
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(out)                               :: meanB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMeanMergedOld_D0_CK2(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D0_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(out)                               :: varB
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(out)                               :: meanB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMeanMergedOld_D0_CK1(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D0_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(out)                               :: varB
        complex(TKG)    , intent(in)                                :: meanA
        complex(TKG)    , intent(out)                               :: meanB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMeanMergedOld_D0_RK5(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D0_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(out)                               :: varB
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(out)                               :: meanB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMeanMergedOld_D0_RK4(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D0_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(out)                               :: varB
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(out)                               :: meanB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMeanMergedOld_D0_RK3(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D0_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(out)                               :: varB
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(out)                               :: meanB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMeanMergedOld_D0_RK2(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D0_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(out)                               :: varB
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(out)                               :: meanB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMeanMergedOld_D0_RK1(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D0_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)                                :: varA
        real(TKG)       , intent(out)                               :: varB
        real(TKG)       , intent(in)                                :: meanA
        real(TKG)       , intent(out)                               :: meanB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Old_D1

    interface setVarMeanMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setVarMeanMergedOld_D1_CK5(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D1_CK5
#endif
        use pm_kind, only: TKG => CK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varB(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:)
        complex(TKG)    , intent(out)   , contiguous                :: meanB(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setVarMeanMergedOld_D1_CK4(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D1_CK4
#endif
        use pm_kind, only: TKG => CK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varB(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:)
        complex(TKG)    , intent(out)   , contiguous                :: meanB(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setVarMeanMergedOld_D1_CK3(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D1_CK3
#endif
        use pm_kind, only: TKG => CK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varB(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:)
        complex(TKG)    , intent(out)   , contiguous                :: meanB(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setVarMeanMergedOld_D1_CK2(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D1_CK2
#endif
        use pm_kind, only: TKG => CK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varB(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:)
        complex(TKG)    , intent(out)   , contiguous                :: meanB(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setVarMeanMergedOld_D1_CK1(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D1_CK1
#endif
        use pm_kind, only: TKG => CK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varB(:)
        complex(TKG)    , intent(in)    , contiguous                :: meanA(:)
        complex(TKG)    , intent(out)   , contiguous                :: meanB(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setVarMeanMergedOld_D1_RK5(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D1_RK5
#endif
        use pm_kind, only: TKG => RK5
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varB(:)
        real(TKG)       , intent(in)    , contiguous                :: meanA(:)
        real(TKG)       , intent(out)   , contiguous                :: meanB(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setVarMeanMergedOld_D1_RK4(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D1_RK4
#endif
        use pm_kind, only: TKG => RK4
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varB(:)
        real(TKG)       , intent(in)    , contiguous                :: meanA(:)
        real(TKG)       , intent(out)   , contiguous                :: meanB(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setVarMeanMergedOld_D1_RK3(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D1_RK3
#endif
        use pm_kind, only: TKG => RK3
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varB(:)
        real(TKG)       , intent(in)    , contiguous                :: meanA(:)
        real(TKG)       , intent(out)   , contiguous                :: meanB(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setVarMeanMergedOld_D1_RK2(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D1_RK2
#endif
        use pm_kind, only: TKG => RK2
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varB(:)
        real(TKG)       , intent(in)    , contiguous                :: meanA(:)
        real(TKG)       , intent(out)   , contiguous                :: meanB(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setVarMeanMergedOld_D1_RK1(varB, meanB, varA, meanA, fracA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setVarMeanMergedOld_D1_RK1
#endif
        use pm_kind, only: TKG => RK1
        real(TKG)       , intent(in)                                :: fracA
        real(TKG)       , intent(in)    , contiguous                :: varA(:)
        real(TKG)       , intent(out)   , contiguous                :: varB(:)
        real(TKG)       , intent(in)    , contiguous                :: meanA(:)
        real(TKG)       , intent(out)   , contiguous                :: meanB(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_sampleVar ! LCOV_EXCL_LINE