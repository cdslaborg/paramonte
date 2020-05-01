"targetAcceptanceRate sets an optimal target for the ratio of the number of accepted objective function calls to the &
&total number of function calls by the " // methodName // " sampler. It is a real-valued array of length 2, whose elements &
&determine the upper and lower bounds of the desired acceptance rate. When the acceptance rate of the sampler is outside the &
&specified limits, the sampler's settings will be automatically adjusted to bring the overall acceptance rate to within the &
&specified limits by the input variable targetAcceptanceRate. When assigned from within a dynamic-language programming &
&environment, such as MATLAB or Python, or from within an input file, targetAcceptanceRate can also be a single real number &
&between 0 and 1. In such case, the " // methodName // " sampler will constantly attempt (with no guarantee of success) &
&to bring the average acceptance ratio of the sampler as close to the user-provided target ratio as possible. The success &
&of " // methodName // " in keeping the average acceptance ratio close to the requested target value depends heavily on:\n&
&    1) the value of adaptiveUpdatePeriod; the larger, the easier.\n&
&    2) the value of adaptiveUpdateCount; the larger, the easier.\n&
&Note that the acceptance ratio adjustments will only occur every adaptiveUpdatePeriod sampling steps for a total &
&number of adaptiveUpdateCount. &
&There is no default value for targetAcceptanceRate, as the acceptance ratio is not directly adjusted during sampling."