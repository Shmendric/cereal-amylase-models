# Global certificates for the supplied nominal maltose curves

This is a self-contained specialization of classical uniform approximation and isotonic regression. No priority is claimed for the general mathematics. The project contribution is the executed, independently checkable analysis of the four specified technical datasets.

## Model and target

At ordered fixed times t_i >= 0 and nominal values y_i, consider

f(t) = P0 + Delta * (1 - exp(-d*t)), with P0 >= 0, Delta >= 0, d > 0.

The target epsilon_E is the infimum over this whole family of max_i abs(y_i-f(t_i)). There is no upper numerical bound on d in the definition. No independent chemical substrate-capacity constraint is imposed. Bounds for this enlarged model remain lower bounds if such a constraint is later added, but the witness need not be admissible for the smaller family.

## At most two intersections

The difference of two members is a + b*exp(-d1*t) + c*exp(-d2*t). If d1 differs from d2, its derivative has at most one zero: after division by a strictly positive exponential, a nonzero constant plus one strictly monotone exponential has at most one zero. Zero coefficients only reduce the number of possible zeros. Rolle's theorem therefore gives at most two zeros of the difference, unless the functions coincide. If d1=d2, the difference is constant plus one exponential and has at most one zero, unless it vanishes identically. The same two-zero conclusion holds when one member is a linear limit, P0+b*t, because its difference from an exponential has derivative constant plus one exponential.

## Four-point lower certificate

Let f_c be an admissible candidate. Suppose its residuals y_i-f_c(t_i) have alternating signs at four ordered indices, and every one of those four residuals has absolute value at least m>0. A different curve with uniform error strictly less than m would have to be strictly above f_c at each positive residual and strictly below f_c at each negative residual. Its difference from f_c would have at least three zeros by continuity, which is impossible by the previous paragraph. Thus epsilon_E >= m. If M bounds the absolute residuals of f_c at every observed point, epsilon_E <= M. Therefore m <= epsilon_E <= M. The proof does not require trusting the search grid, optimizer termination, or a local Hessian. It does not give a confidence region for the parameters or prove the chemical model true.

## Rational enclosure of the exponential

For rational z>=0 let S_N=sum_{k=0}^N z^k/k!, T=z^(N+1)/(N+1)!, q=z/(N+2). If q<1, the omitted positive terms have successive ratios at most q. Therefore S_N <= exp(z) <= S_N+T/(1-q). Taking positive reciprocals gives 1/(S_N+T/(1-q)) <= exp(-z) <= 1/S_N. verify_certificates.py evaluates these expressions with exact rational arithmetic until the exp(-z) interval width is <= 10^-42. It propagates the interval through the candidate curve and rounds all printed bounds outwards. Decimal data entries are treated as exact rationals for this technical calculation only.

## Monotone floor and possible improvement

For any nondecreasing sequence u_i and any i<j with y_i>y_j, max(abs(y_i-u_i),abs(y_j-u_j)) >= (y_i-y_j)/2. Let e_M=max_{i<j}(y_i-y_j)_+/2. The construction u_i=max(0,max_{k<=i}y_k-e_M) is nonnegative, nondecreasing and has uniform error at most e_M for the present nonnegative data. Thus e_M is exactly the minimum uniform error among all nonnegative nondecreasing sampled sequences. This is the unweighted linear-order case of known L-infinity isotonic regression.

Any additional monotone model family, such as a physically constrained Weibull saturation family, has error at least e_M. Hence its absolute improvement over the best exponential is at most M-e_M. When epsilon_E>0, its relative improvement is at most 1-e_M/M. These are upper bounds on a possible benefit, not estimates of the achieved improvement of a fitted Weibull model. The result concerns the uniform norm and does not compare sums of squares or out-of-sample errors.

## Dependence on uncertain nominal values

For fixed time points and a fixed model family, distance to that family in the uniform norm is 1-Lipschitz in the values. If max_i abs(y_i-y'_i)<=delta, then max(0,m-delta)<=epsilon_E(y')<=M+delta. This follows by the triangle inequality. It does not cover time-axis errors, arbitrary correlated calibration transformations, or an unknown experimental error distribution. Consequently a positive nominal lower bound alone is not experimental model invalidation.

## Prior work

Classical alternation: C. B. Dunham (1979), DOI 10.1016/0771-050X(79)90027-5, and G. A. Watson (1990), DOI 10.1093/imanum/10.4.569. The publisher abstracts establish the general prior-art concepts; the full 1979 and 1990 papers were not accessible in this pass. The family-specific proof above is given in full rather than assuming identical admissibility conditions.

Bounded-error methods: L. Jaulin and E. Walter (1993), DOI 10.1016/0005-1098(93)90106-4; J. Garloff, I. Idriss and A. P. Smith (2007), DOI 10.1007/s11155-007-9034-9. The 1993/2007 publisher abstracts were checked; an author-hosted Jaulin/Walter chapter was also read. No equivalence of every model constraint is asserted.

Isotonic regression: Q. F. Stout (2018), DOI 10.1016/j.jcss.2017.09.001; author-hosted published text, especially sections 2–3.
