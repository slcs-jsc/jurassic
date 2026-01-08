# Retrieval theory

This page describes the theoretical background of the inverse modelling
capabilities implemented in JURASSIC. The retrieval framework follows
the principles of **optimal estimation** and provides a consistent
approach to deriving atmospheric state variables from infrared radiance
measurements.

The implementation closely follows the formalism described by Rodgers
(2000) and is fully integrated with the JURASSIC forward radiative
transfer model.

---

## Forward and inverse problems

In atmospheric remote sensing, the **forward problem** consists of
computing radiances from a given atmospheric state using a radiative
transfer model. In JURASSIC, this is performed by the forward model
described in the [Theory](theory.md) section.

The **inverse problem** aims to estimate the atmospheric state from a
set of measured radiances. This problem is generally ill-posed because:

- the measurements are finite and noisy,
- the atmospheric state is high-dimensional,
- different atmospheric states can produce similar radiances.

Optimal estimation provides a statistically well-defined framework to
address these challenges.

---

## State vector and measurement vector

### State vector

The atmospheric state to be retrieved is represented by a **state
vector** \(\mathbf{x}\). Depending on the application, the state
vector may include:

- temperature at selected altitude levels,
- volume mixing ratios of one or more trace gases,
- additional parameters such as background continua or scaling factors.

The mapping between the physical atmospheric representation used by
JURASSIC and the state vector is handled internally by state-mapping
operators.

### Measurement vector

The **measurement vector** \(\mathbf{y}\) contains the observed
radiances (or brightness temperatures) for the selected spectral
channels and observation geometries.

Measurement noise and forward-model uncertainties are represented by a
measurement error covariance matrix.

---

## Forward model operator

The forward model is represented by an operator \(\mathbf{F}\) such
that:

\[
\mathbf{y} = \mathbf{F}(\mathbf{x})
\]

In JURASSIC, \(\mathbf{F}\) corresponds to the radiative transfer
model, including ray tracing, spectral approximations, and emissivity
lookup tables.

For retrieval applications, the forward model must be differentiable
with respect to the state vector. JURASSIC therefore provides analytic
or semi-analytic **Jacobians** (also referred to as weighting functions
or kernels).

---

## Jacobians and sensitivity matrices

The **Jacobian matrix** \(\mathbf{K}\) describes the sensitivity of
the measurements to changes in the state vector:

\[
K_{ij} = \frac{\partial y_i}{\partial x_j}
\]

In JURASSIC, Jacobians are computed consistently with the forward model
using the same numerical discretization and spectral approximations.
This consistency is essential for stable and physically meaningful
retrievals.

Jacobians can be computed for individual atmospheric quantities, such
as temperature or trace gas concentrations, and are reused within the
iterative retrieval algorithm.

---

## Optimal estimation framework

Optimal estimation formulates the retrieval as the minimization of a
cost function that combines measurement information and prior
knowledge:

\[
J(\mathbf{x}) =
(\mathbf{y} - \mathbf{F}(\mathbf{x}))^T
\mathbf{S}_\epsilon^{-1}
(\mathbf{y} - \mathbf{F}(\mathbf{x}))
+
(\mathbf{x} - \mathbf{x}_a)^T
\mathbf{S}_a^{-1}
(\mathbf{x} - \mathbf{x}_a)
\]

where:

- \(\mathbf{x}_a\) is the **a priori** state vector,
- \(\mathbf{S}_a\) is the a priori covariance matrix,
- \(\mathbf{S}_\epsilon\) is the measurement error covariance matrix.

The solution represents a compromise between fitting the measurements
and remaining consistent with prior information.

---

## Iterative solution and linearization

Because the forward model is generally non-linear, the retrieval
problem is solved iteratively. At each iteration, the forward model is
linearized around the current state estimate:

\[
\mathbf{F}(\mathbf{x}) \approx
\mathbf{F}(\mathbf{x}_i) +
\mathbf{K}_i (\mathbf{x} - \mathbf{x}_i)
\]

JURASSIC employs a Levenberg–Marquardt-type scheme to ensure numerical
stability and convergence. The damping parameter balances between
Gauss–Newton and gradient-descent behaviour.

Iterations continue until convergence criteria based on state updates
and cost-function changes are satisfied.

---

## Averaging kernels and information content

The **averaging kernel matrix** \(\mathbf{A}\) describes the
sensitivity of the retrieved state to the true atmospheric state:

\[
\mathbf{A} = \frac{\partial \hat{\mathbf{x}}}{\partial \mathbf{x}}
\]

Averaging kernels provide insight into:

- vertical resolution of the retrieval,
- sensitivity to true atmospheric variations,
- influence of the a priori information.

JURASSIC computes averaging kernels and related diagnostic quantities
that can be used to assess retrieval performance and information
content.

---

## Error characterization

The retrieval error covariance matrix \(\mathbf{S}_x\) quantifies
the expected uncertainty of the retrieved state:

\[
\mathbf{S}_x =
(\mathbf{K}^T \mathbf{S}_\epsilon^{-1} \mathbf{K}
+ \mathbf{S}_a^{-1})^{-1}
\]

This matrix accounts for both measurement noise and the imposed prior
constraints. JURASSIC provides access to retrieval error estimates that
can be used for scientific interpretation and validation.

---

## Practical considerations

When performing retrievals with JURASSIC, users should consider:

- appropriate choice of a priori state and covariance,
- realistic measurement error estimates,
- consistency between forward-model configuration and retrieval setup,
- potential correlations between retrieved parameters.

Poorly chosen prior information or inconsistent configurations can
lead to biased or unstable retrievals.

---

## Summary

JURASSIC implements a physically consistent and well-established
optimal estimation retrieval framework that is tightly integrated with
its forward radiative transfer model. The use of analytic Jacobians and
shared numerical infrastructure ensures stable and efficient retrieval
calculations.

For a detailed theoretical treatment of optimal estimation, users are
referred to Rodgers (2000) and the references listed in the
[References](references.md) section.
