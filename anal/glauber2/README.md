# Glauber Monte Carlo ROOT Macro

This is a simple ROOT macro for Monte Carlo Glauber collisions of two identical
hard-sphere nuclei. By default it is configured close to Pb+Pb at LHC energies:

- `A = 208`
- hard-sphere radius `R = 6.62 fm`
- inelastic nucleon-nucleon cross section `sigmaNN = 70 mb`
- negative-binomial particle production per participant nucleon
- mean `20.0` particles per participant
- negative-binomial width parameter `k = 1.0`
- `useNegativeBinomial = true`

Run with ROOT:

```sh
root -l -q 'glauber_mc.C(10000)'
```

or choose parameters explicitly:

```sh
root -l -q 'glauber_mc.C(50000, 208, 6.62, 70.0, 20.0, 20.0, 1.0, true, "pbpb_nbd.root", 12345)'
```

For constant particles per participant, set `useNegativeBinomial` to `false`:

```sh
root -l -q 'glauber_mc.C(50000, 208, 6.62, 70.0, 20.0, 20.0, 1.0, false, "pbpb_constant.root", 12345)'
```

The macro writes:

- `events` tree with `event`, `b`, `npart`, `ncoll`, `nparticles`, and `multiplicityModel`
- `hB` impact-parameter histogram
- `hNpart` participant histogram
- `hNcoll` binary-collision histogram
- `hNparticles` produced-particle histogram
- `hNpartNcoll` correlation histogram
- `hNpartNparticles` correlation histogram
- `hBNparticles` produced particles versus impact parameter

The particle histograms currently cover `Nparticles = 0..11999`, which is large
enough for the default `20 * Npart` model in Pb+Pb.

The collision condition is the usual black-disk approximation:

```text
d_T < sqrt(sigmaNN / pi)
```

with `sigmaNN` converted from mb to `fm^2` using `1 mb = 0.1 fm^2`.

Nucleons are sampled uniformly in volume inside the sphere:

```text
r = R * u^(1/3)
cos(theta) uniform in [-1, 1]
phi uniform in [0, 2pi]
```

The particle-production model can be selected with `useNegativeBinomial`.

With `useNegativeBinomial = true`, particles are negative-binomial per participant:

```text
for each participant:
  n_i ~ NBD(meanParticlesPerParticipant, nbdK)

Nparticles = sum_i n_i
```

The macro samples the negative binomial using the Gamma-Poisson construction:

```text
lambda ~ Gamma(k, mean / k)
n_i ~ Poisson(lambda)
```

With `useNegativeBinomial = false`, particles are constant per participant:

```text
Nparticles = round(meanParticlesPerParticipant * Npart)
```

In the tree, `multiplicityModel = 1` means NBD and `multiplicityModel = 0` means constant.
