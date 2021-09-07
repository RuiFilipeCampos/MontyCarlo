---
sidebar_position: 3
---

# Simulation Schemes

import useBaseUrl from '@docusaurus/useBaseUrl';


## Detailed History Simulation

In this section, an attempt is made to build an **intuition** for what is the *differential cross section* is and how it can be used in particle simulations. A bottom up approach is taken to walk the reader to the most important algorithm in the monte carlo techniques for the simulation of particle transport.

### Interaction Cross Section


The **differential cross section** (DCS) is the quantity that bridges the analytic solutions of scattering theory to the probability realm. 
For a given collision, it is a function of the final result of the interaction and can be seen as the unnormalized probability distribution of those values.

Consider the classical physics experiment shown in figure 1.
Any given point on the plane that contains all possible initial conditions (plane of incidence $\Sigma$) is injectively mapped to a tuple $(\theta, \phi) \in [0, \pi]\times[0, 2\pi]$ that specifies the final result of the interaction. That is, there is a one-to-one map between $\Sigma$ and the unit sphere ($S$).

<figure>
<img alt="Docusaurus with Keytar" src={useBaseUrl('/img/DCS.jpg')} />
<figcaption>
<b>Figure 1:</b> A mono-energetic beam of particles are thrown at a normal incidence to the plane $\Sigma$. Their trajectory is deflected by the scatterer and their assymptotic behaviour is described by the tuple $(\theta, \phi)$. Note how the physical details of the interaction between the projectile and the target are hidden. This is a blackbox representation of a colision.
</figcaption>
</figure>


Consider now a solid angle $d\Omega$ located at some point $(\theta_0, \phi_0)$ in the unit sphere. Since each point in this solid angle can be uniquely traced back to a point in the plane of incidence, there is a corresponding cross sectional area $d\sigma$ containing the initial conditions that will reproduce the set of final conditions represented by $d \Omega$ at $(\theta_0, \phi_0)$. By maintaining the amplitude of the solid angle constant, but varying its position on the unit sphere, the size of $d\sigma$ will change, thus quantifying the probability of a given final condition occurring.

Intuitively, if particles are randomly thrown at the plane of incidence, in the conditions described in figure \ref{fig:DCS}, $d\sigma$ is the size of the target that the particle must hit in order to reproduce one of the final conditions contained within $d \Omega$ at $(\theta_0, \phi_0)$. Certain final conditions correspond to larger targets, and therefore are more likely to occur. The differential cross section is then defined as the multiplication factor that scales $d\Omega$ placed at a position $(\theta, \phi)$, to the size of this target,

$$
  d\sigma = \frac{d\sigma_E}{d\Omega} |\_{(\theta, \phi)} d\Omega.
$$

The subscript $E$ indicates the dependence of this function on the energy of the incident beam. Furthermore, to account for the energy loss of the projectile, the unit sphere is extended to a more general configuration space $S \times [0, W[$. Finally, the the probability of the interaction occurring at all can be quantified by integrating the differential cross section over the entire configuration space,

$$
            \sigma(E)  = \int_0^W \int_{4\pi} \frac{\partial^2 \sigma_E}{\partial\Omega \partial W}d\Omega dW,
            %\sigma(E)=\int_0^EdT\int_0^\pi2\pi sin\theta d\theta \frac{d^2 \sigma(E,\theta T,)}{d\Omega dT}
$$
where $\sigma$ is known as the integrated cross section or the total cross section of the interaction. Note that, for a *classical* central potential collision, this integration will diverge to infinity. It reflects the fact that, no matter the initial condition, the trajectory of the particle will always be affected due to the infinite range of these potentials. That is, every $d\sigma \in \Sigma$ is being summed over. 

In these cases, the usual procedure is to establish a cut off in the integration limits. For the example in figure \ref{fig:DCS}, it can seen that a small $\delta \Omega$ located at $\theta = 0$ and $\phi = 0$ corresponds to an infinitely large set of initial conditions. A small solid angle can then be excluded from the domain of integration. The size of the excluded $\delta \Omega$ is context dependent, but the usual reasoning is that this amplitude corresponds to the resolution of the available measuring apparatus. In the context of monte carlo simulations, this ad-hoc procedure is acceptable as long as it is consistent throughout all the implemented interactions and it can be freely changed (that is, it can be made as small as possible).  


### Integrated Cross Section

As previously mentioned, the result of the integration shown in equation (\ref{eqn:def:DCS}) quantifies the probability of a given interaction occurring. Since $d\sigma$ is the size of the target that a particle must hit to reproduce some final condition, $\sigma$ is then the size of the target that the particle must hit in order to interact at all.

When several interactions are considered, their total cross sections can be added,

$$
            \sigma_{tot} = \sum_{i} \sigma_{i},
$$


this follows directly from the additive nature of equation (\ref{eqn:def:DCS}). The actual probabilities are obtained by normalizing the quantities being discussed. The probability of the $i^{\textrm{ith}}$ interaction occurring is given by $\sigma_i / \sigma_{tot}$, while the probability density function of the final state of the particle is given by

$$
            p_E^{(i)}(W, \theta, \phi)  = \frac{1}{\sigma_i} \frac{\partial^2 \sigma_E}{\partial\Omega \partial W}.
            %\sigma(E)=\int_0^EdT\int_0^\pi2\pi sin\theta d\theta \frac{d^2 \sigma(E,\theta T,)}{d\Omega dT}
$$

Random sampling of (\ref{eqn:PDF_DCS}) is done via the evaluation of analytic formulas upon the generation of a single uniform random number or through the use of numerical methods on tabulated values of the differential cross section.
 

### Mean Free Path

### The Algorithm


## Condensed History Simulation