---
sidebar_position: 6
---

# Electrons


![Electrons](/img/electrons.png)


## Elastic Scattering

An [elastic scattering](https://en.wikipedia.org/wiki/Elastic_scattering) of an electron with an atom is defined as an interaction where the internal quantum states of the atom are left untouched. Meaning that no excitation/ionization occurs in the shell structure of the atom. When considering the elastic scattering, a reference frame is chosen where the target (atom) is at rest while the projectile (electron) is deflected by the electrostatic field of the atom. Since the rest mass of the atom is much larger than the mass of the electron, the recoil energy and momentum of the atom is negligible. However, the electron may be strongly deflected. This change in movement causes the emission of electromagnetic radiation, the so called Bremsstrahlung emission or Breaking Radiation.

In the theoretical treatment and consequently, in the implementation, this interaction is divided into two independent events:

- **Angular Deflection**: it accounts for the deflection of the projectile, implemented in the `._elastic` method (section ??);

- **Bremsstrahlung Emission**: it accounts for the photonic emission due to the acceleration(change of direction) of the electron, implemented in the `._brem` method (section ??) 

In the following subsections, these two interactions are given a brief theoretical overview. The angular deflection is entirely based off of the PENELOPE code system while the Bremsstrahlung is taken from EGSnrc.

### Angular Deflection

The electric potential of an atom at rest is formally given by

$$
    V(\textbf{r}) = z_0 e \psi(\textbf{r}) +  V_{\textrm{ex}}(\textbf{r})
$$

where $\psi(\textbf{r})$ is the electrostatic potential of the atom, $z_0$ is the charge sign of the projectile (electron or positron), $e$ the elementary charge and $V_{\textrm{ex}}(\textbf{r})$ the ... The electrostatic potential is obtained by a model choice of charge density of the nucleus ($\rho_n(\textbf{r})$) and that of the electron cloud($\rho_e(\textbf{r})$). When the electron cloud is not symmetrical, it is averaged over all directions


$$
    \rho_e(r) = \frac{1}{4\pi}\int_{4\pi} \rho_e(\textbf{r})d\Omega 
$$

to obtain a spherical symmetrical charge density. The proton density is chosen to be the fermi distribution

$$
    \rho_n(r) = \frac{\rho_0}{\exp((r - R_n)(4\ln 3/t)) + 1} 
$$

The total charge density is given by $\rho(r) = \rho_n(r) + \rho_e(r)$ and the electrostatic potential of the atom by

$$
    \psi(r) = \frac{Q(r)}{r} + \int_r^\infty e\rho(r') 4\pi r' dr' 
$$

where 

$$
    Q(r) = 4\pi \int_0^r \rho(r')r'^2 dr'    
$$

is the total charge within a radius of $r$. 
Since the potential is symmetric under any azymuthal transformation, the result of the interaction will also have azymuthal symmetry, that is $\phi = 2 \pi \xi$. The differential cross section is given by 

$$
    \frac{d\sigma^{(\textrm{el})}}{d\Omega} = |f(\theta)|^2 + |g(\theta)|^2
$$

where

$$
    f(\theta) = \frac{1}{2ik} \sum_{l=0}^\infty
   \left \{
        (l + 1) [
                \exp(2i\delta_{l+}) - 1
                ]
                + 
                l [
                \exp(2i\delta_{l-}) - 1
                ] 
    \right \} 
    P_l(\cos \theta),
$$

is the direct scattering amplitude and

$$
    g(\theta) = \frac{1}{2ik} \sum_{l=0}^\infty 
    \{\exp(2i\delta_{l-}) - \exp(2i\delta_{l+})\} P^1_{l}(\cos \theta)
$$

is the spin-flip scattering amplitude. $k$ is the electrons wave number, $P_l(x)$ the legendre polynomials, $P_l^1(x)$ the associated legendre polynomials and $\delta_x$ are phase shifts associated with the assymptotic behaviour of the Dirac radial functions.


Evidently, conjuring up an algorithm for sampling from (\ref{electron:elastic:DCS}) is not a trivial matter. Numerical tables for this DCS was generated with the ELSEPA program \cite{ELSEPA} for the energy range of 1keV to 100MeV and sampling is done via the numerical inverse transform (with hash lookup) as described in section \ref{sec:sampling:numerical_IT}. For the 100MeV to 1GeV range the modified Wentzel model was used.

The Wentzel model for the probability distribution in $\mu = (1-\cos \theta)/2$ is given by

$$
    p^{(\textrm{W})}_{A_0}(\mu) = \frac{A_0 (1+A_0) }{(\mu + A_0)^2}
$$


whose first and second moment are

$$
    \avg{\mu}_W = A_0 \left [ 
    (1 + A_0) \ln \left (
                    \frac{1+A_0}{A_0} - 1
                 \right )
    \right ]
$$

and


$$
    \avg{\mu^2}_W = A_0 \left [1 - 2 \avg{\mu}_W \right].
$$

The value of $A_0$ is chosen such that this distribution reproduces the first moment of the distribution implicitly defined in (\ref{eqn:DCS:elastic}). Whose value is given by

$$
    \langle \mu \rangle = \frac{1}{2} \frac{\sigma_{ \textrm{el,1}} }
                                            {\sigma_{ \textrm{el} }}
$$

Where $\sigma_{\textrm{el}}$ is the total cross section as defined by equation (\ref{}) and $\sigma_{\textrm{el}, 1}$ is the first transport cross section defined in section \ref{sec:CHII}. However, this choice will not necessarily yield the correct second moment: 


$$
    \langle \mu^2 \rangle = \frac{1}{2} \frac{\sigma_{ \textrm{el,1}} }
                                            {\sigma_{ \textrm{el} }}
                                            -
                            \frac{1}{6} \frac{\sigma_{ \textrm{el,2}} }
                                            {\sigma_{ \textrm{el} }}
$$

where $\sigma_{\textrm{el,2}}$ is the second transport cross section. To accommodate for this, Savat introduced a modified Wentzel model:

$$ \label{eqn:MWdcs}
    p^{(\textrm{MW})}(\mu) = \left\{\begin{matrix}
(1-B_1)p^{(\textrm{W})}_{A_1}(\mu) + B_1 \delta(\mu - \langle\mu\rangle) & \langle\mu^2\rangle_W < \langle \mu^2 \rangle \\ 
(1-B_2)p^{(\textrm{W})}_{A_2}(\mu) + B_2 8 (\mu - 1/2) \Theta(\mu-1/2) &  \langle\mu^2\rangle_W > \langle \mu^2 \rangle
\end{matrix}\right.
$$

All these parameters have been defined such that the whole distribution reproduces the correct first and second moments. By equating them, the following expressions are obtained:

\begin{align}
    A_1 = & A_0 \notag \\ 
    B_1 = & \frac{\avg{\mu^2}_{W,A_1} - \avg{\mu^2}}
               {\avg{\mu^2}_{W,A_1} - \avg{\mu}^2}
\end{align}

In the case of $A_2$ the following equation is obtained

$$
    \left ( 
    \frac{17}{24} - \avg{\mu^2}
    \right ) \avg{\mu}_{\textrm{W,}A_2}
    -
    \left ( 
    \frac{5}{6} - \avg{\mu}
    \right )\avg{\mu^2}_{\textrm{W,}A_2}
    =
    \frac{17}{24} \avg{\mu}
    -
    \frac{5}{6} \avg{\mu^2}
$$

which is solved numerically. The value of $B_2$ is then given by

$$
    B = \frac{\avg{\mu} - \avg{\mu}_{\textrm{W,}A_2}}{5/6 - \avg{\mu}_{\textrm{W,}A_2}}.
$$

Sampling from (\ref{eqn:MWdcs}) is done on the basis of the composition method (see section \ref{sec:sampling:composition}): the distribution to be sampled is chosen with probabilities $B_i$ and the chosen distribution is sampled. The W distribution and the trangle distribution are sampled via the inverse transform method (section \ref{rand:IT}). 

### Bremsstrahlung Emission

The differential cross section for the energy loss of the electron upon a Bremsstrahlung event is usually expressed in the form

$$
    \frac{d\sigma^{\textrm{(brem)}}}{dW} =
    \frac{Z^2}{\beta^2}
    \frac{1}{W}
    \chi(E, k)
$$

where $Z$ is the atomic number, $\beta = v/c$, $W$ the electrons energy loss, $E$ its initial energy, $k = W/E$ and $\chi$ is known as the scaled bremsstrahlung DCS. The use of numerical tables to sample from $\chi$ is ubiqutous. MontyCarlo uses the tables prepared at NRC for the EGSnrc code system. They were calculated using the electron-nucleous bremsstrahlung (relevant at high $Z$) provided and recommended by NIST and the improved electron-electron bremsstrahlung tables as calculated by \cite{Frédéri-2008}. The random sampling procedure is done via the analytic inverse transform on $1/W$ (section \ref{rand:IT}) with a rejection step using linearly interpolated values of $\chi$ (section \ref{sec:rejection_sampling}).

Following \cite{bielajew1989improved}, after sampling the energy loss $k$, the differential cross section as defined in equation 2BS of \cite{KOCH-1959} is used to sample the direction of emission of the photon 


$$
%\begin{align}\label{eqn:DCS:angle_brem}
    d\sigma^{(\textrm{brem})} & = 
    \frac{4Z^2r_0^2}{137}
    \frac{dk}{k}
    ydy
    \left \{
        \frac{16y^2E}{(y^2 + 1)^4E_0}
        -
        \frac{(E+E_0)^2}{(y^2+1)E_0^2}
        +
        \left  [ 
            \frac{E_0^2 + E^2}{(y^2+1)^2E_0^2}
            -
            \frac{4y^2E}{(y^2+1)^4 E_0}
        \right ]
        \ln M(y)
    \right \} \notag  \\
    & \propto f_{E_0}(y^2) N_r g_r(y^2) d(y^2) 
%\end{align}
$$


where $E_0$ is the energy of the incident electron, $\theta$ the polar deflection angle of the emitted photon, $y = E_0 \theta $, $E = E_0 - W$ is the energy of the outgoing electron, $W$ is the photons energy loss,

$$
    \frac{1}{M(y)} = 
    \left ( 
    \frac{k}{2E_0E}
    \right )^2
    +
    \left (
    \frac{Z^{1/2}}{111(y^2+1)}
    \right )^2,
$$

$$
    f_{E_0}(x) = \frac{1 + 1/(\pi E_0)^2}{(x+1)^2} ,
$$

$$
    g_r(x) = 2r - 3(1+r^2)
    - 
    \left (
    4 + \ln m(x)
    \right )
    \left [ 
    (1 + r^2)
    -
    \frac{4xr}{(x + 1)^2}
    \right ]
$$

$r = E/E_0$ and

$$
    m(x) = 
    \left ( 
    \frac{1-r}{2E_0r}
    \right )^2
    +
    \left (
    \frac{Z^{1/3}}{111(x+1)}
    \right )^2
$$

The proportionality in equation  (\ref{eqn:DCS:angle_brem}) hides the constant leading factors, which include both $k$ and $dk$. This form was chosen such that $f_{E_0}$ is a valid distribution in $x = y^2$, which can be sampled using the analytical inverse transform (section \ref{rand:IT}). while $N_r g_r(x)$ is then used for a rejection step(section \ref{sec:rejection_sampling}). For maximum efficiency, $N_r$ is chosen such that the maximum of $g_r$ is, at least approximately, normalized:

$$
    \frac{1}{N_r} = \max(g(0), g(1), g(\pi^2E_0^2))
$$




## Inelastic Scattering


### Kinematics

Collision between a charged particle and the medium, which is at rest:



\newcommand{\po}{\mathbf{p}_0}
\newcommand{\pf}{\mathbf{p}_f}
\newcommand{\Eo}{\mathbf{E}_0}
\newcommand{\Ef}{\mathbf{E}_f}
\newcommand{\q}{\mathbf{q}}


- \textbf{Projectile}
    - $\po$: linear momentum before collision 
    - $\pf$: linear momentum after collision
    - $\Eo$: energy before collision
    - $\Ef$: energy after collision
    - $W$: energy loss

The momentum transfer in a collision is

$$
    \q = \pf - \po
$$

$$
    Q = \q^2/(2m_ec^2)
$$


### DCS

 What follows is a short theoretical treatment of the interaction of an high energy electron with the medium it is travelling in, leading up to the construction of a simplified model, appropriate for random sampling during a monte carlo simulation.

Reference \cite{} will be followed very closely and the model itself was constructed for the PENELOPE code system. It is, however, an outdated model which does not account for post ionization relaxation effects. In the older version of PENELOPE, this was accounted for by introducing an artificial impact ionization event. In MontyCarlo, an ad-hoc modification was made so that the model is able to account for this without introducing an artificial event.


The differential cross section for the interaction of the electron with the bound electrons of an atom is given by \cite{} (cite fano)

$$ \label{eqn:inel:DCS:at}
    \frac{d \sigma^{\textrm{(inel)}}}{d W dQ} = 
    \frac{2\pi z_0^2e^4}{m_ev^2}     % leading constants
    \frac{df(Q, W)}{dW}              % GOS
    \left \{ 
        \frac{2m_ec^2}{WQ(Q + 2m_ec^2)}
        + 
        \frac{\beta^2\sin^2\theta_rW2m_ec^2}{[Q(Q+2m_ec^2)-W^2]^2}
    \right \}
$$




which is corrected to include the E-Field screening 

$$ \label{eqn:inel:DCS:cm}
    \frac{d \sigma_k^{\textrm{(inel)}}}{dW dQ} = 
    \frac{2\pi z_0^2e^4}{m_ev^2}     % leading constants
    \underbrace{
    \frac{df(Q, W)}{dW}              % GOS
     }_{\textrm{GOS}}
    \left \{ 
        \underbrace{
        \frac{2m_ec^2}{WQ(Q + 2m_ec^2)}
        }_{\textrm{Transverse Interaction}}
        +
        \underbrace{
        \left (
            \frac{\beta^2\sin^2\theta_rW2m_ec^2}{[Q(Q+2m_ec^2)-W^2]^2}
            -
            D(Q, W)
        \right ) 
        }_{\textrm{Longitudinal Interaction}}
    \right \}
$$

where $\frac{df(Q, W)}{dW}$ is known as the generalized oscillator strength(GOS) and $D(Q, W)$ as the density correction. Equations (\ref{eqn:inel:DCS:at}) and (\ref{eqn:inel:DCS:cm}) are very much the same apart from the $D$ subtraction.

### The Generalised Oscillator Strength (GOS)

The generalised oscillator strenght for an $mn$ electronic transition was first defined by \cite{Bethe-1997} as 

$$
    f_{mn}(q) = \frac{E_n - E_m}{Rh}\psi_{mn}(q)
$$

where $R$ is the Rydburg energy,

$$
    \psi_{mn}(q) = \frac{a^2}{q^3} |\varepsilon_{mn} (q)|^2
$$

is the known as the generalized transition probability from $m$ to $n$ and $|\varepsilon_{mn} (q)|^2$ gives the conditional probability of the $mn$ transition given that a magnitude of momentum transfer of $q$ will occur. The generalized oscillator strength (for discrete excitations) is then defined as $f_m(q) = f_{0m}(q)$.

In \cite{INOKUTI-1971} this quantitiy is generalized to include, not only the discrete excitations, but the transitions to continuum (i.e. ionization), which is expressed as "number of oscillators" per "unit energy loss":

$$ \label{GOS_def}
    \frac{df(W, Q)}{dW} = \sum_n \frac{W}{Q} | \varepsilon_n (q(Q))|^2 \delta(W-E_n)
$$

where $n$ runs over transitions from ground state to discrete and the continuum. Calculations of this quantity are only known for two systems: the free-electron gas and the hydrogen atom. The GOS for inelastic colisions with the hydrogenic K-Shell is shown in figure \ref{}.



### The Liljequist GOS model

The Liljequist GOS model for the $k$ shell is defined as

$$ %\label{GOSmodel}
     F_k(Q, W) = f_k \left\{\begin{matrix}
 \delta(W-W_k) & Q < W_k  & \textrm{(distant/resonant collision)} \\ 
 \delta(W-Q) &  Q>W_k     & \textrm{(close/binary collision)}
\end{matrix}\right.
$$

where $W_k$ is taken as the energy loss on a resonant interaction with the $k$-shell, known as resonance energy and $f_k$ is the oscillator strength of the $k$-shell, identified as the number of electrons in the shell. When the the energy of the momentum transfer is larger than the resonance energy, the collision is taken to be a binary collision between two free electrons, wherein the energy loss is the same as the recoil energy. 



By summing over all the shells in the atom, the GOS is modelled:

$$
     \frac{df(Q,W)}{dW} = \sum_k F_k(Q, W).
$$

However the $(f_k, W_k)$ parameters are defined, they must obey the following rules:

$$
\begin{align} 
& \sum_k f_k = Z \label{eqn:bethe_sum} \\
& \sum_k f_k \ln W_k = Z \ln I  \label{eqn:ionization}
\end{align}
$$

where $I$ is the mean ionization energy of the medium. Equation (\ref{eqn:bethe_sum}) is the Bethe sum rule and equation (\ref{eqn:ionization}) is a requirement that comes from ?? and implies that $W_k$ are the partial mean ionization energies. 

The conduction band of the medium is first modelled as a free electron gas, whose plasmon energy is given by 

$$
    W = \hbar w_p = \sqrt{\frac{4\pi e^2 \hbar^2}{m_e}} n_0^{1/2}
$$

where $n_0$ is the electron number density of the gas. Secondly, in the Liljequist GOS model, the conduction band is modelled as its own shell. Then, the resonance energy in equation (\ref{}) is identified with the plasmon energy in equation (\ref{}), resulting in the following definition

$$
\begin{align}
    W_{cb} & = \sqrt{\frac{4\pi N e^2 \hbar^2}{m_e}} f_{cb}^{1/2} \\
           & = W_m \left (\frac{f_{cb}}{Z} \right )^{1/2}
\end{align}
$$

where $W_m$ is the plasmon energy of an electron gas with the same density as the medium under consideration. The oscilator strenght of this conduction band model is taken to be the number of electrons (per atom) that are present in the conduction band. The resonance energy of the bound electrons are then determined by

$$
    W_k = \sqrt{(aU_i)^2 + \frac{2W_m^2}{3}  \frac{f_{k}}{Z}  }
$$

where $U_i$ are the ionization energies and $a$ is known as the correction parameter (???). Note that $a$ is the only unknown parameter and equation (\ref{}) is automatically satisified by identifying the $f_k$ as the number of electrons in each shell. Equation (\ref{}) is used to numerically determine the value of $a$, thus completely defining the model.


Inserting equation \ref{} into \ref{} results in 

$$
    \frac{d \sigma_k^{\textrm{(inel)}}}{d W d Q} =
     \frac{2\pi z_0^2e^4 f_k }{m_ev^2} % leading constants
    \left\{\begin{matrix}
         \delta(W-W_k)\frac{2m_ec^2}{WQ(Q + 2m_ec^2)} +   \delta(W-W_k)      \left (
            \frac{\beta^2\sin^2\theta_rW2m_ec^2}{[Q(Q+2m_ec^2)-W^2]^2}
            -
            D(Q, W)
        \right )   & Q < W_k   \\ 
        \delta(W-Q)     \left \{ 
        \frac{2m_ec^2}{WQ(Q + 2m_ec^2)}
        +
        \left (
            \frac{\beta^2\sin^2\theta_rW2m_ec^2}{[Q(Q+2m_ec^2)-W^2]^2}
            -
            D(Q, W)
        \right ) 
    \right \}  &  Q>W_k
    \end{matrix}\right.
$$

which defines the three different interaction regimes of the electron projectile with the $k$-shell:

- Distant Transverse
- Distant Longitudinal
- Close 


\newcommand{\farT}{^{\textrm{(far.t)}}}
\newcommand{\farL}{^{\textrm{(far.l)}}}
\newcommand{\close}{^{\textrm{(close)}}}


\subsubsection{Distant Transverse Interaction}

The DCS for the inelastic interaction in the distant-transverse regime is given by (cf. eq. \ref{})

$$
    \frac{d \sigma_k^{\textrm{(inel)}}}{d W d Q} = \frac{2\pi z_0^2e^4 f_k }{m_ev^2}        \left (
            \frac{\beta^2\sin^2\theta_rW2m_ec^2}{[Q(Q+2m_ec^2)-W^2]^2}
            -
            D(Q, W)
        \right )  \delta(W-W_k)
$$

for $W < W_k$. Following \cite{}, the expression in brackets is substituted by a simpler expression

$$
    \frac{d \sigma_k^{\textrm{(inel)}}}{d W d Q} = \frac{2\pi z_0^2e^4 f_k }{m_ev^2}
    \frac{1}{W}
    \left (
        \ln \left ( 
            \frac{1}{1 - \beta^2}
        \right )
        -
        \beta^2
        -
        \delta_F
        \right )  \delta(W-W_k) \delta(Q-Q_-)
$$

where

$$
    Q_- = 
$$

is the minimum kinematically allowed recoil energy and $\delta_F$ is known as the Fermi density effect correction on the stopping power defined as 

$$
    \delta_F = 
    \left \langle \ln \left ( 
    1 + \frac{L^2}{W^2}
    \right ) \right \rangle_{\textrm{OOS}}
    -
    \frac{L^2}{W_m^2}(1-\beta^2)
$$

where the average performed using the normalized optical oscillator strength ($df/dW|_{Q=0}$) and $L$ can be obtained by numerically solving

$$
    1-\beta^2 = \frac{W_m^2}{Z} \sum_k  \frac{f_k}{W_k^2 + L^2}
$$

In this interaction, the angular deflection impinged on the projectile is negligeble. When chosen, the particle's direction is not changed. Equation (\ref{}) is used to calculate the total cross section of the event, the stopping power and the straggling parameter. Discussion of how these calculations were done is being postponed to section \ref{}, since they are performed within within the framework of the class II condensed history scheme.

\subsubsection{Distant Longitudinal Interaction}

The DCS is: 

$$
    \frac{d \sigma_k^{\textrm{(inel)}}}{d W d Q} =
    \delta(W-W_k)\frac{2m_ec^2}{WQ(Q + 2m_ec^2)}.
$$

In simulations of the inelastic interaction, when it is known that a distant interaction with a particular $k$-shell is to occur, it implies that the energy loss of the projectile is known and given by  $W = W_k$. As such, the distribution implicitly defined by the DCS in (\ref{}) will be integrated over $W$ to obtain the partial probability distribution. The result of which is a simple enough distribution whose cumulative function can be inverted. The recoil energy is then sampled with the analytical inverse transform. 


\subsubsection{Close Interaction}

The DCS is:

$$
    \frac{d \sigma_k^{\textrm{(inel.close)}}}{d W d Q} = 
    \frac{2\pi z_0^2e^4 f_k }{m_ev^2}  \delta(W-Q)     \left \{ 
        \frac{2m_ec^2}{WQ(Q + 2m_ec^2)}
        +
        \left (
            \frac{\beta^2\sin^2\theta_rW2m_ec^2}{[Q(Q+2m_ec^2)-W^2]^2}
            -
            D(Q, W)
        \right ) 
    \right \} 
$$

Note that $W=Q$, meaning that 

\subsection{Adaptation of Electrons to Positrons}