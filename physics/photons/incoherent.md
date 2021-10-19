# Compton Scattering

:::caution

This stuff is still being written.

:::


Traditionally, most Monte Carlo codes have used the Klein-Nishina(KN) differential cross section for sampling the incoherent interaction,

$$
    \frac{d\sigma^{(\textrm{KN})}}{d\Omega} 
    =
    \pi r_0^2 Z \left  (
                        \frac{w_C}{w_0}
                \right )^2 
    \left ( 
    \frac{w_C}{w_0} + \frac{w_0}{w_C} - \sin^2 \theta
    \right ),
$$

where
$$
    w_C = \frac{w_0}{1 + (w_0/m_e c^2) (1 - \cos \theta)}
$$

is known as the compton line, $w_0$ is the incident photon energy and $\theta$ its polar deflection angle after the interaction, $r_0$ the classical electron radius, and $Z$ the atomic number of the atom.  This DCS is derived by assuming that the electron is free and at rest \cite{book:quantum}, which ends up being a good a approximation at high energies, but neglects a series of low energy considerations that results in inconsistencies.

The neglected effects include \cite{P.M-1997}:


- **Binding**: Since the electron is bound, the photon must provide the necessary energy for ionization. Otherwise the interaction cannot occur.
- **Doppler Broadening**: The bound electrons are not only moving, their momentum is undetermined. Thus causing a broadening of the energy distribution for a fixed polar angle.
- **Compton Defect**: The peak of the energy distribution is shifted in relation to its value for the free at rest case.
- **Compton-Ramman Resonance**: 
- **Infrared Divergence**:

Currently, monte carlo codes have proceeded to implement random sampling algorithms based on doubly differential cross sections (DDCS) that account for some of these effects \cite{R.Mar-2004, Longo-2008, Y-1994, Ye-2006}.
These models are based on the impulse approximation(IA) in conjunction with the independent particle approximation(IPA):


- **IPA**: Electron-Electron interactions are assumed to be negligible, resulting in the fact that the electron that will interact with the photon can be equiprobably chosen.
- **IA**: Binding energy of the electrons is small. It is assumed that the electron is scattered into plane waves.

The implementation in Monty Carlo draws from both the current Geant4 and EGSnrc implementations. 
 
The DDCS for the collision between an unpolarized photon and a bound electron in the i'th shell can be given by


$$
    \frac{d^2\sigma^{(i)}}{d\Omega dw} =
        \frac{r_0^2}{2} 
        \left(\frac{w_C w}{w_0^2} \right) 
        \frac{dp_z}{dw}
        \left ( 
        \frac{w_C}{w_0} + \frac{w_0}{w_C} - \sin^2 \theta
        \right )
        J_i(p_z),
$$


where $w_0$ is the incident photon energy, $w$ the outgoing photon energy, $\theta$ the polar deflection angle, $p_z$ the z projection of the electrons momentum given by

$$
    p_z = -137 \frac{w_0 - w - w_0w(1 -\cos \theta)/(m_ec^2)}{|\textbf{k} - \textbf{k}_0 |},
$$


$\mathbf{k}$ and $\mathbf{k}_0$ the outgoing and incident photon momentum vectors,


$$
    \frac{dp_z}{dw} = \frac{137w_0}{| \textbf{k} - \textbf{k}_0 | w_C} - \frac{p_z(w-w_0 \cos \theta)}{ |\textbf{k} - \textbf{k}_0 |},
$$

$J_i(p_z)$ is the probability distribution of the bound electrons momentum, known as the compton profile for the i'th shell, given by

$$
    J_i(p_z) = \int_{\textbf{R}^2} dp_x dp_y  | \psi_i(\mathbf{p}) | ^ 2,
$$

and $\psi_i$ is the electron wavefunction projected onto momentum space. 
The compton profiles have been derived and tabulated by Biggs \cite{biggs1975hartree} in the non-relativistic impulse approximation for atoms with $Z < 37$ and in the relativistic impulse approximation for atoms with $Z > 35$. 
The data itself can be found in Geant4's repositories \cite{site:geant4data}.
The DCS in solid angle is obtained by integrating (\ref{compton:bound_DDCS}) with respect to the outgoing photon energy while assuming that $k = k_C$

$$
    \frac{d\sigma^{(i)}}{d\Omega} =
        \frac{r_0^2}{2} 
        \left(\frac{k_C^2}{k_0^2} \right) 
        \left ( 
        \frac{k_C}{k_0} + \frac{k_0}{k_C} - \sin^2 \theta
        \right )
        S_i(k_0, \theta)
$$


where $S_i(k_0, \theta)$ is the incoherent scattering function for the i'th shell as given by

$$
    S_i(p_z) = \int_{-\infty}^{p_{z, max}^{(i)}}  J_i(u)du.
$$

The upper bound of this integration is obtained from (\ref{compton:p_z}) by substituting $w$ by its maximum value of $k_0 - U_i$, where $U_i$ is the binding energy of the ith shell. 
Finally, the DCS for the whole atom is obtained by summing the contribution from each shell

$$
    \frac{d\sigma}{d\Omega} =
        \frac{r_0^2}{2} 
        \left(\frac{w_C^2}{w_0^2} \right) 
        \left ( 
        \frac{w_C}{w_0} + \frac{w_0}{w_C} - \sin^2 \theta
        \right )
        S(k_0, \theta)
$$

where $S(k_0, \theta)$ is the incoherent scattering function for the atom and is given by

$$
S(k_0, \theta) = \sum_i S_i(k_0, \theta).
$$

Note how in equation (\ref{comp:atom_dcs}), the KN DCS has been recovered. 
Then, to sample this interaction, the IPA is used  by randomly selecting an active atom based with discrete probabilities based on the normalized stoichiometric coefficients of the molecules formula. 
This sampling is done with the Wilkersons Aliasing. 

Then, the polar deflection from the KN DCS is sampled. For this, the well known sampling algorithmn developed by (namito again?) was used to generate a proposal, which is then rejected against $S(k_0, \theta)$.

Monty Carlo uses the incoherent scattering function tabulated in the EPDL. 
The tabulation of the incoherent scattering function present in the EPDL is made in function of

$$
    x = \frac{\sin(\theta/2)}{\lambda}  \textrm{    [cm]}^{-1}
$$

Where $\lambda$ and $\theta$ are the wavelenght and polar angle deflection of the \textit{scattered} photon. In order to save computational time, the following transformations were made on the array containing the $x$ values.

First, note that since $\theta$ is the polar angle, we have that $\theta \in [0, \pi]$. Thus, we can write $\sin (\theta/2)$ in terms of $\cos(\theta)$ without fear of losing information about a possible negative sign.

$$
\sin \frac{\theta}{2} = \sqrt{\frac{1  - \cos \theta}{2}} 
$$

Secondly, we can write the wavelenght in terms of the photons energy, and put it in a form that is readily available from the sampling algorithm,

$$
    \frac{k_0}{m_e c^2} = \frac{h}{m_e c} \frac{1}{\lambda} = k'_0
$$

Putting equation (\ref{eqn:sinthetaovertwo}) and (\ref{eqn:k}) into (\ref{eqn:S_argument}) we get,

$$
    \frac{h}{m_e c} x \cdot 10^{-2} =  \sqrt{\frac{1  - \cos \theta}{2}} k_0
$$

where $x$ is multiplied by $10^{-2}$ to change its units to SI. We can further manipulate to obtain an even more suitable form for the sampling algorithm,

$$
2 \cdot 10^{-4} \left ( \frac{h}{m_e c} x \right ) ^ 2  
    =  
    (1  - \cos \theta) k_0  ^ 2
$$

Both $1 - \cos \theta$ and $ \varepsilon k_C$ are direct outputs of the algorithm. 
The transformation of $x$ indicated in the LHS of (\ref{eqn:TRANSFORMx}) is made at initialization time. 

When $\theta$ has been sampled, a subshell is randomly chosen based on occupancy numbers. 
The energy of the scattered photon is sampled from (\ref{compton:dcs_bound_electron}). 
The upper bound for the integration in equation (\ref{eqn:S_i(p_z)}) is calculated and $p_z$ is sampled using 

$$
p_z = C^{-1}(\xi C(p_{z, max})) 
$$

Where

$$
C(x) = \int_0^x J_i(u)du 
$$

and $\xi$ is a random number uniformly distributed in $[0, 1[$. Both the cumulutative and inverse cumulutative tables have been interpolated using the method described in section \ref{sec:interp}, which does greatly increase the ammount of numerical data present in memory during simulation, but note that this interaction is already quite slow. 
A proposal for the new energy value is calculated by introducing the sampled $p_z$ value into equation (\ref{compton:p_z}) and solving for $k$. 
A rejection step, $k/k_C<1$, is performed to account for the second term in equation (\ref{compton:bound_DDCS}).


