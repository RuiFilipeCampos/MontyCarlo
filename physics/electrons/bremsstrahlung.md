# Bremsstrahlung

The differential cross section for the **energy loss** of the electron upon a Bremsstrahlung event is frequently expressed in the form

$$
    \left (
        \frac{d\sigma}{dW}
    \right )^{\textrm{(brem)}}
        
         =
    \frac{Z^2}{\beta^2}
    \frac{1}{W}
    \chi(E_0, k)
$$

where $Z$ is the atomic number of the atom, $\beta = v/c$, $W$ the electrons energy loss, $E_0$ its initial energy, $k = W/E$ the reduced energy and $\chi$ is known as the scaled bremsstrahlung DCS.

The use of numerical tables to sample from $\chi$ is ubiqutous. MontyCarlo uses the tables prepared at NRC for the EGSnrc code system. They were calculated using the electron-nucleous bremsstrahlung (relevant at high $Z$) provided and recommended by NIST and the improved electron-electron bremsstrahlung tables as calculated by \cite{Frédéri-2008}. The random sampling procedure is done via the analytic inverse transform on $1/W$ (section \ref{rand:IT}) with a rejection step using linearly interpolated values of $\chi$ (section \ref{sec:rejection_sampling}).

Following \cite{bielajew1989improved}, after sampling the energy loss $k$, the differential cross section as defined in equation 2BS of \cite{KOCH-1959} is used to sample the direction of emission of the photon 


$$
    d\sigma^{(\textrm{brem})}  = 
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
    \right \}   \\
$$

where $\theta$ the polar deflection angle of the emitted photon, $ y = E_0 \theta$ , $E = E_0 - W$ is the energy of the outgoing electron, $W$ is the photons energy loss.

This equation can be neatly packed into the following form


$$
    d\sigma^{(\textrm{brem})}  \propto f_{E_0}(y^2) N_r g_r(y^2) d(y^2) 
$$

where


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
