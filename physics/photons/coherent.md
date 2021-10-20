# Coherent Scattering

## Thomson Scattering

Consider an incident electromagnetic plane wave on a free charged particle that is initially at rest. The EM wave can be written in versor notation

$$
\tilde \mathbf{E} = \hat \mathbf{ s} E_0 \exp \left ( 
i (\mathbf{k} \cdot \mathbf{r} - wt)
\right )
$$

where $\hat \mathbf{ s}$ is the polarization of the field (which is linear), $E_0$ is the amplitude of the field, $\mathbf{k}$ is the wave vector (indicates the direction of propagation) and $w$ is the frequency of oscillationç. The actual value of the electric field can be extracted by taking the real part of $\tilde \mathbf{E}$, that is,

$$
\mathbf{E} = \textrm{Re}(
    \tilde \mathbf{E}
).
$$

Since the field is always pointing in the s-direction, which can be chosen to be any direction without loss of generality, the charge will be set in motion in that direction. Meaning that this is a 1d problem and only the s-coordinate needs to be considered. Using Newtons Law,

$$
m  \mathbf{a} \cdot \hat \mathbf{s} = e \textrm{Re}(
    \tilde \mathbf{E} \cdot \hat \mathbf{s}
).
$$

The radiated power per unit solid angle of an electric dipole (i.e, the oscillating charge), is given by 

> (not to self) need to either re-derive this eq., or cite it

$$
\frac{dP}{d\Omega}
=
\frac{
    e^2 \left < 
       \left (  \mathbf{a} \cdot \hat \mathbf{s} \right ) ^2
    \right >
}{

    16 \pi ^2 \epsilon_0 c^3
}

\sin ^2 \theta       

$$

Substituting the previous equation here results in

$$

\frac{dP}{d\Omega}
=
\frac{
    e^2 \left < 
          \left (  
                e \textrm{Re}(
                    \tilde \mathbf{E} \cdot \hat \mathbf{s} ) / m         
           \right ) ^2
    \right >
}{

    16 \pi ^2 \epsilon_0 c^3
}

\sin ^2 \theta       

$$

which can be simplified to

$$

\frac{dP}{d\Omega}
=
\frac{
    e^4 \left < 
          
                 \textrm{Re}(
                    \tilde \mathbf{E} \cdot \hat \mathbf{s}  )          
            ^2
    \right >
}{

    16 \pi ^2 m^2 \epsilon_0 c^3
}

\sin ^2 \theta       

$$

substituting in the first equation

$$

\frac{dP}{d\Omega}
=
\frac{
    e^4 \left < 
          
                 \textrm{Re}(
                    E_0 \exp \left ( 
i (\mathbf{k} \cdot \mathbf{r} - wt)
\right ) )          
            ^2
    \right >
}{

    16 \pi ^2 m^2 \epsilon_0 c^3
}

\sin ^2 \theta       

$$

which can be better seen as


$$

\frac{dP}{d\Omega}
=

\textrm{Re} \left (


    
\frac{e^4 E_0^2 \sin ^2 \theta   }{    16 \pi ^2 m^2 \epsilon_0 c^3}

          \left <
                 
                     \exp \left ( 
i (\mathbf{k} \cdot \mathbf{r} - wt)
\right )           
            ^2
    \right >


     

\right )
$$

finally resulting in

$$

\frac{dP}{d\Omega}
=

\textrm{Re} \left (


    
\frac{e^4 E_0^2 \sin ^2 \theta   }{    16 \pi ^2 m^2 \epsilon_0 c^3}

          \frac{1}{2}


     

\right )
$$

that can be trivially simplified to


$$

\frac{dP}{d\Omega}
=




    
\frac{e^4 E_0^2 \sin ^2 \theta   }{    32 \pi ^2 m^2 \epsilon_0 c^3}

      


     
$$

## Rayleigh Scattering

The DCS:

$$
        \frac{d\sigma}{d\Omega} 
        =
        r_e^2 \frac{1 + \cos^2 \theta}{2} F_Z(q)^2,
$$