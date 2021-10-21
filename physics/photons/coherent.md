# Coherent Scattering

## Thomson Scattering

### Setting up the problem

Consider an plane electromagnetic wave incident on a charged particle that is free at rest on the position $\mathbf{r}$. The most general form of the electric field component of the plane EM wave can be written in phasor notation as 

$$
\tilde \mathbf{E} = \hat \mathbf{ s} E_0 \exp \left ( 
i (\mathbf{k} \cdot \mathbf{r} - wt + \phi_0)
\right )
$$

where $\hat \mathbf{ s}$ is the polarization of the field, which is linear and can be chosen at will without loss of generality, $E_0$ is the amplitude of the E-field, $\mathbf{k}$ is the wave vector (indicates the direction of propagation) and $w$ is the frequency of oscillation. The actual value of the electric field can be extracted by taking the real part of $\tilde \mathbf{E}$, that is,

$$
\mathbf{E} = \textrm{Re}(
    \tilde \mathbf{E}
).
$$

Since the field is always pointing in the s-direction, which can be chosen to be any direction without loss of generality, the charge will be set in motion in that direction. Meaning that this is a 1d problem and only the s-coordinate needs to be considered. Using Newtons Law,

$$
m  \mathbf{a} \cdot \hat \mathbf{s} = -e \textrm{Re}(
    \tilde \mathbf{E} \cdot \hat \mathbf{s}
).
$$

This equation can be rewritten as

$$
\ 

\textrm{Re}\left ( \frac{d^2 s}{dt^2} \right ) 
=
\textrm{Re}\left( -\frac{e}{m} 
    \tilde \mathbf{E} \cdot \hat \mathbf{s}
\right).
$$

Laying all out in phasor notation,

$$
\frac{d^2 \tilde s}{dt^2}
=
- \left [ 
\frac{e}{m} 
E_0 
\exp  ( 
        i \mathbf{k} \cdot \mathbf{r} 
      + \phi_0) \right ]
\exp  ( 
    - iwt
 )
$$

Which can be easily solved by repeated integration:

$$
\frac{d \tilde s}{dt}
=
- \left [ 
\frac{e}{-iwm} 
E_0 
\exp  ( 
        i (\mathbf{k} \cdot \mathbf{r}  + \phi_0)
      ) \right ]  \exp  ( 
    - iwt
 ) + C_1
$$

and

$$
\tilde s =
- \left [ 
\frac{e}{w^2m} 
E_0 
\exp  ( 
        i (\mathbf{k} \cdot \mathbf{r}  + \phi_0)
      ) \right ] \exp  ( 
    - iwt
 ) + C_1t + C_2

$$

Applying the boundary conditions of the problem ($v_0 = 0$ and $s_0 = 0$) results in 

$$
C_1
=

- \frac{e}{iwm} 
E_0 
\exp  ( 
        i (\mathbf{k} \cdot \mathbf{r}  + \phi_0)
      ) 
$$

and

$$
C_2 =
 \frac{e}{w^2m} 
E_0 
\exp  ( 
        i (\mathbf{k} \cdot \mathbf{r}  + \phi_0)
      ) 

$$



The movement of the electron is then described by 

$$
\tilde s =
- \left [ 
\frac{e}{w^2m} 
E_0 
\exp  ( 
        i (\mathbf{k} \cdot \mathbf{r}  + \phi_0)
      ) \right ] \exp  ( 
    - iwt
 ) 
 \\ - \left [ \frac{e}{iwm} 
E_0 
\exp  ( 
        i (\mathbf{k} \cdot \mathbf{r}  + \phi_0)
      ) \right ] t 
  \\ + \frac{e}{w^2m} 
E_0 
\exp  ( 
        i (\mathbf{k} \cdot \mathbf{r}  + \phi_0)
      ) 

$$

which can be further manipulated to a form where taking the real part is trivial,

$$
\tilde s = \frac{E_0e}{wm}
\left \{
\frac{

\exp  ( 
        i (\mathbf{k} \cdot \mathbf{r}  + \phi_0)
      )
      -
\exp  ( 
        i (\mathbf{k} \cdot \mathbf{r} - wt + \phi_0)
      ) 
      
 }{w} 
 \\ + t 

\exp  ( 
        i (\mathbf{k} \cdot \mathbf{r}  + \phi_0 + \pi/2)
      )  


\right \}
$$

where taking the real part is then

$$
 s = \frac{E_0e}{wm}
\left \{
\frac{

\cos  (\mathbf{k} \cdot \mathbf{r}  + \phi_0)
-
\cos (
    \mathbf{k} \cdot \mathbf{r} - wt + \phi_0
)

      }{w} 
 \\ + t 

\sin   (\mathbf{k} \cdot \mathbf{r}  + \phi_0) 


\right \}
$$



- [ ] use that to get to dP/dOmega
- [ ] motivate the way in which dP/dOmega is used to get to the DCS

### The (differential) re-radiated power



> (notE to self) need to either re-derive this eq., or cite it

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


### The incident power


### The differential cross section

## Rayleigh Scattering

The DCS:

$$
        \frac{d\sigma}{d\Omega} 
        =
        r_e^2 \frac{1 + \cos^2 \theta}{2} F_Z(q)^2,
$$