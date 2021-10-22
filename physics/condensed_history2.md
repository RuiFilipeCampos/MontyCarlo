# Class II Condensed History

:::caution

This stuff is still being written.

:::


Consider an isotropic medium of atomic number density $N$.
Without loss of generality, any point in this medium can be chosen as the origin for the trajectory of some charged particle. 
Let $s$ be the displacement of the particle without interacting with the medium, $f(s; \mathbf{r}, \hat \mathbf{d})$ be the conditional probability density function of the random variable $(\mathbf{r}, \hat \mathbf{d})$ given that a displacement $s$ has occurred, where $\mathbf{r}  \in \mathbf{R}^3$ locates a point in space and $\hat \mathbf{d} \in S$ is the direction of the particle. The diffusion equation for this distribution can be written as

$$
    (\partial_s + \partial_{\hat d} ) f(s; \mathbf{r}, \hat d)
    =
    N \sigma_{el} \int_{S^2} 
    \left [  
    f(s; \mathbf{r}, \hat d) - f(s; \mathbf{r}, \hat d') 
    \right ] 
    p(\theta, \phi)
    d\Omega'
$$

where $p(\theta, \phi)$ is the probability density function of a polar angle deflection of $\theta$ with respect with the original direction of movement. 

Consider the harmonical expansion in $\hat d$,
 

$$
    f(s; \mathbf{r}, \hat d) = \sum_{l=0}^{\infty} \sum_{m = -l}^l f_l^m(s; \mathbf{r}) Y_l^m(\hat d)
$$


inserting it into equation ():

$$
    Y_l^m(\hat d) \partial_s f_l^m(s; \mathbf{r})
    +
    f_l^m(s; \mathbf{r}) \partial_{\hat d} Y_l^m(\hat d)
    = 
    N \sigma_{el} 
    f_l^m(s; \mathbf{r})
    \int_{S^2}
    \left [  
    Y_l^m(\hat d) -  Y_l^m(\hat d') 
    \right ] 
    p(\theta, \phi)
    d\Omega'
$$



$$
    Y_\lambda^\mu(\hat d)^* Y_l^m(\hat d) \partial_s f_l^m(s; \mathbf{r})
    +
    \nabla_{\mathbf{r}} f_l^m(s; \mathbf{r}) \cdot Y_\lambda^\mu(\hat d)^* \hat d Y_l^m(\hat d)
    \\
   =  \\
    N \sigma_{el} 
    f_l^m(s; \mathbf{r})
    \int_{S^2}
    \left [  
    Y_\lambda^\mu(\hat d)^* Y_l^m(\hat d) -  Y_\lambda^\mu(\hat d)^* Y_l^m(\hat d') 
    \right ] 
    p(\theta, \phi)
    d\Omega'
$$

$$
    \delta_{\mu,m}\delta_{\lambda,l} \partial_s f_l^m(s; \mathbf{r})
    +
    f_l^m(s; \mathbf{r}) \mathbf{Q}_{\lambda \mu}^{lm}
   =  
    N \sigma_{el} 
    f_l^m(s; \mathbf{r})
    \int_{(S^2)^2}
    \left [  
    Y_\lambda^\mu(\hat d)^* Y_l^m(\hat d) -  Y_\lambda^\mu(\hat d)^* Y_l^m(\hat d') 
    \right ] 
    p(\theta, \phi)
    d\Omega' d\Omega
$$

Consider the expansion of $p(\theta)$ in legendre polynomials


$$
    p(\theta) = \sum_n a_n P_n(\cos \theta)
$$

and the addition theorem for spherical harmonics

$$    
P_l(\cos \theta) = \sum_{l} Y()T()
$$


$$
    f_l^m(s; \mathbf{r})
    \int_{(S^2)^2}
    \left [  
    Y_\lambda^\mu(\hat d)^* Y_l^m(\hat d) -  Y_\lambda^\mu(\hat d)^* Y_l^m(\hat d') 
    \right ] 
    P_n(\cos \theta)
    d\Omega' d\Omega
$$

$$
    f_l^m(s; \mathbf{r})
    \left [ 
    \int_{(S^2)^2}
    Y_\lambda^\mu(\hat d)^* Y_l^m(\hat d) -  Y_\lambda^\mu(\hat d)^* Y_l^m(\hat d') 
    P_n(\cos \theta)
    d\Omega' d\Omega
    \right ] 

$$

---

$$
    f_l^m(s; \mathbf{r})
    \int_{(S^2)^2}
    \left [  
    Y_\lambda^\mu(\hat d)^* Y_l^m(\hat d) -  Y_\lambda^\mu(\hat d)^* Y_l^m(\hat d') 
    \right ] 
    Y_a^b(\hat d)^*Y_a^b(\hat d') 
    d\Omega' d\Omega
    =
$$

$$
    f_l^m(s; \mathbf{r})
    \int_{(S^2)^2}
    \left [  
    Y_\lambda^\mu(\hat d)^* Y_l^m(\hat d)Y_a^b(\hat d)^*Y_a^b(\hat d')  -  Y_\lambda^\mu(\hat d)^* Y_l^m(\hat d')Y_a^b(\hat d)^*Y_a^b(\hat d')  
    \right ] 
    d\Omega' d\Omega
    =
$$

$$
    f_l^m(s; \mathbf{r})
    \int_{(S^2)^2}
    \left [  
    Y_\lambda^\mu(\hat d)^* Y_l^m(\hat d)Y_a^b(\hat d)^*Y_a^b(\hat d')  -  Y_\lambda^\mu(\hat d)^* Y_l^m(\hat d')Y_a^b(\hat d)^*Y_a^b(\hat d')  
    \right ] 
    d\Omega' d\Omega
    =
$$