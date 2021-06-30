



from ..._init import eax
import numpy as np
from ..._init cimport EAX
from scipy.integrate import quad

from libc.math cimport sqrt, log, fmin, pi
cdef double  ELECTRON_REST_MASS      = 0.51099895000e6            
cdef double  _2ELECTRON_REST_MASS    = 2 *ELECTRON_REST_MASS
cdef double  _4ELECTRON_REST_MASS_2  = 4*ELECTRON_REST_MASS**2
cdef double  SPEED_OF_LIGHT          = 2.99792458e10
cdef double  H_BAR                   = 6.5821e-16 
cdef double  MASS_ELECTRON           = 9.1094e-28 
cdef double  ELECTRON_CHARGE         = 4.8032e-10

cdef double two_mec2 = _2ELECTRON_REST_MASS


cdef double CONST = 2*pi *ELECTRON_CHARGE**4 / MASS_ELECTRON  *6.242e+11






def do(formula, shells):
	"""
	Calculate contribution of soft inelastic events
	to sIMFP1 and sIMFP2. 

	formula -> dynamic_dict that holds pretty much everything
	shells  -> memview or array -> [[Wk, fk], [Wk, Fk], ...] 


	The DCS of the inelastic event is split into three parts:
		(1) Close interactions
		(2) Far longitudinal interaction
		(3) Far transverse interaction

	The contribution of (3) to the angular defflection of the
	projectile is negligeable. So it is not considered.

	The contribution of (1) has been calculated analyticaly
	using sympy. Those values are calculated numerically in this
	function.

	It was not possible to solve the contribution of (2) anlytically,
	so I decided, at least for now, use numerical integration.

	Each contribution is completly specified by:
		Oscillator Parameters
			Wk, fk
		Projectile Parameters
			E
		Cut off values
			Wcc
			Wmax -> maximum value of energy transfer kinematically allowed
	
	The contributions are integrated in mu (angular deflection on the projectile)
	in their corresponding domains

	(1) 0 < mu < mu_k -> 0 < W < Wk
	(2) mu_k < mu < mu_cc -> Wk < W < min(Wcc, Wmax)

	**The shells provided to this function are assumed to have
	resonance energies below Wcc.**


	References:
	See PENELOPE manual 2016 version chapter (??).
	"""

	formula.log.add_header("Calculating inelastic sIMFP1 and sIMFP2")
	print("CALCULATING INELASTIC sIMFP's...")



	cdef double CONST0 = formula.N*CONST

	shells = np.array(shells)
	cdef double Wcc = formula.Wcc

	FK = shells[:, 0]
	WK = shells[:, 1]
	cdef double[:] avgMU = np.zeros(len(eax))
	cdef double[:] avgMU2 = np.zeros(len(eax))

	cdef double[:] sIMFP1 = np.zeros(len(eax))
	cdef double[:] sIMFP2 = np.zeros(len(eax))
	cdef double X, beta2, v2


	cdef double E
	cdef int i, j
	cdef double upper
	cdef int Ni = len(eax)
	cdef int Nj = len(WK)

	cdef double cp, cpk, cp0sq, Rk, arg
	cdef double Dcpcpk
	
	for i in range(Ni):
		E = EAX[i]
		cp = sqrt(E*(E + two_mec2))
		Wmax = fmin(Wcc, E/2)
		for j in range(Nj):
			if WK[j] > E/2 or WK[j] == 0 or FK[j] == 0: continue 
			Wk = WK[j]


			 
			#print(E, Wk, FK[j])

			# CONTRIBUTION FROM CLOSE
			sc_avgMU  = quad(muDCSB,  Wk, Wmax, args=(eax[i],) )[0]
			sc_avgMU2 = quad(mu2DCSB, Wk, Wmax, args=(eax[i],) )[0]

			# CONTRIBUTION FROM FAR LONGITUDINAL (FAR TRANSVERSE IGNORED)
			cpk = sqrt((E - Wk)*(E - Wk + two_mec2))
			Dcpcpk = (cp - cpk)**2
			cp0sq = Wk*(Wk + two_mec2)

			sc_avgMU += FK[j]*(Dcpcpk*(1 - log(Dcpcpk/cp0sq)) + cp0sq)/Wk/cp/cpk/4


			sc_avgMU2 +=  FK[j]*( 2*Dcpcpk*Dcpcpk* (log(cp0sq/Dcpcpk) - cp0sq + Dcpcpk) + (cp0sq - Dcpcpk)**2 )/32/Wk/cp**2/cpk**2

			avgMU[i]  += FK[j]*sc_avgMU
			avgMU2[i] += FK[j]*sc_avgMU2

		X = (E/ELECTRON_REST_MASS + 1)**2
		beta2 = (X-1)/X
		v2 = SPEED_OF_LIGHT**2 * beta2

		# multiplying missing constant and it's up and ready! 
		sIMFP1[i] += CONST0/v2 * 2 * avgMU[i]
		sIMFP2[i] += CONST0/v2 * 6 * (avgMU[i] - avgMU2[i])

	# logging the graph for debug purposes

	fig = formula.log.new_plot()

	formula.log.add_to_plot(fig, EAX, sIMFP1)
	formula.log.add_to_plot(fig, EAX, sIMFP2)
	formula.log.finish_plot(fig)

	return sIMFP1, sIMFP2


	




cdef double mu2DCSB(double W, double E):
	return mu(W, E)**2 *DCSB(W, E)

cdef double muDCSB(double W, double E):
	return mu(W, E)*DCSB(W, E)

cdef inline double DCSB(double W, double E):
	cdef double a = (E/(E + two_mec2))
	return 1/(E-W)**2 - 1/W*1/(E-W) + a*(1/W*1/(E - W) + 1/(E*E)) + 1/W**2


cdef inline double mu(double W, double E):
	cdef double cpE = sqrt(E*(E + two_mec2))
	cdef double cpW = sqrt((E - W)*(E - W + two_mec2))
	return ( W*(W+ two_mec2) - (cpE - cpW)**2 )/4/cpW/cpE




