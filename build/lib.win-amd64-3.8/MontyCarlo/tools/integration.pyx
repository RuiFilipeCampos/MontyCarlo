

#roots of the 20th legendre polynomials
cdef tuple roots = (-0.9931285991850949,
-0.9639719272779137,
-0.912234428251326,
-0.8391169718222189,
-0.7463319064601508,
-0.6360536807265149,
-0.510867001950827,
-0.37370608871541955,
-0.22778585114164504,
-0.0765265211334973,
0.0765265211334973,
0.22778585114164504,
0.37370608871541955,
0.510867001950827,
0.6360536807265149,
0.7463319064601508,
0.8391169718222189,
0.912234428251326,
0.9639719272779137,
0.9931285991850949)

#corresponding weights
cdef tuple w = (0.017614007139152694 ,
0.04060142980038749 ,
0.06267204833410933 ,
0.08327674157670428 ,
0.10193011981724028 ,
0.11819453196151843 ,
0.1316886384491765 ,
0.1420961093183818 ,
0.1491729864726036 ,
0.15275338713072562 ,
0.15275338713072562 ,
0.1491729864726036 ,
0.1420961093183818 ,
0.1316886384491765 ,
0.11819453196151843 ,
0.10193011981724028 ,
0.08327674157670428 ,
0.06267204833410933 ,
0.04060142980038749 ,
0.017614007139152694 )


cdef class Integrator:

	@staticmethod #public Integrator 
	def create(object f, long double x0, long double xf):
		self = <Integrator>Integrator.__new__(Integrator)

		self.tol = 1e-16
		self.f = f

		self.unstable_intervals = (Interval.create(f, x0, xf),)
		self.stable_intervals = ()
		self.result = self.sumResult()

		cdef tuple temp
		cdef Interval interval
		cdef long double previous_result
		cdef int k = 0

		while self.unstable_intervals:
			temp = self.unstable_intervals
			self.unstable_intervals = ()

			for interval in temp:
				I1, I2 = interval.split()

				if abs(I1.result + I2.result - interval.result) < self.tol:
					self.stable_intervals += (interval,)

				else: self.unstable_intervals += (interval,)

			previous_result = self.result

			self.result = self.sumResult()

			if abs(self.result - previous_result) < 1e-16 and k>2:
				break
			k += 1
			if k > 10:
				print("###########################################")
				print("I_(n-1) = ", previous_result)
				print("I_n = ", self.result)
				print("error =", self.error)
				raise ValueError("Integral is not converging!")

		self.error = abs(self.result - previous_result)
		return self
		


	cdef long double sumResult(self):
		cdef Interval interval
		cdef long double result = 0.

		if self.unstable_intervals:

			for interval in self.unstable_intervals:
				result += interval.result
		if self.stable_intervals:

			for interval in self.stable_intervals:
				result += interval.result

		if not self.unstable_intervals and not self.stable_intervals:
			raise RuntimeError("No intervals?")

		return result




	cdef splitAll(self):
		cdef tuple temp = ()
		cdef Interval interval, I1, I2
		for interval in self.unstable_intervals:
			I1, I2 = interval.split()
			temp += (I1, I2)
		self.unstable_intervals = temp

cdef class Interval(Integrator):

	@staticmethod
	cdef public Interval create(object f, long double x0, long double xf):
		self = <Interval>Interval.__new__(Interval)
		self.x0, self.xf = x0, xf
		self.f = f
		self.result = self.integrate(x0, xf)
		self.stable = False
		return self

	cdef public tuple split(self):
		cdef long double mid_point = (self.x0 + self.xf)*.5

		cdef Interval I1 = Interval.create(self.f, self.x0,   mid_point)
		cdef Interval I2 = Interval.create(self.f, mid_point, self.xf)

		return I1, I2


	cdef long double integrate(self, long double a, long double b):
			#normalizing the interval [a, b]
		cdef long double z, m, c
		m = .5*(b - a)
		c = .5*(b + a)

		cdef long double I = 0.

		for i in range(20):
			z = m*roots[i] + c
			I += self.f(z)*w[i]
		return I*(b-a)*.5

	def __repr__(self):
		return f"[{self.x0}, {self.xf}]"

	def __str__(self):
		return f"[{self.x0}, {self.xf}]"


