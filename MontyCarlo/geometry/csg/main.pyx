





# the kind of thing that will mess up multiprocessing/multithreading
cdef double displacement


cdef class CSGproxy(BVH):
	cdef intIterator iterator
	
	cdef void set_iterator(self, intIterator iterator):
		self.iterator = iterator
	

	cdef void _set_safest_distance(self, double3& pos):
		# note: assumes pos inside current volume
		self.distance = self.iterator.current() - displacement

	cdef void set_safest_distance(self, double3& pos):
		self.distance = self.iterator.current() - displacement

	# signals that it already is a proxy, no need for intersecting
	cdef bint main_intersect(self, double3& origin, double3& dire):
		return False



cdef class CSGvolume(BVH):
	cdef Proxy proxy

	def __init__(self, *args, **kwargs):
		# Opening lock, volume can be modified
		super(CSGvol, self).__init__(*args, **kwargs)


	cdef bint move(self, STATE& state, double SP):
		cdef double3 origin = state.pos
		displacement = 0

		cdef Closest first
		cdef Closest second
		cdef int i

		IF DEBUG_MODE:
			input(string(state) + "Entering move method.")
			input(string(state) + "How does the workspace look like?")
			input("\t" + self.print_ws(state)) 
			input(string(state) + "Starting event loop:")

		
		
		while True: 

			# gets the safest KNOWN distance 
			self._set_safest_distance(state.pos)


			# find the closest surface
			first.index = 0
			first.distance = self.distance

			for i in range(1, self.Nws):
				(<V> self.ws[i]).set_safest_distance(state)
				if (<V> self.ws[i]).distance < first.distance:
					first.distance = (<V> self.ws[i]).distance
					first.index = i

			# WHERE WE GONAN GO NOW


			# safe to just advance the particle 
			if state.L < first.distance:
				state.pos.x += state.dire.x*state.L
				state.pos.y += state.dire.y*state.L
				state.pos.z += state.dire.z*state.L

				# state.L = 0
				return False


			if first.distance < .1:

				if (<V> self.ws[first.index]).main_intersect(origin, state.dire):
					self.ws[first.index] = (<V> self.ws[first.index]).proxy

				first.distance  = (<V> self.ws[first.index])._get_safest_distance()
				second.distance = INF

				for i in range(0, first.index):
					IF DEBUG_MODE: print(i, (<V> self.ws[i]).distance, (<V> self.ws[i]))
					if (<V> self.ws[i]).distance < second.distance:
						second.distance = (<V> self.ws[i]).distance
						second.index = i

				for i in range(first.index+1, self.Nws):
					IF DEBUG_MODE: print(i, (<V> self.ws[i]).distance, (<V> self.ws[i]))
					if (<V> self.ws[i]).distance < second.distance:
						second.distance = (<V> self.ws[i]).distance
						second.index = i


				if first.distance == INF:
					if state.L < second.distance: # < first.distance

					
						state.pos.x += state.dire.x*state.L
						state.pos.y += state.dire.y*state.L
						state.pos.z += state.dire.z*state.L
						return False
					
					# second.distance < state.L < first.distance   
					# self.virtual_event(state, second.distance)

					state.pos.x += state.dire.x*second.distance
					state.pos.y += state.dire.y*second.distance
					state.pos.z += state.dire.z*second.distance
					state.L -= second.distance
					displacement += second.distance
					continue


				if first.distance < second.distance:
					if state.L < first.distance: 
						state.pos.x += state.dire.x*state.L
						state.pos.y += state.dire.y*state.L
						state.pos.z += state.dire.z*state.L
						return False

					#  self.virtual_event(state, first.distance)
					
					state.pos.x += state.dire.x*first.distance
					state.pos.y += state.dire.y*first.distance
					state.pos.z += state.dire.z*first.distance
					state.L -= first.distance
					displacement += first.distance

					(<V> self.ws[first.index]).iterator.inc()


					if first.index == 0:

						for i in range(0, self.index_in_outer):
							if ((<V> self.outer).ws[i]).is_inside(state.pos):
								state.current_region = self.ws[i]

								# staying in outer, must keep cached intersections
								if state.current_region == <void*> self.outer:
									self.keep = True
									self.exitINNER_TO_OUTER()
									return True
							
								# entering some adjacent volume, must intersect it then
								(<V> state.current_region).main_intersect(state)
								(<V> state.current_region).keep = True
								self.exitINNER_TO_INNER()
								return True

						for i in range(self.index_in_outer + 1, self.Nws):
							if ((<V> self.outer).ws[i]).is_inside(state.pos):
								state.current_region = self.ws[i]

								# staying in outer, must keep cached intersections
								if state.current_region == <void*> self.outer:
									self.keep = True
									self.exitINNER_TO_OUTER()
									return True

								# entering some adjacent volume, must intersect it then
								(<V> state.current_region).main_intersect(state)
								(<V> state.current_region).keep = True
								self.exitINNER_TO_INNER()
								return True


					# from outer to inner
					state.current_region = self.ws[first.index]
					self.exitOUTER_TO_INNER()
					return True

				# min() == L
				if state.L < second_nearest:
					IF VERBOSE: print("min() == L 222")
					self.final(state)
					self.exit()
					return False ? 

			self.virtual_event(state, first.distance)


	cdef void _set_safest_distance(self, double3& pos):
		# note: assumes pos inside current volume
		self.distance = -self.SDF(pos)

	cdef void set_safest_distance(self, double3& pos):
		self.distance = self.SDF(pos)

	cdef bint main_intersect(self, double3& origin, double3& dire):
		
		self.proxy.set_iterator(intIterator(self.intersect(
			pos, state.dire
		)))
		
		return True











	cdef str print_ws(self, STATE& state):
		cdef int i
		to_print = "["
		for i in range(self.Nws):
			to_print += f"(is_inside = {(<V> self.ws[i]).is_inside(state.pos)}, keep = {(<V> self.ws[i]).keep}, cache = {(<V> self.ws[i]).cache}) ,"
		
		to_print += "]"
		return to_print




	cdef void set_safest_distance(self, double3& pos):
		pass

	cdef bint is_inside(self, double3& pos):
		raise RuntimeError("`is_inside` was called from virtual (in CSGvol)")

	cdef void depositUNIFORM(self, STATE& state, double SP):
		pass

	cdef void depositLOCAL(self, double3& pos, double E):
		pass


	cdef void depositRANDOM(self, STATE& state, double E, double tau):
		pass


	cdef void depositLocaly(self, double3& pos, double E):
		pass

	#@lock("Modifiying volume after being closed")
	def rotate(self, axis, angle):
		return NotImplemented

	#@lock("Modifiying volume after being closed")
	def translate(self, direction, displacement):
		return NotImplemented



	cdef intLIST intersect(self, double3& pos, double3& dire):
		raise RuntimeError(".intersect called from virtual")


	cdef double SDF(self, double3 pos):
		raise RuntimeError("`.SDF` called from virtual (in CSGvol)")
















