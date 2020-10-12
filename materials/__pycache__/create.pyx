

ctypedef (double)(*func)(double)


cdef create(object material):
	#CONSTRUCTING COHERENT########################################################
cdef struct CoherentType:
	RationalInterpolation FF
	func				  thomsonDCS

ctypedef CoherentType coherent


#the actual thomson DCS
cdef double f(double cos):
	return (1 + cos)*.5

coherent.thomsonDCS = &f #giving its adress to the coherent struct

#note that FF is a callable, RITA should be able to instantiate
coherent.FF = RationalInterpolation.create(material.photon.coherent.FF, 
										   0, 300**2, True)
########################################################


#CONSTRUCTING PHOTON 
#############################################################
cdef struct PhotonType:
	CoherentType coherent

ctypedef PhotonType photon
photon.coherent = coherent
################################################################



#####################################################33
cdef struct MaterialType:
	PhotonType photon

ctypedef MaterialType material
material.photon = photon
#####################################################








cdef struct photon
	coherent

cdef struct material
	photon

material
	photon
		coherent
			FF
			thomsonDCS
		incoherent

