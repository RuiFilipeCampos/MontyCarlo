cdef extern from "mixmax_release_200final/mixmax.hpp":
    ctypedef unsigned int myID_t;
    cdef cppclass mixmax_engine:
        mixmax_engine(myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID );	   #Constructor with four 32-bit seeds
        mixmax_engine(); # Constructor, no seeds
        double get_next_float();
