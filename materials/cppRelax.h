#ifndef AAA_HEADER
#define AAA_HEADER


#include <random>
#include <vector>
#include <string>
#include <deque>
#include <stack>
#include <list>
#include <forward_list>
#include "../_random/mixmax/mixmax.hpp"


using namespace std;

namespace model{
	

	/*
	Represents a transition from the shell where it is stored to the shells it is pointing to.
	- Radiative Transitions.
	- Non Radiative Transitions.
	*/
	//struct Transition;
	
	
	
	
	struct Transition ;


	
	/*A memoryview analog in C++ -> array of Transitions.*/


	class Shell{
			public:
				/* INITIALIZING SHELL*/
				Shell ();
				
				bool LAST;
				
				Shell(double* frac, bool dontSIMULATE);
				void offset_trSTART(int size);
				
				double *original_frac;
				bool dontSIMULATE;
				double *frac; //this will be the most dangerous part of the program

				

				/* SETTING TRANSITIONS/WILKERSONS ALIASING*/
				void setTRANSITIONS(Transition* arr, int N);
				Transition *trSTART;
				int Nt;
				
	

				
				
				Transition* sample_transition();




				
				
				/*FOR DEBUGGING*/
				string list_transitions();

				
		};
		

			

	Transition constructRadiative(double E, Shell* j);
	Transition constructnonRadiative(double E, Shell* j, Shell* k);

	

	
	/*A memoryview analog in C++ -> array of shells.*/
	struct arraySHELL{
		Shell *start;
		Shell *end;
	};	
	
	struct PARTICLES{
		vector<double> PHOTONS; 
		vector<double> ELECTRONS;
	};
	
	
	void setADRESS_RAD(Transition * tr);
	void setADRESS_NONRAD(Transition * tr);


	
	
	
	class Atom{

		
		
		public:
			/* INITIALIZING ATOM*/
			
			
			double R;
			int i ;
			Shell* shSTART ;
			Transition *trSTART;
			Atom();

			Atom(Shell* arr, int Nsh, int Nt);
			arraySHELL SHELLS;
			int Nsh;            //neded for bound checking, probably will remove it
			

			Shell* fetchFI(int index); // helper method for setting shell transitions, probably unnecessary, could declare shell array in .pxd and use that
			
			void run(int shell_index,  PARTICLES *particles, mixmax_engine *genPTR);
			
			
			
			

			
			
			/*FOR DEBUGGING*/
			string list_array();

			
			

			
			


			};
			


			
			
			
};

#endif