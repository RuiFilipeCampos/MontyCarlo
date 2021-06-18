#include <vector>
#include <string>
#include <iostream>
#include "cppRelax.h"
#include <time.h>
#include <stdio.h>
#include <assert.h>
#include <stdexcept>
#include <deque> 
#include<time.h>
#include <fstream> 
#include <stack>
#include <list>
#include <forward_list>
#include <stdlib.h> 
#include "../_random/mixmax/mixmax.hpp"

// will eventually clean this up 

using namespace std;



namespace model{
	
	
	//srand(1212321);
	struct Transition {
		double p;
		Transition* t;
		double E;
		Shell* j;
		Shell* k;
		void (*perform)(Transition *, PARTICLES *);
	};
	
	





		
	

	Shell::Shell () {};
	
	
	Shell::Shell(double* frac, bool dontSIMULATE){
		this->frac = frac;
		this->original_frac = frac;
		this->dontSIMULATE = dontSIMULATE;
		this->LAST = false;
	};
	
	void Shell::setTRANSITIONS(Transition* arr, int N){
		this->trSTART = arr;
		this->Nt = N;
	};
	


	Transition* Shell::sample_transition(){
		/*
		double R = ((rand() % 100)*1e-2)*Nt;
		int i = (int) R;
		Transition *tr = TRANSITIONS.start + i;
		
		// urand() sometimes will be one, which will cause an out of bounds access here
		if (R - i < tr->p){ return tr; };
		return tr->t;*/
		return trSTART;
	};

	
	string Shell::list_transitions(){/*
		string repr = to_string(*frac) + "\n";
		for (Transition *ptr = TRANSITIONS.start; ptr < TRANSITIONS.end; ++ptr){
			repr += "        <Transition@" + to_string((int)ptr) +  " j = " + to_string((int)ptr->j) + ", ";
			repr += "k = " + to_string((int)ptr->k)  + ", ";
			repr += "p = " + to_string(ptr->p)  + ", ";
			repr += "t = " + to_string((int)ptr->t)  + + "> \n";
			
			
			
		};*/
		return "";
		
	};
	


	


	

	

	
	

			
	
	/**------------------------------ATOM--------------------------------------  **/
	/**------------------------------ATOM--------------------------------------  **/
	/**------------------------------ATOM--------------------------------------  **/

	Atom::Atom(){};
	
	
	
	Atom::Atom(Shell *arr, int Nsh, int Nt){
		

		this->shSTART = arr;
		

		//(this->SHELLS).start = arr;
		//(this->SHELLS).end = arr + N;
		//this->Nsh = N;
		
	
	};
	

	
	Shell* Atom::fetchFI(int index){
		//if (index >= Nsh || index < 0){
		//	throw out_of_range(">>> cppRelax.cpp > Atom::fetchFI > OUT OF BOUNDS - index = " + to_string(index));
		//};
		
		return (this->shSTART + index);
	};
	
	string Atom::list_array(){/*
		string a = "";
		int i = 0;
		for (Shell *ptr = SHELLS.start; ptr < SHELLS.end; ++ptr){
			a += "    <Shell #" + to_string(i) + "  ";
			a += "Nel = " + to_string(ptr->Nel) + "  ";
			a += "Uk = " + to_string(ptr->binding_energy) + " > \n";
			i++;
		};*/
		return "";
	};
	

	

	
	double ONE = 1;
	Shell NULL_SHELL = Shell(&ONE, true);
	Shell *active_shell;

	void performRAD(Transition* _this, PARTICLES *particles){
		(particles->PHOTONS).push_back(_this->E);
		--(_this->j->frac);
		(active_shell->frac)++;
	};

	void performNONRAD(Transition* _this, PARTICLES *particles){
		(particles->ELECTRONS).push_back(_this->E);
		--(_this->j->frac);
		--(_this->k->frac);
		(active_shell->frac)++;
	};
	

	bool EQ;
	
	Transition *tr;

	void Atom::run(int shell_index, PARTICLES *particles, mixmax_engine *genPTR){ // cut off value should be chosen in a smart way, cut_off = min ( E(last M-Shell), simulation cut off ) 



		active_shell = shSTART + shell_index; // the cython class should ensure that this is in bounds

		if (active_shell->dontSIMULATE){ //cython class should take care of this
			return;
		};
		
		
		
		(active_shell->frac)--;
		/* --- CHOOSE THE FIRST TRANSITION ----  */
		R = genPTR->get_next_float()*(active_shell->Nt);
		i = (int) R;
		
		tr = active_shell->trSTART + i;
		
		
		//tr->perform(particles);
		
		
		if (R - i > tr->p){
			tr->t->perform(tr->t,particles);
		}else{
			tr->perform(tr, particles);
		};

		
		
		
		while (1) { /*NOW THE SAME THING, BUT WITH REJECTION SAMPLING*/
		
		
		
			/*  ----- SELECT NEXT ACTIVE SHELL ------ */
			
			


			EQ = active_shell->frac == active_shell->original_frac; //need to use this twice

			while ( EQ || active_shell->dontSIMULATE){ // TRUE = IGNORE ACTIVE SHELL
				
				if (!EQ){ // shell is marked dontSIMULATE, will be ignored, but has vacancies -> remove them
					active_shell->frac = active_shell->original_frac;
				};
				++active_shell;
				if (active_shell->LAST){ // increment shell and check if it's the last one
					active_shell->frac = active_shell->original_frac; // remove any vacancies before terminating 
					return;}; 
				
				EQ = active_shell->frac == active_shell->original_frac; // set up another condition to decide if will be ignored
			};
			
			
			
			/*  ----- CHOOSE TRANSITION ------ */
			do{
				
				// GENERATING PROPOSAL (WILERSONS ALIASING)
				R = genPTR->get_next_float()*active_shell->Nt;
				i = (int) R;

				tr = active_shell->trSTART + i;

				if (R - i > tr->p){
					tr = tr->t;
					};
					
			if (tr->j == tr->k) {
				  if (*(tr->j->frac - 1) == 0){
					  
					  continue;
					  
				  }
				
			}
				

			}while(genPTR->get_next_float() > *(tr->j->frac) * *(tr->k->frac) );

			tr->perform(tr, particles);
			
		};
		
		
		


		
	};
	
	
	
	
	void setADRESS_RAD(Transition * tr){
		tr->perform = &performRAD;
		tr->k = &NULL_SHELL;

	};
	
	
	void setADRESS_NONRAD(Transition * tr){
		tr->perform = &performNONRAD;
	};
	
	
	
	
};




			/*while (*(active_shell->frac) == 1 || active_shell->binding_energy < cut_off || active_shell->Nt == 0){
				
				
				if (*(active_shell->frac) != 1){
					active_shell->frac = active_shell->original_frac;
					};
				
				
				
				++active_shell;
				
				
				
				if (active_shell->shells_ahead == 0){
					//reset();
					active_shell->frac = active_shell->original_frac;
					return;
				};
			};
			
			*/















			/*  ----- CHECK EDGE CASES ------ */
			
			
			


			/*
			switch (active_shell->shells_ahead - *empty_shells){
				
				case 0:
					reset();
					return;
				
				
				case 1: // only one shell left, i.e., only one transition 
					// in all likelihood, this loop will just do about 2 iterations 
					for (Transition *TR = (active_shell->TRANSITIONS).start; TR < (active_shell->TRANSITIONS).end; ++TR){ // search for available transition
						if ( *(TR->j->frac) *  *(TR->k->frac) > 0){ 
							
							if (TR->k->Nt < 0 ){  // radiative, insert photon and end simulation 
								PHOTONS->push_back(tr);
								reset();
								return;
							};
							//if (tr->k == tr->j){continue;};
							ELECTRONS->push_back(tr);
							reset();
							return;
						}; 
					};
					reset();
					return;

			};*/