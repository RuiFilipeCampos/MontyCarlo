// IntervalArithmetics.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "IntervalArithmetics.h"
#include <list>

typedef std::list<Interval> intLIST;


Interval::Interval(double t1, double t2) {
        this->t1 = t1;
        this->t2 = t2;
};

void Interval::print(void) {
    std::cout << "[" << this->t1 << ", " << this->t2 << "] \n" ;
};

Interval::Interval() {};

void Interval::pprint() {
    std::cout << "[" << this->t1;
    double dt = this->t2 - this->t1;
    int N = (int)(dt / .1);
    while (N != 0) { 
        std::cout << "-";
    --N; }
    std::cout << this->t2 << "]";
}

intLIST intIntersect(intLIST& left, intLIST& right) {
    intLIST result;

    intLIST::iterator first;
    intLIST::iterator last;

    for (auto l_it = left.begin(); l_it != left.end(); ++l_it) {
        for (auto r_it = right.begin(); r_it != right.end(); ++r_it) {
            
            // order the intervals -> should be a way to make this branchless
            if (l_it->t1 < r_it->t1) {
                first = l_it;
                last = r_it;
            }
            else {
                first = r_it;
                last = l_it;
            };

            // they don't intersect
            if (first->t2 < last->t1) { continue; };

            // apply normal rules
            result.emplace_back(std::max(first->t1, last->t1), std::min(first->t2, last->t2));
        };
    };

    return result;
};

intLIST intPlus(intLIST& left, intLIST& right) {
    /*
    The result is the list "left".
    */

    auto l_it = left.begin();
    auto r_it = right.begin();

    auto l_end = left.end();
    auto r_end = right.end();

    intLIST::iterator move_this;

    do {

        if (l_it->t2 < r_it->t1) { // TARGET INTERVAL IS AHEAD -> MOVE 

           //  l --[t1, t2]-------------------some other interval----
           //  r ---------------[t1, t2]-------------------------------

            ++l_it;

        } else if (l_it->t1 > r_it->t2) { // TARGET IS BEHIND -> PLACE BEHIND 

            //  l --------------[t1, t2]----
            //  r --[t1, t2]----------------

            move_this = r_it;
            ++r_it;
            left.splice(l_it, right, move_this);

        }else{  // INTERVALS INTERSECT
            
            //  l -------[t1------t2]-------- //
            //  r --[t1-------t2]------------ //

            l_it->t1 = std::min(r_it->t1, l_it->t1);

            //  l --[t1------t2]------???????????
            //  r --[t1-------------------t2]-----
            
            intLIST::iterator probe = l_it;

            while (probe->t2 < r_it->t2) {

                //  l --[t1------t2]---[t1 ----- t2]----
                //  r --[t1--------------------------------t2]---------

                ++probe;

                if (probe == l_end) {
                    l_it->t2 = r_it->t2;
                    return left;
                };
            };
            if (r_it->t2 > probe->t1) {
            //  l --[t1------t2]--------[t1 ----- t2]----
            //  r --[t1----------------------t2]--????????--------
                l_it->t2 = probe->t2;
            } else {
                l_it->t2 = r_it->t2;
                //  l --[t1------t2]-----------[t1 ----- t2]----
                //  r --[t1--------------t2]----------
                ++l_it;
            };
            ++r_it;
        };

        if (l_it == l_end) {
            if (r_it != r_end) {
                left.splice(l_end, right, r_it, r_end);
            }
            return left;
        };
    } while (r_it != r_end);
    return left;
};

intLIST intMinus(intLIST& left, intLIST& right) {
    auto l_it = left.begin();
    auto r_it = right.begin();

    auto l_end = left.end();
    auto r_end = right.end();

    intLIST::iterator move_this;

    do {

        if (l_it->t2 < r_it->t1) { // TARGET INTERVAL IS AHEAD -> MOVE 

           //  l --[t1, t2]-------------------some other interval----
           //  r ---------------[t1, t2]-------------------------------

            ++l_it;

        }
        else if (l_it->t1 > r_it->t2) { // TARGET IS BEHIND -> IGNORE
         //  l --------------[t1, t2]----
         //  r --[t1, t2]----------------

            //move_this = r_it;
            ++r_it;
            //left.splice(l_it, right, move_this);

        }
        else {  // INTERVALS INTERSECT

                //  l -------[t1-------t2]-------- //
                //  r -------[t1-------t2]------------ //

            if (l_it->t1 < r_it->t1) {

                if (l_it->t2 > r_it->t2) {
                    //  l -------[t1-----------t2]-------- //
                    //  r ----------[t1----t2]------------ //

                    //Interval I1 = { l_it->t1, r_it->t1 };
                    left.insert(l_it, { l_it->t1, r_it->t1 });
                    l_it->t1 = r_it->t2;

                    ++r_it;


                }
                else {
                    //  l --[t1------t2]-------- //
                    //  r -------[t1-------t2]------------ //

                    l_it->t2 = r_it->t1;

                    ++l_it;


                    if (l_it == l_end) { return left; };

                    while (l_it->t2 < r_it->t2) {
                        //  l --[t1--t2]----[t1---t2]------------- //
                        //  r -------[t1--------------t2]------------ //
                        l_it = left.erase(l_it);
                        if (l_it == l_end) { return left; };
                    };

                    if (l_it->t1 < r_it->t2) {
                        //  l --[t1--t2]-----------[t1---t2]------------- //
                        //  r -------[t1--------------t2]------------ //
                        l_it->t1 = r_it->t2;
                    };
                    //  l --[t1--t2]-------------------[t1---t2]------------- //
                    //  r -------[t1------------t2]------------ //
                    ++r_it;
                };

            }
            else {
                if (l_it->t2 < r_it->t2) {
                    //  l ---------[t1-----t2]-------- //
                    //  r ------[t1-------------t2]------------ //
                    l_it = left.erase(l_it);

                }
                else {
                    //  l --------------[t1---------t2]-------- //
                    //  r -------[t1----------t2]------------ //
                    l_it->t1 = r_it->t2;
                    ++r_it;
                };
            };


        }


        if (l_it == l_end) {
            return left;
        };

    } while (r_it != r_end);
    return left;
};


void print_intLIST(intLIST& the_list, double start) {
    


    auto first = the_list.begin();
    double dt = first->t1 - start;
    int N = (int)(dt / .1);

    while (N != 0) {
        std::cout << "*";
        --N;
    }



    first->pprint();


    auto it = first;
    ++it;

    for (it; it != the_list.end(); ++it) {

        double dt = it->t1 - first->t2;

        int N = (int)(dt / .1);

        while (N != 0) {
            std::cout << "*";
            --N;
        }
        it->pprint();

        ++first;
    }
};

intIterator::intIterator() {};


intIterator::intIterator(intLIST crosses) {
    this->crosses = crosses;
    this->it = this->crosses.begin();
    this->end = this->crosses.end();

    this->second = this->it->t1 < 0;


};

double intIterator::current() {
    if (this->it == this->end) {
        return std::numeric_limits<double>::infinity();
    };

    if (this->second) {
        return this->it->t2;
    };

    return this->it->t1;

};

void intIterator::inc() {

    if (this->second) {
        
        ++(this->it);
        this->second = false;
        return;
    };

    this->second = true;
};






int main() { return 0; };
/*
int main()
{

    std::list<double> L = { 1, 2, 3, 4 };
    std::list<double> R = { -1, -2, -3, -4 };

    std::list<double>::iterator insert_here = L.begin();
    std::list<double>::iterator take_this = R.begin();
    std::list<double>::iterator itt = L.begin();

    ++itt;
    ++itt;
    ++insert_here;
    ++insert_here;

    ++take_this;
    std::cout << *insert_here << "\n"; 
    L.splice(insert_here, R, take_this);
    
    for (auto it = L.begin(); it != L.end(); ++it) {
        std::cout << *it << "\n";
    
    };

    std::cout << "\n itt = " << *itt << "\n";
    std::cout << *insert_here << "\n";





    std::list<Interval> crosses1;
    std::list<Interval> crosses2;
    intLIST result;
    crosses1.emplace_back(-4, -3);
    crosses1.emplace_back(-1., 2.); crosses2.emplace_back(0.5, 2.5);
    crosses1.emplace_back(3., 4.2); crosses2.emplace_back(3.1, 3.9);

    crosses2.emplace_back(4., 5);




    std::cout << "\n";
    print_intLIST(crosses1, -5);
    std::cout << "\n";

    print_intLIST(crosses2, -5);
    std::cout << "\n";





    std::cout << "Intersecting... \n";


    std::cout << "\n";
    print_intLIST(crosses1, -5);
    std::cout << "\n";

    print_intLIST(crosses2, -5);
    std::cout << "\n";

    result = intersect(crosses1, crosses2);

    print_intLIST(result, -5);
    std::cout << "\n";

    std::cout << "\n";
    print_intLIST(crosses1, -5);
    std::cout << "\n";

    print_intLIST(crosses2, -5);
    std::cout << "\n";

    minus(crosses1, crosses2);

    print_intLIST(crosses1, -5);
    std::cout << "\n";



    std::cout << "Hello World!\n";
};
*/
// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
