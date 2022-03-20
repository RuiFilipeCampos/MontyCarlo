#pragma once
#ifndef __IntervalArithmeticsHH__
#define __IntervalArithmeticsHH__

#include <list>



class Interval {
    public:
        double t1, t2; 
        Interval();
        Interval(double t1, double t2);
        void print(void);
        void pprint();
};


typedef std::list<Interval> intLIST;

class intIterator {
    bool second;
    intLIST crosses;
    intLIST::iterator it;
    intLIST::iterator end;

    public:
        intIterator();
        intIterator(intLIST crosses);
        void inc();
        double deref();
        double current();

   
};

intLIST intIntersect(intLIST& left, intLIST& right);
intLIST intPlus(intLIST& left, intLIST& right);
intLIST intMinus(intLIST& left, intLIST& right);


#endif