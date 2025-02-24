#ifndef __Orbital_H__
#define __Orbital_H__

class Orbital{

 private:

   int l,m;
   double Power, Exponent;

 public:

   Orbital() {};
   ~Orbital(){}; 
   void initialize(int il, int im, double ipow, double iexp){
    l=il; m=im; Power=ipow; Exponent=iexp; 
   } 
   int get_l() { return l;}
   int get_m() { return m;}
   double get_power() { return Power;}
   double get_exponent() { return Exponent;}
};
#endif 
