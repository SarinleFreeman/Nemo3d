#ifndef __Random_H__
#define __Random_H__

#include <math.h>

const double Pi=M_PI, TwoPi=2.0*Pi, RandomBase=1.0/(0xFFFFFFFF+1.0),
  TwoPiRandomBase=TwoPi*RandomBase, TwoRandomBase=2.0*RandomBase;

class Random
 {
  unsigned int n, v, w, x, y, z, c;
  public:
  Random(void) {Init(0);}
  Random(unsigned int Seed) {Init(Seed);}
  void Init(unsigned int Seed);
  unsigned int Next(void);
    // uniform random integer, 0 <= i < 2^32
  signed int SignedNext(void) {return (signed int)Next();}
    // uniform random integer, -2^31 <= i < 2^31
  double Uniform(void) {return (Next()*RandomBase+Next())*RandomBase;}
    // returns uniform deviate, 0 <= x < 1
  double UniformPhase(void) {return (Next()*RandomBase+Next())*TwoPiRandomBase;}
    // returns uniform deviate, 0 <= x < 2*Pi
  double SignedUniform(void)
   {
    return (Next()*RandomBase+SignedNext())*TwoRandomBase;
   }
    // returns uniform deviate, -1 <= x < 1
  double Normal(void);
    // returns normal deviate, mean = 0, variance = 1
 };

#endif

