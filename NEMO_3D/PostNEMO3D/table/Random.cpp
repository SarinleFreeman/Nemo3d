#include "Random.h"

void Random::Init(unsigned int Seed)
 {
  n=Seed; v=w=x=y=z=c=0; for(int i=0; i<8; i++) z=Next();
 }

unsigned int Random::Next(void) // period 2^192
// G. Marsaglia, A. Zaman, Comp. in Phys. 8, 117 (1994)
 {
  unsigned int s=v+c; c=(y<s); s=y-s; if(c) s-=10U;
  v=w; w=x; x=y; y=z; z=s; n=n*69069U+1013904243U;
  return s+n;
 }

double Random::Normal(void)
 {
  static bool b=true;
  static double old=0.0;
  if(b)
   {
    b=false;
    double v1, v2, vv;
    do {v1=SignedUniform(); v2=SignedUniform(); vv=v1*v1+v2*v2;}
      while(vv>=1.0 || vv==0.0);
    vv=sqrt(-2.0/vv*log(vv));
    old=v1*vv;
    return v2*vv;
   }
  else
   {
    b=true;
    return old;
   }
 }

