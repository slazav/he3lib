#include "../he3.h"
#include <unistd.h>
#include <string.h>

main(){
  double ttc=0.2; /* T/Tc */
  double p=10;    /* P, bar */
  double nuB;

  // Leggett frequency
  nuB = he3_nu_b_(&ttc, &p);
  printf("nu_b = %.3f kHz\n", nuB/1000);
  printf("m_3  = %.3e g\n", he3_amass_);
}
