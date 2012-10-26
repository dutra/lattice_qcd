// gcc -lm ho_gswv.c -o ho_gswv


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int N = 8;
int npaths = 200000;
int nrepeat = 5;
double a = 0.5;

void update(double *, int);
double S(int, double *);
double compute_O(double *);
double MC(double *, double *);
double randr();

double randr()
{
  double scaled = (double)rand()/(RAND_MAX+1.0);
  //  printf("%e\n", scaled);
  return scaled;
}

void update(double *x, int repeat) {
  int i;
  for(i=0; i<repeat; i++) {
    int j;
    for(j=1; j<N-1; j++) {
      double old_x = *(x+j);
      double old_S = S(j,x);
      double dS = 0;
    
      *(x+j) += 2.0*(1.0-2.0*randr());
      //      printf("%e\n", *(x+j));
      dS = S(j,x) - old_S;
      //      printf("%e\n", dS);

      if((dS > 0) && (exp(-dS) < randr())) {
	*(x+j) = old_x;
      }
    }
  }
  return;
}

double S(int j, double *x) {
  return a* ( 0.5*((x[j+1]-x[j])/a)*((x[j+1]-x[j])/a) + 0.5*x[j]*x[j]);
}

double compute_O(double *x) {
  double g;
  double s=0;
  int i;
  for(i=0; i<N; i++)
    s += S(i,x);
  //  printf("%e\n",s);
  return exp(-s);

}

double MC(double *x, double *G) {
  int i;
  double avg=0;

  update(x, 10*nrepeat);

  for(i=0; i<npaths; i++) {
    update(x,nrepeat);
    *(G+i) = compute_O(x);
  }
  for (i=0; i<npaths; i++)
    avg += *(G+i);
  return avg;
}

int main() {

  int i,j;
  srand ( time(NULL) );
  for(j=0; j<100; j++) {
    double x[N];
    double G[npaths];
  
    x[0]=x[N-1]=((double)j)/50.0;

    for(i=1; i<N-1; i++)
      x[i] = 0; 
    printf("%e, %e\n", j/50.0, MC(x,G));
  }

  return 0;
}
