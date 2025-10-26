#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<sys/time.h>

#include"../include/geometry_st.h"
#include"../include/random.h"

// when using C99 M_PI is not defined in math.h header!
#ifndef M_PI
  #define M_PI  3.141592653589793238462643383279502884
#endif

#define STDIM 2  // space-time dimensionality
#define STRING_LENGTH 50

static inline double wall_time_sec(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

void staple(double complex ** restrict lattice, 
            long int const * const restrict nnp,
            long int const * const restrict nnm,
            long int r,
            int dir,
            long int stvol,
            double complex *ris)
  {
//               ^ dir
//         (4)   |   (1)
//     +----<----+---->----+
//     |         |         |
//  (5)|         |         |
//     V         ^         V (2)
//     |         |         |
//     +---->----+----<----+---> i
//     r1  (6)   r   (3)
//
  
  int i, k;
  long int r1;
  double complex aux;
  
  *ris=0.0+I*0.0;
  
  for(k=1; k<STDIM; k++)
     {
     i=(dir+k)%STDIM;
 
     // forward
     aux=lattice[nnp[dirgeo(r,dir,stvol)]][i];    // 1
     aux*=conj(lattice[nnp[dirgeo(r,i,stvol)]][dir]);  // 2
     aux*=conj(lattice[r][i]);  //3
     *ris+=aux;

     // backward
     r1=nnm[dirgeo(r,i,stvol)];
     aux=conj(lattice[nnp[dirgeo(r1,dir,stvol)]][i]);  // 4
     aux*=conj(lattice[r1][dir]);
     aux*=lattice[r1][i];
     *ris+=aux; 
     }
  }


// return 1 if the update is accepted, 0 otherwise
int metropolis(double complex ** restrict lattice, 
               long int const * const restrict nnp,
               long int const * const restrict nnm,
               long int r,
               int dir,
               double epsilon,
               double beta,
               long int stvol)
  {
  int acc=0;
  double oldS, newS, deltaS, aux;
  double complex stap, change;

  aux=epsilon*(2.0*myrand()-1.0);
  change=1.0/sqrt(1.0+aux*aux) + I*aux/sqrt(1.0+aux*aux);  // complex phase that becomes 1 if epsilon=0
  
  // compute staple
  staple(lattice, nnp, nnm, r, dir, stvol, &stap);

  oldS=-beta*creal(lattice[r][dir]*stap);
  newS=-beta*creal(lattice[r][dir]*stap*change);
  deltaS=newS-oldS;

  if(deltaS<0) 
    {
    lattice[r][dir]*=change;
    acc=1;    
    }
  else
    {
    if(myrand()<exp(-deltaS))
      {
      lattice[r][dir]*=change;
      acc=1;    
      }
    }

  return acc;
  }


void overrelaxation(double complex ** restrict lattice, 
                    long int const * const restrict nnp,
                    long int const * const restrict nnm,
                    long int r,
                    int dir,
                    long int stvol)
  {
  double complex stap, aux;

  // compute staple
  staple(lattice, nnp, nnm, r, dir, stvol, &stap);

  if(cabs(stap)>1.0e-10)
    {
    aux=conj(lattice[r][dir]*stap*stap);
    aux/=(cabs(stap)*cabs(stap));

    lattice[r][dir]=aux;
    }
  }


// plaquette 
double plaquette(double complex ** restrict lattice, 
                 long int const * const restrict nnp,
                 long int stvol)
  {
  int i, j;
  long int r;
  double ris=0.0;
  double complex aux;

//     ^i
//     |  (2)  
//     +--->--+
//     |      |
// (1) ^      V (3)
//     |      |
//     +--<---+--->j
//     r  (4)

  for(r=0; r<stvol; r++)
     {
     for(i=0; i<STDIM; i++)
        {
        for(j=i+1; j<STDIM; j++)
           {
           aux=lattice[r][i];
           aux*=lattice[nnp[dirgeo(r,i,stvol)]][j];
           aux*=conj(lattice[nnp[dirgeo(r,j,stvol)]][i]);
           aux*=conj(lattice[r][j]);
 
           ris+=creal(aux);
           }
        }
     }

  ris/=(double)stvol;
  ris/=(0.5*(double)STDIM*(double)(STDIM-1));

  return ris;
  }


// topological charge 
double topch(double complex ** restrict lattice, 
             long int const * const restrict nnp,
             long int stvol)
  {
  #if STDIM==2
    long int r;
    double ris=0.0;
    double complex aux;
  
//     ^1
//     |  (2)  
//     +--->--+
//     |      |
// (1) ^      V (3)
//     |      |
//     +--<---+--->0
//     r  (4)

    for(r=0; r<stvol; r++)
       {
       aux=lattice[r][1];
       aux*=lattice[nnp[dirgeo(r,1,stvol)]][0];
       aux*=conj(lattice[nnp[dirgeo(r,0,stvol)]][1]);
       aux*=conj(lattice[r][0]);
   
       ris+=atan2(cimag(aux),creal(aux))/(2.0*M_PI);
       }
  
    return ris;
  #else
    return 0.0;
  #endif
  }

double compute_hamiltonian(double complex ** restrict lattice,
                           double *momentum,
                           long int const * const restrict nnp,
                           long int stvol,
                           double beta)
{
    double kinetic_energy = 0.0;
    double potential_energy = 0.0;

    for (long int r = 0; r < stvol; r++) {
        for (int i = 0; i < STDIM; i++) {
            kinetic_energy += 0.5 * momentum[dirgeo(r, i, stvol)] * momentum[dirgeo(r, i, stvol)];
        }
    }

    double nplaq_per_site = 0.5 * (double)STDIM * (double)(STDIM - 1);
    double total_plaquettes = nplaq_per_site * (double)stvol;
    potential_energy = -beta * plaquette(lattice, nnp, stvol)* total_plaquettes;

    return kinetic_energy + potential_energy;
}

void leapfrog(double complex ** restrict lattice,
              double *momentum,
              long int const * const restrict nnp,
              long int const * const restrict nnm,
              long int stvol,
              double beta,
              double epsilon,
              int nsteps)
{
    double complex staple_val;
    double force;

    // Half-step update for momenta
    for (long int r = 0; r < stvol; r++) {
        for (int i = 0; i < STDIM; i++) {
            staple(lattice, nnp, nnm, r, i, stvol, &staple_val);
            force = -beta * cimag(lattice[r][i] * staple_val);
            momentum[dirgeo(r, i, stvol)] += 0.5 * epsilon * force;
        }
    }

    // Full-step update for lattice
    for (int step = 0; step < nsteps; step++) {
        for (long int r = 0; r < stvol; r++) {
            for (int i = 0; i < STDIM; i++) {
                lattice[r][i] *= cexp(I * epsilon * momentum[dirgeo(r, i, stvol)]);
            }
        }

        // Full-step update for momenta (except the last step)
        if (step < nsteps - 1) {
            for (long int r = 0; r < stvol; r++) {
                for (int i = 0; i < STDIM; i++) {
                    staple(lattice, nnp, nnm, r, i, stvol, &staple_val);
                    force = -beta * cimag(lattice[r][i] * staple_val);
                    momentum[dirgeo(r, i, stvol)] += epsilon * force;
                }
            }
        }
    }

    // Final half-step update for momenta
    for (long int r = 0; r < stvol; r++) {
        for (int i = 0; i < STDIM; i++) {
            staple(lattice, nnp, nnm, r, i, stvol, &staple_val);
            force = -beta * cimag(lattice[r][i] * (staple_val));
            momentum[dirgeo(r, i, stvol)] += 0.5 * epsilon * force;
        }
    }
}


void hmc(double complex ** restrict lattice,
         long int const * const restrict nnp,
         long int const * const restrict nnm,
         long int stvol,
         double beta,
         double epsilon,
         int nsteps)
{
    double *momentum = (double *)malloc((unsigned long int)(STDIM * stvol) * sizeof(double));
    if (momentum == NULL) {
        fprintf(stderr, "Allocation problem for momenta (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    double complex **lattice_backup = (double complex **)malloc((unsigned long int)(stvol) * sizeof(double complex *));
    if (lattice_backup == NULL) {
        fprintf(stderr, "Allocation problem for lattice backup (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    for (long int r = 0; r < stvol; r++) {
        lattice_backup[r] = (double complex *)malloc((unsigned long int)(STDIM) * sizeof(double complex));
        if (lattice_backup[r] == NULL) {
            fprintf(stderr, "Allocation problem for lattice backup (%s, %d)\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        memcpy(lattice_backup[r], lattice[r], STDIM * sizeof(double complex));
    }

    for (long int r = 0; r < stvol; r++) {
        for (int i = 0; i < STDIM; i++) {
            double u = myrand();
            if (u < 1e-12) u = 1e-12;
            momentum[dirgeo(r, i, stvol)] = sqrt(-2.0 * log(u)) * cos(2.0 * M_PI * myrand());
        }
    }

    double initial_hamiltonian = compute_hamiltonian(lattice, momentum, nnp,stvol, beta);

    // leapfrog integration
    //printf("Debug Inside HMC, before leapfrog\n");
    leapfrog(lattice, momentum, nnp, nnm, stvol, beta, epsilon, nsteps);
    //printf("Debug Inside HMC, after leapfrog\n");

    double final_hamiltonian = compute_hamiltonian(lattice, momentum, nnp,stvol, beta);

    if (myrand() >= exp(initial_hamiltonian - final_hamiltonian)) {
        // Reject: restore the original lattice
        for (long int r = 0; r < stvol; r++) {
            memcpy(lattice[r], lattice_backup[r], STDIM * sizeof(double complex));
        }
    }

    for (long int r = 0; r < stvol; r++) {
        free(lattice_backup[r]);
    }
    free(lattice_backup);
    free(momentum);
}



// main
int main(int argc, char **argv)
   {
   int i, Nt, Ns;
   double beta, rand, plaq, Q;
   double complex **lattice;
   long int *nnp, *nnm;
   long int iter, sample, r, stvolume, acc, count;

   char datafile[STRING_LENGTH];
   FILE *fp;

   const int overrelax=5;
   const int measevery=1;
   const int unitarizeevery=10;
   const double epsilon=1;
   
   const unsigned long int seed1=(unsigned long int) time(NULL);
   const unsigned long int seed2=seed1+127;
   int algorithm = 0; // 0 default for Metro+OR std, 1 now do just OR as placeholder but i want to implement HB, 2 for HMC
   if (argc == 7) algorithm = atoi(argv[6]);
    
   if(argc != 7)
     {
     fprintf(stdout, "How to use this program:\n");
     fprintf(stdout, "  %s Nt Ns beta sample datafile alg\n\n", argv[0]);
     fprintf(stdout, "  Nt = temporal size of the lattice\n");
     fprintf(stdout, "  Ns = spatial size of the lattice (space-time dimension defined by macro STDIM)\n");
     fprintf(stdout, "  beta = coupling\n");
     fprintf(stdout, "  sample = number of drawn to be extracted\n");
     fprintf(stdout, "  datafile = name of the file on which to write the data\n\n");
     fprintf(stdout, "Compiled for:\n");
     fprintf(stdout, "  dimensionality = %d\n\n", STDIM);
     fprintf(stdout, "Output:\n");
     fprintf(stdout, "  plaquette  (topological charge)\n");
     fprintf(stdout, "  one line for each configuration\n");

     return EXIT_SUCCESS;
     }
   else
     {  
     // read input values 
     Nt=atoi(argv[1]);
     Ns=atoi(argv[2]);
     beta=atof(argv[3]);
     sample=atol(argv[4]);

     if(strlen(argv[5]) >= STRING_LENGTH)
       {
       fprintf(stderr, "File name too long. Increse STRING_LENGTH or shorten the name (%s, %d)\n", __FILE__, __LINE__);
       return EXIT_FAILURE;
       }
     else
       {
       strcpy(datafile, argv[5]);
       }
     }

   // initialize random number generator
   myrand_init(seed1, seed2);

   // compute the spacetime volume
   stvolume=Nt;
   for(i=1; i<STDIM; i++)
      {
      stvolume*=Ns;
      }

   // allocate the lattice
   // and next neighbors: nnp[dirgeo(r, i, volume)]= next neighbor in positive "i" direction of site r 
   lattice=(double complex **)malloc((unsigned long int)(stvolume)*sizeof(double complex *));
   if(lattice == NULL)
     {
     fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
     return EXIT_FAILURE;
     }
   for(r=0; r<stvolume; r++)
      {
      lattice[r]=(double complex *)malloc((unsigned long int)(STDIM)*sizeof(double complex));
      if(lattice[r] == NULL)
        {
        fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
        }
      }
   nnp=(long int *)malloc((unsigned long int)(STDIM*stvolume)*sizeof(long int));
   if(nnp == NULL){
     fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
     return EXIT_FAILURE;
     }
   nnm=(long int *)malloc((unsigned long int)(STDIM*stvolume)*sizeof(long int));
   if(nnm == NULL){
     fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
     return EXIT_FAILURE;
     }

   // initialize nnp and nnm
   init_neighbors_st(nnp, nnm, Nt, Ns, STDIM);

   // initialize lattice to ordered start
   for(r=0; r<stvolume; r++)
      {
      for(i=0; i<STDIM; i++)
         {
         lattice[r][i]=1.0+I*0.0;
         }
      }

   // open data file
   fp=fopen(datafile, "w");
   if(fp==NULL)
     {
     fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
     return EXIT_FAILURE;
     }

   acc=0.0;
   count=0;
   int prog = 0;

  // time stats for prog to new prog
   double t_start = wall_time_sec();
   double t_last_prog = t_start;
   long int it_last_prog = 0;

   for (iter = 0; iter < sample; iter++) {
        int new_prog = (int)((double)(iter + 1) / (double)sample * 100.0);
        if (new_prog > prog) {
            prog = new_prog;
            double now = wall_time_sec();
            long int iters_since = (iter + 1) - it_last_prog;
            double dt = now - t_last_prog;
            double per_it = iters_since > 0 ? dt / (double)iters_since : 0.0;
            double elapsed = now - t_start;
            double eta = ((double)sample - (iter + 1)) * per_it;

            printf("Progress: %d%% | %ld iters in %.3fs (%.6fs/iter) | elapsed %.1fs | ETA %.1fs\n", prog, iters_since, dt, per_it, elapsed, eta);
            fflush(stdout);

            t_last_prog = now;
            it_last_prog = iter + 1;
        }
    if (algorithm == 0) {
        // Original behavior: Metropolis and Overrelaxation
        rand = myrand();
        if (rand < 1.0 / (double)overrelax) {
            count++;
            // Metropolis
            for (r = 0; r < stvolume; r++) {
                for (i = 0; i < STDIM; i++) {
                    acc += metropolis(lattice, nnp, nnm, r, i, epsilon, beta, stvolume);
                }
            }
        } else {
            // Overrelaxation
            for (r = 0; r < stvolume; r++) {
                for (i = 0; i < STDIM; i++) {
                    overrelaxation(lattice, nnp, nnm, r, i, stvolume);
                }
            }
        }
    } else if (algorithm == 1) {
        // Overrelaxation only
        for (r = 0; r < stvolume; r++) {
            for (i = 0; i < STDIM; i++) {
                overrelaxation(lattice, nnp, nnm, r, i, stvolume);
            }
        }
    } else if (algorithm == 2) {
        // HMC
        //printf("Debug: Before HMC call\n");
        hmc(lattice, nnp, nnm, stvolume, beta, 0.01, 20); // 10 leapfrog steps
        //printf("Debug: After HMC call\n");
      }

    // Reunitarize
    if (iter % unitarizeevery == 0) {
        for (r = 0; r < stvolume; r++) {
            for (i = 0; i < STDIM; i++) {
                lattice[r][i] /= cabs(lattice[r][i]);
            }
        }
    }

    // Perform measurements
    if (iter % measevery == 0) {
        plaq = plaquette(lattice, nnp, stvolume);
        Q = topch(lattice, nnp, stvolume);

        fprintf(fp, "%.12f %.12f\n", plaq, Q);
    }
  }

   printf("Acceptance rate %f\n", (double)acc/(double)count/(double)stvolume/(double)STDIM);

   // close datafile
   fclose(fp);

   // free memory
   for(r=0; r<stvolume; r++)
      {
      free(lattice[r]);
      }
   free(lattice);
   free(nnp);
   free(nnm);

   return EXIT_SUCCESS;
   }


