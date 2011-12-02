#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PROGNAME "ising"


#define ITERMAX 100000
#define I0 100
#define T 4.4
#define H 0.0
//note that here J=1

/* Periodic boundary conditions */
#define plus(i) ((i<L)?i+1:1)
#define minus(i) ((i>1)?i-1:L)

main(int argc,char **argv)
{
    short int **S;
    int i,j,n,iter,sum,m,en,L,Eold,Enew;
    FILE * outfile;
    float w[5],x,z[2];
    float T0,T1,dT;
    long seed=-1;

    short int **smatrix(long nrl, long nrh, long ncl, long nch);
    float ran1(long *idum);

    outfile = fopen ("m_E.dat","w");
    setbuf ( outfile , NULL );

    
    if(argc!=2) {
        fprintf(stderr,"usage: %s Nspins \n",PROGNAME);
        exit(1);
    }

     if((L=atol(argv[1]))<=0) {
        fprintf(stderr,"%s: error: T0 %d should be positive\n",
            PROGNAME,L);
        exit(1);
    }

       

    S=smatrix(1,L,1,L);

    


    /* Initialization at infinite temperature */
    for(i=1;i<=L;i++) {
        for(j=1;j<=L;j++) {
            S[i][j]=(ran1(&seed)<0.5?+1:-1);
        }
    }

    
    for(iter=0;iter<=ITERMAX;iter++) {
        for(n=0;n<L*L;n++) {
            /* Choose a spin at random */
            i=ran1(&seed)*L+1;
            j=ran1(&seed)*L+1;
            Eold=-(S[plus(i)][j]+S[minus(i)][j]+S[i][plus(j)]+S[i][minus(j)]+H)
                *S[i][j];
            Enew= (S[plus(i)][j]+S[minus(i)][j]+S[i][plus(j)]+S[i][minus(j)]+H)
                *S[i][j];//you don't really need that
            
            x=exp(-2*Eold/T);// exp(DeltaE/T)
            /* Monte Carlo update */
            if(x>1)S[i][j]=-S[i][j];//if DeltaE<0
            else {
                if(ran1(&seed)<x)S[i][j]=-S[i][j];
            }
        }
        
        /* Output magnetization and internal energy per spin every I0 steps */
        if(iter%I0==0) {
            m=en=0;
            for(i=1;i<=L;i++) {
                for(j=1;j<=L;j++) {
                    m+=S[i][j];             en-=(S[plus(i)][j]+S[minus(i)][j]+S[i][plus(j)]+S[i][minus(j)]) *S[i][j];
                }
            }
            fprintf(outfile,"%d %g %g\n",iter,(float)m/(L*L),(float)en/(2*L*L));
        }
    }

}


/*-------------------- smatrix -----------------------------*/
#define NR_END 1

short int **smatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a short int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    short int **m;

    /* allocate pointers to rows */
    m=(short int **) malloc((size_t)((nrow+NR_END)*sizeof(short int*)));
    if (!m) {
        fprintf(stderr,"%s: error: allocation failure 1 in smatrix()\n",
            PROGNAME);
        exit(1);
    }
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl]=(short int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(short int)));
    if (!m[nrl]) {
        fprintf(stderr,"%s: error: allocation failure 2 in smatrix()\n", PROGNAME);
        exit(1);
    }
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

    /* return pointer to array of pointers to rows */
    return m;
}


#undef NR_END

/*-------------------------------------------------------------*/

/*------------------------ran1.c-------------------------------*/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0) *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}


#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/*-------------------------------------------------------------------*/

