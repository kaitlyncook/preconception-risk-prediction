/*
 TO COMPILE USE THE CODE:
 
 R CMD SHLIB loop.c -lgsl -lgslcblas
 
 */

#include <stdio.h>
#include <math.h>
#include <time.h>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_sort_vector.h"
#include "gsl/gsl_sf.h"

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#include "R.h"
#include "Rmath.h"

/*
#include "loop.h"
*/
#define Pi 3.141592653589793238462643383280


void ej_func(gsl_vector *e_j,
                 int j)
{
    gsl_vector_set_zero(e_j);
    gsl_vector_set(e_j, j, 1);
    
    return;
}

double L2_norm(gsl_vector *x1,
               gsl_vector *x2)
{
    int p = x1 -> size;
    int i;
    double val = 0;
    
    for(i = 0; i < p; i++)
    {
        val += pow(gsl_vector_get(x1, i)-gsl_vector_get(x2, i), 2);
    }
    return pow(val, 0.5);
}

/*
Code credit: https://stackoverflow.com/questions/14003487/round-in-c-with-n-digits-after-decimal-point
*/

double round1(double num, int N) {
      double p10 = pow(10,N);
      return round(num* p10) / p10;
}

/* */
void loop(double Data[],
		  double Prob[],
          double Seq[],
          double Cat1[],
          double Cat2[],
          double Cat3[],
          double Cat4[],
          int *n,
          int *p,
          int *I,
          int *J,
          int *K,
          int *L,
          double *CRval,
          double *Denom)
{
    GetRNGstate();
    
    int i, j, k, l, m, s, t, indexMin, ties;
    double val, rval, pi, pj, pk, pl, lossMin;
    
    gsl_matrix *riskMat     = gsl_matrix_calloc(*n, *p);
    gsl_matrix *seq     	= gsl_matrix_calloc(24, 4);
    gsl_vector *indCat1     = gsl_vector_calloc(*I);
    gsl_vector *indCat2     = gsl_vector_calloc(*J);
    gsl_vector *indCat3     = gsl_vector_calloc(*K);
    gsl_vector *indCat4     = gsl_vector_calloc(*L);
    
    gsl_vector *p1     = gsl_vector_calloc(4);
    gsl_vector *p2     = gsl_vector_calloc(4);
    gsl_vector *p3     = gsl_vector_calloc(4);
    gsl_vector *p4     = gsl_vector_calloc(4);
    gsl_vector *e_j     = gsl_vector_calloc(4);
    
    gsl_vector *min_vec = gsl_vector_calloc(24);
    
    for(i = 0; i < *n; i++)
    {
        for(j = 0; j < *p; j++)
        {
            gsl_matrix_set(riskMat, i, j, Data[(j* *n) + i]);
        }
    }
    
    for(i = 0; i < 24; i++)
    {
        for(j = 0; j < 4; j++)
        {
            gsl_matrix_set(seq, i, j, Seq[(j* 24) + i]-1);
        }
    }
    
    for(i = 0; i < *I; i++)
    {
        gsl_vector_set(indCat1, i, Cat1[i]-1);
    }
    for(i = 0; i < *J; i++)
    {
        gsl_vector_set(indCat2, i, Cat2[i]-1);
    }
    for(i = 0; i < *K; i++)
    {
        gsl_vector_set(indCat3, i, Cat3[i]-1);
    }
    for(i = 0; i < *L; i++)
    {
        gsl_vector_set(indCat4, i, Cat4[i]-1);
    }
    
    
    /*
    for(i = 0; i < 180; i++)
    {
        for(j = 0; j < *p; j++)
        {
            Rprintf("Mat_%d_%d = %.10f\n", i+1, j+1, gsl_matrix_get(riskMat, i, j));
        }
    }
    
    for(i = 0; i < 24; i++)
    {
        for(j = 0; j < 4; j++)
        {
            Rprintf("Seq_%d_%d = %.10f\n", i+1, j+1, gsl_matrix_get(seq, i, j));
        }
    }
    */
    
    for(i = 0; i < *I; i++)
    {
        for(j = 0; j < *J; j++)
        {
            for(k = 0; k < *K; k++)
            {
                for(l = 0; l < *L; l++)
                {
                    /*if((gsl_matrix_get(riskMat, i, 4) < gsl_matrix_get(riskMat, j, 4)) && (gsl_matrix_get(riskMat, j, 4) < gsl_matrix_get(riskMat, k, 4)) && (gsl_matrix_get(riskMat, k, 4) < gsl_matrix_get(riskMat, l, 4)))
                    { */
                        for(m = 0; m < 4; m++)
                        {
                            gsl_vector_set(p1, m, gsl_matrix_get(riskMat, (int) gsl_vector_get(indCat1, i), m));
                            gsl_vector_set(p2, m, gsl_matrix_get(riskMat, (int) gsl_vector_get(indCat2, j), m));
                            gsl_vector_set(p3, m, gsl_matrix_get(riskMat, (int) gsl_vector_get(indCat3, k), m));
                            gsl_vector_set(p4, m, gsl_matrix_get(riskMat, (int) gsl_vector_get(indCat4, l), m));
                        }
                        
                        pi = Prob[(int) gsl_vector_get(indCat1, i)];
                        pj = Prob[(int) gsl_vector_get(indCat2, j)];
                        pk = Prob[(int) gsl_vector_get(indCat3, k)];
                        pl = Prob[(int) gsl_vector_get(indCat4, l)];
                        
                        gsl_vector_set_zero(min_vec);
                        for(s = 0; s < 24; s++)
                        {
                            ej_func(e_j, (int) gsl_matrix_get(seq, s, 0));
                            val = L2_norm(p1, e_j);
                            
                            ej_func(e_j, (int) gsl_matrix_get(seq, s, 1));
                            val += L2_norm(p2, e_j);
                            
                            ej_func(e_j, (int) gsl_matrix_get(seq, s, 2));
                            val += L2_norm(p3, e_j);
                            
                            ej_func(e_j, (int) gsl_matrix_get(seq, s, 3));
                            val += L2_norm(p4, e_j);
                            
                            rval = round1(val, 10);
                            
                            gsl_vector_set(min_vec, s, rval);
                        }
                        
                        indexMin = (int) gsl_vector_min_index(min_vec);
                        lossMin = gsl_vector_min(min_vec);
                        
                        ties = 0;
                        
                        for(t=0; t < 24; t++)
                        {
                        	if(gsl_vector_get(min_vec, t)==lossMin) ties++;	
                        }

                        *CRval = (indexMin == 0) ? *CRval+1/(ties*pi*pj*pk*pl) : *CRval;
                        *Denom += 1/(pi*pj*pk*pl);
                        
                        /*
                        Rprintf("CRval = %.3f\n", *CRval);
                         */
                    /*}*/
                    
                }
            }
        }
    }
    
    
    PutRNGstate();
    return;
    
    
    
    
}


