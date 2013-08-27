#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <R.h>

// counts the occurrence of each row in a data frame
void tabulate (int *data, int nrow, int ncol, int *dims, int *Freq) {
 
  int i;
  
  // int *cumDims  = malloc (ncol * sizeof(int));
  int *cumDims = (int *) R_alloc (ncol, sizeof(int));

  if (cumDims == NULL) 
    error ("Memory allocation error.");

  cumDims[ncol - 1] = 1;

  for (i = ncol - 2; i >= 0; i--)
    cumDims[i] = cumDims[i+1] * dims[i+1];    
 
  int index = 0;

  for (i = 0; i < nrow * ncol; i++) {
    index = index + data[i] * cumDims [i % ncol];   
    if ((i+1) % ncol == 0) {
      Freq[index]++;
      index = 0;      
    }
  }  
   
}

// collapses a contingency table over the response variable
void collapse (int *Freq, int lFreq, int *margFreq, int nRespCats) {

  int i,j;
  int k = 0;

  for (i = 0; i < lFreq; i+=nRespCats) {
    for (j = 0; j < nRespCats; j++)
      margFreq[k] += Freq[i+j]; 
    k++;
  }

}

// finds the log marginal likelihood of a regression
// function uses equation 17 in Reference [2]
void findLogMargLik (int *data, int *nrow, int *ncol, int *dims, double *alpha, double *logMargLik) {

  int i;
  int lFreq = 1;
  
  for (i = 0; i < ncol[0]; i++)
    lFreq *= dims[i];  
  
  // int *Freq = calloc (lFreq, sizeof(int));
  int *Freq = (int*) S_alloc (lFreq, sizeof(int));

  if (Freq == NULL) 
    error ("Memory allocation error");

  tabulate (data, nrow[0], ncol[0], dims, Freq); 
     
  double term1, term2;
  term1 = term2 = 0;

  double priorValue = alpha[0] / lFreq;

  for (i = 0; i < lFreq; i++)
    term1 += lgamma(Freq[i] + priorValue);
  term2 = lFreq * lgamma (priorValue);

  int lMargFreq = lFreq / dims[ncol[0]-1];
  // int *margFreq = calloc (lMargFreq, sizeof(int));
  int *margFreq = (int*) S_alloc (lMargFreq, sizeof(int));

  if (margFreq == NULL)
    error ("Memory allocation error.");

  collapse (Freq, lFreq, margFreq, dims[ncol[0]-1]);

  double term3, term4;
  term3 = term4 = 0; 

  priorValue = alpha[0] / lMargFreq;

  for (i = 0; i < lMargFreq; i++)
    term3 += lgamma(margFreq[i] + priorValue);  
  term4 = lMargFreq * lgamma(priorValue);

  logMargLik[0] = term1 - term2 + term4 - term3;
  
}
