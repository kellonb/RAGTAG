/*******************************************||********************************************
                           Genetic algorithm optimizer of RAGTAG                               *
                                      genA.cu                                            *
Runs iterations of genetic algoirthm to optimize molecular mechanics dihedral parameters * 
                                                                                         *
              @author James Maier, Kellon Belfon, Chuan Tian                             *
              @lab Carlos Simmerling lab, Stony Brook University                         *
              @version 3.0 2019 Feb                                                      *
********************************************||*******************************************/
/*****************************************************************************************
* 	                ---------------LOAD LIBRARIES-------------                       *  
*****************************************************************************************/
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/generate.h>
#include <thrust/device_ptr.h>
#include <list>
#include <map>
#include "load.h"
#include "parse.h"
using namespace std;

/******** Number of threads for a given block, 256 block threads (index 0 to 255) *******/
const int BLOCK_SIZE=256;

#define HANDLE_ERROR(x) x;


/*****************************************************************************************
*                  Defining the six pivotal functions for the genetic algorithm          *
*  (1) mateIt, (2) mutateIt, (3) scoreIt, (4) calcAreas, (5) moveEm, (6) getSumAreas     *
* note: getSumAreas uses two other functions sumEm and sumEmIndex                        *
*****************************************************************************************/

/*****************************************************************************************
*                                | function1: mateIt |                                   *
*                                                                                        *
* @purpose creates offspring from a population, generating crossovers according to pCross*
* @param Vs a global array of all the parent and child genomes (Amplitude parameters)    *
* @param ptrs array of pointers from logical indices to actual indices into Vs for       * 
*        each individual                                                                 *
* @param areas the probabilities for choosing each individual for mating                 *
* @param sumArea pointer to the sum of all the individual areas                          *
* @param rands array of random numbers for crossover                                     *         
* @param pCross probability that crossover occurs                                        *               
* @param pSize number of individuals in the population (possible amplitudes solutions)   *
* @param genomeSize number of genes in a genome (number of dihedral * periodicity)       *
*****************************************************************************************/

__global__ void mateIt(float *Vs, int *ptrs, const float *areas, const float *sumArea, 
        const float *rands, const float pCross, const int pSize, const int genomeSize)
{
  /* figure out index for threads  blockId.x is the index for blocks, 
     blockDIM.x is the elements per blocks (# of threads in a block)
     threadIdx is the index for threads */
  int i=blockIdx.x * blockDim.x + threadIdx.x;

  /* random numbers for first parent, second parent, crossover */
  int randi=i*3;

  /* multiply i by 2, as we will have 2 parents and 2 offspring using a left bitwise 
  (<<) by 1*/
  i<<=1;

  /* if we're in the population (sometimes warps may go past) */ 
  if (i<pSize) {
    int parent[2];
    int j;
  /* figure out parents */
    parent[0]=parent[1]=-1;
  /* find parent where cumulative (cum) area (A) is less than random target (tgt) area
    selection of parents depends on cumulative probability being less than the 
    random probabilities (random numbers). rands[] is a function that returns a 
    pseudo-random number (pseudo random number generator)
    The pseudo random numbers range from 0 to RAND_MAX(default is 32767) 
    The random probabilities (tgtA) is random numbers multiply by sum of all the 
    individual probabilities*/
   
    float cumA=0.0f, tgtA=rands[randi++]* *sumArea; //tgtA random number from 0 to the sumArea
    while(cumA<=tgtA){
      ++parent[0];
      cumA+=areas[ptrs[parent[0]]/genomeSize]; // areas is e^-(scorei/score0)
      /* rands[randi-1] is the index back to zero since it is the first set of parents */
    }

#if DEBUG>2
    printf("rands[%d] ; %f ; %f=%f * %f\n",randi, cumA, tgtA, rands[randi-1], *sumArea);
    printf("first parent\n");
#endif

    /* This substract 1st parent area from sum of area  */
    cumA=0.0f; tgtA=rands[randi++]* (*sumArea-areas[ptrs[parent[0]]/genomeSize]); 
    while (cumA<=tgtA){
      ++parent[1];
      if (parent[1]==parent[0])  //Ensure you don't pick the same parents
        ++parent[1];
      cumA+=areas[ptrs[parent[1]]/genomeSize];
    }
#if DEBUG>2
    printf("Make offspring %d from %d and %d (%f=%f*(%f-%f)) %d\n", i, parent[0], 
       parent[1], tgtA, rands[randi-1], *sumArea, areas[ptrs[parent[0]]/genomeSize], randi);
#endif
    /* add offset of pSize to i because it is a child (next population) */
    i+=pSize;
    /* use ptrs to get indices into Vs */
    int i0=ptrs[i], i1=ptrs[i+1];
    parent[0]=ptrs[parent[0]];
    parent[1]=ptrs[parent[1]];
    /* set j to index for the next set of Vs */
    j=i0+genomeSize;
    /* put parent[0], parent[1], and i1 relative to i0, so we can just add i0 for index */
    parent[0]-=i0;
    parent[1]-=i0;
    i1-=i0;
    /* start with crossover pt at the end (no crossover) */
    int crossPt=j;
    /* check if we need to do crossover, 
       only do crossover if random number is less than pCross */
    if(rands[randi]<pCross){
      crossPt=i0+1+(int)(rands[randi]/pCross*(float)(genomeSize-1));
    }
    while(i0<crossPt){
      /* load next bit from parent and increment i */
      Vs[i0]=Vs[parent[0]+i0];
      Vs[i1+i0]=Vs[parent[1]+i0];
      ++i0;
    }
    while(i0<j){
      Vs[i0]=Vs[parent[1]+i0];
      Vs[i1+i0]=Vs[parent[0]+i0];
      ++i0;
    }  //end of while loop
  } // end of if i<pSize loop 
}

/*****************************************************************************************
                                | function 2: mutateIt |

 * @brief introduces mutations to the genomes in Vs, according to probability pMut, 
    with a max perturbation of max
 *
 * @param Vs a global array of all the parent and child genomes
 * @param ptrs array of pointers from logical indices to actual indices into Vs for
     each individual
   @param rands array of random numbers
 * @param pSize number of individuals in the population
 * @param pMut probability that a mutation occurs, evaluated for each gene
 * @param max maximum perturbation to an allele
 * @param genomeSize number of genes in a genome
*******************************************************************************************/

__global__ void mutateIt(float *Vs, int *ptrs, const float *rands, const int pSize, const float pMut, const float max, const int genomeSize, const float *AZP )
{
  /* figure out index */
  int i=blockIdx.x * blockDim.x + threadIdx.x;
  if(i<pSize){
    // get index into random number array
    int r=i*genomeSize;
    i=ptrs[i];
    int j=i+genomeSize;
    // want random numbers from [-max, max). will subtract max later
    float scale=2.0f*max/pMut;
    int az = -1;
    int mv;
    int indx;
    // iterate through genome
    while(i<j){
      az++;
      //printf("i is %d, j is %d, az is %d, azp is %f\n", i,j,az,AZP[az]);
      //printf("Vs before random mutation is %f\n", Vs[i]);
      if (AZP[az] < 0.0){
           mv = i + (-AZP[az]);
           Vs[i] =  Vs[mv];
           indx = i;
           //printf("in if: Vs[i] %f, Vs[mv] %f\n", Vs[i],Vs[mv]);
      }
      Vs[mv] = Vs[indx];
      //printf("in if: Vs[mv] %f, Vs[indx] %f\n", Vs[mv],Vs[indx]);
      if(rands[r]<pMut){
         // if Vs equal to zero (if amplitude is zero, leave it as zero)
        if (Vs[i] == 0.0) {
          // Ensure if Vs start off as zero it remain zero]
          //printf("Vs before zero is %f\n", Vs[i]);
          Vs[i] = 0.0;
          //printf("Vs after zero is %f\n", Vs[i]);
        }
        else {
          // mutate the amplitude by adding perturbation based on max, random number and pMut
          //printf("Vs before mutation is %f\n", Vs[i]);
          Vs[i]+=rands[r]*scale-max;
          // if AZP is negative then make Vs for that periodicity equal to the Vs plus the number
          // This is a way to ensure two Vs be the same value
          if (AZP[az] < 0.0){
            mv = i + (-AZP[az]);
            //printf("BEFORE: AZP is %f, mv is %d, i is %d, Vs[i] %f =  Vs[i+mv] %f\n", AZP[az],mv,i,Vs[i],Vs[i+mv]); 
            Vs[i] =  Vs[mv];
            indx = i;
            //printf("AZP is %f, mv is %d, i is %d, Vs[i] %f =  Vs[i+mv] %f\n", AZP[az],mv,i,Vs[i],Vs[i+mv]); 
            //printf("Vs after mutation is %f\n", Vs[i]);
          }
          Vs[mv] = Vs[indx];
          //printf(" and Vs after mutation is %f\n", Vs[i]);
        }
      }
      ++i;
      ++r;
    } // end of while loop
  } 
}

/************************************************************************************************
                                | function 3: scoreIt | 

 * @brief calculates a score indicating the closeness of fit for each individual/chromosome
   (set of parameters) against the training set
 * @param scores score for each conformation, calculated here, output array
 * @param areas weighting for each conformation, no longer need
 * @param Vs a global array of all the parent and child genomes (amplitudes)
 * @param ptrs array of pointers from logical indices to actual indices into Vs for each individual
 * @param tset training set
 * @param tgts targets for training
 * @param wts weights of each point in the training set
 * @param breaks breaks in training set, where different data should not be compared across breaks
 * @param nConf number of conformations in training set
 * @param pSize number of individuals in the population
 * @param genomeSize number of genes in a genome
 * @param xx space to store energy differences for each conformation with test parameters
************************************************************************************************/

__global__ void scoreIt(float *scores, float *areas, const float *Vs, const int *ptrs, const float *tset, const float *tgts, const float *wts, const int *breaks, const int nConf, const int pSize, const int genomeSize, float *xx )
{
  // i represent a chromosome , a set of amplitude parameters, this function will be done for each i (chromosome) at the same time
  int i=blockIdx.x * blockDim.x + threadIdx.x;
  if(i<pSize){
    float *x=xx+i*nConf;  // for the error of each conformation
    // get reference to score, S is the AAE 
    float *S=scores+i;
    // set score to 0
    *S=0.0f;
    // accumulate little s for each set
    float s;
    // get first index in genome
    int i0=ptrs[i];
    // get index of next genome space for looping bounds
    int j=i0+genomeSize;
    // start with the first element in the training set
    int t=0;
    /* start at break 0 */
    int b=0;
    /* loop over conformations c */
    int c=0;
    while(c<nConf){
      //s is the sum of REE 
      s=0.0f;
      /* loop only in units without break points */
      while(c<breaks[b+1]){
      /* start with delta E (tgts) for a given conformation (c) within a break; see load.cpp 
           conf (c) goes through until it reach a break. the loop will set delta E */
        x[c]=tgts[c];
        /* subtract contributions from each parameter for conformation c for each conformation 
          e.g deltaE - cos (dihedral * periodicity) * parameter generated from chromosomes 
          Therefore, it is delta E - sum of cosines for each dihedral */
        for(i=i0;i<j;i++,t++){
          x[c]-=Vs[i] * tset[t]; // tset is cos(n * dih)
#if DEBUG>2
          printf("scoreIt: Azp for c%d is %f with az index of %d and target is %f with V of %f\n",c,AZP[az],az,tset[t],Vs[i]);
#endif
          // this one uses V *(1 + cos(n*dih) but it is the same as above
          //x[c]-= (Vs[i] * (1 + tset[t]));
        }
        /* add differences in this error from all other errors */
        for(int c2=breaks[b];c2<c;c2++){
#if DEBUG>2
          printf("%d - %d\n",c,c2); //print the pairs index
#endif
          // calculate the absolute error for each pairs 
          float err=x[c]-x[c2];
          // sum the absolute of the errors (err) - -err = + err ; +err = +err
          //s+=(err<0.0f?-err:err); //ternary operator, condition is err < 0.0; if true err is negative, if false error is positive 
          s+=abs(err); 
        }
        /* next conformation */
        ++c;
      } 
      /* add little error to big error S, weighted by number of pairs, wt  is 2 / nconf*(nconf-1) */
      *S+=s*wts[b];
      /* go to next breakpoint (data set) */
      ++b;
    }
  } //end if in Psize
}

/**************************************************************************************************
*                                 | function 4: calcAreas |                                       *
*                                                                                                 *
*     calculates the areas (the probability) each individual has of mating                        *
*___________________________________Parameters____________________________________________________*
* @param scores scores for each individual (set of parameters)                                    *
* @param areas fitness for each individual, in terms of probability of mating                     *
* @param ptrs array of pointers from logical indices to actual indices into Vs for each individual*
* @param pSize number of individuals in the population                                            *
* @param genomeSize number of genes in a genome                                                   *
**************************************************************************************************/

__global__ void calcAreas(float *scores, float *areas, const int *ptrs, const int pSize, const int genomeSize) {
  int i=blockIdx.x * blockDim.x + threadIdx.x;
  //if((i<<1)<(pSize-1)*pSize){
  if(i<pSize){
    areas[ptrs[i]/genomeSize]=__expf(-scores[i]/scores[0]);
  }
}
/*****************************************************************************************
*                                | function 5: moveEm |
*
* @brief simple helper function for copying data from oldF, oldI to neWF, newI
*
* @param newF pointer to new float array
* @param newI pointer to new int array
* @param oldF pointer to old float array
* @param oldI pointer to old int array
* @param N number of floats/ints to copy
*****************************************************************************************/

__global__ void moveEm(float * newF, int *newI, float *oldF, int *oldI, int N) {
  int i=blockIdx.x * blockDim.x + threadIdx.x;
  if(i<N){
    newF[i]=oldF[i];
    newI[i]=oldI[i];
  }
}
/******************************| function 5 ends |***************************************/

/*****************************************************************************************
                   | sumEm and sumEmIndex : helper function for getSumAreas |

* @brief performs a sum of each successive pair of N numbers in source and stores the sums 
         in sums. intended to be run multiple times to sum over a whole array. if N is odd, 
         the last sum index will be N/2-1 and contain the sum of the last 3 numbers
*
* @param sums where to store the sums
* @param source where to get the numbers to sum together
* @param N the dimension of source
*
* @return                        ********************************************************/

__global__ void sumEm(float *sums, float *source, int N){
  int i=blockIdx.x*blockDim.x+threadIdx.x;
  int j=(i<<1);
  if(j+3<N)sums[i]=source[j]+source[j+1];
  else if(j+3==N) sums[i]=source[j]+source[j+1]+source[j+2];
  else if(j+2==N) sums[i]=source[j]+source[j+1];
}

/*
* @brief performs a sum of pairs of N numbers in source, using locations indicated 
         by pointers. pointers has indices multiplied by genomeSize. intended to be
         run multiple times to sum over a whole array. if N is odd, the last sum index 
         will be N/2-1 and contain the sum of the last 3 numbers
*
* @param sums where to store the sums
* @param source an array where to get the numbers to sum together
* @param N the dimension of source
* @param ptrs the indices to use when gathering pairs for summation
* @param genomeSize the number by which the indices in ptrs are scaled
*
* @return 
*/
__global__ void sumEmIndex(float *sums, float *source, int N, const int *ptrs, const int genomeSize){
  int i=blockIdx.x*blockDim.x+threadIdx.x;
  int j=(i<<1); // j = i*2 (mutiplication using a left bitwise shift)
  if(j+3<N)sums[i]=source[ptrs[j]/genomeSize]+source[ptrs[j+1]/genomeSize];
  else if(j+3==N) sums[i]=source[ptrs[j]/genomeSize]+source[ptrs[j+1]/genomeSize]+source[ptrs[j+2]/genomeSize];
  else if(j+2==N) sums[i]=source[ptrs[j]/genomeSize]+source[ptrs[j+1]/genomeSize];
#if DEBUG>1
  if(j+2<=N)printf(" %d:%f",i,sums[i]);
#endif
}
/*******************************| end of helper function |*******************************/
/*****************************************************************************************
*                                | function 6: getSumAreas |                             * 
*                        ---------uses sumEmIndex and sumEM--------                      *
*                                                                                        *
* @brief get sum of all areas                                                            *
* @param areas_d pointer to areas on device                                              *
* @param ptrs_d pointer to indices for each individual in population                     *
* @param pSize population size                                                           *
* @param temp_d pointer to temporary array on device                                     *
* @param genomeSize number of alleles in genome                                          *
*****************************************************************************************/

float *getSumAreas(float *areas_d, int *ptrs_d, int pSize, float *temp_d, const int & genomeSize){
  int dim=pSize; //Set dim to pSize
  int offset=0;

/* The triple chevron below describes an execution configuration the first
    argument(((dim>>1)+BLOCK_SIZE-1)/BLOCK_SIZE) in the execution configuration specifies 
    the  number of thread blocks in the grid, and the second specifies (BLOCK_SIZE) 
     the number of threads in a thread block */

  // return an array of sums (temp_d), sum up the probabilities in areas_d array
  sumEmIndex <<<((dim>>1)+BLOCK_SIZE-1)/BLOCK_SIZE, BLOCK_SIZE>>> (temp_d, areas_d, dim, ptrs_d, genomeSize);

#if DEBUG>1
  std::cout << std::endl;
#endif

  pSize >>= 1; // divide pSize by 2 using a right bitwise shift 
  while((dim>>=1)>1){  // while pSize/2 is greater than 1: Keep dividing (1/2 psize) by 2  
    offset^=pSize;  //bitwise XOR offest is 1/2 pSize then 0, then 1/2 pSize, then 0...
    // doing this switch the source to be (temp+pSize/2) then the source changes to (temp_d+0), then back and forth
    sumEm <<<((dim>>1)+BLOCK_SIZE-1)/BLOCK_SIZE, BLOCK_SIZE>>> (temp_d+offset, temp_d+(offset^pSize), dim);

#if DEBUG>1
  std::cout << std::endl;
#endif
  }
  return temp_d+offset;
}


/*
///////////////////////////////////////////////////////                 `
//////////////////////////////////                                       `
/////////////////////                                                  |   | 
/////////////                                                     ~ ~ ~ ~ ~ ~ ~
////////                                                         |              |
/////                                                        ____|              |____  
///                                                         |                        | 
//                                                       ___|          J.M           |___
/                                                       |              K.B               |
/                              PROGRAM BEGINS HERE      |              C.T               |
*****************************************************************************************/

/*****************************************************************************************
argc is a vairable with the number of arguments passed to GenA
argv is a vector of strings representing the the arguments the GenA takes
input file: parametersfitting data using the following format:
 _____________________________________________________________________        
|-<dihedral> <AMBER atom type for dihedral 1>                         |
|-<dihedral> <AMBER atom type for dihedral 2>                         |
|<name of data set> <dihedral 1> <dihedral 2>                         |
| <dihedral 1 value> <dihedral 2 value> <E_QM> <E_MM>                 |
| <dihedral 1 value> <dihedral 2 value> <E_QM> <E_MM>                 |
|                    ...                                              |
|/                                                                    | 
|<name of data set> <dihedral 1> <dihedral 2>                         |
| <dihedral 1 value> <dihedral 2 value> <E_QM> <E_MM>                 |
| <dihedral 1 value> <dihedral 2 value> <E_QM> <E_MM>                 |  
|                   ...                                               |
|/                                                                    |  
|_____________________________________________________________________|

<dihedral> is the name of dihedral e.g phi, psi, chi1, chi2, chi3, etc
<AMBER atom type for dihedral 1> e.g chi1 is N -CX-2C-2C for Met, get from frcmod file
<name of data set> is any name, e.g Metalpha, Metbeta, Metcharge
<dihedral 1 value> this is the dihedral value (deg) of the optimized QM structures 
     e.g 105.62
<E_QM> the QM energy of conformation i with restraint dihedral
<E_MM> the MM energy of conformation i with with zeroed dihedral parameters in the 
       frcmod
... repeat for all conformations within a break 
/ (refer to as break (brk))
a break seperate conformations that are different database
    e.g alpha backbone, beta backbone, charge amino acids
                                  GOODLUCK!!!
                                  [ O    O ]
                                  [    b ' ]
                                  [  ----- ]
contact: kellonbelfon@gmail.com with genA title for help
*****************************************************************************************/

int main(int argc, char *argv[]){

  /* start the timer */
  auto t1=std::chrono::high_resolution_clock::now();

  /*specify the string name of the savefile, scorefile, loadfile etc */
  std::string saveFile, loadFile, scoreFile, logFile, frcmodFile, inputFile, fitFile, azpFile;

  /* genetic algorithm parameters initiated */
  int pSize, nGen, rseed, peng, ncp, nCos, nChrom;
  float pMut, max, pCross, keep, nDataset;

  /* getting the filenames from the commands -r, -c, -s, -o, -f -y -a */
  for (int i=1;i<argc;i++){
    if(i+1<argc){
      if(argv[i][0]=='-'&&argv[i][1]=='r')saveFile=argv[++i];  //file that save amplitudes parameter (Vs)
      else if(argv[i][0]=='-'&&argv[i][1]=='c')loadFile=argv[++i]; //file with Vs for restart or from other forcefields
      else if(argv[i][0]=='-'&&argv[i][1]=='s')scoreFile=argv[++i]; // file that save the scores
      else if(argv[i][0]=='-'&&argv[i][1]=='f')frcmodFile=argv[++i]; //file that save frcmod file
      else if(argv[i][0]=='-'&&argv[i][1]=='o')logFile=argv[++i]; //file that save outputs 
      else if(argv[i][0]=='-'&&argv[i][1]=='i')inputFile=argv[++i]; // input file with dihedral info
      else if(argv[i][0]=='-'&&argv[i][1]=='y')fitFile=argv[++i]; // file with and idea of how your target energy change
      else if(argv[i][0]=='-'&&argv[i][1]=='a')azpFile=argv[++i]; // file with the periodicty that amplpitude will be zero
    }
  }

  /* open the output file which is the log file */
  std::ofstream logfile;
  logfile.open (logFile.c_str(), ios::out);
  
  /* open the score file to store scores */
  std::ofstream scorefile;
  scorefile.open (scoreFile.c_str(), ios::out); 
  scorefile << "#Generation" << std::setw(14) << "Chromosomes" << std::setw(12) << "Scores" << std::setw(14) << "areas\n";
 
  /* Now load genA parameters, from the parmfile -p  */
  for (int i=1;i<argc;i++){
    if(i+1<argc){
      if(argv[i][0]=='-'&&argv[i][1]=='p'){
      ConfigFile cfg(argv[++i]); //file that has the genetic algorithm parameters
  
      // check if keys exixt
      if (!(cfg.keyExists("pSize"))) std::cout << "pSize was not specified, using default of 2000\n";  
      if (!(cfg.keyExists("nGen"))) std::cout << "nGen was not specified, using default of 1000\n";  

      // Retreive the value of keys 
      pSize = cfg.getValueOfKey<int>("pSize", 2000);
      logfile << "Population Size (pSize): " << pSize << "\n\n";
      nGen = cfg.getValueOfKey<int>("nGen", 1000);
      logfile << "Number of Generations (nGen): " << nGen << "\n\n";
      pMut = cfg.getValueOfKey<float>("pMut", 1);
      logfile << "Probability of Mutations (pMut): " << pMut << "\n\n";
      max = cfg.getValueOfKey<float>("max", 0.5);
      logfile << "Maximal permissible mutation (max): " << max << "\n\n";
      pCross = cfg.getValueOfKey<float>("pCross", 0.8);
      logfile << "Probability of crossover (pCross): " << pCross << "\n\n";
      rseed = cfg.getValueOfKey<int>("rseed", 314245);
      logfile << "Random seed (rseed): " << rseed << "\n\n";
      peng  = cfg.getValueOfKey<int>("peng", 5);
      logfile << "Print scores every  " << peng << "generations (peng)\n\n";
      ncp  = cfg.getValueOfKey<int>("ncp", 1);
      logfile << "Print scores of only " << ncp << " chromosomes every peng \n\n";
      nCos = cfg.getValueOfKey<int>("nCos", 4);
      logfile << "Periodicity (nCos): " << nCos << "\n\n";
      keep = cfg.getValueOfKey<float>("keep", 0.1);
      logfile << "We will use " << keep << " for the elitist regime\n\n"; 
      nDataset = cfg.getValueOfKey<int>("nDataset", 1);
      logfile << "Number of Dataset is : " << nDataset << "\n\n";
      if(!loadFile.empty()) {
        nChrom = cfg.getValueOfKey<int>("nChrom", 1);
        logfile << "Number of chromosome reported is : " << nChrom << "\n\n";
        }
      }
    } 
  }

 
/* initializing GPU (_d) and CPU arrays */ 
  cudaError_t error;
  size_t nRands;
  curandGenerator_t gen;
  float *Vs, *Vs_d, *rands, *rands_d, *tset, *tset_d, *tgts, *tgts_d, *wts, *wts_d, *xx_d;
  float *AZP_d, *AZP, *scores, *scores_d, *areas, *areas_d;
  int genomeSize, g, *ptrs_d, *ptrs, N, nConf=0, *breaks, *breaks_d, nBreaks; 
  int save=pSize*keep; //save is number of chromosome we will keep as elitist


/***************************| load data from load.cpp |***********************************
*  check load.cpp for this section                                                       *
*  map is a way to create a dictionary, correction map is an array with key              * 
*****************************************************************************************/

/* initiating container with key and values name correctionMap */
  std::map<std::string,DihCorrection> correctionMap;

/* input file open, with dihedral info */ 
  std::ifstream inputfile;
  inputfile.open (inputFile.c_str(), std::ios::in);
 
/* load in arrays generated from load.cpp, check it out for further comments */
  load(inputfile, &tset, &tgts, &wts, &nConf, &breaks, &nBreaks, &genomeSize, 
       correctionMap, nCos);
  logfile << "Input file loaded ('_')" << "\n\n";
/****************************************************************************************/

/*************************| memory allocation |*******************************************
*   Declare and allocate host and device memory, copy data arrays from CPU host 
       (breaks,tset,                                                                 
*     tgts,wts) to device GPU (breaks_d, etc)                                        
*****************************************************************************************/

#if DEBUG && 0
  for(int i=0;i<nConf;i++){
    for(int j=0;j<genomeSize;j++)
      std::cerr << ' ' << tset[i*genomeSize+j];
    std::cerr << std::endl;
  }
  std::cerr << tgts[0] << ' ' << tgts[1] << ' ' << tgts[2] << ' ' << tgts[3] << std::endl;
  std::cerr << "first cudaMalloc, " << nBreaks << " breaks" << std::endl;
#endif

/* Allocate memory on GPU */
  cudaMalloc((void **)&breaks_d, nBreaks*sizeof(int));
  cudaMalloc((void **)&tgts_d, (nBreaks-1+nConf*(1+genomeSize))*sizeof(float));
  wts_d=tgts_d+nConf;
  tset_d=wts_d+nBreaks-1;

#if DEBUG
  std::cerr << "COPY" << std::endl;
#endif

/* Copying over the arrays from the CPU to GPU
nbreaks is the # of dataset + 1. e.g if you are doing alpha and beta backbone set then nbreaks=3
genomesize is the # of fitting dihedral * periodicity, e.g 3 set of dihedral * 4 periodicity = 12
nconf is the # of conformations you are fitting
tgts is (E_QMi-E_MMi) + (E_MMref-E_QMref) for each conformation, which = nconf, see load.cpp
tset is the cos(dih*periodicity) for 4 periodicity for a dihedral for each conformation
so 20 conf will give tgts of 20 (nconf) * 12 (# of dih * periodicity) = 120 
*/
  cudaMemcpy(breaks_d, breaks, nBreaks*sizeof(breaks[0]), cudaMemcpyHostToDevice);
  if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s\n", cudaGetErrorString(error));}
  cudaMemcpy(tset_d, tset, nConf*genomeSize*sizeof(float), cudaMemcpyHostToDevice);
  if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s\n", cudaGetErrorString(error));}
  cudaMemcpy(tgts_d, tgts, nConf*sizeof(float), cudaMemcpyHostToDevice);
  if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s\n", cudaGetErrorString(error));}
  cudaMemcpy(wts_d, wts, (nBreaks-1)*sizeof(*wts), cudaMemcpyHostToDevice);
  if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s\n", cudaGetErrorString(error));}

/**********************| initiate GPU blocks and # of random variable |*************************** 
*          we need randoms, new pop 3xcrossover, genomeSizexmut                                  *    
*        genome size is the number of genes which is all the parameters,                         *
*   e.g for 4 periodicity and three dihedral fitting, then genomesize will be 4 * 3 = 12         *
*   nRands is number of randoms we need for each set of parameters                               *
*   e.g if psize (population size) is 10, then number of random number we will need is           *
*                   (3+(# of periodicity x # of dihedral)) * psize                               *
* so for 4 periodicity and 3 dihedral fitting (chi1 chi2 chi3), then nRands = 3+12 * 10 = 150    *
*________________________________________________________________________________________________*  
*  nBlocks is dependent on the population size, it is use to figure out how many GPU blocks      *
*  we need to initialize the arrays for calculations. Each block has 256 threads.                *
*  one thread represent one individual (chromosome with soln parameters) from the population     *
*   e.g population size of 2000 will require (2000+256-1)/256 = 8.81 => 8 blocks                 *
*                                                                                                *
*************************************************************************************************/
  nRands=(3+genomeSize)*pSize;
  int nBlocks=(pSize+BLOCK_SIZE-1)/BLOCK_SIZE;

#ifdef DEBUG
  std::cerr << nRands << "nRands\n";
  std::cerr << nBlocks << " blocks\n";
#endif

/*******************************| initializing host and device variables|************************
*         N (bitwise operation below) is the pSize (1st input) multiply by 2;                   *
*       initiating the chromosomes  which have the solns                                        *
************************************************************************************************/
#if DEBUG
  printf("Allocate memory\n");
#endif

  rands=(float *)malloc(nRands*sizeof(float));
  N=(pSize<<1);
  HANDLE_ERROR(cudaMalloc((void **)&Vs_d, (N*(genomeSize+4)+pSize*nConf+nRands)*sizeof(float)));
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s\n", cudaGetErrorString(error));}
  rands_d=Vs_d+N*genomeSize;
  scores_d=rands_d+nRands;
  areas_d=scores_d+(N<<1);
  xx_d=areas_d+(N<<1);
  scores=(float *)malloc(sizeof(*scores)*N);
  float *scores_ds[2];
  scores_ds[0]=scores_d;
  scores_ds[1]=scores_d+N;

  // allocate memory to host Vs (amplitudes or barrier height for the cosine function)
  Vs=(float *)malloc(N*genomeSize*sizeof(float));
  areas=(float *)malloc(N*sizeof(float));
  /* allocate the memory space to hold array of pointers (prts) of size N (2*pSize)
  these pointers point to the individuals (chromosome) in the population */
  ptrs=(int *)malloc(sizeof(int)*N);
  ptrs[0]=0;
  for(g=1;g<N;g++)ptrs[g]=ptrs[g-1]+genomeSize;
  HANDLE_ERROR(cudaMalloc((void **)&ptrs_d, N*2*sizeof(int)));
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s\n", cudaGetErrorString(error));}
  int *ptrs_ds[2];
  ptrs_ds[0]=ptrs_d;
  ptrs_ds[1]=ptrs_d+N;
  cudaMemcpy(ptrs_d, ptrs, sizeof(int)*N, cudaMemcpyHostToDevice);
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s\n", cudaGetErrorString(error));}
  int curList=0;

#if 0
  HANDLE_ERROR(cudaMalloc((void **)&scores_d, N*sizeof(float)));

  HANDLE_ERROR(cudaMalloc((void **)&xx_d, nOffspring*nConf*sizeof(float)));
#endif

/* thrust is a c++ template library for CUDA similar to STL it have two containers: 
       thrust::host_vector<type> and thrust::device_vector<type>
  The containers make common operations such as cudaMalloc, cudaFree, cudaMemcpy, more concise
  e.g thrust::host_vector<int> vec_h(2) will allocate host vector with 2 elements
    thrust::device_vectore<int> vec_d = vec_h will copy host vector to device
  This will allow you to directly manipulate device values from the host
    so vec_d[0] = 5; can be done from host  and once you output vector memory is 
    automatically released 
   it have a few algorithms, we use thrust::sort(), */
  thrust::device_ptr<int> dPtrs(ptrs_d), dPtrs_save(ptrs_d+save);
  thrust::device_ptr<float> dScores(scores_d), dVs(Vs_d);
  thrust::device_ptr<float> dScores_save(scores_d+save),
                            dScores_pSize(scores_d+pSize),
                            dScores_N(scores_d+N);


/**************************| Create a random generator |********************************************
* curandCreateGenerator takes two parameters: pointer to generator (*gen), type of generator       *
Once created,random number generators can be defined using the general options seed, offset,& order*
When rng_type is CURAND_RNG_PSEUDO_DEFAULT, the type chosen is CURAND_RNG_PSEUDO_XORWOW            *
*__________________________________________________________________________________________________*
*curandSetPseudoRandomGeneratorSeed takes two parameters (1) the generator (gen) & (2) seed value  *
* seed value # is used to initialize the generator and control the set of random numbers;          *
* same seed will the give same set of random numbers of the psuedorandom generator                 *
* rseed is the random number specified from the 6th input)                                         *
*__________________________________________________________________________________________________*
*    curandGenerateNormal take 5 parameters:                                                       * 
*  (1) generator - Generator to use                                                                *
*  (2) outputPtr - Pointer to device memory to store CUDA-generated results,                       *
                or Pointer to host memory to store CPU-generated resluts                           *
*  (3) num - Number of floats to generate                                                          *
*  (4) mean - Mean of normal distribution                                                          *
*  (5) stddev - Standard deviation of normal distribution                                          *
* Results are 32-bit floating point values with mean and standard deviation.                       * 
***************************************************************************************************/

#if DEBUG
  printf("Create random generator\n");
#endif

  // create the generator name gen
  curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);

#if DEBUG
  printf("Seed random generator\n");
#endif
  // initiate the generator with the random seed (rseed)
  curandSetPseudoRandomGeneratorSeed(gen, rseed);
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s (seed)\n", cudaGetErrorString(error));}
#if DEBUG
   std::cerr << "GenerateNormal" << std::endl;
#endif
  // Vs_d is the amplitudes which is random numbers
  curandGenerateNormal(gen, Vs_d, N*genomeSize, 0, 1);
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s (normal)\n", cudaGetErrorString(error));}

#if DEBUG
  cudaMemcpy(Vs, Vs_d, sizeof(float)*genomeSize*N, cudaMemcpyDeviceToHost);
  /// print the three Vs from the first two chromosomes. 
  std::cout << "random Vs, created on GPU" << std::endl;
  for(int i=0;i<1;i++){
    std::cout <<  Vs[ptrs[i]] << " " << Vs[ptrs[i]+1] << " " << Vs[ptrs[i]+2] << std::endl;  
  }    
#endif

  /*****  if we have a load file copy Vs (amplitude parameters) from the loaded file and populate Vs ***********/
  if(!loadFile.empty()) {
    std::ifstream loadfile;
    loadfile.open (loadFile.c_str(), std::ios::in);
    // copy the random Vs to add previous chromosome of nChrom
    cudaMemcpy(Vs, Vs_d, sizeof(float)*genomeSize*N, cudaMemcpyDeviceToHost);
    if (loadfile.is_open()) {
      for (int i=0;i<nChrom;i++) {
        for (int j=0;j<genomeSize;j++) {
          loadfile >> Vs[ptrs[i]+j]; 
        }
      }
    }
    // print the three Vs from the first two chromosomes, to ensure your Vs were loaded. 
    logfile << "Here is your loaded Vs(amplitudes) for first two chromosomes: \n\n" << std::endl;
    for(int i=0;i<2;i++){
      for(int j=0;j<genomeSize;j++){
        logfile <<  Vs[ptrs[i]+j] << " ";  
      }
      logfile << "\n";
    }      

    // copy loaded Vs  to the GPU and overwrite random Vs. If user only create two chromosomes or 
    // previous Vs then the rest of the chromosome will be random 
    cudaMemcpy(Vs_d, Vs, N*genomeSize*sizeof(*Vs), cudaMemcpyHostToDevice);// copy to GPU 
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s\n", cudaGetErrorString(error));}
  }

#if DEBUG
  // check to see if Vs was transfer to gpu successful 
  /// print the three Vs from the first two chromosomes. 
  std::cout << "loaded Vs" << std::endl;
  for(int i=0;i<1;i++){
    std::cout <<  Vs[ptrs[i]] << " " << Vs[ptrs[i]+1] << " " << Vs[ptrs[i]+2] << std::endl;  
  }    
  cudaMemcpy(Vs_d, Vs, pSize*genomeSize*sizeof(*Vs), cudaMemcpyHostToDevice);// copy to GPU 
  cudaMemcpy(Vs, Vs_d, sizeof(float)*genomeSize*N, cudaMemcpyDeviceToHost); // copy back to CPU
  /// print the three Vs from the first two chromosomes. 
  std::cout << "After transfer of loaded Vs to GPU" << std::endl;
  for(int i=0;i<1;i++){
    std::cout <<  Vs[ptrs[i]] << " " << Vs[ptrs[i]+1] << " " << Vs[ptrs[i]+2] << std::endl;  
  }    
#endif
  
   /*  This section zero amplitudes that we want based on the AZP file  */
  // if we have AZP file we will read in the values and store it in AZP arrays
  // AZPfile has to have a number of elements equal to genomsize 
  if(!azpFile.empty()) {
    // allocate memory for AZP to zero amplitudes, This array is multiply by Vs   
    AZP=(float *)malloc(genomeSize*sizeof(float));
    std::ifstream azpfile;
    azpfile.open (azpFile.c_str(), std::ios::in);
    if (azpfile.is_open()) {
      for(int j=0;j<genomeSize;j++){
        azpfile >> AZP[j]; 
      }
    }
    HANDLE_ERROR(cudaMalloc((void **)&AZP_d, (genomeSize)*sizeof(float)));
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s\n", cudaGetErrorString(error));}
 
    // now apply AZP, first copy over Vs from GPU
    cudaMemcpy(Vs, Vs_d, sizeof(float)*genomeSize*N, cudaMemcpyDeviceToHost);
    for(int i=0;i<N;i++){
      int az = 0;
      for(int j=0;j<genomeSize;j++,az++){
        // Because some of the Vs is negative let us ensure Vs become 0. If az is 1 then nothing Vs remain
#ifdef DEBUG
        printf("Before Vs is (%f) and azp is %f\n", Vs[ptrs[i]+j], AZP[az]);
#endif
        // if AZP is 0 then the amplitude will be zero
        if (AZP[az] == 0.0){
          Vs[ptrs[i]+j] = 0.00000; 
        }
        // if AZP is negative then make Vs for that periodicity equal to the Vs plus the number
        // This is a way to ensure two Vs be the same value 
        if (AZP[az] < 0.0){
          int mv = j + (-AZP[az]);
          Vs[ptrs[i]+j] =  Vs[ptrs[i]+mv];
        }
#ifdef DEBUG
        printf("az index is %d, AZP is %f, mv is %d, i is %d, Vs[i] %f =  Vs[i+mv] %f\n", az,AZP[az],mv,i,Vs[ptrs[i]+j],Vs[ptrs[i]+mv]); 
        printf("After Vs is (%f)\n", Vs[ptrs[i]+j]);
        printf("Vs is %f for az%d of AZP value %f for genome %d and psize %d\n",  Vs[ptrs[i]+j],az,AZP[az],j,i);
#endif
      }
    }      
    // copy new Vs to GPU, need to clean up the code (too much copying, only copy Vs once) 
    cudaMemcpy(Vs_d, Vs, sizeof(float)*genomeSize*N, cudaMemcpyHostToDevice);// copy to GPU 
    cudaMemcpy(AZP_d, AZP, sizeof(float)*genomeSize, cudaMemcpyHostToDevice);// copy to GPU 
  } //end AZP section  

#ifdef DEBUG
  for(int j=0;j<genomeSize;j++){
    printf("aZp is %f\n",AZP[j]);
  }
#endif

 
/***************************| score of the first set of chromosomes |*******************************
* Here we score the two arrays of parents with solution parameters in the initial population       * 
*****************************************nm  *******************************************************/
#if DEBUG
    std::cerr << "1stscore" << std::endl;
#endif

    /* lauch first kernel to score the initial set of chromsomes (Vs_d) and output scores in scores_ds
      betweem the triple chervon is called the execution configuration that takes two parts
      1st part takes the number of thread blocks and the second part take the number of threads in a block */
    scoreIt <<<(N+BLOCK_SIZE-1)/BLOCK_SIZE, BLOCK_SIZE>>> (scores_ds[curList], areas_d, Vs_d, ptrs_ds[curList], tset_d, tgts_d, wts_d, breaks_d, nConf, pSize, genomeSize, xx_d);
    /* score of chromosomes out of psize since we initiated 2 times psize */
    scoreIt <<<(N+BLOCK_SIZE-1)/BLOCK_SIZE, BLOCK_SIZE>>> (scores_ds[curList]+pSize, areas_d, Vs_d, ptrs_ds[curList]+pSize, tset_d, tgts_d, wts_d, breaks_d, nConf, pSize, genomeSize, xx_d);
  
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s (1stscore)\n", cudaGetErrorString(error));}

#if DEBUG
    std::cerr << "1stsort" << std::endl;
#endif

    /* sort the scores from each chromosome of the initial population */
    thrust::sort_by_key(thrust::device_pointer_cast(scores_ds[curList]), thrust::device_pointer_cast(scores_ds[curList]+N), thrust::device_pointer_cast(ptrs_ds[curList]));
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s (1stsort)\n", cudaGetErrorString(error));}

    // print the initial scores based on ncp as -1 
    cudaMemcpy(scores, scores_ds[curList], sizeof(*scores)*ncp, cudaMemcpyDeviceToHost);
    for(int m=0;m<ncp;m++){
      scorefile << std::setw(6) << "-1" << std::setw(14) << m << std::setw(18) << scores[m]/nDataset << "\n";
      }


#if DEBUG>2
    cudaMemcpy(scores, scores_ds[curList], sizeof(*scores)*N, cudaMemcpyDeviceToHost);
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s\n", cudaGetErrorString(error));}
    cudaMemcpy(Vs, Vs_d, sizeof(*Vs)*N*genomeSize, cudaMemcpyDeviceToHost);
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s\n", cudaGetErrorString(error));}
    cudaMemcpy(ptrs, ptrs_ds[curList], sizeof(*ptrs)*N, cudaMemcpyDeviceToHost);
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s\n", cudaGetErrorString(error));}
       /* i is each chromosome, scores[i] is scores, Vs[ptrs[i]] is the amplitude parameters;
         Vs[ptrs[i]]+n specifies the next n amplitude. e.g chromosome i have genomesize amplitude parms 
         e.g  Vs[ptrs[i]]+1 is the amplitude term when the periodicity is 3 for the 1st dihedral being
        fitted, and  Vs[ptrs[i]]+4, the amplitude term when the periodicity is 4 for the 2nd dihedral */
    for(int i=0;i<N;i++){
      std::cerr << i << ": [" << ptrs[i] << "] = " << scores[i] << " {"<<Vs[ptrs[i]]<<" "<<Vs[ptrs[i]+1]<<" "<<Vs[ptrs[i]+2]<<" "<<Vs[ptrs[i]+3]<<"}\n";
    }
#endif

/****************************| Let us begin the iterations through generations |********************

 Genetic algorithm iterations through the number of generations (nGen: 2nd input) 

****************************************************************************************************/

  /* for loop for the generation */
  for(g=0;g<nGen;g++){
  
/*************************| Step1: Generate random numbers |****************************************/
#if DEBUG>1
  printf("Generate random numbers\n");
  printf(" %d",g);fflush(stdout);
#endif
    
    // create an array of random numbers (rands_d) used for mutations and crossover where the number of random #s is nRands 
    curandGenerateUniform(gen, rands_d, nRands);
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s\n", cudaGetErrorString(error));}

/***| Step2: calculate the probabilities (areas) each individual (chromosome) has of mating |******/
#if DEBUG>2
    std::cerr << "Mate" << std::endl;
#endif
    calcAreas <<<nBlocks, BLOCK_SIZE>>> (scores_ds[curList], areas_d, ptrs_d, pSize, genomeSize);

/***| Step3:  mate the individuals (chromosomes,Parent[0],[1]) selected for the next generation |***/
    mateIt <<<nBlocks, BLOCK_SIZE>>> (Vs_d, ptrs_ds[curList], areas_d, 
  getSumAreas(areas_d, ptrs_ds[curList], pSize, areas_d+N, genomeSize),
 rands_d, pCross, pSize, genomeSize);
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s (mate)\n", cudaGetErrorString(error));}

/*****************| Step4: mutate individuals generated after mating |*****************************/
#if DEBUG>2
    std::cerr << "Mutate" << std::endl;
#endif
    mutateIt <<<nBlocks, BLOCK_SIZE>>> (Vs_d, ptrs_ds[curList]+pSize, rands_d+pSize*3, pSize, pMut, max, genomeSize, AZP_d);
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s (mutate)\n", cudaGetErrorString(error));}

/**************| Step5: Score the individuals to select for the next generation |*******************/
#if DEBUG>2
    std::cerr << "Score" << std::endl;
#endif
    scoreIt <<<nBlocks, BLOCK_SIZE>>> (scores_ds[curList]+pSize, areas_d, Vs_d, ptrs_ds[curList]+pSize, tset_d, tgts_d, wts_d, breaks_d, nConf, pSize, genomeSize, xx_d);
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s (score)\n", cudaGetErrorString(error));}

#if DEBUG>2
    //std::cerr << "Display em:\n\tCopy scores" << std::endl;
    cudaMemcpy(scores, scores_ds[curList], sizeof(*scores)*N, cudaMemcpyDeviceToHost);
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s\n", cudaGetErrorString(error));}
    //std::cerr << "\tCopy Vs" << std::endl;
    cudaMemcpy(Vs, Vs_d, sizeof(*Vs)*N*genomeSize, cudaMemcpyDeviceToHost);
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s\n", cudaGetErrorString(error));}
    //std::cerr << "\tCopy ptrs" << std::endl;
    cudaMemcpy(ptrs, ptrs_ds[curList], sizeof(*ptrs)*N, cudaMemcpyDeviceToHost);
    if((error=cudaGetLastError())!=cudaSuccess){fprintf(stderr, "Cuda error: %s\n", cudaGetErrorString(error));}
    for(int i=0;i<N;i++){
      /* below you can print the scores for a chromosomes every generation */
      std::cout << "This is Generation: " << g << " and Chromosome (set of parameters): " << i << std::endl;
      std::cout << "Score: " << scores[i] << std::endl;
    /* below you can print out the scores and the first four barrier parameters,since we are using 
     4 periodicity, the first 4 barrier parameters are for the 1st dihedral in the input file */
     std::cerr << i << ": [" << ptrs[i] << "] = " << scores[i] << " {"<<Vs[ptrs[i]]<<" "<<Vs[ptrs[i]+1]<<" "<<Vs[ptrs[i]+2]<<" "<<Vs[ptrs[i]+3]<<"}\n";
    }
#endif

/*****| Step6: Sort the scored chromosomes (individuals) & select for mating for next generation |**/
#if DEBUG>2
    std::cerr << "Move 1" << std::endl;
#endif
    moveEm <<<(save+BLOCK_SIZE-1)/BLOCK_SIZE, BLOCK_SIZE>>> (scores_ds[curList^1], ptrs_ds[curList^1], scores_ds[curList], ptrs_ds[curList], save);
#if DEBUG>2
    std::cerr << "Move 2" << std::endl;
#endif
    moveEm <<<(pSize+BLOCK_SIZE-1)/BLOCK_SIZE, BLOCK_SIZE>>> (scores_ds[curList^1]+save, ptrs_ds[curList^1]+save, scores_ds[curList]+pSize, ptrs_ds[curList]+pSize, pSize);//nOffspring);
#if DEBUG>2
    std::cerr << "Move 3" << std::endl;
#endif
    moveEm <<<(pSize-save+BLOCK_SIZE-1)/BLOCK_SIZE, BLOCK_SIZE>>> (scores_ds[curList^1]+save+pSize, ptrs_ds[curList^1]+save+pSize, scores_ds[curList]+save, ptrs_ds[curList]+save, pSize-save);
    curList^=1;

    /* first sort only the non elitist with the offspring */
#if DEBUG>1
    std::cerr << "Selection sort (" << N << " items, less " << save << ")" << std::endl;
#endif
    thrust::sort_by_key(thrust::device_pointer_cast(scores_ds[curList]+save), thrust::device_pointer_cast(scores_ds[curList]+pSize+save), thrust::device_pointer_cast(ptrs_ds[curList]+save));

    /* second sort  is to order the first sort within psize by fitness */
#if DEBUG>1
    std::cerr << "Rank sort" << std::endl;
#endif
    thrust::sort_by_key(thrust::device_pointer_cast(scores_ds[curList]), thrust::device_pointer_cast(scores_ds[curList]+pSize), thrust::device_pointer_cast(ptrs_ds[curList]));

/****************************************************************************************************
* Here you can print the score of chromosomes (total is 2 x population size) to score file (-s)     *
****************************************************************************************************/
    //peng --> print every n generation
    //ncp --> number of chromosomes to print
    //if generation is divisable by peng
    if(g%peng==0) {
      //scorefile << "#Generation" << std::setw(14) << "Chromosomes" << std::setw(12) << "Scores\n";
      cudaMemcpy(scores, scores_ds[curList], sizeof(*scores)*ncp, cudaMemcpyDeviceToHost); //copy over ncp scores
      cudaMemcpy(areas, areas_d, sizeof(*areas)*ncp, cudaMemcpyDeviceToHost); //copy over ncp areas
      // divide score by the number of datasets to print the average of the datasets since score is sum of each dataset score
      for(int m=0;m<ncp;m++){
        scorefile << std::setw(6) << g << std::setw(14) << m << std::setw(18) << scores[m]/nDataset << std::setw(18) << areas[m] << "\n";
      }
    }

/* END GENETIC ALGORITM */
  } 

  scorefile.close();
/****************************************************************************************************
*    TERMINATION, LAST RESULTS < SCORES AND PARAMETERS FOR EACH INDIVIDUAL
****************************************************************************************************/
 
/***************************************************************************************************/
  /*  copy over the results from GPU to the CPU to save the scores and parameters */
  cudaMemcpy(Vs, Vs_d, sizeof(float)*genomeSize*N, cudaMemcpyDeviceToHost);
  cudaMemcpy(ptrs, ptrs_ds[curList], sizeof(int)*N, cudaMemcpyDeviceToHost);
  cudaMemcpy(scores, scores_ds[curList], sizeof(float)*N, cudaMemcpyDeviceToHost);
  cudaMemcpy(tgts, tgts_d, sizeof(float)*nConf, cudaMemcpyDeviceToHost);
  cudaMemcpy(tset, tset_d, nConf*genomeSize*sizeof(float), cudaMemcpyDeviceToHost);

/****************************************************************************************************/

// Here I am writing out the initial dE and the final dE, see load.cpp for description 
  /* file that stores initial dE */
  std::ofstream fitfile;
  fitfile.open (fitFile.c_str(), ios::out); 
  fitfile << "#dE = (E (QMi) - E (QMref) ) - ( E (MMi) - E (MMref) \n\n";
  fitfile << "#ref is the first structure in each dataset  \n\n";
  fitfile << "#Below is the target relative energy and the best relative energy after the fit \n\n";
  fitfile << "#Conformations" << "      "<< "dE" <<  "      " << "dE_from_fit " <<   "      " << "dE - dE_fromfit" << "\n";
  
  /*  save the last dE  */
  float new_dih[nConf]; // hold the new dihedral(s) contribution 
  for(int i=0;i<nConf;i++){ //loop through the conformations 
    float sum_cos = 0.0f; //array that hold the sum the cosine for each conformation 
    for(int j=0;j<genomeSize;j++){
       // Vs[ptrs[0]+j] this is indexing through best chromosome amplitudes (Vs[ptrs[0]])
       // multiply the amplitude (Vs) by (1+ cos(n*dih) which is (tset[i*genomeSize+j] )
      //printf("Vs : %f\n", Vs[ptrs[i]+j]);
      sum_cos += (Vs[ptrs[0]+j] * (1 + tset[i*genomeSize+j]));
      //printf("cos(n*dih) : %f\n", tset[i*genomeSize+j]);
    }
    new_dih[i] = sum_cos; //Effnew new dihedral energy for each conf
  }

  int b = 0; // b is number of dataset
  int c = 0; // c is conformations
  while(c<nConf){ // check through the conformations 
    float conf0 =  new_dih[breaks[b]]; // conf0 is the 1st conf in the dataset (separated by breaks)
    while(c<breaks[b+1]){ // loop through the dataset
      // printf("1st conf: %f\n", conf0);
      new_dih[c] -= conf0; //subtract to 1st conf to get relative energy
      // printf("newdih: %f\n", new_dih[c]);
      ++c; // go to next conformation
      }
  ++b; // go to next data set
  }  
  for(int i=0;i<nConf;i++){ // loop to save to the fit file
    float temp = tgts[i] - new_dih[i];
    fitfile << std::setw(6) << i << std::setw(18) << tgts[i] << std::setw(16) << new_dih[i] << std::setw(18) << temp << "\n";
  }
  fitfile.close();
/****************************************************************************************************/

  /* saving all of the scores, with dihedral parameters to the logfile */
  logfile << "\n";
  logfile << "Printing all of the final dihedral parameters, check your -f file for the best one \n\n";
  logfile << "The first one is the best score, best parameters\n\n";
  /* loop through the population */
  for(int i=0;i<pSize;i++){
    /* these are the final scores for each individual in the population, print in the output file  */
    // divide score by the number of datasets to print the average of the datasets since score is sum of each dataset score
    logfile << std::fixed << scores[i]/nDataset << std::endl;
    for(std::map<std::string,DihCorrection>::iterator it=correctionMap.begin(); it!=correctionMap.end(); ++it){
    /* second.setGenome(Vs+ptrs[i]) is the dihedral parameters for each individual in the population 
      print in the output file                                                                     */
      logfile << it->second.setGenome(Vs+ptrs[i]);
    }
  }

/****************************************************************************************************/
  /* Save a frcmod file to use in Amber */
  if(!frcmodFile.empty()){
    std::ofstream frcmodfile;
    frcmodfile.open (frcmodFile.c_str(), ios::out);
    frcmodfile << "frcmod from GenA.cu \n";
    frcmodfile << "DIHE\n";
    for(std::map<std::string,DihCorrection>::iterator it=correctionMap.begin(); it!=correctionMap.end(); ++it){
      frcmodfile << it->second.setGenome(Vs+ptrs[0]); //the best parameters 
      }
    frcmodfile.close();
  }

/****************************************************************************************************/
  /* Save the amplitudes to a restart file  */
  if(!saveFile.empty()){
    std::ofstream savefile;
    savefile.open (saveFile.c_str(), ios::out);
    for(int i=0;i<N;i++){
      for(int j=0;j<genomeSize;j++){
        savefile << std::setw(9) << Vs[ptrs[i]+j] << " ";
      }
      savefile <<"\n";
    }
  savefile.close();
  }

/****************************************************************************************************/
#if 0
  std::cout << scores[pSize] << std::endl;
  for(std::map<std::string,DihCorrection>::iterator it=correctionMap.begin(); it!=correctionMap.end(); ++it){
    std::cout << it->second.setGenome(Vs+ptrs[pSize]);
    //std::cout << it->second;
  }
#endif


#if 0
  printf("Copy random numbers\n");
  cudaMemcpy(rands, rands_d, nRands*sizeof(unsigned int), cudaMemcpyDeviceToHost);
  printf("Print random numbers\n");
  printf("%d",rands[0]);
  for(i=1;i<nRands;i++){
    printf(" %d",rands[i]);
  }
  putchar('\n');
#endif

  auto t2=std::chrono::high_resolution_clock::now();

  logfile <<"\n\n";
  logfile << "GenA took " 
          << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() 
          << " milli seconds to obtain your parameters" << "\n";

  logfile.close(); //close log file

/*****************| Free up Memory |*******************************************************/
  free(ptrs);
  curandDestroyGenerator(gen);
  //cudaFree(xx_d);
  cudaFree(Vs_d);
  cudaFree(ptrs_d);
  cudaFree(breaks_d);
  cudaFree(tgts_d);
  free(Vs);
  free(scores);
  //cudaFree(rands_d);
  free(rands);
  return 0;
}
