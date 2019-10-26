#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <list>
#include <map>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>

/* dihedral parameters in AMBER are:barrier (amplitudes), phase and periodicity
the format for these parameters in AMBER are:
atom_type        divider barrier_term  phase       periodicity
X -CA-CC-X          1       2.800       180.0         -2.000
The phase can be 0, meaning that a maximum energy is encountered at zero degrees.
Or it can be 180 degrees, meaning there is a minimum at 180 degrees
Equation: Vbarrier/divider × [1+cos(periodicity×φ − phase)]
A "negative" periodicity (-2 in the case above), tells the programs reading the parameter
file that additional terms are present for that particular connectivity
private only class DihCorrection is aware and public is where everything
has access so label and genome is accessible from all other function or class */

class DihCorrection {
private:
  /* maps periodicity to column
   label, is the string from input file denoting the dihedral name(e.g chi1, chi2, etc)
  */
  std::map<int,int> corrs;
  std::map<int,int> periodicity;
  std::map<int,int> fitgroup;
  std::string label;  //dihedral name
  const float *genome;
  
public:
  DihCorrection() {}
  DihCorrection(std::string label) : label(label) {
  }

  /* n is the key of the periodicity (e.g 0 1 2 3)*/
  void addFitgroup(int f, int indx) {
    fitgroup[f]=indx;
  }

  int fitgrpindx;
  void fitgroupindex(int x){
    fitgrpindx=x;
  }
  /*col is the value of (periodicity x # of dihedral) we are fitting
  so for three dihedral (chi1,2,3) then col is index 0 to 11 */
  void addCorr(int n, int col) {
    corrs[n]=col;
  }

  void addPer(int pindx, int per) {
    periodicity[pindx]=per;
  }


  std::map<int,int> getPer() { return periodicity; }
  std::map<int,int> getCorrs() { return corrs; }
  std::map<int,int> getFitgroup() { return fitgroup; }
  const DihCorrection & setGenome(const float *genome) {
    this->genome=genome;
    return *this;
  }
  friend std::ostream & operator<<(std::ostream &os, const DihCorrection &dc);
}; //END of DihCorrection classs 


// This is the output frcmod format 
std::ostream & operator<<(std::ostream &os, const DihCorrection &dc) {
  /* iterating through dc.corrs,  */
  // for loop that iterate through the number of periodicity for a given atom type 
  for(std::map<int,int>::const_iterator it=dc.corrs.begin(), next=it;it!=dc.corrs.end();it=next){
    ++next;
    os << dc.label << std::fixed; // dc.label is the dihedral atomtypes 
    if(dc.label.length()<=11)os << std::setw(8) << 1; //atomtypes and divider of 1 
    float val=dc.genome[it->second]; //This is the Vs fitted 
    /* it->first is the key in dc.corrs() */
    os << std::setw(12) << (val<0?-val:val) << std::setw(8) << (val<0?180:0) << std::setw(5) << (next==dc.corrs.end()?it->first:-it->first) << std::endl;
  }
  return os;
}

/* load function */
void load(std::istream & in, float **tset, int **ptrsV, int **ptrsT, int **ptrsD, int **allFgsinDsi, int **nDperD, float **tgts, float **wts, int *nConfs, int *nDataset, int **brks, int *nBrks, int *trainingSize, int *genomeSize, std::map<std::string,DihCorrection> & correctionMap, int **nVperFgi, int **nCosperFgperDs, int nCos, int nFg, int nDih, int *totdih, int **DihFgindx)
{
  
  int ch;
  int col=0;
  std::vector<int> nVperFg; // Vector that holds the needed Vs for a given Fg (not per Dataset)
  std::vector<int> checkV;
  std::vector<int> holdFg4dih;
  checkV.push_back(-99); // Added to do the checks, I remove the -99 after 
  // Let us read the dihedral name, atomtype read the lines only if it starts with + or -
  // populating the dc datatype
  while(in.peek()=='+'||in.peek()=='-'){
    int pindx=0;
    std::string label;
    std::string dih;
    ch=in.get();
    int strLen=(ch=='-'?11:35);
    in >> label;
#if LOAD_DEBUG
    std::cout << "Label is " << label  << "" << std::endl;
#endif
    in >> std::ws;
    for(int i=0;i<strLen;i++){
      dih.push_back(in.get());
    }
    /* dih is the atom types for the dihedral (e.g CX-2C-2C-S ) 
    label is the dihedral name (e.g chi2) so Dihedral[chi2]=CX-2C-2C-S */
#if LOAD_DEBUG
    std::cout << "Dihedral[" << label  << "]=" << dih << std::endl;
#endif

    DihCorrection dc(dih);
    // Now get the fitting groups and periodicity if defined
    while((ch=in.get())==' ');
    // if lines end and there is no dihedralfitting group info then use nCos from the genA parmfile 
    if(ch=='\n'){
      // nCos is the number of cosine (periodicity) 
      for(int n=nCos;n>0;n--) {
      /* print out n:col which is 4:0, 3:1, ....0;11 for three dihedral, see above */
#if LOAD_DEBUG
          std::cout << n << ": " << col << std::endl;
#endif
          dc.addCorr(n, col++);
          dc.addPer(pindx, n);
          pindx++;
      }
#if LOAD_DEBUG
      std::cout << "getCorrs().size = " << dc.getCorrs().size() << std::endl;
#endif
    } else {
      // While character is not '\n'
      while(ch!='\n'){
        //if(ch<'0'||ch>'9'){
        // IF character is negative sign then this is a fitting group index
        if(ch=='-'){
          std::string Fgind;
          // can only two fittig groups with 2 digits so from 0 to 99
          for(int i=0;i<2;i++){
            Fgind.push_back(in.get());
          //dih.push_back(in.get());
          } 
          // Get the fitting group index from the genA_input file it is the  negative number after atomtypes  
          //std::cout << "Fg index: " << Fgind << std::endl;
          dc.addFitgroup(0, stoi(Fgind));
          dc.fitgroupindex(stoi(Fgind));
          holdFg4dih.push_back(stoi(Fgind));
        // ELSe it is not a fittig group index but periodicity 
        } else {
#if LOAD_DEBUG
           std::cout <<  "ELSE ch > 0 or ch < 9 "<< std::endl;
	   std::cout << "Periodicity n " << ch-'0' << " col: " << col << std::endl;
#endif
           // Get the periodicity, positive numbers after the fitting group index, convert character to integer with -'0'
           dc.addPer(pindx, ch-'0');
           dc.addCorr(ch-'0', col++);
#if LOAD_DEBUG
	   std::cout << "Periodicity n " << ch-'0' << " col: " << col << " pindx " << pindx << std::endl;
           std::cout << "getCorrs().size = " << dc.getCorrs().size() << std::endl;
#endif
           pindx++;
        }
       do ch=in.get();while(ch==' ');
      }
    } // end of the atomtype lines in the input 
    correctionMap[label]=dc;
    if (checkV.back() != (dc.getFitgroup()[0]) ) {
      nVperFg.push_back(pindx);
      checkV.push_back(dc.getFitgroup()[0]);
#if LOAD_DEBUG
      std::cout << "nVperFg = " << pindx << " checkV = " << checkV.back() << " Fg index = " << dc.getFitgroup()[0] << std::endl;
#endif
    }
#if LOAD_DEBUG
    //std::map<std::string,DihCorrection>::iterator it=correctionMap.find(label);
    //if(it!=)
#endif
  } // end of while loop, Finish with the first part of the GenA_input file 

 
  /*  Moving to the second part of the file (dihedral information) */
  *trainingSize=col;  //This is nDih * nCos or the total # of Cosines
  std::vector<float> fgindx; // Vector that holds the number of dihedreal per Fg
  std::vector<int> nCosperFg; // Vector that holds the number of cosines per Fg
  std::vector<int> nDihperDs; // Vector that holds the number of dihedrals in a dataset
  std::vector<std::vector<float> > data; // holds dihedral data 
  std::vector<std::vector<int> > tsetptDs; // holds index into tset for a given dataset
  std::vector<int> breaks;    // holds the number of conformation per dataset 
  std::vector<int> equaldihs; // this will hold the index into tset for a given set of Vs 
  std::vector<float> weights; // this will hold the weights for a give dataset 
  std::vector<int> FgsinDs; //This vector holds only Fitting groups index in each dataset
  std::vector<int> allFgsinDs; //This vector holds all the Fitting groups index in each dataset
  std::vector<std::vector<int> > Fitgrpindx;  // this is a 2D array that holds the index for the fitting grps 
  std::vector<int> nFgperDs; //this will hold the total type of fitting group in each dataset 
  *nConfs=0;  // will enventuall hold the number of conformations  
  *nDataset=0;
  std::string line;
  double off; // use to calculate dE
  int stop_fill=0;
  int kill=-99;
  int nFgsit;
  FgsinDs.push_back(-99); // Added to do the checks -99 is then removed 

  // initialize row in the Fitgrpindx vector
  for(int tf = 0; tf < nFg; tf++) {  
    std::vector<int> rowF;
    Fitgrpindx.push_back(rowF);
  }

  while(in.good()&&std::getline(in,line)){
    /* First assign nConfs to breaks, then will loop back and place the nConf for each dataset 
     e.g if I have 2 dataset of 10 and 13 confs then breaks = {0,10,13} */
    breaks.push_back(*nConfs);
    std::list<DihCorrection*> cols; // will hold the dihedral columns 
    std::string label;
    std::istringstream input(line);
    input >> label; // This label is the dataset name 
#if LOAD_DEBUG
    std::cout << "DatasetName=" << label << std::endl;
#endif
    weights.push_back(1.0f); // Add 1.0 to the weight array , it will be used to divide 2/N(N-1) 
    // if the input is good 
    if(input.good()){
      // Check if the 1st label after the dataset name has a < sign in front the name 
      input >> label;
      /* atof converts string to double */
      // If the label has a < signs then the number that comes after is used  as a weight for the dataset 
      // Example <10 will give a weight of 10 * (2/(Nconf(Nconf-1))) for each conformation in the dataset
      if(label[0]=='<') {
        weights.back()=atof(label.data()+1);
#if LOAD_DEBUG
        std::cout << "weights is " << atof(label.data()+1) << " * 2/(N(N-1))" <<std::endl;
#endif 
      }
      // else save the label in correctionMap
      else { 
        cols.push_back(&correctionMap[label]);
      }
      // continue reading the next label (dihedral name)
      while(input.good()){
        input >> label;
        if(label[0]=='>') {
          /* get the number of different fitting group in this dataset */
          nDihperDs.push_back(atoi(label.data()+1));
#if LOAD_DEBUG
          std::cout << "nDihperDs = " << atoi(label.data()+1) << " for Dataset" << *nDataset <<std::endl;
#endif 
        } 
        else {
          cols.push_back(&correctionMap[label]);
        }
      }
    } // END of if(input.good())

    // Now we have the dihedral names (label) let us get the dihedrals and some important arrays 
    std::vector<float> dataRow;
    std::vector<int> tsetpt;
    int dihcounter=-1;
    /* Loop through each dataset */
    // forwardslash is the end of a given dataset
    
    /* flag is set to one at the begining of each dataset and set to 0 after the first conformation (row) in a dataset is read */
    int flag = 1;

    while(std::getline(in,line)&&line[0]!='/'){  // Loop through the datasets
      input.clear();
      input.str(line);
      /* assignment of the dataRow to be (1 + size of the genome (genome is number of dihedral * # of cosines)
       so for 3 dihedrals with 4 cosines it will be 13 to 0 */
      dataRow.assign(1+*trainingSize, 0);
      tsetpt.assign(*trainingSize, 0); //Using trainingSize here but don't need this amount 
      int tcount=0;
      double dih;
      // loop through columns 
      int fcounter=0;
      int stop_fill1=0; 
      /* loop through the rows of dihedrals for a given dataset */
      for(std::list<DihCorrection*>::iterator it=cols.begin();it!=cols.end();++it){
        input >> dih; //dihedral values in the input file dih1 dih2 ...
#if LOAD_DEBUG
        /* dih is the dihedral values located in input value */
        std::cout << "Dihedral Value (deg): " << dih;
#endif
       /* convert degrees from input file (-180 to 180) to radians(-pi to pi) */
        dih*=3.141592653589793238/180.0;
#if LOAD_DEBUG
        std::cout << " (rad): " << dih << std::endl;
        std::cout << "number of cosine for dihedral " << (*it)->getCorrs().size() << std::endl;
#endif
        // Get nCos from the size of the dihedral [This is the number of periodicity for a dihedral  
        nCos = (*it)->getCorrs().size();

        // Only want the number of cosines for a given fitting group
        // so if stop_fill is less than number of dihedrals and the fitting group id matches stop_fill 
        if (flag == 1){ // only looking at the first row in a given dataset to populate these vectors 
          stop_fill = (*it)->getFitgroup()[0]; //get the fitting group of this dihedral 
          //std::cout << "stop_fill = " << stop_fill << std::endl;
          allFgsinDs.push_back(stop_fill);   // populate the fitting group into allFgsinDs
          //std::cout << "ALL Fgs index: " << allFgsinDs.back() << std::endl;
          //std::cout << "Last element of FgsinDs: " <<FgsinDs.back() << std::endl; 
          if (*nDataset != kill) { // Only true when we move to the next dataset
            //nFgperDs.push_back(nFgsit); 
            //std::cout << "nFgsit " << nFgsit << "kill " << kill << std::endl; 
            nFgsit=0;
          } 
          if (FgsinDs.back() == stop_fill) { // the last Fg in the previous dataset is equal to the first Fg in this dataset
	    //std::cout << "FgsinDs.back() == stop_fill" << std::endl;
            if (*nDataset != kill) {  
	      //std::cout << "kill " << kill << " nCos= " << nCos << std::endl;
              nCosperFg.push_back(nCos);
              FgsinDs.push_back(stop_fill);
              nFgsit++;
              kill=*nDataset;
	      //std::cout << "Populate FginDs with Fg = " << stop_fill <<" nFgsit = " << nFgsit << std::endl;
            }
          }
          if (FgsinDs.back() != stop_fill) { //Only count the types of fitting group)
            nCosperFg.push_back(nCos);
            FgsinDs.push_back(stop_fill);
            //std::cout << "Populate nCosperFg with nCos = " << nCos << std::endl;
            nFgsit++;
            kill=*nDataset;
	    //std::cout << "Populate FginDs with Fg = " << stop_fill <<" nFgsit = " << nFgsit << " nCos = " << nCos << std::endl;
          } 
        }
        // only want the number of dihedral in a fitting group 
        if (stop_fill1 < nDih) {
          for (int f=0;f<nFg;f++) {
            if ((*it)->getFitgroup()[0] == f) {
              fgindx.push_back(f);
              stop_fill1++;
            }
          }
        }
        //std::cout << "number of cosine for dihedral (nCos) \n" << nCos << std::endl;

        /* nCos is the number of periodicity for a given dihedral */
        //for(int n=nCos;n>0;--n){
        for(int pt=0;pt<nCos;pt++){
       /* *it is the pointer to the correction map 
        n[col] = cos(dihedral* periodicity), where n is the periodicty and 
        col is an index for the periodicity, and dihedral in radians; 
        so for a given dihedral the cosine term is evaluated for the number of
        periodicity (n), so n cos terms for a given dihedral*/
          int n=(*it)->getPer()[pt]; // Get n value 

#if LOAD_DEBUG
          //std::cout << "col: (n value) " << (*it)->getCorrs()[n] << " " << std::endl;
          //std::cout << "fitGroup for label " << label << ": " << (*it)->getFitgroup()[0] << "  " << std::endl;

          //std::cout << "dih= " << dih << std::endl; 
          std::cout << "n[col] : " << n << "[" << (*it)->getCorrs()[n] << "]+=" << cos(dih*(float)n) << std::endl;
#endif 
          if (pt == 0) { //only want the largest ns to get the col index into tset
            tsetpt[tcount]=((*it)->getCorrs()[n]); 
            tcount++;
          }
          dihcounter++;
          //std::cout << "dihCounter :" << dihcounter << " trainingSize= " << trainingSize << std::endl;
          if (dihcounter < *trainingSize) {   
            for (int f=0;f<nFg;f++) {
                //std::cout << "NFG PART f getfitGroup for label " << label << ": " << (*it)->getFitgroup()[0] << " == " << f<< std::endl;
               //if f is equal to the fitting group number in the genA_input flag then populate the array 
              if ( (*it)->getFitgroup()[0] == f){
                // Only want this based on the first set ofdihedrl 
                //std::cout << "F loop " << fcounter << std::endl;
                Fitgrpindx[f].push_back(fcounter);
                fcounter++;
              }
            }
          }
          dataRow[(*it)->getCorrs()[n]]+=cos(dih*(double)n);
        } // END of ncos for loop 
#if LOAD_DEBUG
          /* (*it)->getCorrs().size() return the number of elements in *it which is 4 */
        std::cout << ' ' << (*it)->getCorrs().size() << std::endl;
#endif
      } // End of for loop that loops through a given row of dihedrals, now put the dE in the last cell of this row 
      flag = 0;
      /*  This section calculates deltaE btwn QM and MM for a given structure using the first 
          structure in the break(sets) as reference, E_QM is qm energies, E_MM is mm energies
          Only done to make the numbers smaller, pairwise stuff is done on the GPU
          off = E (MM ref) - E (QM ref), where ref is the first QM & MM energy within the break
          delta E = 0 for the first structure 
          delta E = (E (QMi) - E (MMi) ) + off where i is every structure after the break
          delta E = (E (QMi) - E (MMi) ) + ( E (MMref) - E (QMref) )
          delta E = (E (QMi) - E (QMref) ) - ( E (MMi) - E (MMref) )    */
      double E_QM,E_MM,dE;
      /* input E_QM and E_MM, from the inputstream take column E (QM energies) and E0 (MMenergies) */
      input >> E_QM >> E_MM;
      //std::cout << E_QM << E_MM << std::endl;
      
      /* the if statement evaluate the first conf after the break, else statement everything else */
      if(*nConfs==breaks.back()){
        off=E_MM-E_QM; // Emmref - Eqmref
        std::cout.precision(14);
        std::cout << std::setw(14) << "off: " << off << std::endl;
        std::cout.precision(14);
        std::cout << std::setw(14) << "EQM 0: " << E_QM << std::endl;
        std::cout << "EMM 0: " << E_MM << std::endl;
        dE=0; // set the reference pair to zero
      }else{
        dE=E_QM-E_MM+off;
      }
      /*  *trainingSize in dataRow = delta E */
#if LOAD_DEBUG
      std::cout << dE << std::endl;
#endif 
      // make the last column in trainingSize be the target (dE)
      dataRow[*trainingSize]=(float)dE;

#if LOAD_DEBUG
      std::cout << " deltaE="<<dataRow[*trainingSize]<<std::endl;
#endif

      /* for every conformation, push_back will add dataRow to the array data
         dataRow contains the deltaE and the cosine term  periodicity per dihedral */
      ++*nConfs;
      data.push_back(dataRow); //populate the row of cos(n*dih) for a given conformation 
      tsetptDs.push_back(tsetpt);// populate the row of tset index for a given conformation 
    }

    /* weights here is 1 / (nconf * ((nconf-1)/2)) for a given population in a break 
       similar to 2/(nconf)*(nconf-1) in 14SB paper when a weight value is added then 
       multiply 2/(nconf)*(nconf-1) * weight given in the input file  */
    weights.back()/=(float)((*nConfs-breaks.back())*(*nConfs-breaks.back()-1)/2);
    //std::cout << "END of Dataset " << *nDataset << "\n\n";
    ++*nDataset; // iterate to next dataset
    nFgperDs.push_back(nFgsit);
    //std::cout << "nFgsit " << nFgsit << "kill " << kill << std::endl; 
 
  } // End of data reading, now let us build the arrays and pointers for GenA RAGTAG
  breaks.push_back(*nConfs);
  *nBrks=breaks.size();

#if LOAD_DEBUG
  /* print the number of conformations (structures) and breaks */
  std::cout << *nConfs << " confs and " << *nBrks << " breaks" << std::endl;
#endif

  /* sizeof(x), number of bytes to represent type x 
     *tgts is the delta E for a given conformation (structure) 
     *tset is the cos(periodicity * dihedrals) terms for a given angle and periodicity
       e.g three dihedral with 4 periodicity will give 12 tset for a give structure that
           has the three dihedrals
      The idea if you have a structure that have a delta E of 0.4 kcal/mol. Then using 
      amber energy Vbarrier/divider × [1+cos(periodicity×dihedral − phase)]
      parameters Vbarrier and phase will be searched with the genA.cu to give a energy
      (sum of zeroed dihedral) as close to the target energy of 0.4 kcal/mol 

      *brks is the conformation position for the program to identify the breaks
        e.g 20 conf, 10 alpha, 10 beta, will have 3 breaks (break seperate alpha and beta)
            so break will be 0, 10, 20
      *wts is 1 / nconf * (nconf-1)/2 for a given population in a break 
        e.g 20 conf, 10 alpha, 10 beta will have two dataset alpha and beta
            so for alpha wts = 1/(10)*((10-1)/2) = 1/45 = 0.0222
            so for beta wts = 1/(10)*((10-1)/2) = 1/45 = 0.0222 */


  // Allocation of some arrays we need, see explanation above  
  *brks=(int *)malloc(sizeof(**brks)*breaks.size());
  *wts=(float *)malloc(sizeof(**wts)*weights.size());
  *tset=(float *)malloc(sizeof(**tset)* *nConfs* *trainingSize);
  *tgts=(float *)malloc(sizeof(**tgts)* *nConfs);
  *nCosperFgperDs=(int *)malloc(sizeof(**nCosperFgperDs)*nCosperFg.size());
  *ptrsT=(int *)malloc(sizeof(**ptrsT)* *nConfs* *trainingSize);
  *ptrsV=(int *)malloc(sizeof(**ptrsV)* nFg);
  *ptrsD=(int *)malloc(sizeof(**ptrsD)* *nDataset+1);
  *nDperD=(int *)malloc(sizeof(**nDperD)* *nDataset+1);
  *DihFgindx=(int *)malloc(sizeof(**DihFgindx)*nDih);
  *nVperFgi=(int *)malloc(sizeof(**nVperFgi)*nFg);

  /* Get the tset (nConfs * trainingSize) and tgts (nConfs)  */  
  for(int i=0;i<*nConfs;i++){
    (*tgts)[i]=data[i][*trainingSize];
    std::cout << "Target energy for conf " << i << " = " << data[i][*trainingSize] << std::endl;
    std::cout << "Training values for conf " << i << std::endl;
    for(int j=0;j<*trainingSize;j++){
      (*tset)[i* *trainingSize+j]=data[i][j];
        std::cout.precision(4);
        std::cout << std::setw(4) << (*tset)[i* *trainingSize+j] << " ";
    }
    std::cout << "\n" << std::endl;
  }
  
  /* Populate the fitting group index for printing out frdcmod etd*/
  std::cout << "Fitting group index for dihedrals: " << std::endl;
  for(int i=0;i<nDih;i++){
    (*DihFgindx)[i]=holdFg4dih[i];
    std::cout << (*DihFgindx)[i] << " ";
  }
  std::cout << "\n" << std::endl;
 
  /* populate the breaks and weights pointers which points to the first element of each dataset */
  std::cout << "brks(index for 1st conf in a Ds): " << std::endl;
  for(int i=0;i<weights.size();i++){
    (*brks)[i]=breaks[i];
    (*wts)[i]=weights[i];
     std::cout << (*brks)[i] << " "; 
  }
  std::cout << "\n" << std::endl;

  std::cout << "weights (index for Ds weights ): " << std::endl;
  for(int i=0;i<weights.size();i++){
    (*wts)[i]=weights[i];
     std::cout << (*wts)[i] << " ";
  }
  std::cout << "\n" << std::endl;

  std::cout << "Total conf (nConf): " << std::endl;
  for(int i=weights.size();i<*nBrks;i++){
    (*brks)[i]=breaks[i];
     std::cout << (*brks)[i] << " "; 
  }
  std::cout << "\n" << std::endl;

  /* Number of dihedrals in a given dataset, size is nDataset */
  int sumD=0;
  std::cout << "nDihperDs: " << std::endl;
  for (int i=0;i<*nDataset;i++) {
    (*nDperD)[i]=nDihperDs.at(i);
    sumD+=nDihperDs.at(i);
    std::cout << (*nDperD)[i] << " ";
  }
  std::cout << "\n" ;
  *totdih=sumD; 
  (*nDperD)[*nDataset]=sumD; //The last index is the sum maybe use this instead of totdih TODO
   
  /* All fitting groups in  given dataset, size is nDataset */
  *allFgsinDsi=(int *)malloc(sizeof(**allFgsinDsi)* *totdih); // allocate here cus I need the sum of all columns in Ds 
  std::cout << "allFginDs: " << std::endl;
  int indx=0;
  for (int i=0;i<*nDataset;i++) {
    for (int j=0;j<nDihperDs.at(i);j++) {
      std::cout << allFgsinDs.at(indx) << " ";
      (*allFgsinDsi)[indx]= allFgsinDs.at(indx);
      indx++;
    }
    std::cout << "\n";
  }

  /* total type of fitting groups in  given dataset, size is nDataset */
  int nFperD[*nDataset];
  std::cout << "nFgperDs: " << std::endl;
  for (int i=0;i<*nDataset;i++) {
    std::cout << nFgperDs.at(i) << " ";
    nFperD[i]=nFgperDs.at(i);
  }
  std::cout << "\n";
 
  /* Fitting groups in a DS  */
  indx=1; //index from 1 because index 0 is -99 
  std::cout << "FgsinDs: " << std::endl;
  for (int i=0;i<*nDataset;i++) {
    for (int f=0;f<nFperD[i];f++){
      std::cout << FgsinDs.at(indx) << " ";
      indx++;
    }
    std::cout << "\n";
  }
  indx=0;

  /* Number of Cosines per fitting group per Dataset */
  std::cout << "nCosperFgsperDs: " << std::endl;
  for (int i=0;i<*nDataset;i++) {
    for(int j=0;j<nFgperDs[i];j++){
      (*nCosperFgperDs)[indx]= nCosperFg[indx];
      std::cout << (*nCosperFgperDs)[indx] << " ";
      indx++;
    }
    std::cout << "\n";
  }
  indx=0;

  /* Get the Number of Vs needed per fitting group (nVperFg) it is nCosperFg * ndihperFg */ 
  int sum=0;
  std::cout << "nVperFg: " << std::endl;
  for(int i=0;i<nFg;i++){
    (*nVperFgi)[i]= nVperFg.at(i);
    sum+=nVperFg.at(i);
    std::cout << (*nVperFgi)[i] << " ";
  }
  std::cout << "\n";
  /* Get the genomesize which is the number of Vs */
  *genomeSize=sum;
  printf("genomeSize is %d \n", *genomeSize);

 
  /* tset ptrs index   */
  std::cout << "ptrsT: " << std::endl;
  for(int i=0;i<*nConfs;i++){
    for(int j=0;j<*trainingSize;j++){
      (*ptrsT)[i* *trainingSize+j]=tsetptDs[i][j];
    //for (int j=0;j<nDihperDs.at(i);j++) {
      std::cout << tsetptDs[i][j] << " ";
      indx++;
    }
    std::cout << "\n";
  }

  /* Vs ptrs index based on nVperFg, only need 1 copy since Vs indexing are parallelize*/
  std::cout << "ptrsV: " << std::endl;
  (*ptrsV)[0]=0;
  std::cout << (*ptrsV)[0] << " ";
  for (int i=1;i<nFg;i++) {
    (*ptrsV)[i]=(*ptrsV)[i-1] + nVperFg[i-1];
    std::cout << (*ptrsV)[i] << " ";
  }
  //std::cout << (*ptrsV)[nFg] << " "; //the last one 
  std::cout << "\n";
  
  /* dihedral index pointers */
  std::cout << "ptrsD: " << std::endl;
  (*ptrsD)[0]=0;
  std::cout << (*ptrsD)[0] << " ";
  for (int i=1;i<(*nDataset+1);i++) {
    //std::cout << "i is " << i << " (*ptrsD)[i-1] " << (*ptrsD)[i-1] <<" (*nDperD)[i-1] = " << (*nDperD)[i-1] <<std::endl;
    (*ptrsD)[i]=(*ptrsD)[i-1] + (*nDperD)[i-1];
    std::cout << (*ptrsD)[i] << " ";
  }
  std::cout << "\n";
  
}// end of program 
