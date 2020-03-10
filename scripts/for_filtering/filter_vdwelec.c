#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/******************************************************************************************************
* This program checks if a pair of atoms is too close to each other and if a side chain atom has a strong
   electrostatic interaction with backbone atoms during the dihedral fitting stage
    written by James Maier & Kellon Belfon 
*******************************************************************************************************/

/*******************************************************************************
 * Usage
./filter_vdwelec top [crd1 crd2 ...] 
it will print a VDW.log, ELEC.log, bond_info.log & Suspect_structures.dat

*******************************************************************************/

/**********************************
 * CONSTANTS
 **********************************/
// array of elements by atomic number
enum elements_t {H=1,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Ge,As,Se,Br,Kr,Rb};

// If you add anything here, change the size of the arrray and change BB_ATOM and SS_ATOM size as well
// array of backbone atom name, 9 atom name with 4 space for each
char backbone[12][4] = {"HH31","HH33","HH32","CH3","C","CA","O","N","H", "H2", "H1", "H3"};
// array of sidechain heavy atom name and polar hydrogens, atom name with 4 space for each
char sidechain[77][4] = {"CB","CD","CD1","CD2","CH2","CG","CG1","CG2","CE","CE1","CE2","CE3",
                         "CZ","CZ2","CZ3","OD1","OD2","OE1","OE2","OG","OG1","OH","ND","ND1",
                         "ND2","NE","NE1","NE2","NH2","NZ","SD","SED","SG","O1P","O2P","O3P", 
                         "P", "H1P", "HO","HH21","HH22","HH11","HH12","HD2","HD21","HD22",
                         "HG","HE2", "HIA", "CIA", "CI", "OT", "NT", "CHB", "CHA", "CH",
                         "C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17",
                         "C18","C19","C110","C111","C112","C113","C114","C115","C116","O3"};

// atom type structure. xyz is the coordinates
struct atom_t {
    enum elements_t element;
    int numBonds;
    int curBond;
    int bonds[4];
    char atomName[5];
    double charge;
    double vdw;
    double x;
    double y;
    double z;
};

/*************************************************
 * FUNCTIONS
atom_t contain a structure of elements, numBonds
curBonf, bonds, x , y ,z
 *************************************************/

/* function calculate distance between a pair of atoms */

double distance(struct atom_t *a, struct atom_t *b) {
/*assigning x, y and z coordinates to an array to calculate distance between atoms a and atom b
  distance formula: sqrt( (x2-x1)^2 + ((y2-y1)^2 + (z2-z1)^2 )*/
    double t,r,u;
    r = a->x - b->x;   //r = x2 - x1
    r *= r;          //r = r*r = (x2 - x1)^2
    t = a->y - b->y;   //t = y2 - y1
    r += t * t;        //r = r + t*t = (x2 - x1)^2 + (y2 - y1)^2   
    u =a->z-b->z;   //u = z2 - z1
    r += u*u;        //r = r + u*u = (x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2
    return sqrt(r); //sqrt(r) = sqrt( (x2-x1)^2 + ((y2-y1)^2 + (z2-z1)^2)  
}

/* function calculate electrostic energy using E = q1*q2/distance  
TODO: Add 1-4 electrostatics scaling */
double Elec_E(double q1, double q2, double dist){
    double E; 
    E = (q1 * q2) / dist;
    return E;
}

/* function to extract information */
void maskWithin(struct atom_t *atoms, int from, int numBonds, int *mask) {
    int a;
    //printf("  masking %d\n",from+1);
    mask[from]=1;
    if(numBonds>0) {
        numBonds--;
        for(a=0;a<atoms[from].numBonds;a++) {
            if(!mask[atoms[from].bonds[a]]) {
                maskWithin(atoms,atoms[from].bonds[a],numBonds,mask);
            }
        }
    }
}

/**********************************
 * MAIN PROGRAM
 **********************************/
 
#define BUFF_SIZE 256
#define BB_ATOM 12
#define SC_ATOM 77

int main(int argc, char *argv[]) {
    /* declarations */
    char buffer[BUFF_SIZE];
    struct atom_t *atoms; 
    int numAtoms=0;
    int numBondsWithoutH, numBondsIncH, comp, comp1;
    FILE *in;
    int *mask;
    int a,b,c,f,sc,bb;
    double dist;
    double elec;  
    FILE *vdwlog;
    FILE *eleclog;
    FILE *logfile;
    FILE *badvdw;
    FILE *badelec;
    float vdw_filter=1.3; // if distance less than sum of vdw pairs/vdw_filter, james used 1.3
    float elec_filter=42.0;  //if electrostic energy is greater than elec_filter do not keep structure, james used 42


    /* open the output file which is the log file */
    vdwlog = fopen("VDW.log", "w");
    eleclog = fopen("ELEC.log", "w");
    logfile = fopen("Bond_info.log", "w");
    badelec = fopen("Suspect_structures.elec.dat", "w");
    badvdw = fopen("Suspect_structures.vdw.dat", "w");


    /* load an amber topology file */
    in=fopen(argv[1],"r");
    //fputs("LOading top.\n",stderr);

    /* look for the %FLAG POINTERS in AMBER parm7 file, which have 14 spaces
       read this for amber file format http://ambermd.org/formats.html */
    while(fgets(buffer,BUFF_SIZE,in)&&strncmp(buffer,"%FLAG POINTERS",14)) ;
    fgets(buffer,BUFF_SIZE,in); //skip line
    fgets(buffer,BUFF_SIZE,in); 
    /*allocates space for destination of parsed lines from fgets function 
        buffer has 8 spaces this first number on the line is the number of atoms */
    buffer[8]=0; 
    numAtoms=atoi(buffer); // total number of atoms
    fprintf(logfile, "Total number of atoms in topology file: %d atoms\n\n", numAtoms);
    /* 32 spaces into the line is number of bonds without H */
    buffer[32]=0;
    numBondsWithoutH=atoi(buffer+24);
    fprintf(logfile, "Number of bonds without Hydrogen: %d \n\n", numBondsWithoutH);
    /* 24 spaces into the line is the number of bonds with H */
    buffer[24]=0;
    numBondsIncH=atoi(buffer+16);
    fprintf(logfile, "Number of bonds with Hydrogen: %d \n\n", numBondsIncH);
    // done with the first Flag, look at a parm7 file if you want to understand

    /*allocate mememory for atoms and setting up the array atoms */
    atoms=malloc(sizeof(struct atom_t)*numAtoms);

    // loop through atom structure and assign 0 to numBonds 
    for(c=0; c<numAtoms; c++) {
	atoms[c].numBonds=0;
    }
    mask=malloc(sizeof(int)*numAtoms);

    /* continue parsing the parm7 file, next flag we want is ATOM_NAME */
    while(fgets(buffer,BUFF_SIZE,in)&&strncmp(buffer,"%FLAG ATOM_NAME",15)) ;
    fgets(buffer,BUFF_SIZE,in); // next line
    for(c=0; c<numAtoms; c++) {
        fscanf(in,"%4s",atoms[c].atomName);
        //printf("%s\n", atoms[c].atomName);
    }
   
    /* continue parsing the parm file, next flag we want is the CHARGE */
    while(fgets(buffer,256,in)&&strncmp(buffer,"%FLAG CHARGE",12)) ;
    fgets(buffer,128,in);
    for(c=0; c<numAtoms; c++) {
        fscanf(in,"%lf", &atoms[c].charge);
        /* amber uses E = q1*q2/r, where E is kcal/mol, r in angstroms and q is what is in the parmfile, so leave as is
         but if u want actually charges like in the lib file uncomment below */
         //atom[c].charge/=18.222615;  //parm file evaluate charges by * scaling factor 18.222615
         //printf("%lf\n", atoms[c].charge);
    }
	
    /* continue parsing the parm file, next flag we want is the ATOMIC NUMBER */
    while(fgets(buffer,128,in)&&strncmp(buffer,"%FLAG ATOMIC_NUMBER",19)) ;
    fgets(buffer,128,in);
    /* for every atoms get its atomic number by scanning the line and matching
           e.g when read 1 match this to H, when read 6 match to C */
    for(c=0; c<numAtoms; c++) {
        fscanf(in,"%d", &atoms[c].element);
        //printf("%d\n", atoms[c].element);
    }
	
         
    /* continue parsing the parm file, next we will get the bonds information  */
    /* first, bonds with hydrogen */
    while(fgets(buffer,128,in)&&strncmp(buffer,"%FLAG BONDS_INC_HYDROGEN",24)) ;
    fgets(buffer,128,in);
    for(c=0; c<numBondsIncH; c++) {
        fscanf(in,"%d %d %*d",&a,&b);
	a/=3; b/=3; 
        //printf("a=%d, b=%d\n", a, b);
	//a+=1; b+=1; // atom index from parm file is number/3 + 1
        fprintf(logfile, "Bonds with Hydrogens: %d and %d\n", a+1, b+1);
        atoms[a].bonds[atoms[a].numBonds++]=b;
        atoms[b].bonds[atoms[b].numBonds++]=a;
    }
    /* second, bonds with hydrogen */
    while(fgets(buffer,128,in)&&strncmp(buffer,"%FLAG BONDS_WITHOUT_HYDROGEN",28)) ;
    fgets(buffer,128,in);
    for(c=0; c<numBondsWithoutH; c++) {
	fscanf(in,"%d %d %*d",&a,&b);
	a/=3; b/=3;
	//a+=1; b+=1; // atom index from parm file is number/3 + 1
        fprintf(logfile, "Bonds without hydrogens: %d and %d\n", a+1, b+1);
	atoms[a].bonds[atoms[a].numBonds++]=b;
	atoms[b].bonds[atoms[b].numBonds++]=a;
    }

    /* get the VDW radius from the topology file,  */
    while(fgets(buffer,128,in)&&strncmp(buffer,"%FLAG RADII",11));                                                                     
    fgets(buffer,128,in);
    for(c=0; c<numAtoms; c++) {
        fscanf(in,"%lf", &atoms[c].vdw);
        //fscanf(in,"%lf", vdw);
        //printf("%lf\n", atoms[c].vdw);
    } 

    /* Now open coordinate files and loop through to propose structures to reject based on two criteria
       (1) Reject structures if the distance between two atoms (not bonded or within an angle) is less than the
            sum of the Van der Waals divided by 1.3. 1.3 can be adjusted based on your liking
       (2) Reject structures if the coulombic energy between a side chain particle (Doing heavy atoms only) and
           the back bone exceeded 42 kcal/mol. 42 is picked based on Adams, 1979 (the answer to the universe).
            Any number between 40 and 43 will work 
    */
     
    /* loop through each file with coordinates */ 
    for(f=2;f<argc;f++) {
        if((in=(fopen(argv[f],"r")))) {
            fputs("Loading ",stderr);fputs(argv[f],stderr);fputc('\n',stderr);
	    fgets(buffer,BUFF_SIZE,in); //skip first line
	    fgets(buffer,BUFF_SIZE,in); //skip second line
            fprintf(vdwlog, "\n");
            fprintf(eleclog, "\n");
            fprintf(vdwlog, "# This is structure %s: \n\n", argv[f]);
            fprintf(vdwlog, "# If the distance between atom pairs is less than the sum of their VDW radius divided by 1.3, Reject_flag will be 1 \n");
            fprintf(vdwlog, "# If there is a 1 in  Reject_flag column, reject structure %s \n\n", argv[f]);
            fprintf(vdwlog, "#   Atom 1         Atom 2    VDW 1      VDW 2     Distance  Sum of VDW  (Sum of VDW)/%2.1f  Reject_flag\n", vdw_filter);
            fprintf(eleclog, "# This is structure %s: \n\n", argv[f]);
            fprintf(eleclog, "# If a give sidechain atom-backbone charge interaction is greater than 42 kcal/mol, Reject_flag will be 1\n"); 
            fprintf(eleclog, "# If there is a 1 in  Reject_flag column, reject structure %s \n\n", argv[f]);
            fprintf(eleclog, "#   Atom 1        Atom 2      charge 1(iu)  charge 2(iu)    charge 1      charge 2      distance     Elec energy    Reject_flag\n");
            /* dealing with vdw and elec cutoff, loading coordinates and calculating distance */
            /* Loop through the atoms */
	    for(c=0;c<numAtoms;c++) {
                atoms[c].x = 0.0; atoms[c].y = 0.0; atoms[c].z = 0.0;
                /* read in atom xyz coordinates and save to atoms array*/
	        fscanf(in,"%lf %lf %lf",&atoms[c].x, &atoms[c].y, &atoms[c].z);
		//printf("%lf %lf %lf\n", atoms[c].x, atoms[c].y, atoms[c].z);
            }

	    for(c=0;c<numAtoms;c++) {
                /* electrostatics part */
                /* for each atoms now loop through sidechains atom name to find the atom c that is a side chain atom */
                for(sc=0;sc<SC_ATOM;sc++) {
                    /* If sidechain atom name equals the atom c atom name then */
                    comp = strncmp(sidechain[sc], atoms[c].atomName, 4); // return 0, if strings are the same
                    if(comp == 0) { // means the two strings are identical
                        // printf("sidechain atomname: %s == current atomname: %s\n", sidechain[sc], atoms[c].atomName);
                        /* If yes, then loop through backbone atom's name and find an atom name that is a backbone and calculate the electrostatic interaction */
	                for(a=0;a<numAtoms;a++) {
                            for(bb=0;bb<BB_ATOM;bb++) {
                                comp1 = strncmp(backbone[bb], atoms[a].atomName, 4); // return 0, if strings are the same
                                if(comp1 == 0) { // means the two strings are identical
                                    // printf("backbone atomname: %s == current atomname: %s\n", backbone[bb], atoms[a].atomName);
                                    //printf("Distance will be calculated between atomname: %s and atomname: %s\n", atoms[c].atomName, atoms[a].atomName);
                                    /* get the distance between sidechain atoms[c] and backbone atoms[bb]*/
                                    //printf("xyz of atomname %s is %f %f %f \n", atoms[c].atomName, atoms[c].x, atoms[c].y, atoms[c].z);
                                    //printf("xyz of atomname %s is %f %f %f \n", atoms[a].atomName, atoms[a].x, atoms[a].y, atoms[a].z);
                                    dist = distance(atoms+c, atoms+a);
                                    /* get the electrostatic energy */
                                    elec = Elec_E(atoms[c].charge, atoms[a].charge, dist);
                                    // printf("elec energy of %s and %s is %lf:", atoms[c].atomName, atoms[a].atomName, elec); 
                                    if(elec > elec_filter){ // if electrostic energy is greater than 42 do not keep structure
                                         fprintf(eleclog, "%4d ( %4s ) %4d ( %4s )%14lf%14lf%14lf%14lf%14lf%14lf%12d\n", a+1, atoms[a].atomName, c+1, atoms[c].atomName, atoms[a].charge, atoms[c].charge, atoms[a].charge/18.222615, atoms[c].charge/18.222615, dist, elec, 1);
                                         fprintf(badelec, "%s   ELEC\n", argv[f]);
                                    }
                                    else {
                                         fprintf(eleclog, "%4d ( %4s ) %4d ( %4s )%14lf%14lf%14lf%14lf%14lf%14lf%12d\n", a+1, atoms[a].atomName, c+1, atoms[c].atomName, atoms[a].charge, atoms[c].charge, atoms[a].charge/18.222615, atoms[c].charge/18.222615, dist, elec, 0);
    
                                    }
                                }
                            }
			} 
                    }
                }      
                // VDW part                        
		for(a=0; a<numAtoms; a++) {
                    mask[a]=0;
		}
	        //printf("From %d:\n",c+1);
		maskWithin(atoms,c,2,mask);
	        for(a=0; a<c; a++) {
		    if(!mask[a]) {
                        /* If distance between pairs < sum of vdw/vdw_filter(1.3) */
                        //if(distance(atoms+a,atoms+c)<(vdWr[atoms[a].element]+vdWr[atoms[c].element])/1.3) {
                        if(distance(atoms+a,atoms+c)<(atoms[a].vdw+atoms[c].vdw)/vdw_filter) {
                           
                            fprintf(vdwlog, "%4d ( %4s ) %4d ( %4s )%10lf%10lf%12lf%12lf%14lf%12d\n", a+1, atoms[a].atomName, c+1, atoms[c].atomName, atoms[a].vdw, atoms[c].vdw, distance(atoms+a,atoms+c), atoms[a].vdw+atoms[c].vdw, (atoms[a].vdw+atoms[c].vdw)/vdw_filter, 1);
                            fprintf(badvdw, "%s   VDW\n", argv[f]);
		        }
			else {
                            fprintf(vdwlog, "%4d ( %4s ) %4d ( %4s )%10lf%10lf%12lf%12lf%14lf%12d\n", a+1, atoms[a].atomName, c+1, atoms[c].atomName, atoms[a].vdw, atoms[c].vdw, distance(atoms+a,atoms+c), atoms[a].vdw+atoms[c].vdw, (atoms[a].vdw+atoms[c].vdw)/vdw_filter, 0);
			}
		    }  
	        }
            }
            fprintf(vdwlog, "#END of FILE for Structure %s\n",argv[f]);
            fprintf(eleclog, "#END of FILE for Structure %s\n",argv[f]);
            fprintf(vdwlog, "\n");
            fclose(in);
	}
	else fprintf(stderr,"Unable to open %s\n",argv[f]);
    }
    
    fclose(vdwlog);
    fclose(eleclog);    
    free(atoms);
    free(mask);
    return 0;
}


/*
if(hydrogen connected to electronegative) {
	find bonded atom, is it oxygen, nitrogen, etc.?
	for(every other electronegative atom) {
		if(distance < sum of van der Waals radii) {
			discount this structure
		}
	}
}*/
