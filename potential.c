/*************************************************************************

   Program:    
   File:       
   
   Version:    
   Date:       
   Function:   
   
   Copyright:  (c) University of Reading / Dr. Andrew C. R. Martin 2001
   Author:     Dr. Andrew C. R. Martin
   Address:    School of Animal and Microbial Sciences,
               The University of Reading,
               Whiteknights,
               P.O. Box 228,
               Reading RG6 6AJ.
               England.
   Phone:      +44 (0)118 987 5123 Extn. 7022
   Fax:        +44 (0)118 931 0180
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  02.08.02  Original
   V1.1  06.08.02  Corrected to ignore column and row zero in the matrices
   V1.2  23.10.02  Fixed -d option and reports error if unable to read
                   RAM file

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/SysDefs.h"

/************************************************************************/
/* Defines and macros
*/
#define NMATRIX  20
#define MAXLINES 10000
#define MAXBUFF  512
#define RAMBUFF  2000
#define MATSIZE  168
#define RAMDATA  "RAM.gen.par"
#define DATAENV  "DATADIR"

/************************************************************************/
/* Globals
*/
REAL gPot[NMATRIX][MATSIZE][MATSIZE];

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
REAL CalcEnergy(FILE *fp, PDB *pdb, BOOL verbose);
REAL FindEnergy(PDB *p, PDB *q);
BOOL ReadDatafile(FILE *pf);
int Convert(PDB *p);
int get_distance_bin_one(REAL distance);
void PatchOffsets(PDB *pdb);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *verbose, char *ramfile); 
void Usage(void);


/************************************************************************/
int main(int argc, char **argv)
{
   int counta;
   REAL TotalEnergy;
   PDB *pdb;
   int natoms;
   char  infile[MAXBUFF],
      outfile[MAXBUFF],
      ramfile[MAXBUFF];
   FILE  *in      = stdin,
      *out     = stdout,
      *datafp = NULL;
   BOOL verbose, noenv;
   
   
   strcpy(ramfile, RAMDATA);
   counta = 0;
   TotalEnergy = 0.0;

   if(ParseCmdLine(argc, argv, infile, outfile, &verbose, ramfile))
   {
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         if((datafp=OpenFile(ramfile, DATAENV, "r", &noenv))!=NULL)
         {
            /* Opens the potentials datafile and reads it into an array */
            if(!ReadDatafile(datafp))
            {
               fprintf(stderr,"Can't read datafile\n");
               return(1);
            }
            fclose(datafp);
   
            pdb=ReadPDB(in, &natoms);
            PatchOffsets(pdb);
   
            TotalEnergy = CalcEnergy(out, pdb, verbose);
   
            fprintf(out, "Total Energy = %f\n", TotalEnergy);
         }
         
      }
      else
      {
         fprintf(stderr, "Unable to read RAM potential file: %s\n",
                 ramfile);
         if(noenv)
         {
            fprintf(stderr, "Environment variable not set: %s\n",
                    DATAENV);
         }
      }
   }
   else
   {
      Usage();
   }

   return(0);
}


/************************************************************************/
void PatchOffsets(PDB *pdb)
{
   PDB *p;
   for(p=pdb; p!=NULL; NEXT(p))
   {
      p->occ = (REAL)Convert(p);
   }
}


/************************************************************************/
REAL CalcEnergy(FILE *out, PDB *pdb, BOOL verbose)
{
   PDB *p, *q, *start, *stop, *start2;
   REAL totalEnergy = (REAL)0.0, energy;
   
   /* Run through the protein, one residue at a time                    */
   for(start=pdb; start!=NULL; start=stop)
   {
      /* Find the start of the next residue                             */
      stop = FindNextResidue(start);
      /* And the start of the residue after that                        */
      start2 = FindNextResidue(stop);
      
      /* If there was a residue at position n+2                         */
      if(start2!=NULL)
      {
         /* Step through all the atoms in this residue                  */
         for(p=start; p!=stop; NEXT(p))
         {
            /* and all the atoms from start of residue n+2 onwards      */
            for(q=start2; q!=NULL; NEXT(q))
            {
               energy = FindEnergy(p,q);
               if(verbose)
               {
                  fprintf(out, "%c%d%c (%s %s) - %c%d%c (%s %s) : %.3f\n",
                          p->chain[0], p->resnum, p->insert[0],
                          p->resnam, p->atnam,
                          q->chain[0], q->resnum, q->insert[0],
                          q->resnam, q->atnam,
                          energy);
               }
               totalEnergy += energy;
            }
         }
      }
   }

   return(totalEnergy);
}


/************************************************************************/
REAL FindEnergy(PDB *p, PDB *q)
{
   REAL distance, 
        potential = 0.0;
   int  matrixnum, pos1, pos2;
   
   
   distance = DIST(p, q);
   
   if((matrixnum = get_distance_bin_one(distance))>=0)
   {
      pos1 = (int)p->occ;
      pos2 = (int)q->occ;
      
      if((pos1 < 0) || (pos2 < 0))
      {
         if((pos1 < 0) && (pos1 > (-10)))
         {
            fprintf(stderr,"Atom %s %d (%s %c %d) ignored\n",
                    p->atnam, p->atnum, p->resnam, p->chain[0], p->resnum);
            p->occ = (-10);
         }
         if((pos2 < 0) && (pos2 > (-10)))
         {
            fprintf(stderr,"Atom %s %d (%s %c %d) ignored\n",
                    q->atnam, q->atnum, q->resnam, q->chain[0], q->resnum);
            q->occ = (-10);
         }
         potential = (REAL)0.0;
      }
      else
      {
         potential = gPot[matrixnum][pos1][pos2];
      }
      
   }
   
   return(potential);
}


/************************************************************************/
BOOL ReadDatafile(FILE *fp)
{
   char *chp;
   int linenum=0, matrixnum, matrixline, fieldnum;
   char buffer[RAMBUFF], field[MAXBUFF];
   
   
   while(fgets(buffer, RAMBUFF, fp))
   {
      TERMINATE(buffer);
      chp = buffer;
      while(*chp==' ')chp++;
      matrixnum = (int)(linenum / MATSIZE);
      matrixline = linenum%MATSIZE;
      fieldnum = 0;
      /* Throw away the first field */
      chp = GetWordNC(chp, field, MAXBUFF);
      
      while(chp && (fieldnum<MATSIZE))
      {
         chp = GetWordNC(chp, field, MAXBUFF);
         sscanf(field, "%lf", &(gPot[matrixnum][matrixline][fieldnum++]));
      }

      linenum++;
   }

   return(TRUE);
}
   

/************************************************************************/
int Convert(PDB *p)
{
   int pos = (-1);
   char atnam[8];
   
   strcpy(atnam, p->atnam);

   if(!strncmp(atnam, "OXT ", 4))
   {
      strcpy(atnam, "O   ");
   }
   
   if (!strncmp(p->resnam, "ALA", 3))
   {
      if (!strncmp(atnam, "N  ", 3))
      {
         pos = 1;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 2;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 3;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 4;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 5;
      }
   }
   else if (!strncmp(p->resnam, "CYS", 3))
   {
      if (!strncmp(atnam, "N  ", 3))
      {	
         pos = 6;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 7;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 8;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 9;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 10;
      }
      else if (!strncmp(atnam, "SG ", 3))
      {
         pos = 11;
      }
   }
   else if (!strncmp(p->resnam, "ASP", 3))
   {
      if (!strncmp(atnam, "N  ", 3))
      {
         pos = 12;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 13;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 14;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 15;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 16;
      }
      else if (!strncmp(atnam, "CG ", 3))
      {
         pos = 17;
      }
      else if (!strncmp(atnam, "OD1", 3))
      {
         pos = 18;
      }
      else if (!strncmp(atnam, "OD2", 3))
      {
         pos = 19;
      }
   }
   else if (!strncmp(p->resnam, "GLU", 3))
   {
      if (!strncmp(atnam, "N  ", 3))
      {
         pos = 20;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 21;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 22;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 23;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 24;
      }
      else if (!strncmp(atnam, "CG ", 3))
      {
         pos = 25;
      }
      else if (!strncmp(atnam, "CD ", 3))
      {
         pos = 26;
      }
      else if (!strncmp(atnam, "OE1", 3))
      {
         pos = 27;
      }
      else if (!strncmp(atnam, "OE2", 3))
      {
         pos = 28;
      }
   }
   else if (!strncmp(p->resnam, "PHE", 3))
   {
      if (!strncmp(atnam, "N  ", 3))
      {
         pos = 29;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 30;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 31;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 32;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 33;
      }
      else if (!strncmp(atnam, "CG ", 3))
      {
         pos = 34;
      }
      else if (!strncmp(atnam, "CD1", 3))
      {
         pos = 35;
      }
      else if (!strncmp(atnam, "CD2", 3))
      {
         pos = 36;
      }
      else if (!strncmp(atnam, "CE1", 3))
      {
         pos = 37;
      }
      else if (!strncmp(atnam, "CE2", 3))
      {
         pos = 38;
      }
      else if (!strncmp(atnam, "CZ ", 3))
      {
         pos = 39;
      }
   }
   else if (!strncmp(p->resnam, "GLY", 3))
   {
      if (!strncmp(atnam, "N  ", 3))
      {
         pos = 40;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 41;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 42;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 43;
      }
   }
   else if (!strncmp(p->resnam, "HIS", 3))
   {
      if (!strncmp(atnam, "N  ", 3))
      {
         pos = 44;
      } 
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 45;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 46;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 47;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 48;
      }
      else if (!strncmp(atnam, "CG ", 3))
      {
         pos = 49;
      }
      else if (!strncmp(atnam, "ND1", 3))
      {
         pos = 50;
      }
      else if (!strncmp(atnam, "CD2", 3))
      {
         pos = 51;
      }
      else if (!strncmp(atnam, "CE1", 3))
      {
         pos = 52;
      }
      else if (!strncmp(atnam, "NE2", 3))
      {
         pos = 53;
      }
   }
   else if (!strncmp(p->resnam, "ILE", 3))
   {
      if (!strncmp(atnam, "N  ", 3))
      {
         pos = 54;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 55;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 56;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 57;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 58;
      }
      else if (!strncmp(atnam, "CG1", 3))
      {
         pos = 59;
      }
      else if (!strncmp(atnam, "CG2", 3))
      {
         pos = 60;
      }
      else if (!strncmp(atnam, "CD1", 3))
      {
         pos = 61;
      }
   }
   else if (!strncmp(p->resnam, "LYS", 3))
   {
      if (!strncmp(atnam, "N  ", 3))
      {
         pos = 62;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 63;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 64;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 65;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 66;
      }
      else if (!strncmp(atnam, "CG ", 3))
      {
         pos = 67;
      }
      else if (!strncmp(atnam, "CD ", 3))
      {
         pos = 68;
      }
      else if (!strncmp(atnam, "CE ", 3))
      {
         pos = 69;
      }	
      else if (!strncmp(atnam, "NZ ", 3))
      {
         pos = 70;
      }
   }
   else if (!strncmp(p->resnam, "LEU", 3))
   {
      if (!strncmp(atnam, "N  ", 3))
      {
         pos = 71;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 72;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 73;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 74;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 75;
      }
      else if (!strncmp(atnam, "CG ", 3))
      {
         pos = 76;
      }
      else if (!strncmp(atnam, "CD1", 3))
      {
         pos = 77;
      }
      else if (!strncmp(atnam, "CD2", 3))
      {
         pos = 78;
      }
   }
   else if (!strncmp(p->resnam, "MET", 3))
   {
      if (!strncmp(atnam, "N  ", 3))
      {
         pos = 79;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 80;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 81;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 82;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 83;
      }
      else if (!strncmp(atnam, "CG ", 3))
      {
         pos = 84;
      }
      else if (!strncmp(atnam, "SD ", 3))
      {
         pos = 85;
      }
      else if (!strncmp(atnam, "CE ", 3))
      {
         pos = 86;
      }
   }
   else if (!strncmp(p->resnam, "ASN", 3))
   {
      if (!strncmp(atnam, "N  ", 3))
      {
         pos = 87;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 88;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 89;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 90;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 91;
      }
      else if (!strncmp(atnam, "CG ", 3))
      {
         pos = 92;
      }
      else if (!strncmp(atnam, "OD1", 3))
      {
         pos = 93;
      }
      else if (!strncmp(atnam, "ND2", 3))
      {
         pos = 94;
      }
   }
   else if (!strncmp(p->resnam, "PRO", 3))
   {
      if (!strncmp(atnam, "N  ", 3))
      {
         pos = 95;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 96;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 97;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 98;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 99;
      }
      else if (!strncmp(atnam, "CG ", 3))
      {
         pos = 100;
      }
      else if (!strncmp(atnam, "CD ", 3))
      {
         pos = 101;
      }
   }
   else if (!strncmp(p->resnam, "GLN", 3))
   {
      if (!strncmp(atnam, "N  ", 3))
      {
         pos = 102;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 103;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 104;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 105;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 106;
      }
      else if (!strncmp(atnam, "CG ", 3))
      {
         pos = 107;
      }
      else if (!strncmp(atnam, "CD ", 3))
      {
         pos = 108;
      }
      else if (!strncmp(atnam, "OE1", 3))
      {
         pos = 109;
      }
      else if (!strncmp(atnam, "NE2", 3))
      {
         pos = 110;
      }
   }
   else if (!strncmp(p->resnam, "ARG", 3))
   {
      if(!strncmp(atnam, "N  ", 3))
      {
         pos = 111;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 112;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 113;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 114;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 115;
      }
      else if (!strncmp(atnam, "CG ", 3))
      {
         pos = 116;
      }
      else if (!strncmp(atnam, "CD ", 3))
      {
         pos = 117;
      }
      else if (!strncmp(atnam, "NE ", 3))
      {
         pos = 118;
      }
      else if (!strncmp(atnam, "CZ ", 3))
      {
         pos = 119;
      }
      else if (!strncmp(atnam, "NH1", 3))
      {
         pos = 120;
      }
      else if (!strncmp(atnam, "NH2", 3))
      {
         pos = 121;
      }
   }
   else if (!strncmp(p->resnam, "SER", 3))
   {
      if(!strncmp(atnam, "N  ", 3))
      {
         pos = 122;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 123;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 124;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 125;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 126;
      }
      else if (!strncmp(atnam, "OG ", 3))
      {
         pos = 127;
      }
   }
   else if (!strncmp(p->resnam, "THR", 3))
   {
      if(!strncmp(atnam, "N  ", 3))
      {
         pos = 128;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 129;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 130;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 131;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 132;
      }
      else if (!strncmp(atnam, "OG1", 3))
      {
         pos = 133;
      }
      else if (!strncmp(atnam, "CG2", 3))
      {
         pos = 134;
      }
   }
   else if (!strncmp(p->resnam, "VAL", 3))
   {
      if(!strncmp(atnam, "N  ", 3))
      {
         pos = 135;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 136;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 137;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 138;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 139;
      }
      else if (!strncmp(atnam, "CG1", 3))
      {
         pos = 140;
      }
      else if (!strncmp(atnam, "CG2", 3))
      {
         pos = 141;
      }
   }
   else if (!strncmp(p->resnam, "TRP", 3))
   {
      if(!strncmp(atnam, "N  ", 3))
      {
         pos = 142;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 143;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 144;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 145;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 146;
      }
      else if (!strncmp(atnam, "CG ", 3))
      {
         pos = 147;
      }
      else if (!strncmp(atnam, "CD1", 3))
      {
         pos = 148;
      }
      else if (!strncmp(atnam, "CD2", 3))
      {
         pos = 149;
      }
      else if (!strncmp(atnam, "NE1", 3))
      {
         pos = 150;
      }
      else if (!strncmp(atnam, "CE2", 3))
      {
         pos = 151;
      }
      else if (!strncmp(atnam, "CE3", 3))
      {
         pos = 152;
      }
      else if (!strncmp(atnam, "CZ2", 3))
      {
         pos = 153;
      }
      else if (!strncmp(atnam, "CZ3", 3))
      {
         pos = 154;
      }
      else if (!strncmp(atnam, "CH2", 3))
      {
         pos = 155;
      }
   }
   else if (!strncmp(p->resnam, "TYR", 3))
   {
      if(!strncmp(atnam, "N  ", 3))
      {
         pos = 156;
      }
      else if (!strncmp(atnam, "CA ", 3))
      {
         pos = 157;
      }
      else if (!strncmp(atnam, "C  ", 3))
      {
         pos = 158;
      }
      else if (!strncmp(atnam, "O  ", 3))
      {
         pos = 159;
      }
      else if (!strncmp(atnam, "CB ", 3))
      {
         pos = 160;
      }
      else if (!strncmp(atnam, "CG ", 3))
      {
         pos = 161;
      }
      else if (!strncmp(atnam, "CD1", 3))
      {
         pos = 162;
      }
      else if (!strncmp(atnam, "CD2", 3))
      {
         pos = 163;
      }
      else if (!strncmp(atnam, "CE1", 3))
      {
         pos = 164;
      }
      else if (!strncmp(atnam, "CE2", 3))
      {
         pos = 165;
      }
      else if (!strncmp(atnam, "CZ ", 3))
      {
         pos = 166;
      }
      else if (!strncmp(atnam, "OH ", 3))
      {
         pos = 167;
      }
   }

   return(pos);
}


/************************************************************************/
int get_distance_bin_one(REAL distance)
{
  if (distance < 3.0) return 0;
  if (distance < 4.0) return 1;
  if (distance < 5.0) return 2;
  if (distance < 6.0) return 3;
  if (distance < 7.0) return 4;
  if (distance < 8.0) return 5;
  if (distance < 9.0) return 6;
  if (distance < 10.0) return 7;
  if (distance < 11.0) return 8;
  if (distance < 12.0) return 9;
  if (distance < 13.0) return 10;
  if (distance < 14.0) return 11;
  if (distance < 15.0) return 12;
  if (distance < 16.0) return 13;
  if (distance < 17.0) return 14;
  if (distance < 18.0) return 15;
  if (distance < 19.0) return 16;
  if (distance < 20.0) return 17;
  
  return(-1);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     BOOL *verbose, char *ramfile)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            char   *ramfile     RAM data file
   Returns: BOOL                Success?

   Parse the command line
   
   09.08.95 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *verbose, char *ramfile)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   *verbose = FALSE;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'v':
            *verbose = TRUE;
            break;
         case 'd':
            argc--;
            argv++;
            strncpy(ramfile, argv[0], MAXBUFF);
            break;
         case 'h':
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 0--2 arguments left                    */
         if(argc > 2)
            return(FALSE);
         
         if(argc)
         {
            /* Copy the next to infile                                  */
            strcpy(infile, argv[0]);
            argc--;
            argv++;
         
            /* If there's another, copy it to outfile                   */
            if(argc)
               strcpy(outfile, argv[0]);
         }
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
void Usage(void)
{
   fprintf(stderr,"\npotential V1.2 (c) 2002, Dr. Andrew C.R. Martin, \
University of Reading\n");
   fprintf(stderr,"\nUsage: potential [-v] [-d datafile] [in.pdb \
[out.txt]]\n");
   fprintf(stderr,"       -v   Verbose - print all atom pairs and their \
energies\n");
   fprintf(stderr,"       -d   Specify the datafile (default: %s)\n",
           RAMDATA);

   fprintf(stderr,"\nCalculates a pseudo-energy term for a protein \
structure using Ram\n");
   fprintf(stderr,"Samudrala's empirical potential (see \
http://prostar.carb.nist.gov/).\n");
   fprintf(stderr,"The location of the datafile should either be \
specified explicitly,\n");
   fprintf(stderr,"be in the current directory, or in the directory \
defined by the\n");
   fprintf(stderr,"%s environment variable.\n\n", DATAENV);
}
