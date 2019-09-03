#ifndef __GO_MODEL_BUILDER_IO_H__
#define __GO_MODEL_BUILDER_IO_H__
#include "DataStructure.h"
#include "DataType.h"
#include "stdio.h"
//#define testWriteKBGoTop
//#define testWriteKBGoPDB
//#define testWriteKBGoParams
//write topology file

int writeKBGoTop(const Molecule *InMol,const char *filename);
int writeKBGoPDB(const Molecule *InMol,const char *filename);
int writeKBGoParams(const Molecule *InMol,const char *filename);

int writeKBGoTop(const Molecule *InMol,const char *filename){
    FILE *ptr;
#ifdef testWriteKBGoTop
    ptr=stdout;
#else
    ptr=fopen(filename,"w");
#endif
    char currentChainId=InMol->ChainId[0];
    fprintf(ptr,"* This CHARMM .top file describes a Go model\n");
    fprintf(ptr,"*\n\n");
    fprintf(ptr,"read rtf card\n");
    fprintf(ptr,"* Topology for Go model of\n");
    fprintf(ptr,"*\n");
    fprintf(ptr,"   20   1\n");
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        char chainId;
        char resname[5];
        chainId=InMol->ChainId[ii];
        sprintf(resname,"%c%d",chainId,InMol->ResNum[ii]);
        fprintf(ptr,"MASS %-3d %-4s    %f\n",-1,resname,InMol->Mass[ii]);
    }
    fprintf(ptr,"\nDECL +CA\n\n");
    fprintf(ptr,"AUTOGENERATE ANGLES DIHEDRAL\n\n");
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        char chainId;
        char resname[5];
        chainId=InMol->ChainId[ii];
        sprintf(resname,"%c%d",chainId,InMol->ResNum[ii]);
        fprintf(ptr,"RESI %-4s          %-4f\n",resname,InMol->Charge[ii]);
        fprintf(ptr,"GROU\n");
        fprintf(ptr,"Atom  %-4s%-4s     %-4f\n",InMol->AtomType[ii],resname,InMol->Charge[ii]);
        fprintf(ptr,"Bond %-4s +%-4s\n\n",InMol->AtomType[ii],InMol->AtomType[ii]);
    }
    fprintf(ptr,"END\n\n");
#ifndef testWriteKBGoTop
    fclose(ptr);
#endif
    return 0;
}
//write PDB file
int writeKBGoPDB(const Molecule *InMol,const char *filename){
    FILE *ptr;
#ifdef testWriteKBGoPDB
    ptr=stdout;
#else
    ptr=fopen(filename,"w");
#endif
    char currentChainId=InMol->ChainId[0];
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        char chainId;
        char resname[5];
        chainId=InMol->ChainId[ii];
        if(currentChainId!=chainId){
            fprintf(ptr,"TER\n");
            currentChainId=chainId;
        }
        sprintf(resname,"%c%d",chainId,InMol->ResNum[ii]);
        if(mystrcmp(InMol->AtomType[ii],"ZN")==0){
            fprintf(ptr,"HETATM%5d  %-4s%-4s%c%4d",ii+1,InMol->AtomType[ii],resname,
                    chainId,ii+1);
        }else{
            fprintf(ptr,"ATOM  %5d  %-4s%-4s%c%4d",ii+1,InMol->AtomType[ii],resname,
                    chainId,ii+1);
        }
        fprintf(ptr,"%12.3f %7.3f %7.3f",InMol->X[ii],InMol->Y[ii],InMol->Z[ii]);
        fprintf(ptr,"  1.00  1.00      PRO%c\n", InMol->ChainId[ii]);
    }
    fprintf(ptr,"END\n");
#ifndef testWriteKBGoPDB
    fclose(ptr);
#endif
    return 0;
}
//write parameters file
int writeKBGoParams(const Molecule *InMol,const char *filename){
    FILE *ptr;
#ifdef testWriteKBGoParams
    ptr=stdout;
#else
    ptr=fopen(filename,"w");
#endif
    fprintf(ptr,"*This CHARMM .param file describes a Go model\n");
    fprintf(ptr,"*Scaling factor = %f\n",InMol->nativeScalingFactor);
    fprintf(ptr,"*\n");
    fprintf(ptr,"\n");
    fprintf(ptr,"read param card\n");
    fprintf(ptr,"*\n");
    //bond
    fprintf(ptr,"\nBOND\n");
    for(int ii=0; ii<InMol->NumOfBonds; ++ii){
        int atom1,atom2;
        atom1=InMol->BOND[ii].i;
        atom2=InMol->BOND[ii].j;
        float k,r;
        k=InMol->BOND[ii].k;
        r=InMol->BOND[ii].r;
        fprintf(ptr,"%c%-3d    %c%-3d      %f  %f\n",InMol->ChainId[atom1],atom1,
                InMol->ChainId[atom2],atom2,k,r);
    }
    fprintf(ptr,"\nANGLE\n");
    for(int ii=0; ii<InMol->NumOfAngles; ++ii){
        int atom1,atom2,atom3;
        atom1=InMol->ANGLE[ii].a;
        atom2=InMol->ANGLE[ii].b;
        atom3=InMol->ANGLE[ii].c;
        float k,theta;
        k=InMol->ANGLE[ii].k;
        theta=InMol->ANGLE[ii].theta;
        fprintf(ptr,"%c%-3d    %c%-3d    %c%-3d      %f  %f\n",InMol->ChainId[atom1],atom1,
                InMol->ChainId[atom2],atom2,InMol->ChainId[atom3],atom3,k,theta);
    }
    for(int ii=0; ii<InMol->NumOfZincAngles; ++ii){
        Angle angle=InMol->ZINCANGLE[ii];
        int atom1,atom2,atom3;
        atom1=angle.a;
        atom2=angle.b;
        atom3=angle.c;
        float k,theta;
        k=angle.k;
        theta=angle.theta;
        fprintf(ptr,"%c%-3d    %c%-3d    %c%-3d      %f  %f\n",InMol->ChainId[atom1],atom1,
                InMol->ChainId[atom2],atom2,InMol->ChainId[atom3],atom3,k,theta);
    }
    fprintf(ptr,"\nDIHEDRAL\n");
    for(int ii=0; ii<InMol->NumOfDihedrals; ++ii){
        int atom1,atom2,atom3,atom4,mult;
        atom1=InMol->DIHEDRAL[ii].a;
        atom2=InMol->DIHEDRAL[ii].b;
        atom3=InMol->DIHEDRAL[ii].c;
        atom4=InMol->DIHEDRAL[ii].d;
        float k,theta;
        k=InMol->DIHEDRAL[ii].k;
        mult=InMol->DIHEDRAL[ii].multiplicity;
        theta=InMol->DIHEDRAL[ii].theta;
        fprintf(ptr,"%c%-3d %c%-3d %c%-3d %c%-3d   %f  %d  %f\n",InMol->ChainId[atom1],
                atom1,InMol->ChainId[atom2],atom2,InMol->ChainId[atom3],atom3,
                InMol->ChainId[atom4],atom4,k,mult,theta);
    }
    fprintf(ptr,"\nNONBONDED NBXMOD 3 ATOM CDIEL SHIFT VATOM VDISTANCE VSWITCH -\n");
    fprintf(ptr,"  CUTNB 399.0 CTOFNB 398.5 CTONNB 395.5 EPS 1.0 WMIN 1.5\n\n");
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        int atom1;
        atom1=ii;
        float eps,vdwr;
        eps=InMol->EPS[ii];
        vdwr=InMol->VDWR[ii];
        fprintf(ptr,"%c%-3d     0.0  %f  %f\n",InMol->ChainId[atom1],atom1,eps,vdwr);
    }
    fprintf(ptr,"\nNBFIX\n");
    for(int ii=0; ii<InMol->NumOfNatives; ++ii){
        int atom1,atom2;
        float k,r;
        atom1=InMol->NATIVE[ii].i;
        atom2=InMol->NATIVE[ii].j;
        k=InMol->NATIVE[ii].k;
        r=InMol->NATIVE[ii].r;
        fprintf(ptr,"%c%-3d    %c%-3d       %f    %f\n",InMol->ChainId[atom1],atom1,
                InMol->ChainId[atom2],atom2,k,r);
    }
    if(InMol->NumOfZincBonds!=0){
        fprintf(ptr,"\nZINCBOND\n");
    }
    for(int ii=0; ii<InMol->NumOfZincBonds; ++ii){
        int atom1,atom2;
        float k,r;
        atom1=InMol->ZINCBOND[ii].i;
        atom2=InMol->ZINCBOND[ii].j;
        k=InMol->ZINCBOND[ii].k;
        r=InMol->ZINCBOND[ii].r;
        fprintf(ptr,"%c%-3d    %c%-3d       %f    %f\n",InMol->ChainId[atom1],atom1,
                InMol->ChainId[atom2],atom2,k,r);
    }
    fprintf(ptr,"\nEND\n");
#ifndef testWriteKBGoParams
    fclose(ptr);
#endif
    return 0;
}
#endif
