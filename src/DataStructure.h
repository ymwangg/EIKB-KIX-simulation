#ifndef __DATA_STRUCTURE_H__
#define __DATA_STRUCTURE_H__
#include "DataType.h"
struct Molecule;
struct Bond;
struct Angle;
struct Dihedral;
struct Native;
struct State;
struct Molecule{
    Molecule(): nativeScalingFactor(0),NumOfBonds(0),
    NumOfAngles(0),NumOfDihedrals(0),NumOfNatives(0),
    NumOfZincBonds(0),NumOfZincAngles(0){};
    unsigned int NumOfAtoms;
    unsigned int NumOfResidues;
    float nativeScalingFactor;
    float *X;
    float *Y;
    float *Z;
    float *tmpX;
    float *tmpY;
    float *tmpZ;
    float *EPS;
    float *VDWR;
    float *Charge;
    float *Mass;
    int *AtomNum; //7-11 col in PDB
    char **AtomType; //13-16 col in PDB 5*sizeof(char)
    char **ResType; //17-20 col in PDB 4*sizeof(char)
    char *ChainId; //22 col in PDB 1*sizeof(char)
    char **ChainCHARMMId; //73-76 in PDB 5*sizeof(char)
    int *ResNum; //23-26 col in PDB
    unsigned int NumOfBonds;
    Bond *BOND;
    unsigned int NumOfAngles;
    Angle *ANGLE;
    unsigned int NumOfDihedrals;
    Dihedral *DIHEDRAL;
    unsigned int NumOfNatives;
    Native *NATIVE;
    unsigned int NumOfZincBonds;
    Bond *ZINCBOND;
    unsigned int NumOfZincAngles;
    Angle *ZINCANGLE;
};
struct Bond{
    int i,j;
    float k;
    float r;
};
struct Angle{
    int a,b,c;
    float k;
    float theta;
};
struct Dihedral{
    int a,b,c,d;
    float k;
    int multiplicity;
    float theta;
};
struct Native{
    int i,j;
    float k;
    float r;
};
struct State{
    unsigned int NumOfAtoms;
    float *X;
    float *Y;
    float *Z;
    float *VX;
    float *VY;
    float *VZ;
    float *FX;
    float *FY;
    float *FZ;
};
#endif
