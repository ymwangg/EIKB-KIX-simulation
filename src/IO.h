#ifndef __IO_H__
#define __IO_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "DataType.h"
#include "DataStructure.h"
#include "StringUtility.h"
#define DEBUG_IO true
int readPDB(const char* filename, Molecule *Mol);
int renumberMolecule(Molecule *InMol);
int printMolecule(Molecule *InMol);
int calcBonds(Molecule *InMol);
/* Read the Go model topology file *.top
 */
int readCHARMMTopology(const char* filename, Molecule *OutMol);

/* Read the Go model parameter file *.param 
 */
int readCHARMMParameter(const char* filename, Molecule *OutMol);

/* Compare two trimmed strings, return 0 if two strings match each other
 * Example: " ABC   " matches "ABC   " or "     ABC"
 */

/*

 */
int initializeState(const Molecule *InMol, State *OutState);
int printVelocity(const State *InState);
int printForce(const State *InState);

int readPDB(const char* filename,Molecule *OutMol,const int model = -1){
    printf("model=%d\n",model);
    printf("Loading PDB file %s\n",filename);
    FILE *fp;
    fp=fopen(filename,"r");
    char str[MAXCHAR];
    int InputNumOfAtoms=0;
    int InputNumOfModels=0;
    if(fp==NULL){
        printf("Error: could not open file %s\n",filename);
        return 1;
    }
    //get number of atoms
    while(fgets(str,MAXCHAR,fp)!=NULL){
        if((!strncmp("ATOM  ",str,6) || !strncmp("HETATM",str,6))
                &&InputNumOfModels==1){
            InputNumOfAtoms++;
        }
        if(!strncmp("MODEL ",str,6)){
            InputNumOfModels++;
        }
    }
    if(InputNumOfModels==0){
        rewind(fp);
        while(fgets(str,MAXCHAR,fp)!=NULL){
            if(!strncmp("ATOM  ",str,6) || !strncmp("HETATM",str,6)){
                InputNumOfAtoms++;
            }
        }
    }
    printf("Number of Atoms Read = %d\n",InputNumOfAtoms);
    printf("Number of Models = %d\n",InputNumOfModels);
    OutMol->NumOfAtoms=InputNumOfAtoms;
    rewind(fp);
    //allocate memory
    OutMol->X=(float*)malloc(InputNumOfAtoms*sizeof(float));
    OutMol->Y=(float*)malloc(InputNumOfAtoms*sizeof(float));
    OutMol->Z=(float*)malloc(InputNumOfAtoms*sizeof(float));
    OutMol->tmpX=(float*)malloc(InputNumOfAtoms*sizeof(float));
    OutMol->tmpY=(float*)malloc(InputNumOfAtoms*sizeof(float));
    OutMol->tmpZ=(float*)malloc(InputNumOfAtoms*sizeof(float));
    OutMol->EPS=(float*)malloc(InputNumOfAtoms*sizeof(float));
    OutMol->Charge=(float*)malloc(InputNumOfAtoms*sizeof(float));
    OutMol->VDWR=(float*)malloc(InputNumOfAtoms*sizeof(float));
    OutMol->Mass=(float*)malloc(InputNumOfAtoms*sizeof(float));
    OutMol->AtomNum=(int*)malloc(InputNumOfAtoms*sizeof(int*));
    OutMol->ResNum=(int*)malloc(InputNumOfAtoms*sizeof(int*));
    OutMol->AtomType=(char**)malloc(InputNumOfAtoms*sizeof(char*));
    OutMol->ResType=(char**)malloc(InputNumOfAtoms*sizeof(char*));
    OutMol->ChainId=(char*)malloc(InputNumOfAtoms*sizeof(char));
    OutMol->ChainCHARMMId=(char**)malloc(InputNumOfAtoms*sizeof(char*));
    //read in properties
    //need to allocate an extra char to store terminating \0
    char x[9],y[9],z[9];
    x[8]='\0';y[8]='\0';z[8]='\0';
    char atomNum[6];
    atomNum[5]='\0';
    char resNum[5];
    resNum[4]='\0';
    char atomType[5];
    atomType[4]='\0';
    char resType[5];
    resType[4]='\0';
    char chainCHARMMId[5];
    chainCHARMMId[4]='\0';
    char chainId;
    int i=0;
    int modelNum=0;
    if((model!=-1)||InputNumOfModels==0){
        printf("Reading a single model\n");
        while(fgets(str,MAXCHAR,fp)!=NULL){
            if(!strncmp("MODEL",str,5)){
                modelNum++;
                i=0;
            }
            //if only need to read 1 model
            if(((model!=-1)&&(model==modelNum-1))||InputNumOfModels==0){
                if(!strncmp("ATOM  ",str,6) or !strncmp("HETATM",str,6)){
                    //printf("%s",str);
                    strncpy(atomNum,&str[6],5);
                    strncpy(resNum,&str[22],4);
                    strncpy(atomType,&str[12],4);
                    strncpy(resType,&str[17],4);
                    strncpy(&chainId,&str[21],1);
                    strncpy(chainCHARMMId,&str[72],4);
                    strncpy(x,&str[30],8);
                    strncpy(y,&str[38],8);
                    strncpy(z,&str[46],8);
                    OutMol->X[i]=atof(x);
                    OutMol->Y[i]=atof(y);
                    OutMol->Z[i]=atof(z);
                    OutMol->AtomNum[i]=atoi(atomNum);
                    OutMol->ResNum[i]=atoi(resNum);
                    OutMol->ChainId[i]=chainId;
                    OutMol->AtomType[i]=(char*)malloc(5*sizeof(char));
                    memcpy(OutMol->AtomType[i],atomType,4);
                    OutMol->ResType[i]=(char*)malloc(5*sizeof(char));
                    memcpy(OutMol->ResType[i],resType,4);
                    OutMol->ChainCHARMMId[i]=(char*)malloc(5*sizeof(char));
                    memcpy(OutMol->ChainCHARMMId[i],chainCHARMMId,4);
                    i++;
                }
            } 
        }
    }
    else{
        printf("Averaging over multiple models\n");
        while(fgets(str,MAXCHAR,fp)!=NULL){
            if(!strncmp("MODEL",str,5)){
                modelNum++;
                i=0;
            }
            //if only need to read 1 model
            //average over all models
            if(!strncmp("ATOM  ",str,6) or !strncmp("HETATM",str,6)){
                if(modelNum==1){
                    //printf("%s",str);
                    memcpy(atomNum,&str[6],5);
                    memcpy(resNum,&str[22],4);
                    memcpy(atomType,&str[12],4);
                    memcpy(resType,&str[17],4);
                    memcpy(&chainId,&str[21],1);
                    memcpy(chainCHARMMId,&str[72],4);
                    memcpy(x,&str[30],8);
                    memcpy(y,&str[38],8);
                    memcpy(z,&str[46],8);
                    OutMol->X[i]=atof(x)/InputNumOfModels;
                    OutMol->Y[i]=atof(y)/InputNumOfModels;
                    OutMol->Z[i]=atof(z)/InputNumOfModels;
#if DEBUG_IO
                    if(i==0){
                        printf("(%f,%f,%f)/%d\n",OutMol->X[i],OutMol->Y[i],OutMol->Z[i],InputNumOfModels);
                    }
#endif
                    OutMol->AtomNum[i]=atoi(atomNum);
                    OutMol->ResNum[i]=atoi(resNum);
                    OutMol->ChainId[i]=chainId;
                    OutMol->AtomType[i]=(char*)malloc(5*sizeof(char));
                    memcpy(OutMol->AtomType[i],atomType,4);
                    OutMol->ResType[i]=(char*)malloc(5*sizeof(char));
                    memcpy(OutMol->ResType[i],resType,4);
                    OutMol->ChainCHARMMId[i]=(char*)malloc(5*sizeof(char));
                    memcpy(OutMol->ChainCHARMMId[i],chainCHARMMId,4);
                    i++;
                }else{
                    memcpy(x,&str[30],8);
                    memcpy(y,&str[38],8);
                    memcpy(z,&str[46],8);
                    OutMol->X[i]+=atof(x)/InputNumOfModels;
                    OutMol->Y[i]+=atof(y)/InputNumOfModels;
                    OutMol->Z[i]+=atof(z)/InputNumOfModels;
#if DEBUG_IO
                    if(i==0){
                        printf("(%f,%f,%f)/%d\n",OutMol->X[i],
                                OutMol->Y[i],OutMol->Z[i],InputNumOfModels);
                    }
#endif
                    i++;
                }
            }

        }
    }
    printf("Finish reading PDB file %s\n",filename);
    fclose(fp);
    return 0;
}
int renumberMolecule(Molecule *InMol){
    int atomNum=0;
    int resNum=0;
    int currentResNum=InMol->ResNum[0];
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        InMol->AtomNum[ii]=atomNum;
        atomNum++;
        if(currentResNum!=InMol->ResNum[ii]){
            resNum++;
            currentResNum=InMol->ResNum[ii];
        }
        InMol->ResNum[ii]=resNum;
    }
    InMol->NumOfResidues=resNum+1;
    if(isblank(InMol->ChainId[0])){
        char chainId='A';
        char *currentChain=InMol->ChainCHARMMId[0];
        for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
            if(mystrcmp(InMol->ChainCHARMMId[ii],currentChain)!=0){
                currentChain=InMol->ChainCHARMMId[ii];
                chainId++;
            }
            InMol->ChainId[ii]=chainId;
        }
    }
    return 0;
}
int printMolecule(Molecule *InMol){
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        printf("Atom%d %f %f %f\n",ii,InMol->X[ii],InMol->Y[ii],InMol->Z[ii]);
        printf("Atom%d atomNum=%d resNum=%d\n",ii,InMol->AtomNum[ii],
                InMol->ResNum[ii]);
        printf("Atom%d atomType='%s' resType='%s' chainId=%c chainCHARMMId='%s'\n",
                ii,InMol->AtomType[ii],InMol->ResType[ii],InMol->ChainId[ii],
                InMol->ChainCHARMMId[ii]);
        //printf("Atom%d %f %f %f %d %d %s %s %s %s\n",ii,InMol->X[ii],InMol->Y[ii],InMol->Z[ii],
        //        InMol->AtomNum[ii],InMol->ResNum[ii],InMol->AtomType[ii],InMol->ResType[ii],
        //        InMol->ChainId[ii],InMol->ChainCHARMMId[ii]);
    }
    return 0;
}
int readCHARMMTopology(const char* filename, Molecule *OutMol){
    printf("Reading topology file %s\n",filename);
    FILE *fp;
    fp=fopen(filename,"r");
    char str[MAXCHAR];
    int InputNumOfAtoms=0;
    if(fp==NULL){
        printf("Error: could not open file %s\n",filename);
        return 1;
    }
    //get number of atoms
    while(fgets(str,MAXCHAR,fp)!=NULL){
        if(!strncmp("MASS",str,4)){
            InputNumOfAtoms++;
        }
    }
    if(InputNumOfAtoms!=OutMol->NumOfAtoms){
        printf("Warning: the number of atoms in topology file didn't match PDB file\n");
    }
    rewind(fp);
    char *tmp,*atomNum,*resType,*mass,*charge;
    char buffer[4];
    int i=0;
    int j=0;
    while(fgets(str,MAXCHAR,fp)!=NULL){
        if(!strncmp("MASS",str,4)){
           tmp=strtok(str," \n");
           atomNum=strtok(NULL," \n");
           resType=strtok(NULL," \n");
           mass=strtok(NULL," \n");
           sprintf(buffer,"%s",resType);
           memcpy(OutMol->ResType[i],buffer,4);
           OutMol->Mass[i]=atof(mass);
#if DEBUG_IO
           printf("restype=%s mass=%f\n",OutMol->ResType[i],OutMol->Mass[i]);
#endif
           i++;
        }
        if(!strncmp("RESI",str,4)){
            tmp=strtok(str," \n");
            resType=strtok(NULL," \n");
            charge=strtok(NULL," \n");
            OutMol->Charge[j]=atof(charge);
#if DEBUG_IO
            printf("restype=%s charge=%f\n",resType,OutMol->Charge[j]);
#endif
            j++;
        }
    }
    printf("Number of atoms = %d\n",InputNumOfAtoms);
    printf("Finish reading topology file %s\n",filename);
    fclose(fp);
    return 0;
}
int readCHARMMParameter(const char* filename, Molecule *OutMol){
    printf("Reading parameter file %s\n",filename);
    FILE *fp;
    fp=fopen(filename,"r");
    char str[MAXCHAR];
    if(fp==NULL){
        printf("Error: could not open file %s\n",filename);
        return 1;
    }
    //get number of bonds
    bool rightField=false;
    int numOfBonds=0,numOfAngles=0,numOfDihedrals=0,numOfNatives=0,numOfZincBonds=0;
    int numOfAtoms=0;
    while(fgets(str,MAXCHAR,fp)!=NULL){
        //if(regex_match(str,regex("[b|B][o|O][n|N][d|D]"))){
        if(rightField && (strncmp(str," ",1)==0 || strncmp(str,"\n",1)==0)){
            break;
        }
        if(rightField){
            numOfBonds++;
        }
        if(!strncmp("BOND",str,4)){
            rightField=true;
        }
    }
    //get number of angles
    rightField=false;
    while(fgets(str,MAXCHAR,fp)!=NULL){
        //if(regex_match(str,regex("[b|B][o|O][n|N][d|D]"))){
        if(rightField && (strncmp(str," ",1)==0 || strncmp(str,"\n",1)==0)){
            break;
        }
        if(rightField){
            numOfAngles++;
        }
        if(!strncmp("ANGL",str,4)){
            rightField=true;
        }
    }
    //get number of dihedrals
    rightField=false;
    while(fgets(str,MAXCHAR,fp)!=NULL){
        //if(regex_match(str,regex("[b|B][o|O][n|N][d|D]"))){
        if(rightField && (strncmp(str," ",1)==0 || strncmp(str,"\n",1)==0)){
            break;
        }
        if(rightField){
            numOfDihedrals++;
        }
        if(!strncmp("DIHE",str,4)){
            rightField=true;
        }
    }
    //get number of atoms
    rightField=false;
    while(fgets(str,MAXCHAR,fp)!=NULL){
        //if(regex_match(str,regex("[b|B][o|O][n|N][d|D]"))){
        if(rightField && (strncmp(str," ",1)==0 || strncmp(str,"\n",1)==0)){
            break;
        }
        if(rightField){
            numOfAtoms++;
        }
        if(!strncmp("NONB",str,4)){
            //printf("%s",str);
            rightField=true;
            fgets(str,MAXCHAR,fp);
            //printf("%s",str);
            fgets(str,MAXCHAR,fp);
            //printf("%s",str);
        }
    }
    if(numOfAtoms!=OutMol->NumOfAtoms){
        printf("Warning: the number of atoms in topology file didn't match PDB file\n");
    }
    //get number of nbfix
    rightField=false;
    while(fgets(str,MAXCHAR,fp)!=NULL){
        //if(regex_match(str,regex("[b|B][o|O][n|N][d|D]"))){
        if(rightField && (strncmp(str," ",1)==0 || strncmp(str,"\n",1)==0)){
            break;
        }
        if(rightField){
            numOfNatives++;
        }
        if(!strncmp("NBFI",str,4)){
            rightField=true;
            //printf("%s",str);
        }
    }
    //get number of zincBond
    rightField=false;
    while(fgets(str,MAXCHAR,fp)!=NULL){
        //if(regex_match(str,regex("[b|B][o|O][n|N][d|D]"))){
        if(rightField && (strncmp(str," ",1)==0 || strncmp(str,"\n",1)==0)){
            break;
        }
        if(rightField){
            numOfZincBonds++;
        }
        if(!strncmp("ZINC",str,4)){
            rightField=true;
        }
    }
    printf("Num of bonds=%d\n",numOfBonds);
    printf("Num of angles=%d\n",numOfAngles);
    printf("Num of dihedrals=%d\n",numOfDihedrals);
    printf("Num of atoms=%d\n",numOfAtoms);
    printf("Num of nbfixs=%d\n",numOfNatives);
    printf("Num of zincBonds=%d\n",numOfZincBonds);
    //read parameters
    rewind(fp);
    OutMol->NumOfBonds=numOfBonds;
    OutMol->BOND=(Bond*)malloc(numOfBonds*sizeof(Bond));
    OutMol->NumOfAngles=numOfAngles;
    OutMol->ANGLE=(Angle*)malloc(numOfAngles*sizeof(Angle));
    OutMol->NumOfDihedrals=numOfDihedrals;
    OutMol->DIHEDRAL=(Dihedral*)malloc(numOfDihedrals*sizeof(Dihedral));
    OutMol->NumOfNatives=numOfNatives;
    OutMol->NATIVE=(Native*)malloc(numOfNatives*sizeof(Native));
    OutMol->NumOfZincBonds=numOfZincBonds;
    OutMol->ZINCBOND=(Bond*)malloc(numOfZincBonds*sizeof(Bond));
    //read bonds
    int i=0;
    rightField=false;
    while(fgets(str,MAXCHAR,fp)!=NULL){
        //if(regex_match(str,regex("[b|B][o|O][n|N][d|D]"))){
        char *atom1,*atom2,*k,*r;
        if(rightField && (strncmp(str," ",1)==0 || strncmp(str,"\n",1)==0)){
            break;
        }
        if(rightField){
           atom1=strtok(str," \n");
           atom2=strtok(NULL," \n");
           k=strtok(NULL," \n");
           r=strtok(NULL," \n");
           for(int ii=0; ii<OutMol->NumOfAtoms; ++ii){
               if(mystrcmp(atom1,OutMol->ResType[ii])==0){
                   OutMol->BOND[i].i=ii;
               }
               if(mystrcmp(atom2,OutMol->ResType[ii])==0){
                   OutMol->BOND[i].j=ii;
               }
           }
           OutMol->BOND[i].k=atof(k);
           OutMol->BOND[i].r=atof(r);
#if DEBUG_IO
           printf("Bond%d (%d,%d) k=%f r=%f\n",i,OutMol->BOND[i].i,OutMol->BOND[i].j,
                   OutMol->BOND[i].k,OutMol->BOND[i].r);
#endif
           //printf("bond%d %s %s %f %f\n",i,atom1,atom2,
           //        atof(k),atof(r));
           i++;
        }
        if(!strncmp("BOND",str,4)){
            rightField=true;
        }
    }
    //get number of angles
    i=0;
    rightField=false;
    while(fgets(str,MAXCHAR,fp)!=NULL){
        //if(regex_match(str,regex("[b|B][o|O][n|N][d|D]"))){
        char *atom1,*atom2,*atom3,*k,*theta;
        if(rightField && (strncmp(str," ",1)==0 || strncmp(str,"\n",1)==0)){
            break;
        }
        if(rightField){
            atom1=strtok(str," \n");
            atom2=strtok(NULL," \n");
            atom3=strtok(NULL," \n");
            k=strtok(NULL," \n");
            theta=strtok(NULL," \n");
            for(int ii=0; ii<OutMol->NumOfAtoms; ++ii){
                if(mystrcmp(atom1,OutMol->ResType[ii])==0){
                    OutMol->ANGLE[i].a=ii;
                }
                if(mystrcmp(atom2,OutMol->ResType[ii])==0){
                    OutMol->ANGLE[i].b=ii;
                }
                if(mystrcmp(atom3,OutMol->ResType[ii])==0){
                    OutMol->ANGLE[i].c=ii;
                }
            }
            OutMol->ANGLE[i].k=atof(k);
            OutMol->ANGLE[i].theta=atof(theta);
#if DEBUG_IO
            printf("Angle%d (%d,%d,%d) k=%f theta=%f\n",i,OutMol->ANGLE[i].a,
                    OutMol->ANGLE[i].b,OutMol->ANGLE[i].c,OutMol->ANGLE[i].k,
                    OutMol->ANGLE[i].theta);
#endif
            i++;
        }
        if(!strncmp("ANGL",str,4)){
            rightField=true;
        }
    }
    //get number of dihedrals
    i=0;
    rightField=false;
    while(fgets(str,MAXCHAR,fp)!=NULL){
        //if(regex_match(str,regex("[b|B][o|O][n|N][d|D]"))){
        char *atom1,*atom2,*atom3,*atom4,*k,*mult,*theta;
        if(rightField && (strncmp(str," ",1)==0 || strncmp(str,"\n",1)==0)){
            break;
        }
        if(rightField){
            atom1=strtok(str," \n");
            atom2=strtok(NULL," \n");
            atom3=strtok(NULL," \n");
            atom4=strtok(NULL," \n");
            k=strtok(NULL," \n");
            mult=strtok(NULL," \n");
            theta=strtok(NULL," \n");
            for(int ii=0; ii<OutMol->NumOfAtoms; ++ii){
                //printf("%s\n",OutMol->ResType[ii]);
                if(mystrcmp(atom1,OutMol->ResType[ii])==0){
                    OutMol->DIHEDRAL[i].a=ii;
                }
                if(mystrcmp(atom2,OutMol->ResType[ii])==0){
                    OutMol->DIHEDRAL[i].b=ii;
                }
                if(mystrcmp(atom3,OutMol->ResType[ii])==0){
                    OutMol->DIHEDRAL[i].c=ii;
                }
                if(mystrcmp(atom4,OutMol->ResType[ii])==0){
                    OutMol->DIHEDRAL[i].d=ii;
                }
            }
            OutMol->DIHEDRAL[i].k=atof(k);
            OutMol->DIHEDRAL[i].multiplicity=atoi(mult);
            OutMol->DIHEDRAL[i].theta=atof(theta);
#if DEBUG_IO
            printf("Dihedral%d (%d,%d,%d,%d) k=%f mult=%d theta=%f\n",i,
                    OutMol->DIHEDRAL[i].a,
                    OutMol->DIHEDRAL[i].b,OutMol->DIHEDRAL[i].c,
                    OutMol->DIHEDRAL[i].d,OutMol->DIHEDRAL[i].k,
                    OutMol->DIHEDRAL[i].multiplicity,OutMol->DIHEDRAL[i].theta);
#endif
            i++;
        }
        if(!strncmp("DIHE",str,4)){
            rightField=true;
        }
    }
    //get number of atoms
    i=0;
    rightField=false;
    while(fgets(str,MAXCHAR,fp)!=NULL){
        //if(regex_match(str,regex("[b|B][o|O][n|N][d|D]"))){
        char *atom, *tmp, *eps, *vdwr;
        if(rightField && (strncmp(str," ",1)==0 || strncmp(str,"\n",1)==0)){
            break;
        }
        if(rightField){
           atom=strtok(str," \n");
           tmp=strtok(NULL," \n");
           eps=strtok(NULL," \n");
           vdwr=strtok(NULL," \n");
           OutMol->EPS[i]=atof(eps);
           OutMol->VDWR[i]=atof(vdwr);
#if DEBUG_IO
           printf("Atom%d EPS=%f VDWR=%f\n",i,OutMol->EPS[i],OutMol->VDWR[i]);
#endif
           //printf("bond%d %s %s %f %f\n",i,atom1,atom2,
           //        atof(k),atof(r));
           i++;
        }
        if(!strncmp("NONB",str,4)){
            //printf("%s",str);
            rightField=true;
            fgets(str,MAXCHAR,fp);
            //printf("%s",str);
            fgets(str,MAXCHAR,fp);
            //printf("%s",str);
        }
    }
    //get number of nbfix
    i=0;
    rightField=false;
    while(fgets(str,MAXCHAR,fp)!=NULL){
        //if(regex_match(str,regex("[b|B][o|O][n|N][d|D]"))){
        char *atom1,*atom2,*k,*r;
        if(rightField && (strncmp(str," ",1)==0 || strncmp(str,"\n",1)==0)){
            break;
        }
        if(rightField){
           atom1=strtok(str," \n");
           atom2=strtok(NULL," \n");
           k=strtok(NULL," \n");
           r=strtok(NULL," \n");
           for(int ii=0; ii<OutMol->NumOfAtoms; ++ii){
               if(mystrcmp(atom1,OutMol->ResType[ii])==0){
                   OutMol->NATIVE[i].i=ii;
               }
               if(mystrcmp(atom2,OutMol->ResType[ii])==0){
                   OutMol->NATIVE[i].j=ii;
               }
           }
           OutMol->NATIVE[i].k=atof(k);
           OutMol->NATIVE[i].r=atof(r);
#if DEBUG_IO
           printf("Native%d (%d,%d) k=%f r=%f\n",i,OutMol->NATIVE[i].i,OutMol->NATIVE[i].j,
                   OutMol->NATIVE[i].k,OutMol->NATIVE[i].r);
#endif
           //printf("bond%d %s %s %f %f\n",i,atom1,atom2,
           //        atof(k),atof(r));
           i++;
        }
        if(!strncmp("NBFI",str,4)){
            rightField=true;
            //printf("%s",str);
        }
    }

    //get number of zincbonds
    i=0;
    rightField=false;
    while(fgets(str,MAXCHAR,fp)!=NULL){
        //if(regex_match(str,regex("[b|B][o|O][n|N][d|D]"))){
        char *atom1,*atom2,*k,*r;
        if(rightField && (strncmp(str," ",1)==0 || strncmp(str,"\n",1)==0)){
            break;
        }
        if(rightField){
           atom1=strtok(str," \n");
           atom2=strtok(NULL," \n");
           k=strtok(NULL," \n");
           r=strtok(NULL," \n");
           for(int ii=0; ii<OutMol->NumOfAtoms; ++ii){
               if(mystrcmp(atom1,OutMol->ResType[ii])==0){
                   OutMol->ZINCBOND[i].i=ii;
               }
               if(mystrcmp(atom2,OutMol->ResType[ii])==0){
                   OutMol->ZINCBOND[i].j=ii;
               }
           }
           OutMol->ZINCBOND[i].k=atof(k);
           OutMol->ZINCBOND[i].r=atof(r);
#if DEBUG_IO
           printf("ZincBond%d (%d,%d) k=%f r=%f\n",i,OutMol->ZINCBOND[i].i,OutMol->ZINCBOND[i].j,
                   OutMol->ZINCBOND[i].k,OutMol->ZINCBOND[i].r);
#endif
           //printf("bond%d %s %s %f %f\n",i,atom1,atom2,
           //        atof(k),atof(r));
           i++;
        }
        if(!strncmp("ZINC",str,4)){
            rightField=true;
            //printf("%s",str);
        }
    }

    fclose(fp);
    return 0;
}
int initializeState(const Molecule *InMol, State *OutState){
    if(InMol->NumOfAtoms==0){
        return 1;
    }   
    int numOfAtoms=InMol->NumOfAtoms;
    OutState->NumOfAtoms=numOfAtoms;
    OutState->X=(float*)malloc(numOfAtoms*sizeof(float));
    OutState->Y=(float*)malloc(numOfAtoms*sizeof(float));
    OutState->Z=(float*)malloc(numOfAtoms*sizeof(float));
    OutState->VX=(float*)malloc(numOfAtoms*sizeof(float));
    OutState->VY=(float*)malloc(numOfAtoms*sizeof(float));
    OutState->VZ=(float*)malloc(numOfAtoms*sizeof(float));
    OutState->FX=(float*)malloc(numOfAtoms*sizeof(float));
    OutState->FY=(float*)malloc(numOfAtoms*sizeof(float));
    OutState->FZ=(float*)malloc(numOfAtoms*sizeof(float));
    for(int ii=0; ii<numOfAtoms; ++ii){
        OutState->X[ii]=InMol->X[ii];
        OutState->Y[ii]=InMol->Y[ii];
        OutState->Z[ii]=InMol->Z[ii];
        OutState->VX[ii]=0.0;
        OutState->VY[ii]=0.0;
        OutState->VZ[ii]=0.0;
        OutState->FX[ii]=0.0;
        OutState->FY[ii]=0.0;
        OutState->FZ[ii]=0.0;
        printf("%f %f %f\n",OutState->X[ii],OutState->Y[ii],OutState->Z[ii]);
    }   
    return 0;
}
int printForce(const State *InState){
    int numOfAtoms=InState->NumOfAtoms;
    printf("********Force********\n");
    for(int ii=0; ii<numOfAtoms; ++ii){
        printf("Atom%d (%f,%f,%f)\n",ii,InState->FX[ii],
                InState->FY[ii],InState->FZ[ii]);
    }   
    printf("******END-Force******\n");
    return 0;
}
int printVelocity(const State *InState){
    int numOfAtoms=InState->NumOfAtoms;
    printf("********Velocity********\n");
    for(int ii=0; ii<numOfAtoms; ++ii){
        printf("Atom%d (%f,%f,%f)\n",ii,InState->VX[ii],
                InState->VY[ii],InState->VZ[ii]);
    }   
    printf("******END-Velocity******\n");
    return 0;
}
int printChargePerChain(const Molecule *InMol){
    std::vector<float> charges;
    char chainId=InMol->ChainId[0];
    float sum=0;
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        if(chainId!=InMol->ChainId[ii]){
            charges.push_back(sum);
            sum=0;
            chainId=InMol->ChainId[ii];
        }
        sum+=InMol->Charge[ii];
    }
    charges.push_back(sum);
    int i=0;
    for(float charge : charges){
        printf("Chain%d charge=%f\n",i,charge);
        i++;
    }
    return 0;
}
#endif
