#ifndef __COORDINATE_H__
#define __COORDINATE_H__
#include "DataStructure.h"
#include "DataType.h"
int getCenterOfMass(const Molecule *InMol,float *coor);
int getCenterOfGeometry(const Molecule *InMol,float *coor);
int setMoleculeToPosition(Molecule *InMol,float bx, float by, float bz);
int setChainToPosition(Molecule *InMol,const char chainId, float bx, float by,
        float bz);
int setMultChainToPosition(Molecule *InMol,const char chainId, float bx, float by,
        float bz);
int copyCoordinate(const Molecule *InMol,Molecule *OutMol);
int copyCoordinateByChain(const Molecule *InMol, Molecule *OutMol, const char chainId);
int getCenterOfMassByChain(const Molecule *InMol,const char chainId,float *coor);
int getCenterOfGeometryByChain(const Molecule *InMol,const char chainId, float *coor);
int translateByChain(Molecule *InMol,const char chainId, const float dx, const float dy, const float dz);
int moveBackToOriginalBox(Molecule *InMol, const float bx, const float by, const float bz);

int getCenterOfMass(const Molecule *InMol,float *coor){
    int numOfAtoms=InMol->NumOfAtoms;
    float x,y,z;
    x=0;
    y=0;
    z=0;
    for(int ii=0; ii<numOfAtoms; ++ii){
        x+=InMol->X[ii]*InMol->Mass[ii];
        y+=InMol->Y[ii]*InMol->Mass[ii];
        z+=InMol->Z[ii]*InMol->Mass[ii];
    }
    float totalMass=0;
    for(int ii=0; ii<numOfAtoms; ++ii){
        totalMass+=InMol->Mass[ii];
    }
    x=x/(totalMass);
    y=y/(totalMass);
    z=z/(totalMass);
    coor[0]=x;
    coor[1]=y;
    coor[2]=z;
    return 0;
}
int getCenterOfGeometry(const Molecule *InMol,float *coor){
    int numOfAtoms=InMol->NumOfAtoms;
    float x,y,z;
    x=0;
    y=0;
    z=0;
    for(int ii=0; ii<numOfAtoms; ++ii){
        x+=InMol->X[ii];
        y+=InMol->Y[ii];
        z+=InMol->Z[ii];
    }
    x=x/(numOfAtoms);
    y=y/(numOfAtoms);
    z=z/(numOfAtoms);
    coor[0]=x;
    coor[1]=y;
    coor[2]=z;
    return 0;
}
int setMoleculeToPosition(Molecule *InMol,float bx, float by, float bz){
    int numOfAtoms=InMol->NumOfAtoms;
    float offset[3];
    getCenterOfMass(InMol,offset);
    for(int ii=0; ii<numOfAtoms; ++ii){
        InMol->X[ii]+=bx-offset[0];
        InMol->Y[ii]+=by-offset[1];
        InMol->Z[ii]+=bz-offset[2];

    }
    return 0;
}
int setChainToPosition(Molecule *InMol,const char chainId, float bx, float by,
        float bz){
    int numOfAtoms=InMol->NumOfAtoms;
    float offset[3];
    getCenterOfMassByChain(InMol,chainId,offset);
    for(int ii=0; ii<numOfAtoms; ++ii){
        if(InMol->ChainId[ii]==chainId){
            InMol->X[ii]+=bx-offset[0];
            InMol->Y[ii]+=by-offset[1];
            InMol->Z[ii]+=bz-offset[2];
        }
    }
    return 0;
}
int setMultChainToPosition(Molecule *InMol,const std::vector<char> &chains, float bx, float by,float bz){
    int numOfAtoms=InMol->NumOfAtoms;
    float offset[3];
    getCenterOfMassByChain(InMol,chains[0],offset);
    for(int ii=0; ii<numOfAtoms; ++ii){
        for(char chainId : chains){
            if(InMol->ChainId[ii]==chainId){
                InMol->X[ii]+=bx-offset[0];
                InMol->Y[ii]+=by-offset[1];
                InMol->Z[ii]+=bz-offset[2];
            }
        }
    }
    return 0;
}
int copyCoordinate(const Molecule *InMol,Molecule *OutMol){
    if(InMol->NumOfAtoms!=OutMol->NumOfAtoms){
        return 1;
    }else{
        for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
            OutMol->X[ii]=InMol->X[ii];
            OutMol->Y[ii]=InMol->Y[ii];
            OutMol->Z[ii]=InMol->Z[ii];
        }
    }
    return 0;
}

int copyCoordinateByChain(const Molecule *InMol, const char chainId1, Molecule *OutMol, const char chainId2){
    int index1,index2;
    int numOfAtoms1, numOfAtoms2;
    numOfAtoms1=0;
    numOfAtoms2=0;
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        if(InMol->ChainId[ii]==chainId1){
            index1=ii;
            break;
        }
    }
    for(int ii=0; ii<OutMol->NumOfAtoms; ++ii){
        if(OutMol->ChainId[ii]==chainId2){
            index2=ii;
            break;
        }
    }
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        if(InMol->ChainId[ii]==chainId1){
            numOfAtoms1++;
        }
    }
    for(int ii=0; ii<OutMol->NumOfAtoms; ++ii){
        if(OutMol->ChainId[ii]==chainId2){
            numOfAtoms2++;
        }
    }
    if(numOfAtoms1==numOfAtoms2){
        for(int ii=0; ii<numOfAtoms1; ++ii){
            OutMol->X[index2+ii]=InMol->X[index1+ii];
            OutMol->Y[index2+ii]=InMol->Y[index1+ii];
            OutMol->Z[index2+ii]=InMol->Z[index1+ii];
        }
    }else{
        return 1;
    }
    return 0;
}
int getCenterOfMassByChain(const Molecule *InMol,const char chainId,float *coor){
    int numOfAtoms=InMol->NumOfAtoms;
    float x,y,z;
    x=0;
    y=0;
    z=0;
    float totalMass=0;
    for(int ii=0; ii<numOfAtoms; ++ii){
        if(InMol->ChainId[ii]==chainId){
            x+=InMol->X[ii]*InMol->Mass[ii];
            y+=InMol->Y[ii]*InMol->Mass[ii];
            z+=InMol->Z[ii]*InMol->Mass[ii];
            totalMass+=InMol->Mass[ii];
        }
    }
    x=x/(totalMass);
    y=y/(totalMass);
    z=z/(totalMass);
    coor[0]=x;
    coor[1]=y;
    coor[2]=z;
    return 0;
}
int getCenterOfGeometryByChain(const Molecule *InMol,const char chainId, float *coor){
    int numOfAtoms=InMol->NumOfAtoms;
    float x,y,z;
    x=0;
    y=0;
    z=0;
    int numOfSelectedAtoms=0;
    for(int ii=0; ii<numOfAtoms; ++ii){
        if(InMol->ChainId[ii]==chainId){
            x+=InMol->X[ii];
            y+=InMol->Y[ii];
            z+=InMol->Z[ii];
            numOfSelectedAtoms++;
        }
    }
    x=x/(numOfSelectedAtoms);
    y=y/(numOfSelectedAtoms);
    z=z/(numOfSelectedAtoms);
    coor[0]=x;
    coor[1]=y;
    coor[2]=z;
    return 0;
}
int translateByChain(Molecule *InMol,const char chainId, const float dx, const float dy, const float dz){
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        if(InMol->ChainId[ii]==chainId){
           InMol->X[ii]+=dx;
           InMol->Y[ii]+=dy;
           InMol->Z[ii]+=dz;
        }
    }
    return 0;
}
int moveBackToOriginalBox(Molecule *InMol, const float bx, const float by, const float bz){
    float coor[3];
    float newcoor[3];
    std::vector<char> chainIds;
    char currentChainId=InMol->ChainId[0];
    chainIds.push_back(currentChainId);
    for(int ii=1; ii<InMol->NumOfAtoms; ++ii){
        if(currentChainId!=InMol->ChainId[ii]){
            chainIds.push_back(InMol->ChainId[ii]);
            currentChainId=InMol->ChainId[ii];
        }
    }
    for(char chainId : chainIds){
        //printf("coor fixing chain %c\n",chainId);
        getCenterOfMassByChain(InMol,chainId,coor);
        //printf("old COM = (%f,%f,%f)\n",coor[0],coor[1],coor[2]);
        //fix negative
        while(coor[0]<0){
            coor[0]+=bx;
        }
        while(coor[1]<0){
            coor[1]+=by;
        }
        while(coor[2]<0){
            coor[2]+=bz;
        }
        newcoor[0]= fmod(coor[0],bx);
        newcoor[1]= fmod(coor[1],by);
        newcoor[2]= fmod(coor[2],bz);
        //printf("new COM = (%f,%f,%f)\n",newcoor[0],newcoor[1],newcoor[2]);
        setChainToPosition(InMol,chainId,newcoor[0],newcoor[1],newcoor[2]);
    }
    return 0;
}
#endif
