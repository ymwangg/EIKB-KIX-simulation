#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <regex>
#include "DataType.h"
#include "DataStructure.h"
#include "IO.h"
#include "Bonded.h"
#include "Vector.h"
#include "OpenMMDataStructure.h"
#include "OpenMMKernel.h"
#include "Coordinate.h"
#include "GoModelBuilder.h"
#include "GoModelBuilderIO.h"
#include "OpenMMRepd.h"
using namespace::std;

int doMD(Molecule *InMol, OpenMMParameters* InParams, OpenMMData* InOpenMM);
int doTemperatureRepd(Molecule *InMol,OpenMMParameters* InParams, 
        OpenMMData* InOpenMM);
int doHamiltonianRepd(Molecule *InMol,OpenMMParameters* InParams, 
        OpenMMData* InOpenMM);
int doHamiltonianRepd3D(Molecule *InMol,OpenMMParameters* InParams, 
        OpenMMData* InOpenMM);

int main(int argc, char** argv){
    Molecule mol;
    for(int ii=0; ii<argc; ++ii){
        printf("args[%d] = %s\n",ii,argv[ii]);
    }
    readPDB(argv[1],&mol,-1);
    renumberMolecule(&mol);
    Molecule gomol;
    //buildGoModel(&mol,&gomol);
    readCHARMMTopology(argv[2],&mol);
    readCHARMMParameter(argv[3],&mol);
    //float scaleAB = atof(argv[4]);

    copyGoMolecule(&mol,&gomol);
    printChargePerChain(&gomol);

    OpenMMData ommData;
    OpenMMParameters ommParams;
    //ommParams.timestep=0.022;

    /*TAZ1*/
    //chainA:TAZ1, chainB:HIF, chainC:CITED2
    //TAZ1-HIF
    //scaleNative(&gomol,{{'B','B',0.1}}); //HIF
    //TAZ1-CITED2
    //scaleNative(&gomol,{{'B','B',0.1}}); //CITED2
    //
    /*TAZ1-ternary*/
    /*KIX*/
    //KIX-MLL
    scaleNative(&gomol,{{'B','B',0.1}}); //MLL
    //KIX-pKID
    //scaleNativepKID1(&gomol,{{'B','B',0.7}}); //pKID alpha1
    //scaleNativepKID2(&gomol,{{'B','B',0.3}}); //pKID alpha2
    /*KIX-ternary*/
    /*
    scaleNative(&gomol,{{'B','B',0.05}}); //mll
    scaleNativepKID2(&gomol,{{'C','C',0.6}}); //pKID
    scaleNative(&gomol,{{'A','B',1.360604}}); //KIX-mll
    scaleNative(&gomol,{{'A','C',1.891507}}); //KIX-pKID
    */
    //hif beta=1.516053
    //cited2 beta=1.200999
    scaleNative(&gomol,{{'B','B',0.1}});
    scaleNative(&gomol,{{'C','C',0.1}});
    scaleNative(&gomol,{{'A','B',1.516053}});
    scaleNative(&gomol,{{'A','C',1.200999}});
    scaleNative(&gomol,{{'B','B',0.4}}); //c-myb
    scaleNative(&gomol,{{'C','C',0.1}}); //mll
    scaleNative(&gomol,{{'A','B',1.097403}}); //KIX-c-myb
    //scaleNative(&gomol,{{'A','C',1.394810}}); //KIX-MLL

    MPI_Init(NULL,NULL);
    setupOpenMMSystem(&gomol,&ommParams,&ommData);

    //do simulations
    doMD(&gomol, &ommParams, &ommData);
    //doTemperatureRepd(&gomol, &ommParams, &ommData);
    //doHamiltonianRepd(&gomol, &ommParams, &ommData);

    return 0;
}

int doMD(Molecule *InMol, OpenMMParameters* InParams, OpenMMData* InOpenMM){
    InParams->usesPeriodicBoundaryConditions=true;
    InParams->bx=150;
    InParams->by=150;
    InParams->bz=150;
    InParams->outDCDName="1.dcd";
    //InParams->steps=400000000;
    InParams->steps=40000000;
    InParams->printFreq=10000;
    InParams->temperature=300;
    //setMultChainToPosition(InMol,{'A','B'},InParams->bx/2.0,InParams->by/2.0,InParams->bz/2.0);
    setMultChainToPosition(InMol,{'A','B','C'}, InParams->bx/2.0, InParams->by/2.0, InParams->bz/2.0);
    //setChainToPosition(InMol,'A',InParams->bx/2.0,InParams->by/2.0,InParams->bz/2.0);
    //setChainToPosition(InMol,'B',InParams->bx/2.0+40,InParams->by/2.0,InParams->bz/2.0);
    setChainToPosition(InMol,'C',InParams->bx/2.0+60,InParams->by/2.0,InParams->bz/2.0);
    minimize(InOpenMM->context);

    simulateOpenMM(InMol,InParams,InOpenMM);
    return 0;
}
int doTemperatureRepd(Molecule *InMol,OpenMMParameters* InParams, 
        OpenMMData* InOpenMM){
    InParams->usesPeriodicBoundaryConditions=true;
    InParams->bx=150;
    InParams->by=150;
    InParams->bz=150;
    //InParams->outDCDName="1.dcd";
    InParams->steps=3000000000;
    InParams->printFreq=10000;
    InParams->temperature=300;

    setChainToPosition(InMol,'A',InParams->bx/2.0,InParams->by/2.0,InParams->bz/2.0);
    setChainToPosition(InMol,'B',InParams->bx/2.0-50,InParams->by/2.0,InParams->bz/2.0);
    setChainToPosition(InMol,'C',InParams->bx/2.0+50,InParams->by/2.0,InParams->bz/2.0);
    //minimize(InOpenMM->context);

    RepdParameters repdParams;
    repdParams.exchangeFreq=-1;
    //repdParams.temperatures
    //for(int ii=0; ii<5; ++ii){
        //repdParams.temperatures.push_back(300+ii*10);
        //repdParams.temperatures.push_back(300);
    //}
    simulateOpenMMRepdTemperature(InMol,InParams,InOpenMM,&repdParams);
    return 0;
}
int doHamiltonianRepd(Molecule *InMol,OpenMMParameters* InParams, 
        OpenMMData* InOpenMM){
    InParams->steps=400000000;
    InParams->printFreq=10000;
    InParams->temperature=300;

    //KIX-MLL
    //setChainToPosition(InMol,'A',InParams->bx/2.0,InParams->by/2.0,InParams->bz/2.0);
    //setChainToPosition(InMol,'B',InParams->bx/2.0-40,InParams->by/2.0,InParams->bz/2.0);
    //KIX-Myb
    //setChainToPosition(InMol,'A',InParams->bx/2.0,InParams->by/2.0,InParams->bz/2.0);
    //setChainToPosition(InMol,'B',InParams->bx/2.0,InParams->by/2.0,InParams->bz/2.0);
    //setChainToPosition(InMol,'C',InParams->bx/2.0,InParams->by/2.0,InParams->bz/2.0);
    //setMultChainToPosition(InMol,{'B','A','C'}, InParams->bx/2.0, InParams->by/2.0, InParams->bz/2.0);
    //setMultChainToPosition(InMol,{'A','B','C'}, InParams->bx/2.0, InParams->by/2.0, InParams->bz/2.0);
    setMultChainToPosition(InMol,{'A','B'}, InParams->bx/2.0, InParams->by/2.0, InParams->bz/2.0);
    setMultChainToPosition(InMol,{'A','B','C'}, InParams->bx/2.0, InParams->by/2.0, InParams->bz/2.0);
    //setMultChainToPosition(InMol,{'A'}, InParams->bx/2.0, InParams->by/2.0, InParams->bz/2.0);

    RepdParameters repdParams;
    repdParams.exchangeFreq=10000;
    //repdParams.exchangeFreq=-1;
    repdParams.dim1=8;
    repdParams.dim2=3; //20k per replica
    simulateOpenMMRepdHamiltonian(InMol,InParams,InOpenMM,&repdParams);
    return 0;
}

int doHamiltonianRepd3D(Molecule *InMol,OpenMMParameters* InParams, 
        OpenMMData* InOpenMM){
    InParams->steps=400000000;
    InParams->printFreq=10000;
    InParams->temperature=300;

    //setChainToPosition(InMol,'A',InParams->bx/2.0,InParams->by/2.0,InParams->bz/2.0);
    //setChainToPosition(InMol,'B',InParams->bx/2.0-20,InParams->by/2.0,InParams->bz/2.0);

    RepdParameters repdParams;
    repdParams.exchangeFreq=10000;
    repdParams.dim1=16;
    repdParams.dim2=4; //15k per replica
    repdParams.dim3=4; //15k per replica
    simulateOpenMMRepdHamiltonian3D(InMol,InParams,InOpenMM,&repdParams);
    return 0;
}
