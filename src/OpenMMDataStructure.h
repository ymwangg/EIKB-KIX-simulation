#ifndef __OPENMM_DATA_STRUCTURE_H__
#define __OPENMM_DATA_STRUCTURE_H__
#include "OpenMM.h"

struct OpenMMData;
struct OpenMMParameters;
struct OpenMMData {
    //OpenMMData() : system(0), context(0), integrator(0) {}
    //~OpenMMData() {delete context; delete integrator; delete system;}
    OpenMM::System*         system;
    OpenMM::Integrator*     integrator;
    OpenMM::Context*  context;
    OpenMM::Platform* platform;
    std::map<const std::string,int> ForceIdInSystem;
};
struct OpenMMParameters {
    OpenMMParameters(): 
        usesBondConstraint(true),usesPeriodicBoundaryConditions(true),
        bx(150),by(150),bz(150),temperature(300),timestep(0.022),fbeta(0.1),constraintTol(1e-6){}
    //~OpenMMParameters(){};
    //system
    bool usesBondConstraint;
    //periodic boundary conditions
    bool usesPeriodicBoundaryConditions;
    float bx,by,bz; //angstrom
    float temperature; //kelvin
    float timestep; //picoseconds
    float fbeta; //picoseconds^(-1)
    float constraintTol;
    std::string outDCDName;
    unsigned long int steps;
    unsigned long int printFreq;
    //repd
    void *logFileHandle;
    void *repdFileHandle;
    unsigned long int currentTime;
    int idx1;
    int idx2;
    int idx3;

};
struct RepdParameters {
    RepdParameters(): dim2(1){};
    unsigned long int exchangeFreq;
    std::vector<float> temperatures;
    int dim1;
    int dim2;
    int dim3;
};

#endif
