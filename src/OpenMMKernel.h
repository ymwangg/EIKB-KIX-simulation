#ifndef __OPENMM_KERNEL_H__
#define __OPENMM_KERNEL_H__
#include "OpenMM.h"
#include "io/dcdplugin.cpp"
#include "DataStructure.h"
#include "OpenMMDataStructure.h"
#include "DataType.h"
#include "Coordinate.h"
#include <mpi.h>
using OpenMM::Vec3;

#define BOND_CUTOFF 3
#define NONBOND_CUTOFF 25
#define USE_ELEC true
#define USE_BACKBONE_CORRECTION false
#define USE_CENTROID_RESTRAINT true
#define USE_MPI_CUDA
//#define testOpenMMRepd
//#define testOMMNonbondedExclusion
#define testOMMZincBond


OpenMM::System* createSystem(const Molecule *InMol, 
        const OpenMMParameters* InParams,
        OpenMMData* InOpenMM,std::vector<Vec3> &initialPosInNm);

/*forces*/
OpenMM::HarmonicBondForce* addBondForce(const Molecule *InMol, 
        const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int forceId);

OpenMM::HarmonicAngleForce* addAngleForce(const Molecule *InMol, 
        const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int forceId);

OpenMM::PeriodicTorsionForce* addDihedralForce(const Molecule *InMol, 
        const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int forceId);

OpenMM::CustomTorsionForce* addDihedralCorrectionForce(const Molecule *InMol, 
        const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int forceId);

OpenMM::CustomNonbondedForce* addNonbondedForce(const Molecule *InMol, 
        const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int forceId);

OpenMM::CustomBondForce* addNativeForce(const Molecule *InMol,
        const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int forceId, 
        OpenMM::CustomNonbondedForce *nonBondedForce);

int addNonbondedExclusionAndCutoff(const Molecule *InMol,
        const OpenMMParameters* InParams,
        OpenMM::CustomNonbondedForce *nonBondedForce,
        OpenMM::CustomBondForce* nativeForce, const int nonBondedExclusionCutoff, 
        const int distCutoffInAnstrom);

OpenMM::CustomNonbondedForce* addElecForce(const Molecule *InMol, 
        const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int forceId);

int addElecExclusionAndCutoff(const Molecule *InMol,const OpenMMParameters* InParams,
        OpenMM::CustomNonbondedForce *elecForce, 
        OpenMM::CustomBondForce* nativeForce, const int elecExclusionCutoff, 
        const int distCutoffInAngstrom);

int addCentroidForce(const Molecule *InMol, OpenMMData* InOpenMM, const char chainId,
const double k, const double x0, const double y0, const double z0);

int addInterChainCentroidBondForce(const Molecule *InMol, OpenMMData* InOpenMM,const char chainA, const char chainB,const double k, const double r0);

int addZincBondForce(const Molecule *InMol, 
        const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int forceId);

OpenMM::CustomBondForce* addElecBondForce(const Molecule *InMol,
        const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int elecExclusionCutoff, 
        const int forceId);

int setupOpenMMSystem(const Molecule *InMol,const OpenMMParameters* InParams,
        OpenMMData* InOpenMM);
int Molecule2Molfile(const Molecule *InMol,const OpenMMParameters* InParams,
        molfile_timestep_t *ts);
int simulateOpenMM(const Molecule* InMol, const OpenMMParameters* InParams,
        OpenMMData* InOpenMM);

int resetPositionsVelocities(const Molecule *InMol,const OpenMMParameters* InParams,
        OpenMMData* InOpenMM);
int resetPositions(const Molecule *InMol, OpenMMData* InOpenMM);
int resetVelocities(const OpenMMParameters* InParams, OpenMMData* InOpenMM);
int dumpOpenMMEnergy(const OpenMMData* InOpenMM);
int minimize(OpenMM::Context *context,const double tolerance,int maxIterations);
int resetNativeForce(const Molecule *InMol, OpenMMData* InOpenMM);


int setupOpenMMSystem(const Molecule *InMol,const OpenMMParameters* InParams,
        OpenMMData* InOpenMM){
    //load OpenMM plugins
    //std::string pluginsDirectory="/home/ymwang/apps/openmm_src/openmm_bin/lib/plugins";
    //std::string pluginsDirectory="/home/ymwang/apps/anaconda2/pkgs/openmm-7.1.1-py27_0/lib/plugins/";
    std::string pluginsDirectory="/home/ymwang/apps/anaconda2/pkgs/openmm-7.2.1-py27_1/lib/plugins/";
    OpenMM::Platform::loadPluginsFromDirectory(pluginsDirectory);
    //OpenMM::Platform::loadPluginsFromDirectory
    //    (OpenMM::Platform::getDefaultPluginsDirectory());
    int numOfPlatforms=OpenMM::Platform::getNumPlatforms();
    printf("Number of registered platforms=%d\n",numOfPlatforms);
    OpenMM::Platform* platform;
    float maxspeed=0;
    for(int ii=0; ii<numOfPlatforms;++ii){
        OpenMM::Platform *plt= &OpenMM::Platform::getPlatform(ii);
        std::string name=plt->getName();
        float speed=plt->getSpeed();
        printf("Platform%d = %s speed=%f\n",ii,name.c_str(),speed);
        if(speed > maxspeed){
            maxspeed = speed;
            platform=plt;
        }
    }
    printf("Platform %s was chosen as the fastest\n",platform->getName().c_str());
#ifdef USE_MPI_CUDA
    const char* deviceId = std::getenv("CUDA_VISIBLE_DEVICES");
    char newEnv[128] = "CUDA_VISIBLE_DEVICES=0,1,2,3";
    putenv(newEnv);
    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD,&pid);
    char buff[2];
    sprintf(buff,"%d",pid%2);
    printf("pid%d deviceid='%s'\n",pid,deviceId);
    platform->setPropertyDefaultValue("DeviceIndex",deviceId);
#endif
    InOpenMM->platform = platform;

    /*create system*/
    std::vector<Vec3> initialPosInNm;
    OpenMM::System *system;
    system = createSystem(InMol,InParams,InOpenMM,initialPosInNm);
    InOpenMM->system = system;

    /*create forces*/
    /*bond*/
    OpenMM::HarmonicBondForce*  bondForce;
    bondForce=addBondForce(InMol,InParams,InOpenMM,1);

    /*angle*/
    OpenMM::HarmonicAngleForce *angleForce;
    angleForce=addAngleForce(InMol,InParams,InOpenMM,2);

    /*dihedral*/
    OpenMM::PeriodicTorsionForce* dihedralForce;
    dihedralForce=addDihedralForce(InMol,InParams,InOpenMM,3);

    /*non-native nonbond*/
    OpenMM::CustomNonbondedForce *nonBondedForce;
    nonBondedForce=addNonbondedForce(InMol,InParams,InOpenMM,4);

    /*native nonbond*/
    OpenMM::CustomBondForce* nativeForce;
    nativeForce=addNativeForce(InMol,InParams,InOpenMM,5,nonBondedForce);

    /*electrostatics using nonbond*/
    //OpenMM::CustomNonbondedForce *elecForce;
    //elecForce=addElecForce(InMol,InParams,InOpenMM,6);

    /*electrostatics using bond*/
#if USE_ELEC
    OpenMM::CustomBondForce* elecBondForce;
    elecBondForce=addElecBondForce(InMol,InParams,InOpenMM,BOND_CUTOFF,7);
    printf("Using electrostatics\n");
#endif

    /*add nonbonded exclusions*/
    /*nonbond residue cutoff=3, distance cutoff=25.0 Anstrom*/
    addNonbondedExclusionAndCutoff(InMol,InParams,nonBondedForce,nativeForce,
            BOND_CUTOFF, NONBOND_CUTOFF);

    /*add elec exclusions*/
    //addElecExclusionAndCutoff(InMol,InParams,elecForce,nativeForce, 3,24.0);

    /*backbone correction*/
#if USE_BACKBONE_CORRECTION
    OpenMM::CustomTorsionForce* dihedralCorrectionForce;
    dihedralCorrectionForce=addDihedralCorrectionForce(InMol,InParams,InOpenMM,8);
    printf("Using backbone correction\n");
#endif

    /*zincBond*/
    addZincBondForce(InMol,InParams,InOpenMM,9);

#if USE_CENTROID_RESTRAINT
    //TAZ1
    //addCentroidForce(InMol, InOpenMM, 'A', 50, InParams->bx/2.0, InParams->by/2.0, InParams->bz/2.0);
    //KIX
    addCentroidForce(InMol, InOpenMM, 'A', 50, InParams->bx/2.0, InParams->by/2.0, InParams->bz/2.0);
#endif

    /*create integrator*/
    InOpenMM->integrator = new OpenMM::LangevinIntegrator(InParams->temperature,
            InParams->fbeta,InParams->timestep);
    /*create context*/
#ifdef USE_MPI_CUDA
    InOpenMM->context = new OpenMM::Context(*InOpenMM->system,*InOpenMM->integrator,*platform);
#else
    InOpenMM->context = new OpenMM::Context(*InOpenMM->system,*InOpenMM->integrator);
#endif
    /*initialize positions*/
    InOpenMM->context->setPositions(initialPosInNm);
    const std::map<std::string, double> params=InOpenMM->context->getParameters();
    for(auto const& x : params){
        printf("%s=%f\n",x.first.c_str(),x.second);
    }
    /*check system*/
    printf("OpenMMKernel: system NumParticles=%d\n",system->getNumParticles());
    printf("OpenMMKernel: system NumOfForces=%d\n",system->getNumForces());
    printf("OpenMMKernel: system NumOfConstraints=%d\n",system->getNumConstraints());
    //printf("OpenMMKernel: nonBondedForce NumExclusions=%d\n",nonBondedForce->getNumExclusions());
    printf("OpenMMKernel: numOfNatives=%d\n",InMol->NumOfNatives);
    dumpOpenMMEnergy(InOpenMM);
    printf("OpenMMKernel: finish context setup\n");

    printf("*******************************************\n");
    printf("Platform properties\n");
    OpenMM::Platform *platformInContext = &InOpenMM->context->getPlatform();
    std::vector<std::string> property=platformInContext->getPropertyNames();
    for(std::string s : property){
        printf("%s\n",s.c_str());
        std::string propertyValue=platformInContext->getPropertyValue(
                *InOpenMM->context,s);
        printf("value=%s\n",propertyValue.c_str());
    }
    return 0;
}

int dumpOpenMMEnergy(const OpenMMData* InOpenMM){
    float nonbond=0;
    OpenMM::State state = InOpenMM->context->getState(OpenMM::State::Energy,false,1<<1);
    printf("OpenMMKernel: bond=%f kcal/mol\n",state.getPotentialEnergy()*OpenMM::KcalPerKJ);
    state = InOpenMM->context->getState(OpenMM::State::Energy,false,1<<2);
    printf("OpenMMKernel: angle=%f kcal/mol\n",state.getPotentialEnergy()*OpenMM::KcalPerKJ);
    state = InOpenMM->context->getState(OpenMM::State::Energy,false,1<<3);
    printf("OpenMMKernel: dihedral=%f kcal/mol\n",state.getPotentialEnergy()*OpenMM::KcalPerKJ);
    state = InOpenMM->context->getState(OpenMM::State::Energy,false,1<<4);
    nonbond+=state.getPotentialEnergy()*OpenMM::KcalPerKJ;
    printf("OpenMMKernel: nonbond=%f kcal/mol\n",state.getPotentialEnergy()*OpenMM::KcalPerKJ);
    state = InOpenMM->context->getState(OpenMM::State::Energy,false,1<<5);
    nonbond+=state.getPotentialEnergy()*OpenMM::KcalPerKJ;
    printf("OpenMMKernel: native=%f kcal/mol\n",state.getPotentialEnergy()*OpenMM::KcalPerKJ);
    printf("OpenMMKernel: nonbond=%f kcal/mol\n",nonbond);
    state = InOpenMM->context->getState(OpenMM::State::Energy,false,1<<6);
    printf("OpenMMKernel: elec=%f kcal/mol\n",state.getPotentialEnergy()*OpenMM::KcalPerKJ);
    state = InOpenMM->context->getState(OpenMM::State::Energy,false,1<<7);
    printf("OpenMMKernel: elecBond=%f kcal/mol\n",state.getPotentialEnergy()*OpenMM::KcalPerKJ);
    state = InOpenMM->context->getState(OpenMM::State::Energy,false,1<<8);
    printf("OpenMMKernel: dihedralCorrection=%f kcal/mol\n",state.getPotentialEnergy()*OpenMM::KcalPerKJ);
    state = InOpenMM->context->getState(OpenMM::State::Energy,false,1<<9);
    printf("OpenMMKernel: zincBond=%f kcal/mol\n",state.getPotentialEnergy()*OpenMM::KcalPerKJ);
    state = InOpenMM->context->getState(OpenMM::State::Energy,false,1<<10);
    printf("OpenMMKernel: resCentroid=%f kcal/mol\n",state.getPotentialEnergy()*OpenMM::KcalPerKJ);
    state = InOpenMM->context->getState(OpenMM::State::Energy,false,1<<11);
    printf("OpenMMKernel: centroidBond=%f kcal/mol\n",state.getPotentialEnergy()*OpenMM::KcalPerKJ);
    return 0;
}

//run a simulation
int simulateOpenMM(Molecule* InMol, const OpenMMParameters* InParams,
        OpenMMData* InOpenMM){
    /*initialize positions*/
    std::vector<Vec3> initialPosInNm;
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        const Vec3 posInNm(InMol->X[ii]*OpenMM::NmPerAngstrom,
                InMol->Y[ii]*OpenMM::NmPerAngstrom,
                InMol->Z[ii]*OpenMM::NmPerAngstrom);
        initialPosInNm.push_back(posInNm);
    }
    InOpenMM->context->setPositions(initialPosInNm);

    /*initialize velocities*/
    InOpenMM->context->setVelocitiesToTemperature(InParams->temperature);

    /*apply constraint*/
    InOpenMM->context->applyConstraints(InParams->constraintTol);
    InOpenMM->context->applyVelocityConstraints(InParams->constraintTol);

    /*dcd file handle*/
    /*ts is needed as an intermediate file between Molecule and dcdfile*/
    void *v = open_dcd_write(InParams->outDCDName.c_str(), "dcd", InMol->NumOfAtoms);
    dcdhandle *dcd;
    molfile_timestep_t ts;
    ts.coords = (float *)malloc(3*sizeof(float)*InMol->NumOfAtoms);
    dcd = (dcdhandle *)v;
    /*set what information should be retrived from each timestep*/
    int infomask = OpenMM::State::Energy + OpenMM::State::Positions;

    /*set integrator temperature*/
    dynamic_cast<OpenMM::LangevinIntegrator*>(InOpenMM->integrator)->setTemperature(InParams->temperature);
    float temp = dynamic_cast<OpenMM::LangevinIntegrator*>(InOpenMM->integrator)->getTemperature();
    //printf("temperature=%f\n",temp);
    OpenMM::State state = InOpenMM->context->getState(infomask,true);
    //calculate the number of steps to save as dcd file
    int steps=ceil(float(InParams->steps)/float(InParams->printFreq));
    //print energy decomposition
    dumpOpenMMEnergy(InOpenMM);

    for(int ii=0; ii<steps; ++ii){
        InOpenMM->integrator->step(InParams->printFreq);
        state=InOpenMM->context->getState(infomask,true);
        printf("step=%d V=%f EK=%f\n",ii,state.getPotentialEnergy()*OpenMM::KcalPerKJ, state.getKineticEnergy()*OpenMM::KcalPerKJ);
        std::vector<Vec3> positions = state.getPositions();
        for(int n=0; n<InMol->NumOfAtoms; ++n){
            //printf("atom%d ",n);
            //printf("%f %f %f\n",positions[n][0]*OpenMM::AngstromsPerNm,
            //        positions[n][1]*OpenMM::AngstromsPerNm,
            //        positions[n][2]*OpenMM::AngstromsPerNm);
            InMol->X[n]=positions[n][0]*OpenMM::AngstromsPerNm;
            InMol->Y[n]=positions[n][1]*OpenMM::AngstromsPerNm;
            InMol->Z[n]=positions[n][2]*OpenMM::AngstromsPerNm;
        }
        //move back to original box
        moveBackToOriginalBox(InMol, InParams->bx, InParams->by, InParams->bz);
        resetPositions(InMol,InOpenMM);
        Molecule2Molfile(InMol,InParams,&ts);
        write_timestep(v,&ts);
    }
    close_file_write(v);
    printf("Done simulation\n");
    return 0;
}

int resetPositionsVelocities(const Molecule *InMol,const OpenMMParameters* InParams,
        OpenMMData* InOpenMM){
    std::vector<Vec3> newPosInNm;
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        const Vec3 posInNm(InMol->X[ii]*OpenMM::NmPerAngstrom,
                InMol->Y[ii]*OpenMM::NmPerAngstrom,
                InMol->Z[ii]*OpenMM::NmPerAngstrom);
        newPosInNm.push_back(posInNm);
    }
    InOpenMM->context->setPositions(newPosInNm);
    /*initialize velocities*/
    InOpenMM->context->setVelocitiesToTemperature(InParams->temperature);
    /*apply constraint*/
    InOpenMM->context->applyConstraints(InParams->constraintTol);
    InOpenMM->context->applyVelocityConstraints(InParams->constraintTol);
    dynamic_cast<OpenMM::LangevinIntegrator*>(InOpenMM->integrator)->
        setTemperature(InParams->temperature);
    return 0;
}

int resetPositions(const Molecule *InMol, OpenMMData* InOpenMM){
    std::vector<Vec3> newPosInNm;
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        const Vec3 posInNm(InMol->X[ii]*OpenMM::NmPerAngstrom,
                InMol->Y[ii]*OpenMM::NmPerAngstrom,
                InMol->Z[ii]*OpenMM::NmPerAngstrom);
        newPosInNm.push_back(posInNm);
    }
    InOpenMM->context->setPositions(newPosInNm);
    return 0;
}

int resetVelocities(const OpenMMParameters* InParams,
        OpenMMData* InOpenMM){
    /*initialize velocities*/
    InOpenMM->context->setVelocitiesToTemperature(InParams->temperature);
    /*apply constraint*/
    InOpenMM->context->applyConstraints(InParams->constraintTol);
    InOpenMM->context->applyVelocityConstraints(InParams->constraintTol);
    dynamic_cast<OpenMM::LangevinIntegrator*>(InOpenMM->integrator)->
        setTemperature(InParams->temperature);
    return 0;
}

OpenMM::System* createSystem(const Molecule *InMol,const OpenMMParameters* InParams,OpenMMData* InOpenMM,std::vector<Vec3> &initialPosInNm){
    OpenMM::System* system = new OpenMM::System();
    //add particles and initialize positions
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        system->addParticle(InMol->Mass[ii]);
        const Vec3 posInNm(InMol->X[ii]*OpenMM::NmPerAngstrom,
                InMol->Y[ii]*OpenMM::NmPerAngstrom,
                InMol->Z[ii]*OpenMM::NmPerAngstrom);
        initialPosInNm.push_back(posInNm);
    }
    //set up periodic boundary condition
    if(InParams->usesPeriodicBoundaryConditions==true && InParams->bx>0 &&InParams->by>0
            &&InParams->bz>0){
        system->setDefaultPeriodicBoxVectors(Vec3(InParams->bx*OpenMM::NmPerAngstrom,0,0),
                Vec3(0,InParams->by*OpenMM::NmPerAngstrom,0),Vec3(0,0,InParams->bz*OpenMM::NmPerAngstrom));
        printf("OpenMMKernel: uses periodic boundary condition dx=%f dy=%f dz=%f\n",
                InParams->bx,InParams->by,InParams->bz);
    } else {
        printf("OpenMMKernel: periodic boundary condition is not used\n");
    }
    return system;
}
//add zinc bond force
int addZincBondForce(const Molecule *InMol, 
        const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int forceId){
    OpenMM::CustomBondForce* zincBondForce = new OpenMM::CustomBondForce("k*(r-r0)^2");
    zincBondForce->addPerBondParameter("r0");
    zincBondForce->addPerBondParameter("k");
    for(int ii=0;ii<InMol->NumOfZincBonds;++ii){
        int atom1=InMol->ZINCBOND[ii].i;
        int atom2=InMol->ZINCBOND[ii].j;
        float k=InMol->ZINCBOND[ii].k;
        float r=InMol->ZINCBOND[ii].r;

        //if k=0, no need to add this interaction
        std::vector<double> Params={r*OpenMM::NmPerAngstrom,fabs(k*OpenMM::KJPerKcal)};
        int bondId = zincBondForce->addBond(atom1,atom2,Params);
#ifdef testOMMZincBond
        printf("OpenMMKernel: add zincBond (%d,%d) k=%f r=%f\n",atom1,atom2,k,r);
#endif
    }
    zincBondForce->setForceGroup(forceId);
    int IdInSystem=InOpenMM->system->addForce(zincBondForce);
    printf("zincBondForce IdInSystem = %d\n",IdInSystem);
    InOpenMM->ForceIdInSystem.insert(std::pair<const std::string,int>("zincBondForce",IdInSystem));
    return 0;
};

//add bond force
OpenMM::HarmonicBondForce* addBondForce(const Molecule *InMol, const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int forceId){
    OpenMM::HarmonicBondForce*  bondForce;
    //if bond constraint is applied, no need to add bond forces
    if(InParams->usesBondConstraint==false){
        OpenMM::HarmonicBondForce*  bondForce = new OpenMM::HarmonicBondForce();
    }
    //add bond pair to bondForce
    std::vector<std::pair<int,int>> bondPairs;
    for(int ii=0; ii<InMol->NumOfBonds;++ii){
        int atom1=InMol->BOND[ii].i;
        int atom2=InMol->BOND[ii].j;
        float k=InMol->BOND[ii].k;
        float r=InMol->BOND[ii].r;
        bondPairs.push_back(std::make_pair(atom1,atom2));
        if(InParams->usesBondConstraint==false){
            bondForce->addBond(atom1,atom2,r*OpenMM::NmPerAngstrom,2*k*
                    OpenMM::KJPerKcal*OpenMM::AngstromsPerNm*OpenMM::AngstromsPerNm);
#ifdef testOMMBond
            printf("OpenMMKernel: add bond (%d,%d) k=%f r=%f\n",atom1,atom2,k,r);
#endif
        } else {
            InOpenMM->system->addConstraint(atom1,atom2,r*OpenMM::NmPerAngstrom);
#ifdef testOMMBond
            printf("OpenMMKernel: add bond (%d,%d) k=%f r=%f\n",atom1,atom2,k,r);
#endif
        }
    }
    if(InParams->usesBondConstraint==false){
        bondForce->setForceGroup(forceId);
        InOpenMM->system->addForce(bondForce);
    }
    return bondForce;
}
OpenMM::HarmonicAngleForce* addAngleForce(const Molecule *InMol, const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int forceId){
    OpenMM::HarmonicAngleForce* angleForce = new OpenMM::HarmonicAngleForce();
    for(int ii=0; ii<InMol->NumOfAngles;++ii){
        int atom1=InMol->ANGLE[ii].a;
        int atom2=InMol->ANGLE[ii].b;
        int atom3=InMol->ANGLE[ii].c;
        float k=InMol->ANGLE[ii].k;
        float theta=InMol->ANGLE[ii].theta;
        angleForce->addAngle(atom1,atom2,atom3,theta*OpenMM::RadiansPerDegree,
                2*k*OpenMM::KJPerKcal);
#ifdef testOMMAngle
        printf("OpenMMKernel: add angle (%d,%d,%d) k=%f theta=%f\n",atom1,atom2,atom3,k,theta);
#endif
    }
    angleForce->setForceGroup(forceId);
    InOpenMM->system->addForce(angleForce);
    return angleForce;
}

OpenMM::PeriodicTorsionForce* addDihedralForce(const Molecule *InMol, const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int forceId){
    OpenMM::PeriodicTorsionForce* dihedralForce = new OpenMM::PeriodicTorsionForce();
    for(int ii=0;ii<InMol->NumOfDihedrals;++ii){
        int atom1=InMol->DIHEDRAL[ii].a;
        int atom2=InMol->DIHEDRAL[ii].b;
        int atom3=InMol->DIHEDRAL[ii].c;
        int atom4=InMol->DIHEDRAL[ii].d;
        float k=InMol->DIHEDRAL[ii].k;
        int mult=InMol->DIHEDRAL[ii].multiplicity;
        float theta=InMol->DIHEDRAL[ii].theta;
        dihedralForce->addTorsion(atom1,atom2,atom3,atom4,mult,theta*OpenMM::RadiansPerDegree,
                k*OpenMM::KJPerKcal);
#ifdef testOMMDihedral
        printf("OpenMMKernel: add dihedral (%d,%d,%d,%d) k=%f mult=%d theta=%f\n",
               atom1,atom2,atom3,atom4,k,mult,theta);
#endif
    }
    dihedralForce->setForceGroup(forceId);
    InOpenMM->system->addForce(dihedralForce);
    return dihedralForce;
}

//backbone dihedral correction
OpenMM::CustomTorsionForce* addDihedralCorrectionForce(const Molecule *InMol, const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int forceId){
    OpenMM::CustomTorsionForce* dihedralForce = new OpenMM::CustomTorsionForce("kDC*cos(theta-theta0)");
    float k=-1.16; //kcal/mol
    float theta0=297.35; //degree
    dihedralForce->addGlobalParameter("kDC",k*OpenMM::KJPerKcal);
    dihedralForce->addGlobalParameter("theta0",theta0*OpenMM::RadiansPerDegree);
    for(int ii=0;ii<InMol->NumOfDihedrals;++ii){
        int atom1=InMol->DIHEDRAL[ii].a;
        int atom2=InMol->DIHEDRAL[ii].b;
        int atom3=InMol->DIHEDRAL[ii].c;
        int atom4=InMol->DIHEDRAL[ii].d;
        dihedralForce->addTorsion(atom1,atom2,atom3,atom4);
#ifdef testOMMDihedralCorrection
        printf("OpenMMKernel: add dihedral (%d,%d,%d,%d) k=%f theta=%f\n",
               atom1,atom2,atom3,atom4,k,theta0);
#endif
    }
    dihedralForce->setForceGroup(forceId);
    InOpenMM->system->addForce(dihedralForce);
    return dihedralForce;
}
OpenMM::CustomNonbondedForce* addNonbondedForce(const Molecule *InMol, const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int forceId){
    OpenMM::CustomNonbondedForce* nonBondedForce = new OpenMM::CustomNonbondedForce(
            "sqrt(epsilon1 * epsilon2) * (13*sr^12 - 18*sr^10 + 4*sr^6); sr = (rmin1 + rmin2)/(2*r)");
            //"sqrt(epsilon1 * epsilon2) * (sr^12); sr = (rmin1 + rmin2)/(2*r)");
    nonBondedForce->addPerParticleParameter("rmin");
    nonBondedForce->addPerParticleParameter("epsilon");
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        //actually the VDWR(ii) is rmin(ii)/2
        std::vector<double> nonbondedParams={2*InMol->VDWR[ii]*OpenMM::NmPerAngstrom,
            InMol->EPS[ii]*OpenMM::KJPerKcal};
        nonBondedForce->addParticle(nonbondedParams);
    }
    nonBondedForce->setForceGroup(forceId);
    InOpenMM->system->addForce(nonBondedForce);
    return nonBondedForce;
}
OpenMM::CustomBondForce* addNativeForce(const Molecule *InMol,const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int forceId, OpenMM::CustomNonbondedForce *nonBondedForce){
    /*
    OpenMM::CustomBondForce* nativeForce = new OpenMM::CustomBondForce(
            "step(cutoff^2-pr2) * epsilon * (13*sr^6 - 18*sr^5 + 4*sr^3);\\
            sr = rmin^2 / pr2; pr2 = dx^2 + dy^2 + dz^2;\\
            dx = dx - bx*ceil(dx/bx-0.5); dx = x2 - x1;\\
            dy = dy - by*ceil(dy/by-0.5); dy = y2 - y1;\\
            dz = dz - bz*ceil(dz/bz-0.5); dz = z2 - z1");
            */
    /*non-native nonbond*/
    OpenMM::CustomBondForce* nativeForce = new OpenMM::CustomBondForce(
            //"step(cutoff - r) * epsilon * (13 * sr^12 - 18 * sr^10 + 4 * sr^6);sr = rmin / r");
            "epsilon * (13 * sr^12 - 18 * sr^10 + 4 * sr^6);sr = rmin / r");
    nativeForce->addPerBondParameter("rmin");
    nativeForce->addPerBondParameter("epsilon");
    //nativeForce->addGlobalParameter("cutoff",99.0);
    //add exclusions
    //nonBondedForce->createExclusionsFromBonds(bondPairs,nonBondedExclusionCutoff);
    for(int ii=0;ii<InMol->NumOfNatives;++ii){
        int atom1=InMol->NATIVE[ii].i;
        int atom2=InMol->NATIVE[ii].j;
#if USE_ELEC
        //if uses elec force
        //the force constants of the salt bridge interactions is decreased by 40%
        float k=InMol->NATIVE[ii].k;
        if(InMol->Charge[atom1]*InMol->Charge[atom2]<0){
            k=k*0.6;
            printf("Salt bridge interaction(%d(%f),%d(%f)) deacrease from %f to %f\n",
                    atom1,InMol->Charge[atom1],atom2,InMol->Charge[atom2],
                    InMol->NATIVE[ii].k,k);
        }
#else
        float k=InMol->NATIVE[ii].k;
#endif

        float r=InMol->NATIVE[ii].r;
        //if k=0, no need to add this interaction
        std::vector<double> nativeParams={r*OpenMM::NmPerAngstrom,fabs(k*OpenMM::KJPerKcal)};
        int bondId = nativeForce->addBond(atom1,atom2,nativeParams);
#ifdef testOMMNative
        printf("OpenMMKernel: add native (%d,%d) k=%f r=%f\n",atom1,atom2,k,r);
        printf("OpenMMKernel: add native exclusion (%d,%d)\n",atom1,atom2);
#endif
        nonBondedForce->addExclusion(atom1,atom2);
    }
    nativeForce->setForceGroup(forceId);
    int IdInSystem=InOpenMM->system->addForce(nativeForce);
    printf("nativeForce IdInSystem = %d\n",IdInSystem);
    InOpenMM->ForceIdInSystem.insert(std::pair<const std::string,int>("nativeForce",IdInSystem));
    return nativeForce;
}
int resetNativeForce(const Molecule *InMol, OpenMMData* InOpenMM){
    int IdInSystem=InOpenMM->ForceIdInSystem["nativeForce"];
    OpenMM::CustomBondForce* nativeForce = dynamic_cast<OpenMM::CustomBondForce*>(&InOpenMM->system->getForce(IdInSystem));
    printf("Updating nativeForce\n");
    printf("Num of natives = %d\n",nativeForce->getNumBonds());
    for(int ii=0;ii<InMol->NumOfNatives;++ii){
        int atom1=InMol->NATIVE[ii].i;
        int atom2=InMol->NATIVE[ii].j;
        float k=InMol->NATIVE[ii].k;
        float r=InMol->NATIVE[ii].r;
        if(k==0){
            k= -0.0001;
        }
        std::vector<double> nativeParams={r*OpenMM::NmPerAngstrom,fabs(k*OpenMM::KJPerKcal)};
        nativeForce->setBondParameters(ii,atom1,atom2,nativeParams);
    }
    nativeForce->updateParametersInContext(*InOpenMM->context);
    return 0;
};

int addNonbondedExclusionAndCutoff(const Molecule *InMol,const OpenMMParameters* InParams,
        OpenMM::CustomNonbondedForce *nonBondedForce,
        OpenMM::CustomBondForce* nativeForce, const int nonBondedExclusionCutoff, 
        const int distCutoffInAngstrom)
{
    //zinc bond exclusion
    for(int ii=0;ii<InMol->NumOfZincBonds;++ii){
        int atom1=InMol->ZINCBOND[ii].i;
        int atom2=InMol->ZINCBOND[ii].j;
        if(InMol->ZINCBOND[ii].k>0){
            nonBondedForce->addExclusion(atom1,atom2);
#ifdef testOMMNonbondedExclusion
                    printf("OpenMMKernel: add zincBond exclusion (%d,%d)\n",atom1,atom2);
#endif
        }
    }
    //nonbond exclusion
    for(int ii=0;ii<InMol->NumOfAtoms;++ii){
        for(int jj=1;jj<nonBondedExclusionCutoff;++jj){
            bool flag=true;
            int atom1=ii;
            int atom2=ii+jj;
            if(atom2<InMol->NumOfAtoms && mystrcmp(InMol->AtomType[atom1],"ZN")!=0 &&
                    mystrcmp(InMol->AtomType[atom2],"ZN")!=0){
                for(int kk=0;kk<InMol->NumOfNatives;++kk){
                    if(atom1==InMol->NATIVE[kk].i&&atom2==InMol->NATIVE[kk].j ||
                            atom2==InMol->NATIVE[kk].i&&atom1==InMol->NATIVE[kk].j ){
                        //if(InMol->NATIVE[kk].k!=0){
                            flag=false;
                        //}
                    }
                }
                //make sure do not double exclude native contact
                //and do not exclude inter-chain residues
                if(flag && InMol->ChainId[atom1]==InMol->ChainId[atom2]){
                //if(flag){
                    nonBondedForce->addExclusion(atom1,atom2);
#ifdef testOMMNonbondedExclusion
                    printf("OpenMMKernel: add nonbond exclusion (%d,%d)\n",atom1,atom2);
#endif
                }
            }
        }
    }
    if(InParams->usesPeriodicBoundaryConditions==true){
        nativeForce->setUsesPeriodicBoundaryConditions(true);
        nonBondedForce->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
        nonBondedForce->setCutoffDistance(distCutoffInAngstrom*OpenMM::NmPerAngstrom);
    }else{
        nonBondedForce->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
        nonBondedForce->setCutoffDistance(distCutoffInAngstrom*OpenMM::NmPerAngstrom);
    }
    return 0;
}

OpenMM::CustomNonbondedForce* addElecForce(const Molecule *InMol, const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int forceId){
    const double DIELEC=138.93875744; //(kJ)*(mol-1)*(nm)*(unitcharge-2)
    OpenMM::CustomNonbondedForce* elecForce = new OpenMM::CustomNonbondedForce(
            "DIELEC*q1*q2*exp(-1.0*r/zeta)/(D*r)");
            //"DIELEC*q1*q2/(r)");
    elecForce->addPerParticleParameter("q");
    elecForce->addGlobalParameter("zeta",1.0); //nm
    elecForce->addGlobalParameter("D",80); //epsilon0
    elecForce->addGlobalParameter("DIELEC",DIELEC); //(kJ)*(mol-1)*(nm)*(unitcharge-2)
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        std::vector<double> nonbondedParams={InMol->Charge[ii]};
        elecForce->addParticle(nonbondedParams);
    }
    elecForce->setForceGroup(forceId);
    InOpenMM->system->addForce(elecForce);
    return elecForce;
}

int addElecExclusionAndCutoff(const Molecule *InMol,const OpenMMParameters* InParams,
        OpenMM::CustomNonbondedForce *elecForce, 
        OpenMM::CustomBondForce* nativeForce, const int elecExclusionCutoff, 
        const int distCutoffInAngstrom)
{
    for(int ii=0;ii<InMol->NumOfNatives;++ii){
        elecForce->addExclusion(InMol->NATIVE[ii].i,InMol->NATIVE[ii].j);
    }
    for(int ii=0;ii<InMol->NumOfAtoms;++ii){
        for(int jj=1;jj<elecExclusionCutoff;++jj){
            int atom1=ii;
            int atom2=ii+jj;
            bool foundInNative=false;
            for(int kk=0; kk<InMol->NumOfNatives; ++kk){
                int atomi,atomj;
                atomi=InMol->NATIVE[kk].i;
                atomj=InMol->NATIVE[kk].j;
                if((atom1==atomi&&atom2==atomj) || 
                        (atom1==atomj&&atom2==atomi)){
                    foundInNative=true;
                    break;
                }
            }
            if(atom2<InMol->NumOfAtoms&&!foundInNative){
                if(InMol->ChainId[atom1]==InMol->ChainId[atom2]){
                    elecForce->addExclusion(atom1,atom2);
                    printf("OpenMMKernel: add elec exclusion (%d,%d)\n",atom1,atom2);
                }
            }
        }
    }
    if(InParams->usesPeriodicBoundaryConditions==true){
        elecForce->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
        elecForce->setCutoffDistance(distCutoffInAngstrom*OpenMM::NmPerAngstrom);
    }else{
        elecForce->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
        elecForce->setCutoffDistance(distCutoffInAngstrom*OpenMM::NmPerAngstrom);
    }
    return 0;
}

OpenMM::CustomBondForce* addElecBondForce(const Molecule *InMol,const OpenMMParameters* InParams, OpenMMData* InOpenMM, const int elecExclusionCutoff, const int forceId){
    OpenMM::CustomBondForce* elecBondForce = new OpenMM::CustomBondForce(
            "DIELEC*q1*q2*exp(-1.0*r/zeta)/(D*r)");
            //"DIELEC*q1*q2/(D*r)");
            //"DIELEC*q1*q2/(r)");
    const double DIELEC=138.93875744; //(kJ)*(mol-1)*(nm)*(unitcharge-2)
    elecBondForce->addGlobalParameter("zeta",1.0); //nm
    elecBondForce->addGlobalParameter("D",80); //epsilon0
    elecBondForce->addGlobalParameter("DIELEC",DIELEC); //(kJ)*(mol-1)*(nm)*(unitcharge-2)
    //nonBondedForce->createExclusionsFromBonds(bondPairs,nonBondedExclusionCutoff);
    elecBondForce->addPerBondParameter("q1");
    elecBondForce->addPerBondParameter("q2");
    std::vector<int> chargedAtoms;
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        if(InMol->Charge[ii]!=0){
            chargedAtoms.push_back(ii);
        }
    }
    int nChargedAtoms=chargedAtoms.size();
    for(int ii=0; ii<nChargedAtoms; ++ii){
        int atom1=chargedAtoms[ii];
        float q1=InMol->Charge[atom1];
        for(int jj=ii+1; jj<nChargedAtoms; ++jj){
            int atom2=chargedAtoms[jj];
            if(abs(atom1-atom2)>=elecExclusionCutoff||
                    (InMol->ChainId[atom1]!=InMol->ChainId[atom2])||
                    mystrcmp(InMol->AtomType[atom1],"ZN")==0||
                    mystrcmp(InMol->AtomType[atom2],"ZN")==0){
                float q2=InMol->Charge[atom2];
                elecBondForce->addBond(atom1,atom2,{q1,q2});
                //                    printf("ElecBond (%d,%d) charge=(%f,%f)\n",atom1,atom2,
                //                            InMol->Charge[atom1],InMol->Charge[atom2]);
                for(int kk=0; kk<InMol->NumOfNatives; ++kk){
                    if((atom1==InMol->NATIVE[kk].i&&atom2==InMol->NATIVE[kk].j)||
                            (atom1==InMol->NATIVE[kk].j&&atom2==InMol->NATIVE[kk].i)){
                        //                            printf("ElecBond (%d,%d) charge=(%f,%f) is also native\n",atom1,atom2,
                        //                                    InMol->Charge[atom1],InMol->Charge[atom2]);
                    }
                }
            }
        }
    }
    //elecBondForce->addBond();
    if(InParams->usesPeriodicBoundaryConditions==true){
        elecBondForce->setUsesPeriodicBoundaryConditions(true);
    }
    elecBondForce->setForceGroup(forceId);
    InOpenMM->system->addForce(elecBondForce);
    return elecBondForce;
}

int addCentroidForce(const Molecule *InMol, OpenMMData* InOpenMM, const char chainId,
const double k, const double x0, const double y0, const double z0){
    /*add restraint*/
    //char buff[256];
    //sprintf(buff,"%f*( (x1-%f)^2 + (y1-%f)^2 + (z1-%f)^2 )",5.0*OpenMM::KJPerKcal,InParams->bx*OpenMM::NmPerAngstrom/2.0,InParams->by*OpenMM::NmPerAngstrom/2.0,InParams->bz*OpenMM::NmPerAngstrom/2.0);
    //sprintf(buff,"%f*( (x1-%f)^2 + (y1-%f)^2 + (z1-%f)^2 )",5.0,0,0,0);
    //OpenMM::CustomCentroidBondForce* resCentroidForce = new OpenMM::CustomCentroidBondForce(2,"5.0*((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)");
    OpenMM::CustomCentroidBondForce* resCentroidForce = new OpenMM::CustomCentroidBondForce(1,"k*((x1-x0)^2+(y1-y0)^2+(z1-z0)^2)");
    resCentroidForce->addPerBondParameter("k");
    resCentroidForce->addPerBondParameter("x0");
    resCentroidForce->addPerBondParameter("y0");
    resCentroidForce->addPerBondParameter("z0");

    std::vector<int> particles1;
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        if(InMol->ChainId[ii]==chainId){
            particles1.push_back(ii);
        }
    }
    if(particles1.size()==0){
        printf("Warning: addCentroidForce chain %c to (%f,%f,%f) k=%f failed\n",
                chainId,x0,y0,z0,k);
        return 1;
    }
    resCentroidForce->addGroup(particles1);

    std::vector<double> bondParameters;
    bondParameters.push_back(k*OpenMM::KJPerKcal*
            OpenMM::AngstromsPerNm*OpenMM::AngstromsPerNm);
    bondParameters.push_back(x0*OpenMM::NmPerAngstrom);
    bondParameters.push_back(y0*OpenMM::NmPerAngstrom);
    bondParameters.push_back(z0*OpenMM::NmPerAngstrom);

    std::vector<int> bondGroups;
    bondGroups.push_back(0);
    resCentroidForce->addBond(bondGroups,bondParameters);

    resCentroidForce->setForceGroup(10);
    InOpenMM->system->addForce(resCentroidForce);
    printf("addCentroidForce %s\n",resCentroidForce->getEnergyFunction().c_str());
    return 0;
}

int addInterChainCentroidBondForce(const Molecule *InMol, OpenMMData* InOpenMM,const char chainA, const char chainB,const double k, const double r0){
    OpenMM::CustomCentroidBondForce* force = new OpenMM::CustomCentroidBondForce(2,"k*(distance(g1,g2)-r0)^2");
    force->addPerBondParameter("k");
    force->addPerBondParameter("r0");
    std::vector<int> particles1;
    std::vector<int> particles2;
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        if(InMol->ChainId[ii]==chainA){
            particles1.push_back(ii);
        }
        else if(InMol->ChainId[ii]==chainB){
            particles2.push_back(ii);
        }
    }
    if(particles1.size()==0 || particles2.size()==0){
        printf("Warning: addInterChainCentroidBondForce (%c,%c) r=%f k=%f failed\n",
                chainA,chainB,r0,k);
        return 1;
    }
    force->addGroup(particles1);
    force->addGroup(particles2);

    std::vector<int> bondGroups;
    bondGroups.push_back(0);
    bondGroups.push_back(1);

    std::vector<double> bondParameters;
    bondParameters.push_back(k*OpenMM::KJPerKcal*OpenMM::AngstromsPerNm*
            OpenMM::AngstromsPerNm);
    bondParameters.push_back(r0*OpenMM::NmPerAngstrom);

    force->addBond(bondGroups,bondParameters);
    force->setForceGroup(11);
    if(InOpenMM->system->usesPeriodicBoundaryConditions()){
        force->setUsesPeriodicBoundaryConditions(true);    
    } 
    InOpenMM->system->addForce(force);
    //reinitialize context because a new force was added
    InOpenMM->context->reinitialize();
    printf("addInterChainCentroidBondForce %s with k=%f r0=%f\n",force->getEnergyFunction().c_str(),k,r0);
    return 0;
}

int minimize(OpenMM::Context *context,const double tolerance=10,int maxIterations=0){
    OpenMM::LocalEnergyMinimizer minimizer=OpenMM::LocalEnergyMinimizer();
    minimizer.minimize(*context,tolerance,maxIterations);
    return 0;
}

int Molecule2Molfile(const Molecule *InMol,const OpenMMParameters* InParams,
        molfile_timestep_t *ts){
    int natoms=InMol->NumOfAtoms;
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        ts->coords[3*ii+0]=InMol->X[ii];
        ts->coords[3*ii+1]=InMol->Y[ii];
        ts->coords[3*ii+2]=InMol->Z[ii];
    }   
    //get box lengths in angstrom
    ts->A=InParams->bx;
    ts->B=InParams->by;
    ts->C=InParams->bz;
    //currently only cubix box was supported
    ts->alpha=90.0;
    ts->beta=90.0;
    ts->gamma=90.0;
    return 0;
}

#endif
