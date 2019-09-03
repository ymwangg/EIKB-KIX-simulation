#ifndef __OPENMM_REPD_H__
#define __OPENMM_REPD_H__
#include "OpenMM.h"
#include "DataStructure.h"
#include "OpenMMDataStructure.h"
#include "DataType.h"
#include "OpenMMKernel.h"
using OpenMM::Vec3;
#include <mpi.h>

//#define testOpenMMRepd

int primaryProcSendRecv(const int receiver, const double potentialE, 
        Molecule *InMol, OpenMMData* InOpenMM, const OpenMMParameters* InParams);

int secondaryProcSendRecv(const int receiver, const double potentialE, 
        Molecule *InMol, OpenMMData* InOpenMM, const OpenMMParameters* InParams);

bool queryExchange(float T1, float T2, double E1, double E2);

bool queryExchangeHamiltonian(double E1, double E2);

int setupOpenMMRepdTemperature(Molecule* InMol, OpenMMParameters* InParams,
        OpenMMData* InOpenMM,const RepdParameters *InRepdParams);

int setupOpenMMRepdHamiltonian(Molecule* InMol, OpenMMParameters* InParams,
        OpenMMData* InOpenMM,const RepdParameters *InRepdParams);

int simulateOpenMMRepdTemperature(Molecule* InMol, OpenMMParameters* InParams,
        OpenMMData* InOpenMM,const RepdParameters *InRepdParams){
    /*initialize mpi*/
    int pid;
    int nprocs;
    double startTime,endTime;
    MPI_Comm_rank(MPI_COMM_WORLD,&pid);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    printf("pid=%d was initialized and using platform %s\n",pid,
            InOpenMM->context->getPlatform().getName().c_str());
    srand(time(NULL)+pid);

    setupOpenMMRepdTemperature(InMol, InParams, InOpenMM, InRepdParams);

    /*dcd file handle*/
    /*ts is needed as an intermediate file between Molecule and dcdfile*/
    char outname[128];
    sprintf(outname,"%d.dcd",pid);
    //void *v = open_dcd_write(InParams->outDCDName.c_str(), "dcd", InMol->NumOfAtoms);
    void *v = open_dcd_write(outname, "dcd", InMol->NumOfAtoms);
    dcdhandle *dcd;
    molfile_timestep_t ts;
    ts.coords = (float *)malloc(3*sizeof(float)*InMol->NumOfAtoms);
    dcd = (dcdhandle *)v;

    /*simulation log file handle*/
    char outlogname[128];
    sprintf(outlogname,"%d.log",pid);
    FILE *logfp = fopen(outlogname,"w");
    //write log files to stdout
    //logfp = stdout; 
    InParams->logFileHandle = logfp; 

    /*repd log file handle*/
    char outrepdlogname[128];
    sprintf(outrepdlogname,"%d.repd",pid);
    FILE *repdlogfp = fopen(outrepdlogname,"w");
    //write log files to stdout
    //logfp = stdout; 
    InParams->repdFileHandle = repdlogfp; 

    /*set what information should be retrived from each timestep*/
    int infomask = OpenMM::State::Energy + OpenMM::State::Positions;

    /*set integrator temperature*/
    //float temp = dynamic_cast<OpenMM::LangevinIntegrator*>(InOpenMM->integrator)->getTemperature();
    //printf("temperature=%f\n",temp);
    OpenMM::State state = InOpenMM->context->getState(infomask,true);

    //calculate the number of steps to save as dcd file
    unsigned long int querySteps=0;
    unsigned long int stepsPerQuery=0;
    //check if do exchange or not
    bool doexchange=true;
    if(InRepdParams->exchangeFreq<=0 || 
            InRepdParams->exchangeFreq>=InParams->steps){
        doexchange=false;
    }
    if(InParams->printFreq >= InRepdParams->exchangeFreq){
        querySteps=ceil(float(InParams->steps)/float(InRepdParams->exchangeFreq));
        stepsPerQuery=InRepdParams->exchangeFreq;
    }else{
        querySteps=ceil(float(InParams->steps)/float(InParams->printFreq));
        stepsPerQuery=InParams->printFreq;
    }

    unsigned long int currentTime=0;
    unsigned long int timestep=0;
    unsigned long int exchstep=0;
    InParams->currentTime=0;

    double potentialE;
    double kineticE;

    printf("***************************************************\n");
    dumpOpenMMEnergy(InOpenMM);
    printf("***************************************************\n");
    startTime=MPI_Wtime();

    //main dynamics loop
    for(int ii=0; ii<querySteps; ++ii){
        InOpenMM->integrator->step(stepsPerQuery);
        state=InOpenMM->context->getState(infomask,true);
        std::vector<Vec3> positions = state.getPositions();
        for(int n=0; n<InMol->NumOfAtoms; ++n){
            InMol->X[n]=positions[n][0]*OpenMM::AngstromsPerNm;
            InMol->Y[n]=positions[n][1]*OpenMM::AngstromsPerNm;
            InMol->Z[n]=positions[n][2]*OpenMM::AngstromsPerNm;
        }
        potentialE=state.getPotentialEnergy()*OpenMM::KcalPerKJ;
        kineticE=state.getKineticEnergy()*OpenMM::KcalPerKJ;
        /*save this time frame to dcd file*/
        if(currentTime%InParams->printFreq==0){
            fprintf((FILE*)InParams->logFileHandle,"step=%d V=%f Ek=%f\n",
                    timestep,potentialE,kineticE);
            printf("step=%d V=%f Ek=%f\n",timestep,potentialE,kineticE);
            //move back to original box
            moveBackToOriginalBox(InMol, InParams->bx, InParams->by, InParams->bz);
            resetPositions(InMol,InOpenMM);
            Molecule2Molfile(InMol,InParams,&ts);
            write_timestep(v,&ts); 
            timestep++;
        }
        /*do replica exchange*/
        if(currentTime%InRepdParams->exchangeFreq==0 && doexchange){
            double neighborE=0;
            int exchange;
            //exchange even
            if(exchstep%2==0){
                //even nodes and not the last one
                if(pid%2==0 && pid!=(nprocs-1)){
                    //get energy from pid+1 neighbor
                    MPI_Recv(&neighborE, 1, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    //decide if the exchange succeeds
                    exchange=queryExchange(InRepdParams->temperatures[pid],
                            InRepdParams->temperatures[pid+1],potentialE,neighborE);
                    //send the exchange status to pid+1 neighbor
                    MPI_Send(&exchange,1, MPI_INT, pid+1, 1, MPI_COMM_WORLD);
                    if(exchange){
                        fprintf((FILE*)InParams->repdFileHandle,"step%d exchange(%d,%d)\n",
                                InParams->currentTime,pid,pid+1);
#ifdef testOpenMMRepd
                        printf("T1=%f T2=%f E1=%f E2=%f %d exchange with %d\n",
                                InRepdParams->temperatures[pid],
                                InRepdParams->temperatures[pid+1],
                                potentialE,neighborE,
                                pid,pid+1);
                        float x0,y0,z0;
                        int last=InMol->NumOfAtoms-1;
                        x0=InMol->X[last];
                        y0=InMol->Y[last];
                        z0=InMol->Z[last];
#endif
                        MPI_Send(InMol->X, InMol->NumOfAtoms, MPI_FLOAT, pid+1, 
                                2, MPI_COMM_WORLD);
                        MPI_Send(InMol->Y, InMol->NumOfAtoms, MPI_FLOAT, pid+1,
                                3, MPI_COMM_WORLD);
                        MPI_Send(InMol->Z, InMol->NumOfAtoms, MPI_FLOAT, pid+1,
                                4, MPI_COMM_WORLD);
                        MPI_Recv(InMol->X, InMol->NumOfAtoms, MPI_FLOAT, pid+1,
                                2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(InMol->Y, InMol->NumOfAtoms, MPI_FLOAT, pid+1,
                                3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(InMol->Z, InMol->NumOfAtoms, MPI_FLOAT, pid+1,
                                4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        resetPositionsVelocities(InMol,InParams,InOpenMM);
#ifdef testOpenMMRepd
                        printf("pid%d (%f,%f,%f) -> (%f,%f,%f)\n",pid,x0,y0,z0,InMol->X[last],
                                InMol->Y[last],InMol->Z[last]);
#endif
                    }else if(pid!=(nprocs-1)){
#ifdef printFail
                        printf("T1=%f T2=%f E1=%f E2=%f %d failed exchange with %d\n",
                                InRepdParams->temperatures[pid],
                                InRepdParams->temperatures[pid+1],
                                potentialE,neighborE,
                                pid,pid+1);
#endif
                    }
                } //end even nodes
                //odd nodes
                else if(pid%2==1){
                    //send energy to pid-1 neighbor
                    MPI_Send(&potentialE, 1, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD);
                    //receive the exchange status from pid-1 neighbor
                    MPI_Recv(&exchange, 1, MPI_INT, pid-1, 1 ,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    if(exchange){
#ifdef testOpenMMRepd
                        float x0,y0,z0;
                        int last=InMol->NumOfAtoms-1;
                        x0=InMol->X[last];
                        y0=InMol->Y[last];
                        z0=InMol->Z[last];
#endif
                        MPI_Send(InMol->X, InMol->NumOfAtoms, MPI_FLOAT, pid-1,
                                2, MPI_COMM_WORLD);
                        MPI_Send(InMol->Y, InMol->NumOfAtoms, MPI_FLOAT, pid-1,
                                3, MPI_COMM_WORLD);
                        MPI_Send(InMol->Z, InMol->NumOfAtoms, MPI_FLOAT, pid-1,
                                4, MPI_COMM_WORLD);
                        MPI_Recv(InMol->X, InMol->NumOfAtoms, MPI_FLOAT, pid-1,
                                2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(InMol->Y, InMol->NumOfAtoms, MPI_FLOAT, pid-1,
                                3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(InMol->Z, InMol->NumOfAtoms, MPI_FLOAT, pid-1,
                                4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        resetPositionsVelocities(InMol,InParams,InOpenMM);
#ifdef testOpenMMRepd
                        printf("pid%d (%f,%f,%f) -> (%f,%f,%f)\n",pid,x0,y0,z0,InMol->X[last],
                                InMol->Y[last],InMol->Z[last]);
#endif
                    }
                } //end off nodes
            } //endif even exchange steps

            else{
                //odd nodes and not the last one
                if(pid%2==1 && pid!=(nprocs-1)){
                    //get energy from pid+1 neighbor
                    MPI_Recv(&neighborE, 1, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    //decide if the exchange succeeds
                    exchange=queryExchange(InRepdParams->temperatures[pid],
                            InRepdParams->temperatures[pid+1],potentialE,neighborE);
                    //send the exchange status to pid+1 neighbor
                    MPI_Send(&exchange,1, MPI_INT, pid+1, 1, MPI_COMM_WORLD);
                    if(exchange){
                        fprintf((FILE*)InParams->repdFileHandle,"step%d exchange(%d,%d)\n",
                                InParams->currentTime,pid,pid+1);
#ifdef testOpenMMRepd
                        printf("T1=%f T2=%f E1=%f E2=%f %d exchange with %d\n",
                                InRepdParams->temperatures[pid],
                                InRepdParams->temperatures[pid+1],
                                potentialE,neighborE,
                                pid,pid+1);
                        float x0,y0,z0;
                        int last=InMol->NumOfAtoms-1;
                        x0=InMol->X[last];
                        y0=InMol->Y[last];
                        z0=InMol->Z[last];
#endif
                        MPI_Send(InMol->X, InMol->NumOfAtoms, MPI_FLOAT, pid+1, 
                                2, MPI_COMM_WORLD);
                        MPI_Send(InMol->Y, InMol->NumOfAtoms, MPI_FLOAT, pid+1,
                                3, MPI_COMM_WORLD);
                        MPI_Send(InMol->Z, InMol->NumOfAtoms, MPI_FLOAT, pid+1,
                                4, MPI_COMM_WORLD);
                        MPI_Recv(InMol->X, InMol->NumOfAtoms, MPI_FLOAT, pid+1,
                                2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(InMol->Y, InMol->NumOfAtoms, MPI_FLOAT, pid+1,
                                3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(InMol->Z, InMol->NumOfAtoms, MPI_FLOAT, pid+1,
                                4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        resetPositionsVelocities(InMol,InParams,InOpenMM);
#ifdef testOpenMMRepd
                        printf("pid%d (%f,%f,%f) -> (%f,%f,%f)\n",pid,x0,y0,z0,InMol->X[last],
                                InMol->Y[last],InMol->Z[last]);
#endif
                    } else {
#ifdef printFail
                        printf("T1=%f T2=%f E1=%f E2=%f %d failed exchange with %d\n",
                                InRepdParams->temperatures[pid],
                                InRepdParams->temperatures[pid+1],
                                potentialE,neighborE,
                                pid,pid+1);
#endif
                    }
                }
                //even node
                else if(pid%2==0 && pid!=0){
                    //send energy to pid-1 neighbor
                    MPI_Send(&potentialE, 1, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD);
                    //receive the exchange status from pid-1 neighbor
                    MPI_Recv(&exchange, 1, MPI_INT, pid-1, 1 ,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    if(exchange){
#ifdef testOpenMMRepd
                        float x0,y0,z0;
                        int last=InMol->NumOfAtoms-1;
                        x0=InMol->X[last];
                        y0=InMol->Y[last];
                        z0=InMol->Z[last];
#endif
                        MPI_Send(InMol->X, InMol->NumOfAtoms, MPI_FLOAT, pid-1,
                                2, MPI_COMM_WORLD);
                        MPI_Send(InMol->Y, InMol->NumOfAtoms, MPI_FLOAT, pid-1,
                                3, MPI_COMM_WORLD);
                        MPI_Send(InMol->Z, InMol->NumOfAtoms, MPI_FLOAT, pid-1,
                                4, MPI_COMM_WORLD);
                        MPI_Recv(InMol->X, InMol->NumOfAtoms, MPI_FLOAT, pid-1,
                                2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(InMol->Y, InMol->NumOfAtoms, MPI_FLOAT, pid-1,
                                3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(InMol->Z, InMol->NumOfAtoms, MPI_FLOAT, pid-1,
                                4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        resetPositionsVelocities(InMol,InParams,InOpenMM);
#ifdef testOpenMMRepd
                        printf("pid%d (%f,%f,%f) -> (%f,%f,%f)\n",pid,x0,y0,z0,InMol->X[last],
                                InMol->Y[last],InMol->Z[last]);
#endif
                    }
                }
            } //end if odd exchange steps
            exchstep++;
        } //endif exchange

        currentTime += stepsPerQuery;
        InParams->currentTime = currentTime;
        fflush((FILE*)InParams->logFileHandle);
        fflush((FILE*)InParams->repdFileHandle);
        
    } //end dynamics for loop

    //close log file handles
    if((FILE*)InParams->logFileHandle!=stdout){
        fclose((FILE*)InParams->logFileHandle);
    }
    if((FILE*)InParams->repdFileHandle!=stdout){
        fclose((FILE*)InParams->repdFileHandle);
    }

    //close dcd file handle
    close_file_write(v);

    //close MPI
    endTime=MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    if(pid==0){
        printf("Done simulation\n");
        printf("Elapsed time = %f seconds\n",endTime-startTime);
    }
    MPI_Finalize();
    return 0;
}

int simulateOpenMMRepdHamiltonian(Molecule* InMol, OpenMMParameters* InParams,
        OpenMMData* InOpenMM,const RepdParameters *InRepdParams){
    /*initialize mpi*/
    int pid;
    int nprocs;
    double startTime,endTime;
    MPI_Comm_rank(MPI_COMM_WORLD,&pid);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    int dim1=InRepdParams->dim1;
    int dim2=InRepdParams->dim2;
    if(dim1*dim2!=nprocs){
        printf("Warning: Replica Exchange Dim1=%d Dim2=%d didn't "
                "match %d processors\n",dim1,dim2,nprocs);
    }
    printf("pid=%d was initialized and using platform %s\n",pid,
            InOpenMM->context->getPlatform().getName().c_str());
    setupOpenMMRepdHamiltonian(InMol,InParams, InOpenMM,InRepdParams);
    srand(time(NULL)+pid);
    //calculate the replica coordinate (idx1,idx2) in (coordinate1,coordinate2)
    int idx1 = pid % dim1;
    int idx2 = pid / dim1;
    InParams->idx1=idx1;
    InParams->idx2=idx2;


    /*dcd file handle*/
    /*ts is needed as an intermediate file between Molecule and dcdfile*/
    char outname[128];
    sprintf(outname,"%d-%d.dcd",idx1,idx2);
    void *v = open_dcd_write(outname, "dcd", InMol->NumOfAtoms);
    dcdhandle *dcd;
    molfile_timestep_t ts;
    ts.coords = (float *)malloc(3*sizeof(float)*InMol->NumOfAtoms);
    dcd = (dcdhandle *)v;

    /*simulation log file handle*/
    char outlogname[128];
    sprintf(outlogname,"%d-%d.log",idx1,idx2);
    FILE *logfp = fopen(outlogname,"w");
    //write log files to stdout
    //logfp = stdout; 
    InParams->logFileHandle = logfp; 

    /*repd log file handle*/
    char outrepdlogname[128];
    sprintf(outrepdlogname,"%d-%d.repd",idx1,idx2);
    FILE *repdlogfp = fopen(outrepdlogname,"w");
    //write log files to stdout
    //logfp = stdout; 
    InParams->repdFileHandle = repdlogfp; 

    /*set what information should be retrived from each timestep*/
    int infomask = OpenMM::State::Energy + OpenMM::State::Positions;

    /*set integrator temperature*/
    //float temp = dynamic_cast<OpenMM::LangevinIntegrator*>(InOpenMM->integrator)->getTemperature();
    //printf("temperature=%f\n",temp);
    OpenMM::State state = InOpenMM->context->getState(infomask,true);
    //calculate the number of steps to save as dcd file
    unsigned long int querySteps=0;
    unsigned long int stepsPerQuery=0;
    bool doexchange=true;
    if(InRepdParams->exchangeFreq<=0 || 
            InRepdParams->exchangeFreq>=InParams->steps){
        doexchange=false;
    }
    if(InParams->printFreq >= InRepdParams->exchangeFreq){
        querySteps=ceil(float(InParams->steps)/float(InRepdParams->exchangeFreq));
        stepsPerQuery=InRepdParams->exchangeFreq;
    }else{
        querySteps=ceil(float(InParams->steps)/float(InParams->printFreq));
        stepsPerQuery=InParams->printFreq;
    }
    unsigned long int currentTime=0;
    unsigned long int timestep=0;
    unsigned long int exchstep=0;
    InParams->currentTime=0;

    double potentialE;
    double kineticE;

    printf("***************************************************\n");
    printf("pid=%d was running at %fK\n",pid,
            dynamic_cast<OpenMM::LangevinIntegrator*>(
                InOpenMM->integrator)->getTemperature());
    dumpOpenMMEnergy(InOpenMM);
    printf("***************************************************\n");
    startTime=MPI_Wtime();

    //main dynamics loop
    for(int ii=0; ii<querySteps; ++ii){
        InOpenMM->integrator->step(stepsPerQuery);
        state=InOpenMM->context->getState(infomask,true);
        std::vector<Vec3> positions = state.getPositions();
        for(int n=0; n<InMol->NumOfAtoms; ++n){
            InMol->X[n]=positions[n][0]*OpenMM::AngstromsPerNm;
            InMol->Y[n]=positions[n][1]*OpenMM::AngstromsPerNm;
            InMol->Z[n]=positions[n][2]*OpenMM::AngstromsPerNm;
        }
        potentialE=state.getPotentialEnergy()*OpenMM::KcalPerKJ;
        kineticE=state.getKineticEnergy()*OpenMM::KcalPerKJ;
        /*save this time frame to dcd file*/
        if(currentTime%InParams->printFreq==0){
            fprintf((FILE*)InParams->logFileHandle,"step=%d V=%f Ek=%f\n",
                    timestep,potentialE,kineticE);
            printf("step=%d V=%f Ek=%f\n",timestep,potentialE,kineticE);
            //move back to original box
            moveBackToOriginalBox(InMol, InParams->bx, InParams->by, InParams->bz);
            resetPositions(InMol,InOpenMM);
            Molecule2Molfile(InMol,InParams,&ts);
            write_timestep(v,&ts); 
            timestep++;
        }
        /*do replica exchange*/
        if(currentTime%InRepdParams->exchangeFreq==0 && doexchange){
            //exchange 1d even
            if(exchstep%4==0){
                //even nodes and not the last one
                if(idx1%2==0 && idx1!=(dim1-1)){
                    primaryProcSendRecv(pid+1,potentialE,InMol,InOpenMM,InParams);
                } //end even nodes
                //odd nodes
                else if(idx1%2==1){
                    secondaryProcSendRecv(pid-1,potentialE,InMol,InOpenMM,InParams);
                } //end off nodes
            } //endif even exchange steps

            //exchange 1d odd
            else if (exchstep%4==1) {
                //odd nodes and not the last one
                if(idx1%2==1 && idx1!=(dim1-1)){
                    primaryProcSendRecv(pid+1,potentialE,InMol,InOpenMM,InParams);
                }
                //even node
                else if(idx1%2==0 && idx1!=0){
                    secondaryProcSendRecv(pid-1,potentialE,InMol,InOpenMM,InParams);
                }
            } //end if odd exchange steps

            //exchange 2d even
            else if(exchstep%4==2){
                //even nodes and not the last one
                if(idx2%2==0 && idx2!=(dim2-1)){
                    primaryProcSendRecv(pid+dim1,potentialE,InMol,InOpenMM,InParams);
                } //end even nodes
                //odd nodes
                else if(idx2%2==1){
                    secondaryProcSendRecv(pid-dim1,potentialE,InMol,InOpenMM,InParams);
                } //end off nodes
            } //endif even exchange steps

            //exchange 2d odd
            else if (exchstep%4==3) {
                //odd nodes and not the last one
                if(idx2%2==1 && idx2!=(dim2-1)){
                    primaryProcSendRecv(pid+dim1,potentialE,InMol,InOpenMM,InParams);
                }
                //even node
                else if(idx2%2==0 && idx2!=0){
                    secondaryProcSendRecv(pid-dim1,potentialE,InMol,InOpenMM,InParams);
                }
            } //end if odd exchange steps

            exchstep++;
        } //endif exchange

        currentTime += stepsPerQuery;
        InParams->currentTime = currentTime;
        fflush((FILE*)InParams->logFileHandle);
        fflush((FILE*)InParams->repdFileHandle);
    } //end dynamics for loop

    //close log file handles
    if((FILE*)InParams->logFileHandle!=stdout){
        fclose((FILE*)InParams->logFileHandle);
    }
    if((FILE*)InParams->repdFileHandle!=stdout){
        fclose((FILE*)InParams->repdFileHandle);
    }

    //close dcd file handle
    close_file_write(v);

    //close MPI
    endTime=MPI_Wtime();
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    if(pid==0){
        printf("Done simulation\n");
        printf("Elapsed time = %f seconds\n",endTime-startTime);
    }
    MPI_Finalize();
    return 0;
}

int simulateOpenMMRepdHamiltonian3D(Molecule* InMol, OpenMMParameters* InParams,
        OpenMMData* InOpenMM,const RepdParameters *InRepdParams){
    /*initialize mpi*/
    int pid;
    int nprocs;
    double startTime,endTime;
    MPI_Comm_rank(MPI_COMM_WORLD,&pid);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    int dim1=InRepdParams->dim1;
    int dim2=InRepdParams->dim2;
    int dim3=InRepdParams->dim3;
    if(dim1*dim2*dim3!=nprocs){
        printf("Warning: Replica Exchange Dim1=%d Dim2=%d didn't "
                "match %d processors\n",dim1,dim2,nprocs);
    }
    printf("pid=%d was initialized and using platform %s\n",pid,
            InOpenMM->context->getPlatform().getName().c_str());
    setupOpenMMRepdHamiltonian(InMol,InParams, InOpenMM,InRepdParams);
    srand(time(NULL)+pid);
    //calculate the replica coordinate (idx1,idx2,idx3) in (coordinate1,coordinate2)
    int idx1 = pid % (dim1*dim2) % dim1;
    int idx2 = pid % (dim1*dim2) / dim1;
    int idx3 = pid / (dim1*dim2);
    InParams->idx1=idx1;
    InParams->idx2=idx2;
    InParams->idx3=idx3;


    /*dcd file handle*/
    /*ts is needed as an intermediate file between Molecule and dcdfile*/
    char outname[128];
    sprintf(outname,"%d-%d-%d.dcd",idx1,idx2,idx3);
    void *v = open_dcd_write(outname, "dcd", InMol->NumOfAtoms);
    dcdhandle *dcd;
    molfile_timestep_t ts;
    ts.coords = (float *)malloc(3*sizeof(float)*InMol->NumOfAtoms);
    dcd = (dcdhandle *)v;

    /*simulation log file handle*/
    char outlogname[128];
    sprintf(outlogname,"%d-%d-%d.log",idx1,idx2,idx3);
    FILE *logfp = fopen(outlogname,"w");
    //write log files to stdout
    //logfp = stdout; 
    InParams->logFileHandle = logfp; 

    /*repd log file handle*/
    char outrepdlogname[128];
    sprintf(outrepdlogname,"%d-%d-%d.repd",idx1,idx2,idx3);
    FILE *repdlogfp = fopen(outrepdlogname,"w");
    //write log files to stdout
    //logfp = stdout; 
    InParams->repdFileHandle = repdlogfp; 

    /*set what information should be retrived from each timestep*/
    int infomask = OpenMM::State::Energy + OpenMM::State::Positions;

    /*set integrator temperature*/
    //float temp = dynamic_cast<OpenMM::LangevinIntegrator*>(InOpenMM->integrator)->getTemperature();
    //printf("temperature=%f\n",temp);
    OpenMM::State state = InOpenMM->context->getState(infomask,true);
    //calculate the number of steps to save as dcd file
    unsigned long int querySteps=0;
    unsigned long int stepsPerQuery=0;
    bool doexchange=true;
    if(InRepdParams->exchangeFreq<=0 || 
            InRepdParams->exchangeFreq>=InParams->steps){
        doexchange=false;
    }
    if(InParams->printFreq >= InRepdParams->exchangeFreq){
        querySteps=ceil(float(InParams->steps)/float(InRepdParams->exchangeFreq));
        stepsPerQuery=InRepdParams->exchangeFreq;
    }else{
        querySteps=ceil(float(InParams->steps)/float(InParams->printFreq));
        stepsPerQuery=InParams->printFreq;
    }
    unsigned long int currentTime=0;
    unsigned long int timestep=0;
    unsigned long int exchstep=0;
    InParams->currentTime=0;

    double potentialE;
    double kineticE;

    printf("***************************************************\n");
    dumpOpenMMEnergy(InOpenMM);
    printf("***************************************************\n");
    startTime=MPI_Wtime();

    //main dynamics loop
    for(int ii=0; ii<querySteps; ++ii){
        InOpenMM->integrator->step(stepsPerQuery);
        state=InOpenMM->context->getState(infomask,true);
        std::vector<Vec3> positions = state.getPositions();
        for(int n=0; n<InMol->NumOfAtoms; ++n){
            InMol->X[n]=positions[n][0]*OpenMM::AngstromsPerNm;
            InMol->Y[n]=positions[n][1]*OpenMM::AngstromsPerNm;
            InMol->Z[n]=positions[n][2]*OpenMM::AngstromsPerNm;
        }
        potentialE=state.getPotentialEnergy()*OpenMM::KcalPerKJ;
        kineticE=state.getKineticEnergy()*OpenMM::KcalPerKJ;
        /*save this time frame to dcd file*/
        if(currentTime%InParams->printFreq==0){
            fprintf((FILE*)InParams->logFileHandle,"step=%d V=%f Ek=%f\n",
                    timestep,potentialE,kineticE);
            printf("step=%d V=%f Ek=%f\n",timestep,potentialE,kineticE);
            //move back to original box
            moveBackToOriginalBox(InMol, InParams->bx, InParams->by, InParams->bz);
            resetPositions(InMol,InOpenMM);
            Molecule2Molfile(InMol,InParams,&ts);
            write_timestep(v,&ts); 
            timestep++;
        }
        /*do replica exchange*/
        if(currentTime%InRepdParams->exchangeFreq==0 && doexchange){
            //exchange 1d even
            if(exchstep%6==0){
                //even nodes and not the last one
                if(idx1%2==0 && idx1!=(dim1-1)){
                    primaryProcSendRecv(pid+1,potentialE,InMol,InOpenMM,InParams);
                } //end even nodes
                //odd nodes
                else if(idx1%2==1){
                    secondaryProcSendRecv(pid-1,potentialE,InMol,InOpenMM,InParams);
                } //end off nodes
            } //endif even exchange steps

            //exchange 1d odd
            else if (exchstep%6==1) {
                //odd nodes and not the last one
                if(idx1%2==1 && idx1!=(dim1-1)){
                    primaryProcSendRecv(pid+1,potentialE,InMol,InOpenMM,InParams);
                }
                //even node
                else if(idx1%2==0 && idx1!=0){
                    secondaryProcSendRecv(pid-1,potentialE,InMol,InOpenMM,InParams);
                }
            } //end if odd exchange steps

            //exchange 2d even
            else if(exchstep%6==2){
                //even nodes and not the last one
                if(idx2%2==0 && idx2!=(dim2-1)){
                    primaryProcSendRecv(pid+dim1,potentialE,InMol,InOpenMM,InParams);
                } //end even nodes
                //odd nodes
                else if(idx2%2==1){
                    secondaryProcSendRecv(pid-dim1,potentialE,InMol,InOpenMM,InParams);
                } //end off nodes
            } //endif even exchange steps

            //exchange 2d odd
            else if (exchstep%6==3) {
                //odd nodes and not the last one
                if(idx2%2==1 && idx2!=(dim2-1)){
                    primaryProcSendRecv(pid+dim1,potentialE,InMol,InOpenMM,InParams);
                }
                //even node
                else if(idx2%2==0 && idx2!=0){
                    secondaryProcSendRecv(pid-dim1,potentialE,InMol,InOpenMM,InParams);
                }
            } //end if odd exchange steps

            //exchange 3d even
            else if(exchstep%6==2){
                //even nodes and not the last one
                if(idx3%2==0 && idx3!=(dim3-1)){
                    primaryProcSendRecv(pid+dim1*dim2,potentialE,InMol,InOpenMM,InParams);
                } //end even nodes
                //odd nodes
                else if(idx3%2==1){
                    secondaryProcSendRecv(pid-dim1*dim2,potentialE,InMol,InOpenMM,InParams);
                } //end off nodes
            } //endif even exchange steps

            //exchange 3d odd
            else if (exchstep%6==3) {
                //odd nodes and not the last one
                if(idx3%2==1 && idx3!=(dim3-1)){
                    primaryProcSendRecv(pid+dim1*dim2,potentialE,InMol,InOpenMM,InParams);
                }
                //even node
                else if(idx3%2==0 && idx3!=0){
                    secondaryProcSendRecv(pid-dim1*dim2,potentialE,InMol,InOpenMM,InParams);
                }
            } //end if odd exchange steps

            exchstep++;
        } //endif exchange

        currentTime += stepsPerQuery;
        InParams->currentTime = currentTime;
        fflush((FILE*)InParams->logFileHandle);
        fflush((FILE*)InParams->repdFileHandle);
    } //end dynamics for loop

    //close log file handles
    if((FILE*)InParams->logFileHandle!=stdout){
        fclose((FILE*)InParams->logFileHandle);
    }
    if((FILE*)InParams->repdFileHandle!=stdout){
        fclose((FILE*)InParams->repdFileHandle);
    }

    //close dcd file handle
    close_file_write(v);

    //close MPI
    endTime=MPI_Wtime();
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    if(pid==0){
        printf("Done simulation\n");
        printf("Elapsed time = %f seconds\n",endTime-startTime);
    }
    MPI_Finalize();
    return 0;
}
int primaryProcSendRecv(const int receiver, const double potentialE, 
        Molecule *InMol, OpenMMData* InOpenMM, 
        const OpenMMParameters* InParams){ 
    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD,&pid);
    float *swap;
    double E1=potentialE;
    double E2=0;
    double deltaE1=0;
    double deltaE2=0;
    float T=InParams->temperature;
    bool exchange;
    MPI_Send(InMol->X, InMol->NumOfAtoms, MPI_FLOAT, receiver, 
            1, MPI_COMM_WORLD);
    MPI_Send(InMol->Y, InMol->NumOfAtoms, MPI_FLOAT, receiver,
            2, MPI_COMM_WORLD);
    MPI_Send(InMol->Z, InMol->NumOfAtoms, MPI_FLOAT, receiver,
            3, MPI_COMM_WORLD);
    //get coordinates from idx+1 neighbor
    MPI_Recv(InMol->tmpX, InMol->NumOfAtoms, MPI_FLOAT, receiver,
            1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(InMol->tmpY, InMol->NumOfAtoms, MPI_FLOAT, receiver,
            2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(InMol->tmpZ, InMol->NumOfAtoms, MPI_FLOAT, receiver,
            3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    swap=InMol->X; InMol->X=InMol->tmpX; InMol->tmpX=swap;
    swap=InMol->Y; InMol->Y=InMol->tmpY; InMol->tmpY=swap;
    swap=InMol->Z; InMol->Z=InMol->tmpZ; InMol->tmpZ=swap;

    resetPositions(InMol,InOpenMM);
    OpenMM::State state=InOpenMM->context->getState(OpenMM::State::Energy,true);
    E2=state.getPotentialEnergy()*OpenMM::KcalPerKJ;
    deltaE1 = (E2-E1)/T;
    MPI_Recv(&deltaE2, 1, MPI_DOUBLE, receiver, 
            4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    exchange=queryExchangeHamiltonian(deltaE1,deltaE2);
    MPI_Send(&exchange, 1, MPI_INT, receiver, 0, MPI_COMM_WORLD);

    if(exchange){
#ifdef testOpenMMRepd
        printf("E1=%f E2=%f dE1=%f dE2=%f %d exchange with %d\n",
                E1,E2,
                deltaE1,deltaE2,pid,receiver);
#endif
        fprintf((FILE*)InParams->repdFileHandle,"step=%d exchange(%d,%d)\n",
                InParams->currentTime,pid,receiver);

        //reset velocities
        resetVelocities(InParams,InOpenMM);
#ifdef testOpenMMRepd
        float x0,y0,z0;
        int last=InMol->NumOfAtoms-1;
        last=0;
        x0=InMol->tmpX[last];
        y0=InMol->tmpY[last];
        z0=InMol->tmpZ[last];
        printf("pid%d (%f,%f,%f) -> (%f,%f,%f)\n",pid,x0,y0,z0,
                InMol->X[last],InMol->Y[last],InMol->Z[last]);
#endif
    }else{
        swap=InMol->X; InMol->X=InMol->tmpX; InMol->tmpX=swap;
        swap=InMol->Y; InMol->Y=InMol->tmpY; InMol->tmpY=swap;
        swap=InMol->Z; InMol->Z=InMol->tmpZ; InMol->tmpZ=swap;
        resetPositions(InMol,InOpenMM);
#ifdef printFail
        printf("E1=%f E2=%f dE1=%f dE2=%f %d failed exchange with %d\n",
                E1,E2,
                deltaE1,deltaE2,pid,receiver);
#endif
    }
    return 0;
}

int secondaryProcSendRecv(const int receiver, const double potentialE, 
        Molecule *InMol, OpenMMData* InOpenMM, const OpenMMParameters* InParams){
    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD,&pid);
    float *swap;
    double E1=potentialE;
    double E2=0;
    double deltaE2=0;
    float T=InParams->temperature;
    bool exchange;
    MPI_Send(InMol->X, InMol->NumOfAtoms, MPI_FLOAT, receiver, 
            1, MPI_COMM_WORLD);
    MPI_Send(InMol->Y, InMol->NumOfAtoms, MPI_FLOAT, receiver,
            2, MPI_COMM_WORLD);
    MPI_Send(InMol->Z, InMol->NumOfAtoms, MPI_FLOAT, receiver,
            3, MPI_COMM_WORLD);
    //get coordinates from idx+1 neighbor
    MPI_Recv(InMol->tmpX, InMol->NumOfAtoms, MPI_FLOAT, receiver,
            1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(InMol->tmpY, InMol->NumOfAtoms, MPI_FLOAT, receiver,
            2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(InMol->tmpZ, InMol->NumOfAtoms, MPI_FLOAT, receiver,
            3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    swap=InMol->X; InMol->X=InMol->tmpX; InMol->tmpX=swap;
    swap=InMol->Y; InMol->Y=InMol->tmpY; InMol->tmpY=swap;
    swap=InMol->Z; InMol->Z=InMol->tmpZ; InMol->tmpZ=swap;

    resetPositions(InMol,InOpenMM);
    OpenMM::State state=InOpenMM->context->getState(OpenMM::State::Energy,true);
    E2=state.getPotentialEnergy()*OpenMM::KcalPerKJ;
    deltaE2 = (E2-E1)/T;
    MPI_Send(&deltaE2, 1, MPI_DOUBLE, receiver,
            4, MPI_COMM_WORLD);
    MPI_Recv(&exchange, 1, MPI_INT, receiver, 0, MPI_COMM_WORLD,
            MPI_STATUS_IGNORE);

    if(exchange){
#ifdef testOpenMMRepd
        float x0,y0,z0;
        int last=InMol->NumOfAtoms-1;
        last=0;
        x0=InMol->tmpX[last];
        y0=InMol->tmpY[last];
        z0=InMol->tmpZ[last];
#endif
        resetVelocities(InParams,InOpenMM);
#ifdef testOpenMMRepd
        printf("E1=%f E2=%f pid%d (%f,%f,%f) -> (%f,%f,%f)\n",
                E1,E2,
                pid,x0,y0,z0,InMol->X[last],
                InMol->Y[last],InMol->Z[last]);
#endif
    }else{
#ifdef printFail
        printf("E1=%f E2=%f pid%d failed exchange with pid%d\n",
                E1,E2,pid,receiver);
#endif
        swap=InMol->X; InMol->X=InMol->tmpX; InMol->tmpX=swap;
        swap=InMol->Y; InMol->Y=InMol->tmpY; InMol->tmpY=swap;
        swap=InMol->Z; InMol->Z=InMol->tmpZ; InMol->tmpZ=swap;
        resetPositions(InMol,InOpenMM);
    }
    return 0;
}

int setupOpenMMRepdTemperature(Molecule* InMol, OpenMMParameters* InParams,
        OpenMMData* InOpenMM,const RepdParameters *InRepdParams){
    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD,&pid);
    /*set each replica's temperature*/
    //InParams->temperature=InRepdParams->temperatures[pid];
    InParams->temperature=300;
    //reset positions and scale velocities to InParams->temperature
    resetPositionsVelocities(InMol,InParams,InOpenMM);
    return 0;
}

int setupOpenMMRepdHamiltonian(Molecule* InMol, OpenMMParameters* InParams,
        OpenMMData* InOpenMM,const RepdParameters *InRepdParams){
    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD,&pid);
    int dim1=InRepdParams->dim1;
    int dim2=InRepdParams->dim2;
    int idx1 = pid % dim1;
    int idx2 = pid / dim1;
    /*modify coordinate*/

    /*set each replica's temperature*/
    InParams->temperature=300+idx2*20;
    //InParams->temperature=300+idx1*1;
    //InParams->temperature=300;
    
    /*add forces*/

    //force constant
    //double k=5; //kcal/mol/A
    //distance interval
    //double deltaR=0.5; //angstrom
    
    //do something with idx1
    /*addInterChainCentroidBondForce(InMol,InOpenMM,'A','B',k,(idx1)*deltaR);
    if(idx1*deltaR<20){
        setChainToPosition(InMol,'B',InParams->bx/2.0+20,InParams->by/2.0,
                InParams->bz/2.0);
    } else{
        setChainToPosition(InMol,'B',InParams->bx/2.0+(idx1)*deltaR,
                InParams->by/2.0,
                InParams->bz/2.0);
    }*/
    //float scaleFactor=idx1*0.02+1.38; //HIF init
    float scaleFactor=idx1*0.02+1.10; //CITED2 init
    //float scaleFactor=idx1*0.02+1.7;
    //float scaleFactor=idx1*0.02+1.3; //MLL first time
    float scaleFactor=idx1*0.02+1.24; //MLL
    //float scaleFactor=idx1*0.02+1.76; //pKID first time
    //float scaleFactor=idx1*0.02+1.82; //pKID
    scaleNative(InMol,{{'A','B',scaleFactor}});
    resetNativeForce(InMol,InOpenMM);

    //do someting with idx2
    /*
    addInterChainCentroidBondForce(InMol,InOpenMM,'A','C',k,(idx2)*deltaR);
    if(idx2*deltaR<20){
        setChainToPosition(InMol,'C',InParams->bx/2.0-20,InParams->by/2.0,
                InParams->bz/2.0);
    }else{
        setChainToPosition(InMol,'C',InParams->bx/2.0-(idx2)*deltaR,
                InParams->by/2.0,
                InParams->bz/2.0);
    }*/

    resetPositionsVelocities(InMol,InParams,InOpenMM);

    /*minimize*/
    OpenMM::LocalEnergyMinimizer minimizer=OpenMM::LocalEnergyMinimizer();
    minimizer.minimize(*InOpenMM->context);
    return 0;
}

bool queryExchangeHamiltonian(double E1, double E2){
    //T1,T2 kelvein
    //E1,E2 kcal/mol
    double KBOLTZ=1.987191e-3;
    //E1=deltaE1=(Ei(ri)-Ei(rj))
    //E2=deltaE2=(Ej(rj)-Ej(ri))
    //if E1+E2 < 0, probability = 1
    //else probability = exp(-1.0*(E1+E2))
    double delta = 1.0 / KBOLTZ  * (E1 + E2);
    double probability = fmin(1.0, exp(-1.0*delta));
    if(probability==1){
        return true;
    }else{
        double random=rand()/double(RAND_MAX);
        if(random <= probability){
            return true;
        }
        else{
            return false;
        }
    }
}


bool queryExchange(float T1, float T2, double E1, double E2){
    //T1,T2 kelvein
    //E1,E2 kcal/mol
    double KBOLTZ=1.987191e-3;
    double delta = -1.0 / KBOLTZ * (1.0/T1 - 1.0/T2) * (E1 - E2);
    double probability = fmin(1.0, exp(-1.0*delta));
    if(probability==1){
        return true;
    }else{
        double random=rand()/double(RAND_MAX);
        if(random <= probability){
            return true;
        }
        else{
            return false;
        }
    }
}

int setupOpenMMRepdHamiltonian3D(Molecule* InMol, OpenMMParameters* InParams,
        OpenMMData* InOpenMM,const RepdParameters *InRepdParams){
    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD,&pid);
    int dim1=InRepdParams->dim1;
    int dim2=InRepdParams->dim2;
    int dim3=InRepdParams->dim3;
    int idx1 = pid % (dim1*dim2) % dim1;
    int idx2 = pid % (dim1*dim2) / dim1;
    int idx3 = pid / (dim1*dim2);
    /*modify coordinate*/

    /*set each replica's temperature*/
    InParams->temperature=300+idx1*15;
    //InParams->temperature=300;
    
    /*add forces*/

    //force constant
    //double k=5; //kcal/mol/A
    //distance interval
    //double deltaR=0.5; //angstrom
    
    //do something with idx1
    /*
    addInterChainCentroidBondForce(InMol,InOpenMM,'A','B',k,(idx1)*deltaR);
    if(idx1*deltaR<20){
        setChainToPosition(InMol,'B',InParams->bx/2.0+20,InParams->by/2.0,
                InParams->bz/2.0);
    } else{
        setChainToPosition(InMol,'B',InParams->bx/2.0+(idx1)*deltaR,
                InParams->by/2.0,
                InParams->bz/2.0);
    }*/
    float scaleFactor=idx1*0.02+0.8;
    scaleNative(InMol,{{'A','B',scaleFactor}});
    resetNativeForce(InMol,InOpenMM);

    //do someting with idx2
    /*
    addInterChainCentroidBondForce(InMol,InOpenMM,'A','C',k,(idx2)*deltaR);
    if(idx2*deltaR<20){
        setChainToPosition(InMol,'C',InParams->bx/2.0-20,InParams->by/2.0,
                InParams->bz/2.0);
    }else{
        setChainToPosition(InMol,'C',InParams->bx/2.0-(idx2)*deltaR,
                InParams->by/2.0,
                InParams->bz/2.0);
    }*/

    resetPositionsVelocities(InMol,InParams,InOpenMM);

    /*minimize*/
    OpenMM::LocalEnergyMinimizer minimizer=OpenMM::LocalEnergyMinimizer();
    minimizer.minimize(*InOpenMM->context);
    return 0;
}
#endif
