#ifndef __BONDED_H__
#define __BONDED_H__
#include "DataStructure.h"
#include "DataType.h"
#include "Vector.h"
#define DEBUG_BOND true
int computeBonds(const Molecule *InMolecule, State *OutState){
    float energy=0.0;
    int numOfBonds=InMolecule->NumOfBonds;
    for(int ii=0; ii<numOfBonds; ++ii){
        int atomIdx1=InMolecule->BOND[ii].i;
        int atomIdx2=InMolecule->BOND[ii].j;
        float r0=InMolecule->BOND[ii].r;
        float k=InMolecule->BOND[ii].k;
        float x1,y1,z1;
        x1=OutState->X[atomIdx1];
        y1=OutState->Y[atomIdx1];
        z1=OutState->Z[atomIdx1];
        float x2,y2,z2;
        x2=OutState->X[atomIdx2];
        y2=OutState->Y[atomIdx2];
        z2=OutState->Z[atomIdx2];
        Vector r12(x2-x1,y2-y1,z2-z1);
        float r=r12.length();
        float deltaR= r-r0;
        Vector f12=-2.0*r12*deltaR*k/r;
        OutState->FX[atomIdx1]+=f12.x;
        OutState->FY[atomIdx1]+=f12.y;
        OutState->FZ[atomIdx1]+=f12.z;
        OutState->FX[atomIdx2]-=f12.x;
        OutState->FY[atomIdx2]-=f12.y;
        OutState->FZ[atomIdx2]-=f12.z;
        float e;
        e=k*deltaR*deltaR;
        energy+=e;
#if DEBUG_BOND
        printf("Bond%d (%d,%d) k=%f r0=%f r=%f dr=%f e=%f\n",ii,atomIdx1,
                atomIdx2,k,r0,r,deltaR,e);
#endif
    }
    printf("EBond=%f\n",energy);
    return 0;
}
int computeAngles(const Molecule *InMolecule, State *OutState){
    const float PI=4.0*atan(1.0);
    const float RAD2DEG=180.0/PI;
    float energy=0.0;
    const int numOfAngles=InMolecule->NumOfAngles;
    for(int ii=0; ii<numOfAngles; ++ii){
        const int atomIdx1=InMolecule->ANGLE[ii].a;
        const int atomIdx2=InMolecule->ANGLE[ii].b;
        const int atomIdx3=InMolecule->ANGLE[ii].c;
        float theta0=InMolecule->ANGLE[ii].theta;
        float k=InMolecule->ANGLE[ii].k;
        float x1,y1,z1;
        x1=OutState->X[atomIdx1];
        y1=OutState->Y[atomIdx1];
        z1=OutState->Z[atomIdx1];
        float x2,y2,z2;
        x2=OutState->X[atomIdx2];
        y2=OutState->Y[atomIdx2];
        z2=OutState->Z[atomIdx2];
        float x3,y3,z3;
        x3=OutState->X[atomIdx3];
        y3=OutState->Y[atomIdx3];
        z3=OutState->Z[atomIdx3];
        Vector r12(x1-x2,y1-y2,z1-z2);
        Vector r32(x3-x2,y3-y2,z3-z2);
        float d12=r12.length();
        float d32=r32.length();
        float cosTheta=(r12*r32)/(d12*d32);
        if(cosTheta>1.0){
            cosTheta=1.0;
        }else if(cosTheta<-1.0){
            cosTheta=-1.0;
        }
        float sinTheta=sqrt(1.0-cosTheta*cosTheta);
        float theta=acos(cosTheta);
        float deltaTheta=theta-theta0/RAD2DEG;
        //use k instead of 0.5k to sue CHARMM convention
        float e=k*deltaTheta*deltaTheta;
        energy+=e;
        float d12inv=1.0/d12;
        float d32inv=1.0/d32;
        //use a factor of 2 to sue CHARMM convention
        float factor= -2.0*k*deltaTheta/sinTheta;
        float c1=factor*d12inv;
        float c2=factor*d32inv;
        /*
         * f12=-{k*deltaTheta}/{sin(theta)*|d12|}*
         * {r12/|r12|*cos(theta)-r32/|r32|}
         * f32=-{k*deltaTheta}/{sin(theta)*|d32|}*
         * {r32/|r32|*cos(theta)-r12/|r12|}
         */
        //it seems that mindy is wrong
        Vector f12=-1.0*c1*(r12*(d12inv*cosTheta)-r32*d32inv);
        Vector f1=f12;
        Vector f32=-1.0*c2*(r32*(d32inv*cosTheta)-r12*d12inv);
        Vector f3=f32;
        Vector f2=-f12-f32;
        OutState->FX[atomIdx1]+=f1.x;
        OutState->FY[atomIdx1]+=f1.y;
        OutState->FZ[atomIdx1]+=f1.z;

        OutState->FX[atomIdx2]+=f2.x;
        OutState->FY[atomIdx2]+=f2.y;
        OutState->FZ[atomIdx2]+=f2.z;

        OutState->FX[atomIdx3]+=f3.x;
        OutState->FY[atomIdx3]+=f3.y;
        OutState->FZ[atomIdx3]+=f3.z;
#if DEBUG_BOND
        printf("Angle%d (%d,%d,%d) k=%f theta0=%f theta=%f dTheta=%f e=%f\n",
                ii,atomIdx1,atomIdx2,atomIdx3,k,theta0,theta*RAD2DEG,deltaTheta*RAD2DEG,e);
        //printf("factor=%f k=%f deltaTheta=%f sinTheta=%f d12inv=%f\n",factor,k,deltaTheta,sinTheta,d12inv);
#endif
    }
    printf("EAngle=%f\n",energy);
    printf("RAD2DEG=%f\n",RAD2DEG);
    return 0;
}
int computeDihedrals(const Molecule *InMolecule, State *OutState){
    const float PI=4.0*atan(1.0);
    const float RAD2DEG=180.0/PI;
    float energy=0.0;
    const int numOfDihedrals=InMolecule->NumOfDihedrals;
    for(int ii=0; ii<numOfDihedrals; ++ii){
        const int atomIdx1=InMolecule->DIHEDRAL[ii].a;
        const int atomIdx2=InMolecule->DIHEDRAL[ii].b;
        const int atomIdx3=InMolecule->DIHEDRAL[ii].c;
        const int atomIdx4=InMolecule->DIHEDRAL[ii].d;
        const float theta0=InMolecule->DIHEDRAL[ii].theta;
        const float k=InMolecule->DIHEDRAL[ii].k;
        const float mult=InMolecule->DIHEDRAL[ii].multiplicity;
        float e=0.0;
        float theta;
        float deltaTheta;

        float x1,y1,z1;
        x1=OutState->X[atomIdx1];
        y1=OutState->Y[atomIdx1];
        z1=OutState->Z[atomIdx1];
        float x2,y2,z2;
        x2=OutState->X[atomIdx2];
        y2=OutState->Y[atomIdx2];
        z2=OutState->Z[atomIdx2];
        float x3,y3,z3;
        x3=OutState->X[atomIdx3];
        y3=OutState->Y[atomIdx3];
        z3=OutState->Z[atomIdx3];
        float x4,y4,z4;
        x4=OutState->X[atomIdx4];
        y4=OutState->Y[atomIdx4];
        z4=OutState->Z[atomIdx4];
        Vector r12(x1-x2,y1-y2,z1-z2);
        Vector r23(x2-x3,y2-y3,z2-z3);
        Vector r34(x3-x4,y3-y4,z3-z4);

        Vector dcosdA;
        Vector dcosdB;
        Vector dsindC;
        Vector dsindB;
        Vector f1, f2, f3; 

        Vector A, B, C;
        A.cross(r12, r23);
        B.cross(r23, r34);
        C.cross(r23, A); 

        float rA = A.length();
        float rB = B.length(); 
        float rC = C.length();

        float cos_phi = (A*B)/(rA*rB);
        float sin_phi = (C*B)/(rC*rB);
        // Normalize B
        rB = 1.0/rB;
        B *= rB; 

        float phi = -atan2(sin_phi, cos_phi);
        theta=phi*RAD2DEG;
        //mindy has bugs
        deltaTheta=theta-theta0;

        if (fabs(sin_phi) > 0.1) {
            // Normalize A
            rA = 1.0/rA;
            A *= rA; 
            dcosdA = rA*(cos_phi*A-B);
            dcosdB = rB*(cos_phi*B-A);
        }   
        else {
            // Normalize C
            rC = 1.0/rC;
            C *= rC;
            dsindC = rC*(sin_phi*C-B);
            dsindB = rB*(sin_phi*B-C);
        }
        for (int j=0; j<1; j++) {
            float n = mult;
            float delta = (-1.0*theta0)/RAD2DEG;
            float K, K1;
            if (n) {
                K = k * (1.0+cos(n*phi + delta));
                K1 = -n*k*sin(n*phi + delta);
            }
            else {
                float diff = phi-delta;
                if (diff < -PI) diff += 2.0*PI;
                else if (diff > PI) diff -= 2.0*PI;
                K = k*diff*diff;
                K1 = 2.0*k*diff;
            }
            e+=K;
            energy += K;
            // forces
            if (fabs(sin_phi) > 0.1) {
                //mindy bug
                //K1 = K1/sin_phi; 
                K1 = -1.0*K1/sin_phi;
                f1.x += K1*(r23.y*dcosdA.z - r23.z*dcosdA.y);
                f1.y += K1*(r23.z*dcosdA.x - r23.x*dcosdA.z);
                f1.z += K1*(r23.x*dcosdA.y - r23.y*dcosdA.x);

                f3.x += K1*(r23.z*dcosdB.y - r23.y*dcosdB.z);
                f3.y += K1*(r23.x*dcosdB.z - r23.z*dcosdB.x);
                f3.z += K1*(r23.y*dcosdB.x - r23.x*dcosdB.y);

                f2.x += K1*(r12.z*dcosdA.y - r12.y*dcosdA.z
                        + r34.y*dcosdB.z - r34.z*dcosdB.y);
                f2.y += K1*(r12.x*dcosdA.z - r12.z*dcosdA.x
                        + r34.z*dcosdB.x - r34.x*dcosdB.z);
                f2.z += K1*(r12.y*dcosdA.x - r12.x*dcosdA.y
                        + r34.x*dcosdB.y - r34.y*dcosdB.x);
            }
            else {
                //  This angle is closer to 0 or 180 than it is to
                //  90, so use the cos version to avoid 1/sin terms
                //mindy bug
                //K1 = -K1/cos_phi;
                K1=K1/cos_phi;

                f1.x += K1*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
                        - r23.x*r23.y*dsindC.y
                        - r23.x*r23.z*dsindC.z);
                f1.y += K1*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
                        - r23.y*r23.z*dsindC.z
                        - r23.y*r23.x*dsindC.x);
                f1.z += K1*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
                        - r23.z*r23.x*dsindC.x
                        - r23.z*r23.y*dsindC.y);

                f3 += cross(K1,dsindB,r23);

                f2.x += K1*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
                        +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
                        +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
                        +dsindB.z*r34.y - dsindB.y*r34.z);
                f2.y += K1*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
                        +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
                        +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
                        +dsindB.x*r34.z - dsindB.z*r34.x);
                f2.z += K1*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
                        +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
                        +(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y
                        +dsindB.y*r34.x - dsindB.x*r34.y);
            }
        }    // end loop over multiplicity
        //f[dihedral->atom1] += f1;
        OutState->FX[atomIdx1]+=f1.x;
        OutState->FY[atomIdx1]+=f1.y;
        OutState->FZ[atomIdx1]+=f1.z;
        //f[dihedral->atom2] += f2-f1;
        OutState->FX[atomIdx2]+=f2.x-f1.x;
        OutState->FY[atomIdx2]+=f2.y-f1.y;
        OutState->FZ[atomIdx2]+=f2.z-f1.z;
        //f[dihedral->atom3] += f3-f2;
        OutState->FX[atomIdx3]+=f3.x-f2.x;
        OutState->FY[atomIdx3]+=f3.y-f2.y;
        OutState->FZ[atomIdx3]+=f3.z-f2.z;
        //f[dihedral->atom4] += -f3;
        OutState->FX[atomIdx4]+=-f3.x;
        OutState->FY[atomIdx4]+=-f3.y;
        OutState->FZ[atomIdx4]+=-f3.z;
#if DEBUG_BOND
        printf("Dihedral%d (%d,%d,%d,%d) k=%f theta0=%f theta=%f dTheta=%f e=%f\n",
                ii,atomIdx1,atomIdx2,atomIdx3,atomIdx4,k,theta0,theta,deltaTheta,e);
#endif
    }
    printf("EDiheral=%f\n",energy);
    printf("RAD2DEG=%f\n",RAD2DEG);
}
#endif
