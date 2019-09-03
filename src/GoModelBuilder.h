#ifndef __GO_MODEL_BUILDER_H__
#define __GO_MODEL_BUILDER_H__
#include "DataType.h"
#include "DataStructure.h"
#include "IO.h"
#include <map>
#include <tuple>
#include <set>
#include "StringUtility.h"

//#define testBuildGoModel
//#define testCalcBonds
//#define testCalcAngles
//#define testCalcDihedrals
//#define testGetNumOfChains
//#define testCalcHBonds
//#define testCalcSidechainContacts
//#define testCalcNativeContacts
//#define testBond
//#define testAngle
//#define testDihedral
//#define testNative
//#define testVDW
//#define testScaleNative
//#define testMassCharge

namespace KBGoParams{
    //residue charge map
    std::map<const std::string,const float> ResidueCharge = {
        {"ALA",0},{"ARG",1},{"ASN",0},{"ASP",-1},
        {"CYS",0},{"GLU",-1},{"GLN",0},{"GLY",0},
        {"HIS",0.5},{"ILE",0},{"LEU",0},{"LYS",1},
        {"MET",0},{"PHE",0},{"PRO",0}, {"SER",0},
        {"THR",0},{"TRP",0},{"TYR",0},{"VAL",0},
        //patches
        {"HSD",0.5},{"HSE",0.5},{"HSP",0.5},{"HSC",0.5},
        {"ZN",-1}
    };
    //residue mass map
    std::map<const std::string,const float> ResidueMass = {
        {"ALA",71.0788}, {"ARG",156.1875},{"ASN",114.1038},{"ASP",115.0886},
        {"CYS",103.1388},{"GLU",129.1155},{"GLN",128.1307},{"GLY",57.0519},
        {"HIS",137.1411},{"ILE",113.1594},{"LEU",113.1594},{"LYS",128.1741},
        {"MET",131.1926},{"PHE",147.1766},{"PRO",97.1167}, {"SER",87.0782},
        {"THR",101.1051},{"TRP",186.2132},{"TYR",163.1760},{"VAL",99.1326},
        {"HSD",137.1411},{"HSE",137.1411},{"HSP",138.1491},{"HSC",138},{"ZN",65.38}
    };
    //residue 3 letter to 1 letter name map
    std::map<const std::string,const char> AminoThreeToOne = {
        {"GLY",'G'},{"PRO",'P'},{"ALA",'A'},{"VAL",'V'},
        {"LEU",'L'},{"ILE",'I'},{"MET",'M'},{"CYS",'C'},
        {"PHE",'F'},{"TYR",'Y'},{"TRP",'W'},{"HIS",'H'},
        {"HSD",'H'},{"HSE",'H'},{"LYS",'K'},{"ARG",'R'},
        {"GLN",'Q'},{"ASN",'N'},{"GLU",'E'},{"ASP",'D'},
        {"SER",'S'},{"THR",'T'},{"HSC",'H'},{"HSP",'H'}
    };
    //residue 1 letter to 3 letter name map
    //H is HSD only
    std::map<const char,const std::string> AminoOneToThree = {
        {'G',"GLY"},{'P',"PRO"},{'A',"ALA"},{'V',"VAL"},
        {'L',"LEU"},{'I',"ILE"},{'M',"MET"},{'C',"CYS"},
        {'F',"PHE"},{'Y',"TYR"},{'W',"TRP"},
        {'H',"HSD"},{'K',"LYS"},{'R',"ARG"},
        {'Q',"GLN"},{'N',"ASN"},{'E',"GLU"},{'D',"ASP"},
        {'S',"SER"},{'T',"THR"}
    };

    /*hydrogen bond related*/
    static const float KSfac=27.888; //0.42 * 0.2 * 332
    static const float HBondEcut=-0.5; //energy cutoff
    static const float HBondDistCut=5.22; //distance cutoff
    static const float HFactor=1.0;
    /*side chain contact related*/
    static const float SidechainDistCut=4.5; //sidechain-sidechain distance cutoff
    static const float SidechainCADistCut=14.0; //CA-CA distance cutoff
    /*MJ potential*/
    static const float MJFactor=0.5; //
    static const float MJAve=-3.1757;
    /*Scaling Native Nonbond*/
    static const float TfScaling=0.0054;
    /*Bonded*/
    static const float BondFactor=200;
    static const float AngleFactor=40;
    static const float DihedralFactor=0.4;
    /*Folding Temperature*/
    static const float Tf=350;
    /*nonbond nonnative nonbond*/
    static const float VDWEpsilonFactor=0.00007;
    static const float VDWSigmaFactor=0.561231;

#include "data/GoPseudoDihedral.dat"
    //defined PseudoDihedral[1600]
#include "data/GoMJPotential96.dat"
    //defined MJPotential[400]
};

struct nativeScaler;
struct nativeScaler{
    char first;
    char second;
    float factor;
};

//data structure for hbond and sidechain native pair
struct contactPair;
struct contactPair{
    int first;
    int second;
    float E;
    float r;
    bool operator==(contactPair const &rhs) const {
        return (this->first==rhs.first)&&(this->second==rhs.second);
    }
};
//sorting function for contactPair
bool sortContactPair(const contactPair &left, const contactPair &right){
    if(left.first<right.first){
        return true;
    }
    else if (left.first>right.first){
        return false;
    }
    else{
        if(left.second<right.second){
            return true;
        }else{
            return false;
        }
    }
}

int calcNativeContacts(std::vector<contactPair> &hbondList, std::vector<contactPair> &sideChainContact, std::vector<contactPair> &nativeContact, Molecule *OutMol);
int calcSidechainContacts(const Molecule *InMol, Molecule *OutMol, std::vector<contactPair> &sideChainContact);
int isblankString(const char *str);
int allocateGoMolecule(const int NumOfAtoms, Molecule *OutMol);
int getNumOfChains(const Molecule *InMol, int *NumOfChains, int **NumOfAtomsPerChain);
int calcBonds(const Molecule *InMol, Molecule *OutMol,const int NumOfChains, 
        const int *NumOfAtomsPerChain);
int calcAngles(const Molecule *InMol, Molecule *OutMol,const int NumOfChains,
        const int *NumOfAtomsPerChain);
int calcDihedrals(const Molecule *InMol, Molecule *OutMol,const int NumOfChains,
        const int *NumOfAtomsPerChain);
int calcHBonds(const Molecule *InMol, Molecule *OutMol, const int NumOfChains, const int *NumOfAtomsPerChain, std::vector<contactPair> &hbondList);
int _isSameChain(const int NumOfChains, const int *NumOfAtomsPerChain, int atomi, int atomj);
int calcZincFinger(const Molecule *InMol, Molecule *OutMol);

int buildGoModel(const Molecule *InMol, Molecule *OutMol){
    int numOfCGAtoms=0;
    char alphaC[3]="CA";
    char zinc[3]="ZN"; //zinc finger
    int NumOfChains;
    int *NumOfAtomsPerChain;
    if(mystrcmp(" CA ","CA")==0){
        printf("mystrcmp is OK\n");
    }
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        if(mystrcmp(InMol->AtomType[ii],alphaC)==0){
#ifdef testBuildGoModel
            printf("Atom%d '%s' is CA\n",ii,InMol->AtomType[ii]);
#endif
            numOfCGAtoms++;
        }
        else if (mystrcmp(InMol->AtomType[ii],zinc)==0){
#ifdef testBuildGoModel
            printf("Atom%d '%s' is ZN\n",ii,InMol->AtomType[ii]);
#endif
            numOfCGAtoms++;
        }
    }
    allocateGoMolecule(numOfCGAtoms,OutMol);
    //extract CG atoms
    int indexCG=0;
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        if(mystrcmp(InMol->AtomType[ii],alphaC)==0 ||
                mystrcmp(InMol->AtomType[ii],zinc)==0){
            char *tmp;
#ifdef testBuildGoModel
            printf("Atom%d '%s' is extracted\n",ii,InMol->AtomType[ii]);
#endif
            OutMol->X[indexCG]=InMol->X[ii];

            OutMol->Y[indexCG]=InMol->Y[ii];

            OutMol->Z[indexCG]=InMol->Z[ii];

            OutMol->AtomNum[indexCG]=InMol->AtomNum[ii];

            mytrim(InMol->AtomType[ii],&tmp);
            OutMol->AtomType[indexCG]=tmp;

            mytrim(InMol->ResType[ii],&tmp);
            OutMol->ResType[indexCG]=tmp;

            OutMol->ChainId[indexCG]=InMol->ChainId[ii];

            mytrim(InMol->ChainCHARMMId[ii],&tmp);
            OutMol->ChainCHARMMId[indexCG]=tmp;

            OutMol->ResNum[indexCG]=InMol->ResNum[ii];

            //get residue mass
            char *resType = OutMol->ResType[indexCG];
            char *trimResType;
            mytrim(resType,&trimResType);
            std::string resTypeString = trimResType;
            OutMol->Mass[indexCG]=KBGoParams::ResidueMass[trimResType];
            OutMol->Charge[indexCG]=KBGoParams::ResidueCharge[trimResType];

            indexCG++;
        }
    }
    //get number of chains and number of atoms per chain
    //int NumOfChains int NumOfAtomsPerChain[int]
    getNumOfChains(OutMol, &NumOfChains, &NumOfAtomsPerChain);
#ifdef testGetNumOfChains
    printf("Testing getNumOfChains()...\n");
    for(int ii=0; ii<NumOfChains; ++ii){
        printf("Chain%d has %d atoms\n",ii,NumOfAtomsPerChain[ii]);
    }
#endif
    //calculate zinc finger
    calcZincFinger(InMol, OutMol);
    //calculate bond
    calcBonds(InMol,OutMol,NumOfChains,NumOfAtomsPerChain);
    //calculate angle
    calcAngles(InMol,OutMol,NumOfChains,NumOfAtomsPerChain);
    //calculate dihedrals
    calcDihedrals(InMol, OutMol,NumOfChains,NumOfAtomsPerChain);
    //hydrogen bond list
    std::vector<contactPair> hbondList;
    //sidechain contact list
    std::vector<contactPair> sideChainContact;
    //native contact list
    std::vector<contactPair> nativeContact;
    std::vector<contactPair> nativeContact2;
    //calculate hydrogen bond
    calcHBonds(InMol,OutMol,NumOfChains,NumOfAtomsPerChain,hbondList);
#ifdef testCalcHBonds
    printf("Testing calcHBond()...\n");
    for(contactPair tmp : hbondList){
        printf("HBond (%d,%d) r=%f\n",tmp.first,tmp.second,tmp.r);
    }
#endif
    //calculate sidechain contacts
    calcSidechainContacts(InMol, OutMol,sideChainContact);
#ifdef testCalcSidechainContacts
    printf("Testing calcSidechainContacts()...\n");
    for(contactPair pair : sideChainContact){
        printf("Sidechain contact(%d,%d) r=%f E=%f\n",pair.first,pair.second,
                pair.r,pair.E);
    }
#endif
    //calculate native contacts
    calcNativeContacts(hbondList,sideChainContact,nativeContact2,OutMol);
    for(contactPair pair : nativeContact2){
        if(mystrcmp(OutMol->AtomType[pair.first],"ZN")!=0 &&
                mystrcmp(OutMol->AtomType[pair.second],"ZN")!=0){
#ifdef testCalcNativeContacts
        printf("Native contact(%d,%d) r= %f E=%f\n",pair.first,pair.second,
                pair.r,pair.E);
#endif
            nativeContact.push_back(pair);
        }
    }
    printf("New native contact number=%d\n",nativeContact.size());
    //calcualte total native nonbond energy
    float nonbondEtotal=0;
    for(contactPair pair : nativeContact){
        nonbondEtotal+=pair.E;
    }

    //calculate native nonbond energy scaling factor
    int NumOfZNs=0;
    for(int ii=0; ii<OutMol->NumOfAtoms; ++ii){
        if(mystrcmp(OutMol->AtomType[ii],"ZN")==0){
            NumOfZNs++;
        }
    }
    printf("NumOfZNs=%d\n",NumOfZNs);
    float nonbondScale=-1.0*(KBGoParams::TfScaling)*(KBGoParams::Tf)*
        (OutMol->NumOfResidues-NumOfZNs)/nonbondEtotal;

    printf("NumOfResidues=%d\n",OutMol->NumOfResidues);
    printf("ScalingFactor=%f\n",nonbondScale);
    OutMol->nativeScalingFactor=nonbondScale;
    //calcualte new total native nonbond energy
    nonbondEtotal=nonbondEtotal*nonbondScale;
    //calculate native nonbond energy per residue
    float EperResidue=nonbondEtotal/(OutMol->NumOfResidues-NumOfZNs);
    //calculate scaling factor for bond 
    float bondScale = -1.0*KBGoParams::BondFactor*EperResidue;
    //calculate scaling factor for angle
    float angleScale = -1.0*KBGoParams::AngleFactor*EperResidue;
    //calculate scaling factor for dihedral
    float dihedralScale = -1.0*KBGoParams::DihedralFactor*EperResidue;
    //calculate scaling factor for nonnative nonbond Epsilon
    float epsilon = KBGoParams::VDWEpsilonFactor*EperResidue;
    //scale the bond k
    for(int ii=0;ii<OutMol->NumOfBonds;++ii){
        OutMol->BOND[ii].k = bondScale;
#ifdef testBond
        printf("Bond(%d,%d) k=%f r=%f\n",OutMol->BOND[ii].i,OutMol->BOND[ii].j,
                OutMol->BOND[ii].k,OutMol->BOND[ii].r);
#endif
    }
    //scale the angle k
    for(int ii=0;ii<OutMol->NumOfAngles;++ii){
        OutMol->ANGLE[ii].k=angleScale;
#ifdef testAngle
        printf("Angle(%d,%d,%d) k=%f theta=%f\n",OutMol->ANGLE[ii].a,OutMol->ANGLE[ii].b,
                OutMol->ANGLE[ii].c,OutMol->ANGLE[ii].k,OutMol->ANGLE[ii].theta);
#endif
    }
    //scale the dihedral k
    for(int ii=0;ii<OutMol->NumOfDihedrals;++ii){
        OutMol->DIHEDRAL[ii].k*=dihedralScale;
#ifdef testDihedral
        printf("Dihedral(%d,%d,%d,%d) k=%f theta=%f\n",OutMol->DIHEDRAL[ii].a,
                OutMol->DIHEDRAL[ii].b,OutMol->DIHEDRAL[ii].c,
                OutMol->DIHEDRAL[ii].d, OutMol->DIHEDRAL[ii].k,
                OutMol->DIHEDRAL[ii].theta);
#endif
    }
    //calculate number of native contacts
    OutMol->NumOfNatives=nativeContact.size();
#ifdef testNative
    printf("Num of total natives=%d\n",OutMol->NumOfNatives);
#endif
    //allocate native contact memory
    OutMol->NATIVE=(Native*)malloc(nativeContact.size()*sizeof(Native));
    //map nativeContact vector to OutMol->NATIVE array
    for(int ii=0; ii<OutMol->NumOfNatives; ++ii){
        OutMol->NATIVE[ii].i = nativeContact[ii].first;
        OutMol->NATIVE[ii].j = nativeContact[ii].second;
        //scale native contact energy
        OutMol->NATIVE[ii].k = nativeContact[ii].E * nonbondScale;
        //recalculate r to avoid potential bugs
        int atom1=OutMol->NATIVE[ii].i;
        int atom2=OutMol->NATIVE[ii].j;
        float newDistance=sqrt((OutMol->X[atom1]-OutMol->X[atom2])*
                (OutMol->X[atom1]-OutMol->X[atom2]) + 
                (OutMol->Y[atom1]-OutMol->Y[atom2])*
                (OutMol->Y[atom1]-OutMol->Y[atom2]) +
                (OutMol->Z[atom1]-OutMol->Z[atom2])*
                (OutMol->Z[atom1]-OutMol->Z[atom2]) );
        OutMol->NATIVE[ii].r = newDistance;
#ifdef testNative
        printf("Native(%d,%d) k=%f r=%f\n",OutMol->NATIVE[ii].i,OutMol->NATIVE[ii].j,
                OutMol->NATIVE[ii].k,OutMol->NATIVE[ii].r);
#endif
    }
    //calculate nonnative nonbond parameters VDWR and EPS
    for(int ii=0; ii<OutMol->NumOfAtoms; ++ii){
        //nearest residue id to residue ii
        int nearest;
        //nearest residue distance to residue ii
        float dist=9999;
        for(int jj=0; jj<OutMol->NumOfAtoms; ++jj){
            //nearest residue need to be both more than 3 residues away
            //and not forming native contact with residue ii
            if(abs(ii-jj)>3){
                bool foundInNative=false;
                //test if (ii,jj) is native contact
                for(contactPair pair : nativeContact){
                    if((pair.first==ii && pair.second==jj) ||
                            (pair.first==jj && pair.second==ii)){
                        foundInNative=true;
                        break;
                    }
                }
                //if not, calculate newdist
                if(!foundInNative){
                    float newdist=sqrt((OutMol->X[ii]-OutMol->X[jj])*
                            (OutMol->X[ii]-OutMol->X[jj])+
                            (OutMol->Y[ii]-OutMol->Y[jj])*
                            (OutMol->Y[ii]-OutMol->Y[jj])+
                            (OutMol->Z[ii]-OutMol->Z[jj])*
                            (OutMol->Z[ii]-OutMol->Z[jj]));
                    //if newdist is smaller
                    //set it to be the nearest
                    if(newdist<dist){
                        dist=newdist;
                        nearest=jj;
                    }
                }
            }
        }
        //scale the van der waals radii
        OutMol->VDWR[ii]=dist*KBGoParams::VDWSigmaFactor;
        //set the van der waals well depth
        OutMol->EPS[ii]=epsilon;
#ifdef testVDW
        printf("Atom%d %d eps=%f vdwr=%f\n",ii,nearest,epsilon,OutMol->VDWR[ii]);
#endif
#ifdef testMassCharge
        printf("Atom%d resid=%c mass=%f charge=%f\n",ii,OutMol->ChainId[ii],OutMol->Mass[ii],OutMol->Charge[ii]);
#endif
    }
    printf("Done generating Go Model\n");
    return 0;
}
int calcNativeContacts(std::vector<contactPair> &hbondList, std::vector<contactPair> &sideChainContact, std::vector<contactPair> &nativeContact, Molecule *OutMol){    
    printf("Combining hydrogen bond and sidechain contact to form native contact list\n");
    int NumOfChains;
    int *NumOfAtomsPerChain;
    getNumOfChains(OutMol, &NumOfChains, &NumOfAtomsPerChain);

    std::vector<contactPair> singleHbond;
    std::vector<contactPair> doubleHbond;
    std::vector<contactPair> orientSingleHbond;
    std::vector<contactPair> orientDoubleHbond;
    //hbondList.push_back({111,111});
    //hbondList.push_back({111,111});
    //hbondList.push_back({110,110});
    //find single and double hydrogen bond
    while(!hbondList.empty()){
        contactPair pair=hbondList.front();
        hbondList.erase(hbondList.begin());
        bool isdouble=false;
        bool exist=false;
        for(contactPair tmpPair : hbondList){
            for(contactPair tmpPair2 : doubleHbond){
                if(pair==tmpPair2){
                    exist=true;
                }
            }
            if(tmpPair==pair){
                isdouble=true;
            }
        }
        if(isdouble==true){
            doubleHbond.push_back(pair);
        }
        else if (!exist){
            singleHbond.push_back(pair);
        }
    }
    for(contactPair pair : sideChainContact){
        nativeContact.push_back(pair);
    }
    for(contactPair pair : singleHbond){
        bool foundInSide=false;
        for(contactPair sidepair : sideChainContact){
            if(pair.first==sidepair.first && pair.second==sidepair.second){
                foundInSide=true;
                break;
            }
        }
        if(!foundInSide){
            contactPair newpair;
            newpair.first=pair.first;
            newpair.second=pair.second;
            newpair.r=pair.r;
            newpair.E=-1.0*KBGoParams::HFactor;
            nativeContact.push_back(newpair);
        }
        else{
            orientSingleHbond.push_back(pair);
        }
    }
    for(contactPair pair : doubleHbond){
        bool foundInSide=false;
        for(contactPair sidepair : sideChainContact){
            if(pair.first==sidepair.first && pair.second==sidepair.second){
                foundInSide=true;
                break;
            }
        }
        if(!foundInSide){
            contactPair newpair;
            newpair.first=pair.first;
            newpair.second=pair.second;
            newpair.r=pair.r;
            newpair.E=-1.0*KBGoParams::HFactor;
            nativeContact.push_back(newpair);
            orientSingleHbond.push_back(pair);
        }
        else{
            orientDoubleHbond.push_back(pair);
        }
    }
    //orientation hbond
    for(contactPair pair : orientSingleHbond){
        int atom1,atom2;
        atom1=pair.first;
        atom2=pair.second;
        int numorie=0;
        int numRes=OutMol->NumOfAtoms;
        if(atom1>0){numorie++;};
        if(atom2<numRes-1){numorie++;};
        if(abs(atom1-atom2)>=4){numorie+=2;};
        float Eorie=-1.0*KBGoParams::HFactor*1/numorie;
        if(atom1>0){
            contactPair oriepair;
            oriepair.first=atom1-1;
            oriepair.second=atom2;
            //calculate new r
            oriepair.r=sqrt((OutMol->X[oriepair.first]-OutMol->X[oriepair.second])*
                    (OutMol->X[oriepair.first]-OutMol->X[oriepair.second]) +
                    (OutMol->Y[oriepair.first]-OutMol->Y[oriepair.second])*
                    (OutMol->Y[oriepair.first]-OutMol->Y[oriepair.second]) +
                    (OutMol->Z[oriepair.first]-OutMol->Z[oriepair.second])*
                    (OutMol->Z[oriepair.first]-OutMol->Z[oriepair.second]));
            oriepair.E=Eorie;
            bool foundInNative=false;
            int pos=0;
            for(contactPair nativepair : nativeContact){
                if(nativepair==oriepair){
                    foundInNative=true;
                    break;
                }
                pos++;
            }
            if(foundInNative){
                nativeContact[pos].E+=oriepair.E;
            }
            else{
                nativeContact.push_back(oriepair);
            }
        }
        if(atom2<numRes-1){
            contactPair oriepair;
            oriepair.first=atom1;
            oriepair.second=atom2+1;
            oriepair.r=sqrt((OutMol->X[oriepair.first]-OutMol->X[oriepair.second])*
                    (OutMol->X[oriepair.first]-OutMol->X[oriepair.second]) +
                    (OutMol->Y[oriepair.first]-OutMol->Y[oriepair.second])*
                    (OutMol->Y[oriepair.first]-OutMol->Y[oriepair.second]) +
                    (OutMol->Z[oriepair.first]-OutMol->Z[oriepair.second])*
                    (OutMol->Z[oriepair.first]-OutMol->Z[oriepair.second]));
            oriepair.E=Eorie;
            bool foundInNative=false;
            int pos=0;
            for(contactPair nativepair : nativeContact){
                if(nativepair==oriepair){
                    foundInNative=true;
                    break;
                }
                pos++;
            }
            if(foundInNative){
                nativeContact[pos].E+=oriepair.E;
            }
            else{
                nativeContact.push_back(oriepair);
            }

        }
        if(abs(atom1-atom2)>=4){
            contactPair oriepair1;
            oriepair1.first=atom1;
            oriepair1.second=atom2-1;
            oriepair1.r=sqrt((OutMol->X[oriepair1.first]-OutMol->X[oriepair1.second])*
                    (OutMol->X[oriepair1.first]-OutMol->X[oriepair1.second]) +
                    (OutMol->Y[oriepair1.first]-OutMol->Y[oriepair1.second])*
                    (OutMol->Y[oriepair1.first]-OutMol->Y[oriepair1.second]) +
                    (OutMol->Z[oriepair1.first]-OutMol->Z[oriepair1.second])*
                    (OutMol->Z[oriepair1.first]-OutMol->Z[oriepair1.second]));
            oriepair1.E=Eorie;
            bool foundInNative=false;
            int pos=0;
            for(contactPair nativepair : nativeContact){
                if(nativepair==oriepair1){
                    foundInNative=true;
                    break;
                }
                pos++;
            }
            if(foundInNative){
                nativeContact[pos].E+=oriepair1.E;
            }
            else{
                nativeContact.push_back(oriepair1);
            }

            contactPair oriepair2;
            oriepair2.first=atom1+1;
            oriepair2.second=atom2;
            oriepair2.r=sqrt((OutMol->X[oriepair2.first]-OutMol->X[oriepair2.second])*
                    (OutMol->X[oriepair2.first]-OutMol->X[oriepair2.second]) +
                    (OutMol->Y[oriepair2.first]-OutMol->Y[oriepair2.second])*
                    (OutMol->Y[oriepair2.first]-OutMol->Y[oriepair2.second]) +
                    (OutMol->Z[oriepair2.first]-OutMol->Z[oriepair2.second])*
                    (OutMol->Z[oriepair2.first]-OutMol->Z[oriepair2.second]));
            oriepair2.E=Eorie;
            foundInNative=false;
            pos=0;
            for(contactPair nativepair : nativeContact){
                if(nativepair==oriepair2){
                    foundInNative=true;
                    break;
                }
                pos++;
            }
            if(foundInNative){
                nativeContact[pos].E+=oriepair2.E;
            }
            else{
                nativeContact.push_back(oriepair2);
            }

        }
    }
    for(contactPair pair : orientDoubleHbond){
        int atom1,atom2;
        atom1=pair.first;
        atom2=pair.second;
        int numorie=0;
        int numRes=OutMol->NumOfAtoms;
        if(atom1>0){numorie++;};
        if(atom2<numRes-1){numorie++;};
        if(abs(atom1-atom2)>=4){numorie+=2;}; 
        //2 for double hbond
        float Eorie=-1.0*KBGoParams::HFactor*2/numorie;
        if(atom1>0){
            contactPair oriepair;
            oriepair.first=atom1-1;
            oriepair.second=atom2;
            oriepair.r=sqrt((OutMol->X[oriepair.first]-OutMol->X[oriepair.second])*
                    (OutMol->X[oriepair.first]-OutMol->X[oriepair.second]) +
                    (OutMol->Y[oriepair.first]-OutMol->Y[oriepair.second])*
                    (OutMol->Y[oriepair.first]-OutMol->Y[oriepair.second]) +
                    (OutMol->Z[oriepair.first]-OutMol->Z[oriepair.second])*
                    (OutMol->Z[oriepair.first]-OutMol->Z[oriepair.second]));
            oriepair.E=Eorie;
            bool foundInNative=false;
            int pos=0;
            for(contactPair nativepair : nativeContact){
                if(nativepair==oriepair){
                    foundInNative=true;
                    break;
                }
                pos++;
            }
            if(foundInNative){
                nativeContact[pos].E+=oriepair.E;
            }
            else{
                nativeContact.push_back(oriepair);
            }
        }
        if(atom2<numRes-1){
            contactPair oriepair;
            oriepair.first=atom1;
            oriepair.second=atom2+1;
            oriepair.r=sqrt((OutMol->X[oriepair.first]-OutMol->X[oriepair.second])*
                    (OutMol->X[oriepair.first]-OutMol->X[oriepair.second]) +
                    (OutMol->Y[oriepair.first]-OutMol->Y[oriepair.second])*
                    (OutMol->Y[oriepair.first]-OutMol->Y[oriepair.second]) +
                    (OutMol->Z[oriepair.first]-OutMol->Z[oriepair.second])*
                    (OutMol->Z[oriepair.first]-OutMol->Z[oriepair.second]));
            oriepair.E=Eorie;
            bool foundInNative=false;
            int pos=0;
            for(contactPair nativepair : nativeContact){
                if(nativepair==oriepair){
                    foundInNative=true;
                    break;
                }
                pos++;
            }
            if(foundInNative){
                nativeContact[pos].E+=oriepair.E;
            }
            else{
                nativeContact.push_back(oriepair);
            }

        }
        if(abs(atom1-atom2)>=4){
            contactPair oriepair1;
            oriepair1.first=atom1;
            oriepair1.second=atom2-1;
            oriepair1.r=sqrt((OutMol->X[oriepair1.first]-OutMol->X[oriepair1.second])*
                    (OutMol->X[oriepair1.first]-OutMol->X[oriepair1.second]) +
                    (OutMol->Y[oriepair1.first]-OutMol->Y[oriepair1.second])*
                    (OutMol->Y[oriepair1.first]-OutMol->Y[oriepair1.second]) +
                    (OutMol->Z[oriepair1.first]-OutMol->Z[oriepair1.second])*
                    (OutMol->Z[oriepair1.first]-OutMol->Z[oriepair1.second]));
            oriepair1.E=Eorie;
            bool foundInNative=false;
            int pos=0;
            for(contactPair nativepair : nativeContact){
                if(nativepair==oriepair1){
                    foundInNative=true;
                    break;
                }
                pos++;
            }
            if(foundInNative){
                nativeContact[pos].E+=oriepair1.E;
            }
            else{
                nativeContact.push_back(oriepair1);
            }

            contactPair oriepair2;
            oriepair2.first=atom1+1;
            oriepair2.second=atom2;
            oriepair2.r=sqrt((OutMol->X[oriepair2.first]-OutMol->X[oriepair2.second])*
                    (OutMol->X[oriepair2.first]-OutMol->X[oriepair2.second]) +
                    (OutMol->Y[oriepair2.first]-OutMol->Y[oriepair2.second])*
                    (OutMol->Y[oriepair2.first]-OutMol->Y[oriepair2.second]) +
                    (OutMol->Z[oriepair2.first]-OutMol->Z[oriepair2.second])*
                    (OutMol->Z[oriepair2.first]-OutMol->Z[oriepair2.second]));
            oriepair2.E=Eorie;
            foundInNative=false;
            pos=0;
            for(contactPair nativepair : nativeContact){
                if(nativepair==oriepair2){
                    foundInNative=true;
                    break;
                }
                pos++;
            }
            if(foundInNative){
                nativeContact[pos].E+=oriepair2.E;
            }
            else{
                nativeContact.push_back(oriepair2);
            }

        }
    }
    std::sort(nativeContact.begin(),nativeContact.end(),sortContactPair);
    return 0;
}
int calcSidechainContacts(const Molecule *InMol, Molecule *OutMol, std::vector<contactPair> &sideChainContact){
    int atomIdx=0;
    int hbondIdx=0;
    int numOfCA=OutMol->NumOfAtoms;
    printf("Calculating sidechain contacts: NumOfAtoms=%d\n",InMol->NumOfAtoms);
    int numOfAtoms=InMol->NumOfAtoms;
    std::vector<std::vector<int>> residue;
    std::vector<int> sidechainAtoms;
    std::vector<int> alphaCarbon;
    int resNum=InMol->ResNum[0];
    for(int ii=0; ii<numOfAtoms; ++ii){
        if(resNum!=InMol->ResNum[ii]){
            //new residue
            if(sidechainAtoms.size()>0){
                residue.push_back(sidechainAtoms);
                sidechainAtoms.clear();
            }
            resNum=InMol->ResNum[ii];
        }
        char *atomType=InMol->AtomType[ii];
        if(mystrcmp(atomType,"CA")==0){
            alphaCarbon.push_back(ii);
        }else if (mystrcmp(atomType,"H")!=0&&mystrcmp(atomType,"N")!=0
                &&mystrcmp(atomType,"C")!=0&&mystrcmp(atomType,"O")!=0){
            char *trimAtomType;
            mytrim(atomType,&trimAtomType);
            if(trimAtomType[0]!='H'){
                sidechainAtoms.push_back(ii);
            }
        }
    }
    //the last sidechainAtoms array was not added in the previous loop
    if(!sidechainAtoms.empty()){
        residue.push_back(sidechainAtoms);
        sidechainAtoms.clear();
    }
    for(std::vector<int> residue1 : residue){
        for(std::vector<int> residue2 : residue){
            int resNum1,resNum2;
            resNum1=InMol->ResNum[residue1[0]];
            resNum2=InMol->ResNum[residue2[0]];
            if(resNum2>=resNum1+3){
                float minDist=KBGoParams::SidechainDistCut+1;
                float caDist;
                for(int atom1 : residue1){
                    for(int atom2 : residue2){
                        float dist=sqrt(
                                (InMol->X[atom1]-InMol->X[atom2])*
                                (InMol->X[atom1]-InMol->X[atom2])+
                                (InMol->Y[atom1]-InMol->Y[atom2])*
                                (InMol->Y[atom1]-InMol->Y[atom2])+
                                (InMol->Z[atom1]-InMol->Z[atom2])*
                                (InMol->Z[atom1]-InMol->Z[atom2]) );
                        if(dist<minDist){
                            minDist=dist;
                        }
                    }
                }
                int CA1,CA2;
                CA1=alphaCarbon[resNum1];
                CA2=alphaCarbon[resNum2];
                caDist=sqrt(
                        (InMol->X[CA1]-InMol->X[CA2])*
                        (InMol->X[CA1]-InMol->X[CA2])+
                        (InMol->Y[CA1]-InMol->Y[CA2])*
                        (InMol->Y[CA1]-InMol->Y[CA2])+
                        (InMol->Z[CA1]-InMol->Z[CA2])*
                        (InMol->Z[CA1]-InMol->Z[CA2]) );
                if(minDist<KBGoParams::SidechainDistCut && caDist<KBGoParams::SidechainCADistCut){
                    char *resType1=InMol->ResType[CA1];
                    char *resType2=InMol->ResType[CA2];
                    char *trimResType1,*trimResType2;
                    mytrim(resType1,&trimResType1);
                    mytrim(resType2,&trimResType2);
                    //convert three letters to one letter and then convert back
                    //to get a unique three letters code of amino acid
                    std::string resTypeString1=KBGoParams::AminoOneToThree[
                        KBGoParams::AminoThreeToOne[trimResType1]];
                    std::string resTypeString2=KBGoParams::AminoOneToThree[
                        KBGoParams::AminoThreeToOne[trimResType2]];
                    //printf("residue%d %s and residue%d %s are native\n",
                    //resNum1,resTypeString1.c_str(),resNum2,resTypeString2.c_str());
                    contactPair pair;
                    pair.first=resNum1;
                    pair.second=resNum2;
                    pair.r=caDist;
                    for(int ii=0;ii<KBGoParams::NumOfMJPairs; ++ii){
                        if((resTypeString1==KBGoParams::MJPotential[ii].residue1 &&
                                    resTypeString2==KBGoParams::MJPotential[ii].residue2)||
                                (resTypeString1==KBGoParams::MJPotential[ii].residue2 &&
                                 resTypeString2==KBGoParams::MJPotential[ii].residue1))
                        {
                            /*printf("%s %s E=%f\n",resTypeString1.c_str(),resTypeString2.c_str(),
                                    KBGoParams::MJPotential[ii].value);*/
                            pair.E=-1.0*KBGoParams::MJPotential[ii].value*
                                KBGoParams::MJFactor/KBGoParams::MJAve;
                            //pair.E=KBGoParams::MJPotential[ii].value;
                            break;
                        }
                    }
                    sideChainContact.push_back(pair);
                }
            }
        }
    }
    return 0;
}
int calcHBonds(const Molecule *InMol, Molecule *OutMol, const int NumOfChains, const int *NumOfAtomsPerChain, std::vector<contactPair> &hbondList){
    int atomIdx=0;
    int hbondIdx=0;
    int numOfCA=OutMol->NumOfAtoms;
    //multiplicity = 4
    printf("Calculating hbonds: NumOfChains=%d NumOfAtoms=%d NumOfCA=%d\n",NumOfChains,OutMol->NumOfAtoms,numOfCA);
    int numOfAtoms=InMol->NumOfAtoms;
    //residue stores the CNOH array 
    std::vector<int*> residue;
    //CNOH stores the C,N,O,H,CA index and the residue index
    int CNOH[6]={-1,-1,-1,-1,-1,-1};
    int resNum=InMol->ResNum[0];
    for(int ii=0; ii<numOfAtoms;++ii){
        if(resNum!=InMol->ResNum[ii]){
            CNOH[5]=resNum;
            resNum=InMol->ResNum[ii];
            //printf("residue%d\n",resNum);
            bool isCorrectResidue=true;
            if(isCorrectResidue){
                //printf("residue%d was added\n",resNum);
                int *tmp;
                tmp=(int*)malloc(5*sizeof(int));
                tmp[0]=CNOH[0];
                tmp[1]=CNOH[1];
                tmp[2]=CNOH[2];
                tmp[3]=CNOH[3];
                tmp[4]=CNOH[4];
                tmp[5]=CNOH[5];
                residue.push_back(tmp);
            }
            CNOH[0]=-1;
            CNOH[1]=-1;
            CNOH[2]=-1;
            CNOH[3]=-1;
            CNOH[4]=-1;
            CNOH[5]=-1;
        }
        char *atomType=InMol->AtomType[ii];
        if(mystrcmp(atomType,"C")==0){
            CNOH[0]=ii;
            //printf("atom%d C\n",ii);
        }
        else if(mystrcmp(atomType,"N")==0){
            CNOH[1]=ii;
            //printf("atom%d N\n",ii);
        }
        else if(mystrcmp(atomType,"O")==0){
            CNOH[2]=ii;
            //printf("atom%d O\n",ii);
        }
        else if(mystrcmp(atomType,"H")==0){
            CNOH[3]=ii;
            //printf("atom%d H\n",ii);
        }
        else if(mystrcmp(atomType,"HN")==0){
            CNOH[3]=ii;
            printf("Warning: found CHARMM HN equivalent to H atom%d\n",ii);
        }
        else if(mystrcmp(atomType,"CA")==0){
            CNOH[4]=ii;
            //printf("atom%d H\n",ii);
        }
        
    }

    //handle the last residue
    CNOH[5]=resNum;
    //printf("residue%d\n",resNum);
    bool isCorrectResidue=true;
    if(isCorrectResidue){
        //printf("residue%d was added as last\n",resNum);
        int *tmp;
        tmp=(int*)malloc(5*sizeof(int));
        tmp[0]=CNOH[0];
        tmp[1]=CNOH[1];
        tmp[2]=CNOH[2];
        tmp[3]=CNOH[3];
        tmp[4]=CNOH[4];
        tmp[5]=CNOH[5];
        residue.push_back(tmp);
    }

    for(int *residue1 : residue){
        int resNum1=InMol->ResNum[residue1[0]];
        //printf("residue%d in chain%c\n",resNum1,InMol->ChainId[residue1[0]]);
        //printf("*********************\n");
    }

    for(int *residue1 : residue){
        //printf("residue%d C=%d N=%d O=%d H=%d\n",tmp[4],
        //        tmp[0],tmp[1],tmp[2],tmp[3]);
        for(int *residue2 : residue){
            int resNum1=residue1[5];
            int resNum2=residue2[5];
            char chainId1=InMol->ChainId[residue1[0]];
            char chainId2=InMol->ChainId[residue2[0]];
            if(((chainId1==chainId2)&&((resNum2-resNum1)>=3))||
                    (chainId1!=chainId2&&(resNum2>resNum1))){
                //printf("res1=%d res2=%d\n",residue1[4],residue2[4]);
                int C1,N1,O1,H1;
                int C2,N2,O2,H2;
                C1=residue1[0];N1=residue1[1];
                O1=residue1[2];H1=residue1[3];
                C2=residue2[0];N2=residue2[1];
                O2=residue2[2];H2=residue2[3];
                float rON,rCH,rOH,rCN;
                //N1H1 to C2O2
                if(N1!=-1 && H1!=-1 && C2!=-1 && O2!=-1){
                    rON=sqrt((InMol->X[N1]-InMol->X[O2])*(InMol->X[N1]-InMol->X[O2])+
                            (InMol->Y[N1]-InMol->Y[O2])*(InMol->Y[N1]-InMol->Y[O2])+
                            (InMol->Z[N1]-InMol->Z[O2])*(InMol->Z[N1]-InMol->Z[O2]));
                    rCH=sqrt((InMol->X[H1]-InMol->X[C2])*(InMol->X[H1]-InMol->X[C2])+
                            (InMol->Y[H1]-InMol->Y[C2])*(InMol->Y[H1]-InMol->Y[C2])+
                            (InMol->Z[H1]-InMol->Z[C2])*(InMol->Z[H1]-InMol->Z[C2]));
                    rOH=sqrt((InMol->X[H1]-InMol->X[O2])*(InMol->X[H1]-InMol->X[O2])+
                            (InMol->Y[H1]-InMol->Y[O2])*(InMol->Y[H1]-InMol->Y[O2])+
                            (InMol->Z[H1]-InMol->Z[O2])*(InMol->Z[H1]-InMol->Z[O2]));
                    rCN=sqrt((InMol->X[N1]-InMol->X[C2])*(InMol->X[N1]-InMol->X[C2])+
                            (InMol->Y[N1]-InMol->Y[C2])*(InMol->Y[N1]-InMol->Y[C2])+
                            (InMol->Z[N1]-InMol->Z[C2])*(InMol->Z[N1]-InMol->Z[C2]));
                    if(rON<KBGoParams::HBondDistCut){    
                        float HbondE = 1.0/rON + 1.0/rCH - 1.0/rOH - 1.0/rCN;
                        HbondE=HbondE*KBGoParams::KSfac;
                        if(HbondE<KBGoParams::HBondEcut){
                            //add hbond
                            //printf("%d %d hbond (ON=%f,CH=%f,OH=%f,CN=%f)\n",
                            //        resNum1,resNum2,rON,rCH,rOH,rCN);
                            int CA1=residue1[4];
                            int CA2=residue2[4];
                            float rCACA=sqrt((InMol->X[CA1]-InMol->X[CA2])*(InMol->X[CA1]-InMol->X[CA2])+
                                    (InMol->Y[CA1]-InMol->Y[CA2])*(InMol->Y[CA1]-InMol->Y[CA2])+
                                    (InMol->Z[CA1]-InMol->Z[CA2])*(InMol->Z[CA1]-InMol->Z[CA2]));
                            contactPair pair;
                            pair.first=resNum1;
                            pair.second=resNum2;
                            pair.E=0;
                            pair.r=rCACA;
                            hbondList.push_back(pair);
                        } //endif
                    } //endif
                } //endif

                if(N2!=-1 && H2!=-1 && C1!=-1 && O1!=-1){
                    //N2H2 to C1O1
                    rON=sqrt((InMol->X[N2]-InMol->X[O1])*(InMol->X[N2]-InMol->X[O1])+
                            (InMol->Y[N2]-InMol->Y[O1])*(InMol->Y[N2]-InMol->Y[O1])+
                            (InMol->Z[N2]-InMol->Z[O1])*(InMol->Z[N2]-InMol->Z[O1]));
                    rCH=sqrt((InMol->X[H2]-InMol->X[C1])*(InMol->X[H2]-InMol->X[C1])+
                            (InMol->Y[H2]-InMol->Y[C1])*(InMol->Y[H2]-InMol->Y[C1])+
                            (InMol->Z[H2]-InMol->Z[C1])*(InMol->Z[H2]-InMol->Z[C1]));
                    rOH=sqrt((InMol->X[H2]-InMol->X[O1])*(InMol->X[H2]-InMol->X[O1])+
                            (InMol->Y[H2]-InMol->Y[O1])*(InMol->Y[H2]-InMol->Y[O1])+
                            (InMol->Z[H2]-InMol->Z[O1])*(InMol->Z[H2]-InMol->Z[O1]));
                    rCN=sqrt((InMol->X[N2]-InMol->X[C1])*(InMol->X[N2]-InMol->X[C1])+
                            (InMol->Y[N2]-InMol->Y[C1])*(InMol->Y[N2]-InMol->Y[C1])+
                            (InMol->Z[N2]-InMol->Z[C1])*(InMol->Z[N2]-InMol->Z[C1]));
                    if(rON<KBGoParams::HBondDistCut){    
                        float HbondE = 1.0/rON + 1.0/rCH - 1.0/rOH - 1.0/rCN;
                        HbondE=HbondE*KBGoParams::KSfac;
                        if(HbondE<KBGoParams::HBondEcut){
                            //printf("%d %d hbond (ON=%f,CH=%f,OH=%f,CN=%f)\n",
                            //        resNum1,resNum2,rON,rCH,rOH,rCN);
                            //add hbond
                            int CA1=residue1[4];
                            int CA2=residue2[4];
                            float rCACA=sqrt((InMol->X[CA1]-InMol->X[CA2])*(InMol->X[CA1]-InMol->X[CA2])+
                                    (InMol->Y[CA1]-InMol->Y[CA2])*(InMol->Y[CA1]-InMol->Y[CA2])+
                                    (InMol->Z[CA1]-InMol->Z[CA2])*(InMol->Z[CA1]-InMol->Z[CA2]));
                            contactPair pair;
                            pair.first=resNum1;
                            pair.second=resNum2;
                            pair.E=0;
                            pair.r=rCACA;
                            hbondList.push_back(pair);
                        } //endif
                    } //endif
                } //endif
            } //endif
        }//end for
    } //end for
    /*
    for(int ii=0; ii<numOfAtoms-3; ++ii){
        resNum1=InMol->ResNum[ii];
        for(int jj=ii+1; jj<numOfAtoms;++jj){
            resNum2=InMol->ResNum[jj];
            //within the same chain
            //j>=i+3
            if(InMol->ChainId[ii]==InMol->ChainId[jj]){
            }
            //in different chains not i,j distance constraint
            else{
                //printf("(%d,%d) (%c,%c) is not in the same chain\n",ii,jj,
                //        InMol->ChainId[ii],InMol->ChainId[jj]);
            }
        }
    }*/
    return 0;
}
int _isSameChain(const int NumOfChains, const int *NumOfAtomsPerChain, int atomi, int atomj){
    int chaina=0,chainb=0;
    for(int ii=0; ii<NumOfChains; ++ii){
        if(atomi<NumOfAtomsPerChain[ii]&&(ii!=0&&(atomi>NumOfAtomsPerChain[ii-1]))){
            chaina=ii;
        }
        if(atomj<NumOfAtomsPerChain[ii]&&(ii!=0&&(atomj>NumOfAtomsPerChain[ii-1]))){
            chainb=ii;
        }
    }
    printf("Atom%d is in chain %d, atom%d is in chain %d\n",atomi,chaina,atomj,chainb);
    if(chaina==chainb){
        return 1;
    }else{
        return 0;
    }
}
int calcDihedrals(const Molecule *InMol, Molecule *OutMol,const int NumOfChains,
        const int *NumOfAtomsPerChain){
    const float PI = 4*atan(1);
    const float RAD2DEG = 180.0/PI;
    int atomIdx=0;
    int diheIdx=0;
    int numOfCG=OutMol->NumOfAtoms;
    int numOfCA=0;
    for(int ii=0; ii<numOfCG; ++ii){
        if(mystrcmp(OutMol->AtomType[ii],"CA")==0){
            numOfCA++;
        }
    }
    //multiplicity = 4
    const int numOfDihedral=4*(numOfCA-3*NumOfChains);
    OutMol->NumOfDihedrals=numOfDihedral;
    OutMol->DIHEDRAL=(Dihedral*)malloc(numOfDihedral*sizeof(Dihedral));
    printf("Calculating dihedral: NumOfChains=%d NumOfAtoms=%d NumOfDihedrals=%d\n",NumOfChains,OutMol->NumOfAtoms,numOfDihedral);
    for(int ii=0; ii<NumOfChains; ++ii){
#ifdef testCalcDihedrals
        printf("******CHAIN%d=%d******\n",ii,NumOfAtomsPerChain[ii]);
#endif
        for(int jj=0; jj<NumOfAtomsPerChain[ii]-3; ++jj){
            //printf("atom id = %d\n",atomIdx);
            int atoma=atomIdx;
            int atomb=atomIdx+1;
            int atomc=atomIdx+2;
            int atomd=atomIdx+3;
            if(mystrcmp(OutMol->AtomType[atoma],"CA")==0 &&
                    mystrcmp(OutMol->AtomType[atomb],"CA")==0 &&
                    mystrcmp(OutMol->AtomType[atomc],"CA")==0 &&
                    mystrcmp(OutMol->AtomType[atomd],"CA")==0){
                char *res1,*res2;
                mytrim(OutMol->ResType[atomb],&res1);
                mytrim(OutMol->ResType[atomc],&res2);
                char res1code,res2code;
                res1code=KBGoParams::AminoThreeToOne[res1];
                res2code=KBGoParams::AminoThreeToOne[res2];
                //printf("res1='%c' res2='%c'\n",res1code,res2code);
                int mult[4]={1,2,3,4};
                float k[4];
                float phi[4];
                for(int kk=0; kk<KBGoParams::NumOfPseudoDihedral; ++kk){
                    char tmpRes1=KBGoParams::PseudoDihedral[kk].residue1;
                    char tmpRes2=KBGoParams::PseudoDihedral[kk].residue2;
                    int multValue=KBGoParams::PseudoDihedral[kk].mult;
                    float kValue=KBGoParams::PseudoDihedral[kk].k;
                    float phiValue=KBGoParams::PseudoDihedral[kk].phi;
                    if((res1code==tmpRes1) && (res2code==tmpRes2)){
                        //printf("%c %c %d %f %f\n",tmpRes1,tmpRes2,multValue,kValue,phiValue);
                        if(multValue==mult[0]){
                            k[0]=kValue;
                            phi[0]=phiValue;
                        }
                        else if(multValue==mult[1]){
                            k[1]=kValue;
                            phi[1]=phiValue;
                        }
                        else if(multValue==mult[2]){
                            k[2]=kValue;
                            phi[2]=phiValue;
                        }
                        else if(multValue==mult[3]){
                            k[3]=kValue;
                            phi[3]=phiValue;
                        }
                    }
                }
                for(int n=0; n<4; ++n){
                    Dihedral dihe;
                    dihe.a=atoma;
                    dihe.b=atomb;
                    dihe.c=atomc;
                    dihe.d=atomd;
                    dihe.theta=phi[n];
                    dihe.k=k[n];
                    dihe.multiplicity=mult[n];
                    OutMol->DIHEDRAL[diheIdx*4+n]=dihe;
#ifdef testCalcDihedrals
                    printf("dihe%d (%d,%d,%d,%d) mult=%d theta=%f k=%f\n",(diheIdx*4+n),
                            atoma,atomb,atomc,atomd,dihe.multiplicity,dihe.theta,dihe.k);
#endif
                }
                diheIdx++;
            } //endif
            atomIdx++;
        } //end for
        atomIdx+=3;
    }
}
int calcAngles(const Molecule *InMol, Molecule *OutMol,const int NumOfChains,
        const int *NumOfAtomsPerChain){
    const float PI = 4*atan(1);
    const float RAD2DEG = 180.0/PI;
    printf("Calculating angle: NumOfChains=%d NumOfAtoms=%d\n",NumOfChains,OutMol->NumOfAtoms);
    int atomIdx=0;
    int angleIdx=0;
    int numOfCG=OutMol->NumOfAtoms;
    int numOfCA=0;
    for(int ii=0; ii<numOfCG; ++ii){
        if(mystrcmp(OutMol->AtomType[ii],"CA")==0){
            numOfCA++;
        }
    }
    OutMol->ANGLE=(Angle*)malloc((numOfCA-2*NumOfChains)*sizeof(Angle));
    OutMol->NumOfAngles=numOfCA-2*NumOfChains;
    for(int ii=0; ii<NumOfChains; ++ii){
#ifdef testCalcAngles
        printf("******CHAIN%d=%d******\n",ii,NumOfAtomsPerChain[ii]);
#endif
        for(int jj=0; jj<NumOfAtomsPerChain[ii]-2; ++jj){
            //printf("atom id = %d\n",atomIdx);
            if(mystrcmp(OutMol->AtomType[atomIdx],"CA")==0 &&
                    mystrcmp(OutMol->AtomType[atomIdx+1],"CA")==0 &&
                    mystrcmp(OutMol->AtomType[atomIdx+2],"CA")==0){
                float theta=0;
                float x1,y1,z1,x2,y2,z2,x3,y3,z3;
                x1=OutMol->X[atomIdx];
                y1=OutMol->Y[atomIdx];
                z1=OutMol->Z[atomIdx];
                x2=OutMol->X[atomIdx+1];
                y2=OutMol->Y[atomIdx+1];
                z2=OutMol->Z[atomIdx+1];
                x3=OutMol->X[atomIdx+2];
                y3=OutMol->Y[atomIdx+2];
                z3=OutMol->Z[atomIdx+2];
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
                theta=acos(cosTheta);
                Angle angle;
                angle.a=atomIdx;
                angle.b=atomIdx+1;
                angle.c=atomIdx+2;
                angle.theta=theta*RAD2DEG;
                OutMol->ANGLE[angleIdx]=angle;
#ifdef testCalcAngles
                printf("angle%d (%d,%d,%d) theta=%f k=%f\n",angleIdx,atomIdx,atomIdx+1,atomIdx+2,angle.theta,angle.k);
#endif
                angleIdx++;
            } //endif
            atomIdx++;
        } //endfor
        atomIdx+=2;
    }
}

int calcBonds(const Molecule *InMol, Molecule *OutMol,const int NumOfChains, 
        const int *NumOfAtomsPerChain){
    printf("Calculating bond: NumOfChains=%d NumOfAtoms=%d\n",NumOfChains,OutMol->NumOfAtoms);
    int atomIdx=0;
    int bondIdx=0;
    int numOfCA=0;
    for(int ii=0; ii< OutMol->NumOfAtoms; ++ii){
        if(mystrcmp(OutMol->AtomType[ii],"CA")==0){
            numOfCA++;
        }
    }
    OutMol->NumOfBonds=numOfCA-NumOfChains;
    OutMol->BOND=(Bond*)malloc((numOfCA-NumOfChains)*sizeof(Bond));
    for(int ii=0; ii<NumOfChains; ++ii){
#ifdef testCalcBonds
        printf("******CHAIN%d=%d******\n",ii,NumOfAtomsPerChain[ii]);
#endif
        for(int jj=0; jj<NumOfAtomsPerChain[ii]-1; ++jj){
            //printf("atom id = %d\n",atomIdx);
            if(mystrcmp(OutMol->AtomType[atomIdx],"CA")==0 &&
                    mystrcmp(OutMol->AtomType[atomIdx+1],"CA")==0
              ){
                float bondLength=0;
                bondLength=sqrt((OutMol->X[atomIdx]-OutMol->X[atomIdx+1])*(OutMol->X[atomIdx]-OutMol->X[atomIdx+1])+
                        (OutMol->Y[atomIdx]-OutMol->Y[atomIdx+1])*(OutMol->Y[atomIdx]-OutMol->Y[atomIdx+1])+
                        (OutMol->Z[atomIdx]-OutMol->Z[atomIdx+1])*(OutMol->Z[atomIdx]-OutMol->Z[atomIdx+1]));
                Bond bond;
                bond.i=atomIdx;
                bond.j=atomIdx+1;
                bond.r=bondLength;
                OutMol->BOND[bondIdx]=bond;
#ifdef testCalcBonds
                printf("bond%d (%d,%d) r=%f k=%f\n",bondIdx,atomIdx,atomIdx+1,bond.r,bond.k);
#endif
                bondIdx++;
            }
            atomIdx++;
        }
        atomIdx++;
    }
    return 0;
}
//get the number of chains
int getNumOfChains(const Molecule *InMol, int *NumOfChains, int **NumOfAtomsPerChain){
    std::vector<int> natoms;
    bool usesChainCHARMMId=false;
    *NumOfChains=1;
    if(isblank(InMol->ChainId[0])&&isblankString(InMol->ChainCHARMMId[0])){
        *NumOfChains=1;
    }
    else if(!isblank(InMol->ChainId[0])){
        usesChainCHARMMId=false;
    }
    else if(!isblankString(InMol->ChainCHARMMId[0])){
        usesChainCHARMMId=true;
        printf("ChainCHARMMId is empty\n");
    }
    else{
        usesChainCHARMMId=false;
    }
    if(!usesChainCHARMMId){
        char currentId=InMol->ChainId[0];
        int numOfAtoms=1;
        //printf("%c\n",currentId);
        for(int ii=1; ii < InMol->NumOfAtoms; ++ii){
            if(currentId!=InMol->ChainId[ii]){
                natoms.push_back(numOfAtoms);
                numOfAtoms=0;
                currentId=InMol->ChainId[ii];
                (*NumOfChains)++;
            }
            numOfAtoms++;
        }
        natoms.push_back(numOfAtoms);
    }else{
        char *currentId=InMol->ChainCHARMMId[0];
        int numOfAtoms=1;
        for(int ii=1; ii < InMol->NumOfAtoms; ++ii){
            numOfAtoms++;
            if(strcmp(currentId,InMol->ChainCHARMMId[ii])!=0){
                natoms.push_back(numOfAtoms);
                numOfAtoms=0;
                currentId=InMol->ChainCHARMMId[ii];
                *(NumOfChains)++;
            }
        }
        natoms.push_back(numOfAtoms);
    }
    (*NumOfAtomsPerChain)=(int*)malloc((*NumOfChains)*sizeof(int));
    for(int ii=0; ii<(*NumOfChains); ++ii){
        (*NumOfAtomsPerChain)[ii]=natoms[ii];
    }
    return 0;
}
//check if a c string is blank
int isblankString(const char *str){
    for(int ii=0; ii<strlen(str); ++ii){
        if(isblank(str[ii])){
            return 0;
        }
    }
    return 1;
}
//allocate memory for Molecule
int allocateGoMolecule(const int NumOfAtoms, Molecule *OutMol){
    printf("Allocating new moledule with %d atoms\n",NumOfAtoms);
    OutMol->NumOfAtoms=NumOfAtoms;
    //each residue only contains alpha carbon 
    OutMol->NumOfResidues=NumOfAtoms;
    OutMol->X=(float*)malloc(NumOfAtoms*sizeof(float));
    OutMol->Y=(float*)malloc(NumOfAtoms*sizeof(float));
    OutMol->Z=(float*)malloc(NumOfAtoms*sizeof(float));
    OutMol->tmpX=(float*)malloc(NumOfAtoms*sizeof(float));
    OutMol->tmpY=(float*)malloc(NumOfAtoms*sizeof(float));
    OutMol->tmpZ=(float*)malloc(NumOfAtoms*sizeof(float));
    OutMol->VDWR=(float*)malloc(NumOfAtoms*sizeof(float));
    OutMol->EPS=(float*)malloc(NumOfAtoms*sizeof(float));
    OutMol->Charge=(float*)malloc(NumOfAtoms*sizeof(float));
    OutMol->Mass=(float*)malloc(NumOfAtoms*sizeof(float));
    OutMol->AtomNum=(int*)malloc(NumOfAtoms*sizeof(int));
    OutMol->AtomType=(char**)malloc(NumOfAtoms*sizeof(char**));
    OutMol->ResType=(char**)malloc(NumOfAtoms*sizeof(char**));
    OutMol->ChainId=(char*)malloc(NumOfAtoms*sizeof(char*));
    OutMol->ChainCHARMMId=(char**)malloc(NumOfAtoms*sizeof(char**));
    OutMol->ResNum=(int*)malloc(NumOfAtoms*sizeof(int));
    return 0;
}
//copy a new molecule form a given molecule
int copyGoMolecule(const Molecule *InMol, Molecule *OutMol){
    allocateGoMolecule(InMol->NumOfAtoms,OutMol);
    OutMol->NumOfResidues=InMol->NumOfResidues;
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        OutMol->X[ii]=InMol->X[ii];
        OutMol->Y[ii]=InMol->Y[ii];
        OutMol->Z[ii]=InMol->Z[ii];
        OutMol->tmpX[ii]=InMol->X[ii];
        OutMol->tmpY[ii]=InMol->Y[ii];
        OutMol->tmpZ[ii]=InMol->Z[ii];
        OutMol->EPS[ii]=InMol->EPS[ii];
        OutMol->VDWR[ii]=InMol->VDWR[ii];
        OutMol->Charge[ii]=InMol->Charge[ii];
        OutMol->Mass[ii]=InMol->Mass[ii];
        OutMol->AtomNum[ii]=InMol->AtomNum[ii];
        OutMol->AtomType[ii]=InMol->AtomType[ii];
        OutMol->ResType[ii]=InMol->ResType[ii];
        OutMol->ChainId[ii]=InMol->ChainId[ii];
        OutMol->ChainCHARMMId[ii]=InMol->ChainCHARMMId[ii];
        OutMol->ResNum[ii]=InMol->ResNum[ii];
    }
    OutMol->NumOfBonds=InMol->NumOfBonds;
    OutMol->BOND=(Bond*)malloc(OutMol->NumOfBonds*sizeof(Bond));
    for(int ii=0; ii<InMol->NumOfBonds; ++ii){
        OutMol->BOND[ii]=InMol->BOND[ii];
    }

    OutMol->NumOfAngles=InMol->NumOfAngles;
    OutMol->ANGLE=(Angle*)malloc(OutMol->NumOfAngles*sizeof(Angle));
    for(int ii=0; ii<InMol->NumOfAngles; ++ii){
        OutMol->ANGLE[ii]=InMol->ANGLE[ii];
    }

    OutMol->NumOfDihedrals=InMol->NumOfDihedrals;
    OutMol->DIHEDRAL=(Dihedral*)malloc(OutMol->NumOfDihedrals*sizeof(Dihedral));
    for(int ii=0; ii<InMol->NumOfDihedrals; ++ii){
        OutMol->DIHEDRAL[ii]=InMol->DIHEDRAL[ii];
    }

    OutMol->NumOfNatives=InMol->NumOfNatives;
    OutMol->NATIVE=(Native*)malloc(OutMol->NumOfNatives*sizeof(Native));
    for(int ii=0; ii<InMol->NumOfNatives; ++ii){
        OutMol->NATIVE[ii]=InMol->NATIVE[ii];
    }

    OutMol->NumOfZincBonds=InMol->NumOfZincBonds;
    OutMol->ZINCBOND=(Bond*)malloc(OutMol->NumOfZincBonds*sizeof(Bond));
    for(int ii=0; ii<InMol->NumOfZincBonds; ++ii){
        OutMol->ZINCBOND[ii]=InMol->ZINCBOND[ii];
    }
    OutMol->NumOfZincAngles=InMol->NumOfZincAngles;
    OutMol->ZINCANGLE=(Angle*)malloc(OutMol->NumOfZincAngles*sizeof(Angle));
    for(int ii=0; ii<InMol->NumOfZincAngles; ++ii){
        OutMol->ZINCBOND[ii]=InMol->ZINCBOND[ii];
    }
    return 0;
}
//merge multiple pdbfiles
//this is specially designed for taz1-hif-cited2 system
int mergeGoMolecule(const std::vector<Molecule> &InMols, 
        const std::vector<std::vector<int>> &ChainIds,Molecule *OutMol){
    int i=0;
    std::vector<int*> AtomsPerChain;
    int numOfAtoms=0;
    for(Molecule mol: InMols){
        int numOfChains;
        int *atomsPerChain;
        getNumOfChains(&mol,&numOfChains,&atomsPerChain);
        AtomsPerChain.push_back(atomsPerChain);
        //printMolecule(&mol);
        printf("mol%d nchains=%d\n",i,numOfChains);
        for(int ii=0;ii<numOfChains;++ii){
            printf("chain%d natoms=%d\n",ii,atomsPerChain[ii]);
        }
        for(int chainId: ChainIds[i]){
            numOfAtoms+=atomsPerChain[chainId];
        }
        i++;
    }
    allocateGoMolecule(numOfAtoms,OutMol);
    std::vector<std::vector<int>> atomIds;
    i=0;
    int atomId=0;
    char newChainId='A';
    for(Molecule mol: InMols){
        std::vector<int> chainId=ChainIds[i];
        for(int chain : chainId){
            char chainIdLetter=chain+'A';
            for(int ii=0; ii<mol.NumOfAtoms; ++ii){
                if(chainIdLetter==mol.ChainId[ii]){
                    OutMol->X[atomId]=mol.X[ii];
                    OutMol->Y[atomId]=mol.Y[ii];
                    OutMol->Z[atomId]=mol.Z[ii];
                    OutMol->VDWR[atomId]=mol.VDWR[ii];
                    OutMol->EPS[atomId]=mol.EPS[ii];
                    OutMol->Charge[atomId]=mol.Charge[ii];
                    OutMol->Mass[atomId]=mol.Mass[ii];
                    OutMol->AtomNum[atomId]=atomId;
                    OutMol->AtomType[atomId]=mol.AtomType[ii];
                    OutMol->ResType[atomId]=mol.ResType[ii];
                    OutMol->ChainId[atomId]=newChainId;
                    OutMol->ChainCHARMMId[atomId]=mol.ChainCHARMMId[ii];
                    OutMol->ResNum[atomId]=atomId;
                    atomId++;
                }
            }
            newChainId++;
        }
        i++;
    }
    Molecule taz1=InMols[0]; 
    Molecule taz1_hif=InMols[1];
    Molecule cited2_taz1=InMols[2];
    int ntaz1=AtomsPerChain[0][0];
    int nhif=AtomsPerChain[1][1];
    int ncited2=AtomsPerChain[2][0];
    printf("ntaz1=%d,nhif=%d,ncited2=%d\n",ntaz1,nhif,ncited2);
    std::vector<Bond> bond;
    std::vector<Angle> angle;
    std::vector<Dihedral> dihedral;
    std::vector<Native> native;
    std::vector<Bond> zincBond;
    std::vector<Angle> zincAngle;
    int offset;
    for(int ii=0; ii<taz1.NumOfZincBonds; ++ii){
        //no need for renumbering
        zincBond.push_back(taz1.ZINCBOND[ii]);
    }
    for(int ii=0; ii<taz1.NumOfZincAngles; ++ii){
        //no need for renumbering
        zincAngle.push_back(taz1.ZINCANGLE[ii]);
    }
    for(int ii=0; ii<taz1.NumOfBonds; ++ii){
        //no need for renumbering
        bond.push_back(taz1.BOND[ii]);
    }
    for(int ii=0; ii<taz1.NumOfAngles; ++ii){
        //no need for renumbering
        angle.push_back(taz1.ANGLE[ii]);
    }
    for(int ii=0; ii<taz1.NumOfDihedrals; ++ii){
        //no need for renumbering
        dihedral.push_back(taz1.DIHEDRAL[ii]);
    }
    for(int ii=0; ii<taz1.NumOfNatives; ++ii){
        //no need for renumbering
        native.push_back(taz1.NATIVE[ii]);
    }

    offset=0;
    for(int ii=0; ii<taz1_hif.NumOfBonds; ++ii){
        if(taz1_hif.ChainId[taz1_hif.BOND[ii].i]=='B'
                && taz1_hif.ChainId[taz1_hif.BOND[ii].j]=='B'){
            Bond newbond;
            newbond.i=taz1_hif.BOND[ii].i+offset;
            newbond.j=taz1_hif.BOND[ii].j+offset;
            newbond.k=taz1_hif.BOND[ii].k;
            newbond.r=taz1_hif.BOND[ii].r;
            bond.push_back(newbond);
        }
    }
    for(int ii=0; ii<taz1_hif.NumOfAngles; ++ii){
        if(taz1_hif.ChainId[taz1_hif.ANGLE[ii].a]=='B'
                && taz1_hif.ChainId[taz1_hif.ANGLE[ii].b]=='B'){
            Angle newangle;
            newangle.a=taz1_hif.ANGLE[ii].a+offset;
            newangle.b=taz1_hif.ANGLE[ii].b+offset;
            newangle.c=taz1_hif.ANGLE[ii].c+offset;
            newangle.k=taz1_hif.ANGLE[ii].k;
            newangle.theta=taz1_hif.ANGLE[ii].theta;
            angle.push_back(newangle);
        }
    }
    for(int ii=0; ii<taz1_hif.NumOfDihedrals; ++ii){
        if(taz1_hif.ChainId[taz1_hif.DIHEDRAL[ii].a]=='B'
                && taz1_hif.ChainId[taz1_hif.DIHEDRAL[ii].b]=='B'){
            Dihedral newdihedral;
            newdihedral.a=taz1_hif.DIHEDRAL[ii].a+offset;
            newdihedral.b=taz1_hif.DIHEDRAL[ii].b+offset;
            newdihedral.c=taz1_hif.DIHEDRAL[ii].c+offset;
            newdihedral.d=taz1_hif.DIHEDRAL[ii].d+offset;
            newdihedral.k=taz1_hif.DIHEDRAL[ii].k;
            newdihedral.theta=taz1_hif.DIHEDRAL[ii].theta;
            newdihedral.multiplicity=taz1_hif.DIHEDRAL[ii].multiplicity;
            dihedral.push_back(newdihedral);
        }
    }
    for(int ii=0; ii<taz1_hif.NumOfNatives; ++ii){
        char chaini,chainj;
        chaini=taz1_hif.ChainId[taz1_hif.NATIVE[ii].i];
        chainj=taz1_hif.ChainId[taz1_hif.NATIVE[ii].j];
        if(chaini!=chainj){
            Native newnative;
            int newi;
            int newj;
            if(chaini=='A'){
                newi=taz1_hif.NATIVE[ii].i;
            }else if(chaini=='B'){
                newi=taz1_hif.NATIVE[ii].i+offset;
            }
            if(chainj=='A'){
                newj=taz1_hif.NATIVE[ii].j;
            }else if(chainj=='B'){
                newj=taz1_hif.NATIVE[ii].j+offset;
            }
            newnative.i=newi;
            newnative.j=newj;
            newnative.k=taz1_hif.NATIVE[ii].k;
            newnative.r=taz1_hif.NATIVE[ii].r;
            native.push_back(newnative);
        }
        if(chaini=='B'&&chainj=='B'){
            Native newnative;
            newnative.i=taz1_hif.NATIVE[ii].i+offset;
            newnative.j=taz1_hif.NATIVE[ii].j+offset;
            newnative.k=taz1_hif.NATIVE[ii].k;
            newnative.r=taz1_hif.NATIVE[ii].r;
            native.push_back(newnative);
        }
    }
    //cited2_taz1
    offset=ntaz1+nhif;
    for(int ii=0; ii<cited2_taz1.NumOfBonds; ++ii){
        if(cited2_taz1.ChainId[cited2_taz1.BOND[ii].i]=='A'
                && cited2_taz1.ChainId[cited2_taz1.BOND[ii].j]=='A'){
            Bond newbond;
            newbond.i=cited2_taz1.BOND[ii].i+offset;
            newbond.j=cited2_taz1.BOND[ii].j+offset;
            newbond.k=cited2_taz1.BOND[ii].k;
            newbond.r=cited2_taz1.BOND[ii].r;
            bond.push_back(newbond);
        }
    }
    for(int ii=0; ii<cited2_taz1.NumOfAngles; ++ii){
        if(cited2_taz1.ChainId[cited2_taz1.ANGLE[ii].a]=='A'
                && cited2_taz1.ChainId[cited2_taz1.ANGLE[ii].a]=='A'){
            Angle newangle;
            newangle.a=cited2_taz1.ANGLE[ii].a+offset;
            newangle.b=cited2_taz1.ANGLE[ii].b+offset;
            newangle.c=cited2_taz1.ANGLE[ii].c+offset;
            newangle.k=cited2_taz1.ANGLE[ii].k;
            newangle.theta=cited2_taz1.ANGLE[ii].theta;
            angle.push_back(newangle);
        }
    }
    for(int ii=0; ii<cited2_taz1.NumOfDihedrals; ++ii){
        if(cited2_taz1.ChainId[cited2_taz1.DIHEDRAL[ii].a]=='A'
                && cited2_taz1.ChainId[cited2_taz1.DIHEDRAL[ii].a]=='A'){
            Dihedral newdihedral;
            newdihedral.a=cited2_taz1.DIHEDRAL[ii].a+offset;
            newdihedral.b=cited2_taz1.DIHEDRAL[ii].b+offset;
            newdihedral.c=cited2_taz1.DIHEDRAL[ii].c+offset;
            newdihedral.d=cited2_taz1.DIHEDRAL[ii].d+offset;
            newdihedral.k=cited2_taz1.DIHEDRAL[ii].k;
            newdihedral.theta=cited2_taz1.DIHEDRAL[ii].theta;
            newdihedral.multiplicity=cited2_taz1.DIHEDRAL[ii].multiplicity;
            dihedral.push_back(newdihedral);
        }
    }
    for(int ii=0; ii<cited2_taz1.NumOfNatives; ++ii){
        char chaini,chainj;
        chaini=cited2_taz1.ChainId[cited2_taz1.NATIVE[ii].i];
        chainj=cited2_taz1.ChainId[cited2_taz1.NATIVE[ii].j];
        if(chaini!=chainj){
            Native newnative;
            int newi;
            int newj;
            if(chaini=='A'){
                newi=cited2_taz1.NATIVE[ii].i+offset;
            }else if(chaini=='B'){
                newi=cited2_taz1.NATIVE[ii].i-ncited2;
            }
            if(chainj=='A'){
                newj=cited2_taz1.NATIVE[ii].j+offset;
            }else if(chainj=='B'){
                newj=cited2_taz1.NATIVE[ii].j-ncited2;
            }
            newnative.i=newi;
            newnative.j=newj;
            newnative.k=cited2_taz1.NATIVE[ii].k;
            newnative.r=cited2_taz1.NATIVE[ii].r;
            native.push_back(newnative);
        }
        if(chaini=='A'&&chainj=='A'){
            Native newnative;
            newnative.i=cited2_taz1.NATIVE[ii].i+offset;
            newnative.j=cited2_taz1.NATIVE[ii].j+offset;
            newnative.k=cited2_taz1.NATIVE[ii].k;
            newnative.r=cited2_taz1.NATIVE[ii].r;
            native.push_back(newnative);
        }
    }

    OutMol->ZINCBOND=(Bond*)malloc(zincBond.size()*sizeof(Bond));
    OutMol->NumOfZincBonds=zincBond.size();
    i=0;
    for(Bond tmp : zincBond){
#ifdef testMerge
        printf("zincBond (%d,%d) k=%f r=%f\n",tmp.i,tmp.j,tmp.k,tmp.r);
#endif
        OutMol->ZINCBOND[i]=tmp;
        i++;
    }

    OutMol->ZINCANGLE=(Angle*)malloc(zincAngle.size()*sizeof(Angle));
    OutMol->NumOfZincAngles=zincAngle.size();
    i=0;
    for(Angle tmp : zincAngle){
#ifdef testMerge
        printf("zincAngle (%d,%d,%d) k=%f theta=%f\n",tmp.a,tmp.b,tmp.c,tmp.k,tmp.theta);
#endif
        OutMol->ZINCANGLE[i]=tmp;
        i++;
    }

    OutMol->BOND=(Bond*)malloc(bond.size()*sizeof(Bond));
    OutMol->NumOfBonds=bond.size();
    i=0;
    for(Bond tmp : bond){
#ifdef testMerge
        printf("bond (%d,%d) k=%f r=%f\n",tmp.i,tmp.j,tmp.k,tmp.r);
#endif
        OutMol->BOND[i]=tmp;
        i++;
    }

    OutMol->ANGLE=(Angle*)malloc(angle.size()*sizeof(Angle));
    OutMol->NumOfAngles=angle.size();
    i=0;
    for(Angle tmp : angle){
#ifdef testMerge
        printf("angle (%d,%d,%d) k=%f theta=%f\n",tmp.a,tmp.b,tmp.c,tmp.k,tmp.theta);
#endif
        OutMol->ANGLE[i]=tmp;
        i++;
    }
    
    OutMol->DIHEDRAL=(Dihedral*)malloc(dihedral.size()*sizeof(Dihedral));
    OutMol->NumOfDihedrals=dihedral.size();
    i=0;
    for(Dihedral tmp : dihedral){
#ifdef testMerge
        printf("diheral (%d,%d,%d,%d) k=%f theta=%f\n",tmp.a,tmp.b,tmp.c,tmp.d,tmp.k,tmp.theta);
#endif
        OutMol->DIHEDRAL[i]=tmp;
        i++;
    }

    OutMol->NATIVE=(Native*)malloc(native.size()*sizeof(Native));
    OutMol->NumOfNatives=native.size();
    i=0;
    for(Native tmp : native){
#ifdef testMerge
        printf("native (%d,%d) k=%f r=%f\n",tmp.i,tmp.j,tmp.k,tmp.r);
#endif
        OutMol->NATIVE[i]=tmp;
        i++;
    }
    printf("total number of atoms merged=%d\n",numOfAtoms);
    return 0;
}
//this function pick the molecule with chainId out of InMol
int splitGoMolecule(const Molecule *InMol, Molecule *OutMol,const std::vector<char> &chainIds){
    int numOfChains;
    int *atomsPerChain;
    int numOfAtoms;
    getNumOfChains(InMol,&numOfChains,&atomsPerChain);
    std::vector<int> chainOffset;
    std::vector<int> numOfAtomsPerChain;
    std::vector<int> newChainOffset;
    char currentChainId;
    int atomId=0;
    int newAtomId=0;
    for(int ii=0; ii<numOfChains; ++ii){
        int natoms=atomsPerChain[ii];
        currentChainId=InMol->ChainId[atomId];
        bool found=false;
        for(char chainId : chainIds){
            if(chainId==currentChainId){
                found=true;
                break;
            }
        }
        if(found){
            newChainOffset.push_back(newAtomId);
            chainOffset.push_back(atomId);
            numOfAtomsPerChain.push_back(natoms);
            newAtomId+=natoms;
        }
        atomId+=natoms;
    }
    int i=0;
    for(int offset: chainOffset){
        printf("offset%d=%d\n",i,offset);
        i++;
    }
    i=0;
    for(int newoffset: newChainOffset){
        printf("newoffset%d=%d\n",i,newoffset);
        i++;
    }
    i=0;
    for(int natoms: numOfAtomsPerChain){
        printf("natoms%d=%d\n",i,natoms);
        i++;
    }
    //determine the number of atoms
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        bool found=false;
        for(char chainId : chainIds){
            if(chainId==InMol->ChainId[ii]){
                found=true;
                break;
            }
        }
        if(found){
            numOfAtoms++;
        }
    }
    allocateGoMolecule(numOfAtoms,OutMol);
    //copy peratom properties
    atomId=0;
    currentChainId='A';
    char tmpChainId='0';
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        bool found=false;
        for(char chainId : chainIds){
            if(chainId==InMol->ChainId[ii]){
                found=true;
                if(tmpChainId!=chainId && tmpChainId!='0'){
                    currentChainId++;
                }
                tmpChainId=chainId;
                break;
            }
        }
        if(found){
            OutMol->X[atomId]=InMol->X[ii];
            OutMol->Y[atomId]=InMol->Y[ii];
            OutMol->Z[atomId]=InMol->Z[ii];
            OutMol->VDWR[atomId]=InMol->VDWR[ii];
            OutMol->EPS[atomId]=InMol->EPS[ii];
            OutMol->Charge[atomId]=InMol->Charge[ii];
            OutMol->Mass[atomId]=InMol->Mass[ii];
            OutMol->AtomNum[atomId]=atomId;
            OutMol->AtomType[atomId]=InMol->AtomType[ii];
            OutMol->ResType[atomId]=InMol->ResType[ii];
            OutMol->ChainId[atomId]=currentChainId;
            OutMol->ChainCHARMMId[atomId]=InMol->ChainCHARMMId[ii];
            OutMol->ResNum[atomId]=atomId;
            atomId++;
        }
    }
    //copy parameters info
    std::vector<Bond> bond;
    std::vector<Angle> angle;
    std::vector<Dihedral> dihedral;
    std::vector<Native> native;
    std::vector<Bond> zincBond;
    std::vector<Angle> zincAngle;
    //zinc bond
    for(int ii=0; ii<InMol->NumOfZincBonds; ++ii){
        int chainIndex=-1;
        int atomIndex1=InMol->ZINCBOND[ii].i;
        int atomIndex2=InMol->ZINCBOND[ii].j;
        int currentChainIndex=0;
        for(int chainId: chainIds){
            if(chainId==InMol->ChainId[atomIndex1]){
                chainIndex=currentChainIndex;
            }
            currentChainIndex++;
        }
        if(chainIndex>=0){
            Bond newbond;
            int offset=chainOffset[chainIndex]-newChainOffset[chainIndex];
            newbond.i=atomIndex1-offset;
            newbond.j=atomIndex2-offset;
            newbond.k=InMol->ZINCBOND[ii].k;
            newbond.r=InMol->ZINCBOND[ii].r;
            zincBond.push_back(newbond);
        }
    }
    //zinc angle
    for(int ii=0; ii<InMol->NumOfZincAngles; ++ii){
        int chainIndex=-1;
        int atomIndex1=InMol->ZINCANGLE[ii].a;
        int currentChainIndex=0;
        for(int chainId: chainIds){
            if(chainId==InMol->ChainId[atomIndex1]){
                chainIndex=currentChainIndex;
            }
            currentChainIndex++;
        }
        if(chainIndex>=0){
            Angle newangle;
            int offset=chainOffset[chainIndex]-newChainOffset[chainIndex];
            newangle.a=InMol->ZINCANGLE[ii].a-offset;
            newangle.b=InMol->ZINCANGLE[ii].b-offset;
            newangle.c=InMol->ZINCANGLE[ii].c-offset;
            newangle.k=InMol->ZINCANGLE[ii].k;
            newangle.theta=InMol->ZINCANGLE[ii].theta;
            zincAngle.push_back(newangle);
        }
    }
    //bond
    for(int ii=0; ii<InMol->NumOfBonds; ++ii){
        int chainIndex=-1;
        int atomIndex1=InMol->BOND[ii].i;
        int atomIndex2=InMol->BOND[ii].j;
        int currentChainIndex=0;
        for(int chainId: chainIds){
            if(chainId==InMol->ChainId[atomIndex1]){
                chainIndex=currentChainIndex;
            }
            currentChainIndex++;
        }
        if(chainIndex>=0){
            Bond newbond;
            int offset=chainOffset[chainIndex]-newChainOffset[chainIndex];
            newbond.i=atomIndex1-offset;
            newbond.j=atomIndex2-offset;
            newbond.k=InMol->BOND[ii].k;
            newbond.r=InMol->BOND[ii].r;
            bond.push_back(newbond);
        }
    }
    //angle
    for(int ii=0; ii<InMol->NumOfAngles; ++ii){
        int chainIndex=-1;
        int atomIndex1=InMol->ANGLE[ii].a;
        int currentChainIndex=0;
        for(int chainId: chainIds){
            if(chainId==InMol->ChainId[atomIndex1]){
                chainIndex=currentChainIndex;
            }
            currentChainIndex++;
        }
        if(chainIndex>=0){
            Angle newangle;
            int offset=chainOffset[chainIndex]-newChainOffset[chainIndex];
            newangle.a=InMol->ANGLE[ii].a-offset;
            newangle.b=InMol->ANGLE[ii].b-offset;
            newangle.c=InMol->ANGLE[ii].c-offset;
            newangle.k=InMol->ANGLE[ii].k;
            newangle.theta=InMol->ANGLE[ii].theta;
            angle.push_back(newangle);
        }
    }
    //dihedral
    for(int ii=0; ii<InMol->NumOfDihedrals; ++ii){
        int chainIndex=-1;
        int atomIndex1=InMol->DIHEDRAL[ii].a;
        int currentChainIndex=0;
        for(int chainId: chainIds){
            if(chainId==InMol->ChainId[atomIndex1]){
                chainIndex=currentChainIndex;
            }
            currentChainIndex++;
        }
        if(chainIndex>=0){
            Dihedral newdihedral;
            int offset=chainOffset[chainIndex]-newChainOffset[chainIndex];
            newdihedral.a=InMol->DIHEDRAL[ii].a-offset;
            newdihedral.b=InMol->DIHEDRAL[ii].b-offset;
            newdihedral.c=InMol->DIHEDRAL[ii].c-offset;
            newdihedral.d=InMol->DIHEDRAL[ii].d-offset;
            newdihedral.k=InMol->DIHEDRAL[ii].k;
            newdihedral.multiplicity=InMol->DIHEDRAL[ii].multiplicity;
            newdihedral.theta=InMol->DIHEDRAL[ii].theta;
            dihedral.push_back(newdihedral);
        }
    }
    for(int ii=0; ii<InMol->NumOfNatives; ++ii){
        int chainIndex1=-1;
        int chainIndex2=-1;
        int atomIndex1=InMol->NATIVE[ii].i;
        int atomIndex2=InMol->NATIVE[ii].j;
        int currentChainIndex=0;
        for(int chainId: chainIds){
            if(chainId==InMol->ChainId[atomIndex1]){
                chainIndex1=currentChainIndex;
            }
            if(chainId==InMol->ChainId[atomIndex2]){
                chainIndex2=currentChainIndex;
            }
            currentChainIndex++;
        }
        if(chainIndex1>=0 && chainIndex2>=0){
            Native newnative;
            int offset1=chainOffset[chainIndex1]-newChainOffset[chainIndex1];
            int offset2=chainOffset[chainIndex2]-newChainOffset[chainIndex2];
            newnative.i=InMol->NATIVE[ii].i-offset1;
            newnative.j=InMol->NATIVE[ii].j-offset2;
            newnative.k=InMol->NATIVE[ii].k;
            newnative.r=InMol->NATIVE[ii].r;
            native.push_back(newnative);
        }
    }

    OutMol->ZINCBOND=(Bond*)malloc(zincBond.size()*sizeof(Bond));
    OutMol->NumOfZincBonds=zincBond.size();
    i=0;
    for(Bond tmp : zincBond){
#ifdef testSplit
        printf("zincBond (%d,%d) k=%f r=%f\n",tmp.i,tmp.j,tmp.k,tmp.r);
#endif
        OutMol->ZINCBOND[i]=tmp;
        i++;
    }

    OutMol->ZINCANGLE=(Angle*)malloc(zincAngle.size()*sizeof(Angle));
    OutMol->NumOfZincAngles=zincAngle.size();
    i=0;
    for(Angle tmp : zincAngle){
#ifdef testSplit
        printf("zincAngle (%d,%d,%d) k=%f theta=%f\n",tmp.a,tmp.b,tmp.c,tmp.k,tmp.theta);
#endif
        OutMol->ZINCANGLE[i]=tmp;
        i++;
    }

    OutMol->BOND=(Bond*)malloc(bond.size()*sizeof(Bond));
    OutMol->NumOfBonds=bond.size();
    i=0;
    for(Bond tmp : bond){
#ifdef testSplit
        printf("bond (%d,%d) k=%f r=%f\n",tmp.i,tmp.j,tmp.k,tmp.r);
#endif
        OutMol->BOND[i]=tmp;
        i++;
    }

    OutMol->ANGLE=(Angle*)malloc(angle.size()*sizeof(Angle));
    OutMol->NumOfAngles=angle.size();
    i=0;
    for(Angle tmp : angle){
#ifdef testSplit
        printf("angle (%d,%d,%d) k=%f theta=%f\n",tmp.a,tmp.b,tmp.c,tmp.k,tmp.theta);
#endif
        OutMol->ANGLE[i]=tmp;
        i++;
    }
    
    OutMol->DIHEDRAL=(Dihedral*)malloc(dihedral.size()*sizeof(Dihedral));
    OutMol->NumOfDihedrals=dihedral.size();
    i=0;
    for(Dihedral tmp : dihedral){
#ifdef testSplit
        printf("diheral (%d,%d,%d,%d) k=%f theta=%f\n",tmp.a,tmp.b,tmp.c,tmp.d,tmp.k,tmp.theta);
#endif
        OutMol->DIHEDRAL[i]=tmp;
        i++;
    }

    OutMol->NATIVE=(Native*)malloc(native.size()*sizeof(Native));
    OutMol->NumOfNatives=native.size();
    i=0;
    for(Native tmp : native){
#ifdef testSplit
        printf("native (%d,%d) k=%f r=%f\n",tmp.i,tmp.j,tmp.k,tmp.r);
#endif
        OutMol->NATIVE[i]=tmp;
        i++;
    }
    return 0;
}

//scale the native contact strengths by chainIds
int scaleNative(Molecule *InMol,const std::vector<nativeScaler> &newScales){
    for(nativeScaler newScale : newScales){
        printf("Scale (chain%c,chain%c) interactions by %f\n",
                newScale.first,newScale.second,newScale.factor);
        for(int ii=0; ii<InMol->NumOfNatives; ++ii){
            float oldE;
            int atom1=InMol->NATIVE[ii].i;
            int atom2=InMol->NATIVE[ii].j;
            char chain1=InMol->ChainId[atom1];
            char chain2=InMol->ChainId[atom2];
            if((chain1==newScale.first && chain2==newScale.second) ||
                    (chain2==newScale.first && chain1==newScale.second)){
                oldE=InMol->NATIVE[ii].k;
                InMol->NATIVE[ii].k *= newScale.factor;
#ifdef testScaleNative
                printf("Scale atom%d%c atom%d%c from %f to %f by %f\n",atom1,
                        chain1,atom2,chain2,oldE,InMol->NATIVE[ii].k,newScale.factor);
#endif
            }
        }
    }
    return 0;
}

int scaleNativepKID1(Molecule *InMol,const std::vector<nativeScaler> &newScales){
    const int offset=16; //residues numbered from 16 to the end of pKID
    int NumOfChains;
    int *NumOfAtomsPerChain;
    getNumOfChains(InMol, &NumOfChains, &NumOfAtomsPerChain);
    for(nativeScaler newScale : newScales){
        printf("Scale (chain%c,chain%c) interactions by %f\n",
                newScale.first,newScale.second,newScale.factor);
        int cursor=0;
        int id=newScale.first-'A';
        for(int ii=0; ii<id; ++ii){
            cursor+=NumOfAtomsPerChain[ii];
        }
        for(int ii=0; ii<InMol->NumOfNatives; ++ii){
            float oldE;
            int atom1=InMol->NATIVE[ii].i;
            int atom2=InMol->NATIVE[ii].j;
            char chain1=InMol->ChainId[atom1];
            char chain2=InMol->ChainId[atom2];
            if((chain1==newScale.first && chain2==newScale.second) ||
                    (chain2==newScale.first && chain1==newScale.second)){
                if(atom1 <= (cursor+offset) && (atom2 <= cursor+offset)){
                    oldE=InMol->NATIVE[ii].k;
                    InMol->NATIVE[ii].k *= newScale.factor;
#ifdef testScaleNative
                    printf("Scale atom%d%c atom%d%c from %f to %f by %f\n",atom1,
                            chain1,atom2,chain2,oldE,InMol->NATIVE[ii].k,newScale.factor);
#endif
                }
            }
        }
    }
    return 0;
}
int scaleNativepKID2(Molecule *InMol,const std::vector<nativeScaler> &newScales){
    const int offset=16; //residues numbered from 16 to the end of pKID
    int NumOfChains;
    int *NumOfAtomsPerChain;
    getNumOfChains(InMol, &NumOfChains, &NumOfAtomsPerChain);
    for(nativeScaler newScale : newScales){
        printf("Scale (chain%c,chain%c) interactions by %f\n",
                newScale.first,newScale.second,newScale.factor);
        int cursor=0;
        int id=newScale.first-'A';
        for(int ii=0; ii<id; ++ii){
            cursor+=NumOfAtomsPerChain[ii];
        }
        for(int ii=0; ii<InMol->NumOfNatives; ++ii){
            float oldE;
            int atom1=InMol->NATIVE[ii].i;
            int atom2=InMol->NATIVE[ii].j;
            char chain1=InMol->ChainId[atom1];
            char chain2=InMol->ChainId[atom2];
            if((chain1==newScale.first && chain2==newScale.second) ||
                    (chain2==newScale.first && chain1==newScale.second)){
                if(atom1 > (cursor+offset) && (atom2 > cursor+offset)){
                    oldE=InMol->NATIVE[ii].k;
                    InMol->NATIVE[ii].k *= newScale.factor;
#ifdef testScaleNative
                    printf("Scale atom%d%c atom%d%c from %f to %f by %f\n",atom1,
                            chain1,atom2,chain2,oldE,InMol->NATIVE[ii].k,newScale.factor);
#endif
                }
            }
        }
    }
    return 0;
}

int calcZincFinger(const Molecule *InMol, Molecule *OutMol){
    const float PI = 4*atan(1);
    const float RAD2DEG = 180.0/PI;
    const float ZincBondStrength = 50;
    const float ZincAngleStrength = 30;
    //ZincBond
    const float ZincBondCutoff=2.5;
    int numOfZinc=0;
    std::vector<int> Zinc;
    std::vector<Bond> ZincBond;
    for(int ii=0; ii<OutMol->NumOfAtoms; ++ii){
        if(!mystrcmp(OutMol->AtomType[ii],"ZN")){
            numOfZinc++;
            Zinc.push_back(ii);
        }
    }
    for(int ii=0; ii<InMol->NumOfAtoms; ++ii){
        if(!mystrcmp(InMol->AtomType[ii],"ZN")){
            printf("**************************\n");
            float x0,y0,z0;
            x0=InMol->X[ii];
            y0=InMol->Y[ii];
            z0=InMol->Z[ii];
            int currentResNum=0;
            std::set<int> bondingResNum;
            for(int jj=0; jj<InMol->NumOfAtoms; ++jj){
                char *resname;
                char resletter;
                float x1,y1,z1;
                float mindist;
                x1=InMol->X[jj];
                y1=InMol->Y[jj];
                z1=InMol->Z[jj];
                mindist=sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0));
                mytrim(InMol->ResType[jj],&resname);
                resletter=KBGoParams::AminoThreeToOne[resname];
                if(mindist < ZincBondCutoff && (resletter=='C' || resletter=='H')){
                    bondingResNum.insert(InMol->ResNum[jj]);
                    //printf("Atom%s%d Residue%d %s %f\n",InMol->AtomType[jj],jj,InMol->ResNum[jj],InMol->ResType[jj],mindist);
                }
            }

            for(int resnum : bondingResNum){
                printf("Residue%d is bonding with Zinc%d\n",resnum,InMol->ResNum[ii]);
            }

            for(int jj=0; jj<InMol->NumOfAtoms; ++jj){
                char *resname;
                char resletter;
                float x1,y1,z1;
                float mindist;
                x1=InMol->X[jj];
                y1=InMol->Y[jj];
                z1=InMol->Z[jj];
                mindist=sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0));
                if(mystrcmp(InMol->AtomType[jj],"CA")==0){
                    for(int resnum : bondingResNum){
                        if(InMol->ResNum[jj] == resnum){
                            Bond bond;
                            bond.i=InMol->ResNum[jj];
                            bond.j=InMol->ResNum[ii];
                            bond.r=mindist;
                            bond.k=ZincBondStrength;
                            ZincBond.push_back(bond);
                        }
                    }
                }//endfor
            }//endfor
        } //endif
    }
    OutMol->NumOfZincBonds=ZincBond.size();
    OutMol->ZINCBOND=(Bond*)malloc(ZincBond.size()*sizeof(Bond));
    int i=0;
    for(Bond bond : ZincBond){
        printf("ZincBond(%d,%d) r=%f\n",bond.i,bond.j,bond.r);
        OutMol->ZINCBOND[i]=bond;
        i++;
    }
    //ZincAngle
    std::vector<Angle> ZincAngle;
    for(int id : Zinc){
        std::vector<int> nearbyAtoms;
        for(Bond bond : ZincBond){
            if(id==bond.j){
                nearbyAtoms.push_back(bond.i);
            }
        }
        for(int ii=0; ii < nearbyAtoms.size(); ++ii){
            for(int jj=ii+1; jj< nearbyAtoms.size(); ++jj){
                Angle angle;
                angle.a=nearbyAtoms[ii];
                angle.b=id;
                angle.c=nearbyAtoms[jj];
                angle.k=ZincAngleStrength;
                float theta=0;
                float x1,y1,z1,x2,y2,z2,x3,y3,z3;
                x1=OutMol->X[angle.a];
                y1=OutMol->Y[angle.a];
                z1=OutMol->Z[angle.a];
                x2=OutMol->X[angle.b];
                y2=OutMol->Y[angle.b];
                z2=OutMol->Z[angle.b];
                x3=OutMol->X[angle.c];
                y3=OutMol->Y[angle.c];
                z3=OutMol->Z[angle.c];
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
                theta=acos(cosTheta);
                angle.theta=theta*RAD2DEG;
                ZincAngle.push_back(angle);
                printf("ZincAngle(%d,%d,%d) k=%f theta=%f\n",
                        angle.a,angle.b,angle.c,angle.k,angle.theta);
            }
        }
        //change the histidine charge
        for(int idx : nearbyAtoms){
            char *resname;
            char resletter;
            mytrim(OutMol->ResType[idx],&resname);
            resletter=KBGoParams::AminoThreeToOne[resname];
            if(resletter=='H'){
                printf("Histidine %c%d charge was changed from %f to 0\n",
                        OutMol->ChainId[idx],idx,OutMol->Charge[idx]);
                OutMol->Charge[idx]=0;
            }
        }
    }
    OutMol->NumOfZincAngles=ZincAngle.size();
    OutMol->ZINCANGLE=(Angle*)malloc(ZincAngle.size()*sizeof(Angle));
    int ii=0;
    for(Angle angle: ZincAngle){
       OutMol->ZINCANGLE[ii]=angle; 
       ii++;
    }
    
    return 0;
}
#endif
