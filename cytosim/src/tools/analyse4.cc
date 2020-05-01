//RCS: $Id: analyse4.cc,v 2.1 2004/12/02 16:20:45 nedelec Exp $
//---------------------------analyse4.cc------------------------
//provides general statistics about complex binding and microtubules

#include <cstring>
#include "sim.h"
#include "reader.h"
#include "iowrapper.h"

int count;

int nbasters[4], emptyX, MTVLinks[4], MTVLinksPop[4];
int junk[1000];
int biggestaster[4][3];
 
int MTLinks[4][3], MTLinksPop[4][3];
real MTpop[4][3], MTVPop[4];

int monb[4][3][3], bundlenb[3][4];
real atime=0, meanlength, dp, dm;


void reset()
{
  int ty, i, r, e;

  for (ty=0; ty<4; ++ty)
    {
      for(r=0; r<3; ++r) 
	{
	  bundlenb[r][ty]=0;
	  for(e=0; e<3; ++e) monb[ty][r][e]=0;
	}
    }
  meanlength=0;
  dp=0;
  dm=0;
  
  for(ty=0; ty<4; ++ty) {
    
    for(i=0; i<3; ++i) { 
      MTLinks[ty][i]=0; 
      MTLinksPop[ty][i]=0;
      MTpop[ty][i]=0;
      biggestaster[ty][i]=0;
    }
    
    nbasters[ty]=0;
    MTVLinks[ty]=0;
    MTVLinksPop[ty]=0;
    MTVPop[ty]=0;
    }
  
  emptyX=0;

  count=0;
  atime=0;
}

//=========================================================================
//=========================================================================
//=========================================================================

void CountComplex()
{
  Complex * mo;
  int ty, endh;

  for(mo=sim.firstFreeComplex(); mo ; mo=mo->next() )
    {
      ty=mo->getType();
      monb[ty][0][0]++;
    }

  for(mo=sim.firstBoundComplex(); mo ; mo=mo->next() )
    {
      ty=mo->getType();
      endh=(mo->getEnd1()!=MT_NOT_END) + (mo->getEnd2()!=MT_NOT_END);

      ++( monb[ty][1][endh] );
    }
  
  for(mo=sim.firstBridgeComplex(); mo ; mo=mo->next() )
    {
      ty=mo->getType();
      endh=(mo->getEnd1()!=MT_NOT_END) + (mo->getEnd2()!=MT_NOT_END);

      ++( monb[ty][2][endh] );

      if (mo->cosAngle() > 0.98) ++bundlenb[endh][ty];
    }
}

//=========================================================================
//=========================================================================
//=========================================================================


void Topology2(int type)
{
    Complex * mo;
    int i, j, k, endh;
    int mt1, mt2;
    int a, nbmts, nblinks;
    
    int NbMT=sim.nbMicrotub()+1;
    int * asters = new int[NbMT];
    int *** MTMT = new int ** [3];
    MTMT[0] = new int * [NbMT];
    MTMT[1] = new int * [NbMT];
    MTMT[2] = new int * [NbMT];
    
    for(endh=0; endh<3; ++endh)
        for(i=0; i<NbMT; ++i)
        {
            MTMT[endh][i] = new int[NbMT];
            for(j=0; j<NbMT; ++j) MTMT[endh][i][j]=0;
        }
            
            for(mo=sim.firstBridgeComplex(); mo ; mo=mo->next() )
            {
                if ( mo->getType() != type ) continue;
                
                endh=(mo->getEnd1()!=MT_NOT_END) + (mo->getEnd2()!=MT_NOT_END);
                
                mt1 = mo->getMT1()->getName();
                mt2 = mo->getMT2()->getName();
                
                if (  (mt1>=NbMT) ||  (mt2>=NbMT) ) { printf("d"); continue; }
                ++MTMT[endh][mt1][mt2];
                ++MTMT[endh][mt2][mt1];
            }
            
            for(i=0; i<NbMT; ++i)
                for(j=0; j<NbMT; ++j)
                    for(endh=2; endh>=0; --endh)
                        if (MTMT[endh][i][j])
                        { 
                            ++MTLinks[type][endh];
                            MTLinksPop[type][endh]+=MTMT[0][i][j]+MTMT[1][i][j]+MTMT[2][i][j]; 
                            endh=-1;
                        } 
                            
                            //----------------------find the clusters of connected MTs:
                            for(i=1; i<NbMT; ++i) asters[i]=i;
    
    do  {
        a=0;
        for(j=1; j<NbMT; ++j) 
            for(k=j+1; k<NbMT; ++k)
                if ( ( MTMT[2][j][k] ) && ( asters[j] != asters[k] ) )
                {
                    a = minT( asters[j], asters[k] );
                    asters[j] = a;
                    asters[k] = a;
                }
    } while (a);
        
        //----------------------print-out the asters:
        for(i=1; i<NbMT; ++i)
            if ((a = asters[i]))
            {
                nbmts=0;
                nblinks=0;
                for(j=i; j<NbMT; ++j)
                    if (asters[j]==a)
                    { 
                        ++nbmts; 
                        asters[j]=0;                    //to avoid counting it twice
                        for(k=1; k<NbMT; ++k) nblinks+=MTMT[2][j][k];
                    }
                        
                        if ( nbmts > 1 )
                        {
                            MTVLinks[type]    += 2*(nbmts - 1); 
                            MTVLinksPop[type] += nblinks; 
                        }
                        
                        if ( nbmts > 2 )
                        {
                            ++nbasters[type];
                            //printf("%i : %i\n",type, nbmts);
                            for( k=0; k<3; ++k)
                                if ( nbmts > biggestaster[type][k] )
                                {
                                    for(j=2; j>k; --j)
                                        biggestaster[type][j] = biggestaster[type][j-1];
                                    biggestaster[type][k] = nbmts;
                                    break;
                                }
                        }
            }
                        
                        for(endh=0; endh<3; ++endh)
                            for(i=0; i<NbMT; ++i)
                                delete [] MTMT[endh][i];
                
                delete [] MTMT[0];
                delete [] MTMT[1];
                delete [] MTMT[2];
                delete [] MTMT;
                delete [] asters;
}


//=========================================================================
//=========================================================================
//=========================================================================


void MTDistance()
{
    Microtub * mt1, * mt2;
    Vecteur x;
    
    for (mt1=sim.firstMicrotub(); mt1 ; mt1=mt1->next() )
        for (mt2=mt1, mt2=mt2->next();  mt2 ; mt2=mt2->next() )
        {
            //x = mt1->whereendEnd(MT_PLUS_END) - mt2->whereendEnd(MT_PLUS_END);
            //Space::modulo(x);
            //dp += x.norm();
            
            x = mt1->whereEnd(MT_MINUS_END) - mt2->whereEnd(MT_MINUS_END);
            Space::modulo(x);
            dm += x.norm();
        }
            
    dp /= sqr( sim.nbMicrotub() );
    dm /= sqr( sim.nbMicrotub() );
}


//=========================================================================
//=========================================================================
//=========================================================================

void Report(int ty)
{
    int endh;
    int Tb, Tl;
    real cnt=count;
    
    if (count==0) return;
    
    for(endh=0; endh<3; ++endh)
    {
        MTLinks[ty][endh] /= 2; 
        MTLinksPop[ty][endh] /= 2; 
        if ( MTLinks[ty][endh] ) 
            MTpop[ty][endh]=MTLinksPop[ty][endh]/real(MTLinks[ty][endh]);
    }
    MTVLinks[ty]    /= 2; 
    MTVLinksPop[ty] /= 2;
    if (MTVLinks[ty]) MTVPop[ty]=MTVLinksPop[ty]/real(MTVLinks[ty]);  
    
    printf("MT%i %5.1f %4.1f : ", ty, atime/cnt, meanlength/cnt);
    printf("X %5.1f %4.1f  T %5.1f %4.1f  V %5.1f %4.1f", 
           MTLinks[ty][0]/cnt, MTpop[ty][0],
           MTLinks[ty][1]/cnt, MTpop[ty][1],
           MTLinks[ty][2]/cnt, MTpop[ty][2]);
    printf(" E %7.1f", emptyX/cnt);
    
    printf(" MTV %5.1f %4.1f", MTVLinks[ty]/cnt, MTVPop[ty]);
    printf(" asters %i %i %i ", biggestaster[ty][0], 
           biggestaster[ty][1], biggestaster[ty][2] );
    if (nbasters[ty])
        printf(" %5.3f  %5.3f",
               MTVLinks[ty]/real(nbasters[ty]*cnt),
               MTVLinksPop[ty]/real(nbasters[ty]*cnt) ); 
    else printf(" 0  0 ");
    //printf(" End-End  %8.5f %8.5f", dp/cnt, dm/cnt);
    printf("\n");
    
    
    Tb=monb[ty][1][0]+monb[ty][1][1];
    Tl=monb[ty][2][0]+monb[ty][2][1]+monb[ty][2][2];
    printf("%i   %5.1f : ", ty, atime/cnt);
    printf("%7.1f %7.1f (ME %7.1f %7.1f )",
           monb[ty][0][0]/cnt, 
           Tb/cnt, monb[ty][1][0]/cnt, monb[ty][1][1]/cnt);
    printf(" %7.1f (XTVB %7.1f %7.1f %7.1f  %7.1f )",
           Tl/cnt, monb[ty][2][0]/cnt, monb[ty][2][1]/cnt, 
           monb[ty][2][2]/cnt, 
           (bundlenb[0][ty]/cnt+bundlenb[1][ty]+bundlenb[2][ty])/cnt); 
    printf( "%7.1f",  bundlenb[0][ty]/cnt);
    printf("\n");
}


void Report()
{
    Report(0);
    Report(1);
    reset();
}

//=========================================================================
//=========================================================================
//=========================================================================

int main(int argc, char * argv[])
{
    int begin=0, end=1000;
    MSG.shutUp();
    
    if (argc==2)
        if (1 == sscanf(argv[1],"%u-%u",&begin,&end))
            end=begin;
    
    if (argc==3)
        if (1 == sscanf(argv[2],"%u-%u",&begin,&end))
            end=begin;
    
    SimReader reader;
    if ( NO_ERROR != reader.openFile( RESULT_OUT ) ) {
        printf("analyse4.cc: cannot open %s\n", RESULT_OUT);
        return 0;
    }
    
    if ( NO_ERROR != MP.parseFile(DATA_OUT) )
        if ( NO_ERROR != MP.parseFile(MP.datafile) )
            MSG.error("analyse.cc","Cannot read %s or data.out", MP.datafile);
    MP.massageParameters();
    
    for( int frame = begin; frame <= end ; ++frame )
        if ( NO_ERROR == reader.readFrame( frame ) )
        {
            ++count;
            CountComplex();
            //MTDistance();
            Topology2(0);
            Topology2(1);
            //CountMTCrossings();
            meanlength += sim.meanMTLength();
            atime += sim.simTime();
            if ( argc >= 2 )
            {
                if (argv[1][0]=='a') Report();
                if ( (argv[1][0]=='b') && ( 0 == (count % 2 ) ) ) Report();
                if ( (argv[1][0]=='c') && ( 0 == (count % 4 ) ) ) Report();
                if ( (argv[1][0]=='d') && ( 0 == (count % 8) ) )  Report();
                if ( (argv[1][0]=='e') && ( 0 == (count % 16) ) ) Report();
            }
        }
            
    Report();
    return 0;
}
