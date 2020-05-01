//RCS: $Id: analyse0.cc,v 2.6 2005/04/18 10:25:47 nedelec Exp $
//---------------------------analyse0.cc------------------------
// This is a template for a program to analyse simulation outputs
// reads result.out, calculates some characteristics and print them

#include "sim.h"
#include "reader.h"
#include "iowrapper.h"
#include "iomessages.h"

const unsigned int STRLEN = 512;
char commands[STRLEN];

int first_frame = 0;
int last_frame  = SimReader::MAX_FRAME;

//===================================================================

//Two examples of 'analyse' functions are below,
//write and add your own analyse procedures here

void countPoints()
{
  int result = 0;
  for (Microtub *  mt=sim.firstMicrotub(); mt ; mt=mt->next() )
    result += mt->nbSegments();
  
  printf("F %4i  %.3f s. %5i points\n", sim.frameInBuffer(), sim.simTime(), result);
}


//-------------------------------------------------------------------
void microtubTopology()
{
  Name nn, aa, ii, jj, kk;
  Name mt1, mt2;

  ///rename the names in a continuous manner:
  sim.renameMicrotub();
  Name NbMT = sim.nbMicrotub() + 1;

  
  int * MTMT = new int[NbMT*NbMT];
  for( nn = 0; nn < NbMT*NbMT; ++nn)
    MTMT[ nn ]=0;
  
  //------- set the array MTMT: which MT is connected to which
  for(Complex * mo=sim.firstBridgeComplex(); mo ; mo=mo->next() ) {
        
    mt1 = mo->getMT1()->getName();
    mt2 = mo->getMT2()->getName();
    //this is an error condition:
    if (  ( mt1 >= NbMT ) || ( mt2 >= NbMT ) ) {
      MSG.error("analyse0.cc", "internal error in microtubTopology()");
      continue; 
    }
    ++MTMT[ mt1 + NbMT * mt2 ];
    ++MTMT[ mt2 + NbMT * mt1 ];
  }
  
  //------- find the clusters of connected MTs:
  Name * cluster = new Name[NbMT];

  for( nn = 1; nn < NbMT; ++nn) 
    cluster[ nn ] = 0;
  
  for(  Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() ) {
    if ( mt->getName() >= NbMT ) {
      MSG.error("analyse0.cc", "internal error in microtubTopology()");
      continue;
    }
    cluster[ mt->getName() ] = mt->getName();
  }
  
  ///calculate the clusters:
  do  {
    aa = 0;
    for( jj = 1;   jj < NbMT; ++jj) 
      for( kk = jj+1; kk < NbMT; ++kk)
        if ( ( MTMT[ jj + NbMT * kk ] ) && ( cluster[ jj ] != cluster[ kk ] ) ) {
          aa = minT( cluster[ jj ], cluster[ kk ] );
          cluster[ jj ] = aa;
          cluster[ kk ] = aa;
        }
  } while (aa);
    
  //------------------------print out the clusters:
   int nbclusters=0;
   int nbpairs=0;
   int biggest[3] = { 0, 0, 0 };
   Vecteur center;
    
    for(ii = 1; ii < NbMT; ++ii)
      if ( ( aa = cluster[ii] ) != 0 ) {
        int nbmts=0;
        int nblinks=0;
        center = VZERO;
        //find all the MT in this cluster:
        for( jj = ii ; jj < NbMT; ++jj )
          if ( cluster[jj] == aa ) { 
            ++nbmts; 
            center += sim.findMicrotub(jj)->whereEnd(MT_MINUS_END);
            for(kk=1; kk < NbMT; ++kk)
              nblinks += MTMT[ jj + NbMT * kk ];
            cluster[ jj ]=0;
          }
            
        nblinks /= 2;
        
        for( jj=0; jj<3; ++jj )
          if ( nbmts > biggest[jj] ) {
            for( kk=2; kk>jj; --kk)
              biggest[kk] = biggest[kk-1];
            biggest[jj] = nbmts;
            break;
          }
          
          
        if ( nbmts == 2 )
          ++nbpairs;
        
        if ( nbmts > 2 ) {
          ++nbclusters;
          center /= nbmts;
          printf("M%-4lx  nbmt %4i   links %6i    center %6.1f %6.1f\n", 
              aa, nbmts, nblinks, center.XX, center.YY);
          }
        }

    printf("nbclusters ");
    printf("%4i  nbpairs %4i   biggest %4i %4i %4i\n", 
           nbclusters, nbpairs, biggest[0], biggest[1], biggest[2]);
    
    delete [] MTMT;
    delete [] cluster;
}



//=====================================================================

void showHelp()
{
  printf("Cytosim / analyse:\n");
  printf("Reads [%s] and prints some statistics on the objects\n", RESULT_OUT);
  printf("USAGE: analyse [ -f(index) ] [-f(index1-index2)] [ command_name1 command_name2 ... ]\n");
  printf("\n");
  printf("Commands supported:\n");
  printf(" nbpoints         report number of points in the microtubules\n");
  printf(" topology         report topolgy of microtubule connections\n");
  printf("\n");
  printf("analyse0  DIM=%i, compiled %s, %s\n", DIM, __TIME__, __DATE__ );
}

//=====================================================================

void parseCommandLine(int argc, char * argv[])
{
  for(int ii = 1; ii < argc; ++ii ) {
    
    if ( '-' == argv[ii][0] ) {
      switch( argv[ii][1] ) {
        case 'h': showHelp(); exit(EXIT_SUCCESS);
        case 'v': MSG.setVerboseLevel( argv[ii]+2 ); continue;
        case 'f':   //frame index
          if ( argv[ii][2] )
            if ( 1 == sscanf(argv[ii]+2, "%u-%u", &first_frame, &last_frame ) )
              last_frame = first_frame; //if only one value set, make window square
          continue;
          
        default:
          printf("ignored option [%s] on command line\n", argv[ii] );
          continue;
      }
    }
    
    if ( strlen( commands ) + strlen( argv[ii] ) + 1 < STRLEN ) {
      strcat( commands, " " );
      strcat( commands, argv[ii] );
    }
  }
  strcat( commands, " " );
}

//=====================================================================
int main(int argc, char * argv[])
{
  if ( argc < 2 ) { 
    showHelp(); 
    return EXIT_SUCCESS;
  } 
  
  parseCommandLine(argc, argv);

  MSG.shutUp();
    
  if ( NO_ERROR != MP.parseFile(DATA_OUT) )
    if ( NO_ERROR != MP.parseFile(MP.datafile) )
      MSG.error("analyse.cc","Cannot read %s or data.out", MP.datafile);

  MP.massageParameters();

  SimReader reader;
  if ( NO_ERROR != reader.openFile( RESULT_OUT ) ) {
    fprintf(stderr, "analyse0.cc: cannot open %s\n", RESULT_OUT);
    return EXIT_FAILURE; 
  }

  for( int frame = first_frame; frame <= last_frame ; ++frame )
    if ( NO_ERROR == reader.readFrame( frame ) ) {
      
      if ( strstr(commands, " nbpoints " ))      countPoints();
      if ( strstr(commands, " topology " )) microtubTopology();
      //add a line here to call you analyse procedures
      
    }
  return EXIT_SUCCESS;
}
