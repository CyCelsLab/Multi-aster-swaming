//RCS: $Id: analyse3.cc,v 2.2 2005/04/18 10:26:02 nedelec Exp $
//------------------------------analyse2.cc----------------------------

#include <cstring>
#include "sim.h"
#include "reader.h"
#include "iowrapper.h"

//=====================================================================

int     ends_cnt[2][2];
Vecteur ends_pos[2][2];
real    ends_pos_sqr[2][2];
Vecteur dirs[2];
int nbmts;

//=====================================================================

void reset()
{
  nbmts = 0;
  for(int ii = 0; ii < 2 ; ++ii )
  for(int jj = 0; jj < 2 ; ++jj )
    {
      ends_pos[ii][jj] = VZERO;
      ends_pos_sqr[ii][jj] = 0;
      ends_cnt[ii][jj] = 0;
      dirs[ii] = VZERO;
    }
}


//=====================================================================


void stats()
{
  Vecteur wm, wp;

  for ( Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() )
    {
      ++nbmts;

      wm = mt->whereEnd(MT_MINUS_END);
      ends_cnt[ 0 ][ wm.XX > 0 ]++;
      ends_pos[ 0 ][ wm.XX > 0 ] += wm;
      ends_pos_sqr[ 0 ][ wm.XX > 0 ] += wm.normSquare();
 
      wp = mt->whereEnd(MT_PLUS_END);
      ends_cnt[ 1 ][ wp.XX > 0 ]++;
      ends_pos[ 1 ][ wp.XX > 0 ] += wp;
      ends_pos_sqr[ 1 ][ wp.XX > 0 ] += wp.normSquare();

      dirs[ wm.XX > 0 ] += ( wp - wm ) / mt->length();
    }
}



//=====================================================================
//=====================================================================


void report(int cnt)
{
  printf("analyse3 DIM = %i", DIM );
  printf("  compiled %s, %s\n", __TIME__, __DATE__ );
 
  for(int ii = 0; ii < 2 ; ++ii )
  for(int jj = 0; jj < 2 ; ++jj )
    if ( ends_cnt[ii][jj] ) 
      {
	ends_pos[ii][jj]     /= ends_cnt[ii][jj];  
	ends_pos_sqr[ii][jj] /= ends_cnt[ii][jj]; 
	ends_pos_sqr[ii][jj] -= ends_pos[ii][jj].normSquare();
      } 
    
  for(int ii = 0; ii < 2 ; ++ii )
    if ( ends_cnt[0][ii] ) 
      dirs[ii] /= ends_cnt[ 0 ][ ii ];


  printf("cnt ");
  printf(" L- %7i", ends_cnt[0][0]);
  printf(" R- %7i", ends_cnt[0][1]);
  printf(" L+ %7i", ends_cnt[1][0]);
  printf(" R+ %7i", ends_cnt[1][1]);
  printf("\n");

  printf("pos ");
  printf(" L- "); ends_pos[0][0].print();
  printf(" R- "); ends_pos[0][1].print();
  printf(" L+ "); ends_pos[1][0].print(); 
  printf(" R+ "); ends_pos[1][1].print();
  printf("\n");

  printf("var ");
  printf(" L- %7.2f ", ends_pos_sqr[0][0] );
  printf(" R- %7.2f ", ends_pos_sqr[0][1] );
  printf(" L+ %7.2f ", ends_pos_sqr[1][0] );
  printf(" R+ %7.2f ", ends_pos_sqr[1][1] );
  printf("\n");

  printf("dir ");
  printf(" L  "); dirs[0].print();
  printf(" R  "); dirs[1].print();
  printf("\n");

  printf("nbs frames %8i mts %8i\n", cnt, nbmts);

}


//=====================================================================
//=====================================================================


int main(int argc, char * argv[])
{
  int cnt=0, begin=0, end=2000;

  if (argc==2)
    if (1 == sscanf(argv[1],"%d-%d",&begin,&end))
      end=begin;
  
  MSG.shutUp();
  if ( NO_ERROR != MP.parseFile(DATA_OUT) )
    if ( NO_ERROR != MP.parseFile(MP.datafile) )
      MSG.error("analyse.cc","Cannot read %s or data.out", MP.datafile);
  MP.massageParameters();

  SimReader reader;
  if ( NO_ERROR != reader.openFile( RESULT_OUT ) )
    { printf("analyse3.cc: cannot open [%s]\n", RESULT_OUT); return 0; }

  MSG.setVerboseLevel(1);

  reset();

  for( int frame = begin; frame <= end ; ++frame )
    if ( NO_ERROR == reader.readFrame( frame ) ) {
      ++cnt;
      stats();
    }

  report(cnt);

  return 0;
}


