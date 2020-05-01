//RCS: $Id: analyse2.cc,v 2.4 2005/04/18 10:25:54 nedelec Exp $
//------------------------------analyse2.cc----------------------------


#include <cstring>
#include "sim.h"
#include "reader.h"
#include "iowrapper.h"
#include "map.h"

//=====================================================================

Map<1, real> mtlength;
Map<1, real> haforce[6];
Map<1, real> haload[6];

int mt_state[10];

void setUp()
{
  real limit;
  mtlength.create( 0, 40, 100 ); //minT( MP.mtmaxlength, 2 * MP.mtmonlength) );
  mtlength.clear();

  for( int ty=0; ty<6; ++ty )
    {
      limit = 2.5 * MP.haforce[ty];

      haforce[ty].create( 0, limit, 50 );
      haforce[ty].clear();

      haload[ty].create( -limit, limit, 50 );
      haload[ty].clear();
    }

  mt_state[ MT_SHRINKING ] = 0;
  mt_state[ MT_GROWING ] = 0;
  mt_state[ 2 ] = 0;
  
}


//=====================================================================



void MTStats()
{
  for ( Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() )
    ++mtlength( mt->length() );

  for ( Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() )
    {
      if ( mt->length() < MP.mtminlength + 0.05 )
	++mt_state[2];
      else
	++mt_state[ mt->getDynamicState(MT_PLUS_END) ];
    }
}

//=====================================================================

void handsForce()
{
  Vecteur v;
  int hty;

  for ( Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() )
    for(Hand * h=mt->firstHand(); h ; h=h->next() )
      if (h->otherBound())
	{
	  hty = h->getType();

	  v = MP.km * h->getStretch();
	  //printf("x=%10.4f y=%10.4f z=%10.4f  ", v.XX, v.YY, v.ZZ );

	  ++ haforce[ hty ]( v.norm() );
	  ++ haload[ hty ]( v * h->dirMT() );
	}
}


//=====================================================================
//=====================================================================


void report()
{
  char filename[128];

  FILE * outfile = fopen("mtl.out","w");
  mtlength.write( outfile, 0 );
  fclose( outfile );

  for( int ty=0; ty < 6; ++ty )
    if ( haforce[ty].sum() > 1 ) {
        snprintf(filename, sizeof(filename), "haf%i.out", ty );
        outfile = fopen(filename, "w");
        haforce[ty].write(outfile, 0);
        fclose( outfile );
        
        snprintf(filename, sizeof(filename), "hal%i.out", ty );
        outfile = fopen(filename, "w");
        haload[ty].write(outfile, 0);
        fclose( outfile );
    }

  printf("MT at minlength %8i  MT_GROWING %8i MT_SHRINKING %8i\n",
	 mt_state[2], mt_state[MT_GROWING], mt_state[MT_SHRINKING]);


  real total = mt_state[2] + mt_state[MT_GROWING] + mt_state[MT_SHRINKING];


  printf("MT at minlength %8.2f  MT_GROWING %8.2f MT_SHRINKING %8.2f (%%)\n",
	 mt_state[2] * 100.0 / total,
	 mt_state[MT_GROWING] * 100.0 / total,
	 mt_state[MT_SHRINKING] * 100.0 / total);
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
    { printf("analyse2.cc: cannot open [%s]\n", RESULT_OUT); return 0; }

  MSG.setVerboseLevel(1);

  setUp();

  for( int frame = begin; frame <= end ; ++frame )
    if ( NO_ERROR == reader.readFrame( frame ) ) {
      ++cnt;
      MTStats();
      handsForce();
    }

  report();
  printf("%i images read\n", cnt);

  return 0;
}


