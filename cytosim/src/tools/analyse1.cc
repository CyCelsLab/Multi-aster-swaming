//RCS: $Id: analyse1.cc,v 2.4 2005/04/12 07:24:03 nedelec Exp $
//---------------------------analyse1.cc------------------------
//provides general statistics about complex binding and microtubules

#include <cstring>
#include <cctype>
#include "sim.h"
#include "reader.h"
#include "iowrapper.h"

const int LMAX = 1000;
const int SCR  = 100;

//=====================================================================
//=====================================================================


void Report(int nbc[], int scale=1)
{
  int limit=SCR, pool=1, i, j, tr=0, m=0, last=0, tr2;
  int s;

  for (i=0; i<LMAX; i++)
    { 
      tr+=nbc[i];
      m +=i*nbc[i];
      if (nbc[i]) last=i;
    }
  tr2=tr-nbc[0];
  if (tr2==0) tr2=1;

  MSG(":");
  if (tr == 0) { MSG("empty\n"); return;}

  pool = 1 + int( last / limit ) ;

  for (i=0; i<=last; i+=pool) 
    {
      for(s=0, j=0; j<pool; ++j) s+=nbc[i+j];
      if (scale)
	s  = int( 100*s/real( tr*pool ) );
      else
	s /= pool;
      if (s > 61) { MSG("@"); continue; }
      if (s > 35) { MSG("%c", s-36+'A'); continue; }
      if (s >  9) { MSG("%c", s-10+'a'); continue; }
      if (s == 0) { if (nbc[i]) MSG("."); else MSG(" "); continue; }
      MSG("%c", s+'0');
    }

  for (; i <= pool*limit; i+=pool) MSG(" ");
  MSG("| m %5.2f(%5.2f) mx %3d T %6d",
	real(m)/real(tr), real(m)/real(tr2), last, tr);
  if (scale) MSG(" (1=%.2f)",tr / 100.0);
  if (pool>1) MSG(" pool %i", pool);
  MSG("\n");
}

//=====================================================================
//=====================================================================

void CountRods()
{
  int i,last,nbs[LMAX];
  for (i=0; i<LMAX; i++) nbs[i]=0;  
  for (Microtub *  mt=sim.firstMicrotub(); mt ; mt=mt->next() )
    ++nbs[ mt->nbSegments() ];
   
  last=0;
  for (i=0; i<LMAX; i++) if (nbs[i]) last=i;
  MSG("nb mt[nb rods] (last=%i) :\n",last);
  for (i=0; i<=last; i++) MSG("%d ",nbs[i]);
  MSG("\n");
  MSG("NbR");
  Report(nbs);
}




//=====================================================================
//=====================================================================

void EachMTlength()
{
  for ( Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() )
    MSG("M%lx : %i %i : %7.5f \n", mt->getName(),
	  mt->getDynamicState(MT_MINUS_END), mt->getDynamicState(MT_PLUS_END), mt->length());
}

//---------------------------------------------------------------------------


void MesureMTLength()
{
  static real ml_old=0, time_old=0, mminus_old=0;

  int state[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  int max=0, min=0;
  real mtop, top[LMAX], mtip, tip[LMAX];
  real Vm, Vp, l, mminus, ml=0, mlsqr=0, var=0;
  int t, w;
  int MAX=minT(sim.nbMicrotub()/10, LMAX);
  

  for(t=0; t<LMAX; t++) { top[t]=0; tip[t]=MP.mtmaxlength; }

  for ( Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() )
    {
      ++state[ mt->getDynamicState(MT_PLUS_END) ];
      l       = mt->length();
      mminus += mt->abscissa( MT_MINUS_END );
      ml     += l;
      mlsqr  += sqr(l);
      if ( l < MP.mtminlength+MP.mtrodlength ) ++min;
      if ( l > MP.mtmaxlength-MP.mtrodlength ) ++max;

      for(t=0; t<MAX; t++)
	if ( l > top[t] )
	  {
	    for(w=MAX-1; w>t; --w) top[w]=top[w-1];
	    top[t]=l;
	    break;
	  }

      for(t=0; t<MAX; t++)
	if ( l < tip[t] )
	  {
	    for(w=MAX-1; w>t; --w) tip[w]=tip[w-1];
	    tip[t]=l;
	    break;
	  }
    }
  
  mtop=0;
  mtip=0;

  if ( MAX > 0 ) {
    for(t=0; t<MAX; t++) mtop+=top[t];
    mtop/=MAX;
    for(t=0; t<MAX; t++) mtip+=tip[t];
    mtip/=MAX;
  }

  if ( sim.nbMicrotub() )
    {
      ml     /= real(sim.nbMicrotub());
      mminus /= real(sim.nbMicrotub());
      mlsqr  /= real(sim.nbMicrotub());
      var     = mlsqr - sqr(ml);
      if ( var > 0 ) var = sqrt( var ); else var = 0;
    }

  MSG("F %3i %8.2fs ", sim.frameInBuffer(), sim.simTime() );
  MSG("%4i mt ( %4u %4u %4u %4u )G/S %5.2f +/- %.3f um",
	sim.nbMicrotub(), state[0], state[1], state[3], state[4], ml, var);

  MSG(" (%.2f %.2f)B/T%i ", mtip, mtop, MAX);
  MSG(" (%4i %4i)@m/M", min, max);

  if ( sim.simTime() != time_old )
    {
      real delay = sim.simTime() - time_old;
      time_old = sim.simTime();
      Vm = ( mminus - mminus_old ) / delay;
      Vp = ( ml - ml_old ) / delay - Vm;
      MSG("  %6.2f %6.2f  nm/s\n", 1000*Vm, 1000*Vp);
    }
  else
    MSG("\n");

  ml_old = ml;
  mminus_old = mminus;
}

//=====================================================================
//=====================================================================


void MesureMTDir()
{
  Vecteur a, b;
  real x, y, angle;
  int r;

  for ( Microtub *  mt = sim.firstMicrotub(); mt ; mt=mt->next() )
    for( r = 0; r < mt->nbSegments(); ++r)
      {
	b = a;
	y = x;

	a = mt->dpts( r );
	x = a.norm();

	if ( r > 0 ) 
	  {
	    angle = acos( a * b / ( x * y )  );
	    MSG("     %i : %.3f %.3f %.3f  : norm %10.6f  rod-rod angle %10.2f\n", 
		  r, a.XX, a.YY, a.ZZ, x, angle);
	  }
	else MSG("M%lx  0 : %.3f %.3f %.3f  : norm %10.6f\n", 
		   mt->getName(), a.XX, a.YY, a.ZZ, x);
		   
      }
}



void RodLength()
{
  int i, nbs[LMAX];
  real rl, scale = 0.5 / MP.mtrodlength;
  real dev=0, mx = 0, mn = 1000;

  for (i=0; i<LMAX; i++) nbs[i]=0;
  
  for ( Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() )
    {
      rl = mt->rodLength() / MP.mtrodlength;
      ++nbs[ int( 100 * rl * scale ) ];
      if ( rl > mx ) mx = rl;
      if ( rl < mn ) mn = rl;
      dev += sqr( rl - 1 );
    }
  dev /= sim.nbMicrotub();
  dev = sqrt( dev );

  MSG("mt->rodlength / MP.rodlength : min. %.6f  max. %.6f  std dev. %.6f\n",
	mn, mx, dev);
  MSG("rod L"); Report(nbs);
}



void RodLengthDeviation()
{
  int r;
  Vecteur a;
  real x, mx, mn, amx=0, amn=1000;
  
  for ( Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() )
    {
      mx = 0;
      mn = 1000;
      for(r=0; r<mt->nbSegments(); ++r)
	{
	  a = mt->dpts(r);
	  x = a.norm();
	  if ( x < mn ) mn = x;
	  if ( x > mx ) mx = x;
	}
      mn /= mt->rodLength();
      mx /= mt->rodLength();
      if ( mx > amx ) amx = mx;
      if ( mn < amn ) amn = mn;
    }
  MSG("pts-pts distance / mt->rodlength : min. %.6f  max. %.6f\n", amn, amx);
}

//=====================================================================
//=====================================================================

void MesureMTFold()
{
  Vecteur a, b;
  real x;
  int r, cnt=0;

  for ( Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() )
    {
      b = mt->dirP( 0 );
      for( r = 1; r < mt->nbSegments(); ++r )
	{
	  a = b;
	  b = mt->dirP( r );
	  x = a * b;
	  if ( x < 0 )
	    {
	      if ( x < -1 ) x = -1;
	      MSG("M%x  rod %i angle %.2f\n",
		    mt->getName(), r, 180*acos(x)/PI );
	      ++cnt;
	    }
	}
    }
  MSG("fold %i\n", cnt);
}


//=====================================================================
//=====================================================================


void MesureMTSpeed()
  //only for a small nb of microtubules, in one direction...
{
  Name n;
  static Name mtnames[50];
  static Vecteur vminus[50], vplus[50];
  Vecteur vm, vp;
  real vxm, vxp, vvm, vvp, dt;
  static real time=-999;

  if (time==-999)
    {
      n=0;
      for ( Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() )
	{
	  mtnames[n] = mt->getName();
	  vminus[n]  = mt->whereEnd(MT_MINUS_END);
	  vplus[n]   = mt->whereEnd(MT_PLUS_END);
	  if ( ++n>=10 ) break;
	}
      time=sim.simTime();
      return;
    }

  dt   = sim.simTime() - time;
  time = sim.simTime();

  for ( Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() )
    {
      for(n=0; (n<50) && (mtnames[n] != mt->getName()) ; ++n) ;
      if (n>=10) continue;

      vm        = -vminus[n];
      vminus[n] = mt->whereEnd(MT_MINUS_END);
      vm       += vminus[n];
      vxm       = vm[0]/dt;
      vvm       = vm.norm()/dt;
      

      vp        = -vplus[n];
      vplus[n]  = mt->whereEnd(MT_PLUS_END);
      vp       += vplus[n];
      vxp       = vp[0]/dt;
      vvp       = vp.norm()/dt;
      
      MSG("M%lx +/-ends : |v| = %6.3f %6.3f , vx =  %6.3f  %6.3f um/s\n",
	    mt->getName(), vvp, vvm, vxp, vxm);
    }
}


//=====================================================================
//=====================================================================

void AsterDistance()
{
  Vecteur x;

  Aster * as1 = sim.findAster(1);
  Aster * as2 = sim.findAster(2);

  if ( ( as1 == 0 ) || ( as2 == 0 )) return;

  x  = as1->getPosition() - as2->getPosition();

  MSG("%.2f s, asd= %6.3f um : %6.3f %6.3f %6.3f\n",
	sim.simTime(), x.norm(),
	x[0], x[1], x[2] );
}



void AsterStretch()
{
  Vecteur x, pos;
  real avg;
  Aster * as;
  int m, nbmts;

  for( int name=1; name <= 2; ++name )
    {
      as = sim.findAster(name);
      if ( as == 0 ) return;

      avg   = 0;
      pos   = as->getPosition();
      nbmts = as->nbMicrotub();

      for( m=0; m < nbmts; ++m )	{
        x = pos - as->getMicrotub(m)->whereEnd( as->getFocus() );
        avg += x.norm();
      }
      MSG(" as %i nbmts %i : %.2e um\n", name, nbmts, avg/nbmts );
    }
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


//=========================================================================


void Topology1(int type, int astersonly)
{
  Vecteur center;
  Name nn, i, j, k;
  int ty, endh;
  Name mt1, mt2;
  int a, nbpairs, nbasters, nbmts, nblinks, biggest[3];

  sim.renameMicrotub();
  Name NbMT = sim.nbMicrotub() + 1;
  int * asters = new int[NbMT];
  int * MTMT = new int[NbMT*NbMT];
  
  for( nn = 0; nn < NbMT*NbMT; ++nn)
    MTMT[ nn ]=0;
  
  for(Complex * mo=sim.firstBridgeComplex(); mo ; mo=mo->next() )
    {
      ty=mo->getType();
      if ( ty != type ) continue;
      endh = mo->getEnd1() + mo->getEnd2();
      if ( astersonly && ( endh != 2 ) ) continue;

      mt1 = mo->getMT1()->getName();
      mt2 = mo->getMT2()->getName();
      if (  ( mt1 >= NbMT ) || ( mt2 >= NbMT ) ) { MSG("d"); continue; }
      ++MTMT[ mt1 + NbMT * mt2 ];
      ++MTMT[ mt2 + NbMT * mt1 ];
    }

  //------------------------find the clusters of connected MTs:

  for( nn = 1; nn < NbMT; ++nn) asters[ nn ] = 0;
  for(  Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() )
    if ( mt->getName() < NbMT ) asters[ mt->getName() ] = mt->getName();

  do  {
    a=0;
    for( j=1;   j<NbMT; ++j) 
    for( k=j+1; k<NbMT; ++k)
      if ( ( MTMT[ j + NbMT * k ] ) && ( asters[ j ] != asters[ k ] ) )
	{
	  a = minT( asters[ j ], asters[ k ] );
	  asters[j] = a;
	  asters[k] = a;
	}
  } while (a);

  //------------------------print out the clusters:
  nbasters=0;
  nbpairs=0;
  for(i=0; i<3; ++i) biggest[i] = 0;

  for(i=1; i<NbMT; ++i)
    if ( ( a = asters[i] ) != 0 )
      {
	nbmts=0;
	nblinks=0;
	center=VZERO;
	for( j = i ; j < NbMT; ++j )
	  if ( asters[j] == a )
	    { 
	      ++nbmts; 
	      center += sim.findMicrotub(j)->whereEnd(MT_MINUS_END);
	      for(k=1; k<NbMT; ++k) nblinks += MTMT[j + NbMT * k];
	      asters[j]=0;
	    }
	nblinks /= 2;

	for( j=0; j<3; ++j ) if ( nbmts > biggest[j] )
	  {
	    for( k=2; k>j; --k) biggest[k] = biggest[k-1];
	    biggest[j] = nbmts;
	    break;
	  }


	if ( nbmts == 2 ) ++nbpairs;
	if ( nbmts > 2 ) 
	  {
	    ++nbasters;
	    center /= nbmts;
	    MSG("M%-4lx  nbmt %4i   links %6i    center %6.1f %6.1f\n", 
		  a, nbmts, nblinks, center.XX, center.YY);
	  }
      }
  if ( astersonly )
    MSG("nbasters ");
  else
    MSG("nbclusters ");

  MSG("%4i  nbpairs %4i   biggest %4i %4i %4i\n", 
	nbasters, nbpairs, biggest[0], biggest[1], biggest[2]);

  delete [] MTMT;
  delete [] asters;
}



//=====================================================================


void Topology2(int type)

  //counts every pair of filaments, in one of the three class X,T,V
  //the class is defined by the motors linking this pair, in order V, T, X.
  //i.e. all pairs with at least one V-motor will be counted as V
  //then pairs with T-motors will be counted as T,
  //and last the same with X...

{
  real MTpop[3];
  int MTLinks[3], MTLinksPop[3];
  int i, j, ty, endh;
  int mt1, mt2;

  if ( sim.nbComplexOfType(type) == 0 ) return;

  sim.renameMicrotub();
  int NbMT = sim.nbMicrotub() + 1;
  int * MTMT = new int[ 3*NbMT*NbMT + NbMT + 1 ];
  
  for(i=0; i<3*NbMT*NbMT+NbMT; ++i) MTMT[i] = 0;
  
  for(Complex * mo=sim.firstBridgeComplex(); mo ; mo=mo->next() )
    {
      ty=mo->getType();
      if ( ty != type ) continue;

      endh = mo->getEnd1() + mo->getEnd2();
      mt1 = mo->getMT1()->getName();
      mt2 = mo->getMT2()->getName();

      if (  (mt1>=NbMT) ||  (mt2>=NbMT) ) { MSG("d"); continue; }
      ++MTMT[ mt1 + NbMT * ( mt2 + NbMT * endh )];
      ++MTMT[ mt2 + NbMT * ( mt1 + NbMT * endh )];
    }

  for(endh=0; endh<3; ++endh) { MTLinks[endh]=0; MTLinksPop[endh]=0; }

  for(i=0; i<NbMT; ++i)
    for(j=0; j<NbMT; ++j)
      for(endh = 2; endh >= 0; --endh)
	if ( MTMT[i + NbMT * ( j + NbMT * endh )] )
	  { 
	    ++MTLinks[endh];
	    MTLinksPop[endh] += MTMT[ i + NbMT * ( j + NbMT * 0 ) ];
	    MTLinksPop[endh] += MTMT[ i + NbMT * ( j + NbMT * 1 ) ];
	    MTLinksPop[endh] += MTMT[ i + NbMT * ( j + NbMT * 2 ) ];
	    endh = -1;
	  } 

  for(endh=0; endh<3; ++endh)
    {  MTLinks[endh] /= 2; MTLinksPop[endh] /= 2; }

  for(i=0; i<3; ++i)
    if ( MTLinks[i] ) 
      MTpop[i] = MTLinksPop[i] / real(MTLinks[i]);
    else MTpop[i] = 0;

  MSG("%i  X %8i %4.1f  T %8i %4.1f  V %8i %4.1f ", type,
	MTLinks[0], MTpop[0],
	MTLinks[1], MTpop[1],
	MTLinks[2], MTpop[2]);
  MSG("\n");
  delete [] MTMT;
}

//=====================================================================
//=====================================================================

void EachMtHandsForce()
{
  int ty,i,ch[6];
  Vecteur v, vt[6];

  for ( Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() )
    {
      for (i=0; i<6; i++)  { vt[i]=VZERO; ch[i]=0; }

      for(Hand * h=mt->firstHand(); h ; h=h->next() )
	if (h->otherBound())
	  {
	    ty=h->getType();
	    ++ch[ty];
	    v=h->whereOtherSide()-h->whereHand();
	    Space::modulo(v);
	    vt[ty]+=v;
	  }

      for (i=0; i<6; i++)  vt[i]*=MP.km;

      MSG("M%lx :\n",mt->getName());
      MSG("Linking ha [ty]: %7d  %7d  %7d  %7d  %7d  %7d\n",
	    ch[0],ch[1],ch[2],ch[3],ch[4],ch[5]);
            
      MSG("VT. tension[ty]: %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f pN\n",
	    vt[0].XX,vt[1].XX,vt[2].XX,vt[3].XX,vt[4].XX,vt[5].XX);
    }
}


//=====================================================================
//=====================================================================

void CountHandsByType()
{
  int ty,bh[6],ch[6],th,tc;

  for (ty=0; ty<6; ty++)    { bh[ty]=0; ch[ty]=0; }
 
  for ( Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() )
    for(Hand * h=mt->firstHand(); h ; h=h->next() )
	{
	  ty=h->getType();
	  ++bh[ty];
	  if (h->otherBound())   ++ch[ty];
	}

  tc=0;
  th=0;
  for (ty=0; ty<6; ty++)    { th+=bh[ty];  tc+=ch[ty]; }

  MSG("bound hands  : total %5d, ", th);
  MSG("[ty]: %4d  %4d  %4d  %4d  %4d  %4d\n",
	 bh[0],bh[1],bh[2],bh[3],bh[4],bh[5]);

  MSG("linking hands: total %5d, ", tc);
  MSG("[ty]: %4d  %4d  %4d  %4d  %4d  %4d\n",
	 ch[0],ch[1],ch[2],ch[3],ch[4],ch[5]);
}

//============================================================================
//=====================================================================
//=====================================================================
void MesureHandsForce()
{
  int ty,i,ch[6];
  Vecteur v;
  real vn,t[6],max[6];
  int above[6];

  for (i=0; i<6; i++) 
    { ch[i]=0; t[i]=0; above[i]=0; max[i]=0; }
 
  for ( Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() )
    for( Hand * h=mt->firstHand(); h ; h=h->next() )
      if (h->otherBound())
	{
	  ty=h->getType();
	  ++ch[ty];
	  v=h->whereOtherSide()-h->whereHand();
	  Space::modulo(v);
	  vn = v.norm();
	  t[ty]+=vn;
	  if (vn > MP.habreakstretch[ty]) above[ty]++;
	  if (vn > max[ty]) max[ty] = vn;
	}

  ty=0;
  for (i=0; i<6; i++)
    {
      if (ch[i])
	{ ty=1; t[i] *= MP.km / (real)ch[i]; }
      max[i] *= MP.km;
    }

  if ( ty ) {
    MSG("Linking ha[ty]: %10d %10d %10d %10d %10d %10d\n",
	  ch[0], ch[1], ch[2], ch[3], ch[4], ch[5]);
    
    MSG("Avg. force[ty]: %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f pN\n",
	  t[0],t[1],t[2],t[3],t[4],t[5]);
    
    MSG("Max  force[ty]: %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f pN\n",
	  max[0],max[1],max[2],max[3],max[4],max[5]);
    
    MSG("Above Fmax[ty]: %10i %10i %10i %10i %10i %10i\n",
	  above[0],above[1],above[2],above[3],above[4],above[5]);
  }
}



//===========================================================================
//===========================================================================


void HandForce(real scale, int end=0)
{
  int ty,j;
  int vf, f[6][LMAX];
  Vecteur v;

  for (ty=0; ty<6; ty++) 
    for(j=0; j<LMAX; ++j)
      f[ty][j]=0;
 
  for ( Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() )
    for( Hand * h=mt->firstHand(); h ; h=h->next() )
      if (h->otherBound())
	if ( (end==0) || (h->getEnd() == end) )
	  {
	    ty=h->getType();
	    v=h->whereOtherSide()-h->whereHand();
	    Space::modulo(v);
	    vf=int(scale * v.norm() / MP.hamaxstretch[ty] );
	    f[ty][vf]++;
	  }

  MSG("h0"); Report(f[0]);
  MSG("h1"); Report(f[1]);
}

//-----------------------------------------------------------------------

void CountMultiLinks(int )
  //mesure the nb of links at end of MTs.
{
  int i, nbh[LMAX], cnt;
  int nbd[LMAX], div;
  const Microtub * other;
  Hand * hai, * haj, * ha;

  for (i=0; i<LMAX; i++) { nbh[i]=0; nbd[i]=0; }

  for ( Microtub * mti=sim.firstMicrotub(); mti ; mti=mti->next() )
    {
      div=0;
      for(hai=mti->firstHand(); hai ; hai=hai->next() )
	if ( ( other = hai->otherMicrotub() ) != 0 )
	    {
	      ++div;
	      cnt=1;
	      haj=hai;
	      haj=haj->next();
	      for(ha=haj++; ha ; ha=haj++ )
		if (other == ha->otherMicrotub() )
		  {
		    ++cnt;
		    ha->pop();
		    hai->push_after( ha );
		    hai=hai->next();
		  }
	      ++(nbh[cnt]);	
	    }
	++(nbd[div]);
      }

  MSG("ML  ");
  Report(nbh);
  MSG("div ");
  Report(nbd);
}


//=====================================================================
//=====================================================================

void ComplexSymmetry()
{
  int i, r, monb[2][4];

  for (i=0; i<4; i++) for(r=0; r<2; r++) 
    monb[r][i]=0; 

  for(Complex * mo=sim.firstBoundComplex(); mo ; mo=mo->next() )
    ++monb[ mo->isAttached1() ][ mo->getType() ];

  for (i=0; i<4; i++)
      if (sim.nbComplexOfType(i))
          MSG("cx%d %5d :  %4d %4d\n",
                i, sim.nbComplexOfType(i), monb[0][i], monb[1][i]);
}

//=====================================================================
//=====================================================================



void IntraMTLinks()
{
  int cnt = 0, near = 0,  dis, mindis = 10000;

  for(Complex * mo=sim.firstBridgeComplex(); mo ; mo=mo->next() )
    {
      const PointInterpolated & s1 = mo->getHand1();
      const PointInterpolated & s2 = mo->getHand2();
      if( s1.getPS() == s2.getPS() ) {
        ++cnt;
        dis = abs( s1.getPoint1() - s2.getPoint1() );
        if ( mindis > dis ) mindis = dis;
        near += ( dis < 2 );
      }
    }
  if ( cnt )
    MSG("intra mt = %d   min( rod1-rod2 ) = %d  #near %i\n", cnt, mindis, near);
}


//=====================================================================
//=====================================================================


void CountComplex(int sym = 0)
{
  Complex * mo;
  int chk, ty, state, monb[4][9];
  char stc[]="FBEbXTetV";

  for (ty=0; ty<4; ++ty)
    for(state=0; state<9; ++state)
      monb[ty][state]=0; 

  for(mo=sim.firstFreeComplex(); mo ; mo=mo->next() )
    ++monb[mo->getType()][mo->getState()];

  for(mo=sim.firstBoundComplex(); mo ; mo=mo->next() )
    ++monb[mo->getType()][mo->getState()];

  for(mo=sim.firstBridgeComplex(); mo ; mo=mo->next() )
    ++monb[mo->getType()][mo->getState()];

  if ( sym ) 
      for (ty=0; ty<4; ++ty)
      {
          monb[ty][1] += monb[ty][3];
          monb[ty][2] += monb[ty][6];
          monb[ty][5] += monb[ty][7];
      }
          
  for (ty=0; ty<4; ++ty) 
      if ( sim.nbComplexOfType(ty) )
      {
          chk=0;
          for(state=0; state<9; ++state)
              chk+=monb[ty][state];
          
          MSG("cx%d %5d :", ty, chk);
          for( state=0; state<9; ++state)
              if ( ( sym == 0 ) || isupper( stc[state] ) )
                  MSG("     %c %5d %2i%%", stc[state], 
                        monb[ty][state], 100*monb[ty][state]/chk);
          MSG("\n");
      }
}



//=====================================================================
//=====================================================================

void ComplexForceXTV()
{
  Complex * mo;
  Vecteur v;
  int ty, endh, i, j, cnt[6][3];
  real force[6][3], avgforce[6][3];

  for(i=0; i<6; i++) for(j=0; j<3; ++j)
    { cnt[i][j]=0; force[i][j]=0; }
 
  for(mo=sim.firstBridgeComplex(); mo ; mo=mo->next() )
    {
      ty=mo->getType();
      endh = mo->getEnd1() + mo->getEnd2();
      v = mo->getStretch();

      ++cnt[ty][endh];
      force[ty][endh] += v.norm();
    }

  for (i=0; i<4; i++) for(j=0; j<3; ++j)
    {
      force[i][j] *= MP.km;
      avgforce[i][j] = 0;
      if ( cnt[i][j] ) avgforce[i][j] = force[i][j] / real(cnt[i][j]);
    }

  for (i=0; i<4; i++) if ( cnt[i][0] + cnt[i][1] + cnt[i][2] > 0 )
    {
      MSG("complex[%i] nb. avg.F tot.F:", i);
      MSG("     X %6i %6.1f %6.1f", cnt[i][0], avgforce[i][0], force[i][0]);
      MSG("     T %6i %6.1f %6.1f", cnt[i][1], avgforce[i][1], force[i][1]);
      MSG("     V %6i %6.1f %6.1f", cnt[i][2], avgforce[i][2], force[i][2]);
      MSG(" pN\n");
    }
}

void ComplexForceDist(real scale)
{
  Complex * mo;
  int ty,j;
  int vf, pf[6][LMAX], af[6][LMAX];
  Vecteur v;

  for (ty=0; ty<6; ty++) 
    for(j=0; j<LMAX; ++j)
      { pf[ty][j]=0; af[ty][j]=0; }

  for(mo=sim.firstBridgeComplex(); mo ; mo=mo->next() )
    {
      ty=mo->getType();
      v = mo->getStretch();
      vf = int( scale * v.norm() / MP.hamaxstretch[mo->getType1()] );
      
      if ( mo->cosAngle() < 0)
	++af[ty][vf];
      else
	++pf[ty][vf];
    }
  MSG("anti 0"); Report(af[0]);
  MSG("para 0"); Report(pf[0]);
  MSG("anti 1"); Report(af[1]);
  MSG("para 1"); Report(pf[1]);
}

//=====================================================================
//=====================================================================

void MesureGrafteds()
{
  Grafted * gh;
  int ty,i,ghnb[2][6];
  Vecteur force[6];
  real stretch[6];
  Vecteur s;

  for (i=0; i<6; i++)
    {
      force[i]=VZERO;
      stretch[i]=0; 
      ghnb[0][i]=0; 
      ghnb[1][i]=0; 
    }

  for(gh=sim.firstFreeGrafted(); gh ; gh=gh->next() )
    ghnb[0][gh->getType()]++;

  for(gh=sim.firstBoundGrafted(); gh ; gh=gh->next() )
    {
      ty=gh->getType();
      ghnb[1][ty]++;
      s = gh->getStretch();
      force[ty] += s;
      stretch[ty] += s.norm();
    }

  for (i=0; i<6; i++)
    if ( ghnb[0][i] + ghnb[1][i] )
      {
	if  ( ghnb[0][i] + ghnb[1][i] != sim.nbGraftedOfType(i) )
	  MSG("***** problem while counting grafteds\n");

	MSG("%5d graft <%d>: %4d      %4d     ",
	       ghnb[0][i]+ghnb[1][i], i, ghnb[0][i], ghnb[1][i]);

	if ( ghnb[1][i] )
	  {
	    stretch[i] *= MP.km;
	    force[i] *= MP.km;
	    MSG("   To. %7.2f", stretch[i]);
	    MSG("   Av. %7.3f", stretch[i]/(real)ghnb[1][i]);
	    MSG("   VT. %7.2f", force[i].norm() );
	  }
	MSG(" pN\n");
      }
}

//=====================================================================
//=====================================================================
//=====================================================================
//=====================================================================

void showHelp()
{
  printf("nbrod      CountRods()\n");
  printf("rl         RodLength()\n");
  printf("rld        RodLengthDeviation()\n");
  printf("mtl        MesureMTlength()\n");
  printf("fold       MesureMTFold()\n");
  printf("mtd        MesureMTDir()\n");
  printf("eachmtl    EachMTlength()\n");
  printf("mtspeed    MesureMTSpeed()\n");

  printf("eachmtf    EachMtHandsForce()\n");
  printf("asd        AsterDistance()\n");
  printf("cluster    analysis of percolation\n");
  printf("topo       Topology2() for motor 0 and 1\n");
  printf("topo0      Topology() for motor 0\n");
  printf("topo1      Topology() for motor 1\n\n");
  printf("topo00     AsterTopology() for motor 0\n");
  printf("topo10     AsterTopology() for motor 1\n\n");
  
  printf("ha         CountHandsBygetType()\n");
  printf("hf         MesureHandsForce()\n");
  printf("f#         HandForce(#) # = 1, 5, 10, 50, 100, 1000\n");
  printf("f+,f-      HandForce(PLUS/MINUS END)\n\n");
  
  printf("mtha       CountHandsByMT()\n");
  printf("rodha      CountHandsByRod()\n");
  printf("rodha+     CountHandsByRod(MT_PLUS_END)\n");
  printf("rodha-     CountHandsByRod(MT_MINUS_END)\n");
  printf("sat        RodSaturation()\n\n");

  printf("ha+        CountHandsAtEnd(MT_PLUS_END)\n");
  printf("ha-        CountHandsAtEnd(MT_MINUS_END)\n\n");

  printf("lf         ComplexForceXTV()\n");
  printf("link       CountComplex()\n");
  printf("eachlink   EachMTLinks()\n");
  printf("multi      CountMultiLinks()\n");
  printf("multi-     CountMultiLink at minus end\n");
  printf("multi+     CountMultiLink at plus end\n\n");
  
  printf("cx         CountComplex()\n");
  printf("cxs        ComplexSymmetry()\n");
  printf("cxf        ComplexForceAvg()\n");
  printf("cxf#       ComplexForceDist(#) #=1, 5, 10, 50, 100\n");
  printf("gh         MesureGrafteds()\n\n");

  printf("intra      IntraMTLinks()\n");
  printf("pos        LinkPos (1 micron)\n");
  printf("pos#       LinkPos (# micron) #=05, 02, 01\n");
  printf("rpos       in-rod LinkPos (0.1 micron)\n");
  printf("rpos#      in-rod LinkPos (# micron) #=005, 002, 001\n\n");

  printf("asd        AsterDistance()\n");
  printf("ass        AsterStretch()\n");
  printf("asf        AsterForce() pull / push\n");
  printf("asf1       AsterForce() state\n");
  printf("asf2       AsterForce() para / anti with state\n");
  printf("asf3       AsterForce() pull / push with state\n");
  printf("asl        AsterForceAverage() para / anti\n");
}



int main(int argc, char * argv[])
{
  MSG.shutUp();
  MSG.setAbortOnError(0);
    
  if ( NO_ERROR != MP.parseFile(DATA_OUT) )
    if ( NO_ERROR != MP.parseFile(MP.datafile) )
      MSG.error("analyse.cc","Cannot read %s or data.out", MP.datafile);

  MP.massageParameters();

  printf("analyse1 DIM=%i, compiled %s, %s\n", DIM, __TIME__, __DATE__ );

  char code[4096]=" ";
  int  begin = 0, end = LMAX;

  if (argc<2) { showHelp(); return 0; } 

  MSG.setVerboseLevel(1);

  for(int i = 1; i < minT(argc, 8); ++i )
    if ( ( argv[i][0] == '-' ) && ( strlen( argv[i] ) > 1 ) ) 
      switch( argv[i][1] ) {
      case 'h': showHelp(); return 0;
      case 'v': MSG.setVerboseLevel( argv[i]+2 ); break;
      }
    else
      switch( sscanf(argv[i], "%d-%d", &begin, &end) ) 
        {
        case 0:
          code[0]=' ';
          strncat(code+1, argv[i], 127);
          strcat(code, " ");
          break;
        case 1:
          end = begin;
          break;
        }

  SimReader reader;
  if ( NO_ERROR != reader.openFile( RESULT_OUT ) )
    { printf("analyse1.cc: cannot open %s\n", RESULT_OUT); return 0; }

  for( int frame = begin; frame <= end ; ++frame )
    if ( NO_ERROR == reader.readFrame( frame ) ) {
      MSG("- - - - - - - - - - -F %i : %.3f s\n", sim.frameInBuffer(), sim.simTime());
      
      if ( strstr(code, " nbrod " ))    CountRods();
      if ( strstr(code, " rl " ))       RodLength();
      if ( strstr(code, " rld " ))      RodLengthDeviation();
      if ( strstr(code, " eachmtl " ))  EachMTlength();
      if ( strstr(code, " mtspeed " ))  MesureMTSpeed();
      if ( strstr(code, " mtl " ))      MesureMTLength();
      if ( strstr(code, " fold " ))     MesureMTFold();
      if ( strstr(code, " mtd " ))      MesureMTDir();
      if ( strstr(code, " eachmtf " ))  EachMtHandsForce();
      
      if ( strstr(code, " cx " ))       CountComplex();
      if ( strstr(code, " cx1 " ))      CountComplex(1);
      if ( strstr(code, " cxs " ))      ComplexSymmetry();
      if ( strstr(code, " cxf1 " ))     ComplexForceDist(1);
      if ( strstr(code, " cxf5 " ))     ComplexForceDist(5);
      if ( strstr(code, " cxf10 " ))    ComplexForceDist(10);
      if ( strstr(code, " cxf50 " ))    ComplexForceDist(50);
      if ( strstr(code, " cxf100 " ))   ComplexForceDist(100);
      if ( strstr(code, " gh " ))       MesureGrafteds();
      
      
      if ( strstr(code, " ha " ))       CountHandsByType();
      if ( strstr(code, " hf " ))       MesureHandsForce();
      if ( strstr(code, " f+ " ))       HandForce(5,MT_PLUS_END);
      if ( strstr(code, " f- " ))       HandForce(5,MT_MINUS_END);
      if ( strstr(code, " f1 " ))       HandForce(1);
      if ( strstr(code, " f5 " ))       HandForce(5);
      if ( strstr(code, " f10 " ))      HandForce(10);
      if ( strstr(code, " f50 " ))      HandForce(50);
      if ( strstr(code, " f100 " ))     HandForce(100);
      if ( strstr(code, " f1000 " ))    HandForce(1000);
      
      if ( strstr(code, " link " ))     CountComplex();
      if ( strstr(code, " lf " ))       ComplexForceXTV();
      if ( strstr(code, " multi " ))    CountMultiLinks(0);
      if ( strstr(code, " multi- " ))   CountMultiLinks(1);
      if ( strstr(code, " multi+ " ))   CountMultiLinks(2);
            
      if ( strstr(code, " cluster " ))  { Topology1(1, 0); }
      if ( strstr(code, " topo " ))     { Topology2(0); Topology2(1); }
      if ( strstr(code, " topo0 " ))    { Topology1(0, 1); Topology2(0); }
      if ( strstr(code, " topo1 " ))    { Topology1(1, 1); Topology2(1); }
      if ( strstr(code, " topo00 " ))   { Topology1(0, 1); }
      if ( strstr(code, " topo10 " ))   { Topology1(1, 1); }
      if ( strstr(code, " intra " ))    IntraMTLinks();
      
      
      if ( strstr(code, " asd " ))      AsterDistance();
    }
  return(0);
}
