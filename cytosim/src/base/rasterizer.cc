//RCS: $Id: rasterizer.cc,v 2.5 2005/02/07 13:50:44 nedelec Exp $

#include <cstdio>
#include <cstdlib>
#include "rasterizer.h"
#include "assert_macro.h"
#include "smath.h"

//keyword VISUAL_DEBUG is defined for compiling testpaint.cc
//it adds some calls to glVertex() in this file, for testing
#ifdef VISUAL_DEBUG
  #include "glut.h"
#endif

void Rasterizer::printPolygon( int nb, const real xy[] )
{
	for( int ii = 0; ii < nb; ++ii )
		printf("  %6.3f %6.3f", xy[2*ii], xy[2*ii+1] );
	printf("\n");
}

/*
 // Byte-wise swap two items of size SIZE.
 inline void swap(char * a, char * b, const size_t size)
 {
     register size_t ssize = size;                                        
     register char *aa = (a), *bb = (b);                                   
     do {                                                                     
         char tmp = *aa;                                                  
         *aa++ = *bb;                                                      
         *bb++ = tmp;                                                     
     } while (--ssize > 0);                                               
 }
 
 
 inline void swap( real xy[], const int i, const int j ) 
 {
     if ( i == j ) {
         printf("(%i %i)", i, j);
         return;
     }
     real tmp = xy[ 2*i ];
     xy[ 2*i ] = xy[ 2*j ];
     xy[ 2*j ] = tmp;
     tmp = xy[ 2*i + 1 ];
     xy[ 2*i + 1 ] = xy[ 2*j + 1 ];
     xy[ 2*j + 1 ] = tmp;
 }
 */

//------------------------------------------------------------------
void Rasterizer::convexHull2D( real xy[], int * nbpts, real dxy[] )
//change the count and reorders the points to yield the convex hull
{
    //we use three macro to be able to modify the functions while
    //keeping the code unchanged.
    
#define X(a) xy[(a)*2]
#define Y(a) xy[(a)*2+1]
#define swap(a,b) { \
    real tmp=xy[(a)*2]; xy[(a)*2]  =xy[(b)*2];     xy[(b)*2]=tmp; \
    tmp=xy[(a)*2+1];    xy[(a)*2+1]=xy[(b)*2+1];   xy[(b)*2+1]=tmp; \
    if ( dxy ) { \
    tmp=dxy[(a)*2];    dxy[(a)*2]  =dxy[(b)*2];   dxy[(b)*2]=tmp; \
    tmp=dxy[(a)*2+1];  dxy[(a)*2+1]=dxy[(b)*2+1]; dxy[(b)*2+1]=tmp; }}
        
    if ( *nbpts < 3 ) return;
        
    //---------- find the bottom-left and up-right points:
    int ii, right = 0, next = 0;
    for( ii = 1; ii < *nbpts; ++ii )
        {
        if (( X(ii) < X(next) ) 
            || (( X(ii) == X(next) ) && ( Y(ii) < Y(next) )))
            next = ii;
        
        if (( X(ii) > X(right) )
            || (( X(ii) == X(right) ) && ( Y(ii) > Y(right) )))
            right = ii;
        }
    
    if ( next == right ) 
        { *nbpts = 1; return; }  //all points are equal ?
    
    //put the leftmost point at position zero:
    if ( next )
        swap( 0, next );
    
    //put the rightmost point last:
    if ( right )
        { swap( *nbpts-1, right); }
    else  //since we swaped with zero, what was in zero is now in next!
        swap( *nbpts-1, next );
        
    //printPolygon(*nbpts, xy);
    
    int curr = 0;
    right = *nbpts-1;
    real dx, dy, dxt, dyt;
    // First sweep from left to right under the convex hull:
    do {
        
        dy = Y(right) - Y(curr);
        dx = X(right) - X(curr);
        next = right;
        for( ii = right-1; ii > curr; --ii ) {
            dyt = Y(ii) - Y(curr);
            dxt = X(ii) - X(curr);
            if ( dyt * dx < dy * dxt )
                { next = ii; dx = dxt; dy = dyt; }
        }
        if ( ++curr >= *nbpts-1 ) 
            return;
        if ( curr != next )
            swap( curr, next );
        //printf(" cnt %i ", curr );  printPolygon(*nbpts, xy);
        
    } while( next != right );
    
    //printf("sweep curr=%i\n", curr);
    // Sweep from right to left above to finish the convex hull:
    do {
        
        next = 0;
        dy = Y(0) - Y(curr); 
        dx = X(0) - X(curr);
        for( ii = curr+1; ii < *nbpts; ++ii ) {
            dyt = Y(ii) - Y(curr); 
            dxt = X(ii) - X(curr);
            if ( dyt * dx < dy * dxt ) 
                { next = ii; dx = dxt; dy = dyt; }
        }
        ++curr;
        if ( next == 0 ) 
            { *nbpts = curr; return; }
        if ( curr != next ) 
            swap( curr, next );
        //printf(" cnt %i ", curr );  print(*nbpts, xy);
        
    } while ( X(0) < X(curr) );
    
    *nbpts = curr+1;
    
#undef X
#undef Y
#undef swap
}


//=============================================================================
//                             1D
//=============================================================================


void Rasterizer::paintPolygon1D(void paint(const int, const int, const int), 
                                const int nbpts, const real x[],
                                const int yy, const int zz)
{
    
    //---------- find the extreme points:
    int left  =  (int)ceil( x[0] );
    int right = (int)floor( x[0] );
    
    int xii;
    for( int ii = 1; ii < nbpts; ++ii )
        {
        xii = (int)ceil( x[ ii ] );
        if ( xii < left ) left = xii;
        
        xii = (int)floor( x[ ii ] );
        if ( xii > right ) right = xii;
        }
    
    for( ; left <= right; ++left )
        paint( left, yy, zz );
}


//------------------------------------------------------------------
void Rasterizer::paintFatLine1D(void paint(const int, const int, const int), 
                                const real p[], const real q[], const real width, 
                                const real min[], const real delta[])
{
    real pq;
    if ( q[0] > p[0] ) pq = width; else pq = -width;
    real pts[2];
    pts[ 0 ] = ( p[ 0 ] - pq - min[ 0 ] ) / delta[ 0 ];
    pts[ 1 ] = ( q[ 0 ] + pq - min[ 0 ] ) / delta[ 0 ];
    paintPolygon1D( paint, 2, pts );
}

//=============================================================================
//                               2D
//=============================================================================


//call paintPoint(x,y,zz) for every point x,y inside the polygon
//of nbpts points who's vertex are given in xy[]. The polygon
//should be convex, and points should be ordered in a consistent
//manner (clockwise or counter-CW).

//this version can miss some points on the edge of the polygon due to round-off
//errors, especially when some of the vertices have integral values.
void Rasterizer::paintPolygon2D(void paint(const int, const int, const int), 
                                const int nbpts, const real xy[], 
                                const int zz)
{
    //printf("P ");for(int i=0; i<2*nbpts; ++i) printf("%5.2f ",xy[i]);
    //printf("\n");
    
    //---------- find the bottom point:
    int cnt = 0, yy, indx1 = 1;
    
    for( yy = 2*nbpts-1; yy > 1; yy -= 2 )
        if ( xy[ yy ] < xy[ indx1 ] )
            indx1 = yy;
    
    //printf("bottom %i at %8.2f\n", indx1, xy[indx1] );
    
    yy = (int)ceil( xy[ indx1 ] );
    int indx2 = indx1, xx, xxmx;
    int yymx1 = yy-1, yymx2 = yymx1, yymx;
    real xx1 = 0, yy1 = 0, dxx1 = 0, yyn;
    real xx2 = 0, yy2, dxx2 = 0;
    
    while( true ) {
        //printf("find lines at yy %i :\n", yy );
        
        //---------find the two lines starting from there:
        if (( yymx1 < yy ) && ( xy[ indx1 ] <= yy )) {
            do {
                if ( ++cnt > nbpts ) return;
                xx1   = xy[ indx1 - 1 ];
                yy1   = xy[ indx1 ];
                
                indx1 += 2;
                if ( indx1 >= 2*nbpts ) indx1 = 1;
                
                yyn   = xy[ indx1 ];
            } while ( yyn <= yy );
            yymx1 = (int)floor( yyn );
            dxx1 = ( xy[ indx1 - 1 ] - xx1 ) / ( yyn - yy1 );
            xx1 += dxx1 * ( yy - yy1 );
        }
        
        if (( yymx2 < yy ) && ( xy[ indx2 ] <= yy )) {
            do {
                if ( ++cnt > nbpts ) return;
                xx2   = xy[ indx2 - 1 ];
                yy2   = xy[ indx2 ];
                
                indx2 -= 2;
                if ( indx2 < 0 ) indx2 = 2*nbpts-1;
                
                yyn   = xy[ indx2 ];
            } while ( yyn <= yy );
            yymx2 = (int)floor( yyn );
            dxx2 = ( xy[ indx2 - 1 ] - xx2 ) / ( yyn - yy2 );
            xx2 += dxx2 * ( yy - yy2 );
        }
        
        //----------draw the horizontal line at yy, between xx1 and xx2:
        if ( yymx1 < yymx2 ) yymx = yymx1; else yymx = yymx2;      
        
        if (( xx1 < xx2 ) || (( xx1 == xx2 ) && ( dxx1 < dxx2 )))
            for( ; yy <= yymx; ++yy ) {
                xx   = (int)ceil(xx1);
                xxmx = (int)floor(xx2);
                
                //printf("yy %i : xx %5.2f %5.2f  %i -> %i\n", yy, xx1, xx2, xx, xxmx);
                for( ; xx <= xxmx; ++xx ) 
                    paint( xx, yy, zz );
                
                xx1 += dxx1;
                xx2 += dxx2;
            }
                else
                    for( ; yy <= yymx; ++yy ) {
                        xx   = (int)ceil(xx2);
                        xxmx = (int)floor(xx1);
                        
                        //printf("yy %i : xx %5.2f %5.2f  %i -> %i\n", yy, xx1, xx2, xx, xxmx);
                        for( ; xx <= xxmx; ++xx ) 
                            paint( xx, yy, zz );
                        
                        xx1 += dxx1;
                        xx2 += dxx2;
                    }
    }
}


//------------------------------------------------------------------
void Rasterizer::paintFatLine2D(void paint(const int, const int, const int), 
                                const real p[], const real q[], const real width)
{
    static real pts[ 8 ];
    
    real dx = q[ 0 ] - p[ 0 ];
    real dy = q[ 1 ] - p[ 1 ];
    
    real s = 0.5 * width / sqrt( dx * dx + dy * dy );
    
    dx *= s;
    dy *= s;
    
    pts[ 0 ] = p[ 0 ] + dy;
    pts[ 1 ] = p[ 1 ] - dx;
    
    pts[ 2 ] = p[ 0 ] - dy;
    pts[ 3 ] = p[ 1 ] + dx;
    
    pts[ 4 ] = q[ 0 ] - dy;
    pts[ 5 ] = q[ 1 ] + dx;
    
    pts[ 6 ] = q[ 0 ] + dy;
    pts[ 7 ] = q[ 1 ] - dx;
    
#ifdef VISUAL_DEBUG
    glColor3f(1.0, 1.0, 0.0);
    glBegin( GL_LINE_LOOP );
    for(int ii=0; ii<=6; ii+=2 )
        glVertex2f( pts[ii], pts[ii+1] );
    glColor3f(1.0, 1.0, 1.0);
    glEnd();
#endif
    
    paintPolygon2D( paint, 4, pts );
}



//------------------------------------------------------------------
void Rasterizer::paintFatLine2D(void paint(const int, const int, const int), 
                                const real p[], const real q[], const real width,
                                const real min[], const real delta[],
                                real scaling )
{
    real pqx = q[ 0 ] - p[ 0 ];
    real pqy = q[ 1 ] - p[ 1 ];
    
    //scaling is known beforhand for microtubules, since all rods have the 
    //same length. The optional argument enables a small optimisation here
    if ( scaling == 0 )
        scaling = width / sqrt( pqx * pqx + pqy * pqy );
    
    pqx *= scaling / delta[ 0 ];
    pqy *= scaling / delta[ 1 ];
    
    real endPx, endPy, endQx, endQy,  pts[8];
    
    endPx = ( p[ 0 ] - min[ 0 ] ) / delta[ 0 ] - pqx;
    endPy = ( p[ 1 ] - min[ 1 ] ) / delta[ 1 ] - pqy;
    
    endQx = ( q[ 0 ] - min[ 0 ] ) / delta[ 0 ] + pqx;
    endQy = ( q[ 1 ] - min[ 1 ] ) / delta[ 1 ] + pqy;
    
    pts[ 0 ] = endPx + pqy;
    pts[ 1 ] = endPy - pqx;
    
    pts[ 2 ] = endPx - pqy;
    pts[ 3 ] = endPy + pqx;
    
    pts[ 4 ] = endQx - pqy;
    pts[ 5 ] = endQy + pqx;
    
    pts[ 6 ] = endQx + pqy;
    pts[ 7 ] = endQy - pqx;
    
#ifdef VISUAL_DEBUG
    glColor3f(1.0, 1.0, 0.0);
    glBegin( GL_LINE_LOOP );
    for(int ii=0; ii<=6; ii+=2 )
        glVertex2f( pts[ii]*delta[0]+min[0], pts[ii+1]*delta[1]+min[1] );
    glColor3f(1.0, 1.0, 1.0);
    glEnd();
#endif  
    
    paintPolygon2D( paint, 4, pts );
}


//=============================================================================
//                               3D
//=============================================================================


// function for qsort: compares the Z component of the two points
int Rasterizer::comp_higher( const void * a, const void * b)
{
  if (((real*)a)[2] > ((real*)b)[2]) return  1;
  if (((real*)a)[2] < ((real*)b)[2]) return -1;
  return 0;
}


//------------------------------------------------------------------
void Rasterizer::paintPolygon3D(void paint(const int, const int, const int), 
                                const int nbpts, real xyz[])
// the polygon is the convex hull of the given points
// algorithm: we section at each integral Z, collect the intersection of
// all possible lines connecting two points, 
// and call paintPolygon2D with the convex hull of all these points
// there is a better version of this function below.
{
    if ( nbpts < 1 ) return; //need at least one point
                             // first order the points in increasing Z:
  qsort( xyz, nbpts, 3*sizeof(real), &Rasterizer::comp_higher );
    
    static int allocated_size = 0;
    static real * xy = 0, *cxy = 0, *dxy = 0;
    
    if ( allocated_size < 2 * nbpts * nbpts ) {
        allocated_size = 2 * nbpts * nbpts;
        if ( xy ) { delete[] xy; delete[] cxy; delete[] dxy; }
        xy  = new real[ allocated_size ];
        cxy = new real[ allocated_size ];
        dxy = new real[ allocated_size ];
        
        if (( xy == 0 ) || ( cxy == 0 ) || ( dxy == 0 )) {
          fprintf(stderr, "Rasterizer::paintPolygon3D():: memory allocation failed\n");
          exit(1);
        }
    }
    
    real dzz;
    int ii, jj, above = 0, nbl, nbp;
    
    int zz  = (int) ceil( xyz[ 2 ] );
    int zzn;
    
    while( true ) {
        
        //printf("restart at zz %4i\n", zz );
        
        //we move to the next point, ensuring that we will exit for sure...
        if ( ++above >= nbpts ) return;
        
        //find the first point strictly above the plane Z = zz:
        //the index of this point is (above-1)
        while( xyz[ 3*above + 2 ] <= zz )
            if ( ++above >= nbpts ) return;
        
        //the next time we have to recalculate the lines
        //is when this point will be below the plane Z = zzn:
        zzn = (int)ceil( xyz[ 3 * above + 2 ] );
        
        //set-up all the lines, which join any point below the plane 
        //to any point above the plane:
        
        nbl = 0;
        
        for( ii = 0; ii < above; ++ii )
            for( jj = above; jj < nbpts; ++jj ) 
                {
                xy[ 2 * nbl     ]  = xyz[ 3 * ii     ];
                xy[ 2 * nbl + 1 ]  = xyz[ 3 * ii + 1 ];
                dzz = xyz[ 3 * jj + 2 ] - xyz[ 3 * ii + 2 ];
                assert( dzz > 0 );
                if ( dzz > 0 ) {
                    dxy[ 2 * nbl     ] = ( xyz[ 3 * jj     ] - xyz[ 3 * ii     ] ) / dzz;
                    dxy[ 2 * nbl + 1 ] = ( xyz[ 3 * jj + 1 ] - xyz[ 3 * ii + 1 ] ) / dzz;
                    dzz = zz - xyz[ 3 * ii + 2 ];
                    xy[ 2 * nbl     ] += dxy[ 2 * nbl     ] * dzz;
                    xy[ 2 * nbl + 1 ] += dxy[ 2 * nbl + 1 ] * dzz;
                    ++nbl;
                }
                }
                
                assert( nbl < nbpts * nbpts );
        //printf("Z %4i nxt Z %4i split point above %2i #lines %4i\n",
        //	   zz, zzn, above, nbl);
        
        for( ; zz < zzn; ++zz ) {
            
            for( ii = 0; ii < 2 * nbl; ++ii ) {
                cxy[ ii ] = xy[ ii ];
                xy[ ii ] += dxy[ ii ];
            }
            nbp = nbl;
            convexHull2D( cxy, &nbp );
            
            //printf("zz %3i : nbp = %i\n", zz, nbp);
            
            paintPolygon2D( paint, nbp, cxy, zz );
            
#ifdef VISUAL_DEBUG
            glColor3f(0,0,1.0);
            glBegin(GL_LINE_LOOP);
            for( ii = 0; ii < nbp; ++ii )
                glVertex3f( cxy[2*ii], cxy[2*ii+1], zz );
            glEnd();
#endif
        }
    }
}



//------------------------------------------------------------------
void Rasterizer::paintFatLine3D_old(void paint(const int, const int, const int), 
                                    const real p[], const real q[], const real width,
                                    const real min[], const real delta[], real scaling )
{
    real n, pq[3];
    
    pq[0] = q[0] - p[0];
    pq[1] = q[1] - p[1];
    pq[2] = q[2] - p[2];
    
    //we normalize pq to width  
    //scaling is known beforhand for microtubules, since all rods have the 
    //same length. The optional argument enables a small optimisation here
    if ( scaling == 0 )
        scaling = width / sqrt( pq[0]*pq[0] + pq[1]*pq[1] + pq[2]*pq[2] );
    
    pq[0] *= scaling;
    pq[1] *= scaling;
    pq[2] *= scaling;
    
    // we will set-up two vectors A and B perpendicular to PQ,
    // we first set a perpendicular to PQ and one of the coordiate axis
    // discarding the dimension where PQ vector is smallest:
    real a[3], b[3];
    
    if ( fabs( pq[0] ) < fabs( pq[1] ) ) {
        
        if ( fabs( pq[0] ) < fabs( pq[2] ) )
            { a[0] =  0.0;   a[1] = -pq[2]; a[2] =  pq[1]; } //pq[0] is the smallest
        else
            { a[0] = -pq[1]; a[1] =  pq[0]; a[2] =  0.0;   } //pq[2] is the smallest
        
    } else {
        
        if ( fabs( pq[1] ) < fabs( pq[2] ) )
            { a[0] = -pq[2]; a[1] =  0.0;   a[2] =  pq[0]; } //pq[1] is the smallest
        else
            { a[0] = -pq[1]; a[1] =  pq[0]; a[2] =  0.0;   } //pq[2] is the smallest
        
    }
    
    // vector B is set perpendicular to PQ and A by cross product
    b[0] = pq[1] * a[2] - pq[2] * a[1];
    b[1] = pq[2] * a[0] - pq[0] * a[2];
    b[2] = pq[0] * a[1] - pq[1] * a[0];
    
    //vectors A & B are now perpendicular to PQ, we normalize them:
    n = width / sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] );
    a[0] *= n / delta[0];
    a[1] *= n / delta[1];
    a[2] *= n / delta[2];
    
    n = width / sqrt( b[0]*b[0] + b[1]*b[1] + b[2]*b[2] );
    b[0] *= n / delta[0];
    b[1] *= n / delta[1];
    b[2] *= n / delta[2];
    
    //we set up two points on the line PQ by extending the
    //segment PQ by width on each side,
    //at the same time, we convert to the grid's coordinates:
    // grid coordinates = ( coordinates - min ) / delta
    real ends[6];
    
    ends[ 0 ] = ( p[ 0 ] - pq[0] - min[ 0 ] ) / delta[ 0 ];
    ends[ 1 ] = ( p[ 1 ] - pq[1] - min[ 1 ] ) / delta[ 1 ];
    ends[ 2 ] = ( p[ 2 ] - pq[2] - min[ 2 ] ) / delta[ 2 ];
    
    ends[ 3 ] = ( q[ 0 ] + pq[0] - min[ 0 ] ) / delta[ 0 ];
    ends[ 4 ] = ( q[ 1 ] + pq[1] - min[ 1 ] ) / delta[ 1 ];
    ends[ 5 ] = ( q[ 2 ] + pq[2] - min[ 2 ] ) / delta[ 2 ];
    
    //last, we set 8 corners of a square volume aligned along PQ
    real pts[24];
    
    pts[ 0  ] = ends[ 0 ] + a[0];
    pts[ 1  ] = ends[ 1 ] + a[1];
    pts[ 2  ] = ends[ 2 ] + a[2];
    
    pts[ 3  ] = ends[ 0 ] + b[0];
    pts[ 4  ] = ends[ 1 ] + b[1];
    pts[ 5  ] = ends[ 2 ] + b[2];
    
    pts[ 6  ] = ends[ 0 ] - a[0];
    pts[ 7  ] = ends[ 1 ] - a[1];
    pts[ 8  ] = ends[ 2 ] - a[2];
    
    pts[ 9  ] = ends[ 0 ] - b[0];
    pts[ 10 ] = ends[ 1 ] - b[1];
    pts[ 11 ] = ends[ 2 ] - b[2];
    
    pts[ 12 ] = ends[ 3 ] + a[0];
    pts[ 13 ] = ends[ 4 ] + a[1];
    pts[ 14 ] = ends[ 5 ] + a[2];
    
    pts[ 15 ] = ends[ 3 ] + b[0];
    pts[ 16 ] = ends[ 4 ] + b[1];
    pts[ 17 ] = ends[ 5 ] + b[2];
    
    pts[ 18 ] = ends[ 3 ] - a[0];
    pts[ 19 ] = ends[ 4 ] - a[1];
    pts[ 20 ] = ends[ 5 ] - a[2];
    
    pts[ 21 ] = ends[ 3 ] - b[0];
    pts[ 22 ] = ends[ 4 ] - b[1];
    pts[ 23 ] = ends[ 5 ] - b[2];
    
    paintPolygon3D( paint, 8, pts );
    
#ifdef VISUAL_DEBUG
    glColor3f(1.0, 1.0, 0.0);
    for(int ii=0; ii<=21; ii+=3 )
        glVertex3f( pts[ii], pts[ii+1], pts[ii+2] );
    glColor3f(1.0, 1.0, 1.0);
#endif
}


//=============================================================================
//=============================================================================



void Rasterizer::paintPolygon3D(void paint(const int, const int, const int), 
                                const int nbpts, real xyzi[], 
                                const char sides[])
// the polygon is the convex hull of the 'nbpts' points given in xyzi[],
// containing the three coordinate + the id of the point, a short int
// cast into a 'real' in position 4.
// sides[] gives information on which pair of points should be considered:
// sides[ ii + nbpts * jj ] == 1 if points of index ii and jj are connected,
// i.e. they are one of the edge of the 3D polygon 
//
// algorithm: we section at each integer Z, collect the intersection of
// edges and call paintPolygon2D with the convex hull of all these points
{
    if ( nbpts < 1 ) return; //need at least one point
                             // first order the points in increasing Z:
    qsort( xyzi, nbpts, 4*sizeof(real), &Rasterizer::comp_higher );
    
    //we can normally only cross four sides of a parallelogram in 3D
    //but in some degenerate cases, it can be more
    const int allocated = 8;
    real  xy[2*allocated];
    real dxy[2*allocated];
    
    real dzz;
    int ii, jj, above = 0;
    
    int zz  = (int) ceil( xyzi[ 2 ] );
    
    while( true ) {
        
        //printf("restart at zz %4i\n", zz );
        
        //we move to the next point, ensuring that we will exit for sure...
        if ( ++above >= nbpts ) return;
        
        //find the first point strictly above the plane Z = zz:
        //the index of this point is (above-1)
        while( xyzi[ 4*above + 2 ] <= zz )
            if ( ++above >= nbpts ) return;
        
        //the next time we have to recalculate the lines
        //is when this point will be below the plane Z = zzn:
        int zzn = (int)ceil( xyzi[ 4 * above + 2 ] );
        
        //set-up all the lines, which join any point below the plane
        //to any point above the plane, being a edge of the solid polygon:
        
        int nbl = 0;  //number of edges crossing the plane at Z=zz;
        for( ii = 0; ii < above; ++ii )
            {
            
            short iik = *(short*)(xyzi+4*ii+3);       //the point index in the original array
            
            for( jj = above; jj < nbpts; ++jj ) 
                {
                
                short jjk = *(short*)(xyzi+4*jj+3);   //the point index in the original array
                                                      //printf("%i %i -> %i\n", iik, jjk, sides[iik+nbpts*jjk]);
                
                if ( sides[ iik + nbpts * jjk ] ) {   //if ii and jj are joined by a line:
                    
                    xy[ 2 * nbl     ]  = xyzi[ 4 * ii     ];
                    xy[ 2 * nbl + 1 ]  = xyzi[ 4 * ii + 1 ];
                    
                    dzz = xyzi[ 4 * jj + 2 ] - xyzi[ 4 * ii + 2 ];
                    assert( dzz > 0 );
                    
                    dxy[ 2 * nbl     ] = ( xyzi[ 4 * jj     ] - xyzi[ 4 * ii     ] ) / dzz;
                    dxy[ 2 * nbl + 1 ] = ( xyzi[ 4 * jj + 1 ] - xyzi[ 4 * ii + 1 ] ) / dzz;
                    
                    dzz = zz - xyzi[ 4 * ii + 2 ];
                    
                    xy[ 2 * nbl     ] += dxy[ 2 * nbl     ] * dzz;
                    xy[ 2 * nbl + 1 ] += dxy[ 2 * nbl + 1 ] * dzz;
                    
                    ++nbl;
                }
                }
            }
        
        assert( nbl < allocated );
        //printf("Z %4i nxt Z %4i split point above %2i #lines %4i\n",
        //	   zz, zzn, above, nbl);
        
        // the edges of the convex solid polygon should not intersect,
        // so we can take the convex hull only once here:
        int nbp; //number of points after the hull.
        int require_hull = 1;
        
        for( ; zz < zzn; ++zz ) {
            
            if ( require_hull ) {
                //try to make the convex hull of the points from xy[]:
                nbp = nbl;
                convexHull2D( xy, &nbp, dxy );   
                //in the particular case where some points overlap, we might
                //loose them, in which case we need to redo the hull later
                require_hull = ( nbp != nbl );
            }
            
            paintPolygon2D( paint, nbp, xy, zz );
            
#ifdef VISUAL_DEBUG
            glColor3f(0,0,1.0);
            glBegin(GL_LINE_LOOP);
            for( ii = 0; ii < nbp; ++ii )
                glVertex3f( xy[2*ii  ], xy[2*ii+1], zz );
            glEnd();
#endif           
            //update the coordinates according to the slopes, for the next zz:
            for( ii = 0; ii < 2 * nbl; ++ii )
                xy[ ii ] += dxy[ ii ];
        }
    }
}


//------------------------------------------------------------------
void Rasterizer::paintFatLine3D(void paint(const int, const int, const int),
                                const real p[], const real q[], const real width,
                                const real min[], const real delta[], real scaling )
{
    real pqx = q[0] - p[0];
    real pqy = q[1] - p[1];
    real pqz = q[2] - p[2];
    
    //we normalize pq to width
    //scaling is known beforehand for microtubules, since all rods have the
    //same length. The optional argument enables a small optimisation here
    if ( scaling == 0 )
        scaling = width / sqrt( pqx*pqx + pqy*pqy + pqz*pqz );
    
    pqx *= scaling;
    pqy *= scaling;
    pqz *= scaling;
    
    // we will set-up two vectors A and B perpendicular to PQ,
    // we first set A perpendicular to PQ and one of the coordiate axis
    // discarding the dimension where PQ vector is smallest:
    real ax, ay, az;
    
    if ( fabs( pqx ) < fabs( pqy ) ) {
        
        if ( fabs( pqx ) < fabs( pqz ) )
            { ax =  0.0;  ay = -pqz;  az =  pqy; } //pqx is the smallest
        else
            { ax = -pqy;  ay =  pqx;  az =  0.0; } //pqz is the smallest
        
    } else {
        
        if ( fabs( pqy ) < fabs( pqz ) )
            { ax = -pqz;  ay =  0.0;  az =  pqx; } //pqy is the smallest
        else
            { ax = -pqy;  ay =  pqx;  az =  0.0; } //pqz is the smallest
        
    }
    
    // vector B is set perpendicular to PQ and A by cross product
    // we could save two multiplications here, since one of a[] is null
    real bx = pqy * az - pqz * ay;
    real by = pqz * ax - pqx * az;
    real bz = pqx * ay - pqy * ax;
    
    //vectors A & B are now perpendicular to PQ, we normalize them:
    //we could save some multiplications here, since one of a[] is null
    real n = width / sqrt( ax*ax + ay*ay + az*az );
    ax *= n / delta[0];
    ay *= n / delta[1];
    az *= n / delta[2];
    
    n = width / sqrt( bx*bx + by*by + bz*bz );
    bx *= n / delta[0];
    by *= n / delta[1];
    bz *= n / delta[2];
    
    //we set up two points on the line PQ by extending the
    //segment PQ by width on each side,
    //at the same time, we convert to the grid's coordinates:
    // grid coordinates = ( coordinates - min ) / delta
    
    real endPx = ( p[ 0 ] - pqx - min[ 0 ] ) / delta[ 0 ];
    real endPy = ( p[ 1 ] - pqy - min[ 1 ] ) / delta[ 1 ];
    real endPz = ( p[ 2 ] - pqz - min[ 2 ] ) / delta[ 2 ];
    
    real endQx = ( q[ 0 ] + pqx - min[ 0 ] ) / delta[ 0 ];
    real endQy = ( q[ 1 ] + pqy - min[ 1 ] ) / delta[ 1 ];
    real endQz = ( q[ 2 ] + pqz - min[ 2 ] ) / delta[ 2 ];
    
    //last, we set 8 corners of a square volume aligned along PQ
    real pts[32];
    
    pts[ 0  ] = endPx + ax;
    pts[ 1  ] = endPy + ay;
    pts[ 2  ] = endPz + az;
    *((short*)(pts+3)) = 0;
    
    pts[ 4  ] = endPx + bx;
    pts[ 5  ] = endPy + by;
    pts[ 6  ] = endPz + bz;
    *((short*)(pts+7)) = 1;
    
    pts[ 8  ] = endPx - ax;
    pts[ 9  ] = endPy - ay;
    pts[ 10 ] = endPz - az;
    *((short*)(pts+11)) = 2;
    
    pts[ 12 ] = endPx - bx;
    pts[ 13 ] = endPy - by;
    pts[ 14 ] = endPz - bz;
    *((short*)(pts+15)) = 3;
    
    pts[ 16 ] = endQx + ax;
    pts[ 17 ] = endQy + ay;
    pts[ 18 ] = endQz + az;
    *((short*)(pts+19)) = 4;
    
    pts[ 20 ] = endQx + bx;
    pts[ 21 ] = endQy + by;
    pts[ 22 ] = endQz + bz;
    *((short*)(pts+23)) = 5;
    
    pts[ 24 ] = endQx - ax;
    pts[ 25 ] = endQy - ay;
    pts[ 26 ] = endQz - az;
    *((short*)(pts+27)) = 6;
    
    pts[ 28 ] = endQx - bx;
    pts[ 29 ] = endQy - by;
    pts[ 30 ] = endQz - bz;
    *((short*)(pts+31)) = 7;
    
    //the static matrix of vertices connections:
    //it defines which point are connected to form
    //the edges of the polygon.
    static const char sides[8*8] = {
        0,1,0,1,1,0,0,0,
        1,0,1,0,0,1,0,0,
        0,1,0,1,0,0,1,0,
        1,0,1,0,0,0,0,1,
        1,0,0,0,0,1,0,1,
        0,1,0,0,1,0,1,0,
        0,0,1,0,0,1,0,1,
        0,0,0,1,1,0,1,0 };
    
    paintPolygon3D( paint, 8, pts, sides );
    
#ifdef VISUAL_DEBUG
    glColor3f(1.0, 1.0, 0.0);
    for(int ii=0; ii<=21; ii+=3 )
        glVertex3f( pts[ii  ]*delta[0]+min[0],
                    pts[ii+1]*delta[1]+min[1], 
                    pts[ii+2]*delta[2]+min[2] );
    glColor3f(1.0, 1.0, 1.0);
#endif
}
