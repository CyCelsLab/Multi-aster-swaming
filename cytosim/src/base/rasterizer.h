//RCS: $Id: rasterizer.h,v 2.3 2004/11/11 16:47:01 nedelec Exp $

// Simple 2D and 3D rasterizer. 
// (I have searched the web, but could not find a suitable 3D rasterizer)
//
// The functions paintPolygon?D(function paint(), polygon) calls
// the given function paint() for every point of INTEGER coordinate 
// inside a given polygon. This polygon is usually specified as a 
// list of points,  of which it is the convex hull.
// The points defining the polygon do not need to be integers
//
// paintFatLine?D( points p, q; width ) invoke the rasterizer, for a  
// polygon obtained by inflating a segment specified by two ends p and q, 
// by the given width. It is used to 'paint' all points of integer 
// coordinates closer than a certain distance to the segment [pq].


// F.Nedelec EMBL Nov 2002, nedelec@embl.de

#ifndef RASTERIZER_H
#define RASTERIZER_H

#include "types.h"
/// Simple 2D and 3D rasterizer. 
namespace Rasterizer 
{
    
    //debug printf() function:
    void printPolygon( int nb, const real xy[] );

    ///functions to calculate the convex hull in 2D from a set of points:
    //the array xy[] contains the coordinates of all the 'nbpts' points, 
    //at output, they are ordered anti-clockwised
    //the optional array dxy[] is reordered at the same time as xy[]:
    void convexHull2D( real xy[],             ///< coordinates of the points
                       int * nbpts,           ///< number of points
                       real dxy[] = 0 );      ///< optional array to be re-ordered with xy[]
    
    
    
    /// Rasterizer function in 1D:
    void paintPolygon1D(void paint(const int, const int, const int),    ///< painting function
                        const int nbpts,      ///< number of points
                        const real x[],       ///< coordinates of the points ( x, x...)
                        const int yy = 0,     ///< second coordinate, passed as argument to paint()
                        const int zz = 0);    ///< third coordinate, passed as argument to paint()
    
    /// Rasterizer function in 1D:
    void paintFatLine1D(void paint(const int, const int, const int),    ///< painting function
                        const real p[],        ///< first point at the end of the line [dim=1]
                        const real q[],        ///< second point at the other end of the line [dim=1]
                        const real width,      ///< width by which the line [pq] is extended, to make a round cylinder
                        const real min[],      ///< phase of the grid [dim=1]
                        const real delta[]);   ///< period for the grid [dim=1]
    
    
    
    /// Rasterizer function in 2D:
    // paintPolygon2D() calls paintPoint(x,y,zz) for every point x,y 
    // where x and y are integer, and is inside the polygon
    // of nbpts points who's vertex are given in xy[]. The polygon
    // should be convex, and points should be ordered in a consistent
    // manner (clockwise or counter-CW).    
    void paintPolygon2D(void paint(const int, const int, const int),    ///< painting function
                        const int nbpts,      ///< number of points
                        const real xy[],      ///< coordinates of the points ( x y, x y...)
                        const int zz = 0);    ///< third coordinate, passed as argument to paint()
     
    
    /// invoke the rasterizer, for a polygon obtained by inflating the segment [pq]
    /// by the given width. The polygon is a parallelepiped
    void paintFatLine2D(void paint(const int, const int, const int),    ///< painting function
                        const real p[],        ///< first point at the end of the line [dim=2]
                        const real q[],        ///< second point at the other end of the line [dim=2]
                        const real width);     ///< width by which the line [pq] is extended, to make a round cylinder
    
    /// this version is sligtly faster, and should be used if the function is to 
    /// be called with many segment of similar length, as in the case of microtubules.
    void paintFatLine2D(void paint(const int, const int, const int),    ///< painting function
                        const real p[],        ///< first point at the end of the line [dim=2]
                        const real q[],        ///< second point at the other end of the line [dim=2]
                        const real width,      ///< width by which the line [pq] is extended, to make a round cylinder
                        const real min[],      ///< phase of the grid [dim=2]
                        const real delta[],    ///< period for the grid [dim=2]
                        real scaling = 0 );    ///< scaling factor = width / distance( pq )
    
    
    ///---------------------------------- Rasterizer functions in 3D:
    
    // the polygon is the convex hull of the 'nbpts' points given in xyz[]
    // algorithm: we section at each integral Z, collect the intersection of
    // all possible lines connecting two points, calculate the convex hull of
    // all these points and call paintPolygon2D.
    // this is far from beeing optimal
    ///Rasterizer function in 3D
    void paintPolygon3D(void paint(const int, const int, const int), 
                        const int nbpts, real xyz[]);
    
    ///Rasterizer function in 3D (old/slower)
    void paintFatLine3D_old(void paint(const int, const int, const int),   ///< painting function
                            const real p[],       ///< first point at the end of the line [dim=3]
                            const real q[],       ///< second point at the other end of the line [dim=3]
                            const real width,     ///< width by which the line [pq] is extended, to make a round cylinder
                            const real min[],     ///< phase of the grid [dim=3]
                            const real delta[],   ///< period for the grid [dim=3]
                            real scaling = 0 );   ///< scaling factor = width / distance( pq )

    
    ///compare function for qsort (internal)
    int comp_higher( const void *, const void *);
    
    /// Rasterizer function in 3D, where the sides making the convex hull are already known:
    // the polygon is the convex hull of the 'nbpts' points given in xyzi[],
    // containing the three coordinate + the index of the point, a short int
    // cast into the 'real' type.
    // sides[] gives information on which pair of points should be considered:
    // sides[ ii + nbpts * jj ] == 1 if points of index ii and jj are connected,
    // i.e. they are one of the edge of the 3D polygon
    
    void paintPolygon3D(void paint(const int, const int, const int),  ///< painting function
                        const int nbpts,                              ///< number of points
                        real xyzi[],                                  ///< coordinated of the point + index : ( x,y,z, i)
                        const char sides[]);                          ///< array specifying which points are connected
   
    /// Rasterizer function in 3D
    void paintFatLine3D(void paint(const int, const int, const int),   ///< painting function
                        const real p[],       ///< first point at the end of the line [dim=3]
                        const real q[],       ///< second point at the other end of the line [dim=3]
                        const real width,     ///< width by which the line [pq] is extended, to make a round cylinder
                        const real min[],     ///< phase of the grid [dim=3]
                        const real delta[],   ///< period for the grid [dim=3]
                        real scaling = 0 );   ///< scaling factor = width / distance( pq )
    
};

#endif //RASTERIZER_H

