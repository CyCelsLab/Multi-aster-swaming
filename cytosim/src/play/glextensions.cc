//RCS: $Id: glextensions.cc,v 1.7 2005/03/30 12:29:51 nedelec Exp $

#include <cctype>
#include "glut.h"
#include "glextensions.h"
#include "iomessages.h"

//------------------------------------------------------------------------------
void glExtensions::glPrintErrors(const char * name_of_calling_function)
{
  GLenum glError;
  do {
    glError = glGetError();
    switch( glError ) {
      case GL_NO_ERROR:  
        return;
      case GL_INVALID_ENUM:
        MSG("OpenGL error GL_INVALID_ENUM in %s\n", name_of_calling_function); 
        break;
      case GL_INVALID_VALUE: 
        MSG("OpenGL error GL_INVALID_VALUE in %s\n", name_of_calling_function); 
        break;
      case GL_INVALID_OPERATION: 
        MSG("OpenGL error GL_INVALID_OPERATION in %s\n", name_of_calling_function);
        break;
      case GL_STACK_OVERFLOW:
        MSG("OpenGL error GL_STACK_OVERFLOW in %s\n", name_of_calling_function); 
        break;
      case GL_STACK_UNDERFLOW:
        MSG("OpenGL error GL_STACK_UNDERFLOW in %s\n", name_of_calling_function);
        break;
      case GL_OUT_OF_MEMORY:
        MSG("OpenGL error GL_OUT_OF_MEMORY in %s\n", name_of_calling_function); 
        break;
      default:
        MSG("OpenGL unknown error in %s\n", name_of_calling_function);
        break;
    }
  } while ( glError != GL_NO_ERROR );
}

//------------------------------------------------------------------------------
int glExtensions::saveImageAsPPM(FILE * file, const int xpos, const int ypos, const int width, const int height)
{
  //allocate a chunck of memory to hold the image:
  char * pixels = new char[ width * height * 3 ];
  
  if (pixels == 0) {
    fprintf(stderr, "play_glut.cc::saveImageAsPPM() memory allocation failed\n");
    return 1;
  }
  
  //selecting the right OpenGL buffer for reading:
  //this is illegal for off-screen rendering, and not needed in general
  //glReadBuffer( GL_BACK );
  
  glPrintErrors("saveImageAsPPM1()");
  
  //set the alignment ?
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  
  //read the pixel values, from top-left corner:
  glReadPixels( xpos, ypos, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels );
  
  glPrintErrors("saveImageAsPPM()");
  
  //write the file in Portable pixmap file format:
  //we use the 'raw' binary format starting with P6 
  //The format is described in the man pages of ppm (man ppm)
  fprintf(file, "P6\n");
  fprintf(file, "%i %i\n", width, height);
  fprintf(file, "255\n");
  
  //write the pixels binary, line by line:
  for(int ii = height-1; ii >= 0; --ii )
    fwrite(pixels + 3*ii*width, 1, 3*width, file);
  
  //cleanup
  delete [] pixels;
  fclose( file );
  return 0;
}

//------------------------------------------------------------------------------
//                         OVERLOADS FOR glVertex
//------------------------------------------------------------------------------

#ifdef SINGLE

void glExtensions::glVertex( const real * V )
{
#if ( DIM == 1 )
  glVertex2f( V[0], 0.0 );
#elif ( DIM == 2 )
  glVertex2f( V[0], V[1] );
#elif ( DIM == 3 )
  glVertex3f( V[0], V[1], V[2] );
#endif
}

void glExtensions::glVertex( const real * V, const real pseudoY )
{
#if ( DIM == 1 )
  glVertex2f( V[0], pseudoY );
#elif ( DIM == 2 )
  glVertex2f( V[0], V[1] );
#elif ( DIM == 3 )
  glVertex3f( V[0], V[1], V[2] );
#endif
}

void glExtensions::glVertex2( const real x, const real y) {
  glVertex2f( x, y );
}
void glExtensions::glVertex3( const real x, const real y, const real z ) {
  glVertex3f( x, y, z );
}
void glExtensions::glTranslate( const real x, const real y, const real z ) {
  glTranslatef( x, y, z );
}


void glExtensions::glTranslatev(const real* dist)
{
#if ( DIM == 1 )
  glTranslatef( dist[0], 0, 0);
#elif ( DIM == 2 )
  glTranslatef( dist[0], dist[1], 0);
#elif ( DIM == 3 )
  glTranslatef( dist[0], dist[1], dist[2]);
#endif
}

#else

void glExtensions::glVertex( const real * V )
{
#if ( DIM == 1 )
  glVertex2d( V[0], 0.0 );
#elif ( DIM == 2 )
  glVertex2d( V[0], V[1] );
#elif ( DIM == 3 )
  glVertex3d( V[0], V[1], V[2] );
#endif
}

void glExtensions::glVertex( const real * V, const real pseudoY )
{
#if ( DIM == 1 )
  glVertex2d( V[0], pseudoY );
#elif ( DIM == 2 )
  glVertex2d( V[0], V[1] );
#elif ( DIM == 3 )
  glVertex3d( V[0], V[1], V[2] );
#endif
}


void glExtensions::glVertex2( const real x, const real y ) {
  glVertex2d( x, y );
}
void glExtensions::glVertex3( const real x, const real y, const real z ) {
  glVertex3d( x, y, z );
}
void glExtensions::glTranslate( const real x, const real y, const real z ) {
  glTranslated( x, y, z );
}


void glExtensions::glTranslatev(const real* dist)
{
#if ( DIM == 1 )
  glTranslated( dist[0], 0, 0);
#elif ( DIM == 2 )
  glTranslated( dist[0], dist[1], 0);
#elif ( DIM == 3 )
  glTranslated( dist[0], dist[1], dist[2]);
#endif
}

#endif


//------------------------------------------------------------------------------
//================================  COLORS  ====================================
//------------------------------------------------------------------------------



// r,g,b values are from 0 to 1
// h = [0, 360], s = [0,1], v = [0,1]
// if s == 0, then h = -1 (undefined)

void glExtensions::RGBtoHSV( const real r, const real g, const real b, real *h, real *s, real *v )
{
  real min, max, delta;
  min = minT(minT( r, g ), b );
  max = maxT(maxT( r, g ), b );
  *v = max;                         // v
  delta = max - min;
  if( max != 0 )
    *s = delta / max;               // s
  else {
    // r = g = b = 0                // s = 0, v is undefined
    *s = 0;
    *h = -1;
    return;
  }
  if( r == max )
    *h = ( g - b ) / delta;         // between yellow & magenta
  else if( g == max )
    *h = 2 + ( b - r ) / delta;     // between cyan & yellow
  else
    *h = 4 + ( r - g ) / delta;     // between magenta & cyan
  *h *= 60;                         // degrees
  if( *h < 0 )
    *h += 360;
}


//------------------------------------------------------------------------------
void glExtensions::HSVtoRGB( const real h, const real s, const real v, real *r, real *g, real *b )
{
  int i;
  real f, p, q, t;
  if( s == 0 ) {
    // achromatic (grey)
    *r = *g = *b = v;
    return;
  }
  real hc = h/60;                  // sector 0 to 5
  i = (int)floor( hc );
  f = hc - i;                      // factorial part of h
  p = v * ( 1 - s );
  q = v * ( 1 - s * f );
  t = v * ( 1 - s * ( 1 - f ) );
  switch( i ) {
  case 0: *r = v; *g = t; *b = p; break;
  case 1: *r = q; *g = v; *b = p; break;
  case 2: *r = p; *g = v; *b = t; break;
  case 3: *r = p; *g = q; *b = v; break;
  case 4: *r = t; *g = p; *b = v; break;
  case 5: *r = v; *g = p; *b = q; break;
  default: *r = 1; *g = 1; *b = 1; break;
  }
}


//------------------------------------------------------------------------------
//set a RGB color as a function of a value h in [0, 1]
//the colors are in the order: blue, green, red, white
void glExtensions::setRainbowColor( const real h, real *r, real *g, real *b )
{
  int i = (int)floor( 5 * h );
  
  if ( i < 0 ) {
    *r = 0; *g = 0, *b = 1;
    return;
  }
  
  real f = 5 * h - i;

  switch( i ) {
  case 0:  *r = 0;   *g = f;   *b = 1;   break;
  case 1:  *r = 0;   *g = 1;   *b = 1-f; break;
  case 2:  *r = f;   *g = 1;   *b = 0;   break;
  case 3:  *r = 1;   *g = 1-f; *b = 0;   break;
  default: *r = 1;   *g = 0;   *b = 0;   break;
  }
}

//------------------------------------------------------------------------------
void glExtensions::glColorRainbow( const real h )
{
  real r, g, b;
  setRainbowColor( h, &r, &g, &b );
  glColor3f( r, g, b );
}

//------------------------------------------------------------------------------
void glExtensions::glLightningColorRainbow( const real x )
{
  real r, g, b;
  setRainbowColor( x, &r, &g, &b );
  GLfloat rgba[4] = {r, g, b, 1.0};
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, rgba);
}

//------------------------------------------------------------------------------
void glExtensions::glColor( const long color )
{
  GLubyte  r = ( color >> 24 ) & 0xFF; //((GLubyte*)&color)[0];
  GLubyte  g = ( color >> 16 ) & 0xFF; //((GLubyte*)&color)[1];
  GLubyte  b = ( color >> 8  ) & 0xFF; //((GLubyte*)&color)[2];
  GLubyte  a = ( color       ) & 0xFF; //((GLubyte*)&color)[3];

  glColor4ub( r, g, b, a );
}

//------------------------------------------------------------------------------
void glExtensions::glColor( long color, const int dimm_factor )
{
  GLubyte  r = ( color >> 24 ) & 0xFF;
  GLubyte  g = ( color >> 16 ) & 0xFF;
  GLubyte  b = ( color >> 8  ) & 0xFF;
  GLubyte  a =(( color       ) & 0xFF) >> dimm_factor;
  
  glColor4ub( r, g, b, a );
}

//------------------------------------------------------------------------------
 void glExtensions::glColor( const long color1, const long color2, const int dimm_factor )
 {
   //we combine the two colors half-half, on the byte level
   long combine = (( color1 >> 1 ) & 0x7F7F7F7F) + (( color2 >> 1 ) & 0x7F7F7F7F);

   GLubyte  r = ( combine >> 24 ) & 0xFF;
   GLubyte  g = ( combine >> 16 ) & 0xFF;
   GLubyte  b = ( combine >> 8  ) & 0xFF;
   GLubyte  a =(( combine       ) & 0xFF) >> dimm_factor;
   
   glColor4ub( r, g, b, a );
 }

//------------------------------------------------------------------------------
void glExtensions::glLightningColor( const GLenum face, const long color)
{
  GLubyte  r = ( color >> 24 ) & 0xFF;
  GLubyte  g = ( color >> 16 ) & 0xFF;
  GLubyte  b = ( color >> 8  ) & 0xFF;
  GLubyte  a =   color         & 0xFF;
  
  GLfloat rgba[4] = {(GLfloat)r/255., (GLfloat)g/255., (GLfloat)b/255., (GLfloat)a/255.};
  glMaterialfv(face, GL_AMBIENT_AND_DIFFUSE, rgba);
}

//------------------------------------------------------------------------------
void glExtensions::glLightningColor( const GLenum face, const long color, const int dimm_factor)
{
  GLubyte  r = ( color >> 24 ) & 0xFF;
  GLubyte  g = ( color >> 16 ) & 0xFF;
  GLubyte  b = ( color >> 8  ) & 0xFF;
  GLubyte  a =   color         & 0xFF;
  a >>= dimm_factor;

  GLfloat rgba[4] = {(GLfloat)r/255., (GLfloat)g/255., (GLfloat)b/255., (GLfloat)a/255.};
  glMaterialfv(face, GL_AMBIENT_AND_DIFFUSE, rgba);
}

//------------------------------------------------------------------------------
void glExtensions::glLightningColorOI( const long colorOutside, const long colorInside )
{
  unsigned char  ro = ( colorOutside >> 24 ) & 0xFF;
  unsigned char  go = ( colorOutside >> 16 ) & 0xFF;
  unsigned char  bo = ( colorOutside >> 8  ) & 0xFF;
  unsigned char  ao =   colorOutside         & 0xFF;
  
  unsigned char  ri = ( colorInside >> 24 ) & 0xFF;
  unsigned char  gi = ( colorInside >> 16 ) & 0xFF;
  unsigned char  bi = ( colorInside >> 8  ) & 0xFF;
  unsigned char  ai =   colorInside         & 0xFF;
  
  GLfloat rgbaOutside[4] = {(GLfloat)ro/255., (GLfloat)go/255., (GLfloat)bo/255., (GLfloat)ao/255.};
  GLfloat rgbaInside[4]  = {(GLfloat)ri/255., (GLfloat)gi/255., (GLfloat)bi/255., (GLfloat)ai/255.};
  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, rgbaOutside);
  glMaterialfv(GL_BACK,  GL_AMBIENT_AND_DIFFUSE, rgbaInside);
}

//------------------------------------------------------------------------------
//for some reason, compilation fails if we overload the function glClearColor
void glExtensions::gleClearColor( const long color )
{
  GLfloat  r = GLfloat( ( color >> 24 ) & 0xFF ) / 255.;
  GLfloat  g = GLfloat( ( color >> 16 ) & 0xFF ) / 255.;
  GLfloat  b = GLfloat( ( color >> 8  ) & 0xFF ) / 255.;
  GLfloat  a = GLfloat( ( color       ) & 0xFF ) / 255.;
  glClearColor( r, g, b, a );
}

//------------------------------------------------------------------------------

/*
 ///class to hold RGB component of a color
class Color {
public:
  char * name;
  int    r, g, b;
  long   hex;
};

Color colors[] = {
  {"black",      0,   0,   0,    0x000000FF },
  {"blue",       0,   0,   255,  0x0000FFFF },
  {"blueTint",   175, 214, 255,  0xAFD7FFFF },
  {"brown",      175, 117, 89,   0xAF7559FF },
  {"cyan",       0,   255, 255,  0x00FFFFFF },
  {"gold",       255, 156, 0,    0xFC9C00FF },
  {"grey",       125, 125, 125,  0x7D7D7DFF },
  {"green",      0,   255, 0,    0x00FF00FF },
  {"greenBlue",  46,  139, 87,   0x2E8B57FF },
  {"greenTint",  152, 255, 179,  0x98FFB3FF },
  {"hotPink",    255, 0,   101,  0xFF0065FF },
  {"magenta",    255, 0,   255,  0xFF00FFFF },
  {"orange",     255, 165, 0,    0xFFA500FF },
  {"pink",       255, 101, 117,  0xFF6575FF },
  {"pinkTint",   255, 171, 187,  0xFFABBBFF },
  {"purple",     160, 32,  240,  0xA020F0FF },
  {"red",        255, 0,   0,    0xFF0000FF },
  {"redOrange",  255, 69,  0,    0xFF4500FF },
  {"seaGreen",   0,   250, 109,  0x00FA6DFF },
  {"skyBlue",    58,  144, 255,  0x3A90FFFF },
  {"violet",     238, 130, 238,  0xEE82EEFF },
  {"white",      255, 255, 255,  0xFFFFFFFF },
  {"yellow",     255, 255, 0,    0xFFFF00FF },
  {"yellowTint", 246, 246, 117,  0xF6F675FF }
};
*/

//------------------------------------------------------------------------------
void glExtensions::displayText(const char text[], const long color, 
                               const int window_width, const int window_height, const int position)
{
  //compute the max width of all the lines in the given text
  const char * c = text;
  int cw = 0, maxwidth = 0;
  while( *c != '\0' ) {
    if (( *c == '\n' ) || ( *c == '\r' ))
      cw = 0;
    else {
      if ( ++cw > maxwidth ) 
        maxwidth = cw;
    }
    ++c;
  }
    
  GLint shift, raster_position[2];

  switch( position ) {
    case 0: //bottom-left
      raster_position[0] = 10;
      raster_position[1] = 10;
      shift = +13;
      break;
    case 1: //bottom-right
      raster_position[0] = window_width - 8*maxwidth - 10;
      raster_position[1] = 10;
      shift = +13;
      break;
    case 2: //top-right
      raster_position[0] = window_width - 8*maxwidth - 10;
      raster_position[1] = window_height - 15;
      shift = -13;
      break;
    default:
    case 3: //top-left
      raster_position[0] = 10;
      raster_position[1] = window_height - 15;
      shift = -13;
      break;
  }

  //set the matrices
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho( 0, window_width, 0, window_height, 0, 1 );
  
  int line = 0;    
  // set color and position of text
  glColor(color);
  glLightningColor(GL_FRONT, color);
  glRasterPos2i(raster_position[0], raster_position[1] + shift * line);

  // draw the string character per character in a small font:
  for(const char * p = text; *p; ++p) {
    if ( *p == '\n' ) 
      glRasterPos2i(raster_position[0], raster_position[1] + shift * (++line));
    else if ( isprint( *p ) )
      glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *p);
  }
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}
