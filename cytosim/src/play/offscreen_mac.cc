//RCS: $Id: offscreen_mac.cc,v 1.1 2005/03/13 22:28:17 nedelec Exp $
// ------------------  off-screen rendering on the mac -------------------------
// Apple introduced PBuffers from Mac OS X 10.3, 
// for hardware accelerated offscreen rendering. 
// It is supported on all platforms with OpenGL 1.2, for PBuffer
// this is C++ code, but could be translated to C, just by  putting all the
// variable declaration in the beginning of the function

#include <cstdio>

#ifdef __APPLE__
    #include <OpenGL/gl.h>
    #include <OpenGL/glu.h>
    #include <AGL/agl.h>
#endif

//the context:
AGLContext aglContext;

//the offscreen drawable
AGLPbuffer aglPbuffer;


//----------------------------------------------------------------
// report of AGL errors
OSStatus aglReportError(const char msg[])
{
	GLenum err = aglGetError();
	if (AGL_NO_ERROR != err) {
		fprintf (stderr, "AGL error in %s: %s\n", msg, (char *) aglErrorString(err));
	}
	// ensure we are returning an OSStatus noErr if no error condition
	if (err == AGL_NO_ERROR)
		return noErr;
	else
		return (OSStatus) err;
}


//----------------------------------------------------------------

int OffScreen::open(const int window_width, const int window_height)
{

// GLint minor, major;
//  aglGetVersion( &major, &minor );
//  MSG("AGL version %i.%i detected\n", int(major), int(minor) );

    
  //pixel format:
  AGLPixelFormat aglPixFmt; 

  //attributes for the pixel format:
	GLint attribList[] = { AGL_RGBA, // AGL_OFFSCREEN,
    AGL_RED_SIZE, 4, 
    AGL_GREEN_SIZE, 4, 
    AGL_BLUE_SIZE,  4,
    AGL_ALPHA_SIZE, 4,
#if (DIM==3)
    AGL_DEPTH_SIZE, 8,
#endif
    AGL_CLOSEST_POLICY, AGL_NO_RECOVERY,
    AGL_NONE };
  
  // find pixel format
  aglPixFmt = aglChoosePixelFormat ( NULL, 0, attribList );
  
  if (aglPixFmt == 0) {
    aglReportError("aglChoosePixelFormat()");
    MSG("play: off-screen rendering aborted: aglChoosePixelFormat() returned NULL\n");
    return 1;
  }
  
  //create context:
  aglContext = aglCreateContext( aglPixFmt, NULL );
   
  if ( aglReportError("aglCreateContext()") != noErr ) {
    MSG("play: off-screen rendering aborted due to above error\n");
    return 2;
  }

 
  if ( ! aglSetCurrentContext(aglContext) ) {
		aglReportError("aglSetCurrentContext()");
    MSG("play: off-screen rendering aborted: aglSetCurrentContext() returned NULL\n");
		return 3;
	}
  
  
  //check that the extension for PBuffers is supported
  const GLubyte * glExtensions = glGetString(GL_EXTENSIONS);
  GLubyte pbExtension[] = "GL_APPLE_pixel_buffer";
  if ( ! gluCheckExtension( pbExtension, glExtensions) ) {
    MSG("play: off-screen rendering aborted: GL_APPLE_pixel_buffer is not supported: install OS X 10.3 or better\n");
    return 4;
  }
  
  
  // PBuffer Creation /
  long kPbufferMaxLevels = 0;
  if ( ! aglCreatePBuffer( window_width, window_height, GL_TEXTURE_2D, GL_RGBA, kPbufferMaxLevels, &aglPbuffer)) {
    aglReportError("aglCreatePBuffer()");
    MSG("play: off-screen rendering aborted: failed to allocate a PBuffer %i x %i\n", window_width, window_height);
    MSG("      You may try a smaller display window: eg. play -s256 -m\n");
    return 5;
  }
  
  // do attach drawable to correct virtual screen and context setting at draw time
  // this means we will not set up the GL state until that time
  GLint vs = aglGetVirtualScreen (aglContext);
	aglReportError("aglGetVirtualScreen()");
	if (!aglSetPBuffer (aglContext, aglPbuffer, 0, 0, vs)) {
		aglReportError("aglGetVirtualScreen()");
    MSG("play: off-screen rendering aborted: aglSetPBuffer() returned NULL\n");
		return 6;
	}
  
  //initialize OpenGL
  return 0;
}

//----------------------------------------------------------------
void OffScreen::close()
{
  aglDestroyContext( aglContext );
  aglDestroyPBuffer( aglPbuffer );
}
