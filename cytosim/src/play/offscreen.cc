//CVS: $Id: offscreen.cc,v 1.2 2005/04/21 16:08:00 foethke Exp $

#include "offscreen.h"
#include "iomessages.h"

//-----------------------------------------------------------------------


//get the Apple offscreen rendering routines on Macs
#if ( defined __APPLE__ ) && ( ! defined OFF_SCREEN_RENDERING )
  #include "offscreen_mac.cc"
  #define OFF_SCREEN_RENDERING
#endif

//get the X-windows offscreen rendering routines on Linux
#if ( defined __linux ) && ( ! defined OFF_SCREEN_RENDERING )
  #include "offscreen_linux.cc"
  #define OFF_SCREEN_RENDERING
#endif

//setting windows offscreen rendering-copied from Linux -athale
//!!-- doesn't work!
//set dummy routines otherwise (Windows- cygwin offscreen)
//#ifndef OFF_SCREEN_RENDERING
#if ( defined _WIN32 ) && ( ! defined OFF_SCREEN_RENDERING )
  #include "offscreen_cygwin.cc"
  #define OFF_SCREEN_RENDERING
#endif

/*
int OffScreen::open(const int, const int) {
  MSG("ABORT: this player does not have off-screen rendering capabilities\n");
  return 1;
}

void OffScreen::close() {
}
 */
 


