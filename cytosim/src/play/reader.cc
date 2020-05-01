//RCS: $Id: reader.cc,v 2.13 2005/04/22 18:53:48 nedelec Exp $

#include "reader.h"
#include "sim.h"
#include "exceptions.h"
#include "iowrapper.h"


//-------------------------------------------------------------------------
void SimReader::clearBuffer() 
{
  //MSG(9,"SimReader::clearBuffer()\n");
    eofFrame = MAX_FRAME;
    minFrame = MAX_FRAME;
    maxFrame = -1;
    for( int ii = 0; ii < MAX_FRAME; ++ii )
      framestart[ ii ] = 0;
}

//-------------------------------------------------------------------------
bool SimReader::eof() const {
  return IO.eof(); 
}

//-------------------------------------------------------------------------
int SimReader::openFile(const char * filename)
{
  bool tryGunzip = false;
  int error = NO_ERROR;

    try {
      error = IO.openInputFile( filename );
      tryGunzip = ( error != NO_ERROR );
    } catch( IOException e ) {
      tryGunzip = true;
    }
    
    //try if the file is gzipped:
    if ( tryGunzip  && (strlen(filename) < 200)) {
      //we try adding ".gz" to see if that works:
      char cmd[256];
      snprintf( cmd, sizeof(cmd), "gunzip %s.gz", filename );
      MSG("SimReader:: trying system( %s )\n", cmd);
      if ( system( cmd ) )
        return 1; //the system() call failed
      //try again to open the (gunziped) file:
      try {
        error = IO.openInputFile( filename );
        if ( error != NO_ERROR ) return error;
      } catch( IOException e ) {
        return 1; 
      }
    }

    try {
      IO.checkInputFile();
    } catch( IOException e ) {
      return 1;
    }
    
    clearBuffer();
    return NO_ERROR;
}

//-------------------------------------------------------------------------
void SimReader::closeFile()
{
    IO.closeInputFile();
}

//-------------------------------------------------------------------------
void SimReader::rewindFile() 
{
    if ( IO.getInputFile() ) {
        clearerr( IO.getInputFile() );
        rewind( IO.getInputFile() );
    }
}

//-------------------------------------------------------------------------
FILE * SimReader::getInputFile()
{
    return IO.getInputFile(); 
}

//-------------------------------------------------------------------------

int SimReader::seekFrame( int frame )
//find the starting point of a frame by brute force !
//searching the file from the current position, and
//then again from the beggining, for "#frm "
{
    int found = 0;
    unsigned long pos;
    char line[16];
    
    if ( IO.error() ) 
        return 1;
    
    //by default, we start from the beggining:
    rewindFile();
        
    //go as far possible before the frame:
    int ff = minT(frame, MAX_FRAME-1);
    while((ff > 1)  && (framestart[ff] <= 0))
        --ff;
    
    if ( framestart[ ff ] > 0 )
        IO.seek( framestart[ ff ] );
        
    MSG(9, "SimReader: seekFrame(%i) at pos %li\n", frame, IO.getPos());
    //in the first round, we start searching from there:
    //start searching from the beggining:
    
    do {
        
        IO.skipUntil(FRAME_START_TAG);
        pos = IO.getPos();
        
        IO.readLine( line, 16 );
        if ( line == strstr(line, FRAME_START_TAG )) 
          if ( 1 == sscanf(line+strlen(FRAME_START_TAG), "%u", &found )) {
          MSG(9, "SimReader: seekFrame(%i) found <%s%i> at pos %li\n",
              frame, FRAME_START_TAG, found, pos);
          
          //remember the position:
          if ( found < MAX_FRAME )
            framestart[ found ] = pos;
           
          //we found it!
          if ( found == frame ) {
            IO.seek(pos);
            return NO_ERROR;
          }
            
          //we are too far in the file:
          if ( found > frame ) 
            return 2;
        }
        
    } while ( !IO.eof() );
    
    MSG(9, "SimReader: seekFrame(%i) EOF at pos %li\n", frame, IO.getPos());
    
    if ( found )
      eofFrame = found+1;
    else {
      if ( eofFrame < frame )
        eofFrame = frame;
    }
    
    return 3;
}

//wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww

int SimReader::readFrameNoCatch( int frame ) 
{
  MSG(9, "SimReader: readFrame(%i) eofFrame=%i\n", frame, eofFrame);
  
  if ( IO.error() ) 
    return 2;
  
  //the frame we are looking for might already be in the buffer:
  if ( frame == sim.frameInBuffer() )
    return NO_ERROR; 
  
  //frame below limit, or frame=1:
  if ( frame <= 0 ) {
    frame = 0; 
    rewindFile(); 
  }
  
  //---------------------try to find the frame start by brute force from there:
  
  if ( NO_ERROR != seekFrame( frame ) )
    return 3;
  
  //--------------------read the sim-state from there:
  
  MSG(9, "SimReader: reading file for frame %i at pos %li\n", frame, IO.getPos() );
  
  switch ( sim.readState() ) {
    
    case NO_ERROR: // read was successful:
      
      //the next frame should start there:
      if (( sim.frameInBuffer() + 1 < MAX_FRAME ) && ( framestart[ sim.frameInBuffer() + 1 ] <= 0 ))
        framestart[ sim.frameInBuffer() + 1 ] = IO.getPos();
      
      if ( minFrame > sim.frameInBuffer() )
        minFrame = sim.frameInBuffer();
        
      if ( maxFrame < sim.frameInBuffer() )
        maxFrame = sim.frameInBuffer();
          
      if ( eofFrame <= sim.frameInBuffer() )
        eofFrame = MAX_FRAME;
            
      MSG(9, "SimReader:: readState() found frame %i\n", sim.frameInBuffer());
      
      if ( frame == sim.frameInBuffer() )
        return NO_ERROR; 
      else
        return 4;
      
    case 1: 
      // eof:
      if ( frame < eofFrame )
        eofFrame = frame;
      MSG(9, "SimReader::readNextFrame() EOF at frame %i\n", eofFrame);
      return 1;
      
    default: 
      // internal error: readState() should return 0 or 1
      MSG.error("SimReader:readNextFrame()", "readState() returned an invalid value");
      return 5;	
  }
}

//wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww

int SimReader::readFrame( int frame ) 
{
  try {
    return readFrameNoCatch( frame );
  } catch( IOException e ) {
    MSG(10, "SimReader: ERROR reading frame %i: %s\n", frame, e.getMessage());
    if ( sim.frameInBuffer() + 1 < MAX_FRAME ) 
      framestart[ sim.frameInBuffer() + 1 ] = IO.getPos();
    //we clear the all thing (nicer, but not essential)
    sim.eraseState();
    return 6;
  }
}


//-------------------------------------------------------------------------
int SimReader::readNextFrameNoCatch()
{
  if ( IO.error() ) 
    return 2;

  unsigned long pos = IO.getPos();
  switch ( sim.readState() ) {
            
    case NO_ERROR: {
      // read was successful:
      int frame = sim.frameInBuffer();
      
      if ((frame < MAX_FRAME) && (framestart[ frame ] <= 0))
        framestart[ frame ] = pos;
      
      //the next frame should start from the current position:
      if ((frame+1 < MAX_FRAME) && (framestart[ frame + 1 ] <= 0))
        framestart[ frame + 1 ] = IO.getPos();
      
      if ( minFrame > frame )
        minFrame = frame;
      
      if ( maxFrame < frame )
        maxFrame = frame;
      
      if ( eofFrame <= frame )
        eofFrame = MAX_FRAME;
      
      MSG(9, "SimReader::readNextFrame(): readState() found frame %i\n", frame);
      
    } return NO_ERROR; 
      
    case 1: {
      // eof:      
      MSG(9, "SimReader::readNextFrame(): EOF at frame %i\n", eofFrame);
    } return 1;
      
    default: 
      // internal error: readState() should return 0 or 1
      MSG.error("SimReader:readNextFrame()", "readState() returned invalid value");
      return 5;	
  }
}


//-------------------------------------------------------------------------
int SimReader::readNextFrame()
{
  try {
    return readNextFrameNoCatch();
  } catch( IOException e ) {
    MSG(10, "SimReader: ERROR frame %i: %s\n", sim.frameInBuffer(), e.getMessage());
    if (sim.frameInBuffer() + 1 < MAX_FRAME)
      framestart[ sim.frameInBuffer() + 1 ] = IO.getPos();
    //we clear the all state (nicer, but not essential)
    sim.eraseState();
    return 6;
  }
}
