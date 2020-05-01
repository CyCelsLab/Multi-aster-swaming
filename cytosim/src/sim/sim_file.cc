//RCS: $Id: sim_file.cc,v 2.31 2005/04/22 19:30:48 nedelec Exp $


void SIM::cleanUpRead(const int frameToKeep )
{
  //the order is important, as Asters register microtubules
  cleanUpReadAster( frameToKeep );
  cleanUpReadNucleus( frameToKeep );
  cleanUpReadComplex( frameToKeep );
  cleanUpReadGrafted( frameToKeep );
  cleanUpReadMicrotub( frameToKeep );
  cleanUpReadSolid( frameToKeep );
}

//=================================================================
//=======================INPUT FROM A FILE=========================
//=================================================================

// read simulation state from the file 'filename' for input
int SIM::readState(const char * filename)
{
  try {
    int result = IO.openInputFile(filename);
    if ( result == NO_ERROR ) {
      MSG(1,"reading initial state from file [%s]\n", filename);
      result = readState();
      IO.closeInputFile();
      return result;
    } else {
      MSG(1,"SIM::readState: file [%s] not found or has error\n", filename);
      return result;
    }
  } catch( IOException e ) {
    return 2;
  }
}

//---------------------------------------------------------------------------
//read the first sim-image from file 'infile', at the current position
//return values:
// NO_ERROR   : success
// (1)        : eof
// handles errors by throwing IOexceptions
int SIM::readState()
{
  if ( IO.error() )
    throw IOException("SIM::readState() :: IO.error()");

  if ( false == IO.skipUntil(FRAME_START_TAG)) {
    if ( IO.eof() )
      return 1;
    else
      throw IOException("SIM::readState() could not find FRAME_START_TAG");
  }

  unsigned long start = IO.getPos();

  try {
    int code = readState_new();
    //if an old-state has been detected, we call readState_old()
    if ( code == 2 ) {
      IO.seek(start);
      code = readState_old();
    }
    cleanUpRead( frame_in_buffer );
    return code;
  } catch( IOException e ) {
    cleanUpRead( frame_in_buffer );
    throw e;
  }
}

//---------------------------------------------------------------------------
//read the first sim-image from file 'infile', at the current position
//this reads the old format of coordinate file (prior to Feb 1st 2004),
//return values:
// NO_ERROR   : success
// (1)        : eof
// handles errors by throwing IOexceptions
int SIM::readState_old()
{
  Name name;
  char id[STRING_SIZE] = {'\0','\0'};

  do {

    id[0] = IO.readEntry();
    MSG(81, "parsing old format %c\n", id[0]);

    //check for a comment, that includes tags for start and end of frame
    if( id[0] == '#' ) {
      IO.readLine( id+1, sizeof(id)-1 );

      if ( id == strstr(id, FRAME_START_TAG) )
        sscanf( id+strlen(FRAME_START_TAG), "%i", &frame_in_buffer);

      //that signals the end of the frame
      if ( id == strstr(id, FRAME_END_TAG) )
        return NO_ERROR;

      continue;
    }

    //simulation 'real' time in seconds
    if ( id[0] == 'T' ) {
      starting_time = IO.readFloatAscii();
      MSG(80, "time = %f\n", starting_time);
      continue;
    }

    //file format of result.out / historical tag
    if ( id[0] == 'H' ) {
      IO.readInputFormat();
      continue;
    }

    //binary tag
    if ( id[0] == 'B' ) {
      IO.readBinarySignature();
      continue;
    }

    //dimension in space
    if ( id[0] == 'D' ) {
      IO.readInputDIM();
      if ( IO.getInputDIM() > DIM )
        MSG.warning("SIM::read", "projecting state in lower dimension");
      continue;
    }

    //random seed
    if ( id[0] == 'R' ) {
      IO.skipChar();
      IO.readUInt32();
      continue;
    }

    //box shape and size
    if ( id[0] == 'Z' ) {
      MP.boxshape = IO.readStringInt();
      MSG(80, "boxshape = %i\n", MP.boxshape);
      for( int ii = 0; ii < 6; ++ii )
        MP.boxsize[ ii ] = IO.readFloatAscii();
      setSpace();
      //printf("Z %i %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
      //     MP.boxshape, MP.boxsize[0], MP.boxsize[1], MP.boxsize[2],
      //     MP.boxsize[3], MP.boxsize[4], MP.boxsize[5]);
      continue;
    }

    //microtub:
    if ( id[0] == 'f' ) {
      name = IO.readUInt16();
      if ( name <= 0 ) throw IOException("invalid name value in SIM::readState_old");
      findMicrotub(name, 1) -> read();
      continue;
    }

    //Solid:
    if ( id[0] == 's' ) {
      name = IO.readUInt16();
      if ( name <= 0 ) throw IOException("invalid name value in SIM::readState_old");
      findSolid(name, 1) -> read();
      continue;
    }

    //aster:
    if ( id[0] == 'a' ) {
	  name = IO.readUInt16();  // nk comments this out on 6 oct 2014
      //name = IO.readUInt32();								
      if ( name <= 0 ) throw IOException("invalid name value in SIM::readState_old");
      findAster(name, 1) -> read();
      continue;
    }

    //grafted:
    if ( id[0] == 'g' ) {
      name = IO.readUInt16();
      if ( name <= 0 ) throw IOException("invalid name value in SIM::readState_old");
      findGrafted(name, 1) -> read();
      continue;
    }

    //complex:
    if ( id[0] == 'c' ) {
      name = IO.readUInt16();
      if ( name <= 0 ) throw IOException("invalid name value in SIM::readState_old");
      findComplex(name, 1) -> read();
      continue;
    }

    //complex with a long name:
    if ( id[0] == 'x' ) {
      name = IO.readUInt32();
      if ( name <= 0 ) throw IOException("invalid name value in SIM::readState_old");
      findComplex(name, 1) -> read();
      continue;
    }

    /*
    //map 2008-10-20
    if ( id[0] == 'P' ) {
	   name = IO.readUInt32();
	   if ( name <= 0 ) throw IOException("invalid name value in SIM::readState_old");
	   findMap(name, 1) -> read();
	   continue;
    }
    */


    // we just skip the line for the moment
    IO.readLine(id+1, sizeof(id)-1);
    MSG(70, "readState: skipped %s\n", id);

  } while( !IO.eof() );

  return 1;
}

//---------------------------------------------------------------------------
//read the first sim-image from file 'infile', at the current position
//return values:
// NO_ERROR   : success
// (1)        : eof
// (2)        : detected an old format
// handles errors by throwing an IOException
int SIM::readState_new()
{
  Name name;
  char line[STRING_SIZE] = {'\0','\0','\0'};
  //reset the comment associated with the frame:
  frame_info[0] = '\0';

  do {

    if ( IO.error() )
      throw IOException("SIM::readState() :: IO.error()");

    if ( IO.readRecordTag(line) )
      return 1;  // EOF

    MSG(81, "parsing recordTag %c%c\n", line[0], line[1]);

    // we cannot use a switch with two char, so have many if()
    // we try the most probable possibilities first:

    //complex:
    if (IO.compareRecordTag(line, "cx")) {
      name = IO.readUInt16();
      if ( name <= 0 ) throw IOException("invalid Complex name value in SIM::readState");
      findComplex(name, 1) -> read();
      continue;
    }

    //complex with a long name:
    if (IO.compareRecordTag(line, "cX")) {
      name = IO.readUInt32();
      if ( name <= 0 ) throw IOException("invalid Complex name value in SIM::readState");
      findComplex(name, 1) -> read();
      continue;
    }

    //grafted:
    if (IO.compareRecordTag(line, "gh")) {
      name = IO.readUInt16();
      if ( name <= 0 ) throw IOException("invalid Grafted name value in SIM::readState");
      findGrafted(name, 1) -> read();
      continue;
    }

    //grafted with a long name:
    if (IO.compareRecordTag(line,"gH")) {
      name = IO.readUInt32();
      if ( name <= 0 ) throw IOException("invalid Grafted name value in SIM::readState");
      findGrafted(name, 1) -> read();
      continue;
    }

    //microtub:
    if (IO.compareRecordTag(line, "mt")) {
      name = IO.readUInt16();
      if ( name <= 0 ) throw IOException("invalid Microtub name value in SIM::readState");
      findMicrotub(name, 1) -> read();
      continue;
    }

    //Solid:
    if (IO.compareRecordTag(line, "so")) {
      name = IO.readUInt16();
      if ( name <= 0 ) throw IOException("invalid Solid name value in SIM::readState");
      findSolid(name, 1) -> read();
      continue;
    }

    //Nucleus:
    if (IO.compareRecordTag(line, "nu")) {
      name = IO.readUInt16();
      if ( name <= 0 ) throw IOException("invalid Nucleus name value in SIM::readState");
      findNucleus(name, 1) -> read();
      continue;
    }

    //aster:
    if (IO.compareRecordTag(line, "as")) {
      name = IO.readUInt16();
      //name = IO.readUInt32(); nk 6 october 2014
      if ( name <= 0 ) throw IOException("invalid Aster name value in SIM::readState");
      findAster(name, 1) -> read();
      continue;
    }

    //map: CHAITANYA 2008-10-21
	if (IO.compareRecordTag(line, "mp")) {

		//name = IO.readUInt8();
		readMap(); //this reads and sets the dist, res and cat maps


		if ( name <= 0 ) throw IOException("invalid Map name value in SIM::readState");

	    continue;
	 }



    //check for a comment, that includes tags for start and end of frame
    if( line[0] == '#' ) {
      IO.readLine( line+2, sizeof(line)-2 );

      if ( line == strstr(line, FRAME_START_TAG) ) {
        if (( 1 != sscanf( line+strlen(FRAME_START_TAG), "%i", &frame_in_buffer) )
            || ( frame_in_buffer < 0 ))
          throw IOException("invalid frame number in SIM::readState");
      }

      //detect the signal for the end of the frame
      if ( line == strstr(line, FRAME_END_TAG) )
        return NO_ERROR;

      //copy the comment into frame_info[], (play can display it)
      //todo: we should copy all the information, not only the first line
      snprintf( frame_info, sizeof(frame_info), "%s", line+1 );
      continue;
    }

    //simulation 'real' time in seconds
    if (IO.compareRecordTag(line,"ts")) {
      starting_time = IO.readFloatAscii();
      MSG(80, "time = %f\n", starting_time);
      continue;
    }

    //file format of result.out / historical tag
    if (IO.compareRecordTag(line, "ht")) {
      IO.readInputFormat();
      if( IO.getFileFormat() != fileFormatVersion )
        MSG(80, "reading old format: %d", IO.getFileFormat());
      continue;
    }

    //binary signature
    if (IO.compareRecordTag(line, "bs")) {
      IO.readBinarySignature();
      continue;
    }

    //dimension in space
   if (IO.compareRecordTag(line, "di")) {
     IO.readInputDIM();
     if ( IO.getInputDIM() > DIM )
       MSG.warning("SIM::read", "projecting state in lower dimension");
     continue;
   }

    //random seed
    if (IO.compareRecordTag(line, "rs")) {
      IO.skipLine();
      continue;
    }

    //box shape and size, this may overide the settings in data.in
    if (IO.compareRecordTag(line, "sh")) {
      MP.boxshape = IO.readStringInt();
      MSG(80, "boxshape = %i\n", MP.boxshape);
      for( int ii = 0; ii < 6; ++ii )
        MP.boxsize[ ii ] = IO.readFloatAscii();
      setSpace();
      //printf("Z %i %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
      //     MP.boxshape, MP.boxsize[0], MP.boxsize[1], MP.boxsize[2],
      //     MP.boxsize[3], MP.boxsize[4], MP.boxsize[5]);
      continue;
    }

    //an H is a sign of an old format, prior to Feb 1st 2004:
    if (IO.compareRecordTag(line, "H "))
      return 2;
    if (IO.compareRecordTag(line, "D "))
      return 2;

    // we just skip the line
    IO.readLine(line+2, sizeof(line)-2);
    MSG(11, "readState: skipped line: %s\n", line);

  } while( !IO.eof() );

  return 1;
}


//=============================================================================
//=============================OUTPUT TO A FILE================================
//=============================================================================

//prints a short inventory of the current simulation state
void SIM::printShortDescription()
{
  MSG(1,"%i microtubs (<L>=%.2f)\n", nbMicrotub(), meanMTLength());

  MSG(1,"%i solids\n", nbSolid());



  for (int i = 0; i < MP.MAX; ++i)
    if (nbComplexOfType(i))
        MSG(1,"%i complex <%d>: hands <%i><%i>\n", nbComplexOfType(i), i,
            MP.cxhatype1[i], MP.cxhatype2[i]);

  for (int i = 0; i < MP.MAX; ++i)
    if (nbGraftedOfType(i))
      MSG(1,"%i grafted hands <%d>\n", nbGraftedOfType(i), i);
}



//---------------------------------------------------------------------------
//write the current simulation state to IO / output
void SIM::writeStateAuto()
{
  char filename[256], syscmd[256];

  //skipping frame here is usefull when we append to a file
  if ( nb_frames_to_skip > 0 ) {
    --nb_frames_to_skip;
    return;
  }

  try {
    //calling writeState(), handling the files in different ways:
    switch ( IO.getWriteMode() ) {

      case 0: {
        //the file is kept open until the end
        if ( IO.getOutputFile() == 0 )
          IO.openOutputFile();

        writeState();
      } break;

      default:
      case 1: {
        //the file is closed and reopened for every frame

        //we open the file if that is not done,
        //and append to it if frame_in_buffer > 0
        if ( IO.getOutputFile() == 0 )
          IO.openOutputFile( frame_in_buffer > 0 );

        writeState();
        IO.closeOutputFile();
      } break;

      case 2: {
        //Output is done with one frame per file and gziped
        snprintf(filename, sizeof(filename), "result%02d.out", frame_in_buffer);
        IO.openOutputFile(filename);
        writeState();
        IO.closeOutputFile();
        snprintf(syscmd, sizeof(syscmd), "gzip -f %s&", filename);
        system(syscmd);
      } break;
    }
  } catch( IOException e ) {
    MSG.warning("Exception in writeStateAuto", e.getMessage());
  }

  //increment the frame counter, after the write operation:
  ++frame_in_buffer;
}


//---------------------------------------------------------------------------
//write the simulation state to a file 'filename'
void SIM::writeState(const char * filename)
{
    try {
        IO.openOutputFile(filename);
    } catch( IOException e ) {
        MSG.error("SIM::writeState(filename)", "file [%s] not found or erroneous", filename);
    }
    try {
        writeState();
        IO.closeOutputFile();
    } catch( IOException e ) {
        MSG.error("SIM::writeState(filename)", e.getMessage());
    }
}



//---------------------------------------------------------------------------
//write the simulation state for output using the IO file descriptor
void SIM::writeState()
{
    //this throws an exception if the file is incorrect:
    IO.checkOutputFile();

    time_t now = time(0);
    char * date = asctime( localtime(&now ));
    date[19]='\0';

    //frame_in_buffer is initially -1, we set it to a valid value:
    if ( frame_in_buffer < 0 )
      frame_in_buffer = 0;

    //write a line identifying a new frame:
    IO.print("\n\n%s%4i, %9.3f sec. %s", FRAME_START_TAG, frame_in_buffer, simTime(), date);

    //write the message given in MP.info if valid:
    if ( * MP.fileinfo != '\0'  ) {
      IO.writeRecordTag("##");
      IO.writeString(MP.fileinfo);
    }

    //record the historical tag / current file format version:
    IO.writeFileFormatVersion(fileFormatVersion);

    //record the simulated time:
    IO.writeRecordTag("ts");
    IO.writeRealAscii(simTime());

    //record the dimensionality:
    IO.writeRecordTag("di");
    IO.writeIntAscii(DIM);

    //record the binary ID:
    IO.writeBinarySignature();

    //we output the box I.D. and size, in future, we could output its real shape
    IO.writeRecordTag("sh");
    IO.writeIntAscii(MP.boxshape);
    for( int ii = 0; ii < 6; ++ii )
        IO.writeRealAscii( MP.boxsize[ ii ] );

    //record the random seed:
    IO.print("\nrs 0x%-9lx", MP.randseed);

    //the order is such that reading can be straightforward:
    //asters are build from microtubules and solid, and should be after them
    //same for complex, grafteds...

    writeSolid();
    writeMicrotub();
    writeAster();
    writeNucleus();
    writeComplex();
    writeGrafted();
    writeMap();//MAP defined in set_field.cc- showGrad 0=none, 1=res, 2=cat, 3=dist

    IO.writeString("\n");
    IO.writeString(FRAME_END_TAG);
    IO.writeIntAscii(frame_in_buffer);
    IO.writeString("\n\n");
}



