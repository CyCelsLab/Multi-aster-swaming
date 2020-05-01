//RCS: $Id: sim.h,v 2.22 2005/04/27 13:28:56 clausen Exp $
// modified by Chaitanya Athale 08/02/2006
//-------------------------------------sim.h----------------------------------
// Neha Khetan, 3 September 2016:  Output MT stats with bound motors
// Neha Khetan, 11 September 2017: Output statistics of ALL Grafted motors


#ifndef  SIM_H
#define  SIM_H

#include "microtub_list.h"
#include "grafted_list.h"
#include "complex_list.h"
#include "aster_list.h"
#include "solid_list.h"
#include "nucleus_list.h"
#include "iowrapper.h"
#include "scalar_field.h"

///simulator inherits MicrotubList, GraftedList, ComplexList, AsterList, SolidList, NucleusList
class SIM :
  public MicrotubList,
  public GraftedList,
  public ComplexList,
  public AsterList,
  public SolidList,
  public NucleusList
  {
 private:

    ///The counters for iterations
    Step       stepCnt;


    int       haCnt;    			       // Neha Khetan, September 2016
	  real hadtstepValue;             // Neha Khetan, September 2016

    ///real time at stepCnt = 0;
    real       starting_time;

    ///frame number after readState(); or frame number before writeState()
    int        frame_in_buffer;

    ///copy of the last commented line in the frame read by readState()
    char       frame_info[STRING_SIZE];

    ///to skip the next frame in writeState();
    int        nb_frames_to_skip;

    ///IO wrapper for the output of statistics
    IOWrapper  IOStat;
    IOWrapper  IOGhStat;			// Neha Khetan, 3 September 2016: to output statistics of motors


	/// Global variable for updating the counter for motor stepping if hamodel 9 is used
	

 public:

           SIM()                     { resetCounters(); }
  virtual ~SIM()                     { eraseState(); clearSpace(); }

  const Step iterationCnt()   const  { return stepCnt; }



 
 
  // Functions to access real simulated time:
  real      simTime()         const  { return starting_time + MP.dt*stepCnt; }
  real      simStartingTime() const  { return starting_time; }
  real		  updatehadtstepValue() const { return hadtstepValue +  MP.dthastep; }      // Neha Khetan, September 2016
  

  // Functions to get ready to simulate:
  void      resetCounters();
  void      setSpace();
  void      clearSpace();
  ///should be called before the simulation starts, returns an error-code
  int       getReady();
  ///returns true is the simulation is ready to run (or appears so)
  bool      isReady();

  // A step of the simulation engine:
  void      stepAll();

  // accesory function needed for stepAll():
  void      specialEvent();

  int       nbObjects() const;
  int       nbFreeHands() const;


  // Functions to perform the solving step (see simsolve.cc):
  void      solve();
  void      solveAllocate();
  void      solveSetMatrices();
  void      getBlock(const SimObject *, real * );
  int       computePreconditionner();
  void      copyPoints(int);
  void      addTestForces();
  void      dumpAll();


  // Shortcuts to link in the appropriate list directly:
  Microtub *  link(Microtub * x);
  Complex  *  link(Complex * x);
  Grafted  *  link(Grafted * x);
  Aster    *  link(Aster * x);
  Solid    *  link(Solid * x);
  Nucleus  *  link(Nucleus * x);


  // Functions to create the initial state:
  void      pinch1(Microtub*, const real, const Vecteur, const int);
  void      pinch2(Microtub*, const real, const Vecteur, const int);

  void      populateRandom();
  void      populateXenopus();
  void      populateElegans();
  void      populatePombe();
  void      populatePombe2();
  void 		populateXenbeads();
  void      populateChromArray();
  void      populateChrompatt();
  void		populateDNApatt();
  void      populateSpindle();
  void      populateTest1();
  void      populateTest2();
  void      populateTest3();
  void      populateTest4();
  void      populateTest5();
  void      populateState();
  void      populateMouse();
  //CHAITANYA
  //These methods in set_field.cc setup the stabilisation field: called in sim_initial.cc
  void      setstabilField(int superModel);
  real      fieldModel( int model, real dNucl, real maxDist );
  Vecteur   findMinmax( real arr[], int aSize  );
 
  real		  ghMap(real xx, real yy);			      // Neha, 18 October 2014
  Vecteur   ghNucl();                           // Neha 18 october 2014
  
  real		  mtMap(real xx, real yy);			      // Neha, 2 March 2015
  Vecteur   mtNucl();                           // Neha 2 March 2015
  


  //Methods of set_stabilField.cc
  void       writeMap();
  void		   readMap();
  //void		 findMap();


  /// Functions to manipulate the simulation state:
  void      eraseState();
  void      printShortDescription();

  ///the master modulo calling moduloPosition() for all objects
  void      moduloPosition();

  /// Function to read state from a file:
private:
  ///remove objects which do not belong to a certain frame
  void      cleanUpRead(const int frameToKeep);
  ///read State from IO, old format
  int       readState_old();
  ///read State from IO, new format
  int       readState_new();
public:
  ///read sim-state from IO.
  int       readState();
  ///read sim-state from a named file
  int       readState(const char* filename);
  ///the numerical id of the state in the buffer
  int       frameInBuffer()               { return frame_in_buffer; }
  ///an info string read from the file with the sim-state
  char *    frameInfo()                   { return frame_info; }
  ///set the frame numerical id
  void      setFrameInBuffer(const int f) { frame_in_buffer = f; }
  ///increment the frame numerical id
  void      incFrameInBuffer()            { ++frame_in_buffer; }

  ///write sim-state to file specified by IO
  void      writeState();
  ///write sim-state to a named file
  void      writeState(const char* filename);
  ///write sim-state to a file in IO, setting the name appropriately
  void      writeStateAuto();
  ///next call to writeState() will do nothing
  void      skipNextFrame()               { nb_frames_to_skip = 1; }

  //CHAITANYA: I have declared this static here so it's available.
  //           We need to replace this by a FieldList just as with MicrotubList above/.
  //           For now this will do.
  //ScalarField stabilizField;
  //instead i try using map.h// better suited for coords-index-matrix notation inerconversion
  Map<DIM, float> distMap;
  Map<DIM, float> catMap;
  Map<DIM, float> resMap;


  // Maximal number of nuclei
  int maxN;




  /// Function to record statistics on the fly, i.e. at each time step
  void      recordStatistics(bool closeFile = false);
  void      recordStatistics1();
  void      recordStatistics2();
  void      recordStatistics3();
  void      recordStatistics4();
  void      recordStatistics5();
  void      recordStatistics6();
  void      recordStatistics7();
  void 	    recordStatistics8();
  void 	    recordStatistics9();
  void 	    recordStatistics10();

   /// Neha Khetan, 3 September 2016
   /// Function to record Ghstatistics on the fly
   void      recordGhStatistics(bool closeFile = false);  // Neha Khetan, 3 September 2016: to output grafted statistics of motors
   void      recordGhStatistics1();			  // Neha Khetan, 3 September 2016: to output statistics of Only bound motors
   void      recordGhStatistics2();       // Neha Khetan, 11 September 2017: Output statistics of ALL Grafted motors
   void      recordGhStatistics3();       // Neha Khetan, 4 July 2018: Output statistics for all grafted thats bound with the corresponding MT and aster ID
   void      recordGhStatistics4();       // Neha Khetan, 12 Seot, 2018: Output statistics for all grafted thats bound with the corresponding MT and aster ID
};




extern SIM sim;

#endif





