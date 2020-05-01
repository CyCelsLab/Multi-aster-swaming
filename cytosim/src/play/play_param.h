//RCS: $Id: play_param.h,v 2.21 2005/04/22 11:50:18 nedelec Exp $
//---------------------------------------------------------------------------
//                           play_param.h
//
//          F. Nedelec, Oct 2002, nedelec@embl-heidelberg.de

#ifndef PLAY_PARAM_H
#define PLAY_PARAM_H

#include "parameter_list.h"
#include "sim_param.h"
#include "iowrapper.h"
#include "vecteur3.h"
#include "quaternion.h"

//in play_draw1.cc -> lColor( PP.hacolor[ cx->getType1() ], 1 );   // NK: 31 -03 - 2016
// Added a new param cxcolor
///The list of parameters for Play, the interactive Cytosim
struct ParamPlay : public ParameterList
{
  //=======================PARAMETERS FOR PLAY=======================
public:

  char           filepath[STRING_SIZE];
  char           file[STRING_SIZE];

  int            offscreen;
  int            fullscreen;
  int            buffered;
  int            keyboard;
  int            periodic;
  int            fog;
  int            blend;
  int            style;

  int            frame;
  unsigned int   delay;
  int            live;
  int            dir;
  int            loop;

  int            autozoom;
  real           zoom;
  Vecteur3       focus;
  Quaternion     quat;

  int            width;
  int            height;
  unsigned long  bgcolor;

  char           tag[STRING_SIZE];
  int            showtag;
  int            showinfo;
  int            showframe;
  int            showtime;
  int            showparameters;
  int            scalebar;
  int            showhelp;
  char           flashtext[STRING_SIZE];
  int            showflashtext;   //number of frames the text will be shown

  int			 showGrad;  //show a gradient of stabilization
  real			 alphaDisp; //Alpha transparency level of displayed gradient

  int            boxwidth;
  unsigned long  boxcolor;


  int            mtselect;
  int            mtlines;
  int            mtpoints;
  int            mtspeckles;
  real           mtspeckledist;
  int            mtends[3];
  int            mtarrows;
  int            mtratchets;
  int            mtforces;
  int            mtproject;
  real           mtrainbow;
  unsigned long  mtcolor[10];
  unsigned long  mtfcolor;
  int            mtcolormode;
  real           mtwidth;
  real           mtsize;


  int            cxselect;
  int            cxhands;
  real           cxsize;
  int            cxlinks;
  real           cxrainbow;
  int            cxdisp;
  int            cxflash;
  real           cxwidth;
  unsigned long  cxcolor[10];

  real           hasize;
  unsigned long  hacolor[10];

  int            ghselect;
  int            ghhands;
  int            ghlinks;
  int            ghdisp;
  real           ghrainbow;
  int            ghflash;
  real           ghsize;
  real           ghwidth;
  unsigned long  ghcolor[10];

  int            sopoints;
  int            solinks;
  real           sosize;
  unsigned long  socolor;

  int            nupoints;
  int            nurefpts;
  int            nuenvelope;
  int            nulinks;
  real           nuptsize;
  real           nuenvwidth;
  unsigned long  nuptcolor;
  unsigned long  nurefcolor;
  unsigned long  nuenvcolor;
  unsigned long  nulinkcolor;

  ParamPlay() {

    //---------------------------  options that define program major mode of operation:
    link("offscreen",  &offscreen,  "0", "offscreen rendering (use options -m / -f)");
    link("live",       &live,       "0", "if=1, enter live simulation mode (use option -l)");

    //---------------------------  options to specify the input files
    link("filepath",   filepath,    sizeof(filepath), "\0", "path for directory where files are located");
    link("file",       file,        sizeof(file), RESULT_OUT, "file name for simulation state input");
    link("frame",      &frame,      "0", "index of simulation frame to display initially");

    link("zoom",       &zoom,       "1", "display zoom factor, applies after autozoom if enabled");
    link("autozoom",   &autozoom,   "1", "if=1 autozoom will be done while displaying the first frame");
    link("focus",      (real*)focus,  3, "0 0 0",   "vector specifying the point of focus (X, Y, Z)");
    link("quat",       (real*)quat,   4, "1 0 0 0", "quaternion specifying the modelview rotation");

    link("delay",      &delay,    "100", "time (milli-seconds) between display callbacks");
    link("dir",        &dir,        "0", "direction of play: 1 = forward  or  -1 = backward");
    link("loop",       &loop,       "0", "if=1 loop movie indefinitely during display");

    link("fullscreen", &fullscreen, "0", "if=1, enter fullscreen display mode");
    link("width",      &width,      "768",      "main display window horizontal size");
    link("height",     &height,     "768",      "main display window vertical size");

    link("bgcolor",    &bgcolor,    "0xFFFFFFFF", "background color (hex 0xRGBA, 1 byte per color)");
    link("periodic",   &periodic,   "0", "if=1 produce 9 tiled images for periodic 2D display");

    link("style",      &style,      "1", "drawing style selection (1, 2 or 3)");

    //---------------------------  options that define the display appearance:
    link("fog",        &fog,        "0", "if=1 enable OpenGL FOG based distance cues");
    link("blend",      &blend,      "1", "if=1 enable OpenGL BLEND");
    link("buffered",   &buffered,   "1", "if=1 used double buffering display");


    link("showhelp",   &showhelp,   "0", "display short help on window (press h in live mode)");
    link("tag",        tag,         sizeof(tag), "\0", "string displayed in bottom left of window");
    link("showtag",    &showtag,    "1", "if=1 display the tag string in bottom left corner");
    link("showinfo",   &showinfo,   "0", "if=1 display the frame info in bottom left corner");
    link("showframe",  &showframe,  "1", "if=1 display frame number");
    link("showtime",   &showtime,   "1", "if=1 display time in sec.");
    link("showparameters", &showparameters,   "0", "if=1 display the parameters file");
    link("scalebar",   &scalebar,   "0", "if=1 display a 10 um scale bar");

    link("keyboard",   &keyboard,   "1", "if=0 disable keyboard input");

    link("showGrad",   &showGrad,   "2", "0= show no grad, if=1 f_rescue, if=2 f_catastrophe, if=3 distance");
    link("alphaDisp",  &alphaDisp,   "0.5", "Level of transparency (alpha) in displayed field gradient");


    //---------------------------  display of box edge
    link("boxwidth",   &boxwidth,   "1", "width of line for box edge");
    link("boxcolor",   &boxcolor,   "0x77777777", "color for box edge");

    //---------------------------  display of microtubules
    link("mtselect",   &mtselect,   "255", "MT selection bit-field");
    link("mtlines",    &mtlines,    "1", "if=1 display MT segments");
    link("mtwidth",    &mtwidth,    "2", "MT segments display width");
    link("mtpoints",   &mtpoints,   "0", "if=1 display model MT points");
    link("mtspeckles", &mtspeckles, "0", "if=1 display MT fiduciary marks");
    link("mtspeckledist", &mtspeckledist, "0.2", "distance between MT fiduciary marks");
    link("mtends",     mtends,   3, "0", "if=1 display corresponding MT ends");
    link("mtsize",     &mtsize,     "2", "MT points display size");
    link("mtratchets", &mtratchets, "0", "if=1 display MT polarity ratchets");
    link("mtforces",   &mtforces,   "0", "if=1 display MT forces");
    link("mtarrows",   &mtarrows,   "0", "if=1 display MT arrows");
    link("mtproject",  &mtproject,  "0", "if=1 display projections of MT on box");
    link("mtrainbow",  &mtrainbow,  "0", "if>0 display color-coded MT tensions (Lagrange multipliers)");

    link("mtcolor",    mtcolor,      10, "0x555555FF 0x0066FFFF 0xFF00FFFF 0xF55500FF 0x000051FF", "colors for microtubules");
    link("mtcolormode",&mtcolormode, "0", "set how MT color is choosen: (0):uniform (1):aster (2):state");
    link("mtfcolor",   &mtfcolor,   "0x11FF00FF", "color for MT forces");

    //---------------------------  display of solids
    link("sopoints",   &sopoints,   "0", "if=1 display solid points");
    link("solinks",    &solinks,    "0", "if=1 display solid internal links");
    link("sosize",     &sosize,     "2", "solid points and links size");
    link("socolor",    &socolor,    "0xFF04E3FF", "solid display color");

    //---------------------------  display of nucleus
    link("nupoints",   &nupoints,   "1", "if=1 display nu. points");
    link("nurefpts",   &nurefpts,   "0", "if=1 display reference points");
    link("nuenvelope", &nuenvelope, "1", "if=1 display nu. envelope");
    link("nulinks",    &nulinks,    "1", "if=1 display mtoc-mt links");
    link("nuptsize",   &nuptsize,   "6", "nu. point size");
    link("nuenvwidth", &nuenvwidth, "1", "nu. envelope size");
    link("nuptcolor",  &nuptcolor,  "0x10FF00FF", "nu. point display color");
    link("nurefcolor", &nurefcolor, "0x00CCFFFF", "nu. ref. point display col.");
    link("nuenvcolor", &nuenvcolor, "0x007FFFFF", "nu. envelope display color");
    link("nulinkcolor",&nulinkcolor,"0xFF8000FF", "mtoc-mt link display color");

    //---------------------------  display of complexes
    //---------------------------  display of complexes
    link("cxselect",   &cxselect,   "255", "complexes (cx) selection bit-field");
    link("cxhands",    &cxhands,    "1", "if=1 display cx. hands");
    link("cxsize",     &cxsize,     "3", "cx. hands display size");
    link("cxlinks",    &cxlinks,    "1", "if=1 display cx. links");
    link("cxrainbow",  &cxrainbow,  "0", "if=1 display cx. force");
    link("cxwidth",    &cxwidth,    "3", "cx. links display size");
    link("cxdisp",     &cxdisp,     "0", "if=1 display cx. displacements");
    link("cxflash",    &cxflash,    "0", "draw a flash when motor detaches");
    //link("cxcolor",    cxcolor,   10 , " 0x088A08FF  0xFF0040FF 0x00BB00FF 0x00FF0000 0x111FFFFF  0x0000FFFF  0x001F1FFF" ,"colors for complexes");
        // "cx. colors (unused->hacolor is used)");
    // bluish purple: 0x00BB00FF ; yellow: 0x00BBFFFF ; c yan:0xFF00FBFF ; purple:0xFF0000FF
	link("cxcolor",    cxcolor,   10 , "0xFF0000FF 0x00BB00FF 0xFFDD00FF  0xFF00FFFF 0x00FFFFFF 0x0000FFFF 0x00BB00FF ","colors for complexes");

    //---------------------------   display of grafteds
    link("ghselect",   &ghselect,   "255", "grafted-hands (gh) selection bit-field");
    link("ghhands",    &ghhands,    "1", "if=1 display gh. hands");
    link("ghsize",     &ghsize,     "3", "gh. hands display size");
    link("ghwidth",    &ghwidth,    "3", "gh. hands link display width");
    link("ghlinks",    &ghlinks,    "1", "if=1 display gh. links");
    link("ghrainbow",  &ghrainbow,  "0", "if=1 display gh. force");
    link("ghdisp",     &ghdisp,     "0", "if=1 display gh. displacements");
    link("ghflash",    &ghflash,    "0", "draw a flash when motor detach");
    link("ghcolor",    ghcolor,     10,
         " 0xFF0000FF 0x00BB00FF 0xFF0000FF 0xFF00FFFF 0x0000FFFF 0xFFFF00FF  0x001F1FFF", "gh. display colors (unused->hacolor is used)", PARAM_INVISIBLE);

    link("hasize",     &hasize,     "2", "(not used)", PARAM_INVISIBLE);
    link("hacolor",    hacolor,     10, "0xFF0000FF 0x00BB00FF 0xFFDD00FF  0xFF00FFFF 0x00FFFFFF 0x0000FFFF 0x00BB00FF ", "hands colors");
  }

  virtual ~ParamPlay() { }
};

extern ParamPlay PP;

#endif
