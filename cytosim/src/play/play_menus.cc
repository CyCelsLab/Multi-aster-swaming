//RCS: $Id: play_menus.cc,v 2.11 2005/03/15 19:58:03 nedelec Exp $

//-----------------------------------------------------------------------
//--------------------------  basic menus  ------------------------------
//-----------------------------------------------------------------------

void Player::eraseAllMenus() 
{
  for( int ii=0; ii < playMenus.size(); ++ii) {
    if ( playMenus[ii] )
      glutDestroyMenu( playMenus[ii] );
  }
  playMenus.clearAll();
}


//enums for the GLUT menus: their order is irrelevant, they just need to be unique:
enum MENUS_ID {
  MENU_nothing,  MENU_quit,
  
  MENU_mthide,   MENU_mtpoints,  MENU_mtspeckles,
  MENU_mtlines,  MENU_mtratchets,
  MENU_mtend1,   MENU_mtend2,    MENU_mtarrows,    
  MENU_mtforces, MENU_mtproject,
  MENU_mthide2,  MENU_mtfree,    MENU_mtaster1,   
  MENU_mtaster2, MENU_mtcolormode,
  MENU_mtrainbow,
  
  MENU_cxhide,   MENU_cxhands,   MENU_cxlinks,
  MENU_cxdisp,   MENU_cxrainbow,
  MENU_cxhide2,  MENU_cxfree,    MENU_cxbound, MENU_cxlink,
  
  MENU_ghhide,   MENU_ghhands,   MENU_ghlinks, MENU_ghrainbow, MENU_ghflash,
  MENU_ghhide2,  MENU_ghfree,    MENU_ghlink,
  
  MENU_sohide,   MENU_sopoints,  MENU_sodisp, MENU_solinks,
  
  MENU_sim_stop,     MENU_sim_step,    MENU_sim_start, 
  MENU_sim_newstate, MENU_sim_readdata,
  
  MENU_anim_play,    MENU_anim_stop,   MENU_anim_slower, 
  MENU_anim_first,   MENU_anim_prev,   MENU_anim_next,
  
  MENU_view_reset,   MENU_view_scale,  MENU_view_periodic,
  MENU_full_screen,  MENU_show_info,
  MENU_view_style1,  MENU_view_style2, MENU_view_style3
};

//------------------------------------------------------------------------
//poor-man's functions to make 'contextual' menus

void addMenu(const char case_true[], int menu_code)
{
  glutAddMenuEntry( case_true, menu_code );
}

void addMenu(int status, const char case_true[], const char case_false[], int menu_code)
{
  if ( status )
    glutAddMenuEntry( case_true, menu_code );
  else
    glutAddMenuEntry( case_false, menu_code );
}

void addMenu(int status, const char case_true[], int menu_code_true, const char case_false[], int menu_code_false)
{
  if ( status )
    glutAddMenuEntry( case_true, menu_code_true );
  else
    glutAddMenuEntry( case_false, menu_code_false );
}

//------------------------------------------------------------------------
int buildSubMenu1()
{
  //------------------ DISPLAY OPTIONS:
  int menu = glutCreateMenu(Player::processMenuEvents);
  addMenu( PP.mtlines,              "hide lines",          "show lines",           MENU_mtlines );
  addMenu( PP.mtratchets,           "hide ratchets",       "show ratchets",        MENU_mtratchets );
  addMenu( PP.mtpoints,             "hide points",         "show points",          MENU_mtpoints );
  addMenu( PP.mtspeckles,           "hide speckles",       "show speckles",        MENU_mtspeckles );
  addMenu( PP.mtrainbow!=0,         "hide rainbow colors", "show rainbow color",   MENU_mtrainbow );
  addMenu( PP.mtends[MT_MINUS_END], "hide minus-ends",     "show minus-end",       MENU_mtend1 );
  addMenu( PP.mtends[MT_PLUS_END],  "hide plus-ends",      "show plus-end",        MENU_mtend2 );
  addMenu( PP.mtforces,             "hide forces",         "show forces",          MENU_mtforces );
  //addMenu( PP.arrows,             "hide arrows",         "show arrows",          MENU_mtarrows );
  addMenu( PP.mtproject,            "hide projections",    "show projections",     MENU_mtproject );
  addMenu( PP.mtcolormode==1,       "color from plus-end state", "color by asters",   MENU_mtcolormode );
  addMenu("hide all", MENU_mthide );
  return menu;
}

int buildSubMenu2()
{
  int menu = glutCreateMenu(Player::processMenuEvents);
  addMenu( PP.cxhands,              "hide hands",          "show hands",           MENU_cxhands );
  addMenu( PP.cxlinks,              "hide links",          "show links",           MENU_cxlinks );
  addMenu( PP.cxrainbow!=0,         "hide color-forces",   "show color-forces",    MENU_cxrainbow );
  addMenu( PP.cxdisp,               "hide displacements",  "show displacements",   MENU_cxdisp );
  addMenu("hide all", MENU_cxhide );  
  return menu;
}

int buildSubMenu3()
{
  int menu = glutCreateMenu(Player::processMenuEvents);
  addMenu( PP.ghhands,              "hide hands",          "show hands",           MENU_ghhands );
  addMenu( PP.ghlinks,              "hide links",          "show links",           MENU_ghlinks );
  addMenu( PP.ghrainbow!=0,         "hide color-forces",   "show color-forces",    MENU_ghrainbow );
  addMenu( PP.ghflash,              "no flash for attach/detach", "flash for attach/detach", MENU_ghflash );
  addMenu("hide all", MENU_ghhide );  
  return menu;
}

int buildSubMenu4()
{
  int menu = glutCreateMenu(Player::processMenuEvents);
  addMenu( PP.sopoints,             "hide points",          "show points",         MENU_sopoints );
  addMenu( PP.solinks,              "hide links",           "show links",          MENU_solinks );
  addMenu("hide all", MENU_sohide );  
  return menu;
}

// -- display -- menu
int buildMenu1()
{
  int m1 = Player::storeMenu( buildSubMenu1() );
  int m2 = Player::storeMenu( buildSubMenu2() );
  int m3 = Player::storeMenu( buildSubMenu3() );
  int m4 = Player::storeMenu( buildSubMenu4() );
  int menu = glutCreateMenu(Player::processMenuEvents);
  glutAddSubMenu("Tubes",      m1);
  glutAddSubMenu("Complex",    m2);
  glutAddSubMenu("Grafteds",   m3);
  glutAddSubMenu("Solids",     m4);
  return menu;
}


int buildSubMenu5()
{
  int menu = glutCreateMenu(Player::processMenuEvents);
  addMenu("hide all", MENU_mthide2 );
  addMenu( PP.mtselect&1,           "hide free",            "show free",            MENU_mtfree );
  addMenu( PP.mtselect&2,           "hide aster 1",         "show aster 1",         MENU_mtaster1 );
  addMenu( PP.mtselect&4,           "hide aster 2",         "show aster 1",         MENU_mtaster2 );  
  return menu;
}

int buildSubMenu6()
{
  int menu = glutCreateMenu(Player::processMenuEvents);
  addMenu("hide all",                       MENU_cxhide2 );
  addMenu( PP.cxselect&1,           "hide free",            "show free",            MENU_cxfree );
  addMenu( PP.cxselect&2,           "hide bound",           "show bound",           MENU_cxbound );
  addMenu( PP.cxselect&4,           "hide bridge",          "show bridge",          MENU_cxlink );  
  return menu;
}

int buildSubMenu7()
{
  int menu = glutCreateMenu(Player::processMenuEvents);
  addMenu("hide all",                       MENU_ghhide2 );
  addMenu( PP.ghselect&1,           "hide free",            "show free",            MENU_ghfree );
  addMenu( PP.ghselect&2,           "hide bridge",          "show bridge",          MENU_ghlink );  
  return menu;
}

int buildSubMenu8()
{
  int menu = glutCreateMenu(Player::processMenuEvents);
  addMenu("no option here",                 MENU_nothing);
  return menu;
}

// -- select -- menu
int buildMenu2()
{
  int m1 = Player::storeMenu( buildSubMenu5() );
  int m2 = Player::storeMenu( buildSubMenu6() );
  int m3 = Player::storeMenu( buildSubMenu7() );
  int m4 = Player::storeMenu( buildSubMenu8() );
  int menu = glutCreateMenu(Player::processMenuEvents);
  glutAddSubMenu("Tubes",      m1);
  glutAddSubMenu("Complex",    m2);
  glutAddSubMenu("Grafteds",   m3);
  glutAddSubMenu("Solids",     m4);
  return menu;
}

// -- animation -- menu
int buildMenu3()
{
  int menu = glutCreateMenu(Player::processMenuEvents);
  addMenu("new state",                      MENU_sim_newstate );
  addMenu("stop",                           MENU_sim_stop );
  addMenu("step",                           MENU_sim_step );
  addMenu("start",                          MENU_sim_start );
  addMenu("read data.in",                   MENU_sim_readdata );  
  return menu;
}

// -- simulation -- menu
int buildMenu4()
{
  int menu = glutCreateMenu(Player::processMenuEvents);
  addMenu("play/faster",                    MENU_anim_play );
  addMenu("slower",                         MENU_anim_slower );
  addMenu("stop",                           MENU_anim_stop );
  addMenu("-",                              MENU_nothing );
  addMenu("first frame",                    MENU_anim_first );
  addMenu("previous frame",                 MENU_anim_prev );
  addMenu("next frame",                     MENU_anim_next );  
  return menu;
}

int buildSubMenu9()
{
  int menu = glutCreateMenu(Player::processMenuEvents);
  addMenu("1: wireframe",                   MENU_view_style1);
  addMenu("2: -",                           MENU_view_style2);
  addMenu("3: lighting",                    MENU_view_style3);
  return menu;
}

// -- view -- menu
int buildMenu5()
{
  int menu = glutCreateMenu(Player::processMenuEvents);
  addMenu("reset",                          MENU_view_reset);
  addMenu(PP.scalebar, "hide 10 um scale bar", "show 10 um scale bar",  MENU_view_scale);
  addMenu("periodic",                       MENU_view_periodic);
  addMenu("fullscreen",                     MENU_full_screen);
  addMenu( PP.showinfo, "hide info", "show info", MENU_show_info);  
  //glutAddSubMenu("style",             buildSubMenu9()); 
  return menu;
}

// -- final top menu
void Player::buildMenus()
{
  eraseAllMenus();
  
  int m1 = Player::storeMenu( buildMenu1() );
  int m2 = Player::storeMenu( buildMenu2() );
  int m3 = Player::storeMenu( buildMenu3() );
  int m4 = Player::storeMenu( buildMenu4() );
  int m5 = Player::storeMenu( buildMenu5() );
  int menu = glutCreateMenu(Player::processMenuEvents);
  
  glutCreateMenu(processMenuEvents);
  glutAddSubMenu("Display",                          m1);
  glutAddSubMenu("Select",                           m2);
  glutAddSubMenu("Animation",                        m3);
  glutAddSubMenu("Simulation",                       m4);
  glutAddSubMenu("View",                             m5);
  addMenu("Quit",                             MENU_quit);
  
  storeMenu( menu );
  glutAttachMenu(MENU_BUTTON);
}

//------------------------------------------------------------------------
//                    MENU CALLBACK function
//------------------------------------------------------------------------

void Player::processMenuEvents(const int menu) 
{
  //in the switch,
  //  use break; for glutPostRedisplay() to be called
  //  use return; if redrawing is not necessary.
  
  switch (menu) {

  case MENU_nothing:      return;
  case MENU_quit:         exit( EXIT_SUCCESS );

    //------------------------------

  case MENU_mtlines:      PP.mtlines      = !PP.mtlines;          break;
  case MENU_mtpoints:     PP.mtpoints     = !PP.mtpoints;         break;
  case MENU_mtspeckles:   PP.mtspeckles   = !PP.mtspeckles;      break;
  case MENU_mtend1:       PP.mtends[MT_MINUS_END] = !PP.mtends[MT_MINUS_END]; break;
  case MENU_mtend2:       PP.mtends[MT_PLUS_END]  = !PP.mtends[MT_PLUS_END];  break;
  case MENU_mtratchets:   PP.mtratchets   = !PP.mtratchets;       break;
  case MENU_mtarrows:     PP.mtarrows     = !PP.mtarrows;         break;
  case MENU_mtforces:     PP.mtforces     = !PP.mtforces;         break;
  case MENU_mtproject:    PP.mtproject    = !PP.mtproject;        break;
  case MENU_mtrainbow:    PP.mtrainbow    = PP.mtrainbow?0:1;     break;
  case MENU_mtcolormode:  PP.mtcolormode  = PP.mtcolormode==1?2:1;break;
    
  case MENU_mthide:  //disable all MT display options
    PP.mtpoints              = 0;
    PP.mtspeckles            = 0;
    PP.mtlines               = 0;
    PP.mtarrows              = 0;
    PP.mtratchets            = 0;
    PP.mtrainbow             = 0;
    PP.mtends[MT_MINUS_END]  = 0;
    PP.mtends[MT_PLUS_END]   = 0;
    break;

  case MENU_mthide2:      PP.mtselect  = 0;                break;
  case MENU_mtfree:       PP.mtselect ^= 1;                break;
  case MENU_mtaster1:     PP.mtselect ^= 2;                break;
  case MENU_mtaster2:     PP.mtselect ^= 4;                break;

    //------------------------------

  case MENU_cxhide:       PP.cxhands  = 0; PP.cxlinks = 0; break;
  case MENU_cxhands:      PP.cxhands  = ! PP.cxhands;      break;
  case MENU_cxlinks:      PP.cxlinks  = ! PP.cxlinks;      break;
  case MENU_cxdisp:       PP.cxdisp   = ! PP.cxdisp;       break;
  case MENU_cxrainbow:    PP.cxrainbow= ! PP.cxrainbow / MP.hamaxstretch[0]; break;

    //------------------------------

  case MENU_cxhide2:      PP.cxselect  = 0;                break;
  case MENU_cxfree:       PP.cxselect ^= 1;                break;
  case MENU_cxbound:      PP.cxselect ^= 2;                break;
  case MENU_cxlink:       PP.cxselect ^= 4;                break;

    //------------------------------

  case MENU_ghhide:       PP.ghhands = 0; PP.ghlinks = 0;  break;
  case MENU_ghhands:      PP.ghhands = ! PP.ghhands;       break;
  case MENU_ghlinks:      PP.ghlinks = ! PP.ghlinks;       break;
  case MENU_ghflash:      PP.ghflash = ! PP.ghflash;       break;
  case MENU_ghrainbow:    PP.ghrainbow = ! PP.ghrainbow / MP.hamaxstretch[0];  break;

    //------------------------------

  case MENU_ghhide2:      PP.ghselect  = 0;                break; 
  case MENU_ghfree:       PP.ghselect ^= 1;                break;
  case MENU_ghlink:       PP.ghselect ^= 2;                break;

    //------------------------------
    
  case MENU_sohide:       PP.sopoints = 0; PP.solinks = 0; break; 
  case MENU_sopoints:     PP.sopoints = ! PP.sopoints;     break; 
  case MENU_solinks:      PP.solinks  = ! PP.solinks;      break;

    //------------------------------
    
  case MENU_view_reset:   resetView();                     break;
  case MENU_view_scale:   PP.scalebar = !PP.scalebar;      break;
  case MENU_view_periodic:PP.periodic = !PP.periodic;      break;
  case MENU_view_style1:  PP.style = 1;                    break;
  case MENU_view_style2:  PP.style = 2;                    break;
  case MENU_view_style3:  PP.style = 3;                    break;
  case MENU_full_screen:  swapFullScreen();                break;
  case MENU_show_info:    PP.showinfo = !PP.showinfo;      break;
    
    //------------------------------

  case MENU_sim_stop:     processNormalKey('i');           return;
  case MENU_sim_step:     processNormalKey('s');           return;
  case MENU_sim_start:    processNormalKey('a');           return;
  case MENU_sim_newstate: processNormalKey('z');           return;
  case MENU_sim_readdata: processNormalKey('x');           return;
    
    //------------------------------

  case MENU_anim_stop:    processNormalKey('i');           return; 
  case MENU_anim_slower:  processNormalKey('o');           return; 
  case MENU_anim_play:    processNormalKey('p');           return; 
  case MENU_anim_first:   processNormalKey('k');           return; 
  case MENU_anim_prev:    prevFrame();                     return;
  case MENU_anim_next:    nextFrame();                     return;
    
  default: MSG("play: internal error: unknown menu code %i\n", menu);  return;
  }
  glutPostRedisplay();
  buildMenus();
}
