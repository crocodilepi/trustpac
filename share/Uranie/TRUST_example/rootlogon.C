using namespace URANIE::DataServer;
using namespace URANIE::Launcher;
using namespace URANIE::Sampler;
using namespace URANIE::Optimizer;
using namespace URANIE::Modeler;
using namespace URANIE::UncertModeler;
using namespace URANIE::Sensitivity;
using namespace URANIE::Relauncher;
using namespace URANIE::Reoptimizer;
using namespace URANIE::Reliability;
#ifndef WIN32 
using namespace URANIE::XMLProblem;
#endif
#ifdef WITH_MPI
using namespace URANIE::MpiRelauncher;
#endif

void rootlogon()
{

    gStyle->SetPalette(1);
    gStyle->SetOptDate(21);

    //General graphical style
    // Default colors
    int white = 0;
    int color = 30;

    //Legend
    gStyle->SetLegendBorderSize(0);
    gStyle->SetFillStyle(0);

    // Pads
    gStyle->SetPadColor(white);
    gStyle->SetTitleFillColor(white);
    gStyle->SetStatColor(white);

}

/* ==================== Hint ====================
   
   Might be practical to store this in a convenient place (for instance
   your home directory) and to create an alias to make sure that you use
   only one rootlogon file independently of where you are.
  
   example : alias root="root -l ${HOME}/rootlogon.C"
   
   Many style issue can be set once and for all here.
*/
