//////////////////////////////////////////////////////////////////////////////
//
// File:        PEMFC_ToolBox.h
// Directory:   $PEMFC_ROOT/src/
//
//////////////////////////////////////////////////////////////////////////////

#ifndef PEMFC_ToolBox_included
#define PEMFC_ToolBox_included

#include <Objet_U.h>
#include <Zone_VF.h>
#include <Scatter.h>

/////////////////////////////////////////////////////////////////////////////
//
// .NAME        : PEMFC_ToolBox
// .DESCRIPTION : class PEMFC_ToolBox
//
// <Description of class PEMFC_ToolBox>
//
/////////////////////////////////////////////////////////////////////////////


class PEMFC_ToolBox
{

public :

  static void find_common_faces_indices( const Zone_VF& z1, const Zone_VF& z2, IntTab& face1face2 );
  static void find_faces_indices1in2( const Zone_VF& z1, const Zone_VF& z2, ArrOfInt& index1in2 );

};

#endif /* PEMFC_ToolBox_included */
