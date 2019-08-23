//////////////////////////////////////////////////////////////////////////////
//
// File:        PEMFC_ToolBox.cpp
// Directory:   $PEMFC_ROOT/src/
//
//////////////////////////////////////////////////////////////////////////////

#include <PEMFC_ToolBox.h>

// Fill face1face2
//  number of lines = total number of common faces between zone1 and zone2 (including virtual faces)
//  2 columns
//  For each common face face1face2(i,0) = index of face on zone1
//    face1face2(i,1) = index of face on zone2
void PEMFC_ToolBox::find_common_faces_indices( const Zone_VF& z1, const Zone_VF& z2, IntTab& face1face2 )
{
  const DoubleTab& coord1 = z1.xv();
  const DoubleTab& coord2 = z2.xv();
  const int n2 = coord2.dimension(0);
  ArrOfInt correspondance(n2);
  const double epsilon = 1e-10; // tolerance
  Scatter::Chercher_Correspondance(coord1,coord2,correspondance,epsilon);
  face1face2.resize(0,2);
  face1face2.set_smart_resize(1); // allows efficient appending of lines
  for (int i2 = 0; i2 < n2; i2++)
    {
      int i1 = correspondance[i2];
      if (i1 >= 0)
        // face i1 on zone 1 matches face i2 on zone 2
        face1face2.append_line(i1,i2);
    }
}


// Fill index1in2 with:
//  for each face i2 in z2,
//    if face is present in z1 and has index i1 in z1
//       index1in2[i2] = i1
//    else
//       index1in2[i2] = -1
void PEMFC_ToolBox::find_faces_indices1in2( const Zone_VF& z1, const Zone_VF& z2, ArrOfInt& index1in2 )
{
  const DoubleTab coord1 = z1.xv();
  const DoubleTab coord2 = z2.xv();
  const int n2 = coord2.dimension(0);
  index1in2.resize_array(n2);
  const double epsilon = 1e-10; // tolerance
  Scatter::Chercher_Correspondance(coord1,coord2,index1in2,epsilon);
}
