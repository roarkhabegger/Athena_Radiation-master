//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file frequencygrid.cpp
//  \brief implementation of functions in class Radiation
//======================================================================================


// Athena++ headers
#include "./radiation.hpp"


//--------------------------------------------------------------------------------------
// \!fn void FrequencyGrid()

// \brief function to create the frequency grid

void Radiation::FrequencyGrid()
{
  for(int i=0; i<nfreq; ++i)
    wfreq(i) = 1.0/nfreq;
  
}

