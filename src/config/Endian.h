////////////////////////////////////////////////////////////////////////////////
// $Id: Endian.h 4446 2013-05-15 19:54:48Z mellinger $
// Description: System endianness.
//
// $BEGIN_BCI2000_LICENSE$
//
// This file is part of BCI2000, a platform for real-time bio-signal research.
// [ Copyright (C) 2000-2012: BCI2000 team and many external contributors ]
//
// BCI2000 is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// BCI2000 is distributed in the hope that it will be useful, but
//                         WITHOUT ANY WARRANTY
// - without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <http://www.gnu.org/licenses/>.
//
// $END_BCI2000_LICENSE$
////////////////////////////////////////////////////////////////////////////////
#ifndef ENDIAN_H
#define ENDIAN_H

#include <sys/param.h>

#ifndef BYTE_ORDER
# define LITTLE_ENDIAN __LITTLE_ENDIAN
# define BIG_ENDIAN __BIG_ENDIAN
# define BYTE_ORDER __BYTE_ORDER
#endif

enum
{
  LittleEndian = LITTLE_ENDIAN,
  BigEndian = BIG_ENDIAN,
  HostOrder = BYTE_ORDER,
};

#endif // ENDIAN_H
