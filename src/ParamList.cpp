////////////////////////////////////////////////////////////////////////////////
// $Id: ParamList.cpp 4600 2013-10-09 15:32:54Z mellinger $
// Authors: gschalk@wadsworth.org, juergen.mellinger@uni-tuebingen.de
// Description: ParamList is a BCI2000 type that represents a collection of
//   parameters.
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
#include "PCHIncludes.h"
#pragma hdrstop

#include "ParamList.h"
#include "ParamRef.h"
#include "BCIAssert.h"

#include <sstream>
#include <fstream>
#include <set>
#include <algorithm>

using namespace std;

// **************************************************************************
// Function:   ByName
// Purpose:    Access a parameter by its name.
// Parameters: Parameter name.
// Returns:    Returns a reference to a parameter with a given name.
// **************************************************************************
Param&
ParamList::ByName( const std::string& inName )
{
  ParamContainer::iterator i = mParams.find( inName );
  if( i == mParams.end() )
  {
    mIndex.push_back( &mParams[inName] );
    i = mParams.find( inName );
  }
  return i->second.Param;
}

const Param&
ParamList::ByName( const std::string& inName ) const
{
  static Param defaultParam = Param().SetName( "" ).SetSection( "Default" ).SetType( "int" );
  const Param* result = &defaultParam;
  ParamContainer::const_iterator i = mParams.find( inName );
  if( i != mParams.end() )
    result = &i->second.Param;
  return *result;
}

MutableParamRef
ParamList::operator()( const std::string& inName )
{
  return MutableParamRef( &operator[]( inName ) );
}

ParamRef
ParamList::operator()( const std::string& inName ) const
{
  return ParamRef( &operator[]( inName ) );
}

// **************************************************************************
// Function:   Clear
// Purpose:    deletes all parameters in the parameter list
//             as a result, the list still exists, but does not contain any parameter
// Parameters: N/A
// Returns:    N/A
// **************************************************************************
void
ParamList::Clear()
{
  mIndex.clear();
  mParams.clear();
}

// **************************************************************************
// Function:   Add
// Purpose:    adds a new parameter to the list of parameters
//             if a parameter with the same name already exists,
//             it updates the currently stored parameter with the provided values
// Parameters: inParam       - reference to a Param object representing the
//                             parameter,
//             inSortingHint - float value used as an additional sorting
//                             criterion.
// Returns:    N/A
// **************************************************************************
void
ParamList::Add( const Param& inParam, float inSortingHint )
{
  ByName( inParam.Name() );
  ParamEntry& entry = mParams[inParam.Name()];
  entry.Param = inParam;
  entry.SortingHint = inSortingHint;
}

// **************************************************************************
// Function:   Add
// Purpose:    adds a new parameter to the list of parameters
//             if a parameter with the same name already exists,
//             it updates the currently stored parameter with the provided values
// Parameters: inLine - ASCII string, as defined in project description,
//                           defining this parameter
// Returns:    true if inLine is a correct parameter line, false otherwise
// **************************************************************************
bool
ParamList::Add( const string& inLine )
{
  istringstream linestream( inLine );
  Param param;
  if( linestream >> param )
    ByName( param.Name() ) = param;
  return static_cast<bool>(linestream);
}

// **************************************************************************
// Function:   Delete
// Purpose:    deletes a parameter of a given name in the list of parameters
//             it also frees the memory for this particular parameter
//             it does not do anything, if the parameter does not exist
// Parameters: name - name of the parameter
// Returns:    N/A
// **************************************************************************
void
ParamList::Delete( const std::string& inName )
{
  ParamContainer::iterator i = mParams.find( inName );
  if( i != mParams.end() )
  {
    Index::iterator j = mIndex.begin();
    while( j != mIndex.end() && *j != &i->second )
      ++j;
    bciassert( j != mIndex.end() );
    mIndex.erase( j );
    mParams.erase( i );
  }
}

// **************************************************************************
// Function:   WriteToStream
// Purpose:    Member function for formatted stream output of the entire
//             parameter list.
//             For partial output, use another instance of type ParamList
//             to hold the desired subset as in ParamList::SaveParameterList().
//             All formatted output functions are, for consistency's sake,
//             supposed to use this function.
// Parameters: Output stream to write into.
// Returns:    Output stream.
// **************************************************************************
ostream&
ParamList::WriteToStream( ostream& os ) const
{
  for( int i = 0; i < Size(); ++i )
    os << ByIndex( i ) << '\n';
  return os;
}

// **************************************************************************
// Function:   ReadFromStream
// Purpose:    Member function for formatted stream input of the entire
//             parameter list. The list is cleared before reading.
//             For partial input, use another instance of type ParamList
//             to hold the desired subset as in ParamList::LoadParameterList().
//             All formatted input functions are, for consistency's sake,
//             supposed to use this function.
// Parameters: Input stream to read from.
// Returns:    Input stream.
// **************************************************************************
istream&
ParamList::ReadFromStream( istream& is )
{
  Clear();
  Param param;
  is >> ws;
  while( !is.eof() && is >> param >> ws )
    ByName( param.Name() ) = param;
  return is;
}

// **************************************************************************
// Function:   WriteBinary
// Purpose:    Member function for binary stream output of the entire
//             parameter list.
//             For partial output, use another instance of type ParamList
//             to hold the desired subset as in ParamList::SaveParameterList().
// Parameters: Output stream to write into.
// Returns:    Output stream written into.
// **************************************************************************
ostream&
ParamList::WriteBinary( ostream& os ) const
{
  for( int i = 0; i < Size(); ++i )
    ByIndex( i ).WriteBinary( os );
  return os;
}

// **************************************************************************
// Function:   ReadBinary
// Purpose:    Member function for binary stream input of the entire
//             parameter list. The list is cleared before reading.
//             For partial input, use another instance of type ParamList
//             to hold the desired subset as in ParamList::LoadParameterList().
// Parameters: Input stream to read from.
// Returns:    N/A
// **************************************************************************
istream&
ParamList::ReadBinary( istream& is )
{
  return ReadFromStream( is );
}

// **************************************************************************
// Function:   Save
// Purpose:    Saves the current list of paramters in a parameter file
// Parameters: char *filename - filename to save the list to
//             usetags - if usetags is true, then the "tag" value in each parameter
//                       determines whether the parameter should be saved
//                       if the tag value in the parameter is true, then the
//                       parameter will NOT be saved
//                       if usetags is false, then all parameters are saved
// Returns:    true - successful
//             false - error (disk full, etc.)
// **************************************************************************
bool
ParamList::Save( const string& inFileName ) const
{
  ofstream file( inFileName.c_str() );
  if( !file.is_open() )
    return false;
  file << *this;
  return !file.fail();
}

// **************************************************************************
// Function:   LoadParameterList
// Purpose:    Loads the current list of parameters from a parameter file
//             It does NOT load system critical dynamic parameters (e.g., ports,
//             IP addresses)
// Parameters: char *filename - filename of the parameterlist
//             nonexisting - if true, load parameters, even if they currently do
//                       not exist in the list
// Returns:    true - successful
//             false - error
// **************************************************************************
bool
ParamList::Load( const string& inFileName, bool inImportNonexisting )
{
  ifstream file( inFileName.c_str() );
  ParamList paramsFromFile;
  file >> paramsFromFile;
  if( file.fail() )
    return false;

  // If desired, exclude parameters missing from the main parameter list.
  typedef set<string> NameSet;
  NameSet unwantedParams;
  if( !inImportNonexisting )
    for( int i = 0; i < paramsFromFile.Size(); ++i )
      if( !Exists( paramsFromFile[i].Name() ) )
        unwantedParams.insert( paramsFromFile[i].Name() );

  for( NameSet::const_iterator i = unwantedParams.begin(); i != unwantedParams.end(); ++i )
    paramsFromFile.Delete( *i );

  for( int i = 0; i < paramsFromFile.Size(); ++i )
  {
    Param &p = paramsFromFile[i], &q = ByName( p.Name() );
    if( !q.Readonly() )
      q.AssignValues( p );
  }
  return true;
}

// **************************************************************************
// Function:   Sort
// Purpose:    Sorts parameters by their section fields, and sorting hints.
// Parameters: N/A
// Returns:    N/A
// **************************************************************************
void
ParamList::Sort()
{
  // stable_sort will retain the relative order of parameters with identical
  // sections.
  stable_sort( mIndex.begin(), mIndex.end(), ParamEntry::Compare );
}
