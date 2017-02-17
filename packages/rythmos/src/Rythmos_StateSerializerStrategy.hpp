//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_STATE_SERIALIZER_STRATEGY_H
#define Rythmos_STATE_SERIALIZER_STRATEGY_H

#include "Rythmos_Types.hpp"

#include "Teuchos_Describable.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_ModelEvaluatorBase.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParser.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_StringInputSource.hpp"
#include "Thyra_SpmdMultiVectorSerializer.hpp"
#include <string>

namespace Rythmos {


/** \brief Base class for serializing Rythmos state data
 *
 */
template<class Scalar> 
class StateSerializerStrategy
  : virtual public Teuchos::Describable
{
public:

  virtual void serializeScalar(const Scalar& s, std::ostream& oStream) const = 0;
  virtual void deSerializeScalar(const Ptr<Scalar>& s, std::istream& iStream) const = 0;

  virtual void serializeInt(const int& i, std::ostream& oStream) const = 0;
  virtual void deSerializeInt(const Ptr<int>& i, std::istream& iStream) const = 0;

  virtual void serializeBool(const bool& b, std::ostream& oStream) const = 0;
  virtual void deSerializeBool(const Ptr<bool>& b, std::istream& iStream) const = 0;

  virtual void serializeVectorBase(const VectorBase<Scalar>& vec, std::ostream& oStream) const = 0;
  virtual void deSerializeVectorBase(const Ptr<VectorBase<Scalar> >& vec, std::istream& iStream) const = 0;

  virtual void serializeParameterList(const Teuchos::ParameterList& pl, std::ostream& oStream) const = 0;
  virtual void deSerializeParameterList(const Ptr<Teuchos::ParameterList>& pl, std::istream& iStream) const = 0;

};


template<class Scalar>
class XMLStateSerializerStrategy
  : virtual public StateSerializerStrategy<Scalar>
{
  public:

  XMLStateSerializerStrategy() {}
  virtual ~XMLStateSerializerStrategy() {}

  void serializeScalar(const Scalar& s, std::ostream& oStream) const
  {
    oStream.precision(std::numeric_limits<Scalar>::digits10+4);
    oStream << " " << s << " ";
  }
  void deSerializeScalar(const Ptr<Scalar>& s, std::istream& iStream) const
  {
    TEUCHOS_ASSERT( !Teuchos::is_null(s) );
    iStream >> (*s);
  }

  void serializeInt(const int& i, std::ostream& oStream) const 
  {
    oStream.precision(std::numeric_limits<Scalar>::digits10+4);
    oStream << " " << i << " ";
  }
  void deSerializeInt(const Ptr<int>& i, std::istream& iStream) const 
  {
    TEUCHOS_ASSERT( !Teuchos::is_null(i) );
    iStream >> (*i);
  }

  void serializeBool(const bool& b, std::ostream& oStream) const 
  {
    oStream.precision(std::numeric_limits<Scalar>::digits10+4);
    oStream << " " << b << " ";
  }
  void deSerializeBool(const Ptr<bool>& b, std::istream& iStream) const 
  {
    TEUCHOS_ASSERT( !Teuchos::is_null(b) );
    iStream >> (*b);
  }

  void serializeVectorBase(const VectorBase<Scalar>& vec, std::ostream& oStream) const 
  {
    Thyra::SpmdMultiVectorSerializer<double> vectorSerializer(false); // binaryMode = false
    vectorSerializer.serialize(vec, oStream);
  }
  void deSerializeVectorBase(const Ptr<VectorBase<Scalar> >& vec, std::istream& iStream) const 
  {
    TEUCHOS_ASSERT( !Teuchos::is_null(vec) );
    Thyra::SpmdMultiVectorSerializer<double> vectorSerializer(false); // binaryMode = false
    vectorSerializer.deserialize( iStream, vec.get() );
  }

  void serializeParameterList(const Teuchos::ParameterList& pl, std::ostream& oStream) const 
  {
    Teuchos::XMLParameterListWriter paramWriter;
    Teuchos::XMLObject XMLpl = paramWriter.toXML(pl);
    // Write special key to ostream to mark beginning of parameter list
    oStream << "\nRythmos::StateSerializerStrategy::serializeParameterList begin\n";
    oStream << XMLpl;
    // Write special key to ostream to mark end of parameter list
    oStream << "\nRythmos::StateSerializerStrategy::serializeParameterList end\n";
  }
  void deSerializeParameterList(const Ptr<Teuchos::ParameterList>& pl, std::istream& iStream) const 
  {
    TEUCHOS_ASSERT( !Teuchos::is_null(pl) );
    Teuchos::XMLObject XMLpl;
    std::ostringstream oStringStream;
    // Read in special key from istream to make sure this is a parameter list
    {
      std::string specialKey;
      while (specialKey != "Rythmos::StateSerializerStrategy::serializeParameterList begin" ) {
        std::getline(iStream,specialKey);
        TEUCHOS_ASSERT( !iStream.eof() );
      }
    }
    // Read until special key from istream is found that marks end of parameter list
    while (!iStream.eof()) {
      std::string line;
      std::getline(iStream,line);
      //std::cout << "line = >>" << line << "<<\n";
      if (line == "Rythmos::StateSerializerStrategy::serializeParameterList end") {
        break;
      }
      oStringStream << line;
    }
    Teuchos::StringInputSource src(oStringStream.str());
    Teuchos::XMLParser parser(src.stream());
    XMLpl = parser.parse();
    Teuchos::XMLParameterListReader paramReader;
    pl->setParameters(paramReader.toParameterList(XMLpl));
  }

};

template<class Scalar>
class BinaryStateSerializerStrategy
  : virtual public StateSerializerStrategy<Scalar>
{
  public:

  BinaryStateSerializerStrategy() {}
  virtual ~BinaryStateSerializerStrategy() {}

  void serializeScalar(const Scalar& s, std::ostream& oStream) const
  {
    oStream.precision(std::numeric_limits<Scalar>::digits10+4);
    oStream.write( reinterpret_cast<const char*>(&s), sizeof(Scalar) );

  }
  void deSerializeScalar(const Ptr<Scalar>& s, std::istream& iStream) const
  {
    TEUCHOS_ASSERT( !Teuchos::is_null(s) );
    iStream.read( reinterpret_cast<char*>(&*s), sizeof(Scalar) );
  }

  void serializeInt(const int& i, std::ostream& oStream) const 
  {
    oStream.precision(std::numeric_limits<Scalar>::digits10+4);
    oStream.write( reinterpret_cast<const char*>(&i), sizeof(int) );
  }
  void deSerializeInt(const Ptr<int>& i, std::istream& iStream) const 
  {
    TEUCHOS_ASSERT( !Teuchos::is_null(i) );
    iStream.read( reinterpret_cast<char*>(&*i), sizeof(int) );
  }

  void serializeBool(const bool& b, std::ostream& oStream) const 
  {
    oStream.precision(std::numeric_limits<Scalar>::digits10+4);
    oStream.write( reinterpret_cast<const char*>(&b), sizeof(bool) );
  }
  void deSerializeBool(const Ptr<bool>& b, std::istream& iStream) const 
  {
    TEUCHOS_ASSERT( !Teuchos::is_null(b) );
    iStream.read( reinterpret_cast<char*>(&*b), sizeof(bool) );
  }

  void serializeVectorBase(const VectorBase<Scalar>& vec, std::ostream& oStream) const 
  {
    Thyra::SpmdMultiVectorSerializer<double> vectorSerializer(true); // binaryMode = true
    vectorSerializer.serialize(vec, oStream);
  }
  void deSerializeVectorBase(const Ptr<VectorBase<Scalar> >& vec, std::istream& iStream) const 
  {
    TEUCHOS_ASSERT( !Teuchos::is_null(vec) );
    Thyra::SpmdMultiVectorSerializer<double> vectorSerializer(true); // binaryMode = true
    vectorSerializer.deserialize( iStream, vec.get() );
  }

  void serializeParameterList(const Teuchos::ParameterList& pl, std::ostream& oStream) const 
  {
    Teuchos::XMLParameterListWriter paramWriter;
    Teuchos::XMLObject XMLpl = paramWriter.toXML(pl);
    // Write special key to ostream to mark beginning of parameter list
    oStream << "\nRythmos::StateSerializerStrategy::serializeParameterList begin\n";
    oStream << XMLpl;
    // Write special key to ostream to mark end of parameter list
    oStream << "\nRythmos::StateSerializerStrategy::serializeParameterList end\n";
  }
  void deSerializeParameterList(const Ptr<Teuchos::ParameterList>& pl, std::istream& iStream) const 
  {
    TEUCHOS_ASSERT( !Teuchos::is_null(pl) );
    Teuchos::XMLObject XMLpl;
    std::ostringstream oStringStream;
    // Read in special key from istream to make sure this is a parameter list
    {
      std::string specialKey;
      while (specialKey != "Rythmos::StateSerializerStrategy::serializeParameterList begin" ) {
        std::getline(iStream,specialKey);
        TEUCHOS_ASSERT( !iStream.eof() );
      }
    }
    // Read until special key from istream is found that marks end of parameter list
    while (!iStream.eof()) {
      std::string line;
      std::getline(iStream,line);
      //std::cout << "line = >>" << line << "<<\n";
      if (line == "Rythmos::StateSerializerStrategy::serializeParameterList end") {
        break;
      }
      oStringStream << line;
    }
    Teuchos::StringInputSource src(oStringStream.str());
    Teuchos::XMLParser parser(src.stream());
    XMLpl = parser.parse();
    Teuchos::XMLParameterListReader paramReader;
    pl->setParameters(paramReader.toParameterList(XMLpl));
  }

};

} // namespace Rythmos
#endif // Rythmos_STATE_SERIALIZER_STRATEGY_H
