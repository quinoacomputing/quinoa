/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#ifndef SUNDANCE_HNMESHER3D_H_
#define SUNDANCE_HNMESHER3D_H_

#include "SundanceDefs.hpp"
#include "SundanceMeshSourceBase.hpp"
#include "SundanceHNMesh3D.hpp"
#include "SundanceRefinementBase.hpp"
#include "SundanceRefinementClass.hpp"
#include "SundanceDomainDefinition.hpp"

namespace Sundance
{
  /** forward declaration */
  class RefinementClass;
  class MeshDomainDef;

class HNMesher3D : public MeshSourceBase
  {
  public:
    /**     */
	HNMesher3D(
			   double position_x , double position_y , double position_z,
               double offset_x , double offset_y, double offset_z,
               int resolution_x , int resolution_y , int resolution_z,
               const MeshType& meshType,
               const MPIComm& comm = MPIComm::world())
      :
    MeshSourceBase(meshType, 0, comm),
      _position_x(position_x), _position_y(position_y), _position_z(position_z),
      _offset_x(offset_x), _offset_y(offset_y), _offset_z(offset_z),
      _resolution_x(resolution_x) , _resolution_y(resolution_y) , _resolution_z(resolution_z) ,
    	refineClass_(dummyRefineClass_) ,
    	meshDomain_(dummyMeshDomain_) {;}

    /**     */
	HNMesher3D(
			   double position_x, double position_y, double position_z,
               double offset_x , double offset_y, double offset_z,
               int resolution_x , int resolution_y , int resolution_z ,
               const MeshType& meshType,
               const RefinementClass& refineClass ,
               const MPIComm& comm = MPIComm::world())
      :
    MeshSourceBase(meshType, 0, comm),
      _position_x(position_x), _position_y(position_y),_position_z(position_z),
      _offset_x(offset_x), _offset_y(offset_y), _offset_z(offset_z),
      _resolution_x(resolution_x) , _resolution_y(resolution_y) ,_resolution_z(resolution_z) ,
    	refineClass_(refineClass) ,
    	meshDomain_(dummyMeshDomain_) {;}

    /**     */
	HNMesher3D(
			   double position_x, double position_y, double position_z,
               double offset_x , double offset_y, double offset_z,
               int resolution_x , int resolution_y , int resolution_z ,
               const MeshType& meshType,
               const MeshDomainDef& meshDomain ,
               const MPIComm& comm = MPIComm::world())
      :
    MeshSourceBase(meshType, 0,comm),
      _position_x(position_x), _position_y(position_y), _position_z(position_z),
      _offset_x(offset_x), _offset_y(offset_y), _offset_z(offset_z),
      _resolution_x(resolution_x) , _resolution_y(resolution_y) , _resolution_z(resolution_z) ,
    	refineClass_(dummyRefineClass_) ,
    	meshDomain_(meshDomain) {;}

    /**     */
	HNMesher3D(
			   double position_x, double position_y,  double position_z,
               double offset_x , double offset_y, double offset_z,
               int resolution_x , int resolution_y , int resolution_z ,
               const MeshType& meshType,
               const RefinementClass& refineClass ,
               const MeshDomainDef& meshDomain ,
               const MPIComm& comm = MPIComm::world())
      :
    MeshSourceBase(meshType, 0, comm),
      _position_x(position_x), _position_y(position_y), _position_z(position_z),
      _offset_x(offset_x), _offset_y(offset_y), _offset_z(offset_z),
      _resolution_x(resolution_x) , _resolution_y(resolution_y) , _resolution_z(resolution_z) ,
    	refineClass_(refineClass) ,
    	meshDomain_(meshDomain) {;}

    /** Create a rectangle mesher from a ParameterList */
	HNMesher3D(const ParameterList& params);

    /** */
    virtual ~HNMesher3D() {;}

    /** Print a short descriptive std::string */
    virtual std::string description() const
    {return "HNMesher3D[pos x =" + Teuchos::toString(_position_x)
       + ", pos y=" + Teuchos::toString(_position_y)
       + ", pos z=" + Teuchos::toString(_position_z)
       + ", offset x=" + Teuchos::toString(_offset_x) +
       + ", offset y=" + Teuchos::toString(_offset_y)
       + ", offset z=" + Teuchos::toString(_offset_z)
       + ", resolution_x=" + Teuchos::toString(_resolution_x)
       + ", resolution_y=" + Teuchos::toString(_resolution_y)
       + ", resolution_z=" + Teuchos::toString(_resolution_z)
       +"]";}


    /** Return a ref count pointer to self */
    virtual RCP<MeshSourceBase> getRcp() {return rcp(this);}


  protected:

    /** The method which all Mesher should have */
    virtual Mesh fillMesh() const ;

  private:

    /** X coordinate of the origin point (lower left)*/
    double _position_x;
    /** Y coordinate of the origin point (lower left)*/
    double _position_y;
    /** Z coordinate of the origin point (lower left)*/
    double _position_z;
    /** offset (length) of the grid in the X direction*/
    double _offset_x;
    /** offset (length) of the grid in the Y direction*/
    double _offset_y;
    /** offset (length) of the grid in the Z direction*/
    double _offset_z;
    /** On the coarse level the resolution on the X axis */
    int _resolution_x;
    /** On the coarse level the resolution on the Y axis */
    int _resolution_y;
    /** On the coarse level the resolution on the Z axis */
    int _resolution_z;

    /** refinement class */
    const RefinementClass refineClass_;

    /** mesh domain (which must not coincide with the whole mesh)*/
    const MeshDomainDef meshDomain_;


    /** static dummy class if the user does not provide refinement class */
    static const RefinementClass dummyRefineClass_;

    /** static domain class if the user does not provide one */
    static const MeshDomainDef dummyMeshDomain_;
  };
}
#endif /* SUNDANCE_HNMESHER3D_H_ */
