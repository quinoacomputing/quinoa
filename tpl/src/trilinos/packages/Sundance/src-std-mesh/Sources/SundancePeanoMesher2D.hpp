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


/*
 * SundancePeanoMesher2D.h
 *
 *  Created on: Sep 16, 2009
 *      Author: benk
 */

#ifndef SUNDANCE_PEANOMESHER2D_H_
#define SUNDANCE_PEANOMESHER2D_H_

#include "SundanceDefs.hpp"
#include "SundanceMeshSourceBase.hpp"
#include "SundancePeanoMesh2D.hpp"


namespace Sundance
{
class PeanoMesher2D : public MeshSourceBase
  {
  public:
    /**
     * Set up meshing, for the Peano mesh
     * position_x
     * position_y
     * offset_x
     * offset_y
     * resolution,
     */
	PeanoMesher2D(double position_x, double position_y,
                        double offset_x , double offset_y,
                        double resolution,
                        const MeshType& meshType,
                        const MPIComm& comm = MPIComm::world())
      :
    MeshSourceBase(meshType, 0, comm),
      _position_x(position_x), _position_y(position_y),
      _offset_x(offset_x), _offset_y(offset_y),
      _resolution(resolution) {;}


    /** Create a rectangle mesher from a ParameterList */
	PeanoMesher2D(const ParameterList& params);

    /** */
    virtual ~PeanoMesher2D() {;}

    /** Print a short descriptive std::string */
    virtual std::string description() const
    {return "SundancePeanoMesher[pos x =" + Teuchos::toString(_position_x)
       + ", pos y=" + Teuchos::toString(_position_y)
       + ", offset x=" + Teuchos::toString(_offset_x) +
       + ", offset y=" + Teuchos::toString(_offset_y)
       + ", resolution=" + Teuchos::toString(_resolution)+ "]";}


    /** Return a ref count pointer to self */
    virtual RCP<MeshSourceBase> getRcp() {return rcp(this);}


  protected:

    /** The method which all Mesher should have */
    virtual Mesh fillMesh() const ;

  private:

    /** The X coordinate of the origin point (lower left)*/
    double _position_x;
    /** The Y coordinate of the origin point (lower left)*/
    double _position_y;
    /** The offset (length) of the grid in the X direction*/
    double _offset_x;
    /** The offset (length) of the grid in the Y direction*/
    double _offset_y;
    /** The resolution in each dimension, since we do not want to have stretched elements*/
    double _resolution;

  };
}
#endif /* SUNDANCE_PEANOMESHER2D_H_ */
