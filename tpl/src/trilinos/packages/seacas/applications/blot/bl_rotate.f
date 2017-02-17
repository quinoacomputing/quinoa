C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
C         
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C $Log: rotate.f,v $
C Revision 1.2  2009/03/25 12:36:47  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:10:13  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:56:35  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE ROTATE (NUM, NPROT, ROTMAT, ROTCEN,
     &   XN, YN, ZN, HZ, VT, PD)
C=======================================================================

C   --*** ROTATE *** (MESH) Rotate 3D coordinates
C   --   Written by Amy Gilkey - revised 09/09/87
C   --
C   --ROTATE rotates the 3D coordinates by subtracting the rotation center
C   --and multipling by the rotation matrix.
C   --
C   --Parameters:
C   --   NUM - IN - the number of nodes to rotate
C   --   NPROT - IN - the node numbers of the nodes to rotate
C   --   ROTMAT - IN - the rotation matrix
C   --   ROTCEN - IN - the center of the rotation
C   --   XN, YN, ZN - IN - the original nodal coordinates
C   --   HZ, VT, PD - OUT - the rotated nodal coordinates

      INTEGER NPROT(NUM)
      REAL ROTMAT(3,3), ROTCEN(3)
      REAL XN(*), YN(*), ZN(*)
      REAL HZ(*), VT(*), PD(*)

      DO 100 IX = 1, NUM
         INP = NPROT(IX)
         X = XN(INP) - ROTCEN(1)
         Y = YN(INP) - ROTCEN(2)
         Z = ZN(INP) - ROTCEN(3)
         HZ(INP) = X*ROTMAT(1,1) + Y*ROTMAT(2,1) + Z*ROTMAT(3,1)
         VT(INP) = X*ROTMAT(1,2) + Y*ROTMAT(2,2) + Z*ROTMAT(3,2)
         PD(INP) = X*ROTMAT(1,3) + Y*ROTMAT(2,3) + Z*ROTMAT(3,3)
  100 CONTINUE

      RETURN
      END
