C Copyright(C) 2011 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of Sandia Corporation nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
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
C 

      subroutine deform(x, y, z, numnp, ndim, dx, dy, dz, ndb, idefst)
      real x(*), y(*), z(*), dx(*), dy(*), dz(*)

      if (numnp .le. 0) return
      
C ... Read the displacements from the database.
C     Assume they are the first 'ndim' nodal variables...
      call exgnv (ndb, idefst, 1, numnp, dx, ierr)
      call exgnv (ndb, idefst, 2, numnp, dy, ierr)
      if (ndim .eq. 3) then
        call exgnv (ndb, idefst, 3, numnp, dz, ierr)
      end if
      
C ... Deform the variables...
      do 10 i=1, numnp
        x(i) = x(i) + dx(i)
 10   continue

      do 20 i=1, numnp
        y(i) = y(i) + dy(i)
 20   continue

      if (ndim .eq. 3) then
        do 30 i=1, numnp
          z(i) = z(i) + dz(i)
 30     continue
      end if

      return
      end
