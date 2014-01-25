//******************************************************************************
/*!
  \file      src/MonteCarlo/LayoutPolicy.h
  \author    J. Bakosi
  \date      Sat 25 Jan 2014 03:14:48 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Particle-, and property-major data layout policies
  \details   Particle-, and property-major data layout policies
*/
//******************************************************************************
#ifndef LayoutPolicy_h
#define LayoutPolicy_h

#include<make_unique.h>

namespace quinoa {

//! Tags for selecting particle-, or property-major data layout policies
const bool ParticleMajor = true;
const bool PropertyMajor = false;

//! Zero-runtime-cost data-layout wrappers for particle properties with
//! type-based compile-time dispatch
template< bool Major >
class ParticleProperties {

  private:
    //! Don't permit copy constructor
    ParticleProperties(const ParticleProperties&) = delete;
    //! Don't permit copy assigment
    ParticleProperties& operator=(const ParticleProperties&) = delete;
    //! Don't permit move constructor
    ParticleProperties(ParticleProperties&&) = delete;
    //! Don't permit move assigment
    ParticleProperties& operator=(ParticleProperties&&) = delete;

   //! Transform a compile-time bool into a type
   template< bool m >
   struct int2type {
     enum { value = m };
   };

   // Overloads for particle-, and property-major data accesses
   inline tk::real& access( int particle, int property, int offset,
                            int2type<ParticleMajor> ) const
   {
     return *(m_ptr.get() + particle*m_nprop + offset + property);
   }
   inline tk::real& access( int particle, int property, int offset,
                            int2type<PropertyMajor> ) const
   {
     return *(m_ptr.get() + particle*m_nprop + offset + property);
   }

   // Overloads for particle-, and property-major const ptr accesses
   inline const tk::real*
   cptr_access( int property, int offset, int2type<ParticleMajor> ) const {
     return m_ptr.get() + offset + property;
   }
   inline const tk::real*
   cptr_access( int property, int offset, int2type<PropertyMajor> ) const {
     return m_ptr.get() + offset + property;
   }

   // Overloads for particle-, and property-major major queries
   inline const char* major( int2type<ParticleMajor> ) const {
     return "particle";
   }
   inline const char* major( int2type<PropertyMajor> ) const {
     return "property";
   }

   const std::unique_ptr< tk::real[] > m_ptr; //!< Particle data pointer
   const int m_nprop;                         //!< Number of particle properties

  public:
    //! Constructor
    ParticleProperties( uint64_t npar, int nprop ) :
      m_ptr( tk::make_unique< tk::real[] >( npar*nprop ) ),
      m_nprop(nprop) {}

    //! Data access dispatch
    inline tk::real&
    operator()( int particle, int property, int offset ) const {
      return access( particle, property, offset, int2type<Major>() );
    }

    //! Const ptr access dispatch
    inline const tk::real*
    cptr( int property, int offset ) const {
      return cptr_access( property, offset, int2type<Major>() );
    }

    //! Major access dispatch
    inline const char* major() const { return major( int2type<Major>() ); }
};

} // quinoa::

// Initial version of zero-runtime-cost data-layout wrappers with type-based
// compile-time dispatch:
//
// Tags for selecting particle-, or property-major data layout policies
// const bool ParticleMajor = true;
// const bool PropertyMajor = false;
// 
// template< bool Major >
// class Data {
// 
//   private:
//    // Transform a compile-time bool into a type
//    template< bool m >
//    struct int2type {
//      enum { value = m };
//    };
// 
//    // Overloads for particle-, and property-major accesses
//    inline
//    tk::real& access( int particle, int property, int2type<ParticleMajor> ) {
//      return *(m_ptr + particle*m_nprop + m_offset + property);
//    }
//    inline
//    tk::real& access( int particle, int property, int2type<PropertyMajor> ) {
//      // This is the same for now, not callsed, irrelevant in zero-cost-test
//      return *(m_ptr + particle*m_nprop + m_offset + property);
//    }
// 
//    tk::real* const m_ptr;
//    const int m_nprop;
//    const int m_offset;
// 
//   public:
//     // Constructor
//     Data( tk::real* const ptr, int nprop, int offset ) :
//       m_ptr(ptr), m_nprop(nprop), m_offset(offset) {}
// 
//     // Access dispatch
//     inline tk::real& operator()( int particle, int property ) {
//       return access( particle, property, int2type<Major>() );
//     }
// };
// 
// Test of zero-cost:
// ------------------
//   * Add to Dirichlet constructor:
//
//       ParticleProperties<Layout> d( particles, m_nprop, m_offset );
//       Model::aa = d( 34, 3 );
//       Model::bb = *(m_particles + 34*m_nprop + m_offset + 3);
//
//   * Add to Model:
//
//     Model {
//       ...
//       public:
//         tk::real aa;
//         tk::real bb;
//       ...
//     }
//
//   * Add to Physics constructor after m_mix is instantiated:
//
//     std::cout << m_mix->aa << m_mix->bb;
//
//   All the above so the optimizing compiler cannot entirely optimize the
//   assignments of aa and bb away.
//
// Debug assembly:
// ---------------
// Generated DEBUG assembly code of the assignments of aa (line 42) and bb (line
// 43) in Dirichlet's constructor, with clang -g -S -mllvm --x86-asm-syntax=intel,
// gnu and intel generate very similar code:
//
//      .loc    143 42 20  # /home/jbakosi/code/quinoa/src/SDE/Dirichlet.h:42:20
// .Ltmp27038:
//      lea     RDI, QWORD PTR [RBP - 56]
//      mov     ESI, 34
//      mov     EDX, 3
//      call    _ZN6quinoa4DataILb1EEclEii
// .Ltmp27039:
//      mov     QWORD PTR [RBP - 176], RAX # 8-byte Spill
//      jmp     .LBB2550_7
// .LBB2550_7:
//      mov     RAX, QWORD PTR [RBP - 176] # 8-byte Reload
//      movsd   XMM0, QWORD PTR [RAX]
//      mov     RCX, QWORD PTR [RBP - 64] # 8-byte Reload
//      movsd   QWORD PTR [RCX + 8], XMM0
//      .loc    143 43 0    # /home/jbakosi/code/quinoa/src/SDE/Dirichlet.h:43:0
//      mov     RDX, QWORD PTR [RCX + 32]
//      imul    ESI, DWORD PTR [RCX + 48], 34
//      movsxd  RDI, ESI
//      shl     RDI, 3
//      add     RDX, RDI
//      movsxd  RDI, DWORD PTR [RCX + 52]
//      movsd   XMM0, QWORD PTR [RDX + 8*RDI + 24]
//      movsd   QWORD PTR [RCX + 16], XMM0
//
// Line 42 translate to register loads and a function call into
// quinoa::ParticleProperties, while line 43 translates to some integer
// arithmetic of the address and loads.
//
// -----------------------------------------------------------------------------
// From the Intel® 64 and IA-32 Architectures Software Developer Manual:
//
// The LEA (load effective address) instruction computes the effective address
// in memory (offset within a segment) of a source operand and places it in a
// general-purpose register. This instruction can interpret any of the
// processor’s addressing modes and can perform any indexing or scaling that may
// be needed. It is especially useful for initializing the ESI or EDI registers
// before the execution of string instructions or for initializing the EBX
// register before an XLAT instruction.
//
// The MOVSXD instruction operates on 64-bit data. It sign-extends a 32-bit
// value to 64 bits. This instruction is not encodable in non-64-bit modes.
//
// A common type of operation on packed integers is the conversion by zero- or
// sign-extension of packed integers into wider data types. SSE4.1 adds 12
// instructions that convert from a smaller packed integer type to a larger
// integer type (PMOVSXBW, PMOVZXBW, PMOVSXBD, PMOVZXBD, PMOVSXWD, PMOVZXWD,
// PMOVSXBQ, PMOVZXBQ, PMOVSXWQ, PMOVZXWQ, PMOVSXDQ, PMOVZXDQ). The source
// operand is from either an XMM register or memory; the destination is an XMM
// register.
//
// IMUL Signed multiply. The IMUL instruction multiplies two signed integer
// operands. The result is computed to twice the size of the source operands;
// however, in some cases the result is truncated to the size of the source
// operands.
//
// SAL/SHL Shift arithmetic left/Shift logical left. The SAL (shift arithmetic
// left), SHL (shift logical left), SAR (shift arithmetic right), SHR (shift
// logical right) instructions perform an arithmetic or logical shift of the
// bits in a byte, word, or doubleword. The SAL and SHL instructions perform the
// same operation. They shift the source operand left by from 1 to 31 bit
// positions. Empty bit positions are cleared. The CF flag is loaded with the
// last bit shifted out of the operand.
//
// Optimized assembly:
// ------------------
// Generated RELWITHDEBINFO assembly code of the assignments of aa (line 42) and
// bb (line 43) in Dirichlet's constructor, with clang -O2 -g DNDEBUG -S -mllvm
// --x86-asm-syntax=intel, gnu and intel generate very similar optimized code:
//
//      .loc    144 42 20  # /home/jbakosi/code/quinoa/src/SDE/Dirichlet.h:42:20
//      movsd   XMM0, QWORD PTR [R14 + 8*RAX + 24]
//      movsd   QWORD PTR [R13 + 8], XMM0
//      .loc    144 43 0    # /home/jbakosi/code/quinoa/src/SDE/Dirichlet.h:43:0
//      mov     RCX, QWORD PTR [R13 + 32]
//      movsd   XMM0, QWORD PTR [RCX + 8*RAX + 24]
//      movsd   QWORD PTR [R13 + 16], XMM0
//
// Both lines 42 and 43 translate to very similar SSE loads with pointer
// arithmetic, i.e., line 42 costs the same as line 43.
//
// -----------------------------------------------------------------------------
// From the Intel® 64 and IA-32 Architectures Software Developer Manual:
//
// The MOVS instruction moves the string element addressed by the ESI register
// to the location addressed by the EDI register. The assembler recognizes three
// “short forms” of this instruction, which specify the size of the string to be
// moved: MOVSB (move byte string), MOVSW (move word string), and MOVSD (move
// doubleword string).
//
// The MOVSD (move scalar double-precision floating-point) instruction transfers
// a 64-bit double-precision floating- point operand from memory to the low
// quadword of an XMM register or vice versa, or between XMM registers.
// Alignment of the memory address is not required, unless alignment checking is
// enabled.

#endif // LayoutPolicy_h
