/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
// ---------------------------------------------------------------------
// Author:     H. Carter Edwards
//
// Purpose:    Associative container for allocated data objects
// ---------------------------------------------------------------------
// Acknowledgements:
//
//   Most all of the algorithms in this class were obtained from
// the Hewlett-Packard source for the Standard Template Library,
// thus the inclusion of Hewlett-Packard's copyright notice.
// Some minor modifications were obtained from Silicon Graphics'
// Standard Template Library source.
// ---------------------------------------------------------------------
/*
 * Copyright (c) 1996,1997
 * Silicon Graphics Computer Systems, Inc.
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Silicon Graphics makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Hewlett-Packard Company makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 */
/*
Red-black tree class, designed for use in implementing STL
associative containers (set, multiset, map, and multimap). The
insertion and deletion algorithms are based on those in Cormen,
Leiserson, and Rivest, Introduction to Algorithms (MIT Press, 1990),
except that

(1) the header cell is maintained with links not only to the root
but also to the leftmost node of the tree, to enable constant time
begin(), and to the rightmost node of the tree, to enable linear time
performance when used with the generic set algorithms (set_union,
etc.);

(2) when a node being deleted has two children its successor node is
relinked into its place, rather than copied, so that the only
iterators invalidated are those referring to the deleted node.
*/
// ---------------------------------------------------------------------

// The header

#include <util/Setv.hpp>

#include <stdexcept>

namespace phdmesh {

// ---------------------------------------------------------------------

Setv<void,void,void> *
Setv<void,void,void>::container( const SetvMember<void> * n )
{
  SetvMember<void> * x = const_cast<SetvMember<void>*>( n );

  if ( x && x->parent ) {
    // Search for root node, while loop breaks at root node
    while ( x != x->parent->parent ) x = x->parent ;
    // the root node's parent is the header
    x = x->parent ;
  }
  else {
    x = NULL ;
  }
  return x ? static_cast<Setv<void,void,void>*>(x) :
             reinterpret_cast<Setv<void,void,void>*>( NULL );
}

SetvMember<void>::~SetvMember()
{
  Setv<void,void,void> * const h = Setv<void,void,void>::container(this);
  if ( h ) { h->remove( this ); }
  parent = left = right = 0 ; color = black ;
}

SetvMember<void>::SetvMember()
  : parent(0), left(0), right(0), color(black)
{}

template<>
SetvMember<void> * setv_iterate<true>( const SetvMember<void> * n )
{
  SetvMember<void> * x = const_cast<SetvMember<void>*>( n );
  SetvMember<void> * y ;
  if ( ( y = x->right ) ) {
    while ( y->left ) y = y->left ;
  }
  else if ( ( y = x->parent ) ) {
    while ( x == y->right ) y = ( x = y )->parent ;
    if ( x == y->parent ) y = y->right ;
  }
  else {
    throw std::logic_error(std::string("ERROR: incrementing Setv::end()"));
  }
  return y ;
}

template<>
SetvMember<void> * setv_iterate<false>( const SetvMember<void> * n )
{
  SetvMember<void> * x = const_cast<SetvMember<void>*>( n );
  SetvMember<void> * y ;
  if ( ( y = x->left ) ) {
    while ( y->right ) y = y->right ;
  }
  else if ( ( y = x->parent ) ) {
    while ( x == y->left ) y = ( x = y )->parent ;
    if ( x == y->parent ) y = y->left ;
  }
  else {
    throw std::logic_error(std::string("ERROR: decrementing Setv::rend()"));
  }
  return y ;
}

// ---------------------------------------------------------------------

inline
void setv_rotate_left( SetvMember<void> ** root ,
                       SetvMember<void> * x )
{
  SetvMember<void> * y = x->right ;

  x->right = y->left ;
  if ( y->left ) y->left->parent = x ;
  y->parent = x->parent ;

  if      ( x == *root )           *root = y ;
  else if ( x == x->parent->left ) x->parent->left  = y ;
  else                             x->parent->right = y ;

  y->left   = x ;
  x->parent = y ;
}

inline
void setv_rotate_right( SetvMember<void> ** root ,
                        SetvMember<void> * x )
{
  SetvMember<void> * y = x->left ;

  x->left = y->right ;
  if ( y->right ) y->right->parent = x;
  y->parent = x->parent ;

  if      ( x == *root )            *root = y ;
  else if ( x == x->parent->right ) x->parent->right = y ;
  else                              x->parent->left  = y ;

  y->right  = x;
  x->parent = y;
}

// ---------------------------------------------------------------------

void Setv<void,void,void>::initialize()
{
  m_header.parent = NULL ;
  m_header.left   = NULL ;
  m_header.right  = NULL ;
  m_header.color  = red ;  /* Color the header node red */

  m_left_end.parent = NULL ;
  m_left_end.left   = NULL ;
  m_left_end.right  = NULL ;
  m_left_end.color  = black ;

  m_right_end.parent = NULL ;
  m_right_end.left   = NULL ;
  m_right_end.right  = NULL ;
  m_right_end.color  = black ;

  m_size = 0 ;

  SetvMember<void> ** leftmost  = & m_left_end.right ;
  SetvMember<void> ** rightmost = & m_right_end.left ;

  *leftmost  = m_header.left  = nREnd(); // left end of the tree
  *rightmost = m_header.right = nEnd();  // right end of the tree
}

Setv<void,void,void>::Setv()
  : SetvMember<void>(),
    m_header( static_cast<SetvMember<void>&>( *this ) ),
    m_left_end(),
    m_right_end(),
    m_size(0)
{
  initialize();
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void Setv<void,void,void>::insert(
  SetvMember<void> * y , SetvMember<void> * z , bool z_lt_y )
{
  SetvMember<void> ** root      = & m_header.parent ;
  SetvMember<void> ** leftmost  = & m_left_end.right ;
  SetvMember<void> ** rightmost = & m_right_end.left ;

  { // Remove from its current container
    Setv<void,void,void> * const zh = container(z);
    if ( zh ) { zh->remove( z ); }
  }

  if ( y == nEnd() ) { // First node inserted
    *root = z ;
    *leftmost = z ;
    *rightmost = z ;
    z->parent = & m_header ; // header is 'super-root'
  }
  else {
    if ( z_lt_y ) {
      y->left = z ;
      // maintain *leftmost pointing to minimum node
      if ( y == *leftmost ) *leftmost = z ;
    }
    else {
      y->right = z;
      // maintain *rightmost pointing to maximum node
      if ( y == *rightmost ) *rightmost = z ;
    }
    z->parent = y ;
  }
  z->left  = 0 ;
  z->right = 0 ;
  z->color = red ;
  ++m_size ;

  // -------------------------------------------------------------------
  // Rebalance, 'y' and 'z' are reused as a local variable

  while ( z != *root && z->parent->color == red ) {
    if ( z->parent == z->parent->parent->left ) {
      y = z->parent->parent->right ;
      if ( y && y->color == red ) {
        z->parent->color         = black;
        y->color                 = black;
        z->parent->parent->color = red;
        z = z->parent->parent ;
      }
      else {
        if ( z == z->parent->right ) {
            z = z->parent ;

            setv_rotate_left( root , z );
            // rotate_left(z);
        }
        z->parent->color         = black;
        z->parent->parent->color = red;

        setv_rotate_right( root , z->parent->parent );
        // rotate_right( z->parent->parent );
      }
    }
    else {
      y = z->parent->parent->left ;
      if ( y && y->color == red ) {
        z->parent->color         = black;
        y->color                 = black;
        z->parent->parent->color = red;
        z = z->parent->parent ;
      }
      else {
        if ( z == z->parent->left ) {
            z = z->parent ;

            setv_rotate_right( root , z );
            // rotate_right(z);
        }
        z->parent->color         = black;
        z->parent->parent->color = red;

        setv_rotate_left( root , z->parent->parent );
        // rotate_left(z->parent->parent);
      }
    }
  }
  (*root)->color = black;
}

// ---------------------------------------------------------------------

void Setv<void,void,void>::remove( SetvMember<void> * node )
{
  static const char method_name[] = "Setv::remove" ;

  if ( container( node ) != this ) {
    std::string msg("Error Setv<*>::remove() not a member");
    throw std::logic_error(msg);
  }

  SetvMember<void> ** root      = & m_header.parent ;
  SetvMember<void> ** leftmost  = & m_left_end.right ;
  SetvMember<void> ** rightmost = & m_right_end.left ;

  if ( 1 == m_size ) { // The last node ?

    if ( node != *leftmost ||
         node != *rightmost ||
         node != *root ) {
      std::string msg(method_name);
      msg.append(" internal data structure corrupted" );
      throw std::runtime_error( msg );
    }

    *leftmost = nREnd();
    *rightmost = nEnd();
    *root = NULL ;
    m_size = 0 ;
    m_header.color = red ;
    node->left = node->right = node->parent = 0 ; node->color = 0 ;
    return ;
  }

  SetvMember<void> * z = node ;
  SetvMember<void> * y = node ;
  SetvMember<void> * x = 0 ;
  SetvMember<void> * x_parent = 0 ;

  // Ready to remove

  if ( y->left == 0 ) {       // z has at most one non-null child. y == z
    x = y->right ;            // x might be null
  }
  else if ( y->right == 0 ) { // z has exactly one non-null child. y == z
    x = y->left ;             // z is not null
  }
  else {                      // z has two non-null children.
     y = y->right ;           // Set y to z's successor.
     while ( y->left ) y = y->left ;
     x = y->right ;           // x might be null
  }

  if ( y != z ) { // relink y in place of z. y is z's successor
    z->left->parent = y ; 
    y->left = z->left ;
    if ( y != z->right ) {
      x_parent = y->parent ;
      if ( x ) x->parent = x_parent ;
      y->parent->left = x;   // y must be a left child
      y->right = z->right;
      z->right->parent = y;
    } else {
      x_parent = y;  // needed in case x == 0
    }
    if ( *root == z) {
      *root = y ;
    }
    else if ( z->parent->left == z) {
      z->parent->left = y;
    }
    else {
      z->parent->right = y;
    }
    y->parent = z->parent;
    { int c = y->color; y->color = z->color; z->color = c ; }
    y = z;
    // y points to node to be actually deleted
  }
  else {  // y == z
    x_parent = y->parent ;
    if ( x ) x->parent = x_parent ; // possibly x == 0
    if ( *root == z) {
      *root = x ;
    }
    else if ( z->parent->left == z ) {
      z->parent->left = x;
    }
    else {
      z->parent->right = x;
    }
    if ( *leftmost == z )  {
      if ( z->right == 0 ) { // z->left must be null also
        // makes *leftmost == nEnd() if z == *root
        *leftmost = z->parent ;
      }
      else {
        // Minimum descending from from 'x'
        SetvMember<void> * xm = x ;
        while ( xm->left ) { xm = xm->left ; }
        *leftmost = xm ;
      }
    }
    if ( *rightmost == z )  {
      if ( z->left == 0 ) { // z->right must be null also
        // makes *rightmost == nEnd() if z == *root
        *rightmost = z->parent ;
      }
      else { // x == z->left
        // Maximum descending from 'x'
        SetvMember<void> * xm = x ;
        while ( xm->right ) { xm = xm->right ; }
        *rightmost = xm ;
      }
    }
  }
  if ( y->color != red ) { 
    while ( x != *root && ( x == 0 || x->color == black ) ) {
      if ( x == x_parent->left ) {
        SetvMember<void> * w = x_parent->right ;
        if ( w->color == red ) {
          w->color        = black;
          x_parent->color = red;

          setv_rotate_left( root , x_parent );
          // rotate_left(x_parent);

          w = x_parent->right ;
        }
        if ((w->left  == 0 || w->left->color  == black) &&
            (w->right == 0 || w->right->color == black)) {
          w->color = red ;
          x = x_parent ;
          x_parent = x_parent->parent ;
        }
        else {
          if (w->right == 0 || w->right->color == black) {
              if ( w->left ) w->left->color = black;
              w->color = red;

              setv_rotate_right( root , w );
              // rotate_right(w);

              w = x_parent->right ;
          }
          w->color = x_parent->color ;
          x_parent->color = black;
          if ( w->right ) w->right->color = black;

          setv_rotate_left( root , x_parent );
          // rotate_left(x_parent);

          break;
        }
      }
      else {  // same as then clause with "right" and "left" exchanged
        SetvMember<void> * w = x_parent->left ;
        if ( w->color == red ) {
          w->color = black;
          x_parent->color = red;

          setv_rotate_right( root , x_parent );
          // rotate_right(x_parent);

          w = x_parent->left ;
        }
        if ((w->right == 0 || w->right->color == black) &&
            (w->left  == 0 || w->left->color  == black)) {
          w->color = red;
          x = x_parent ;
          x_parent = x_parent->parent ;
        }
        else {
          if ( w->left == 0 || w->left->color == black ) {
            if ( w->right ) w->right->color = black;
            w->color = red;

            setv_rotate_left( root , w );
            // rotate_left(w);

            w = x_parent->left ;
          }
          w->color = x_parent->color ;
          x_parent->color = black;
          if ( w->left ) w->left->color = black;

          setv_rotate_right( root , x_parent );
          // rotate_right(x_parent);

          break;
        }
      }
    }
    if ( x ) x->color = black;
  }

  y->left = y->right = y->parent = 0 ; y->color = 0 ;

  --m_size ; // Decrement the tree's count
}

// ---------------------------------------------------------------------
// A reverse communicating method for deleting all entries

SetvMember<void> *
Setv<void,void,void>::unbalancing_removal( SetvMember<void> ** n )
{
  SetvMember<void> * t = *n ;

  while ( t != & m_header && t->parent ) {
    if      ( t->left  ) { t = t->left ; }
    else if ( t->right ) { t = t->right ; }
    else { // Move to parent and remove this leaf
      *n = t->parent ; t->parent = 0 ;
      if ( (*n)->left == t ) (*n)->left  = 0 ;
      else                   (*n)->right = 0 ;
    }
  }

  if ( t == & m_header ) { t = NULL ; }

  return t ;
}

// ---------------------------------------------------------------------
// Virtual destructor

Setv<void,void,void>::~Setv()
{
  if ( m_size || m_header.parent != 0 ) {
    std::string msg("Setv destructor, container is not empty");
    throw std::logic_error( msg );
  }
}

}

