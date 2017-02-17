//@HEADER
/*
************************************************************************

              IFPACK: Robust Algebraic Preconditioning Package
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

/* \file Ifpack_HashTable.h
 *
 * \brief HashTable used in Ifpack_ICT and Ifpack_ILUT.
 *
 * \author Marzio Sala, ETHZ/D-INFK.
 *
 * \date Last modified on 30-Jun-06.
 */

#ifndef IFPACK_HASHTABLE_H
#define IFPACK_HASHTABLE_H

#include "Ifpack_ConfigDefs.h"

// ============================================================================
// Hash table with good performances and high level of memory reuse.
// Given a maximum number of keys n_keys, this class allocates chunks of memory
// to store n_keys keys and values. 
//
// Usage:
//
// 1) Instantiate a object,
//
//    Ifpack_HashTable Hash(n_keys);
//
//    n_keys - maximum number of keys (This will be the n_keys with zero 
//             collisons.)
//
// 3) use it, then delete it:
//
//    Hash.get(key, value)       --> returns the value stored on key, or 0.0 
//                                   if not found.
//    Hash.set(key, value)       --> sets the value in the hash table, replace
//                                   existing values.
//    Hash.set(key, value, true) --> to sum into an already inserted value
//    Hash.arrayify(...)
//
// 4) clean memory:
//
//    Hash.reset();
//
// \author Marzio Sala, ETHZ/COLAB
//
// \date 30-Jun-06
// ============================================================================ 

class Ifpack_HashTable 
{
  public:
    //! constructor.
    Ifpack_HashTable(const int n_keys = 1031, const int n_sets = 1)
    {
      n_keys_ = getRecommendedHashSize(n_keys) ;
      n_sets_ = n_sets;
      seed_ = (2654435761U);

      keys_.reserve(50);
      vals_.reserve(50);

      keys_.resize(n_sets_);
      vals_.resize(n_sets_);

      for (int i = 0; i < n_sets_; ++i)
      {
        keys_[i].resize(n_keys_);
        vals_[i].resize(n_keys_);
      }

      counter_.resize(n_keys_);
      for (int i = 0; i < n_keys_; ++i) counter_[i] = 0;
    }

    //! Returns an element from the hash table, or 0.0 if not found.
    inline double get(const int key)
    {
      int hashed_key = doHash(key);

      for (int set_ptr = 0; set_ptr < counter_[hashed_key]; ++set_ptr)
      {
        if (keys_[set_ptr][hashed_key] == key)  
          return(vals_[set_ptr][hashed_key]);
      }

      return(0.0);
    }

    //! Sets an element in the hash table.
    inline void set(const int key, const double value,
                    const bool addToValue = false)
    {
      int hashed_key = doHash(key);
      int& hashed_counter = counter_[hashed_key];

      for (int set_ptr = 0; set_ptr < hashed_counter; ++set_ptr)
      {
        if (keys_[set_ptr][hashed_key] == key)
        {
          if (addToValue)
            vals_[set_ptr][hashed_key] += value;
          else
            vals_[set_ptr][hashed_key] = value;
          return;
        }
      }

      if (hashed_counter < n_sets_)
      {
        keys_[hashed_counter][hashed_key] = key;
        vals_[hashed_counter][hashed_key] = value;
        ++hashed_counter;
        return;
      }

      std::vector<int> new_key;
      std::vector<double> new_val;

      keys_.push_back(new_key);
      vals_.push_back(new_val);
      keys_[n_sets_].resize(n_keys_);
      vals_[n_sets_].resize(n_keys_);

      keys_[n_sets_][hashed_key] = key;
      vals_[n_sets_][hashed_key] = value;
      ++hashed_counter;
      ++n_sets_;
    }

    /*! \brief Resets the entries of the already allocated memory. This
     *  method can be used to clean an array, to be reused without additional
     *  memory allocation/deallocation.
     */
    inline void reset()
    {
      memset(&counter_[0], 0, sizeof(int) * n_keys_);
    }

    //! Returns the number of stored entries.
    inline int getNumEntries() const 
    {
      int n_entries = 0;
      for (int key = 0; key < n_keys_; ++key)
        n_entries += counter_[key];
      return(n_entries);
    }

    //! Converts the contents in array format for both keys and values.
    void arrayify(int* key_array, double* val_array)
    {
      int count = 0;
      for (int key = 0; key < n_keys_; ++key)
        for (int set_ptr = 0; set_ptr < counter_[key]; ++set_ptr)
        {
          key_array[count] = keys_[set_ptr][key];
          val_array[count] = vals_[set_ptr][key];
          ++count;
        }
    }

    //! Basic printing routine.
    void print()
    {
      cout << "n_keys = " << n_keys_ << endl;
      cout << "n_sets = " << n_sets_ << endl;
    }

    int getRecommendedHashSize (int n)
    {
        /* Prime number approximately in the middle of the range [2^x..2^(x+1)]
         * is in primes[x-1]. Every prime number stored is approximately two 
         * times the previous one, so hash table size doubles every time.
         */
        int primes[] = {
        3, 7, 13, 23, 53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593,
        49157, 98317, 196613, 393241, 786433, 1572869, 3145739, 6291469,
        12582917, 25165842, 50331653, 100663319, 201326611, 402653189,
        805306457, 1610612741 } ;
        int i, hsize ;

        /* SRSR : err on the side of performance and choose the next largest 
         * prime number. One can also choose primes[i-1] below to cut the 
         * memory by half.
         */
        hsize = primes[29] ;
        for (i = 6 ; i < 30 ; i++)
        {
            if (n <= primes[i])
            {
                /*hsize = (i == 0 ? n : primes[i-1]) ;*/
                hsize = primes[i] ;
                break ;
            }
        }

        return hsize ;
    }

  private:
    //! Performs the hashing.
    inline int doHash(const int key)
    {
      return (key % n_keys_);
      //return ((seed_ ^ key) % n_keys_);
    }

    int n_keys_;
    int n_sets_;
    std::vector<std::vector<double> > vals_;
    std::vector<std::vector<int> > keys_;
    std::vector<int> counter_;
    unsigned int seed_;
};
#endif
