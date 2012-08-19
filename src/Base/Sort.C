// -----------------------------------------------------------------------------
// \file    src/Base/Sort.h
// \author  jbakosi
// \date    The Aug 14 9:32:00 2012
// \brief   Shorting functions
// \note    Copyright 2012 Jozsef Bakosi
//          All rights reserved.
// -----------------------------------------------------------------------------

#include <cstdio>
#include "Sort.h"

static void q_sort_i1(int numbers[], int ind[], int left, int right)
{
  int l_hold, r_hold;
  int pivotind;
  int pivot;
  l_hold = left;
  r_hold = right;
  pivot = numbers[left];
  pivotind = ind[left];
  while (left < right)
  {
    while ((numbers[right] >= pivot) && (left < right))
      right--;
    if (left != right)
    {
      numbers[left] = numbers[right];
      ind[left] = ind[right];
      left++;
    }
    while ((numbers[left] <= pivot) && (left < right))
      left++;
    if (left != right)
    {
      numbers[right] = numbers[left];
      ind[right] = ind[left];
      right--;
    }
  }
  numbers[left] = pivot;
  ind[left] = pivotind;
  pivot = left;
  left = l_hold;
  right = r_hold;
  if (left < pivot)
    q_sort_i1(numbers, ind, (int)left, (int)(pivot-1));
  if (right > pivot)
    q_sort_i1(numbers, ind, (int)(pivot+1), (int)(right));
}


void quickSort_i1(int numbers[], int ind[], int array_size)
// quicksorts (int)array numbers and drags array ind along
{
  q_sort_i1(numbers, ind, 0, array_size - 1);
}

static void q_sort_d1(double numbers[], int ind[], int left, int right)
{
  int l_hold, r_hold;
  int pivotind;
  double pivot;
  l_hold = left;
  r_hold = right;
  pivot = numbers[left];
  pivotind = ind[left];
  while (left < right)
  {
    while ((numbers[right] >= pivot) && (left < right))
      right--;
    if (left != right)
    {
      numbers[left] = numbers[right];
      ind[left] = ind[right];
      left++;
    }
    while ((numbers[left] <= pivot) && (left < right))
      left++;
    if (left != right)
    {
      numbers[right] = numbers[left];
      ind[right] = ind[left];
      right--;
    }
  }
  numbers[left] = pivot;
  ind[left] = pivotind;
  pivot = left;
  left = l_hold;
  right = r_hold;
  if (left < pivot)
    q_sort_d1(numbers, ind, (int)left, (int)(pivot-1));
  if (right > pivot)
    q_sort_d1(numbers, ind, (int)(pivot+1), (int)(right));
}


void quickSort_d1(double numbers[], int ind[], int array_size)
// quicksorts array numbers and drags array ind along
{
  q_sort_d1(numbers, ind, 0, array_size - 1);
}

static void q_sort_d3(double numbers[], int ind[], double d1[], double d2[], int left, int right)
{
  int l_hold, r_hold;
  int pivotind;
  double pivot, pd1, pd2;
  l_hold = left;
  r_hold = right;
  pivot = numbers[left];
  pivotind = ind[left];
  pd1 = d1[left];
  pd2 = d2[left];
  while (left < right)
  {
    while ((numbers[right] >= pivot) && (left < right))
      right--;
    if (left != right)
    {
      numbers[left] = numbers[right];
      ind[left] = ind[right];
      d1[left] = d1[right];
      d2[left] = d2[right];
      left++;
    }
    while ((numbers[left] <= pivot) && (left < right))
      left++;
    if (left != right)
    {
      numbers[right] = numbers[left];
      ind[right] = ind[left];
      d1[right] = d1[left];
      d2[right] = d2[left];
      right--;
    }
  }
  numbers[left] = pivot;
  ind[left] = pivotind;
  d1[left] = pd1;
  d2[left] = pd2;
  pivot = left;
  left = l_hold;
  right = r_hold;
  if (left < pivot)
    q_sort_d3(numbers, ind, d1, d2, (int)left, (int)(pivot-1));
  if (right > pivot)
    q_sort_d3(numbers, ind, d1, d2, (int)(pivot+1), (int)(right));
}


void quickSort_d3(double numbers[], int ind[], double d1[], double d2[], int array_size)
// quicksorts array numbers and drags arrays ind, d1, d2 along
{
  q_sort_d3(numbers, ind, d1, d2, 0, array_size - 1);
}
