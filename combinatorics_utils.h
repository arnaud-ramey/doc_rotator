#ifndef COMBINATORICS_UTILS_H
#define COMBINATORICS_UTILS_H

#include <string>
#include <vector>
#include <sstream>
#include <set>
#include <algorithm> // for min()
#include<math.h> // for log

//#include "utils/debug/error.h"
//#include "utils/string/StringUtils.h"

#define PI             3.14159265

namespace combinatorics_utils {

typedef std::vector<int> Permutation;
typedef std::vector<int> OrderedCombination;
typedef std::set<int> UnorderedCombination;

/*!
 Populates ans with [0 1 2 ... (end_index-1)]
 \param ans the vector to populate
 \param end_index the number of elements at the end in the vector
*/
static inline void create_listing(Permutation & ans, const int end_index) {
  ans.clear();
  for (int var = 0; var < end_index; ++var)
    ans.push_back(var);
}

////////////////////////////////////////////////////////////////////////////////

/*!
 Implementation of
 http://en.wikipedia.org/wiki/Permutation#Systematic_generation_of_even_permutations
 \param permut the permutation to increase
 \return bool true if it was a success
*/
static inline bool permutation_find_next_lexicographical_ordering
(Permutation & permut) {
  // 1. Find the largest index k such that a[k] < a[k + 1].
  // If no such index exists, the permutation is the last permutation.
  int k = permut.size() - 2;
  Permutation::const_iterator current_pos = permut.end() - 2;
  Permutation::const_iterator current_pos_next = permut.end() - 1;
  while (true) {
    if (*current_pos < *current_pos_next)
      break;
    // go backwards
    --current_pos;
    --current_pos_next;
    --k;
    if (k < 0) // it was the last permutation
      return false;
  }
  //printf("k:%i", k);

  // 2. Find the largest index l such that a[k] < a[l].
  // Since k + 1 is such an index, l is well defined and satisfies k < l.
  int l = permut.size() - 1;
  current_pos_next = permut.end() - 1;
  while (true) {
    if (*current_pos < *current_pos_next)
      break;
    // go backwards
    --current_pos_next;
    --l;
  }
  //printf("l:%i", l);

  // 3. Swap a[k] with a[l].
  std::swap(permut[k], permut[l]);

  // 4. Reverse the sequence from a[k + 1] up to and including the final element a[n].
  std::reverse(permut.begin() + (k + 1), permut.end());

  return true;
}

////////////////////////////////////////////////////////////////////////////////

static inline void combination_init(Permutation & comb, int k, int) {
  create_listing(comb, k);
}

////////////////////////////////////////////////////////////////////////////////

static inline void combination_random(UnorderedCombination & comb, int k, int n) {
  if (k > n)
    printf("Impossible to generate a combination with k=%i > n=%i",
                k, n);

  comb.clear();
  //std::set<int> comb;

  for (int new_value_idx = 0; new_value_idx < k; ++new_value_idx) {
    // find a value that was not used before
    int new_value;
    while (true) {
      new_value = rand() % n;
      if (comb.find(new_value) == comb.end())
        break;
    }

    // add this new value
    comb.insert(new_value);
    //comb.push_back(new_value);
  } // end loop new_value_idx
}

////////////////////////////////////////////////////////////////////////////////

static inline bool combination_incr(OrderedCombination & comb, int k, int n) {
  printf("combination_incr()");
  int i = k - 1;
  ++comb[i];
  while ((i >= 0) && (comb[i] >= n - k + 1 + i)) {
    --i;
    ++comb[i];
  }

  if (comb[0] >= n - k) /* Combination (n-k, n-k+1, ..., n) reached */
    return 0; /* No more combinations can be generated */

  /* comb now looks like (..., x, n, n, n, ..., n).
        Turn it into (..., x, x + 1, x + 2, ...) */
  for (i = i + 1; i < k; ++i)
    comb[i] = comb[i - 1] + 1;

  return 1;
}

////////////////////////////////////////////////////////////////////////////////

/*!
  \example combination_incr_with_position_matters(comb, 2, 3) gives
  (0,1) -> (0,2) -> (1,0) -> (1,2) -> (2,0) -> (2,1)
  */
static inline bool combination_incr_with_position_matters(OrderedCombination & comb, int k, int n) {
  printf("combination_incr_with_position_matters()");

  bool repeat_in_comb;
  int min_index_to_check = k -1;
  do {
    /* increment the vector */
    int index_to_incr = k - 1;
    while(index_to_incr >= 0) {
      ++comb[index_to_incr];
      if (comb[index_to_incr] >= n) {
        comb[index_to_incr] = 0;
        --index_to_incr;
      }
      else // it is over
        break;
    }

    // we have reached the final combination
    if (index_to_incr < 0)
      return false;

    /* check if there is a repeat */
    repeat_in_comb = false;
    if (index_to_incr < min_index_to_check)
      min_index_to_check = index_to_incr;
    for (int changed_idx = min_index_to_check; changed_idx < k; ++changed_idx) {
      for (int unchanged_idx = 0; unchanged_idx < min_index_to_check; ++unchanged_idx) {
        if (comb[changed_idx] == comb[unchanged_idx]) {
          repeat_in_comb = true;
          break;
        }
      } // loop unchanged_idx

      // propagate the break if there was one
      if (repeat_in_comb)
        break;
    } // loop changed_idx
  } while (repeat_in_comb);

  return true;
}

////////////////////////////////////////////////////////////////////////////////

template<class _Container, class _Iterable>
static inline void select_random_sample_with_comb(const _Container & in,
                                                  _Container & out,
                                                  const _Iterable & comb) {
  // add the new eleemnts
  out.clear();
  for(typename _Iterable::const_iterator it = comb.begin();
      it != comb.end() ; ++it)
    out.push_back( in.at( *it ) );
}

////////////////////////////////////////////////////////////////////////////////

template<class _Container>
static inline void select_random_sample(const _Container & in,
                                        _Container & out,
                                        int size) {
  if (size >(int)  in.size())
    printf("Impossible to select a sample of %i elements "
                "> sample size = %i",
                size, in.size());
  // generate a random combination
  UnorderedCombination comb;
  combination_random(comb, size, in.size());
  //printf("comb.size():%u", comb.size());
  // add the new eleemnts
  select_random_sample_with_comb<_Container, UnorderedCombination>
      (in, out, comb);
}

/*!
 *\brief   returns a random number with a gaussian law
 * (standard normal distribution)
 */
inline double rand_gaussian() {
  //  double x1 = 1.f * rand() / RAND_MAX;
  //  double x2 = 1.f * rand() / RAND_MAX;
  //  return sqrt(-2 * log(x1)) * cos(2.f * PI * x2);

  // Box–Muller method
  return sqrt(-2 * log(drand48())) * cos(2.f * PI * drand48());

  // from http://c-faq.com/lib/gaussian.html
  //  static double U, V;
  //  static int phase = 0;
  //  double Z;
  //  if(phase == 0) {
  //    U = (rand() + 1.) / (RAND_MAX + 2.);
  //    V = rand() / (RAND_MAX + 1.);
  //    Z = sqrt(-2 * log(U)) * sin(2 * M_PI * V);
  //  } else
  //    Z = sqrt(-2 * log(U)) * cos(2 * M_PI * V);
  //  phase = 1 - phase;
  //  return Z;
} // end rand_gaussian()



}; // end combinatorics_utils

#endif // COMBINATORICS_UTILS_H
