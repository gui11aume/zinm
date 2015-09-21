// Implementation of a counter of (small) 32 bit integers.
// The implementation is a hash table without hashing. Small numbers
// are stored at the beginning of the array, large numbers could
// be anywhere. Collisions (which happen between a small and a large
// number or between two large numbers) are resolved by quadratic
// probing. The size of the table must be a prime for the probing
// to scan every spot.
//
// A counter is just an array of entries (an 8 byte struct with
// the integer and its count). The counter can hold 'HASH_SIZE'
// entries. Once full, it cannot be incremented.

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#define HASH_SIZE 105943

struct entry_t;
typedef struct entry_t entry_t;

struct entry_t {
   int32_t  key;
   uint32_t val;
};


entry_t *
new_counter
(void)
// Total size of the table: 830 kB.
{

   entry_t *new = calloc(HASH_SIZE, sizeof(entry_t));

   if (new == NULL) {
      fprintf(stderr, "memory error in function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return NULL;
   }

   return new;

}


int
increment
(
    uint32_t itg,
    entry_t *table
)
// Return 1 if the operation succeeds, and 0 otherwise
// (i.e. when the table is full).
{

   for (size_t i = 0 ; i < HASH_SIZE ; i++) {

      // Quadratic probing.
      entry_t e = table[(itg + i*i) % HASH_SIZE];

      if (e.key == itg || e.val == 0) {
         // Found the integer or an empty spot.
         e.key = itg;
         e.val++;
         return 1;
      }

   }

   // The whole table has been scanned and there
   // is nowhere to store this integer. Return the
   // code for failure.
   
   return 0;

}
