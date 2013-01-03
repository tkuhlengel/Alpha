#ifndef _____RGL_SET_H_____
#define _____RGL_SET_H_____

#ifndef ____RGL___DEFS___H___
#define ____RGL___DEFS___H___

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef __cplusplus 
#define BEGIN_C_DECL extern "C" {
#define END_C_DECL }
#else
#define BEGIN_C_DECL 
#define END_C_DECL 
#endif

BEGIN_C_DECL
typedef void * ptr_t;
typedef unsigned char byte_t;
typedef unsigned int word_t;
typedef char * string_t;
typedef byte_t * chunk_t;
typedef int boolean_t_new;

#define TRUE 1
#define FALSE 0

#define NIL 0x0
typedef int * set_t;

/**
 * An enumeration of set properties.
 */

typedef enum
  {
    SP_MAP = 0x01            /* Adds a map to the set. */
  } set_properties_t;
typedef struct 
  {
    int capacity;
    int size;
    
    int table_size;

    unsigned int properties;

    void * header;
  } set_head_t;
END_C_DECL

#endif







#define get_set_head(s) (set_head_t *) (((byte_t *) s) - sizeof(set_head_t))


/**
 * A set type which can also be used as an array.  The array should not 
 * written into.
 */



#define NO_SUCH_ELEMENT -1

BEGIN_C_DECL

/**
 * Allocates a set object.
 *
 * @param properties a union of set properties masks.
 */

set_t alloc_set(unsigned int properties);

/**
 * Deallocates a set object.
 */

void free_set(set_t s);

/**
 * Returns the number of elements placed in the set.
 */

int size_set(set_t s);

/**
 * Gets the set header (a pointer that can be freely assigned).
 */

void* get_header_set(set_t s);

/**
 * Sets the set header (a pointer that can be freely assigned).
 */

void set_header_set(set_t s, void* header);

/**
 * Puts an element into the set.
 * @param s a set.
 * @param id the 'id' of the element to put.
 * @return the new set - might be different.
 */

set_t put_set(set_t s, int id);

/**
 * Removes an element from the set.
 * @param s a set.
 * @param id the 'id' of the element to remove.
 * @return TRUE if an element was removed or FALSE if nothing removed.
 */

int remove_set(set_t s, int id);

/**
 * Swap element.
 * @param s a set.
 * @param i the first index.
 * @param j the second index.
 */

void swap_set(set_t s, int i, int j);

/**
 * Checks to see if a set contains a given element.
 * @param s a set.
 * @param id the 'id' of the element to be contained.
 * @return TRUE if the element in the set or FALSE otherwise.
 */

int contains_set(set_t s, int id);

/**
 * Returns the index of a given element in the set.
 * @param s a set.
 * @param id the 'id' of the element of a given element.
 * @return the index of the element or NO_SUCH_ELEMENT if no such element exists.
 */

int index_of_set(set_t s, int id);

/**
 * Associates an object with an element.  Requires SP_SET_MAP to be set.
 * @param s a set.
 * @param id the element to associate with.
 * @param obj the object to associate.
 * @return the new set - may be different.
 */

set_t associate_set(set_t s, int id, void* obj);

/**
 * Returns the object mapped to a given element in the set.  Requires SP_SET_MAP to be set.
 * @param s a set.
 * @param id the element to check.
 * @return the object associated with element 'id' or NIL if not possible.
 */
 

void* mapsto_set(set_t s, int id);

/**
 * The object stored at the ith index.  Requires SP_SET_MAP to be set.
 * @param s a set.
 * @param index the index of the element to use.
 * @return the object associated with the indexth element.
 */

void* ith_map_set(set_t s, int index);

//void* ith_map_set(set_t s, int index);

/**
 * Returns an array of objects associated in the set.  Requires SP_SET_MAP to be set.
 * @param s a set.
 * @return an array of objects.
 */

void* * maptable_set(set_t s);

/**
 * Writes a set in human readable.
 * @param fp a file pointer.
 * @param s a set.
 */

//void pretty_print_set(fstream *fp, set_t s);

END_C_DECL

#endif


