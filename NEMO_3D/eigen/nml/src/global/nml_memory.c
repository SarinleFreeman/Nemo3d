/*****************************************************************************
The <namespace>.
Copyright (C) 2002-2004 California Institute of Technology (Caltech)

This file is part of
The NanoElectronic MOdeling (NEMO) Math Library.

This library is free software which you can redistribute and/or modify
under the terms of the GNU Library General Public License
as published by the Free Software Foundation;
either version 2, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; see the file COPYING. If not, write to the
Free Software Foundation, Inc.,
59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA

Written by E. Robert Tisdale
Modified by Marek J. Korkusinski, October 2004

*****************************************************************************/

#include <nml_global.h>

#ifdef NML_MEM_EXT


#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#endif

#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif


static LONG_8Byte max_active_memory = 0 ;
static LONG_8Byte active_memory = 0 ;
static LONG_8Byte freed_memory  = 0 ;
static LONG_8Byte total_memory  = 0 ;

typedef struct 
{
  LONG_8Byte size ;
  double start ;
} darray ;

static size_t offsetLength = ((size_t) ((char *) &((darray *) 0)->start)) ;



/* IMPLEMENTATION OF THE FUNCTION nml_malloc() */

void * nml_malloc( nml_extent sizeOfBlock ) 
{


  darray * pointerToDarray ;
  /* void * pointerToNewObject ; */
  LONG_8Byte * pointerCopy ;

  /* allocate the memory for the darray structure MINUS one double
   * - this aligns the structure that follows it to the correct
   * word or dword or whatever (i.e., the padding after the "size"
   * is accounted for)
   *
   */

  pointerToDarray = (darray *) malloc( offsetLength + sizeOfBlock ) ;

  /* alternative:
   *
   * pointerToDarray = malloc( sizeof(darray) - sizeof(double) + sizeOfBlock)
   */

  /* check if allocated */

  if ( NULL == pointerToDarray )
    {
      fprintf( stdout, "In nml_malloc: Could not allocate memory!\n" ) ;
      return( NULL ) ;
    }
  else /* memory allocated OK */
    {

      /* WRITE DOWN THE TOTAL SIZE OF THE STRUCTURE */

      pointerCopy = ( LONG_8Byte * ) pointerToDarray ;

      *pointerCopy = sizeOfBlock + offsetLength ;

      active_memory += sizeOfBlock + offsetLength ;
      max_active_memory = max(max_active_memory,active_memory);
      total_memory  += sizeOfBlock + offsetLength ;

#if 0
      if (sizeOfBlock + offsetLength!=40){
         printf("memory malloc %ld bytes  %g MB\n",sizeOfBlock + offsetLength, (sizeOfBlock + offsetLength)/1000000.0);
         nml_memory_report("in mamlloc" );
      }else{
	int j=0;
      }
#endif

      return( (void *) &pointerToDarray->start ) ;

    }

} /* end function nml_malloc() */


/* IMPLEMENTATION OF THE FUNCTION nml_calloc() */

void * nml_calloc( nml_extent noOfElements , nml_extent sizeOfElement ) 
{

  darray * pointerToDarray ;
  /* void * pointerToNewObject ;  */
  LONG_8Byte * pointerCopy ;

  /* allocate the memory for the darray structure MINUS one double
   * - this aligns the structure that follows it to the correct
   * word or dword or whatever (i.e., the padding after the "size"
   * is accounted for)
   *
   */

  pointerToDarray = (darray *) calloc( offsetLength + noOfElements * sizeOfElement, 1 ) ;

  /* alternative:
   *
   * pointerToDarray = calloc( sizeof(darray) - sizeof(double) + 
   *                           noOfElements * sizeOfElement , 1 )
   */

  /* check if allocated */

  if ( NULL == pointerToDarray )
    {
      fprintf( stdout, "In nml_calloc: Could not allocate memory! (%ld x %ld)\n",(long)noOfElements,(long)sizeOfElement) ;
      return( NULL ) ;
    }
  else /* memory allocated OK */
    {

      /* WRITE DOWN THE TOTAL SIZE OF THE STRUCTURE */

      pointerCopy = ( LONG_8Byte * ) pointerToDarray ;

      *pointerCopy = sizeOfElement * noOfElements + offsetLength ;

      active_memory += sizeOfElement * noOfElements + offsetLength ;
      max_active_memory = max(max_active_memory,active_memory);
      total_memory  += sizeOfElement * noOfElements + offsetLength ;

#if 0
      if ((sizeOfElement * noOfElements + offsetLength)!=40){
      printf("memory calloc %ld bytes  %g MB\n",sizeOfElement * noOfElements + offsetLength, (sizeOfElement * noOfElements + offsetLength)/1000000.0);
      nml_memory_report("in calloc" );
      }else{
	int j=0;
      }
#endif

      return( (void *) &pointerToDarray->start ) ;

    }


} /* end function nml_calloc() */


/* IMPLEMENTATION OF THE FUNCTION nml_free() */

void nml_free( void * pointerToObject )
{

  darray * pointerToDarray ;
  LONG_8Byte size ;

  /* check if pointer is NULL */
  if ( NULL == pointerToObject )
    {
      fprintf( stdout , "In nml_free: attempted to free NULL\n" ) ;
      return ;
    }
  else /* the pointer is not NULL */
    {

      /* first get the size of this thing */

      pointerToDarray = (darray *) ((char *) pointerToObject - offsetLength) ;
      size = pointerToDarray->size ;

      active_memory -= size ;
      freed_memory  += size ;

#if 0
      if (size!=40)
      nml_memory_report("in free" );
#endif

      free( pointerToDarray ) ;


      return ;

    }


} /* end function nml_free() */


/* IMPLEMENTATION OF THE FUNCTION nml_memory_report() */

void nml_memory_report( char *tag)
{

  fprintf( stdout , "Memory usage:" ) ;
  fprintf( stdout , "active: %g MB;" , active_memory/1000000.0 ) ;
  fprintf( stdout , "max active: %g MB;" , max_active_memory/1000000.0 ) ;
  fprintf( stdout , "freed: %g MB;" , freed_memory/1000000.0 ) ;
  fprintf( stdout , "total: %g MB;" , total_memory/1000000.0 ) ;
  /**
  fprintf( stdout , "active: %ld MB;" , active_memory ) ;
  fprintf( stdout , "freed: %ld MB;" , freed_memory ) ;
  fprintf( stdout , "total: %ld MB;" , total_memory ) ;
  */
  if (tag){
    fprintf(stdout," %s\n", tag);
  }else
    fprintf(stdout,"\n");
  fflush(stdout);



} /* end function nml_memory_report  */


#else /* NML_MEM_EXT not defined */

    /* map the memory allocations onto the standard calls */

void * nml_malloc( nml_extent size ) 
{
  return( malloc( size ) ) ;
}

void * nml_calloc( nml_extent nmemb , nml_extent size ) 
{
  return( calloc( nmemb , size ) ) ;
}
void nml_free( void * ptr ) 
{
  do { if ( NULL != ptr ) free( ptr ) ; } while( 0 ) ;
  return ;
}

void nml_memory_report( char *tag) 
{

  /* do nothing */
  return ;

}


#endif /* NML_MEM_EXT  */
