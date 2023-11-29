/*  
 *     SUTIL -- A Stochastic Programming Utility Library
 *
 *     VERSION 0.1
 *
 *     Authors:   Joe Czyzyk
 *                Northwestern University
 *
 *		  Jeff Linderoth and Jierui Shen
 *                Lehigh University
 *
 *     (C)opyright 2005 - J. Czyzyk, J. Linderoth and J. Shen
 *
 * $Id: hash.h,v 1.2 2005/03/04 23:38:22 linderot Exp $
 */

#ifdef SUTIL_DLL_EXPORTS  
#define SUTIL_DLL_API __declspec(dllexport)   
#else  
#define SUTIL_DLL_API __declspec(dllimport)   
#endif 

#ifndef HashFile
#define HashFile

typedef struct node *ListPtr;

typedef struct node {
  int     index;
  char   *entry;
  ListPtr next;
} List;

typedef struct {
  ListPtr *list;
  int      size;
} HashTable;

extern "C" SUTIL_DLL_API HashTable  *NewHashTable(int size);
extern "C" SUTIL_DLL_API void DeleteHash( HashTable * );
extern "C" SUTIL_DLL_API int hash(HashTable *table, char      *string);
extern "C" SUTIL_DLL_API int GetIndex(HashTable  *table, char       *name);
extern "C" SUTIL_DLL_API int Insert(HashTable   *table, char        *name, int          index);
extern "C" SUTIL_DLL_API void PrintHashTable(HashTable   *table);
extern "C" SUTIL_DLL_API void PrintHashTableStats(HashTable   *table);


#endif
