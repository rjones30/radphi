/*
 *  smap.h  --
 *	Define String Table, a String Vector with Hashed lookup
 *	which creates a fixed mapping between a string and an integer
 *
 *  Original: 10-Jan-1996 11:48
 *
 *  Author:   Maarten Ballintijn <Maarten.Ballintijn@cern.ch>
 *
 *  $Id: smap.h,v 1.1 2004/07/01 17:59:21 radphi Exp $
 *
 *  $Log: smap.h,v $
 *  Revision 1.1  2004/07/01 17:59:21  radphi
 *  Makefile
 *      - changed -m486 switch to -march=i486 to make gcc 3.3 happy [rtj]
 *  buildtree.tar
 *      - added a link to piafs in build tree, for clarity [rtj]
 *  npantup.c, qp_execute.c
 *      - fixed a divide-by-zero error in the case where PCNTWK.nchain==0
 *        and the progress dialog attempts to use it as an increment [rtj]
 *  qp_hbook_if.h, piafront.c, cern_types.h, smap.h, piafc.c
 *      - add #ifndef / #endif protections around defines of constants found
 *        in standard POSIX headers, which are being redefined [rtj]
 *
 *  Revision 1.8  2001/09/18 13:55:19  couet
 *  - smap_sort is not used anymore. It has been removed.
 *
 *  Revision 1.7  1996/04/23 18:39:10  maartenb
 *  - Add RCS keywords
 *
 *
 */

#ifndef CERN_SMAP
#define CERN_SMAP


#include	"cern_types.h"
#include	"hash_int_table.h"
#include	"svec.h"

typedef struct _smap_struct_ {
	int		fSize;
	int		fEntries;
	SVec		fV;
	HashIntTable	fT;
} SMapStruct;

typedef SMapStruct	*SMap;


extern SMap
smap_new( int max );

extern SMap
smap_copy( SMap old );

extern int
smap_add( SMap st, String s );

extern int
smap_entries( SMap st );

extern String
smap_get( SMap st, const int i );

extern bool
smap_map( SMap st, String s, int * ip );

extern void
smap_del( SMap st );

#endif	/*	CERN_SMAP	*/
