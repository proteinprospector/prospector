/* fseek.c */
/* Copyright (c) 2004 by Troels K. */

#include <stdio.h>
#include "fseek.h"
#include "bool.h"

int fseek_calc(ZPOS_T offset, int origin, ZPOS_T* position, ZPOS_T size)
{
   BOOL bOK = TRUE;
   switch (origin)
   {
      case SEEK_SET :
         //bOK = (offset >= 0) && (offset <= size);
         if (bOK) *position = offset;
         break;
      case SEEK_CUR :
         bOK = ((offset + *position) >= 0) && (((offset + *position) <= size));
         if (bOK) *position = offset + *position;
         break;
      case SEEK_END:
         bOK = ((size - offset) >= 0) && (((size - offset) <= size));
         if (bOK) *position = offset + size - 0;
         break;
      default:
         bOK = FALSE;
         break;
   }
   return bOK ? 0 : -1;
}
