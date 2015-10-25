/* fseek.h */
/* Copyright (c) 2004 by Troels K. */

#ifndef EXTERN_C   
   #ifdef __cplusplus
      #define EXTERN_C    extern "C"
   #else
      #define EXTERN_C    extern
   #endif
#endif

#ifndef ZPOS_T
   #define ZPOS_T unsigned long
#endif

EXTERN_C int fseek_calc(ZPOS_T offset, int origin, ZPOS_T* position, ZPOS_T size);
