/* Copyright (C) 2000, Dr. Antonio Munjiza
 *
 * This code is provided as part of the book entitled "The Combined
 * Finite Discrete Element Method". It is distributed WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE. Inclusion of a part or whole of this code into any other
 * commercial or research or other purpose code is not granted without author's
 * written explicit permission.
 * When results using whole or any part of this code
 * are published, Y code must be mentioned and acknowledgement to Dr Munjiza must be made.
 * Should you modify this source code, the Copyright (C) on the modified code
 * as a whole belongs to Dr. A. Munjiza regardless of the extent or nature
 * of modifications.
 * Copyright (C) to whole of any code containing any part of this code
 * also belongs to Dr. A.Munjiza. 
 * Any code comprising any part of this source code
 * must be called Y program.
 * If you do not agree with this, you are not allowed to do 
 * any modifications to any part of this source code or included 
 * any part of it in any other program.
 */
/* file Ytypes.h  Y  types of objects */
/* properties */
#define YIPROPMAX  1000 /* maximum possible number of propertiest */

/* NODES */
#define YTN2MEC    -1  /* 2D mechanical x+y d.o.f. node   */
#define YTN2RIG    -2  /* 2D mechanical x+y d.o.f. node   */
#define YTN3MEC   -11  /* 3D mechanical x+y+z d.o.f. node */

/* ELEMENTS */
#define YTE2TRIELS  1 /* plain stress triangle */
#define YTE2TRIRIG  2 /* plain stress triangle */
#define YTE2JOINTS  3 /* joint                 */
#define YTE2TRISOF  4 /* plain stress triangle */
#define YTE2BOUNDS  5 /* boundary              */
#define YTE2AESOUR	6 /* acoustic sources      */
#define YTE3TETELS  11 /* elastic tetrahedra   */
#define YTEG2RAD5   21 /* 2D 0.5 radius grain  */
#define YTEGASSEMB  22 /* 2D 0.5 radius grain  */
#define YTEG2POINT  31 /* 2D 0.5 radius grain  */

/* FIELDS */
#define YFLDDEF 0 /* default = nothing */
#define YFLDSXX 1 /* stress sigma_xx   */
#define YFLDSXY 2 /* stress sigma_xy   */
#define YFLDSYY 3 /* stress sigma_yy   */
#define YFLDSZZ 4 /* stress sigma_yy   */
#define YFLDSZX 5 /* stress sigma_yy   */
#define YFLDSZY 6 /* stress sigma_yy   */
#define YFLDVEL 7 /* velocity          */
#define YFLDVEX 8 /* velocity x        */
#define YFLDVEY 9 /* velocity y        */
#define YFLDVEZ 10 /* velocity z       */
#define YFLDAEL 11 /* isotropic elastic damage       */
#define YFLBORM 12 /* borhole mass                   */
#define YFLBORP 13 /* borhole pressure               */
#define YFLBORV 14 /* borhole spec. volume           */
#define YFLEK   15 /* total kinetic energy           */
#define YFLG2PR 16 /* G2 pressure                    */
#define YFLDJOP 17 /* isotropic elastic damage       */