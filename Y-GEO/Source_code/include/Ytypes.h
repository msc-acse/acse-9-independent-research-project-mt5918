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
#define YIPROPMAX  1000      /* maximum possible number of propertiest  */
/* NODES */
#define YTN2MEC    -1        /* 2D mechanical x+y d.o.f. node           */
#define YTN2RIG    -2        /* 2D mechanical x+y d.o.f. node           */
#define YTN3MEC   -11        /* 3D mechanical x+y+z d.o.f. node         */

/* ELEMENTS */
#define YTE2TRIELS  1        /* plain stress triangle                   */
#define YTE2TRIRIG  2        /* plain stress triangle                   */
#define YTE2JOINTS  3        /* joint                                   */
#define YTE2TRISOF  4        /* plain stress triangle                   */
#define YTE2PLANESTRESS  5   /* plane stress triangle                   */
#define YTE2PLANESTRAIN  6   /* plane strain triangle                   */

#define YTE2RBAR    7        /* line element Y-RC                       */  
#define YTE2JOINTSCON  8     /* joint concrete Y-RC                     */  

#define YTE3TETELS 11        /* elastic tetrahedra                      */

# define YTEG2RAD5     21   /* 2D 0.5 radius grain                      */
# define YTEGASSEMB    22   /* 2D 0.5 radius grain                      */
# define YTEG2POINT    31   /* 2D 0.5 radius grain                      */

/* FIELDS */
#define YFLDDEF 0           /* default = nothing                        */
#define YFLDSXX 1           /* stress sigma_xx                          */
#define YFLDSXY 2           /* stress sigma_xy                          */
#define YFLDSYY 3           /* stress sigma_yy                          */
#define YFLDSZZ 4           /* stress sigma_zz                          */
#define YFLDSZX 5           /* stress sigma_zx                          */
#define YFLDSZY 6           /* stress sigma_zy                          */
#define YFLDVEL 7           /* velocity                                 */
#define YFLDVEX 8           /* velocity x                               */
#define YFLDVEY 9           /* velocity y                               */
#define YFLDVEZ 10          /* velocity z                               */
#define YFLDAEL 11          /* isotropic elastic damage                 */
#define YFLBORM 12          /* borhole mass                             */
#define YFLBORP 13          /* borhole pressure                         */
#define YFLBORV 14          /* borhole spec. volume                     */
#define YFLEK   15          /* total kinetic energy                     */
#define YFLG2PR 16          /* G2 pressure                              */
#define YFLDJOP 17          /* isotropic elastic damage                 */
#define YFLDFX  18          /* nodal force Fx                           */
#define YFLDFY  19          /* nodal force Fy                           */
#define YFLDFZ  20          /* nodal force Fz                           */
#define YFLEE   21          /* total strain energy                     */
#define YFLDCCX  22         /* nodal coordinate x                       */
#define YFLDCCY  23         /* nodal coordinate y                       */
#define YFLDGRC  24         /* granulometric curve                      */
#define YFLDSBAR 25         /* strss sigma_bar                          */
#define YFLDNCCX 26         /* current nodal coordinate x               */
#define YFLDNCCY 27         /* current nodal coordinate y               */
