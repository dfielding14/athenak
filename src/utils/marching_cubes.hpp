#ifndef UTILS_MARCHING_CUBES_HPP_
#define UTILS_MARCHING_CUBES_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file marching_cubes.hpp
//!       Implelmented by Lachlan Lancaster, working off of the Cython implementation in
//!       the scikit-image library (https://github.com/scikit-image/scikit-image/blob/v0.23.2/skimage/measure/_marching_cubes_lewiner.py),
//!       which itself is based on the Thomas Lewiner et al. 2003 implementation
//!       (http://thomas.lewiner.org/pdfs/marching_cubes_jgt.pdf).
//!       The method is topologically correct at all spacings

#include <cmath>
#include <iostream>
#include <list>

#include "mc_luts.hpp"
#include "athena.hpp"
#include "athena_tensor.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"

struct Cube {
  Cube();
  KOKKOS_INLINE_FUNCTION
  Cube(Real v0, Real v1, Real v2, Real v3, Real v4, Real v5, Real v6, Real v7);
  KOKKOS_INLINE_FUNCTION
  ~Cube();

  // Values of cube corners (isovalue subtracted)
  Real v0, v1, v2, v3, v4, v5, v6, v7;
};

Cube::Cube():
         v0(0.0), v1(0.0), v2(0.0), v3(0.0),
         v4(0.0), v5(0.0), v6(0.0), v7(0.0){ }

KOKKOS_INLINE_FUNCTION
Cube::Cube(Real v0, Real v1, Real v2, Real v3, Real v4, Real v5, Real v6, Real v7):
         v0(v0), v1(v1), v2(v2), v3(v3),
         v4(v4), v5(v5), v6(v6), v7(v7){ }

KOKKOS_INLINE_FUNCTION
Cube::~Cube() {
  // Destructor
}

KOKKOS_INLINE_FUNCTION
bool test_face(signed char face, Cube cube) {
  Real A,B,C,D ;

  switch( face ) {
    case -1 :
    case  1 :
      A = cube.v0; B = cube.v4;
      C = cube.v5; D = cube.v1;
      break;
    case -2 :
    case  2 :
      A = cube.v1; B = cube.v5;
      C = cube.v6; D = cube.v2;
      break;
    case -3 :
    case  3 :
      A = cube.v2; B = cube.v6;
      C = cube.v7; D = cube.v3;
      break;
    case -4 :
    case  4 :
      A = cube.v3; B = cube.v7;
      C = cube.v4; D = cube.v0;
      break;
    case -5 :
    case  5 :
      A = cube.v0; B = cube.v3;
      C = cube.v2; D = cube.v1;
      break;
    case -6 :
    case  6 :
      A = cube.v4; B = cube.v7;
      C = cube.v6; D = cube.v5;
      break;
    default:
    A = B = C = D = 0 ;
  };

  if( fabs( A*C - B*D ) < FLT_EPSILON )
    return face >= 0 ;
  return face * A * ( A*C - B*D ) >= 0  ;  // face and A invert signs
}

KOKKOS_INLINE_FUNCTION
bool test_interior(signed char s, Cube cube, signed char _case,
                   signed char _config, signed char _subconfig) {
  Real t, At=0, Bt=0, Ct=0, Dt=0, a, b ;
  char  test =  0;
  char  edge = -1; // reference edge of the triangulation

  switch( _case )
  {
  case  4 :
  case 10 :
    a = ( cube.v4 - cube.v0 ) * ( cube.v6 - cube.v2 ) - ( cube.v7 - cube.v3 ) * ( cube.v5 - cube.v1 ) ;
    b =  cube.v2 * ( cube.v4 - cube.v0 ) + cube.v0 * ( cube.v6 - cube.v2 )
             - cube.v1 * ( cube.v7 - cube.v3 ) - cube.v3 * ( cube.v5 - cube.v1 ) ;
    t = - b / (2*a) ;
    if( t<0 || t>1 ) return s>0 ;

    At = cube.v0 + ( cube.v4 - cube.v0 ) * t ;
    Bt = cube.v3 + ( cube.v7 - cube.v3 ) * t ;
    Ct = cube.v2 + ( cube.v6 - cube.v2 ) * t ;
    Dt = cube.v1 + ( cube.v5 - cube.v1 ) * t ;
    break ;

  case  6 :
  case  7 :
  case 12 :
  case 13 :
    switch( _case )
    {
    case  6 : edge = test6 (_config, 2) ; break ;
    case  7 : edge = test7 (_config, 4) ; break ;
    case 12 : edge = test12(_config, 3) ; break ;
    case 13 : edge = tiling13_5_1(_config, _subconfig)[0] ; break ;
    }
    switch( edge )
    {
    case  0 :
      t  = cube.v0 / ( cube.v0 - cube.v1 );
      At = 0;
      Bt = cube.v3 + ( cube.v2 - cube.v3 ) * t;
      Ct = cube.v7 + ( cube.v6 - cube.v7 ) * t;
      Dt = cube.v4 + ( cube.v5 - cube.v4 ) * t;
      break;
    case  1 :
      t  = cube.v1 / ( cube.v1 - cube.v2 );
      At = 0;
      Bt = cube.v0 + ( cube.v3 - cube.v0 ) * t;
      Ct = cube.v4 + ( cube.v7 - cube.v4 ) * t;
      Dt = cube.v5 + ( cube.v6 - cube.v5 ) * t;
      break;
    case  2 :
      t  = cube.v2 / ( cube.v2 - cube.v3 );
      At = 0;
      Bt = cube.v1 + ( cube.v0 - cube.v1 ) * t;
      Ct = cube.v5 + ( cube.v4 - cube.v5 ) * t;
      Dt = cube.v6 + ( cube.v7 - cube.v6 ) * t;
      break;
    case  3 :
      t  = cube.v3 / ( cube.v3 - cube.v0 );
      At = 0;
      Bt = cube.v2 + ( cube.v1 - cube.v2 ) * t;
      Ct = cube.v6 + ( cube.v5 - cube.v6 ) * t;
      Dt = cube.v7 + ( cube.v4 - cube.v7 ) * t;
      break;
    case  4 :
      t  = cube.v4 / ( cube.v4 - cube.v5 );
      At = 0;
      Bt = cube.v7 + ( cube.v6 - cube.v7 ) * t;
      Ct = cube.v3 + ( cube.v2 - cube.v3 ) * t;
      Dt = cube.v0 + ( cube.v1 - cube.v0 ) * t;
      break;
    case  5 :
      t  = cube.v5 / ( cube.v5 - cube.v6 );
      At = 0;
      Bt = cube.v4 + ( cube.v7 - cube.v4 ) * t;
      Ct = cube.v0 + ( cube.v3 - cube.v0 ) * t;
      Dt = cube.v1 + ( cube.v2 - cube.v1 ) * t;
      break;
    case  6 :
      t  = cube.v6 / ( cube.v6 - cube.v7 );
      At = 0;
      Bt = cube.v5 + ( cube.v4 - cube.v5 ) * t;
      Ct = cube.v1 + ( cube.v0 - cube.v1 ) * t;
      Dt = cube.v2 + ( cube.v3 - cube.v2 ) * t;
      break;
    case  7 :
      t  = cube.v7 / ( cube.v7 - cube.v4 );
      At = 0;
      Bt = cube.v6 + ( cube.v5 - cube.v6 ) * t;
      Ct = cube.v2 + ( cube.v1 - cube.v2 ) * t;
      Dt = cube.v3 + ( cube.v0 - cube.v3 ) * t;
      break;
    case  8 :
      t  = cube.v0 / ( cube.v0 - cube.v4 );
      At = 0;
      Bt = cube.v3 + ( cube.v7 - cube.v3 ) * t;
      Ct = cube.v2 + ( cube.v6 - cube.v2 ) * t;
      Dt = cube.v1 + ( cube.v5 - cube.v1 ) * t;
      break;
    case  9 :
      t  = cube.v1 / ( cube.v1 - cube.v5 );
      At = 0;
      Bt = cube.v0 + ( cube.v4 - cube.v0 ) * t;
      Ct = cube.v3 + ( cube.v7 - cube.v3 ) * t;
      Dt = cube.v2 + ( cube.v6 - cube.v2 ) * t;
      break;
    case 10 :
      t  = cube.v2 / ( cube.v2 - cube.v6 );
      At = 0;
      Bt = cube.v1 + ( cube.v5 - cube.v1 ) * t;
      Ct = cube.v0 + ( cube.v4 - cube.v0 ) * t;
      Dt = cube.v3 + ( cube.v7 - cube.v3 ) * t;
      break;
    case 11 :
      t  = cube.v3 / ( cube.v3 - cube.v7 );
      At = 0;
      Bt = cube.v2 + ( cube.v6 - cube.v2 ) * t;
      Ct = cube.v1 + ( cube.v5 - cube.v1 ) * t;
      Dt = cube.v0 + ( cube.v4 - cube.v0 ) * t;
      break;
    default :
      //printf( "Invalid edge %d\n", edge ) ;  print_cube() ;
      break;
    }
    break ;

  default :
    //printf( "Invalid ambiguous case %d\n", _case ) ;  print_cube() ;
    break ;
  }

  if( At >= 0 ) test ++ ;
  if( Bt >= 0 ) test += 2 ;
  if( Ct >= 0 ) test += 4 ;
  if( Dt >= 0 ) test += 8 ;
  switch( test )
  {
  case  0 : return s>0 ;
  case  1 : return s>0 ;
  case  2 : return s>0 ;
  case  3 : return s>0 ;
  case  4 : return s>0 ;
  case  5 : if( At * Ct - Bt * Dt <  FLT_EPSILON ) return s>0 ; break ;
  case  6 : return s>0 ;
  case  7 : return s<0 ;
  case  8 : return s>0 ;
  case  9 : return s>0 ;
  case 10 : if( At * Ct - Bt * Dt >= FLT_EPSILON ) return s>0 ; break ;
  case 11 : return s<0 ;
  case 12 : return s>0 ;
  case 13 : return s<0 ;
  case 14 : return s<0 ;
  case 15 : return s<0 ;
  }

  return s<0 ;
}

KOKKOS_INLINE_FUNCTION
Real add_triangle(const char *trig, char n, Cube cube) {
  // get triangle vertices
  Real xt[3],yt[3],zt[3];
  // initialize the area
  Real area = 0.0;
  Real xv, yv, zv;
  int nv;

  for( int t = 0 ; t < 3*n ; t++ ){
    switch( trig[t] ) {
      case  0:
        xt[ t % 3 ] = 1./(1-cube.v1/cube.v0);
        yt[ t % 3 ] = 0.0;
        zt[ t % 3 ] = 0.0;
        break;
      case  1:
        xt[ t % 3 ] = 1.0;
        yt[ t % 3 ] = 1./(1-cube.v2/cube.v1);
        zt[ t % 3 ] = 0.0;
        break;
      case  2:
        xt[ t % 3 ] = 1./(1-cube.v2/cube.v3);
        yt[ t % 3 ] = 1.0;
        zt[ t % 3 ] = 0.0;
        break;
      case  3:
        xt[ t % 3 ] = 0.0;
        yt[ t % 3 ] = 1./(1-cube.v3/cube.v0);
        zt[ t % 3 ] = 0.0;
        break;
      case  4:
        xt[ t % 3 ] = 1./(1-cube.v5/cube.v4);
        yt[ t % 3 ] = 0.0;
        zt[ t % 3 ] = 1.0;
        break;
      case  5:
        xt[ t % 3 ] = 1.0;
        yt[ t % 3 ] = 1./(1-cube.v6/cube.v5);
        zt[ t % 3 ] = 1.0;
        break ;
      case  6:
        xt[ t % 3 ] = 1./(1-cube.v6/cube.v7);
        yt[ t % 3 ] = 1.0;
        zt[ t % 3 ] = 1.0;
        break;
      case  7:
        xt[ t % 3 ] = 0.0;
        yt[ t % 3 ] = 1./(1-cube.v7/cube.v4);
        zt[ t % 3 ] = 1.0;
        break;
      case  8:
        xt[ t % 3 ] = 0.0;
        yt[ t % 3 ] = 0.0;
        zt[ t % 3 ] = 1./(1-cube.v4/cube.v0);
        break;
      case  9:
        xt[ t % 3 ] = 1.0;
        yt[ t % 3 ] = 0.0;
        zt[ t % 3 ] = 1./(1-cube.v5/cube.v1);
        break;
      case 10:
        xt[ t % 3 ] = 1.0;
        yt[ t % 3 ] = 1.0;
        zt[ t % 3 ] = 1./(1-cube.v6/cube.v2);
        break;
      case 11:
        xt[ t % 3 ] = 0.0;
        yt[ t % 3 ] = 1.0;
        zt[ t % 3 ] = 1./(1-cube.v7/cube.v3);
        break;
      case 12:
        xv=0.0;
        yv=0.0;
        zv=0.0;
        nv = 0;
        for(int iv=0; iv<12; iv++) {
          switch (iv) {
            case  0:
              if (cube.v1*cube.v0 < 0) {
                xv += 1./(1-cube.v1/cube.v0);
                yv += 0.0;
                zv += 0.0;
                nv++;
                break;
              } else {
                break;
              }
            case  1:
              if (cube.v2*cube.v1 < 0) {
                xv += 1.0;
                yv += 1./(1-cube.v2/cube.v1);
                zv += 0.0;
                nv++;
                break;
              } else {
                break;
              }
            case  2:
              if (cube.v2*cube.v3 < 0) {
                xv += 1./(1-cube.v2/cube.v3);
                yv += 1.0;
                zv += 0.0;
                nv++;
                break;
              } else {
                break;
              }
            case  3:
              if (cube.v3*cube.v0 < 0) {
                xv += 0.0;
                yv += 1./(1-cube.v3/cube.v0);
                zv += 0.0;
                nv++;
                break;
              } else {
                break;
              }
            case  4:
              if (cube.v5*cube.v4 < 0) {
                xv += 1./(1-cube.v5/cube.v4);
                yv += 0.0;
                zv += 1.0;
                nv++;
                break;
              } else {
                break;
              }
              break;
            case  5:
              if (cube.v6*cube.v5 < 0) {
                xv += 1.0;
                yv += 1./(1-cube.v6/cube.v5);
                zv += 1.0;
                nv++;
                break;
              } else {
                break;
              }
            case  6:
              if (cube.v6*cube.v7 < 0) {
                xv += 1./(1-cube.v6/cube.v7);
                yv += 1.0;
                zv += 1.0;
                nv++;
                break;
              } else {
                break;
              }
            case  7:
              if (cube.v7*cube.v4 < 0) {
                xv += 0.0;
                yv += 1./(1-cube.v7/cube.v4);
                zv += 1.0;
                nv++;
                break;
              } else {
                break;
              }
            case  8:
              if (cube.v4*cube.v0 < 0) {
                xv += 0.0;
                yv += 0.0;
                zv += 1./(1-cube.v4/cube.v0);
                nv++;
                break;
              } else {
                break;
              }
            case  9:
              if (cube.v5*cube.v1 < 0) {
                xv += 1.0;
                yv += 0.0;
                zv += 1./(1-cube.v5/cube.v1);
                nv++;
                break;
              } else {
                break;
              }
            case 10:
              if (cube.v6*cube.v2 < 0) {
                xv += 1.0;
                yv += 1.0;
                zv += 1./(1-cube.v6/cube.v2);
                nv++;
                break;
              } else {
                break;
              }
            case 11:
              if (cube.v7*cube.v3 < 0) {
                xv += 0.0;
                yv += 1.0;
                zv += 1./(1-cube.v7/cube.v3);
                nv++;
                break;
              } else {
                break;
              }
          }
        }
        // average the positions of edge vertices
        xt[ t % 3 ] = xv/nv;
        yt[ t % 3 ] = yv/nv;
        zt[ t % 3 ] = zv/nv;
        break;
      default :
        break;
    }

    //if( tv[t%3] == -1 ) {
    //  printf("Marching Cubes: invalid triangle %d\n", _ntrigs+1) ;
    //  //print_cube() ;
    //}

    if( t%3 == 2 ) {
      // calculate area of the triangle and add it to the total area
      // first calculate vectors of the sides
      Real v1,v2,v3,w1,w2,w3,u1,u2,u3;
      v1 = xt[0] - xt[1];
      v2 = yt[0] - yt[1];
      v3 = zt[0] - zt[1];
      w1 = xt[0] - xt[2];
      w2 = yt[0] - yt[2];
      w3 = zt[0] - zt[2];
      // calculate the cross product
      u1 = v2*w3 - v3*w2;
      u2 = v3*w1 - v1*w3;
      u3 = v1*w2 - v2*w1;
      // area is half the norm of the cross product
      area += 0.5*sqrt(u1*u1 + u2*u2 + u3*u3);
    }
  }
  return area;
}

KOKKOS_INLINE_FUNCTION
Real process_cube(Cube cube) {
  unsigned char index = 0;
  if (cube.v0 < 0.0) index |= 1;
  if (cube.v1 < 0.0) index |= 2;
  if (cube.v2 < 0.0) index |= 4;
  if (cube.v3 < 0.0) index |= 8;
  if (cube.v4 < 0.0) index |= 16;
  if (cube.v5 < 0.0) index |= 32;
  if (cube.v6 < 0.0) index |= 64;
  if (cube.v7 < 0.0) index |= 128;

  // get case and configuration
  signed char _case = cases(index, 0);
  signed char _config = cases(index, 1);
  // _subconfig is set in certain cases
  signed char _subconfig = 0;

  // "The Big Switch"
  switch( _case ) {
  case  0 :
    return 0.0;

  case  1 :
    return add_triangle( tiling1(_config), 1, cube);

  case  2 :
    return add_triangle( tiling2(_config), 2, cube);

  case  3 :
    if( test_face(test3(_config), cube) )
      return add_triangle(tiling3_2(_config), 4, cube); // 3.2
    else
      return add_triangle(tiling3_1(_config), 2, cube); // 3.1

  case  4 :
    if( test_interior( test4(_config), cube, _case, _config, _subconfig) )
      return add_triangle( tiling4_1(_config), 2, cube); // 4.1.1
    else
      return add_triangle( tiling4_2(_config), 6, cube); // 4.1.2

  case  5 :
    return add_triangle( tiling5(_config), 3, cube);

  case  6 :
    if( test_face(test6(_config, 0), cube) )
      return add_triangle( tiling6_2(_config), 5, cube); // 6.2
    else
    {
      if( test_interior( test6(_config, 1), cube, _case, _config, _subconfig) )
        return add_triangle( tiling6_1_1(_config), 3, cube); // 6.1.1
      else
        return add_triangle( tiling6_1_2(_config), 9, cube); // 6.1.2
    }
    return 0.0;

  case  7 :
    if( test_face( test7(_config, 0), cube ) ) _subconfig +=  1 ;
    if( test_face( test7(_config, 1), cube ) ) _subconfig +=  2 ;
    if( test_face( test7(_config, 2), cube ) ) _subconfig +=  4 ;
    switch( _subconfig )
      {
      case 0 :
        return add_triangle( tiling7_1(_config), 3, cube);
      case 1 :
        return add_triangle( tiling7_2(_config, 0), 5, cube);
      case 2 :
        return add_triangle( tiling7_2(_config, 1), 5, cube);
      case 3 :
        return add_triangle( tiling7_3(_config, 0), 9, cube);
      case 4 :
        return add_triangle( tiling7_2(_config, 2), 5, cube);
      case 5 :
        return add_triangle( tiling7_3(_config, 1), 9, cube);
      case 6 :
        return add_triangle( tiling7_3(_config, 2), 9, cube);
      case 7 :
        if( test_interior( test7(_config, 3), cube, _case, _config, _subconfig) )
          return add_triangle( tiling7_4_2(_config), 9, cube);
        else
          return add_triangle( tiling7_4_1(_config), 5, cube);
      default :
        return 0.0;
      };

  case  8 :
    return add_triangle( tiling8(_config), 2, cube);

  case  9 :
    return add_triangle( tiling9(_config), 4, cube);

  case 10 :
    if( test_face( test10(_config, 0), cube) )
    {
      if( test_face( test10(_config, 1), cube) )
        return add_triangle( tiling10_1_1_(_config), 4, cube); // 10.1.1
      else
        return add_triangle( tiling10_2(_config), 8, cube) ; // 10.2
    }
    else
    {
      if( test_face( test10(_config, 1), cube) )
        return add_triangle( tiling10_2_(_config), 8, cube) ; // 10.2
      else
      {
        if( test_interior( test10(_config, 2), cube, _case, _config, _subconfig) )
          return add_triangle(tiling10_1_1(_config), 4, cube) ; // 10.1.1
        else
          return add_triangle(tiling10_1_2(_config), 8, cube) ; // 10.1.2
      }
    }
    return 0.0;

  case 11 :
    return add_triangle( tiling11(_config), 4, cube) ;

  case 12 :
    if( test_face( test12(_config, 0), cube) )
    {
      if( test_face( test12(_config, 1), cube) )
        return add_triangle( tiling12_1_1_(_config), 4, cube) ; // 12.1.1
      else
        return add_triangle( tiling12_2(_config), 8, cube) ; // 12.2
    }
    else
    {
      if( test_face( test12(_config, 1), cube) )
        return add_triangle( tiling12_2_(_config), 8, cube) ; // 12.2
      else
      {
        if( test_interior( test12(_config, 2), cube, _case, _config, _subconfig) )
          return add_triangle(tiling12_1_1(_config), 4, cube); // 12.1.1
        else
          return add_triangle(tiling12_1_2(_config), 8, cube); // 12.1.2
      }
    }
    return 0.0;

  case 13 :
    if( test_face( test13(_config, 0), cube ) ) _subconfig +=  1 ;
    if( test_face( test13(_config, 1), cube ) ) _subconfig +=  2 ;
    if( test_face( test13(_config, 2), cube ) ) _subconfig +=  4 ;
    if( test_face( test13(_config, 3), cube ) ) _subconfig +=  8 ;
    if( test_face( test13(_config, 4), cube ) ) _subconfig += 16 ;
    if( test_face( test13(_config, 5), cube ) ) _subconfig += 32 ;
    switch( subconfig13(_subconfig) )
    {
      case 0 :/* 13.1 */
        return add_triangle(tiling13_1(_config), 4, cube);

      case 1 :/* 13.2 */
        return add_triangle(tiling13_2(_config, 0), 6, cube);
      case 2 :/* 13.2 */
        return add_triangle(tiling13_2(_config, 1), 6, cube);
      case 3 :/* 13.2 */
        return add_triangle(tiling13_2(_config, 2), 6, cube);
      case 4 :/* 13.2 */
        return add_triangle(tiling13_2(_config, 3), 6, cube);
      case 5 :/* 13.2 */
        return add_triangle(tiling13_2(_config, 4), 6, cube);
      case 6 :/* 13.2 */
        return add_triangle(tiling13_2(_config, 5), 6, cube);

      case 7 :/* 13.3 */
        return add_triangle(tiling13_3(_config, 0), 10, cube);
      case 8 :/* 13.3 */
        return add_triangle(tiling13_3(_config, 1), 10, cube);
      case 9 :/* 13.3 */
        return add_triangle(tiling13_3(_config, 2), 10, cube);
      case 10 :/* 13.3 */
        return add_triangle(tiling13_3(_config, 3), 10, cube);
      case 11 :/* 13.3 */
        return add_triangle(tiling13_3(_config, 4), 10, cube);
      case 12 :/* 13.3 */
        return add_triangle(tiling13_3(_config, 5), 10, cube);
      case 13 :/* 13.3 */
        return add_triangle(tiling13_3(_config, 6), 10, cube);
      case 14 :/* 13.3 */
        return add_triangle(tiling13_3(_config, 7), 10, cube);
      case 15 :/* 13.3 */
        return add_triangle(tiling13_3(_config, 8), 10, cube);
      case 16 :/* 13.3 */
        return add_triangle(tiling13_3(_config, 9), 10, cube);
      case 17 :/* 13.3 */
        return add_triangle(tiling13_3(_config,10), 10, cube);
      case 18 :/* 13.3 */
        return add_triangle(tiling13_3(_config,11), 10, cube);

      case 19 :/* 13.4 */
        return add_triangle(tiling13_4(_config, 0), 12, cube);
      case 20 :/* 13.4 */
        return add_triangle(tiling13_4(_config, 1), 12, cube);
      case 21 :/* 13.4 */
        return add_triangle(tiling13_4(_config, 2), 12, cube);
      case 22 :/* 13.4 */
        return add_triangle(tiling13_4(_config, 3), 12, cube);

      case 23 :/* 13.5 */
        _subconfig = 0 ;
        if(test_interior( test13(_config, 6) , cube, _case, _config, _subconfig))
          return add_triangle(tiling13_5_1(_config, 0), 6, cube);
        else
          return add_triangle(tiling13_5_2(_config, 0), 10, cube);
      case 24 :/* 13.5 */
        _subconfig = 1 ;
        if(test_interior( test13(_config, 6) , cube, _case, _config, _subconfig))
          return add_triangle(tiling13_5_1(_config, 1), 6, cube);
        else
          return add_triangle(tiling13_5_2(_config, 1), 10, cube);
      case 25 :/* 13.5 */
        _subconfig = 2 ;
        if(test_interior( test13(_config, 6) , cube, _case, _config, _subconfig))
          return add_triangle(tiling13_5_1(_config, 2), 6, cube);
        else
          return add_triangle(tiling13_5_2(_config, 2), 10, cube);
      case 26 :/* 13.5 */
        _subconfig = 3 ;
        if(test_interior( test13(_config, 6) , cube, _case, _config, _subconfig))
          return add_triangle(tiling13_5_1(_config, 3), 6, cube);
        else
          return add_triangle(tiling13_5_2(_config, 3), 10, cube);

      case 27 :/* 13.3 */
        return add_triangle(tiling13_3_(_config, 0), 10, cube);
      case 28 :/* 13.3 */
        return add_triangle(tiling13_3_(_config, 1), 10, cube);
      case 29 :/* 13.3 */
        return add_triangle(tiling13_3_(_config, 2), 10, cube);
      case 30 :/* 13.3 */
        return add_triangle(tiling13_3_(_config, 3), 10, cube);
      case 31 :/* 13.3 */
        return add_triangle(tiling13_3_(_config, 4), 10, cube);
      case 32 :/* 13.3 */
        return add_triangle(tiling13_3_(_config, 5), 10, cube);
      case 33 :/* 13.3 */
        return add_triangle(tiling13_3_(_config, 6), 10, cube);
      case 34 :/* 13.3 */
        return add_triangle(tiling13_3_(_config, 7), 10, cube);
      case 35 :/* 13.3 */
        return add_triangle(tiling13_3_(_config, 8), 10, cube);
      case 36 :/* 13.3 */
        return add_triangle(tiling13_3_(_config, 9), 10, cube);
      case 37 :/* 13.3 */
        return add_triangle(tiling13_3_(_config,10), 10, cube);
      case 38 :/* 13.3 */
        return add_triangle(tiling13_3_(_config,11), 10, cube);

      case 39 :/* 13.2 */
        return add_triangle(tiling13_2_(_config, 0), 6, cube);
      case 40 :/* 13.2 */
        return add_triangle(tiling13_2_(_config, 1), 6, cube);
      case 41 :/* 13.2 */
        return add_triangle(tiling13_2_(_config, 2), 6, cube);
      case 42 :/* 13.2 */
        return add_triangle(tiling13_2_(_config, 3), 6, cube);
      case 43 :/* 13.2 */
        return add_triangle(tiling13_2_(_config, 4), 6, cube);
      case 44 :/* 13.2 */
        return add_triangle(tiling13_2_(_config, 5), 6, cube);

      case 45 :/* 13.1 */
        return add_triangle(tiling13_1_(_config), 4, cube);

      default :
        //std::cout << "Marching Cubes: Impossible case 13?" << std::endl;
        return 0.0;
      }
    return 0.0;
  
  case 14 :
    return add_triangle(tiling14(_config), 4, cube);
  }
  return 0.0;
}

#endif // UTILS_MARCHING_CUBES_HPP_
