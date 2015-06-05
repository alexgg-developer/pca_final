/*
This file is part of ftdock, a program for rigid-body protein-protein docking 
Copyright (C) 1997-2000 Gidon Moont

Biomolecular Modelling Laboratory
Imperial Cancer Research Fund
44 Lincoln's Inn Fields
London WC2A 3PX

+44 (0)20 7269 3348
http://www.bmm.icnet.uk/

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#include "structures.h"

void discretise_structure( struct Structure This_Structure , float grid_span , unsigned int grid_size , fftw_real *grid ) {

/************/

  /* Counters */

  int	residue , atom ;

  /* Co-ordinates */

  int	x , y , z ;
  int	steps , x_step , y_step , z_step ;

  float		x_centre , y_centre , z_centre ;
  float		x_centre1 , y_centre1 , z_centre1 ;
  float		x_centre2 , y_centre2 , z_centre2 ;
  float		x_centre3 , y_centre3 , z_centre3 ;

  /* Variables */

  float         distance , one_span, sqrDistance ;

/************/

  one_span = grid_span / (float)grid_size ;

  distance = 1.8 ;
  sqrDistance = 3.24;

/************/

  for( x = 0 ; x < grid_size ; x ++ ) {
    for( y = 0 ; y < grid_size ; y ++ ) {
      for( z = 0 ; z < grid_size ; z ++ ) {

        grid[gaddress(x,y,z,grid_size)] = (fftw_real)0 ;

      }
    }
  }

/************/

  steps = (int)( ( distance / one_span ) + 1.5 ) ;
  unsigned int forStep = 4;
  //allargar la distancia on es fa a on és necessita.


  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {
      for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {
          x = gord( This_Structure.Residue[residue].Atom[atom].coord[1] , grid_span , grid_size ) ;
          y = gord( This_Structure.Residue[residue].Atom[atom].coord[2] , grid_span , grid_size ) ;
          z = gord( This_Structure.Residue[residue].Atom[atom].coord[3] , grid_span , grid_size ) ;
          int max_x_step = max_int( ( x - steps ) , 0 );
          int min_x_step = min_int( ( x + steps ) , ( grid_size - 1 ) );
          int max_y_step = max_int( ( y - steps ) , 0 );
          int min_y_step = min_int( ( y + steps ) , ( grid_size - 1 ) );
          int max_z_step = max_int( ( z - steps ) , 0 );
          int min_z_step = min_int( ( z + steps ) , ( grid_size - 1 ) );
          for( x_step = max_x_step; x_step <= min_x_step; x_step ++ ) {
              x_centre = gcentre( x_step , grid_span , grid_size ) ;
              for( y_step = max_y_step ; y_step <= min_y_step ; y_step ++ ) {
                  y_centre = gcentre( y_step , grid_span , grid_size ) ;
                  for( z_step = max_z_step ; z_step - forStep <= min_z_step ; z_step += forStep ) {
                    z_centre  = gcentre( z_step , grid_span , grid_size ) ;
                    z_centre1  = gcentre( z_step + 1 , grid_span , grid_size ) ;
                    z_centre2  = gcentre( z_step + 2, grid_span , grid_size ) ;
                    z_centre3  = gcentre( z_step + 3, grid_span , grid_size ) ;

                    float a = This_Structure.Residue[residue].Atom[atom].coord[1] - x_centre;
                    float b = This_Structure.Residue[residue].Atom[atom].coord[2] - y_centre;
                    float c = This_Structure.Residue[residue].Atom[atom].coord[3] - z_centre;
                    float sqrDistanceToCenter = a*a + b*b + c*c;
                    //if(sqrDistanceToCenter < sqrDistance)
                    if((sqrDistanceToCenter < sqrDistance)) {
                        grid[gaddress(x_step,y_step,z_step,grid_size)] = (fftw_real)1 ;
                    }


                    float a1 = This_Structure.Residue[residue].Atom[atom].coord[1] - x_centre;
                    float b1 = This_Structure.Residue[residue].Atom[atom].coord[2] - y_centre;
                    float c1 = This_Structure.Residue[residue].Atom[atom].coord[3] - z_centre1;
                    float sqrDistanceToCenter1 = a1*a1 + b1*b1 + c1*c1;
                    //if(sqrDistanceToCenter < sqrDistance)
                    if((sqrDistanceToCenter1 < sqrDistance)) {
                        grid[gaddress(x_step,y_step,z_step + 1,grid_size)] = (fftw_real)1 ;
                    }


                    float a2 = This_Structure.Residue[residue].Atom[atom].coord[1] - x_centre;
                    float b2 = This_Structure.Residue[residue].Atom[atom].coord[2] - y_centre;
                    float c2 = This_Structure.Residue[residue].Atom[atom].coord[3] - z_centre2;
                    float sqrDistanceToCenter2 = a2*a2 + b2*b2 + c2*c2;
                    //if(sqrDistanceToCenter < sqrDistance)
                    if((sqrDistanceToCenter2 < sqrDistance)) {
                        grid[gaddress(x_step,y_step,z_step+2,grid_size)] = (fftw_real)1 ;
                    }


                    float a3 = This_Structure.Residue[residue].Atom[atom].coord[1] - x_centre;
                    float b3 = This_Structure.Residue[residue].Atom[atom].coord[2] - y_centre;
                    float c3 = This_Structure.Residue[residue].Atom[atom].coord[3] - z_centre3;
                    float sqrDistanceToCenter3 = a3*a3 + b3*b3 + c3*c3;
                    //if(sqrDistanceToCenter < sqrDistance)
                    if((sqrDistanceToCenter3 < sqrDistance)) {
                        grid[gaddress(x_step,y_step,z_step+3,grid_size)] = (fftw_real)1 ;
                    }


                  }

                  for( z_step; z_step <= min_z_step ; z_step ++ ) {
                     z_centre = gcentre( z_step , grid_span , grid_size ) ;
                     float a = This_Structure.Residue[residue].Atom[atom].coord[1] - x_centre;
                     float b = This_Structure.Residue[residue].Atom[atom].coord[2] - y_centre;
                     float c = This_Structure.Residue[residue].Atom[atom].coord[3] - z_centre;
                     float sqrDistanceToCenter = a*a + b*b + c*c;
                     //if(sqrDistanceToCenter < sqrDistance)
                     if((sqrDistanceToCenter < sqrDistance)) {
                        grid[gaddress(x_step,y_step,z_step,grid_size)] = (fftw_real)1 ;
                     }
                 }
              }
          }
      }
  }

  /*for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {
    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {

      x = gord( This_Structure.Residue[residue].Atom[atom].coord[1] , grid_span , grid_size ) ;
      y = gord( This_Structure.Residue[residue].Atom[atom].coord[2] , grid_span , grid_size ) ;
      z = gord( This_Structure.Residue[residue].Atom[atom].coord[3] , grid_span , grid_size ) ;
      int max_x_step = max_int( ( x - steps ) , 0 );
      int min_x_step = min_int( ( x + steps ) , ( grid_size - 1 ) );
      int max_y_step = max_int( ( y - steps ) , 0 );
      int min_y_step = min_int( ( y + steps ) , ( grid_size - 1 ) );
      int max_z_step = max_int( ( z - steps ) , 0 );
      int min_z_step = min_int( ( z + steps ) , ( grid_size - 1 ) );
      unsigned int forStep = 4;
      for( x_step =  max_x_step; x_step - forStep <= min_x_step; x_step += forStep ) {

        x_centre  = gcentre( x_step , grid_span , grid_size ) ;
        x_centre1  = gcentre( x_step , grid_span , grid_size ) ;
        x_centre2  = gcentre( x_step , grid_span , grid_size ) ;
        x_centre3  = gcentre( x_step , grid_span , grid_size ) ;

        for( y_step = max_y_step ; y_step - forStep <= min_y_step ; y_step += forStep ) {

          y_centre  = gcentre( y_step , grid_span , grid_size ) ;
          y_centre1  = gcentre( y_step , grid_span , grid_size ) ;
          y_centre2  = gcentre( y_step , grid_span , grid_size ) ;
          y_centre3  = gcentre( y_step , grid_span , grid_size ) ;

          for( z_step = max_z_step ; z_step - forStep <= min_z_step ; z_step += forStep ) {

            z_centre  = gcentre( z_step , grid_span , grid_size ) ;
            z_centre1  = gcentre( z_step , grid_span , grid_size ) ;
            z_centre2  = gcentre( z_step , grid_span , grid_size ) ;
            z_centre3  = gcentre( z_step , grid_span , grid_size ) ;

            float a = This_Structure.Residue[residue].Atom[atom].coord[1] - x_centre;
            float b = This_Structure.Residue[residue].Atom[atom].coord[2] - y_centre;
            float c = This_Structure.Residue[residue].Atom[atom].coord[3] - z_centre;
            float sqrDistanceToCenter = a*a + b*b + c*c;
            //if(sqrDistanceToCenter < sqrDistance)
            if((sqrDistanceToCenter < sqrDistance)) {
                grid[gaddress(x_step,y_step,z_step,grid_size)] = (fftw_real)1 ;
            }


            float a1 = This_Structure.Residue[residue].Atom[atom].coord[1] - x_centre1;
            float b1 = This_Structure.Residue[residue].Atom[atom].coord[2] - y_centre1;
            float c1 = This_Structure.Residue[residue].Atom[atom].coord[3] - z_centre1;
            float sqrDistanceToCenter1 = a*a + b*b + c*c;
            //if(sqrDistanceToCenter < sqrDistance)
            if((sqrDistanceToCenter1 < sqrDistance)) {
                grid[gaddress(x_step,y_step,z_step,grid_size)] = (fftw_real)1 ;
            }


            float a2 = This_Structure.Residue[residue].Atom[atom].coord[1] - x_centre2;
            float b2 = This_Structure.Residue[residue].Atom[atom].coord[2] - y_centre2;
            float c2 = This_Structure.Residue[residue].Atom[atom].coord[3] - z_centre2;
            float sqrDistanceToCenter2 = a*a + b*b + c*c;
            //if(sqrDistanceToCenter < sqrDistance)
            if((sqrDistanceToCenter2 < sqrDistance)) {
                grid[gaddress(x_step,y_step,z_step,grid_size)] = (fftw_real)1 ;
            }


            float a3 = This_Structure.Residue[residue].Atom[atom].coord[1] - x_centre3;
            float b3 = This_Structure.Residue[residue].Atom[atom].coord[2] - y_centre3;
            float c3 = This_Structure.Residue[residue].Atom[atom].coord[3] - z_centre3;
            float sqrDistanceToCenter3 = a*a + b*b + c*c;
            //if(sqrDistanceToCenter < sqrDistance)
            if((sqrDistanceToCenter3 < sqrDistance)) {
                grid[gaddress(x_step,y_step,z_step,grid_size)] = (fftw_real)1 ;
            }


          }

          for( z_step; z_step <= min_z_step ; z_step ++ ) {
             z_centre = gcentre( z_step , grid_span , grid_size ) ;
             float a = This_Structure.Residue[residue].Atom[atom].coord[1] - x_centre;
             float b = This_Structure.Residue[residue].Atom[atom].coord[2] - y_centre;
             float c = This_Structure.Residue[residue].Atom[atom].coord[3] - z_centre;
             float sqrDistanceToCenter = a*a + b*b + c*c;
             //if(sqrDistanceToCenter < sqrDistance)
             if((sqrDistanceToCenter < sqrDistance)) {
             grid[gaddress(x_step,y_step,z_step,grid_size)] = (fftw_real)1 ;
             }
         }
        }
      }

    }
  }*/

/************/

  return ;

}



/************************/



void surface_grid( float grid_span , int grid_size , fftw_real *grid , float surface , float internal_value ) {


/************/

  /* Counters */

  int	x , y , z ;
  int	steps , x_step , y_step , z_step ;

  /* Variables */

  float		one_span ;

  int	at_surface ;

/************/

  one_span = grid_span / (float)grid_size ;

/************/

  /* Surface grid atoms */

  steps = (int)( ( surface / one_span ) + 1.5 ) ;

  for( x = 0 ; x < grid_size ; x ++ ) {
    for( y = 0 ; y < grid_size ; y ++ ) {
      for( z = 0 ; z < grid_size ; z ++ ) {

        if( (int)grid[gaddress(x,y,z,grid_size)] == 1 ) {

          at_surface = 0 ;

          for( x_step = max( x - steps , 0 ) ; x_step <= min( x + steps , grid_size - 1 ) ; x_step ++ ) {
            for( y_step = max( y - steps , 0 ) ; y_step <= min( y + steps , grid_size - 1 ) ; y_step ++ ) {
              for( z_step = max( z - steps , 0 ) ; z_step <= min( z + steps , grid_size - 1 ) ; z_step ++ ) {

                if( (int)grid[gaddress(x_step,y_step,z_step,grid_size)] == 0 ) {

                  if( ( (float)( ( ( x_step - x ) * ( x_step - x ) ) + ( ( y_step - y ) * ( y_step - y ) ) + ( ( z_step - z ) * ( z_step - z ) ) ) * one_span * one_span ) < ( surface * surface ) ) at_surface = 1 ;

                }

              }
            }
          }

          if( at_surface == 0 ) grid[gaddress(x,y,z,grid_size)] = (fftw_real)internal_value ;

        }

      }
    }
  }

/************/

  return ;

}
