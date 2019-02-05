#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <stdint.h>

using namespace std;
using glm::vec3;
using glm::mat3;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false


/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */
int t;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw(screen* screen);
// void Interpolate(float a, float b, vector<float>& result);
void Interpolate(vec3 a, vec3 b, vector<vec3>& result);

int main( int argc, char* argv[] )
{

  // vector<float> result(10);
  // Interpolate(7,13,result);
  // for (int i = 0; i < result.size(); ++i)
  // {
  //   cout << result[i] << " ";
  // }

  // vector<vec3> result( 4 );
  // vec3 a(1,4,9.2);
  // vec3 b(4,1,9.8);
  // Interpolate( a, b, result );
  // for( int i=0; i<result.size(); ++i )
  // {
  //     cout << "( "
  //          << result[i].x << ", "
  //          << result[i].y << ", "
  //          << result[i].z << " ) ";
  // }


  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  t = SDL_GetTicks();	/*Set start value for timer.*/

  while( NoQuitMessageSDL() )
    {
      Draw(screen);
      Update();
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  //Set up corners
  vec3 topLeft(1,0,0); // red
  vec3 topRight(0,0,1); // blue
  vec3 bottomRight(0,1,0); // green
  vec3 bottomLeft(1,1,0); // yellow

  //Set up left-most and right-most column
  vector<vec3> leftSide( SCREEN_HEIGHT );
  vector<vec3> rightSide( SCREEN_HEIGHT );
  Interpolate( topLeft, bottomLeft, leftSide );
  Interpolate( topRight, bottomRight, rightSide );

  vector<vec3> tempRow( SCREEN_WIDTH );

  // vec3 colour(1.0,0,0);
  // for(int i=0; i<1000; i++)
  //   {
  //     uint32_t x = rand() % screen->width;
  //     uint32_t y = rand() % screen->height;
  //     PutPixelSDL(screen, x, y, colour);
  //   }

  for (int row = 0; row < SCREEN_HEIGHT; row++)
  {
    Interpolate(leftSide[row],rightSide[row],tempRow);

    for (int col = 0; col < SCREEN_WIDTH; col++)
    {
      PutPixelSDL(screen, col, row, tempRow[col]);
    }
  }

}

/*Place updates of parameters here*/
void Update()
{
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;
  /*Good idea to remove this*/
  // std::cout << "Render time: " << dt << " ms." << std::endl;
  /* Update variables*/
}

void Interpolate(vec3 a, vec3 b, vector<vec3>& result)
{
  //Get size
  int size = result.size();

  cout << "size: " << size;

  //Case size = 1 (avoid div. by 0)
  // if(size == 1)
  // {
  //   result[0]= (a+b)/2;
  //   return;
  // }

  //Interpolate:
  float incrementX = (b[0]-a[0])/(size-1);
  float incrementY = (b[1]-a[1])/(size-1);
  float incrementZ = (b[2]-a[2])/(size-1);

  float accumulatorX = a[0];
  float accumulatorY = a[1];
  float accumulatorZ = a[2];

  for(int i = 0; i<size; i++)
  {
    result[i].x = accumulatorX;
    accumulatorX += incrementX;

    result[i].y = accumulatorY;
    accumulatorY += incrementY;

    result[i].z = accumulatorZ;
    accumulatorZ += incrementZ;
  }

  return;
}

// void Interpolate(float a, float b, vector<float>& result)
// {
//   //Get size
//   int size = result.size();
//
//   //Case size = 1 (avoid div. by 0)
//   if(size == 1)
//   {
//     result[0]= (a+b)/2;
//     return;
//   }
//
//   //Interpolate:
//   float increment = (b-a)/(size-1);
//   float accumulator = a;
//
//   for(int i = 0; i<size; i++)
//   {
//     result[i] = accumulator;
//     accumulator += increment;
//   }
//   return;
// }
