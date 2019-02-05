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

void Update(vector<vec3>& stars);
void Draw(screen* screen, vector <vec3>& stars);
// void Interpolate(float a, float b, vector<float>& result);
void Interpolate(vec3 a, vec3 b, vector<vec3>& result);

int main( int argc, char* argv[] )
{
  //Initialise stars
  vector<vec3> stars( 1000 );

  //Randomise position of stars
  for (int i = 0; i < 1000; i ++)
  {
    //Gen random numbers
    float rx = (2*(float(rand()) / float(RAND_MAX))) - 1; //between -1,1
    float ry = (2*(float(rand()) / float(RAND_MAX))) - 1; //between -1,1
    float rz = float(rand()) / float(RAND_MAX); //between 0,1

    stars[i].x = rx;
    stars[i].y = ry;
    stars[i].z = rz;

  }


  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  t = SDL_GetTicks();	/*Set start value for timer.*/

  while( NoQuitMessageSDL() )
    {
      Draw(screen, stars);
      Update(stars);
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen, vector<vec3>& stars)
{
  /* Clear buffer */
  //Sets screen to black
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  float f = SCREEN_HEIGHT/2;
  vec3 white(1,1,1);

  //Loop through stars
  for (size_t s = 0; s < stars.size(); ++s)
  {
    int u = (f * stars[s].x / stars[s].z) + (SCREEN_WIDTH/2);
    int v = (f * stars[s].y / stars[s].z) + (SCREEN_HEIGHT/2);

    vec3 colour = 0.2f * white / (stars[s].z * stars[s].z);

    PutPixelSDL(screen, u, v, colour);
  }

}

/*Place updates of parameters here*/
void Update(vector<vec3>& stars)
{
  /* Compute frame time */
  //returns # miliseconds since sdl was initialised
  static int t = SDL_GetTicks();
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;

  float velocity = 0.001;

  /*Good idea to remove this*/
  // std::cout << "Render time: " << dt << " ms." << std::endl;

  /* Update variables*/
  for (int s = 0; s < stars.size(); ++s )
  {
    // Add code for update of stars
    // x and y stay the same

    stars[s].z = stars[s].z - (velocity * dt);


    if ( stars[s].z <= 0 ) { stars[s].z += 1;}
    if ( stars[s].z >  1 ) { stars[s].z -= 1;}
  }

}

void Interpolate(vec3 a, vec3 b, vector<vec3>& result)
{
  //Get size
  int size = result.size();

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
