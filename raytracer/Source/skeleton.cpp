#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false

// struct Triangle
// {
//   vec4 v0;
//   vec4 v1;
//   vec4 v2;
//   vec4 normal;
//   vec3 color;
// };
/* GLOBAL VARIABLES                                                      */
int t;

struct Intersection
{
    vec4 position;
    float distance;
    int triangleIndex;
};

/* ----------------------------------------------------------------------*/
/* FUNCTIONS                                                             */

void Update();
void Draw(screen* screen, vector<Triangle>& triangles);
bool ClosestIntersection(
  vec4 start,
  vec4 dir,
  const vector<Triangle>& triangles,
  Intersection& closestIntersection);

int main( int argc, char* argv[] )
{

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  t = SDL_GetTicks();	/*Set start value for timer.*/

  //Instantiate vector of triangles
  vector<Triangle> triangles;
  LoadTestModel( triangles );

  //Update and draw
  while( NoQuitMessageSDL() )
    {
      Draw(screen, triangles);
      Update();
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen, vector<Triangle>& triangles)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  float focalLength = 1.0f;
  vec4 cameraPos(0.0f,0.0f,0.0f,1.0f);

  Intersection closestIntersection;

  //For each pixel
  for (int row = 0; row < SCREEN_HEIGHT; row++)
  {
    for (int col = 0; col < SCREEN_WIDTH; col++)
    {
      //Compute ray direction
      vec4 direction = vec4(
        col-SCREEN_WIDTH/2,
        row-SCREEN_HEIGHT/2,
        focalLength,
        1.0f);

      //Compute ClosestIntersection
      bool intersect = ClosestIntersection(
        cameraPos,
        direction,
        triangles,
        closestIntersection);

      //If an intersection occurs
      if (intersect)
      {
        //Get triangle from triangles
        Triangle intersectedTriangle = triangles[closestIntersection.triangleIndex];

        //Get colour of triangle
        vec3 colour = intersectedTriangle.color;

        //set to colour of that triangle
        PutPixelSDL(screen, row, col, colour);
      }
      //Else set to black
      // else
      // {
      //   //set to black
      // }



    }
  }


}

/*Place updates of parameters here*/
void Update()
{
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;
  /*Good idea to remove this*/
  // std::cout << "Render time: " << dt << " ms." << std::endl;
  /* Update variables*/
}


bool ClosestIntersection(
  vec4 start,
  vec4 dir,
  const vector<Triangle>& triangles,
  Intersection& closestIntersection)
{
  //Get closest intersection
  float minDist = numeric_limits<float>::max();
  bool result = false;

  //For each triangle
  for (int i = 0; i < triangles.size(); i++)
  {
    Triangle triangle = triangles[i];

    //3D version of dir (ignore last homogenous component)
    vec3 dir3 = vec3(dir.x, dir.y, dir.z);

    //Define vertices
    vec4 v0 = triangle.v0;
    vec4 v1 = triangle.v1;
    vec4 v2 = triangle.v2;

    //Define axes
    vec3 e1 = vec3(v1.x-v0.x, v1.y-v0.y, v1.z-v0.z);
    vec3 e2 = vec3(v2.x-v0.x, v2.y-v0.y, v2.z-v0.z);

    //b = s - v0
    vec3 b = vec3(start.x - v0.x, start.y - v0.y, start.z - v0.z);

    //A
    mat3 A = mat3( -dir3, e1, e2 );

    //x = (t u v)^T
    //COME BACK HERE todo cramar's rule instead of inbuilt inverse
    vec3 x = glm::inverse( A ) * b;

    float t = x.x;//if neg, can reject immediately
    float u = x.y;
    float v = x.z;

    //Check if is intersection is within triangle
    if ( (u > 0.0f) && (v > 0.0f) && (u + v < 1.0f) && (t > 0.0f) )
    {
      result = true;
      //Compute "intersection" structure
      if (t < minDist)
      {
        //We set w = 1 -- not sure of this
        closestIntersection.position = vec4(t,u,v,1);
        closestIntersection.distance = t;
        closestIntersection.triangleIndex = i;
      }
    }
  }
  return result;
}
