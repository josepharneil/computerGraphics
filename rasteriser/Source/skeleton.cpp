#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;
using glm::ivec2;
using glm::vec2;

/* * * * * * * * * * * * * * * * * * * * * * *
 *              Defines
 * * * * * * * * * * * * * * * * * * * * * * */
#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE true

/* * * * * * * * * * * * * * * * * * * * * * *
 *              Global Variables
 * * * * * * * * * * * * * * * * * * * * * * */
bool quit;

/* * * * * * * * * * * * * * * * * * * * * * *
 *              Structures
 * * * * * * * * * * * * * * * * * * * * * * */

// struct Triangle
// {
//   vec4 v0;
//   vec4 v1;
//   vec4 v2;
//   vec4 normal;
//   vec3 color;
// };

struct Pixel
{
    int x;
    int y;
    float zinv;
};


/* * * * * * * * * * * * * * * * * * * * * * *
 *              FUNCTION DEFS
 * * * * * * * * * * * * * * * * * * * * * * */
void Update(vec4& cameraPos, float& yaw);
void Draw(screen* screen, vector<Triangle>& triangles, vec4& cameraPos, float& focalLength);
// void VertexShader( const vec4& vertex, ivec2& p, vec4& cameraPos, float& focalLength );
// void DrawLineSDL( SDL_Surface* surface, ivec2 a, ivec2 b, vec3 colour );
// void DrawPolygonEdges( screen* screen, const vector<vec4>& vertices, vec4& cameraPos, float& focalLength );
// void DrawPolygon(screen* screen, const vector<vec4>& vertices, vec4& cameraPos, float& focalLength, vec3 colour );
// void ComputePolygonRows(const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels, vector<ivec2>& rightPixels, screen* screen );
// void FindLine( ivec2 a, ivec2 b, vector<ivec2>& lineToDraw);
// void DrawPolygonRows( screen* screen, const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels, vec3 colour );
// void TestComputePolygonRows(screen* screen);
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);

//pixel versions
void InterpolatePixel(Pixel a, Pixel b, vector<Pixel>& result);
void DrawPolygonRows( screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels, vec3 colour, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH]);
void FindLine( Pixel a, Pixel b, vector<Pixel>& lineToDraw);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels, screen* screen );
void DrawPolygon(screen* screen, const vector<vec4>& vertices, vec4& cameraPos, float& focalLength, vec3 colour, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH] );
void VertexShader( const vec4& vertex, Pixel& p, vec4& cameraPos, float& focalLength );

/* * * * * * * * * * * * * * * * * * * * * * *
 *                  Main
 * * * * * * * * * * * * * * * * * * * * * * */
int main( int argc, char* argv[] )
{
  #pragma region Main
  //Initially, do not quit
  quit = false;



  // return  0;


  //TEMPORARY PLACE HOLDER
  float yaw = 0;
  
  // Initialise screen
  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );


  // TestComputePolygonRows(screen);

  // return 0;

  //Instantiate vector of triangles
  vector<Triangle> triangles;
  LoadTestModel( triangles );

  //Camera position (fixed)
  vec4 cameraPos( 0, 0, -3.001,1 );

  //Focal length
  float focalLength = 160;

  while( !quit ) //NoQuitMessageSDL() )
    {
      Update(cameraPos, yaw);
      Draw(screen, triangles, cameraPos, focalLength);
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
  #pragma endregion main
}

/* * * * * * * * * * * * * * * * * * * * * * *
 *                  Draw
 * * * * * * * * * * * * * * * * * * * * * * */
/*Place your drawing here*/
void Draw(screen* screen, vector<Triangle>& triangles, vec4& cameraPos, float& focalLength)
{
  #pragma region Draw
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  //Depth buffer
  float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

  for( int y=0; y<SCREEN_HEIGHT; ++y )
  {
    for( int x=0; x<SCREEN_WIDTH; ++x )
    {
      depthBuffer[y][x] = 0;
    }
  }
    
  
  //For each triangle
  for( uint32_t i=0; i<triangles.size(); ++i )
  {
    //Construct vector of triangle vertices
    vector<vec4> vertices(3);
    vertices[0] = triangles[i].v0;
    vertices[1] = triangles[i].v1;
    vertices[2] = triangles[i].v2;

    // //For each vertex in triangle
    // for(int v=0; v<3; ++v)
    // {
    //   //Projected position on screen
    //   ivec2 projPos;//Integer vector of size 2

    //   //Shade each vertex white
    //   VertexShader( vertices[v], projPos, cameraPos, focalLength );
    //   vec3 color(1,1,1);

    //   //Put pixel in screen
    //   // PutPixelSDL( screen, projPos.x, projPos.y, color );

    //   //Draw all polygon edges
    //   DrawPolygonEdges( screen, vertices, cameraPos, focalLength );
    // }


    //Draw polygon for each triangle
    DrawPolygon( screen, vertices, cameraPos, focalLength, triangles[i].color, depthBuffer );

  }

  #pragma endregion Draw
}

/* * * * * * * * * * * * * * * * * * * * * * *
 *                  Update
 * * * * * * * * * * * * * * * * * * * * * * */
void Update(vec4& cameraPos, float& yaw)
{
  #pragma region Update

  float movementSpeed = 0.1f;

  SDL_Event e;

  while(SDL_PollEvent(&e))
  {

    if( e.type != SDL_KEYDOWN) {continue;}


    //Move camera with arrow keys
    if( e.key.keysym.scancode == SDL_SCANCODE_UP )
    {
      cameraPos.z += movementSpeed;

    }
    if( e.key.keysym.scancode == SDL_SCANCODE_DOWN )
    {
      //Move camera backward
      cameraPos.z -= movementSpeed;

    }
    if( e.key.keysym.scancode == SDL_SCANCODE_LEFT )
    {
      //Move camera left
      cameraPos.x -= movementSpeed;
    }
    if( e.key.keysym.scancode == SDL_SCANCODE_RIGHT )
    {
      //Move camera right
      cameraPos.x += movementSpeed;
    }
    if( e.key.keysym.scancode == SDL_SCANCODE_J )
    {
      //Move camera left
      cameraPos.y -= movementSpeed;
    }
    if( e.key.keysym.scancode == SDL_SCANCODE_M )
    {
      //Move camera right
      cameraPos.y += movementSpeed;
    }

    //Move light with WASD
    // if ( e.key.keysym.scancode == SDL_SCANCODE_W )
    // {
    //   lightPos.z += movementSpeed;
    // }
    // if ( e.key.keysym.scancode == SDL_SCANCODE_S )
    // {
    //   lightPos.z -= movementSpeed;
    // }
    // if ( e.key.keysym.scancode == SDL_SCANCODE_D )
    // {
    //   lightPos.x += movementSpeed;
    // }
    // if ( e.key.keysym.scancode == SDL_SCANCODE_A )
    // {
    //   lightPos.x -= movementSpeed;
    // }

    //Quit trigger
    if( e.type == SDL_QUIT )
    {
      quit = true;
    }
    if( e.type == SDL_KEYDOWN )
    {
      if( e.key.keysym.sym == SDLK_ESCAPE)
      {
        quit = true;
      }
    }
  }

  #pragma endregion Update
}


/* * * * * * * * * * * * * * * * * * * * * * *
 *          Other Functions
 * * * * * * * * * * * * * * * * * * * * * * */
// //Projects scene point (triangle vertex) onto image plane
// void VertexShader( const vec4& vertex, ivec2& p, vec4& cameraPos, float& focalLength )
// {
//   #pragma region VertexShader
//   //Translates vertex so that camera is at origin, relative to the vertex
//   vec4 vertexTransformed = vertex - cameraPos;

//   //vertex[0] is X
//   //vertex[1] is Y
//   //vertex[2] is Z

//   int x = (focalLength * (vertexTransformed[0] / vertexTransformed[2])) + (SCREEN_WIDTH /2);
//   int y = (focalLength * (vertexTransformed[1] / vertexTransformed[2])) + (SCREEN_HEIGHT/2);

//   p.x = x;
//   p.y = y;

//   #pragma endregion VertexShader
// }

//Projects scene point (triangle vertex) onto image plane
void VertexShader( const vec4& vertex, Pixel& p, vec4& cameraPos, float& focalLength )
{
  #pragma region VertexShader
  //Translates vertex so that camera is at origin, relative to the vertex
  vec4 vertexTransformed = vertex - cameraPos;

  //vertex[0] is X
  //vertex[1] is Y
  //vertex[2] is Z

  int x = (focalLength * (vertexTransformed[0] / vertexTransformed[2])) + (SCREEN_WIDTH /2);
  int y = (focalLength * (vertexTransformed[1] / vertexTransformed[2])) + (SCREEN_HEIGHT/2);

  float zinv = 1/vertex.z;

  p.x = x;
  p.y = y;
  p.zinv = zinv;

  #pragma endregion VertexShader
}

//Could be sped up with Bresenham's line algorithm
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result)
{
  #pragma region Interpolate

  int N = result.size();
  vec2 step = vec2(b-a) / float(max(N-1,1));
  vec2 current( a );
  for( int i=0; i<N; ++i )
  {
      result[i].x = round(current.x);
      result[i].y = round(current.y);
      current += step;
  }

  #pragma endregion Interpolate
}

void InterpolatePixel(Pixel a, Pixel b, vector<Pixel>& result)
{
  #pragma region InterpolateP

  int N = result.size();
  Pixel step;
  step.x = (b.x-a.x) / float(max(N-1,1));
  step.y = (b.zinv-a.zinv) / float(max(N-1,1));
  step.zinv = (b.zinv-a.zinv) / float(max(N-1,1));


  Pixel current( a );
  for( int i=0; i<N; ++i )
  {
      result[i].x = round(current.x);
      result[i].y = round(current.y);
      current.x += step.x;
      current.y += step.y;
      current.zinv += step.zinv;
  }

  #pragma endregion InterpolateP
}

#pragma region wireframe
// //Draws a line on screen between a and b of some colour
// void DrawLineSDL( screen* screen, ivec2 a, ivec2 b, vec3 colour )
// {
//   #pragma region DrawLineSDL
//   // Find difference between two vertices.
//   ivec2 delta = glm::abs( a - b );

//   //Find number of pixels needed.
//   int numberOfPixels = glm::max( delta.x, delta.y ) + 1;

//   //Instantiate vector of pixels.
//   vector<ivec2> lineToDraw( numberOfPixels );
  
//   //Interpolate between two vertices
//   Interpolate( a, b, lineToDraw );

//   //Put each pixel in the screen
//   for (int i = 0; i < lineToDraw.size(); i++)
//   {
//     PutPixelSDL(screen, lineToDraw[i].x, lineToDraw[i].y, colour);
//   }

//   #pragma endregion Interpolate
// }

// //Takes vertices of a triangle (or any polygon) and draws the edges
// void DrawPolygonEdges( screen* screen, const vector<vec4>& vertices, vec4& cameraPos, float& focalLength )
// {
//   #pragma region DrawPolygonEdges
//   int V = vertices.size();
  
//   // Initialise vector of length 3 (for triangles)
//   vector<ivec2> projectedVertices( V );

//   //For each vertex
//   for( int i=0; i<V; ++i )
//   {
//     //Calculate the perspective projection
//     VertexShader( vertices[i], projectedVertices[i], cameraPos, focalLength );
//   }
  
//   // Loop over all projected vertices and draw the edge from it to the next vertex:
//   for( int i=0; i<V; ++i )
//   {
//     int j = (i+1)%V; // The next vertex
//     vec3 colour( 1, 1, 1 );//white
//     DrawLineSDL( screen, projectedVertices[i], projectedVertices[j], colour );
//   }

//   #pragma endregion DrawPolygonEdges
// }
#pragma endregion wireframe

//This draws the polygon from some vertices
void DrawPolygon(screen* screen, const vector<vec4>& vertices, vec4& cameraPos, float& focalLength, vec3 colour, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH] )
{
  #pragma region DrawPolygon
  //Find number of vertices of polygon (3 for triangle)
  int V = vertices.size();

  //Initialise 3 long vector of projected vertices (pixel locations)
  vector<Pixel> vertexPixels( V );

  //For each vertex
  for( int i=0; i<V; ++i )
  {
    //Compute projection
    VertexShader( vertices[i], vertexPixels[i], cameraPos, focalLength);

  }

  //Initialise vectors to store left-most and right-most positions of each row of the projected triangle
  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;

  //Calculates leftPixels and rightPixels
  ComputePolygonRows( vertexPixels, leftPixels, rightPixels, screen );

  //Draws the rows
  DrawPolygonRows( screen, leftPixels, rightPixels, colour, depthBuffer );

  #pragma endregion DrawPolygon
}

// void ComputePolygonRows(
//     const vector<ivec2>& vertexPixels,
//     vector<ivec2>& leftPixels,
//     vector<ivec2>& rightPixels, screen* screen )
// {
//   #pragma region ComputePolygonRows
//   // 1. Find max and min y-value of the polygon
//   //    and compute the number of rows it occupies.
//   int yMin;
//   int yMax;

//   //Initialise two bounds
//   yMin = +numeric_limits<int>::max();
//   yMax = -numeric_limits<int>::max();

//   //For each vertex pixel
//   for (int i = 0; i < vertexPixels.size(); i++)
//   {
//     //Update min and max bounds
//     if (vertexPixels[i].y < yMin) { yMin = vertexPixels[i].y;}
//     if (vertexPixels[i].y > yMax) { yMax = vertexPixels[i].y; }
//   }

//   //Number of rows of pixels in the triangle
//   int numberOfRows = (yMax - yMin + 1);


//   // 2. Resize leftPixels and rightPixels
//   //    so that they have an element for each row.
//   leftPixels.resize(numberOfRows);
//   rightPixels.resize(numberOfRows);

//   // 3. Initialize the x-coordinates in leftPixels
//   //    to some really large value and the x-coordinates
//   //    in rightPixels to some really small value.
//   for( int i = 0; i < numberOfRows; ++i )
//   {
//     leftPixels[i].x  = +numeric_limits<int>::max();
//     rightPixels[i].x = -numeric_limits<int>::max();

//     leftPixels[i].y  = yMin + i;
//     rightPixels[i].y = yMin + i;
//   }

//   // 4. Loop through all edges of the polygon and use
//   //    linear interpolation to find the x-coordinate for
//   //    each row it occupies. Update the corresponding
//   //    values in rightPixels and leftPixels.

//   vector<vector<ivec2> > polygonEdges(3);

//   FindLine(vertexPixels[0], vertexPixels[1], polygonEdges[0]);
//   FindLine(vertexPixels[1], vertexPixels[2], polygonEdges[1]);
//   FindLine(vertexPixels[2], vertexPixels[0], polygonEdges[2]);

//   //For each edge
//   for ( int i = 0; i < 3; i++)
//   {
//     vector<ivec2>& edge = polygonEdges[i];
//     //For each pixel on edge
//     for (int j = 0; j < edge.size(); j++)
//     {
//       //wireframe test
//       // PutPixelSDL(screen, edge[j].x, edge[j].y, vec3(1,1,1));

//       ivec2 pixel = edge[j];

//       int rowIndex = (pixel.y - yMin);

//       if(pixel.x < leftPixels[rowIndex].x) { leftPixels[rowIndex].x = pixel.x; }
//       if(pixel.x > rightPixels[rowIndex].x) { rightPixels[rowIndex].x = pixel.x; }
//     }
//   }

//   #pragma endregion ComputePolygonRows
// }

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels, screen* screen )
{
  // 1. Find max and min y-value of the polygon
  //    and compute the number of rows it occupies.
  int yMin;
  int yMax;

  //Initialise two bounds
  yMin = +numeric_limits<int>::max();
  yMax = -numeric_limits<int>::max();

  //For each vertex pixel
  for (int i = 0; i < vertexPixels.size(); i++)
  {
    //Update min and max bounds
    if (vertexPixels[i].y < yMin) { yMin = vertexPixels[i].y;}
    if (vertexPixels[i].y > yMax) { yMax = vertexPixels[i].y; }
  }

  //Number of rows of pixels in the triangle
  int numberOfRows = (yMax - yMin + 1);


  // 2. Resize leftPixels and rightPixels
  //    so that they have an element for each row.
  leftPixels.resize(numberOfRows);
  rightPixels.resize(numberOfRows);

  // 3. Initialize the x-coordinates in leftPixels
  //    to some really large value and the x-coordinates
  //    in rightPixels to some really small value.
  for( int i = 0; i < numberOfRows; ++i )
  {
    leftPixels[i].x  = +numeric_limits<int>::max();
    rightPixels[i].x = -numeric_limits<int>::max();

    leftPixels[i].y  = yMin + i;
    rightPixels[i].y = yMin + i;
  }

  // 4. Loop through all edges of the polygon and use
  //    linear interpolation to find the x-coordinate for
  //    each row it occupies. Update the corresponding
  //    values in rightPixels and leftPixels.

  vector<vector<Pixel> > polygonEdges(3);

  FindLine(vertexPixels[0], vertexPixels[1], polygonEdges[0]);
  FindLine(vertexPixels[1], vertexPixels[2], polygonEdges[1]);
  FindLine(vertexPixels[2], vertexPixels[0], polygonEdges[2]);

  //For each edge
  for ( int i = 0; i < 3; i++)
  {
    vector<Pixel>& edge = polygonEdges[i];
    //For each pixel on edge
    for (int j = 0; j < edge.size(); j++)
    {
      //wireframe test
      // PutPixelSDL(screen, edge[j].x, edge[j].y, vec3(1,1,1));

      Pixel pixel = edge[j];

      int rowIndex = (pixel.y - yMin);

      if(pixel.x < leftPixels[rowIndex].x) { leftPixels[rowIndex].x = pixel.x; }
      if(pixel.x > rightPixels[rowIndex].x) { rightPixels[rowIndex].x = pixel.x; }
    }
  }
}

// void TestComputePolygonRows(screen* screen)
// {
//   vector<ivec2> vertexPixels(3);
//   vertexPixels[0] = ivec2(10, 5);
//   vertexPixels[1] = ivec2( 5,10);
//   vertexPixels[2] = ivec2(15,15);
//   vector<ivec2> leftPixels;
//   vector<ivec2> rightPixels;
//   ComputePolygonRows( vertexPixels, leftPixels, rightPixels, screen );
//   for( int row=0; row<leftPixels.size(); ++row )
//   {
//     cout << "Start: ("
//     << leftPixels[row].x << ","
//     << leftPixels[row].y << "). "
//     << "End: ("
//     << rightPixels[row].x << ","
//     << rightPixels[row].y << "). " << endl;
//   }
// }

// //Finds a line on screen between a and b of some colour
// void FindLine( ivec2 a, ivec2 b, vector<ivec2>& lineToDraw)
// {
//   #pragma region DrawLineSDL
//   // Find difference between two vertices.
//   ivec2 delta = glm::abs( a - b );

//   //Find number of pixels needed.
//   int numberOfPixels = glm::max( delta.x, delta.y ) + 1;

//   //Resize lineToDraw for interpolations
//   lineToDraw.resize(numberOfPixels);
  
//   //Interpolate between two vertices
//   Interpolate( a, b, lineToDraw );

//   #pragma endregion Interpolate
// }

//Finds a line on screen between a and b of some colour
void FindLine( Pixel a, Pixel b, vector<Pixel>& lineToDraw)
{
  #pragma region DrawLineSDL
  // Find difference between two vertices.
  Pixel delta;
  delta.x = glm::abs( a.x - b.x );
  delta.y = glm::abs( a.y - b.y );
  delta.zinv = glm::abs( a.zinv - b.zinv );
   

  //Find number of pixels needed.
  int numberOfPixels = glm::max( delta.x, delta.y ) + 1;

  //Resize lineToDraw for interpolations
  lineToDraw.resize(numberOfPixels);
  
  //Interpolate between two vertices
  InterpolatePixel( a, b, lineToDraw );

  #pragma endregion Interpolate
}

// void DrawPolygonRows( screen* screen, 
//                       const vector<ivec2>& leftPixels, 
//                       const vector<ivec2>& rightPixels, 
//                       vec3 colour)
// {
//   #pragma region DrawPolygonRows

//   //For each entry in left/right pixels
//   for (int i = 0; i < leftPixels.size(); i++)
//   {
//     //For each pixel between the leftmost edgepoint and rightmost edgepoint
//     for (int x = leftPixels[i].x; x < rightPixels[i].x; x++)
//     {
//       //Plot
//       PutPixelSDL(screen, x, leftPixels[i].y,colour);
//     }

//   }

//   #pragma endregion DrawPolygonRows
// }

void DrawPolygonRows( screen* screen, 
                      const vector<Pixel>& leftPixels, 
                      const vector<Pixel>& rightPixels, 
                      vec3 colour, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH])
{
  #pragma region DrawPolygonRows

  //For each entry in left/right pixels
  for (int i = 0; i < leftPixels.size(); i++)
  {
    //For each pixel between the leftmost edgepoint and rightmost edgepoint
    for (int x = leftPixels[i].x; x < rightPixels[i].x; x++)
    {
      //Plot
      PutPixelSDL(screen, x, leftPixels[i].y,colour);
    }

  }

  #pragma endregion DrawPolygonRows
}