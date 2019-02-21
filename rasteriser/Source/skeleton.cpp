#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <math.h>

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
#define FULLSCREEN_MODE false
#define PI 3.14159265

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
  // vec3 illumination;
  vec4 pos3d;
};

struct Vertex
{
  vec4 position;
  // vec4 normal;
  // vec3 reflectance;
};

/* * * * * * * * * * * * * * * * * * * * * * *
 *              FUNCTION DEFS
 * * * * * * * * * * * * * * * * * * * * * * */
void Update(vec4& cameraPos, int& yaw, mat4& cameraMatrix, vec4& lightPos);
void Draw(screen* screen, vector<Triangle>& triangles, vec4& cameraPos, float& focalLength, vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea);
void TestComputePolygonRows(screen* screen);
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);

//pixel versions
void InterpolatePixel(Pixel a, Pixel b, vector<Pixel>& result);
void DrawPolygonRows( screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels, vec3 colour, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH], vec4& currentNormal, vec3& currentReflectance, vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea );
void FindLine( Pixel a, Pixel b, vector<Pixel>& lineToDraw);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels, screen* screen );
void DrawPolygon(screen* screen, const vector<Vertex>& vertices, vec4& cameraPos, float& focalLength, vec3 colour, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH], vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea, vec4& currentNormal, vec3& currentReflectance );
void VertexShader( const Vertex& vertex, Pixel& p, vec4& cameraPos, float& focalLength );
void PixelShader(const Pixel& p, screen* s, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH], vec3 color, vec4& currentNormal, vec3& currentReflectance, vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea );

/* * * * * * * * * * * * * * * * * * * * * * *
 *                  Main
 * * * * * * * * * * * * * * * * * * * * * * */
int main( int argc, char* argv[] )
{
  //Initially, do not quit
  quit = false;
  
  // Initialise screen
  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

  // TestComputePolygonRows(screen);

  //Instantiate vector of triangles
  vector<Triangle> triangles;
  LoadTestModel( triangles );

  vector<Triangle> originalTriangles;
  LoadTestModel( originalTriangles );

  //Camera position (fixed)
  vec4 cameraPos( 0, 0, -3.001,1 );
  mat4 cameraMatrix;
  int yaw = 0;

  //Light
  vec4 lightPos(0,-0.5,-0.7,1);
  vec4 originalLightPos( 0.0f, -0.5f, -0.7f, 1.0f );
  vec3 lightPower = 30.0f*vec3( 1, 1, 1 );
  vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );


  //Focal length
  float focalLength = SCREEN_WIDTH/2;

  while( !quit ) //NoQuitMessageSDL() )
    {
      Update(cameraPos, yaw, cameraMatrix, originalLightPos);

      mat4 invCameraMatrix = glm::inverse(cameraMatrix);
      
      for (int t = 0; t < triangles.size(); t++)
      {
        triangles[t].v0 = invCameraMatrix * originalTriangles[t].v0;
        triangles[t].v1 = invCameraMatrix * originalTriangles[t].v1;
        triangles[t].v2 = invCameraMatrix * originalTriangles[t].v2;
        triangles[t].ComputeNormal();
      }

      lightPos = invCameraMatrix * originalLightPos;

      Draw(screen, triangles, cameraPos, focalLength, lightPos, lightPower, indirectLightPowerPerArea);
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

/* * * * * * * * * * * * * * * * * * * * * * *
 *                  Draw
 * * * * * * * * * * * * * * * * * * * * * * */
/*Place your drawing here*/
void Draw(screen* screen, vector<Triangle>& triangles, vec4& cameraPos, float& focalLength, vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea)
{
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
    vector<Vertex> vertices(3);
    vertices[0].position = triangles[i].v0;
    vertices[1].position = triangles[i].v1;
    vertices[2].position = triangles[i].v2;

    // vertices[0].normal = triangles[i].normal;
    // vertices[1].normal = triangles[i].normal;
    // vertices[2].normal = triangles[i].normal;

    // vertices[0].reflectance = triangles[i].color;
    // vertices[1].reflectance = triangles[i].color;
    // vertices[2].reflectance = triangles[i].color;

    vec4 currentNormal = triangles[i].normal;
    vec3 currentReflectance = triangles[i].color;

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
    DrawPolygon( screen, vertices, cameraPos, focalLength, triangles[i].color, depthBuffer, 
    lightPos, lightPower, indirectLightPowerPerArea, currentNormal, currentReflectance );

  }
}

/* * * * * * * * * * * * * * * * * * * * * * *
 *                  Update
 * * * * * * * * * * * * * * * * * * * * * * */
void Update(vec4& cameraPos, int& yaw, mat4& cameraMatrix, vec4& lightPos)
{

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
    if ( e.key.keysym.scancode == SDL_SCANCODE_W )
    {
      lightPos.z += movementSpeed;
    }
    if ( e.key.keysym.scancode == SDL_SCANCODE_S )
    {
      lightPos.z -= movementSpeed;
    }
    if ( e.key.keysym.scancode == SDL_SCANCODE_D )
    {
      lightPos.x += movementSpeed;
    }
    if ( e.key.keysym.scancode == SDL_SCANCODE_A )
    {
      lightPos.x -= movementSpeed;
    }

    //ROTATE AROUND Y-AXIS
    if ( e.key.keysym.scancode == SDL_SCANCODE_PERIOD )
    {
      yaw -= 5;
      yaw = yaw % 360;
    }
    if ( e.key.keysym.scancode == SDL_SCANCODE_COMMA )
    {
      yaw += 5;
      yaw = yaw % 360;
    }

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

  cameraMatrix[0][0] = cos( yaw * PI / 180 );
  cameraMatrix[0][1] = 0;
  cameraMatrix[0][2] = sin( yaw * PI / 180 );
  cameraMatrix[0][3] = 0;

  cameraMatrix[1][0] = 0;
  cameraMatrix[1][1] = 1;
  cameraMatrix[1][2] = 0;
  cameraMatrix[1][3] = 0;

  cameraMatrix[2][0] = -sin( yaw * PI / 180 );
  cameraMatrix[2][1] = 0;
  cameraMatrix[2][2] = cos( yaw * PI / 180 );
  cameraMatrix[2][3] = 0;

  cameraMatrix[3][0] = cameraPos.x;
  cameraMatrix[3][1] = cameraPos.y;
  cameraMatrix[3][2] = cameraPos.z;
  cameraMatrix[3][3] = 1;

}

/* * * * * * * * * * * * * * * * * * * * * * *
 *          Other Functions
 * * * * * * * * * * * * * * * * * * * * * * */

//Projects scene point (triangle vertex) onto image plane
void VertexShader( const Vertex& vertex, Pixel& p, vec4& cameraPos, float& focalLength)
{
  //Translates vertex so that camera is at origin, relative to the vertex
  vec4 vertexTransformed = vertex.position - cameraPos;

  //vertex[0] is X
  //vertex[1] is Y
  //vertex[2] is Z

  int x = (focalLength * (vertexTransformed[0] / vertexTransformed[2])) + (SCREEN_WIDTH /2);
  int y = (focalLength * (vertexTransformed[1] / vertexTransformed[2])) + (SCREEN_HEIGHT/2);

  //float zinv = 1.0f/glm::distance(cameraPos, vertex);
  float zinv;
  if(vertexTransformed.z == 0)
  {
    zinv = 0;
  }
  else
  {
    zinv = 1.0f/(vertexTransformed.z);
  }

  // vec4 r = lightPos - vertex.position;
  // vec4 n = vertex.normal;

  // vec4 rNorm = normalize(r);
  // vec4 nNorm = normalize(n);
  
  // float A = 4 * M_PI * ( pow(glm::length(r),2) );
  // vec3 B = vec3(lightPower.x/A,lightPower.y/A,lightPower.z/A);

  // vec3 D = B * max(glm::dot (rNorm,nNorm), 0.0f);

  // vec3 illumination = vertex.reflectance * (D + indirectLightPowerPerArea);

  p.x = x;
  p.y = y;
  p.zinv = zinv;
  p.pos3d = vertex.position;
  // p.illumination = illumination;
}

//Could be sped up with Bresenham's line algorithm
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result)
{
  int N = result.size();
  vec2 step = vec2(b-a) / float(max(N-1,1));
  vec2 current( a );
  for( int i=0; i<N; ++i )
  {
      result[i].x = round(current.x);
      result[i].y = round(current.y);
      current += step;
  }
}

void InterpolatePixel(Pixel a, Pixel b, vector<Pixel>& result)
{
  int N = result.size();
  vec3 step;
  step.x = (b.x-a.x) / float(max(N-1,1));
  step.y = (b.y-a.y) / float(max(N-1,1));
  step.z = (b.zinv-a.zinv) / float(max(N-1,1));
  
  // vec3 illuminationStep = (b.illumination - a.illumination) / float(max(N-1,1));
  vec4 pos3dStep = ( (b.pos3d*b.zinv) - (a.pos3d*a.zinv) ) / float(max(N-1,1));

  vec3 current( float(a.x), float(a.y), a.zinv);
  // vec3 currentIllumination(a.illumination);
  vec4 currentpos3d( a.pos3d * a.zinv );
  for( int i=0; i<N; ++i )
  {
    result[i].x = round(current.x);
    result[i].y = round(current.y);
    result[i].zinv = current.z;
    // result[i].illumination = currentIllumination;
    result[i].pos3d = currentpos3d / current.z;
    // cout << currentIllumination.x << currentIllumination.y << currentIllumination.z << "\n";
    current += step;

    currentpos3d += pos3dStep;
    // cout << currentpos3d.x << "," << currentpos3d.y << "," << currentpos3d.z << "\n";
    // currentIllumination += illuminationStep;
  }
}

// void InterpolatePixel(Pixel a, Pixel b, vector<Pixel>& result)
// {
//   int N = result.size();
//   vec3 step;
//   step.x = (b.x-a.x) / float(max(N-1,1));
//   step.y = (b.y-a.y) / float(max(N-1,1));
//   step.z = (b.zinv-a.zinv) / float(max(N-1,1));


//   vec3 current( float(a.x), float(a.y), a.zinv );
//   for( int i=0; i<N; ++i )
//   {
//     result[i].x = round(current.x);
//     result[i].y = round(current.y);
//     result[i].zinv = current.z;

//     current += step;
//   }
// }

//This draws the polygon from some vertices
void DrawPolygon(screen* screen, const vector<Vertex>& vertices, vec4& cameraPos, float& focalLength, vec3 colour, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH], vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea, vec4& currentNormal, vec3& currentReflectance )
{
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
  DrawPolygonRows( screen, leftPixels, rightPixels, colour, depthBuffer, currentNormal, currentReflectance, lightPos, lightPower, indirectLightPowerPerArea );
}

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

      if(pixel.x < leftPixels[rowIndex].x) 
      { 
        leftPixels[rowIndex].x = pixel.x;
        leftPixels[rowIndex].zinv = pixel.zinv;
        // leftPixels[rowIndex].illumination = pixel.illumination; 
        leftPixels[rowIndex].pos3d = pixel.pos3d;
      }
      if(pixel.x > rightPixels[rowIndex].x) 
      { 
        rightPixels[rowIndex].x = pixel.x; 
        rightPixels[rowIndex].zinv = pixel.zinv;
        rightPixels[rowIndex].pos3d = pixel.pos3d;
        // rightPixels[rowIndex].illumination = pixel.illumination;
      }
    }
  }

}

void TestComputePolygonRows(screen* screen)
{
  vector<Pixel> vertexPixels(3);

  vertexPixels[0].x = 10;
  vertexPixels[0].y = 5;
  vertexPixels[0].zinv = 0.0f;

  vertexPixels[1].x = 5;
  vertexPixels[1].y = 10;
  vertexPixels[1].zinv = 0.0f;

  vertexPixels[1].x = 5;
  vertexPixels[1].y = 10;
  vertexPixels[1].zinv = 0.0f;

  // vertexPixels[2] = {15,15, 0.0f};
  vertexPixels[2].x = 15;
  vertexPixels[2].y = 15;
  vertexPixels[2].zinv = 0.0f;

  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;
  ComputePolygonRows( vertexPixels, leftPixels, rightPixels, screen );
  
  for( int row=0; row<leftPixels.size(); ++row )
  {
    cout << "Start: ("
    << leftPixels[row].x << ","
    << leftPixels[row].y << "). "
    << "End: ("
    << rightPixels[row].x << ","
    << rightPixels[row].y << "). " << endl;
  }
}

//Finds a line on screen between a and b of some colour
void FindLine( Pixel a, Pixel b, vector<Pixel>& lineToDraw)
{
  // Find difference between two vertices.
  Pixel delta;
  delta.x = glm::abs( a.x - b.x );
  delta.y = glm::abs( a.y - b.y );
  // delta.zinv = glm::abs( a.zinv - b.zinv );//not used
   
  //Find number of pixels needed.
  int numberOfPixels = glm::max( delta.x, delta.y ) + 1;

  //Resize lineToDraw for interpolations
  lineToDraw.resize(numberOfPixels);
  
  //Interpolate between two vertices
  InterpolatePixel( a, b, lineToDraw );
}

void DrawPolygonRows( screen* screen, 
                      const vector<Pixel>& leftPixels, 
                      const vector<Pixel>& rightPixels, 
                      vec3 colour, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH],
                      vec4& currentNormal, vec3& currentReflectance, 
                      vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea)
{
  vector<Pixel> pixelsBetween(SCREEN_HEIGHT);

  //For each entry in left/right pixels
  for (int i = 0; i < leftPixels.size(); i++)
  {
    int numberOfPixelsBetween = rightPixels[i].x - leftPixels[i].x + 1;

    pixelsBetween.resize(numberOfPixelsBetween);
    // cout << leftPixels[i].illumination.x << leftPixels[i].illumination.y << leftPixels[i].illumination.z << "\n";
    InterpolatePixel(leftPixels[i], rightPixels[i], pixelsBetween);
    // cout << leftPixels[i].illumination.x << leftPixels[i].illumination.y << leftPixels[i].illumination.z << "\n";

    for(int p = 0; p < pixelsBetween.size(); p++)
    {
      if(pixelsBetween[p].y < SCREEN_HEIGHT && pixelsBetween[p].y >= 0 && pixelsBetween[p].x < SCREEN_WIDTH && pixelsBetween[p].x >=0)
      {
        
        PixelShader(pixelsBetween[p], screen, depthBuffer, colour, currentNormal, currentReflectance, lightPos, lightPower, indirectLightPowerPerArea );

        // if(pixelsBetween[p].zinv > depthBuffer[pixelsBetween[p].y][pixelsBetween[p].x])
        // {
        //   PutPixelSDL(screen, pixelsBetween[p].x, pixelsBetween[p].y, colour);

        //   depthBuffer[pixelsBetween[p].y][pixelsBetween[p].x] = pixelsBetween[p].zinv;
        // }
      } 
    }
  }
}

void PixelShader(const Pixel& p, screen* screen, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH], 
                vec3 colour, vec4& currentNormal, vec3& currentReflectance, vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea )
{
  if(p.zinv > depthBuffer[p.y][p.x])
  {
    vec4 r = lightPos - p.pos3d;
    vec4 n = currentNormal;

    vec4 rNorm = normalize(r);
    vec4 nNorm = normalize(n);
    
    float A = 4 * M_PI * ( pow(glm::length(r),2) );
    vec3 B = vec3(lightPower.x/A,lightPower.y/A,lightPower.z/A);

    vec3 D = B * max(glm::dot (rNorm,nNorm), 0.0f);

    vec3 illumination = currentReflectance * (D + indirectLightPowerPerArea);

    // cout << illumination.x << "," << illumination.y << "," << illumination.z << "\n";

    // PutPixelSDL(screen, p.x, p.y, p.illumination);
    PutPixelSDL(screen, p.x, p.y, illumination);

    depthBuffer[p.y][p.x] = p.zinv;
  }
}


