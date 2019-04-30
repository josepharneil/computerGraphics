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

//============= Defines =============//
#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 800//256
#define FULLSCREEN_MODE false
#define PI 3.14159265
#define LIGHT_POWER 5.0f
#define NEAR_CLIP 0.5f
#define FAR_CLIP 5.0f
#define ANGLE_OF_VIEW 3.14159265/2  //field of view is 90 deg as long as focal length is half screen dimension 
#define AMBIENT_POWER 0.3f
#define TEXTURE_SIZE 300


//============= Global Variables =============//
bool quit;
bool isWireframe;
bool showNormalMap;


//============= Structures =============//
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
  vec2 textureCoordinates;
};

struct Vertex
{
  vec4 position;
  vec2 textureCoordinates;
};

struct Image
{
  vec3 pixels[TEXTURE_SIZE][TEXTURE_SIZE];
};

//============= Function Definitions =============//
#pragma region FunctionDefinitions
void Update(vec4& cameraPos, int& yaw, mat4& cameraMatrix, vec4& lightPos);
void Draw(screen* screen, vector<Triangle>& triangles, vec4& cameraPos, float& focalLength, vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea);
void TestComputePolygonRows(screen* screen);
void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result);

//pixel versions
void InterpolatePixel(Pixel a, Pixel b, vector<Pixel>& result);
void DrawPolygonRows( screen* screen, const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH], vec4& currentNormal, vec4& currentTangent, vec4& currentBitangent, vec3& currentReflectance, vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea, const string& textureName );
void FindLine( Pixel a, Pixel b, vector<Pixel>& lineToDraw);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels, screen* screen );
void DrawPolygon(screen* screen, const vector<Vertex>& vertices, vec4& cameraPos, float& focalLength, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH], vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea, vec4& currentNormal, vec4& currentTangent, vec4& currentBitangent, vec3& currentReflectance, const string& textureName );
void PerspectiveProject( const Vertex& vertex, Pixel& p, vec4& cameraPos, float& focalLength );
void PixelShader(const Pixel& p, screen* s, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH], vec4& currentNormal, vec4& currentTangent, vec4& currentBitangent, vec3& currentReflectance, vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea, const string& textureName );
vec4 NormaliseNoHomogenous(vec4 vector4);
void CalculateCameraMatrix(vec4& camPos, int& yaw, mat4& camMatrix);

vector<Triangle> Clip(Triangle& triangle);
vector<Triangle> Triangulate(vector<Vertex> vertices, const vec3 color);
float DotNoHomogenous(const vec4 A, const vec4 B);
float LengthNoHomogenous(const vec4 v);
vec4 ReflectNoHomogenous(vec4 i, vec4 n);
void ClipToPlane(vector<Vertex>& inputVertices, vec4 planePoint, vec4 planeNormal);
int decodePNG(std::vector<unsigned char>& out_image, unsigned long& image_width, unsigned long& image_height, const unsigned char* in_png, size_t in_size, bool convert_to_rgba32);
void loadFile(std::vector<unsigned char>& buffer, const std::string& filename);
vec3 getPixelRGB(vector<unsigned char> image,int x, int y);

vec3 CheckerBoard(const float& x, const float& y);
void LoadTexture(Image& imageStruct, const string& textureNameString, bool isNormalize = true);
void CreateCoordinateSystem(const vec3 &N, vec3 &Nt, vec3 &Nb) ;

Image woodAlbedo;
Image woodSpecular;
Image woodNormal;



#pragma endregion FunctionDefinitions
//============= END Function Definitions =============//

//============= Overrides =============//
#pragma region Overrides
//Print vec4s
std::ostream &operator<<( std::ostream &os, vec4 const &v )
{
  return os << "(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w
            << ")";
}

std::ostream &operator<<( std::ostream &os, vec3 const &v )
{
  return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

std::ostream &operator<<( std::ostream &os, mat4 const &m )
{
  glm::mat4 mt = transpose( m );
  return os << mt[ 0 ] << endl
            << mt[ 1 ] << endl
            << mt[ 2 ] << endl
            << mt[ 3 ];
}
#pragma endregion
//============= END Overrides =============//

vec3 getPixelRGB(vector<unsigned char> image, int x, int y)
{
  if(x > TEXTURE_SIZE || y > TEXTURE_SIZE || x <0 || y<0)
  {
    cout << "pixel value out of range";
    cout << "x" << x << "\n";
    cout << "y" << y << "\n";
    return vec3(0,0,0);
  }
  int i = (x * 4) + (y * TEXTURE_SIZE * 4) ;
  int r = (int) image[i];
  int g = (int) image[i+1];
  int b = (int) image[i+2];
  return vec3(r,g,b);
}

//Will load texture with given name into the given (global) image struct
void LoadTexture(Image& imageStruct, const string& textureNameString, bool isNormalize)
{
  float div;
  if(isNormalize)
  {
    div = 255.0f;
  }
  else
  {
    div = 1.0f;
  }
  
  //load wood albedo into vector<unsigned char> woodAlbedo
  std::vector<unsigned char> picoImage;
  unsigned long dimension = TEXTURE_SIZE;
  std::vector<unsigned char> buffer;
  loadFile(buffer, textureNameString);

  //Check for error and print if so
  int error = decodePNG(picoImage, dimension, dimension, buffer.empty() ? 0 : &buffer[0], (unsigned long)buffer.size(), true);
  if(error != 0) std::cout << "error: " << error << std::endl;
  
  
  // picoImage is a vector of bytes, each byte an RGBA value,
  // so vector: RGBARGBARGBA...

  //vec3[300][300]   300x300 array of (R,G,B)
  for(int r = 0; r < TEXTURE_SIZE; r++)
  {
    for(int c = 0; c < TEXTURE_SIZE; c++)
    {
      imageStruct.pixels[r][c] = getPixelRGB(picoImage,r,c)/div;
    }
  }
}

void CreateCoordinateSystem(const vec3 &N, vec3 &Nt, vec3 &Nb) 
{ 
    if (fabs(N.x) > fabs(N.y)) 
    {
      Nt = vec3(N.z, 0.0f, -N.x) / sqrtf(N.x * N.x + N.z * N.z);
    }
    else
    {
      Nt = vec3(0.0f, -N.z, N.y) / sqrtf(N.y * N.y + N.z * N.z); 
    }

    
    
    Nb = glm::cross(N,Nt);
}


//============= Main =============//
int main( int argc, char* argv[] )
{
  //Load in textures
  LoadTexture(woodAlbedo, "woodAlbedo.png");
  LoadTexture(woodSpecular,"woodSpecular.png");
  LoadTexture(woodNormal,"woodNormal.png",false);

  

  //Initially, do not quit
  quit = false;
  showNormalMap = true;


  //Triangulation testing
  // Vertex A = {.position = vec4(0,0,0,1)};
  // Vertex B = {.position = vec4(1,1,1,1)};
  // Vertex C = {.position = vec4(2,2,2,1)};
  // Vertex D = {.position = vec4(3,3,3,1)};
  // Vertex E = {.position = vec4(4,4,4,1)};
  
  // vector<Vertex> myVertices;
  // myVertices.push_back(A);
  // myVertices.push_back(B);
  // myVertices.push_back(C);
  // myVertices.push_back(D);
  // myVertices.push_back(E);
  // Triangulate(myVertices, vec3(1,1,1));

  isWireframe = false;
  
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
  CalculateCameraMatrix(cameraPos, yaw, cameraMatrix);

  //Light
  vec4 lightPos(0,-0.5,-0.7,1);
  vec4 originalLightPos( 0.0f, -0.5f, -0.7f, 1.0f );
  vec3 lightPower = LIGHT_POWER*vec3( 1, 1, 1 );
  vec3 indirectLightPowerPerArea = AMBIENT_POWER *vec3( 1, 1, 1 );

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

      // w divide
      lightPos = invCameraMatrix * originalLightPos;

      // for(int i = 0; i < SCREEN_WIDTH; i++ )
      // {
      //   for(int j = 0; j < SCREEN_HEIGHT; j++)
      //   {
      //     // float x = SCREEN_WIDTH 
      //     // float i = SCREEN_WIDTH/
      //     PutPixelSDL(screen, i, j, CheckerBoard(0.001f*i,0.001f*j));
      //   }
      // }

      Draw(screen, triangles, cameraPos, focalLength, lightPos, lightPower, indirectLightPowerPerArea);
      SDL_Renderframe(screen);
      //SDL_SaveImage( screen, "screenshot.bmp" );
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}
//============= End Main =============//

//============= Draw =============//
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

  
  //============= Clipping =============//
  //vector to hold all clipped triangles
  vector<Triangle> clippedTriangles;
	// clippedTriangles.clear();
	// clippedTriangles.reserve( 40 );

  //Apply clipping to each triangle
  for( uint32_t i=0; i<triangles.size(); ++i )
  {
    //temp holds the resulting triangles from clipping a single triangle
    vector<Triangle> temp = Clip(triangles[i]);

    //each triangle resultant from clipping
    for( uint32_t j=0; j<temp.size(); ++j )
    {
      //add each triangle from the result vector to clippedTriangles
      clippedTriangles.push_back(temp[j]);
    }
  }
  //============= END Clipping =============//

  //For each clipped triangle
  for( uint32_t i=0; i<clippedTriangles.size(); ++i )
  {
    //Construct vector of triangle vertices
    vector<Vertex> vertices(3);
    
    //Triangle vertices
    vertices[0].position = clippedTriangles[i].v0;
    vertices[1].position = clippedTriangles[i].v1;
    vertices[2].position = clippedTriangles[i].v2;

    //Triangle UV coordinates
    vertices[0].textureCoordinates = clippedTriangles[i].t0;
    vertices[1].textureCoordinates = clippedTriangles[i].t1;
    vertices[2].textureCoordinates = clippedTriangles[i].t2;

    //Triangle local coordinate axes
    vec4 currentNormal = NormaliseNoHomogenous(clippedTriangles[i].normal);
    vec4 currentTangent;// = clippedTriangles[i].normal;//DEBUGGING: eventually calculate this
    vec4 currentBitangent;// = clippedTriangles[i].normal;//DEBUGGING: eventually calculate this

    //============= Tangent Bitangent =============//
    vec3 vertex0V3 = vec3(vertices[0].position.x,vertices[0].position.y,vertices[0].position.z);
    vec3 vertex1V3 = vec3(vertices[1].position.x,vertices[1].position.y,vertices[1].position.z);
    vec3 vertex2V3 = vec3(vertices[2].position.x,vertices[2].position.y,vertices[2].position.z);

    vec3 deltaPos1 = vertex1V3 - vertex0V3;
    vec3 deltaPos2 = vertex2V3 - vertex0V3;

    vec2 deltaUV1 = vertices[1].textureCoordinates - vertices[0].textureCoordinates;
    vec2 deltaUV2 = vertices[2].textureCoordinates - vertices[0].textureCoordinates;

    float r = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV1.y * deltaUV2.x);
    vec3 tangentV3 = (deltaPos1 * deltaUV2.y - deltaPos2 * deltaUV1.y)*r;
    vec3 bitangentV3 = (deltaPos2 * deltaUV1.x - deltaPos1 * deltaUV2.x)*r;

    currentTangent = vec4(tangentV3.x,tangentV3.y,tangentV3.z,1.0f);
    currentBitangent = vec4(bitangentV3.x,bitangentV3.y,bitangentV3.z,1.0f);

    //============= END Tangent Bitangent =============//

    //Triangle reflectance/ texture name
    vec3 currentReflectance = clippedTriangles[i].color;
    string textureName =  clippedTriangles[i].textureName;

    //Draw polygon for each triangle
    DrawPolygon( screen, vertices, cameraPos, focalLength, depthBuffer, 
    lightPos, lightPower, indirectLightPowerPerArea, currentNormal, currentTangent, currentBitangent, currentReflectance, textureName );
  }
}
//============= END Draw =============//

/*
For all triangles
  DrawPolygon -> compute and draw polygon rows
*/

//============= Update =============//
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
    if ( e.key.keysym.scancode == SDL_SCANCODE_R )
    {
      lightPos.y -= movementSpeed;
    }
    if ( e.key.keysym.scancode == SDL_SCANCODE_F )
    {
      lightPos.y += movementSpeed;
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

    if ( e.key.keysym.scancode == SDL_SCANCODE_N )
    {
      isWireframe = !isWireframe;
    }

    if (e.key.keysym.scancode == SDL_SCANCODE_G)
    {
      showNormalMap = !showNormalMap;
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
  
    CalculateCameraMatrix(cameraPos, yaw, cameraMatrix);
  }
}

/////////
//This draws the polygon from some vertices
void DrawPolygon(screen* screen, const vector<Vertex>& vertices, vec4& cameraPos, float& focalLength, 
                 float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH], vec4& lightPos, vec3& lightPower, 
                 vec3& indirectLightPowerPerArea, vec4& currentNormal, vec4& currentTangent, vec4& currentBitangent, vec3& currentReflectance, const string& textureName )
{

  //Find number of vertices of polygon (3 for triangle)
  int V = vertices.size();

  //Initialise 3 long vector of projected vertices (pixel locations)
  vector<Pixel> vertexPixels( V );

  //For each vertex
  for( int i=0; i<V; ++i )
  {
    //Compute projection
    PerspectiveProject( vertices[i], vertexPixels[i], cameraPos, focalLength);
  }
  
  //Initialise vectors to store left-most and right-most positions of each row of the projected triangle
  vector<Pixel> leftPixels;
  vector<Pixel> rightPixels;

  //Calculates leftPixels and rightPixels
  ComputePolygonRows( vertexPixels, leftPixels, rightPixels, screen );

  //Draws the rows
  DrawPolygonRows( screen, leftPixels, rightPixels, depthBuffer, currentNormal, currentTangent, currentBitangent, currentReflectance, lightPos, lightPower, indirectLightPowerPerArea, textureName );
}


//Projects scene point (triangle vertex) onto image plane
void PerspectiveProject( const Vertex& vertex, Pixel& p, vec4& cameraPos, float& focalLength)
{
  //Translates vertex so that camera is at origin, relative to the vertex
  vec4 vertexTransformed = vertex.position; //- cameraPos;

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

  p.x = x;
  p.y = y;
  p.zinv = zinv;
  p.pos3d = vertex.position;
  p.textureCoordinates = vertex.textureCoordinates;
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
        leftPixels[rowIndex].textureCoordinates = pixel.textureCoordinates;
        leftPixels[rowIndex].pos3d = pixel.pos3d;
      }
      if(pixel.x > rightPixels[rowIndex].x) 
      { 
        rightPixels[rowIndex].x = pixel.x; 
        rightPixels[rowIndex].zinv = pixel.zinv;
        rightPixels[rowIndex].textureCoordinates = pixel.textureCoordinates;
        rightPixels[rowIndex].pos3d = pixel.pos3d;
      }
    }
  }
}

void DrawPolygonRows( screen* screen, 
                      const vector<Pixel>& leftPixels, 
                      const vector<Pixel>& rightPixels, 
                      float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH],
                      vec4& currentNormal, vec4& currentTangent, vec4& currentBitangent, vec3& currentReflectance, 
                      vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea, const string& textureName)
{
  if(isWireframe)
  {  
    for(int i = 0; i < leftPixels.size(); i++)
    {
      PutPixelSDL(screen, leftPixels[i].x, leftPixels[i].y, vec3(1,1,1));
    }
    for(int i = 0; i < rightPixels.size(); i++)
    {
      PutPixelSDL(screen, rightPixels[i].x, rightPixels[i].y, vec3(1,1,1));
    }
  }
  else
  {
    vector<Pixel> pixelsBetween(SCREEN_HEIGHT);

    //For each entry in left/right pixels
    for (int i = 0; i < leftPixels.size(); i++)
    {
      int numberOfPixelsBetween = rightPixels[i].x - leftPixels[i].x + 1;

      pixelsBetween.resize(numberOfPixelsBetween);
      InterpolatePixel(leftPixels[i], rightPixels[i], pixelsBetween);

      for(int p = 0; p < pixelsBetween.size(); p++)
      {
        if(pixelsBetween[p].y < SCREEN_HEIGHT && pixelsBetween[p].y >= 0 && pixelsBetween[p].x < SCREEN_WIDTH && pixelsBetween[p].x >=0)
        {
          PixelShader(pixelsBetween[p], screen, depthBuffer, currentNormal, currentTangent, currentBitangent, currentReflectance, lightPos, lightPower, indirectLightPowerPerArea, textureName );
        } 
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


//Could be sped up with Bresenham's line algorithm
void InterpolatePixel(Pixel a, Pixel b, vector<Pixel>& result)
{
  int N = result.size();

  //Calc steps
  float xStep = (b.x-a.x) / float(max(N-1,1));
  float yStep = (b.y-a.y) / float(max(N-1,1));
  float zStep = (b.zinv-a.zinv) / float(max(N-1,1));
  //The quantity we interpolate is pos3dStep over Z - this accounts for projection
  vec4 pos3dStep = ( (b.pos3d*b.zinv) - (a.pos3d*a.zinv) ) / float(max(N-1,1));
  vec2 textureCoordinatesStep = ( (b.textureCoordinates*b.zinv) - (a.textureCoordinates*a.zinv) ) / float(max(N-1,1));

  //Initialise
  float xCurrent = float(a.x);
  float yCurrent = float(a.y);
  float zCurrent = a.zinv;
  vec4 pos3dCurrent( a.pos3d * a.zinv );
  vec2 textureCoordinatesCurrent( a.textureCoordinates * a.zinv );

  //Interpolate
  for( int i=0; i<N; ++i )
  {
    result[i].x = round(xCurrent);
    result[i].y = round(yCurrent);
    result[i].zinv = zCurrent;
    result[i].pos3d = pos3dCurrent / zCurrent;
    result[i].textureCoordinates = textureCoordinatesCurrent / zCurrent;
    
    xCurrent += xStep;
    yCurrent += yStep;
    zCurrent += zStep;
    pos3dCurrent += pos3dStep;
    textureCoordinatesCurrent += textureCoordinatesStep;
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


//implements the phong model
void PixelShader(const Pixel& p, screen* screen, float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH], 
                vec4& currentNormal, vec4& currentTangent, vec4& currentBitangent, vec3& currentReflectance, 
                vec4& lightPos, vec3& lightPower, vec3& indirectLightPowerPerArea, const string& textureName )
{
  //defaults
  //default diffuse is the argument currentReflectance
  float k_s = 0.0f; //default specularity
  vec4 N = NormaliseNoHomogenous(currentNormal); //default normal
  bool isMetallic = false;

  //point to shade is p.pos3d  
  if(p.zinv > depthBuffer[p.y][p.x])
  {
    //Check for texture name
    if(textureName == "checkerBoard")
    {
      currentReflectance = CheckerBoard(p.textureCoordinates.x, p.textureCoordinates.y);
    }
    if(textureName == "wood")
    {
      float u = p.textureCoordinates.x;
      float v = p.textureCoordinates.y;

      if(u < 0 || u >= 1 || v < 0 || v >= 1)
      {
        //currentReflectance = vec3(0,0,0);
        u = 0;
        v = 0;
      }
      //if UV in range
      else
      {
        int x = floor(u * TEXTURE_SIZE);
        int y = floor(v * TEXTURE_SIZE);        
        currentReflectance = woodAlbedo.pixels[x][y];
        k_s = length(woodSpecular.pixels[x][y]) * 0.1;
        
        if (showNormalMap) 
        {
          vec3 normalMapPixels = woodNormal.pixels[x][y];
          vec3 mapNormal = vec3(2.0f*(normalMapPixels.x/255.0f)-1.0f,2.0f*(normalMapPixels.y/255.0f)-1.0f,2.0f*(normalMapPixels.z/255.0f)-1.0f); //normalise using 2*(color/255) -1
          vec3 currentNormalV3 = normalize(vec3(currentNormal.x,currentNormal.y,currentNormal.z)); //the surface normal
          vec3 currentTangentV3 = normalize(vec3(currentTangent.x,currentTangent.y,currentTangent.z));
          vec3 currentBitangentV3 = normalize(vec3(currentBitangent.x,currentBitangent.y,currentBitangent.z));
          mat3 TBN = mat3(currentTangentV3,currentBitangentV3,currentNormalV3);
          //transform the normal from the normal map to be oriented with the surface using the TBN matrix
          vec3 NVec3 = TBN * mapNormal;
          N = vec4(NVec3.x,NVec3.y,NVec3.z,1.0f);
          N = NormaliseNoHomogenous(N);
        }

      }
    }
    
    //Vectors (all normalised)

    //vector from p to light
    vec4 L = NormaliseNoHomogenous(lightPos - p.pos3d);
    
    //perfect reflection direction 
    vec4 R = NormaliseNoHomogenous(ReflectNoHomogenous(-L,N));
    //vector from p to camera
    vec4 V = NormaliseNoHomogenous(-p.pos3d);

    
    //Parameters {note: if you decrease alpha, you should decrease k_s, so things look sensible}
    
    //hacky way to make just blue block shiny
    if(currentReflectance == vec3(0.15f, 0.15f, 0.75f ))
    {
      k_s = 0.3f;
      isMetallic = true;
    }
    // if(textureName == "wood")
    // {
    //   k_s = 0.15f;
    // }
    //diffuse constant
    float k_d = 1.0f;
    //shininess constant - controls size of specular highlight
    float alpha = 10.0f;
    //diffuse falloff constant
    float falloff = 2.0f; //this was 2.0f originally


    //diffuse
    float A = 4 * M_PI * ( pow(LengthNoHomogenous(lightPos - p.pos3d),falloff));
    float B = LIGHT_POWER / A;
    float diffuse = B * max(DotNoHomogenous(L,N), 0.0f);

    //specular
    //phong term must be in range [0,1]
    float phongTerm = pow(DotNoHomogenous(V,R),alpha);
    if(DotNoHomogenous(V,R) <= 0)
    {
      phongTerm = 0;
    }
    //backface cull of phong
    if(DotNoHomogenous(V,N)<0)
    {
      phongTerm = 0;
    }
    float specular = LIGHT_POWER * phongTerm;

    //shading = diffuse + specular + ambient {note: specular should not be affected by material color}
    vec3 diffuseShading = k_d * diffuse * currentReflectance;
    vec3 ambientShading = indirectLightPowerPerArea * currentReflectance;
    vec3 specularShading;// = k_s * specular * vec3(1,1,1);
    if(isMetallic)
    {
      specularShading = k_s * specular * currentReflectance;
    }
    else
    {
      specularShading = k_s * specular * vec3(1,1,1);
    }
    
    vec3 shading = diffuseShading + ambientShading + specularShading;
    PutPixelSDL(screen, p.x, p.y, shading);

    depthBuffer[p.y][p.x] = p.zinv;

  }
}

//Helpful Functions
vec4 NormaliseNoHomogenous(vec4 vector4)
{
  vec3 tempVec3 = vec3(vector4.x,vector4.y,vector4.z);
  tempVec3 = normalize(tempVec3);
  return vec4(tempVec3.x,tempVec3.y,tempVec3.z,1.0f);
}

void CalculateCameraMatrix(vec4& camPos, int& yaw, mat4& camMatrix)
{
  camMatrix[0][0] = cos( yaw * PI / 180 );
  camMatrix[0][1] = 0;
  camMatrix[0][2] = sin( yaw * PI / 180 );
  camMatrix[0][3] = 0;

  camMatrix[1][0] = 0;
  camMatrix[1][1] = 1;
  camMatrix[1][2] = 0;
  camMatrix[1][3] = 0;

  camMatrix[2][0] = -sin( yaw * PI / 180 );
  camMatrix[2][1] = 0;
  camMatrix[2][2] = cos( yaw * PI / 180 );
  camMatrix[2][3] = 0;

  camMatrix[3][0] = camPos.x;
  camMatrix[3][1] = camPos.y;
  camMatrix[3][2] = camPos.z;
  camMatrix[3][3] = 1;
}

//clip triangle, and return a vector of the resulting triangles
vector<Triangle> Clip(Triangle& triangle)
{
  vector<Vertex> vertices(3);
  vertices[0].position = triangle.v0;
  vertices[1].position = triangle.v1;
  vertices[2].position = triangle.v2;
  
  vertices[0].textureCoordinates = vec2(triangle.t0.x,triangle.t0.y);
  vertices[1].textureCoordinates = vec2(triangle.t1.x,triangle.t1.y);
  vertices[2].textureCoordinates = vec2(triangle.t2.x,triangle.t2.y);

  //until triangulation, we work with vec4s instead of Vertex structs, for ease

  //vertices is the working list of vertices. We first add the vertices from the input triangle. Then we clip to each plane of the view frustum in turn. Then we triangulate
  // vector<vec4> vertices;
  
  // //add vertices from input triangle
  // vertices.push_back(triangle.v0);
  // vertices.push_back(triangle.v1);
  // vertices.push_back(triangle.v2);

  float cosHalfAlpha = cosf((ANGLE_OF_VIEW)/2);
  float sinHalfAlpha = sinf((ANGLE_OF_VIEW)/2);

  //clip to near plane
  ClipToPlane(vertices,vec4(0,0,NEAR_CLIP,1),vec4(0,0,1,1));

  //clip to far plane
  ClipToPlane(vertices,vec4(0,0,FAR_CLIP,1),vec4(0,0,-1,1));

  //clip to top
  ClipToPlane(vertices,vec4(0,0,0,1), vec4(0.0f, cosHalfAlpha, sinHalfAlpha, 1.0f ) );

  //clip to right
  ClipToPlane(vertices,vec4(0,0,0,1), vec4( -cosHalfAlpha, 0.0f,sinHalfAlpha, 1.0f) );
  
  //clip to bottom
  ClipToPlane(vertices,vec4(0,0,0,1), vec4( 0.0f,-cosHalfAlpha,sinHalfAlpha, 1.0f) );

  //clip to left
  ClipToPlane(vertices,vec4(0,0,0,1), vec4(cosHalfAlpha, 0.0f,sinHalfAlpha, 1.0f) );

  
  //triangulate the polygon (might already be a triangle)
  vector<Triangle> result = Triangulate(vertices,triangle.color);

  for(int i = 0; i < result.size(); i++)
  {
    result[i].textureName = triangle.textureName;
  }
  

  return result;
}

//use fan triangulation to triangluate the given polygon (specified by ordered list of vertices), return a vector of triangles.
//color is only passed so that we can create triangles
vector<Triangle> Triangulate(vector<Vertex> vertices, const vec3 color)
{
  vector<Triangle> result;
  
  //if there are less than three points then return empty list.
  if(vertices.size() < 3)
  {
    return result;
  }
  
  //if the passed in polygon is a triagle, just place it in the list and return immediately.
  if(vertices.size() == 3)
  {
    result.push_back(Triangle(vertices[0].position,vertices[1].position,vertices[2].position,color,vertices[0].textureCoordinates,vertices[1].textureCoordinates,vertices[2].textureCoordinates));
    return result;
  }

  //otherwise, triangulate the polygon:
  
  //number of triangles to result from triangulation
  int numTriangles = (int)(ceil(((float)vertices.size()) /2.0f));

  for(int i = 0; i<numTriangles; i++ )
  {
    result.push_back(Triangle(vertices[0].position,vertices[i+1].position,vertices[i+2].position,color,vertices[0].textureCoordinates,vertices[i+1].textureCoordinates,vertices[i+2].textureCoordinates));
  }
  
  return result;
}

float DotNoHomogenous(const vec4 A, const vec4 B)
{
  vec3 tempA = vec3(A.x,A.y,A.z);
  vec3 tempB = vec3(B.x,B.y,B.z);
  return glm::dot(tempA, tempB);
}

float LengthNoHomogenous(const vec4 v)
{
  vec3 temp = vec3(v.x,v.y,v.z);
  return glm::length(temp);
}

vec4 ReflectNoHomogenous(vec4 i, vec4 n)
{
  vec3 i3 = vec3(i.x,i.y,i.z);
  vec3 n3 = vec3(n.x,n.y,n.z);
  vec3 r3 = glm::reflect(i3,n3);
  vec4 r4 = vec4(r3.x,r3.y,r3.z,1);
  return r4;
}

//clips a convex polygon to the plane specified by the point planePoint and the normal planeNormal, returns a vector<vec4> specifying the resultant polygon
//note that the order of vertices is important to specifying the polygon
void ClipToPlane(vector<Vertex>& inputVertices, vec4 planePoint, vec4 planeNormal)
{
  if(inputVertices.size() < 3)
  {
    return;
  }

  vector<Vertex> result;

  int n = inputVertices.size();

  float pdot = 0;
  for(int i = 0; i<n;i++)
  {
    float dot = DotNoHomogenous(planeNormal, (inputVertices[i].position - planePoint));

    //if current vertex is in, and previous was out, then find intersection and add
    // OR if current vertex is out, and previous vertex was in, then find intersection and add
    // this cannot be called for the first vertex
    if(dot * pdot < 0)
    {
      float t = pdot/(pdot - dot);
      vec4 I = inputVertices[i-1].position + t * (inputVertices[i].position - inputVertices[i-1].position);//find intersection I
      vec2 uv = inputVertices[i-1].textureCoordinates + t * (inputVertices[i].textureCoordinates - inputVertices[i-1].textureCoordinates);//find intersection I
      Vertex newVertex;
      newVertex.position = I;
      newVertex.textureCoordinates = uv;
      result.push_back(newVertex);
    }
    
    //if current vertex is in (or on) the plane, then add to list
    if(dot >= 0)
    {
      Vertex newVertex;
      newVertex.position = inputVertices[i].position;
      newVertex.textureCoordinates = inputVertices[i].textureCoordinates;
      result.push_back(newVertex);
    }

    pdot = dot;
  }

  float idot = DotNoHomogenous(planeNormal,(inputVertices[0].position-planePoint));
  //check final edge (i.e the edge from )
  if (pdot * idot < 0)
  {
    float t = pdot/(pdot - idot);
    vec4 I = inputVertices[n-1].position + t*(inputVertices[0].position - inputVertices[n-1].position);
    vec2 uv = inputVertices[n-1].textureCoordinates + t*(inputVertices[0].textureCoordinates - inputVertices[n-1].textureCoordinates);
    Vertex newVertex;
    newVertex.position = I;
    newVertex.textureCoordinates = uv;
    result.push_back(newVertex);
  }

  inputVertices = result;

}
//presume (x,y) is a UV coord (i.e. in range [0,1])
vec3 CheckerBoard(const float& x, const float& y)
{ 
  //out of range
  if(x > 1|| y > 1)
  {
    return vec3(1,0,0);
  }

  //top left
  if(x<=0.5f && y<=0.5f)
  {
    return vec3(0.2f,0.2f,0.2f);
  }
  //top right
  if(x>0.5f && y<=0.5f)
  {
    return vec3(0.8f,0.8f,0.8f);
  }
  //bottom left
  if(x<=0.5f && y>0.5f)
  {
    return vec3(0.8f,0.8f,0.8f);
  }
  //bottom right
  if(x>0.5f && y>0.5f)
  {
    return vec3(0.2f,0.2f,0.2f);
  }

  // throw 
  return vec3(1,0,1);

}


//============= picoPNG =============//

#include <vector>

/*
decodePNG: The picoPNG function, decodes a PNG file buffer in memory, into a raw pixel buffer.
out_image: output parameter, this will contain the raw pixels after decoding.
  By default the output is 32-bit RGBA color.
  The std::vector is automatically resized to the correct size.
image_width: output_parameter, this will contain the width of the image in pixels.
image_height: output_parameter, this will contain the height of the image in pixels.
in_png: pointer to the buffer of the PNG file in memory. To get it from a file on
  disk, load it and store it in a memory buffer yourself first.
in_size: size of the input PNG file in bytes.
convert_to_rgba32: optional parameter, true by default.
  Set to true to get the output in RGBA 32-bit (8 bit per channel) color format
  no matter what color type the original PNG image had. This gives predictable,
  useable data from any random input PNG.
  Set to false to do no color conversion at all. The result then has the same data
  type as the PNG image, which can range from 1 bit to 64 bits per pixel.
  Information about the color type or palette colors are not provided. You need
  to know this information yourself to be able to use the data so this only
  works for trusted PNG files. Use LodePNG instead of picoPNG if you need this information.
return: 0 if success, not 0 if some error occured.
*/
int decodePNG(std::vector<unsigned char>& out_image, unsigned long& image_width, unsigned long& image_height, const unsigned char* in_png, size_t in_size, bool convert_to_rgba32 = true)
{
  // picoPNG version 20101224
  // Copyright (c) 2005-2010 Lode Vandevenne
  //
  // This software is provided 'as-is', without any express or implied
  // warranty. In no event will the authors be held liable for any damages
  // arising from the use of this software.
  //
  // Permission is granted to anyone to use this software for any purpose,
  // including commercial applications, and to alter it and redistribute it
  // freely, subject to the following restrictions:
  //
  //     1. The origin of this software must not be misrepresented; you must not
  //     claim that you wrote the original software. If you use this software
  //     in a product, an acknowledgment in the product documentation would be
  //     appreciated but is not required.
  //     2. Altered source versions must be plainly marked as such, and must not be
  //     misrepresented as being the original software.
  //     3. This notice may not be removed or altered from any source distribution.
  
  // picoPNG is a PNG decoder in one C++ function of around 500 lines. Use picoPNG for
  // programs that need only 1 .cpp file. Since it's a single function, it's very limited,
  // it can convert a PNG to raw pixel data either converted to 32-bit RGBA color or
  // with no color conversion at all. For anything more complex, another tiny library
  // is available: LodePNG (lodepng.c(pp)), which is a single source and header file.
  // Apologies for the compact code style, it's to make this tiny.
  
  static const unsigned long LENBASE[29] =  {3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258};
  static const unsigned long LENEXTRA[29] = {0,0,0,0,0,0,0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4,  4,  5,  5,  5,  5,  0};
  static const unsigned long DISTBASE[30] =  {1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577};
  static const unsigned long DISTEXTRA[30] = {0,0,0,0,1,1,2, 2, 3, 3, 4, 4, 5, 5,  6,  6,  7,  7,  8,  8,   9,   9,  10,  10,  11,  11,  12,   12,   13,   13};
  static const unsigned long CLCL[19] = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15}; //code length code lengths
  struct Zlib //nested functions for zlib decompression
  {
    static unsigned long readBitFromStream(size_t& bitp, const unsigned char* bits) { unsigned long result = (bits[bitp >> 3] >> (bitp & 0x7)) & 1; bitp++; return result;}
    static unsigned long readBitsFromStream(size_t& bitp, const unsigned char* bits, size_t nbits)
    {
      unsigned long result = 0;
      for(size_t i = 0; i < nbits; i++) result += (readBitFromStream(bitp, bits)) << i;
      return result;
    }
    struct HuffmanTree
    {
      int makeFromLengths(const std::vector<unsigned long>& bitlen, unsigned long maxbitlen)
      { //make tree given the lengths
        unsigned long numcodes = (unsigned long)(bitlen.size()), treepos = 0, nodefilled = 0;
        std::vector<unsigned long> tree1d(numcodes), blcount(maxbitlen + 1, 0), nextcode(maxbitlen + 1, 0);
        for(unsigned long bits = 0; bits < numcodes; bits++) blcount[bitlen[bits]]++; //count number of instances of each code length
        for(unsigned long bits = 1; bits <= maxbitlen; bits++) nextcode[bits] = (nextcode[bits - 1] + blcount[bits - 1]) << 1;
        for(unsigned long n = 0; n < numcodes; n++) if(bitlen[n] != 0) tree1d[n] = nextcode[bitlen[n]]++; //generate all the codes
        tree2d.clear(); tree2d.resize(numcodes * 2, 32767); //32767 here means the tree2d isn't filled there yet
        for(unsigned long n = 0; n < numcodes; n++) //the codes
        for(unsigned long i = 0; i < bitlen[n]; i++) //the bits for this code
        {
          unsigned long bit = (tree1d[n] >> (bitlen[n] - i - 1)) & 1;
          if(treepos > numcodes - 2) return 55;
          if(tree2d[2 * treepos + bit] == 32767) //not yet filled in
          {
            if(i + 1 == bitlen[n]) { tree2d[2 * treepos + bit] = n; treepos = 0; } //last bit
            else { tree2d[2 * treepos + bit] = ++nodefilled + numcodes; treepos = nodefilled; } //addresses are encoded as values > numcodes
          }
          else treepos = tree2d[2 * treepos + bit] - numcodes; //subtract numcodes from address to get address value
        }
        return 0;
      }
      int decode(bool& decoded, unsigned long& result, size_t& treepos, unsigned long bit) const
      { //Decodes a symbol from the tree
        unsigned long numcodes = (unsigned long)tree2d.size() / 2;
        if(treepos >= numcodes) return 11; //error: you appeared outside the codetree
        result = tree2d[2 * treepos + bit];
        decoded = (result < numcodes);
        treepos = decoded ? 0 : result - numcodes;
        return 0;
      }
      std::vector<unsigned long> tree2d; //2D representation of a huffman tree: The one dimension is "0" or "1", the other contains all nodes and leaves of the tree.
    };
    struct Inflator
    {
      int error;
      void inflate(std::vector<unsigned char>& out, const std::vector<unsigned char>& in, size_t inpos = 0)
      {
        size_t bp = 0, pos = 0; //bit pointer and byte pointer
        error = 0;
        unsigned long BFINAL = 0;
        while(!BFINAL && !error)
        {
          if(bp >> 3 >= in.size()) { error = 52; return; } //error, bit pointer will jump past memory
          BFINAL = readBitFromStream(bp, &in[inpos]);
          unsigned long BTYPE = readBitFromStream(bp, &in[inpos]); BTYPE += 2 * readBitFromStream(bp, &in[inpos]);
          if(BTYPE == 3) { error = 20; return; } //error: invalid BTYPE
          else if(BTYPE == 0) inflateNoCompression(out, &in[inpos], bp, pos, in.size());
          else inflateHuffmanBlock(out, &in[inpos], bp, pos, in.size(), BTYPE);
        }
        if(!error) out.resize(pos); //Only now we know the true size of out, resize it to that
      }
      void generateFixedTrees(HuffmanTree& tree, HuffmanTree& treeD) //get the tree of a deflated block with fixed tree
      {
        std::vector<unsigned long> bitlen(288, 8), bitlenD(32, 5);;
        for(size_t i = 144; i <= 255; i++) bitlen[i] = 9;
        for(size_t i = 256; i <= 279; i++) bitlen[i] = 7;
        tree.makeFromLengths(bitlen, 15);
        treeD.makeFromLengths(bitlenD, 15);
      }
      HuffmanTree codetree, codetreeD, codelengthcodetree; //the code tree for Huffman codes, dist codes, and code length codes
      unsigned long huffmanDecodeSymbol(const unsigned char* in, size_t& bp, const HuffmanTree& codetree, size_t inlength)
      { //decode a single symbol from given list of bits with given code tree. return value is the symbol
        bool decoded; unsigned long ct;
        for(size_t treepos = 0;;)
        {
          if((bp & 0x07) == 0 && (bp >> 3) > inlength) { error = 10; return 0; } //error: end reached without endcode
          error = codetree.decode(decoded, ct, treepos, readBitFromStream(bp, in)); if(error) return 0; //stop, an error happened
          if(decoded) return ct;
        }
      }
      void getTreeInflateDynamic(HuffmanTree& tree, HuffmanTree& treeD, const unsigned char* in, size_t& bp, size_t inlength)
      { //get the tree of a deflated block with dynamic tree, the tree itself is also Huffman compressed with a known tree
        std::vector<unsigned long> bitlen(288, 0), bitlenD(32, 0);
        if(bp >> 3 >= inlength - 2) { error = 49; return; } //the bit pointer is or will go past the memory
        size_t HLIT =  readBitsFromStream(bp, in, 5) + 257; //number of literal/length codes + 257
        size_t HDIST = readBitsFromStream(bp, in, 5) + 1; //number of dist codes + 1
        size_t HCLEN = readBitsFromStream(bp, in, 4) + 4; //number of code length codes + 4
        std::vector<unsigned long> codelengthcode(19); //lengths of tree to decode the lengths of the dynamic tree
        for(size_t i = 0; i < 19; i++) codelengthcode[CLCL[i]] = (i < HCLEN) ? readBitsFromStream(bp, in, 3) : 0;
        error = codelengthcodetree.makeFromLengths(codelengthcode, 7); if(error) return;
        size_t i = 0, replength;
        while(i < HLIT + HDIST)
        {
          unsigned long code = huffmanDecodeSymbol(in, bp, codelengthcodetree, inlength); if(error) return;
          if(code <= 15)  { if(i < HLIT) bitlen[i++] = code; else bitlenD[i++ - HLIT] = code; } //a length code
          else if(code == 16) //repeat previous
          {
            if(bp >> 3 >= inlength) { error = 50; return; } //error, bit pointer jumps past memory
            replength = 3 + readBitsFromStream(bp, in, 2);
            unsigned long value; //set value to the previous code
            if((i - 1) < HLIT) value = bitlen[i - 1];
            else value = bitlenD[i - HLIT - 1];
            for(size_t n = 0; n < replength; n++) //repeat this value in the next lengths
            {
              if(i >= HLIT + HDIST) { error = 13; return; } //error: i is larger than the amount of codes
              if(i < HLIT) bitlen[i++] = value; else bitlenD[i++ - HLIT] = value;
            }
          }
          else if(code == 17) //repeat "0" 3-10 times
          {
            if(bp >> 3 >= inlength) { error = 50; return; } //error, bit pointer jumps past memory
            replength = 3 + readBitsFromStream(bp, in, 3);
            for(size_t n = 0; n < replength; n++) //repeat this value in the next lengths
            {
              if(i >= HLIT + HDIST) { error = 14; return; } //error: i is larger than the amount of codes
              if(i < HLIT) bitlen[i++] = 0; else bitlenD[i++ - HLIT] = 0;
            }
          }
          else if(code == 18) //repeat "0" 11-138 times
          {
            if(bp >> 3 >= inlength) { error = 50; return; } //error, bit pointer jumps past memory
            replength = 11 + readBitsFromStream(bp, in, 7);
            for(size_t n = 0; n < replength; n++) //repeat this value in the next lengths
            {
              if(i >= HLIT + HDIST) { error = 15; return; } //error: i is larger than the amount of codes
              if(i < HLIT) bitlen[i++] = 0; else bitlenD[i++ - HLIT] = 0;
            }
          }
          else { error = 16; return; } //error: somehow an unexisting code appeared. This can never happen.
        }
        if(bitlen[256] == 0) { error = 64; return; } //the length of the end code 256 must be larger than 0
        error = tree.makeFromLengths(bitlen, 15); if(error) return; //now we've finally got HLIT and HDIST, so generate the code trees, and the function is done
        error = treeD.makeFromLengths(bitlenD, 15); if(error) return;
      }
      void inflateHuffmanBlock(std::vector<unsigned char>& out, const unsigned char* in, size_t& bp, size_t& pos, size_t inlength, unsigned long btype) 
      {
        if(btype == 1) { generateFixedTrees(codetree, codetreeD); }
        else if(btype == 2) { getTreeInflateDynamic(codetree, codetreeD, in, bp, inlength); if(error) return; }
        for(;;)
        {
          unsigned long code = huffmanDecodeSymbol(in, bp, codetree, inlength); if(error) return;
          if(code == 256) return; //end code
          else if(code <= 255) //literal symbol
          {
            if(pos >= out.size()) out.resize((pos + 1) * 2); //reserve more room
            out[pos++] = (unsigned char)(code);
          }
          else if(code >= 257 && code <= 285) //length code
          {
            size_t length = LENBASE[code - 257], numextrabits = LENEXTRA[code - 257];
            if((bp >> 3) >= inlength) { error = 51; return; } //error, bit pointer will jump past memory
            length += readBitsFromStream(bp, in, numextrabits);
            unsigned long codeD = huffmanDecodeSymbol(in, bp, codetreeD, inlength); if(error) return;
            if(codeD > 29) { error = 18; return; } //error: invalid dist code (30-31 are never used)
            unsigned long dist = DISTBASE[codeD], numextrabitsD = DISTEXTRA[codeD];
            if((bp >> 3) >= inlength) { error = 51; return; } //error, bit pointer will jump past memory
            dist += readBitsFromStream(bp, in, numextrabitsD);
            size_t start = pos, back = start - dist; //backwards
            if(pos + length >= out.size()) out.resize((pos + length) * 2); //reserve more room
            for(size_t i = 0; i < length; i++) { out[pos++] = out[back++]; if(back >= start) back = start - dist; }
          }
        }
      }
      void inflateNoCompression(std::vector<unsigned char>& out, const unsigned char* in, size_t& bp, size_t& pos, size_t inlength)
      {
        while((bp & 0x7) != 0) bp++; //go to first boundary of byte
        size_t p = bp / 8;
        if(p >= inlength - 4) { error = 52; return; } //error, bit pointer will jump past memory
        unsigned long LEN = in[p] + 256 * in[p + 1], NLEN = in[p + 2] + 256 * in[p + 3]; p += 4;
        if(LEN + NLEN != 65535) { error = 21; return; } //error: NLEN is not one's complement of LEN
        if(pos + LEN >= out.size()) out.resize(pos + LEN);
        if(p + LEN > inlength) { error = 23; return; } //error: reading outside of in buffer
        for(unsigned long n = 0; n < LEN; n++) out[pos++] = in[p++]; //read LEN bytes of literal data
        bp = p * 8;
      }
    };
    int decompress(std::vector<unsigned char>& out, const std::vector<unsigned char>& in) //returns error value
    {
      Inflator inflator;
      if(in.size() < 2) { return 53; } //error, size of zlib data too small
      if((in[0] * 256 + in[1]) % 31 != 0) { return 24; } //error: 256 * in[0] + in[1] must be a multiple of 31, the FCHECK value is supposed to be made that way
      unsigned long CM = in[0] & 15, CINFO = (in[0] >> 4) & 15, FDICT = (in[1] >> 5) & 1;
      if(CM != 8 || CINFO > 7) { return 25; } //error: only compression method 8: inflate with sliding window of 32k is supported by the PNG spec
      if(FDICT != 0) { return 26; } //error: the specification of PNG says about the zlib stream: "The additional flags shall not specify a preset dictionary."
      inflator.inflate(out, in, 2);
      return inflator.error; //note: adler32 checksum was skipped and ignored
    }
  };
  struct PNG //nested functions for PNG decoding
  {
    struct Info
    {
      unsigned long width, height, colorType, bitDepth, compressionMethod, filterMethod, interlaceMethod, key_r, key_g, key_b;
      bool key_defined; //is a transparent color key given?
      std::vector<unsigned char> palette;
    } info;
    int error;
    void decode(std::vector<unsigned char>& out, const unsigned char* in, size_t size, bool convert_to_rgba32)
    {
      error = 0;
      if(size == 0 || in == 0) { error = 48; return; } //the given data is empty
      readPngHeader(&in[0], size); if(error) return;
      size_t pos = 33; //first byte of the first chunk after the header
      std::vector<unsigned char> idat; //the data from idat chunks
      bool IEND = false, known_type = true;
      info.key_defined = false;
      while(!IEND) //loop through the chunks, ignoring unknown chunks and stopping at IEND chunk. IDAT data is put at the start of the in buffer
      {
        if(pos + 8 >= size) { error = 30; return; } //error: size of the in buffer too small to contain next chunk
        size_t chunkLength = read32bitInt(&in[pos]); pos += 4;
        if(chunkLength > 2147483647) { error = 63; return; }
        if(pos + chunkLength >= size) { error = 35; return; } //error: size of the in buffer too small to contain next chunk
        if(in[pos + 0] == 'I' && in[pos + 1] == 'D' && in[pos + 2] == 'A' && in[pos + 3] == 'T') //IDAT chunk, containing compressed image data
        {
          idat.insert(idat.end(), &in[pos + 4], &in[pos + 4 + chunkLength]);
          pos += (4 + chunkLength);
        }
        else if(in[pos + 0] == 'I' && in[pos + 1] == 'E' && in[pos + 2] == 'N' && in[pos + 3] == 'D')  { pos += 4; IEND = true; }
        else if(in[pos + 0] == 'P' && in[pos + 1] == 'L' && in[pos + 2] == 'T' && in[pos + 3] == 'E') //palette chunk (PLTE)
        {
          pos += 4; //go after the 4 letters
          info.palette.resize(4 * (chunkLength / 3));
          if(info.palette.size() > (4 * 256)) { error = 38; return; } //error: palette too big
          for(size_t i = 0; i < info.palette.size(); i += 4)
          {
            for(size_t j = 0; j < 3; j++) info.palette[i + j] = in[pos++]; //RGB
            info.palette[i + 3] = 255; //alpha
          }
        }
        else if(in[pos + 0] == 't' && in[pos + 1] == 'R' && in[pos + 2] == 'N' && in[pos + 3] == 'S') //palette transparency chunk (tRNS)
        {
          pos += 4; //go after the 4 letters
          if(info.colorType == 3)
          {
            if(4 * chunkLength > info.palette.size()) { error = 39; return; } //error: more alpha values given than there are palette entries
            for(size_t i = 0; i < chunkLength; i++) info.palette[4 * i + 3] = in[pos++];
          }
          else if(info.colorType == 0)
          {
            if(chunkLength != 2) { error = 40; return; } //error: this chunk must be 2 bytes for greyscale image
            info.key_defined = 1; info.key_r = info.key_g = info.key_b = 256 * in[pos] + in[pos + 1]; pos += 2;
          }
          else if(info.colorType == 2)
          {
            if(chunkLength != 6) { error = 41; return; } //error: this chunk must be 6 bytes for RGB image
            info.key_defined = 1;
            info.key_r = 256 * in[pos] + in[pos + 1]; pos += 2;
            info.key_g = 256 * in[pos] + in[pos + 1]; pos += 2;
            info.key_b = 256 * in[pos] + in[pos + 1]; pos += 2;
          }
          else { error = 42; return; } //error: tRNS chunk not allowed for other color models
        }
        else //it's not an implemented chunk type, so ignore it: skip over the data
        {
          if(!(in[pos + 0] & 32)) { error = 69; return; } //error: unknown critical chunk (5th bit of first byte of chunk type is 0)
          pos += (chunkLength + 4); //skip 4 letters and uninterpreted data of unimplemented chunk
          known_type = false;
        }
        pos += 4; //step over CRC (which is ignored)
      }
      unsigned long bpp = getBpp(info);
      std::vector<unsigned char> scanlines(((info.width * (info.height * bpp + 7)) / 8) + info.height); //now the out buffer will be filled
      Zlib zlib; //decompress with the Zlib decompressor
      error = zlib.decompress(scanlines, idat); if(error) return; //stop if the zlib decompressor returned an error
      size_t bytewidth = (bpp + 7) / 8, outlength = (info.height * info.width * bpp + 7) / 8;
      out.resize(outlength); //time to fill the out buffer
      unsigned char* out_ = outlength ? &out[0] : 0; //use a regular pointer to the std::vector for faster code if compiled without optimization
      if(info.interlaceMethod == 0) //no interlace, just filter
      {
        size_t linestart = 0, linelength = (info.width * bpp + 7) / 8; //length in bytes of a scanline, excluding the filtertype byte
        if(bpp >= 8) //byte per byte
        for(unsigned long y = 0; y < info.height; y++)
        {
          unsigned long filterType = scanlines[linestart];
          const unsigned char* prevline = (y == 0) ? 0 : &out_[(y - 1) * info.width * bytewidth];
          unFilterScanline(&out_[linestart - y], &scanlines[linestart + 1], prevline, bytewidth, filterType,  linelength); if(error) return;
          linestart += (1 + linelength); //go to start of next scanline
        }
        else //less than 8 bits per pixel, so fill it up bit per bit
        {
          std::vector<unsigned char> templine((info.width * bpp + 7) >> 3); //only used if bpp < 8
          for(size_t y = 0, obp = 0; y < info.height; y++)
          {
            unsigned long filterType = scanlines[linestart];
            const unsigned char* prevline = (y == 0) ? 0 : &out_[(y - 1) * info.width * bytewidth];
            unFilterScanline(&templine[0], &scanlines[linestart + 1], prevline, bytewidth, filterType, linelength); if(error) return;
            for(size_t bp = 0; bp < info.width * bpp;) setBitOfReversedStream(obp, out_, readBitFromReversedStream(bp, &templine[0]));
            linestart += (1 + linelength); //go to start of next scanline
          }
        }
      }
      else //interlaceMethod is 1 (Adam7)
      {
        size_t passw[7] = { (info.width + 7) / 8, (info.width + 3) / 8, (info.width + 3) / 4, (info.width + 1) / 4, (info.width + 1) / 2, (info.width + 0) / 2, (info.width + 0) / 1 };
        size_t passh[7] = { (info.height + 7) / 8, (info.height + 7) / 8, (info.height + 3) / 8, (info.height + 3) / 4, (info.height + 1) / 4, (info.height + 1) / 2, (info.height + 0) / 2 };
        size_t passstart[7] = {0};
        size_t pattern[28] = {0,4,0,2,0,1,0,0,0,4,0,2,0,1,8,8,4,4,2,2,1,8,8,8,4,4,2,2}; //values for the adam7 passes
        for(int i = 0; i < 6; i++) passstart[i + 1] = passstart[i] + passh[i] * ((passw[i] ? 1 : 0) + (passw[i] * bpp + 7) / 8);
        std::vector<unsigned char> scanlineo((info.width * bpp + 7) / 8), scanlinen((info.width * bpp + 7) / 8); //"old" and "new" scanline
        for(int i = 0; i < 7; i++)
          adam7Pass(&out_[0], &scanlinen[0], &scanlineo[0], &scanlines[passstart[i]], info.width, pattern[i], pattern[i + 7], pattern[i + 14], pattern[i + 21], passw[i], passh[i], bpp);
      }
      if(convert_to_rgba32 && (info.colorType != 6 || info.bitDepth != 8)) //conversion needed
      {
        std::vector<unsigned char> data = out;
        error = convert(out, &data[0], info, info.width, info.height);
      }
    }
    void readPngHeader(const unsigned char* in, size_t inlength) //read the information from the header and store it in the Info
    {
      if(inlength < 29) { error = 27; return; } //error: the data length is smaller than the length of the header
      if(in[0] != 137 || in[1] != 80 || in[2] != 78 || in[3] != 71 || in[4] != 13 || in[5] != 10 || in[6] != 26 || in[7] != 10) { error = 28; return; } //no PNG signature
      if(in[12] != 'I' || in[13] != 'H' || in[14] != 'D' || in[15] != 'R') { error = 29; return; } //error: it doesn't start with a IHDR chunk!
      info.width = read32bitInt(&in[16]); info.height = read32bitInt(&in[20]);
      info.bitDepth = in[24]; info.colorType = in[25];
      info.compressionMethod = in[26]; if(in[26] != 0) { error = 32; return; } //error: only compression method 0 is allowed in the specification
      info.filterMethod = in[27]; if(in[27] != 0) { error = 33; return; } //error: only filter method 0 is allowed in the specification
      info.interlaceMethod = in[28]; if(in[28] > 1) { error = 34; return; } //error: only interlace methods 0 and 1 exist in the specification
      error = checkColorValidity(info.colorType, info.bitDepth);
    }
    void unFilterScanline(unsigned char* recon, const unsigned char* scanline, const unsigned char* precon, size_t bytewidth, unsigned long filterType, size_t length)
    {
      switch(filterType)
      {
        case 0: for(size_t i = 0; i < length; i++) recon[i] = scanline[i]; break;
        case 1:
          for(size_t i =         0; i < bytewidth; i++) recon[i] = scanline[i];
          for(size_t i = bytewidth; i <    length; i++) recon[i] = scanline[i] + recon[i - bytewidth];
          break;
        case 2:
          if(precon) for(size_t i = 0; i < length; i++) recon[i] = scanline[i] + precon[i];
          else       for(size_t i = 0; i < length; i++) recon[i] = scanline[i];
          break;
        case 3:
          if(precon)
          {
            for(size_t i =         0; i < bytewidth; i++) recon[i] = scanline[i] + precon[i] / 2;
            for(size_t i = bytewidth; i <    length; i++) recon[i] = scanline[i] + ((recon[i - bytewidth] + precon[i]) / 2);
          }
          else
          {
            for(size_t i =         0; i < bytewidth; i++) recon[i] = scanline[i];
            for(size_t i = bytewidth; i <    length; i++) recon[i] = scanline[i] + recon[i - bytewidth] / 2;
          }
          break;
        case 4:
          if(precon)
          {
            for(size_t i =         0; i < bytewidth; i++) recon[i] = scanline[i] + paethPredictor(0, precon[i], 0);
            for(size_t i = bytewidth; i <    length; i++) recon[i] = scanline[i] + paethPredictor(recon[i - bytewidth], precon[i], precon[i - bytewidth]);
          }
          else
          {
            for(size_t i =         0; i < bytewidth; i++) recon[i] = scanline[i];
            for(size_t i = bytewidth; i <    length; i++) recon[i] = scanline[i] + paethPredictor(recon[i - bytewidth], 0, 0);
          }
          break;
        default: error = 36; return; //error: unexisting filter type given
      }
    }
    void adam7Pass(unsigned char* out, unsigned char* linen, unsigned char* lineo, const unsigned char* in, unsigned long w, size_t passleft, size_t passtop, size_t spacex, size_t spacey, size_t passw, size_t passh, unsigned long bpp)
    { //filter and reposition the pixels into the output when the image is Adam7 interlaced. This function can only do it after the full image is already decoded. The out buffer must have the correct allocated memory size already.
      if(passw == 0) return;
      size_t bytewidth = (bpp + 7) / 8, linelength = 1 + ((bpp * passw + 7) / 8);
      for(unsigned long y = 0; y < passh; y++)
      {
        unsigned char filterType = in[y * linelength], *prevline = (y == 0) ? 0 : lineo;
        unFilterScanline(linen, &in[y * linelength + 1], prevline, bytewidth, filterType, (w * bpp + 7) / 8); if(error) return;
        if(bpp >= 8) for(size_t i = 0; i < passw; i++) for(size_t b = 0; b < bytewidth; b++) //b = current byte of this pixel
          out[bytewidth * w * (passtop + spacey * y) + bytewidth * (passleft + spacex * i) + b] = linen[bytewidth * i + b];
        else for(size_t i = 0; i < passw; i++)
        {
          size_t obp = bpp * w * (passtop + spacey * y) + bpp * (passleft + spacex * i), bp = i * bpp;
          for(size_t b = 0; b < bpp; b++) setBitOfReversedStream(obp, out, readBitFromReversedStream(bp, &linen[0]));
        }
        unsigned char* temp = linen; linen = lineo; lineo = temp; //swap the two buffer pointers "line old" and "line new"
      }
    }
    static unsigned long readBitFromReversedStream(size_t& bitp, const unsigned char* bits) { unsigned long result = (bits[bitp >> 3] >> (7 - (bitp & 0x7))) & 1; bitp++; return result;}
    static unsigned long readBitsFromReversedStream(size_t& bitp, const unsigned char* bits, unsigned long nbits)
    {
      unsigned long result = 0;
      for(size_t i = nbits - 1; i < nbits; i--) result += ((readBitFromReversedStream(bitp, bits)) << i);
      return result;
    }
    void setBitOfReversedStream(size_t& bitp, unsigned char* bits, unsigned long bit) { bits[bitp >> 3] |=  (bit << (7 - (bitp & 0x7))); bitp++; }
    unsigned long read32bitInt(const unsigned char* buffer) { return (buffer[0] << 24) | (buffer[1] << 16) | (buffer[2] << 8) | buffer[3]; }
    int checkColorValidity(unsigned long colorType, unsigned long bd) //return type is a LodePNG error code
    {
      if((colorType == 2 || colorType == 4 || colorType == 6)) { if(!(bd == 8 || bd == 16)) return 37; else return 0; }
      else if(colorType == 0) { if(!(bd == 1 || bd == 2 || bd == 4 || bd == 8 || bd == 16)) return 37; else return 0; }
      else if(colorType == 3) { if(!(bd == 1 || bd == 2 || bd == 4 || bd == 8            )) return 37; else return 0; }
      else return 31; //unexisting color type
    }
    unsigned long getBpp(const Info& info)
    {
      if(info.colorType == 2) return (3 * info.bitDepth);
      else if(info.colorType >= 4) return (info.colorType - 2) * info.bitDepth;
      else return info.bitDepth;
    }
    int convert(std::vector<unsigned char>& out, const unsigned char* in, Info& infoIn, unsigned long w, unsigned long h)
    { //converts from any color type to 32-bit. return value = LodePNG error code
      size_t numpixels = w * h, bp = 0;
      out.resize(numpixels * 4);
      unsigned char* out_ = out.empty() ? 0 : &out[0]; //faster if compiled without optimization
      if(infoIn.bitDepth == 8 && infoIn.colorType == 0) //greyscale
      for(size_t i = 0; i < numpixels; i++)
      {
        out_[4 * i + 0] = out_[4 * i + 1] = out_[4 * i + 2] = in[i];
        out_[4 * i + 3] = (infoIn.key_defined && in[i] == infoIn.key_r) ? 0 : 255;
      }
      else if(infoIn.bitDepth == 8 && infoIn.colorType == 2) //RGB color
      for(size_t i = 0; i < numpixels; i++)
      {
        for(size_t c = 0; c < 3; c++) out_[4 * i + c] = in[3 * i + c];
        out_[4 * i + 3] = (infoIn.key_defined == 1 && in[3 * i + 0] == infoIn.key_r && in[3 * i + 1] == infoIn.key_g && in[3 * i + 2] == infoIn.key_b) ? 0 : 255;
      }
      else if(infoIn.bitDepth == 8 && infoIn.colorType == 3) //indexed color (palette)
      for(size_t i = 0; i < numpixels; i++)
      {
        if(4U * in[i] >= infoIn.palette.size()) return 46;
        for(size_t c = 0; c < 4; c++) out_[4 * i + c] = infoIn.palette[4 * in[i] + c]; //get rgb colors from the palette
      }
      else if(infoIn.bitDepth == 8 && infoIn.colorType == 4) //greyscale with alpha
      for(size_t i = 0; i < numpixels; i++)
      {
        out_[4 * i + 0] = out_[4 * i + 1] = out_[4 * i + 2] = in[2 * i + 0];
        out_[4 * i + 3] = in[2 * i + 1];
      }
      else if(infoIn.bitDepth == 8 && infoIn.colorType == 6) for(size_t i = 0; i < numpixels; i++) for(size_t c = 0; c < 4; c++) out_[4 * i + c] = in[4 * i + c]; //RGB with alpha
      else if(infoIn.bitDepth == 16 && infoIn.colorType == 0) //greyscale
      for(size_t i = 0; i < numpixels; i++)
      {
        out_[4 * i + 0] = out_[4 * i + 1] = out_[4 * i + 2] = in[2 * i];
        out_[4 * i + 3] = (infoIn.key_defined && 256U * in[i] + in[i + 1] == infoIn.key_r) ? 0 : 255;
      }
      else if(infoIn.bitDepth == 16 && infoIn.colorType == 2) //RGB color
      for(size_t i = 0; i < numpixels; i++)
      {
        for(size_t c = 0; c < 3; c++) out_[4 * i + c] = in[6 * i + 2 * c];
        out_[4 * i + 3] = (infoIn.key_defined && 256U*in[6*i+0]+in[6*i+1] == infoIn.key_r && 256U*in[6*i+2]+in[6*i+3] == infoIn.key_g && 256U*in[6*i+4]+in[6*i+5] == infoIn.key_b) ? 0 : 255;
      }
      else if(infoIn.bitDepth == 16 && infoIn.colorType == 4) //greyscale with alpha
      for(size_t i = 0; i < numpixels; i++)
      {
        out_[4 * i + 0] = out_[4 * i + 1] = out_[4 * i + 2] = in[4 * i]; //most significant byte
        out_[4 * i + 3] = in[4 * i + 2];
      }
      else if(infoIn.bitDepth == 16 && infoIn.colorType == 6) for(size_t i = 0; i < numpixels; i++) for(size_t c = 0; c < 4; c++) out_[4 * i + c] = in[8 * i + 2 * c]; //RGB with alpha
      else if(infoIn.bitDepth < 8 && infoIn.colorType == 0) //greyscale
      for(size_t i = 0; i < numpixels; i++)
      {
        unsigned long value = (readBitsFromReversedStream(bp, in, infoIn.bitDepth) * 255) / ((1 << infoIn.bitDepth) - 1); //scale value from 0 to 255
        out_[4 * i + 0] = out_[4 * i + 1] = out_[4 * i + 2] = (unsigned char)(value);
        out_[4 * i + 3] = (infoIn.key_defined && value && ((1U << infoIn.bitDepth) - 1U) == infoIn.key_r && ((1U << infoIn.bitDepth) - 1U)) ? 0 : 255;
      }
      else if(infoIn.bitDepth < 8 && infoIn.colorType == 3) //palette
      for(size_t i = 0; i < numpixels; i++)
      {
        unsigned long value = readBitsFromReversedStream(bp, in, infoIn.bitDepth);
        if(4 * value >= infoIn.palette.size()) return 47;
        for(size_t c = 0; c < 4; c++) out_[4 * i + c] = infoIn.palette[4 * value + c]; //get rgb colors from the palette
      }
      return 0;
    }
    unsigned char paethPredictor(short a, short b, short c) //Paeth predicter, used by PNG filter type 4
    {
      short p = a + b - c, pa = p > a ? (p - a) : (a - p), pb = p > b ? (p - b) : (b - p), pc = p > c ? (p - c) : (c - p);
      return (unsigned char)((pa <= pb && pa <= pc) ? a : pb <= pc ? b : c);
    }
  };
  PNG decoder; decoder.decode(out_image, in_png, in_size, convert_to_rgba32);
  image_width = decoder.info.width; image_height = decoder.info.height;
  return decoder.error;
}





//an example using the PNG loading function:

#include <iostream>
#include <fstream>

void loadFile(std::vector<unsigned char>& buffer, const std::string& filename) //designed for loading files from hard disk in an std::vector
{
  std::ifstream file(filename.c_str(), std::ios::in|std::ios::binary|std::ios::ate);

  //get filesize
  std::streamsize size = 0;
  if(file.seekg(0, std::ios::end).good()) size = file.tellg();
  if(file.seekg(0, std::ios::beg).good()) size -= file.tellg();

  //read contents of the file into the vector
  if(size > 0)
  {
    buffer.resize((size_t)size);
    file.read((char*)(&buffer[0]), size);
  }
  else buffer.clear();
}
/*
int main(int argc, char *argv[])
{
  const char* filename = argc > 1 ? argv[1] : "test.png";
  
  //load and decode
  std::vector<unsigned char> buffer, image;
  loadFile(buffer, filename);
  unsigned long w, h;
  int error = decodePNG(image, w, h, buffer.empty() ? 0 : &buffer[0], (unsigned long)buffer.size());
  
  //if there's an error, display it
  if(error != 0) std::cout << "error: " << error << std::endl;
  
  //the pixels are now in the vector "image", use it as texture, draw it, ...
  
  if(image.size() > 4) std::cout << "width: " << w << " height: " << h << " first pixel: " << std::hex << int(image[0]) << int(image[1]) << int(image[2]) << int(image[3]) << std::endl;
}
*/

/*
  //this is test code, it displays the pixels of a 1 bit PNG. To use it, set the flag convert_to_rgba32 to false and load a 1-bit PNG image with a small size (so that its ASCII representation can fit in a console window)
  for(int y = 0; y < h; y++)
  {
    for(int x = 0; x < w; x++)
    {
      int i = y * h + x;
      std::cout << (((image[i/8] >> (7-i%8)) & 1) ? '.' : '#');
    }
    std::cout << std::endl;
  }
*/
