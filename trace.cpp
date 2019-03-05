#include <stdio.h>
#include <GL/glut.h>
#include <math.h>
#include "global.h"
#include "sphere.h"

//
// Global variables
//
extern int win_width;
extern int win_height;

extern GLfloat frame[WIN_HEIGHT][WIN_WIDTH][3];

extern float image_width;
extern float image_height;

extern Point eye_pos;
extern float image_plane;
extern RGB_float background_clr;
extern RGB_float null_clr;

extern Spheres *scene;

// light 1 position and color
extern Point light1;
extern float light1_ambient[3];
extern float light1_diffuse[3];
extern float light1_specular[3];

// global ambient term
extern float global_ambient[3];

// light decay parameters
extern float decay_a;
extern float decay_b;
extern float decay_c;

extern int shadow_on;
extern int step_max;
extern int reflection_on;
extern int chessboard_on;
extern int refraction_on;
extern int stochastic_on;
extern int supersampling_on;

extern Point chessPoint;
extern Vector chessBoardNormal;

extern float ambient_black[3];
extern float ambient_white[3];
extern float diffuse_black[3];
extern float diffuse_white[3];
extern float specular_black[3];
extern float specular_white[3];
extern float shineness;
extern float reflectance;

//prototypes
RGB_float phongChessboard(Point, Vector, Vector, int);
RGB_float boardColor(Point, Vector);
bool intersectChessboard(Vector, Point, Point*);
Vector transmissionVector(Vector, Vector, Spheres*, Point*);
RGB_float recursive_ray_trace(Vector, Point, int);
RGB_float phong(Point, Vector, Vector, Spheres*);
void ray_trace();

/*********************************************************************
 * Phong ChessBoard - helper function
 *********************************************************************/
RGB_float phongChessboard(Point q, Vector v, Vector surf_norm, int base){
  RGB_float amb, spec, dif, col;
  //initialization
  amb = {0, 0, 0};
  spec = {0, 0, 0};
  dif = {0, 0, 0};
  col = {0, 0, 0};
  normalize(&surf_norm);
  float distanceLight; // distance between light source and point
  float diffuseParam;  // diffuse parameter
  float shininessParam; //specular shininess parameter
  float decay; //coefficient for light decay
  Vector Half;
  Vector LightSource;
  float ambient[3];
  float diffuse[3];
  float specular[3];
  if(base == 0){
    for(int i = 0; i < 3; i++){
      ambient[i] = ambient_black[i];
      diffuse[i] = diffuse_black[i];
      specular[i] = specular_black[i];
    }
  }
  else{
    for(int i = 0; i < 3; i++){
      ambient[i] = ambient_white[i];
      diffuse[i] = diffuse_white[i];
      specular[i] = specular_white[i];
    }
  }
  amb.r = amb.r + ambient[0]*2*global_ambient[0]*reflectance;
  amb.g = amb.g + ambient[1]*2*global_ambient[1]*reflectance;
  amb.b = amb.b + ambient[2]*2*global_ambient[2]*reflectance;
  Vector Shadow = get_vec(q, light1);
  if((shadow_on) && (isLightBlocked(q, Shadow, scene))){
    return amb;
  }
  LightSource = get_vec(q, light1);
  distanceLight = vec_len(LightSource);
  normalize(&LightSource);
  decay = (1.0/ ( decay_a + decay_b * distanceLight + decay_c * pow(distanceLight, 2)));
  diffuseParam = vec_dot(surf_norm, LightSource);
  dif.r = dif.r + light1_diffuse[0] * diffuse[0] * diffuseParam;
  dif.g = dif.g + light1_diffuse[1] * diffuse[1] * diffuseParam;
  dif.b = dif.b + light1_diffuse[2] * diffuse[2] * diffuseParam;
  float ref_scalar = 2.0*(vec_dot(surf_norm, LightSource));
  Half = vec_plus(vec_scale(surf_norm, ref_scalar), vec_scale(LightSource, -1.0));
  shininessParam = pow(vec_dot(v, Half), shineness);
  spec.r = spec.r + light1_specular[0] * specular[0] * shininessParam;
  spec.g = spec.g + light1_specular[1] * specular[1] * shininessParam;
  spec.b = spec.b + light1_specular[2] * specular[2] * shininessParam;
  col = clr_add(dif, spec);
  col = clr_scale(col, decay);
  col = clr_add(col, amb);
  return col;
}

/*********************************************************************
 * Board Color - helper function
 *********************************************************************/
RGB_float boardColor(Point hit, Vector u){
  int mod = 2;
  RGB_float retColor;
  bool determine1 = (int(hit.x) % mod == 0) && (int(hit.z) % mod == 0);
  bool determine2 = (int(hit.x) % mod != 0) && (int(hit.z) % mod != 0);
  if(determine1 || determine2){
    if(hit.x <= 0)
      retColor = phongChessboard(hit, u, chessBoardNormal, 0);
    else
      retColor = phongChessboard(hit, u, chessBoardNormal, 1);
  }
  else{
    if(hit.x <= 0)
      retColor = phongChessboard(hit, u, chessBoardNormal, 1);
    else
      retColor = phongChessboard(hit, u, chessBoardNormal, 0);
  }
  return retColor;
}

/*********************************************************************
 * intersect chessboard - helper function
 *********************************************************************/
bool intersectChessboard(Vector u, Point o, Point* hit){
  normalize(&chessBoardNormal);
  Vector originPoint = get_vec(chessPoint, o);
  float normal = vec_dot(chessBoardNormal, originPoint);
  float normalRay = vec_dot(chessBoardNormal, u);
  if(normalRay == 0)
    return false;
  if(normal == 0)
    return false;
  float r = normal / normalRay;
  if(r <= 0)
    return false;
  else{
    u = vec_scale(u, r);
    o = get_point(o, u);
    hit->x = o.x;
    hit->y = o.y;
    hit->z = o.z;
    if(o.x >= 4.0 || o.x < -4.0)
      return false;
    if(o.z >= -2 || o.z < -10)
      return false;
    return true;
  }
}

/*********************************************************************
 * tranmission vector - helper function
 *********************************************************************/
Vector transmissionVector(Vector i, Vector surf_norm, Spheres *nearest, Point* hit){
  Vector transmission;
  Point temp;
  temp.x = hit->x;
  temp.y = hit->y;
  temp.z = hit->z;
  Vector eyeHit = get_vec(eye_pos, temp);
  Vector eyeCenter = get_vec(eye_pos, nearest->center);
  normalize(&eyeHit);
  float n;
  n = 1.0/(nearest->refraction_index);
  Vector transmissionParallel = vec_scale(vec_minus(i, vec_scale(surf_norm, vec_dot(i, surf_norm))), n);
  if(vec_len(transmissionParallel) > 1.0){
    transmission = {0, 0, 0};
    return transmission;
  }
  Vector transmissionOrthogonal = vec_scale(surf_norm, (-1)*(sqrt(1-pow(vec_dot(transmissionOrthogonal, transmissionParallel), 2))));
  transmission = vec_plus(transmissionParallel, transmissionOrthogonal);
  return transmission;
}


/*********************************************************************
 * Phong illumination - you need to implement this!
 *********************************************************************/
RGB_float phong(Point q, Vector v, Vector surf_norm, Spheres *sph) {
  RGB_float amb, spec, dif, col;
  //initialization
  amb = {0, 0, 0};
  spec = {0, 0, 0};
  dif = {0, 0, 0};
  col = {0, 0, 0};
  float distanceLight; //distance between light source and point
  float diffuseParam; //diffuse parameter
  float shininessParam;//specular shininess parameter
  float reflactanceParam; //reflectance parameter
  float decay; //coefficient for light decay 
  Vector Half;
  Vector LightSource;
  //AMBIENT GLOBAL
  reflactanceParam = sph->reflectance;
  amb.r = amb.r + global_ambient[0] * reflactanceParam;
  amb.g = amb.g + global_ambient[1] * reflactanceParam;
  amb.b = amb.b + global_ambient[2] * reflactanceParam;
  Vector Shadow = get_vec(q, light1);
  if((shadow_on) && (isLightBlocked(q, Shadow, scene))){
    return amb;
  }
  LightSource = get_vec(q, light1);
  distanceLight = vec_len(LightSource);
  normalize(&LightSource);
  decay = (1.0/ ( decay_a + decay_b * distanceLight + decay_c*pow(distanceLight, 2)));
  diffuseParam = vec_dot(surf_norm, LightSource);
  dif.r = dif.r + light1_diffuse[0]*sph->mat_diffuse[0] * diffuseParam;
  dif.g = dif.g + light1_diffuse[1]*sph->mat_diffuse[1] * diffuseParam;
  dif.b = dif.b + light1_diffuse[2]*sph->mat_diffuse[2] * diffuseParam;
  Half = vec_plus(LightSource, vec_scale(surf_norm, -2*vec_dot(LightSource, surf_norm)));
  shininessParam = pow(vec_dot(v, Half), sph->mat_shineness);
  spec.r = spec.r + light1_specular[0]*sph->mat_specular[0] * shininessParam;
  spec.g = spec.g + light1_specular[1]*sph->mat_specular[1] * shininessParam;
  spec.b = spec.b + light1_specular[2]*sph->mat_specular[2] * shininessParam;
  col = clr_add(dif, spec);
  col = clr_scale(col, decay);
  col = clr_add(col, clr_scale(amb, 2));
  return col;
}

/************************************************************************
 * This is the recursive ray tracer - you need to implement this!
 * You should decide what arguments to use.
 ************************************************************************/
RGB_float recursive_ray_trace(Vector ray, Point o, int iteration) {
  RGB_float retColor, reflectionColor, refractionColor;
  retColor = background_clr;
  reflectionColor = {0, 0, 0};
  refractionColor = {0, 0, 0};
  //initialization
  Spheres* nearest;
  Point* hit = new Point;
  nearest = intersect_scene(o, ray, scene, hit);
  if(nearest != NULL){
    Vector view = get_vec(*hit, eye_pos);
    normalize(&view);
    Vector surf_norm = sphere_normal(*hit, nearest);
    normalize(&surf_norm);
    Vector L = get_vec(*hit, o);
    normalize(&L);
    retColor = phong(*hit, view, surf_norm, nearest);
    if((reflection_on) && (iteration < step_max)){
      float ref_scalar = 2.0*(vec_dot(surf_norm, L));
      Vector reflection = vec_plus(vec_scale(surf_norm, ref_scalar), vec_scale(L, -1.0));
      normalize(&reflection);
      iteration+=1;
      reflectionColor = recursive_ray_trace(reflection, *hit, iteration);
      reflectionColor = clr_scale(reflectionColor, nearest->reflectance);
      retColor = clr_add(retColor, reflectionColor);
    }
    if ((refraction_on) && (iteration < step_max))
    {
      Vector transmission = transmissionVector(view, surf_norm, nearest, hit);
      normalize(&transmission);
      iteration+=1;
      transmission.x = hit->x + transmission.x;
      transmission.x = hit->y + transmission.y;
      transmission.x = hit->z + transmission.z;
      refractionColor = recursive_ray_trace(transmission, *hit, iteration);
      refractionColor = clr_scale(refractionColor, nearest->transparency);
      retColor = clr_add(retColor, refractionColor);
    }
    return retColor;
  }
  if(chessboard_on){
    Point* chessboard_hit = new Point;
    bool isBoardHit = intersectChessboard(ray, o, chessboard_hit);
    if(isBoardHit){
      Vector view = get_vec(*chessboard_hit, eye_pos);
      normalize(&view);
      Vector L = get_vec(*chessboard_hit, o);
      normalize(&L);
      retColor = boardColor(*chessboard_hit, view);
      normalize(&chessBoardNormal);
      if((reflection_on) && (iteration < step_max)){
        float ref_scalar = 2.0*(vec_dot(chessBoardNormal, L));
        Vector reflection = vec_plus(vec_scale(chessBoardNormal, ref_scalar), vec_scale(L, -1.0));
        normalize(&reflection);
        iteration += 1;
        reflectionColor = recursive_ray_trace(reflection, *chessboard_hit, iteration);
        reflectionColor = clr_scale(reflectionColor, reflectance);
        retColor = clr_add(retColor, reflectionColor);
      }
      return retColor;
    }
  }
  return retColor;
}

/*********************************************************************
 * This function traverses all the pixels and cast rays. It calls the
 * recursive ray tracer and assign return color to frame
 *
 * You should not need to change it except for the call to the recursive
 * ray tracer. Feel free to change other parts of the function however,
 * if you must.
 *********************************************************************/
void ray_trace() {
  int i, j;
  float x_grid_size = image_width / float(win_width);
  float y_grid_size = image_height / float(win_height);
  float x_start = -0.5 * image_width;
  float y_start = -0.5 * image_height;
  RGB_float ret_color;
  Point cur_pixel_pos;
  Vector ray;

  // ray is cast through center of pixel
  cur_pixel_pos.x = x_start + 0.5 * x_grid_size;
  cur_pixel_pos.y = y_start + 0.5 * y_grid_size;
  cur_pixel_pos.z = image_plane;

  if(supersampling_on){
    Point subPixelPos;
    RGB_float cumulativeColor;
    for (i=0; i<win_height; i++) {
      for (j=0; j<win_width; j++) {
        subPixelPos = cur_pixel_pos;
        ray = get_vec(eye_pos, subPixelPos);
        //
        // You need to change this!!!
        //
        // ret_color = recursive_ray_trace();
        //ret_color = background_clr; // just background for now

        // Parallel rays can be cast instead using below
        //
        // ray.x = ray.y = 0;
        // ray.z = -1.0;
        // ret_color = recursive_ray_trace(cur_pixel_pos, ray, 1);

        // Checkboard for testing
        normalize(&ray);
        ret_color = recursive_ray_trace(ray, eye_pos, 0);
        cumulativeColor = clr_add(ret_color, cumulativeColor);
        subPixelPos.x -= 0.25 * x_grid_size;
        subPixelPos.y -= 0.25 * y_grid_size;
        ray = get_vec(eye_pos, subPixelPos);
        normalize(&ray);
        ret_color = recursive_ray_trace(ray, eye_pos, 0);
        cumulativeColor = clr_add(ret_color, cumulativeColor);
        subPixelPos.x += 0.25 * x_grid_size;
        subPixelPos.y -= 0.25 * y_grid_size;
        ray = get_vec(eye_pos, subPixelPos);
        normalize(&ray);
        ret_color = recursive_ray_trace(ray, eye_pos, 0);
        cumulativeColor = clr_add(ret_color, cumulativeColor);
        subPixelPos.x -= 0.25 * x_grid_size;
        subPixelPos.y += 0.25 * y_grid_size;
        ray = get_vec(eye_pos, subPixelPos);
        normalize(&ray);
        ret_color = recursive_ray_trace(ray, eye_pos, 0);
        cumulativeColor = clr_add(ret_color, cumulativeColor);
        subPixelPos.x += 0.25 * x_grid_size;
        subPixelPos.y += 0.25 * y_grid_size;
        ray = get_vec(eye_pos, subPixelPos);
        normalize(&ray);
        ret_color = recursive_ray_trace(ray, eye_pos, 0);
        cumulativeColor = clr_add(ret_color, cumulativeColor);
        frame[i][j][0] = GLfloat(ret_color.r);
        frame[i][j][1] = GLfloat(ret_color.g);
        frame[i][j][2] = GLfloat(ret_color.b);
        cur_pixel_pos.x += x_grid_size;
      }
      cur_pixel_pos.y += y_grid_size;
      cur_pixel_pos.x = x_start;
    }
  }
  else{
    for (i=0; i<win_height; i++) {
      for (j=0; j<win_width; j++) {
        ray = get_vec(eye_pos, cur_pixel_pos);
        //
        // You need to change this!!!
        //
        // ret_color = recursive_ray_trace();
        //ret_color = background_clr; // just background for now

        // Parallel rays can be cast instead using below
        //
        // ray.x = ray.y = 0;
        // ray.z = -1.0;
        // ret_color = recursive_ray_trace(cur_pixel_pos, ray, 1);

        // Checkboard for testing
        normalize(&ray);
        ret_color = recursive_ray_trace(ray, eye_pos, 0);
        frame[i][j][0] = GLfloat(ret_color.r);
        frame[i][j][1] = GLfloat(ret_color.g);
        frame[i][j][2] = GLfloat(ret_color.b);

        cur_pixel_pos.x += x_grid_size;
      }

      cur_pixel_pos.y += y_grid_size;
      cur_pixel_pos.x = x_start;
    }
  }

  if(stochastic_on){
    for(int i = 0; i < 5; i++){
      int r1, r2, n1, n2;
      r1 = rand() % (int)image_width/2;
      r2 = rand() % (int)image_height/2;
      n1 = rand() % 2;
      n2 = rand() % 2;
      if(n1 == 0) r1 = -r1;
      if(n2 == 0) r2 = -r2;
      cur_pixel_pos.x = (float) r1;
      cur_pixel_pos.y = (float) r2;
      cur_pixel_pos.z = image_plane;
      ray = get_vec(eye_pos, cur_pixel_pos);
      normalize(&ray);
      ret_color = recursive_ray_trace(ray, eye_pos, 0);
      frame[r1 + win_height/2][r2+win_width/2][0] = GLfloat(ret_color.r);
      frame[r1 + win_height/2][r2+win_width/2][1] = GLfloat(ret_color.g);
      frame[r1 + win_height/2][r2+win_width/2][2] = GLfloat(ret_color.b);
    }
  }
}
