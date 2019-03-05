#include "sphere.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <cfloat>

/**********************************************************************
 * This function intersects a ray with a given sphere 'sph'. You should
 * use the parametric representation of a line and do the intersection.
 * The function should return the parameter value for the intersection, 
 * which will be compared with others to determine which intersection
 * is closest. The value -1.0 is returned if there is no intersection
 *
 * If there is an intersection, the point of intersection should be
 * stored in the "hit" variable
 **********************************************************************/

float intersect_sphere(Point o, Vector u, Spheres *sph, Point *hit) {
  Vector projection;
  float distance;
  Point center = sph->center;
  float rad = sph->radius;
  Vector eyeCenter = get_vec(o, center);
  float rayCenter = vec_dot(u, eyeCenter);
    if(rayCenter <= 0)
      return -1.0;
  projection = vec_scale(u, rayCenter);
  distance = vec_len(vec_minus(projection, eyeCenter));
    if(distance > rad)
      return -1.0;
    else if(distance == rad){
      hit->x = projection.x;
      hit->y = projection.y;
      hit->z = projection.z;
      return rayCenter;
    }
    else{
      float param = sqrt(rad*rad - distance * distance);
      param = vec_len(projection) - param;
      Point inter_1 = get_point(o, vec_scale(u, param));
      hit->x = inter_1.x;
      hit->y = inter_1.y;
      hit->z = inter_1.z;
      return param; 
    }
}

bool isLightBlocked(Point o, Vector u, Spheres* sph)
{
  while(sph != NULL)
  {
      Point center = sph->center;
      float rad = sph->radius;
      float a = vec_dot(u,u); 
      Vector centerOrigin = get_vec(center, o);
      float b = vec_dot(centerOrigin, vec_scale(u, 2));
      Vector tempCenter = {center.x, center.y, center.z};
      Vector temp = {o.x, o.y, o.z};
      float c = vec_dot(tempCenter, tempCenter) + vec_dot(temp, temp) - vec_dot(tempCenter, temp)*2 - pow(rad, 2);
      float root = pow(b,2) - 4 * a * c;
      float squareRoot = sqrt(root);
      float r1 = (-b - squareRoot)/2;
      float r2 = (-b + squareRoot)/2;
      if(root > 0){
        if(r1 > 0.0001 || r2 > 0.0001)
        return true;
      }
    sph = sph->next;
  }
  return false;
}

/*********************************************************************
 * This function returns a pointer to the sphere object that the
 * ray intersects first; NULL if no intersection. You should decide
 * which arguments to use for the function. For exmaple, note that you
 * should return the point of intersection to the calling function.
 **********************************************************************/
Spheres *intersect_scene(Point o, Vector u, Spheres *sph, Point *hit) {
  Spheres* sphTemp = sph;
  Spheres* nearest = NULL;
  float nearestDistance, currentDistance;
  nearestDistance = FLT_MAX;
  while(sphTemp != NULL) {
    currentDistance = intersect_sphere(o, u, sphTemp, hit);
    if((nearestDistance > currentDistance) && (currentDistance != -1.0)) {
      nearestDistance = currentDistance;
      nearest = sphTemp;
    }
    sphTemp = sphTemp->next;
  }
  return nearest;
}

/*****************************************************
 * This function adds a sphere into the sphere list
 *
 * You need not change this.
 *****************************************************/
Spheres *add_sphere(Spheres *slist, Point ctr, float rad, float amb[],
		    float dif[], float spe[], float shine, 
		    float refl, float refract, float alpha, int sindex) {
  Spheres *new_sphere;

  new_sphere = (Spheres *)malloc(sizeof(Spheres));
  new_sphere->index = sindex;
  new_sphere->center = ctr;
  new_sphere->radius = rad;
  (new_sphere->mat_ambient)[0] = amb[0];
  (new_sphere->mat_ambient)[1] = amb[1];
  (new_sphere->mat_ambient)[2] = amb[2];
  (new_sphere->mat_diffuse)[0] = dif[0];
  (new_sphere->mat_diffuse)[1] = dif[1];
  (new_sphere->mat_diffuse)[2] = dif[2];
  (new_sphere->mat_specular)[0] = spe[0];
  (new_sphere->mat_specular)[1] = spe[1];
  (new_sphere->mat_specular)[2] = spe[2];
  new_sphere->mat_shineness = shine;
  new_sphere->reflectance = refl;
  new_sphere->refraction_index = refract;
  new_sphere->transparency = alpha;
  new_sphere->next = NULL;

  if (slist == NULL) { // first object
    slist = new_sphere;
  } else { // insert at the beginning
    new_sphere->next = slist;
    slist = new_sphere;
  }

  return slist;
}

/******************************************
 * computes a sphere normal - done for you
 ******************************************/
Vector sphere_normal(Point q, Spheres *sph) {
  Vector rc;

  rc = get_vec(sph->center, q);
  normalize(&rc);
  return rc;
}

