// Name: Michael Romero
// Quarter, Year: Fall, 2013
// Lab:022
//
// This file is to be modified by the student.
// main.cpp
////////////////////////////////////////////////////////////
#include <GL/glut.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <utility>
using namespace std;

const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;
const double PI = 3.14159265;

int vertices = 0;
int faces = 0;



struct point{
  int x;
  int y;
  int z;
 
};

struct vectr{
  double a1;
  double a2;
  double a3;

  vectr(double posA, double posB, double dir) : a1(posA), a2(posB), a3(dir)
  {}
};


struct connections{
  int v0;
  int v1;
  int v2;
};

struct face{
  point a;
  point b;
  point c;

};


struct Ray
{
  point position;
  point direction;
  
  Ray()
  {
    position.x = 0;
    position.y = 0;
    position.z = -1;

    direction.x = 0;
    direction.y = 0;
    direction.z = 1;
  }
  Ray(int x, int y, int z)
  {
    position.x = x;
    position.y = y;
    position.z = z;

    direction.x = 0;
    direction.y = 0;
    direction.z = 1;
  }

  Ray(int x, int y, int z, int q, int r, int s)
  {
    position.x = x;
    position.y = y;
    position.z = z;

    direction.x = q;
    direction.y = r;
    direction.z = s;
  }

};

vector<point>cpy;
void getVertsandFaces(vector<point> &vs, vector<connections> &ed, int &v, 
		      int &f, vector<face>&fs )
{
  ifstream infile("triangle.txt");
  int vertices = 0;
  infile >> vertices;
  v = vertices;
  int faces = 0;
  infile >> faces;
  f = faces;
  // cout << vertices << " " << faces << endl;
  
  for(int i = 0; i < vertices; i++)
    {
      point vert;
      infile >> vert.x;
      infile >> vert.y;
      infile >> vert.z;
      vert.x*=8;
      vert.y*=8;
      vert.z*=8;
      vs.push_back(vert);
      //cout << vs[i].x << " " << vs[i].y << " " << vs[i].z << endl;
    }
  cpy = vs;
  for(int i = 0; i < faces; i++)
    {
      connections fc;
      infile >> fc.v0;
      infile >> fc.v1;
      infile >> fc.v2;
      ed.push_back(fc);
      point tmp = vs[ed[i].v0];
      point tmp1 = vs[ed[i].v1];
      point tmp2 = vs[ed[i].v2];
      face f;
      f.a = tmp;
      f.b = tmp1;
      f.c = tmp2;
      fs.push_back(f);
      //cout << ed[i].v0 << " " << ed[i].v1 << " " << ed[i].v2 << endl;
    }
   infile.close();
}




// Renders a quad at cell (x, y) with dimensions CELL_LENGTH
void renderPixel(int x, int y, float r = 1.0, float g = 1.0, float b = 1.0)
{
	// ...
	// Complete this function
	// ...
  glBegin(GL_POINTS);
  glColor3f(r,g,b);
  glVertex2i(x,y);
  glEnd();
}

void renderLine(int x2, int y2, int x1, int y1)
{
  
  int dx = x2-x1;
  int dy = y2-y1;

  int steps = 0;

  if(abs(dx) > abs(dy))
    steps = abs(dx);
  else
    steps = abs(dy);

  double x = x1;
  double y = y1;
  
  double x_inc = ((double)(dx))/((double)(steps));
  double y_inc = ((double)(dy))/((double)(steps));
  
  renderPixel((int)(x+0.5),(int)(y+0.5));
  
   for(int i = 0; i < steps; i++)
     {
       x += x_inc;
       y += y_inc;
       renderPixel((int)(x+0.5),(int)(y+0.5));
     }
}
int view = 0;


void glutMouseInput(int button, int state, int x, int y)
{
  switch(button)
    {
    case GLUT_LEFT_BUTTON:
      if(state == GLUT_DOWN)
	{
	  view++;
	  view%=3;
	  //cout << view << endl;
	  glutPostRedisplay();
	}
    }
}



double dotProd(vectr a, vectr b)
{
  return a.a1*b.a1 + a.a2*b.a2 + a.a3*b.a3;
}


pair<vectr, int> crossProd(vectr a, vectr b, point p)
{
  double i = a.a2*b.a3 - b.a2*a.a3;
  double j = -1*(a.a1*b.a3 - b.a1*a.a3);
  double k = a.a1*b.a2 - b.a1*a.a2;

  int res = int((i*p.x + j*p.y + k*p.z)+0.5);
  vectr cp(i,j,k);
  return make_pair(cp,res);
}

vectr normalize(vectr a)
{
  double mag = sqrt(a.a1*a.a1+a.a2*a.a2+a.a3*a.a3);
  vectr norm(a.a1/mag, a.a2/mag, a.a3/mag);
  return norm;
}


int **makeT(point a, point b, point c, point d)
{
  int ** T = new int*[3];
  for(int i = 0; i < 3; i++)
    {
      T[i] = new int[3];
    }
  T[0][0] = -(b.x-a.x);
  T[0][1] = -(c.x-a.x);
  T[0][2] = a.x-d.x;

  T[1][0] = -(b.y-a.y);
  T[1][1] = -(c.y-a.y);
  T[1][2] = a.y-d.y;

  T[2][0] = -(b.z-a.z);
  T[2][1] = -(c.z-a.z);
  T[2][2] = a.z-d.z;

  return T;
}

int **makeB(point a, point b, point c, point d)
{
  int ** T = new int*[3];
  for(int i = 0; i < 3; i++)
    {
      T[i] = new int[3];
    }
  //b is r0
  T[0][0] = a.x-b.x;
  T[0][1] = -(c.x-a.x);
  T[0][2] = d.x;

  T[1][0] = a.y-b.y;
  T[1][1] = -(c.y-a.y);
  T[1][2] = d.y;

  T[2][0] = a.z-b.z;
  T[2][1] = -(c.z-a.z);
  T[2][2] = d.z;

  return T;
}

int **makeC(point a, point b, point c, point d)
{
  int ** T = new int*[3];
  for(int i = 0; i < 3; i++)
    {
      T[i] = new int[3];
    }
  //c is r0
  T[0][0] = -(b.x-a.x);
  T[0][1] = a.x-c.x;
  T[0][2] = d.x;

  T[1][0] = -(b.y-a.y);
  T[1][1] = a.y-c.y;
  T[1][2] = d.y;

  T[2][0] = -(b.z-a.z);
  T[2][1] = a.z-c.z;
  T[2][2] = d.z;

  return T;
}

int **makeA(point a, point b, point c, point d)
{
  int ** T = new int*[3];
  for(int i = 0; i < 3; i++)
    {
      T[i] = new int[3];
    }
  //c is r0
  T[0][0] = -(b.x-a.x);
  T[0][1] = -(c.x-a.x);
  T[0][2] = d.x;

  T[1][0] = -(b.y-a.y);
  T[1][1] = -(c.y-a.y);
  T[1][2] = d.y;

  T[2][0] = -(b.z-a.z);
  T[2][1] = -(c.z-a.z);
  T[2][2] = d.z;

  return T;
}


int det3x3(int ** a)
{
  if(a != NULL)
    {
      int a1 = a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1]);
      int a2 = a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0]);
      int a3 = a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
      int determinant = a1 - a2 + a3;
      return determinant;
    }
  else
      return -1;
}

double getT(double t, double B, double g)
{
  if(t > 0)
    {
      if(B > 0 && g > 0)
	{
	  if(B+g < 1)
	    return t;
	}
      //return 1000000;
    }
  return 1000000;
}

bool intersect(point a, point b, point c, point start, vectr dir)
{
  vectr v1(b.x-a.x,b.y-a.y,b.z-a.z);
  vectr v2(c.x-a.x,c.y-a.y,c.z-a.z);

  pair<vectr,int> n = crossProd(v1,v2,a);
  double notDir = dotProd(n.first,dir);
  if(notDir == 0)
    return false;
  vectr av(a.x,a.y,a.z);
  double d = dotProd(n.first,av);
  vectr startv(start.x,start.y,start.z);
  double t = -(dotProd(n.first,startv)+d);
  if(t < 0)
    return false;
  point p;
  double px = start.x + t*dir.a1;
  double py = start.y + t*dir.a2;
  double pz = start.z + t*dir.a3;
  p.x = int(px+0.5);
  p.y = int(py+0.5);
  p.z = int(pz+0.5);

  vectr e0(b.x-a.x,b.y-a.y,b.z-a.z);
  vectr vp0(p.x-a.x,p.y-a.y,p.z-a.z);

  pair<vectr,int> c0 = crossProd(e0,vp0,a);
  if(dotProd(n.first,c0.first) < 0)
    return false;

  vectr e1(c.x-b.x,c.y-b.y,c.z-b.z);
  vectr vp1(p.x-b.x,p.y-b.y,p.z-b.z);

  c0 = crossProd(e1,vp1,b);
  if(dotProd(n.first,c0.first) < 0)
    return false;

  vectr e2(a.x-c.x,a.y-c.y,a.z-c.z);
  vectr vp2(p.x-c.x,p.y-c.y,p.z-c.z);

  c0 = crossProd(e2,vp2,c);
  if(dotProd(n.first,c0.first) < 0)
    return false;
  return true;
}

double illuminate2(vector<point>&v, vector<connections> e, int index, 
		   int x, int y, point lSource)
{
  double aveX = 0;
  double aveY = 0;
  double aveZ = 0;
  for(int i = 0; i < v.size(); i++)
    {
      aveX += v[i].x;
      aveY += v[i].y;
      aveZ += v[i].z;
    }
  //centroid
  double centrX = aveX/(vertices);
  double centrY = aveY/(vertices);
  double centrZ = aveZ/(vertices);

  double vecABx = v[e[index].v1].x - v[e[index].v0].x;
  double vecABy = v[e[index].v1].y - v[e[index].v0].y;
  double vecABz = v[e[index].v1].z - v[e[index].v0].z;

  double vecACx =  v[e[index].v2].x - v[e[index].v0].x;
  double vecACy =  v[e[index].v2].y - v[e[index].v0].y;
  double vecACz =  v[e[index].v2].z - v[e[index].v0].z;

  vectr ba(vecABx,vecABy,vecABz);
  vectr ca(vecACx,vecACy,vecACz);

  pair<vectr, int> n1 = crossProd(ba,ca,v[e[index].v0]);
  pair<vectr, int> n2 = crossProd(ca,ba,v[e[index].v0]);

  //cout << n1.first.a1 << " " << n1.first.a2 << " " << n1.first.a3 << endl;
  //cout << n2.first.a1 << " " << n2.first.a2 << " " << n2.first.a3 << endl;
  //cout << n1.second << endl;

  double z = (n1.second-(n1.first.a1*x+n1.first.a2*y))/(n1.first.a3);
  point pntOfNrml;
  pntOfNrml.x = x;
  pntOfNrml.y = y;
  pntOfNrml.z = int(z+0.5);

  double vecCentrToNrmlx = pntOfNrml.x - centrX;
  double vecCentrToNrmly = pntOfNrml.y - centrY;
  double vecCentrToNrmlz = pntOfNrml.z - centrZ;
  vectr vctn(vecCentrToNrmlx,vecCentrToNrmly,vecCentrToNrmlz);

  double outward1 = dotProd(n1.first,vctn);
  vectr normalToUse(0,0,0);
  if(outward1 < 0)
    normalToUse = n2.first;
  else
    normalToUse = n1.first;

  //cout << normalToUse.a2 << endl;
  vectr ntu = normalize(normalToUse);
  double vectToLightx = lSource.x - pntOfNrml.x;
  double vectToLighty = lSource.y - pntOfNrml.y;
  double vectToLightz = lSource.z - pntOfNrml.z;
  vectr vtl(vectToLightx,vectToLighty,vectToLightz);
  vectr vtlNorm = normalize(vtl);
  
  //cout << vtlNorm.a1 << " " << vtlNorm.a2 << " " << vtlNorm.a3 << endl;

  vectr toViewer(x-pntOfNrml.x,y-pntOfNrml.y,0-z);
  double hx = vectToLightx;
  double hy = vectToLighty;
  double hz = vectToLightz + (-z);
  double magLv = sqrt(hx*hx+hy*hy+hz*hz);
  vectr h(hx/magLv,hy/magLv,hz/magLv);
  vectr hNorm = normalize(h);
  double ill=0.5*1*dotProd(vtlNorm,ntu)+0.7*1*pow(dotProd(ntu,hNorm),5)+0.6*1;
  //cout << ill << endl;

  bool doesIntersect = false;
  for(int i = 0; i < e.size(); i++)
    {
      if(i != index)
	{
	  doesIntersect = intersect(v[e[i].v0],v[e[i].v1],v[e[i].v2],
				   pntOfNrml,vtl);
	  if(doesIntersect)
	    {
	      ill/=2;
	      break;
	    }
	}
    }

  return ill;
}


void rayTraceUsingEdges(vector<point>&vs, vector<connections> ed)
{

  int ** T = NULL;
  int ** B = NULL;
  int ** C = NULL;
  int ** A = NULL;

  double smallest = 1000000;
  int indexOfFront = 0;
  point light;
  light.x = 400;
  light.y = 800;
  light.z = 0;
  for(int q  = 0; q < WINDOW_HEIGHT; q++)
    {
      for(int p = 0; p < WINDOW_WIDTH; p++)
	{
	  Ray ry(q,p,-1);
	  for(int i = 0; i < ed.size(); i++)
	    {
	      point a = vs[ed[i].v0];
	      point b = vs[ed[i].v1];
	      point c = vs[ed[i].v2];
	 
	      //cout << ed[i].v0 << " " <<ed[i].v1 << " " << ed[i].v2 << endl;
	      T = makeT(a,b,c,ry.position);
	      B = makeB(a,ry.position,c,ry.direction);
	      C = makeC(a,b,ry.position,ry.direction);
	      A = makeA(a,b,c,ry.direction);
	      
	      
	      int detT = det3x3(T);
	      int detB = det3x3(B);
	      int detC = det3x3(C);
	      int detA = det3x3(A);
	      double ta = double(detT)/double(detA);
	      double ba = double(detB)/double(detA);
	      double ca = double(detC)/double(detA);
	      
	      double sm = getT(ta,ba,ca);
	      if(sm < smallest)
		{
		  smallest = sm;
		  indexOfFront = i;
		}
	      
	      for(int i = 0; i < 3; i++)
		{
		  delete []T[i];
		  delete []B[i];
		  delete []C[i];
		  delete []A[i];		  
		}
	      delete []T;
	      delete []B;
	      delete []C;
	      delete []A;
	    }

	  double ill = illuminate2(vs,ed,indexOfFront,q,p,light);
	  if(smallest < 1000000)
	    renderPixel(q,p,1.0*ill,1.0*ill,1.0*ill);
	  
	  smallest = 1000000;
	}
    }
}


double illuminate(vector<face>&fs, int index, int x, int y, point lSource)
{
  double aveX = 0;
  double aveY = 0;
  double aveZ = 0;
  for(int i = 0; i < fs.size(); i++)
    {
      aveX += (fs[i].a.x+fs[i].b.x+fs[i].c.x);
      aveY += (fs[i].a.y+fs[i].b.y+fs[i].c.y);
      aveZ += (fs[i].a.z+fs[i].b.z+fs[i].c.z);
    }
  //centroid
  double centrX = aveX/(3.0*fs.size());
  double centrY = aveY/(3.0*fs.size());
  double centrZ = aveZ/(3.0*fs.size());

  double vecABx = fs[index].b.x - fs[index].a.x;
  double vecABy = fs[index].b.y - fs[index].a.y;
  double vecABz = fs[index].b.z - fs[index].a.z;

  double vecACx = fs[index].c.x - fs[index].a.x;
  double vecACy = fs[index].c.y - fs[index].a.y;
  double vecACz = fs[index].c.z - fs[index].a.z;

  vectr ba(vecABx,vecABy,vecABz);
  vectr ca(vecACx,vecACy,vecACz);

  pair<vectr, int> n1 = crossProd(ba,ca,fs[index].a);
  pair<vectr, int> n2 = crossProd(ca,ba,fs[index].a);

  //cout << n1.first.a1 << " " << n1.first.a2 << " " << n1.first.a3 << endl;
  //cout << n2.first.a1 << " " << n2.first.a2 << " " << n2.first.a3 << endl;
  //cout << n1.second << endl;

  double z = (n1.second-(n1.first.a1*x+n1.first.a2*y))/(n1.first.a3);
  point pntOfNrml;
  pntOfNrml.x = x;
  pntOfNrml.y = y;
  pntOfNrml.z = int(z+0.5);

  double vecCentrToNrmlx = pntOfNrml.x - centrX;
  double vecCentrToNrmly = pntOfNrml.y - centrY;
  double vecCentrToNrmlz = pntOfNrml.z - centrZ;
  vectr vctn(vecCentrToNrmlx,vecCentrToNrmly,vecCentrToNrmlz);

  double outward1 = dotProd(n1.first,vctn);
  vectr normalToUse(0,0,0);
  if(outward1 < 0)
    normalToUse = n2.first;
  else
    normalToUse = n1.first;

  //cout << normalToUse.a2 << endl;
  vectr ntu = normalize(normalToUse);
  double vectToLightx = lSource.x - pntOfNrml.x;
  double vectToLighty = lSource.y - pntOfNrml.y;
  double vectToLightz = lSource.z - pntOfNrml.z;
  vectr vtl(vectToLightx,vectToLighty,vectToLightz);
  vectr vtlNorm = normalize(vtl);
  
  //cout << vtlNorm.a1 << " " << vtlNorm.a2 << " " << vtlNorm.a3 << endl;

  double ill = 0.5*1.0 + 0.5*1*dotProd(vtlNorm,ntu);
  //cout << ill << endl;
  return ill;
}


void rayTraceUsingFaces(vector<face>&fs)
{
  
  int ** T = NULL;
  int ** B = NULL;
  int ** C = NULL;
  int ** A = NULL;

  double smallest = 1000000;
  int indexOfFront = 0;
  point light;
  light.x = 400;
  light.y = 800;
  light.z = 600;
  for(int q  = 0; q < WINDOW_HEIGHT; q++)
    {
      for(int p = 0; p < WINDOW_WIDTH; p++)
	{
	  Ray ry(q,p,-1);
	  for(int i = 0; i < fs.size(); i++)
	    {
	      point a = fs[i].a;
	      point b = fs[i].b;
	      point c = fs[i].c;
	 
	      //cout << ed[i].v0 << " " <<ed[i].v1 << " " << ed[i].v2 << endl;
	      T = makeT(a,b,c,ry.position);
	      B = makeB(a,ry.position,c,ry.direction);
	      C = makeC(a,b,ry.position,ry.direction);
	      A = makeA(a,b,c,ry.direction);
	      
	      
	      int detT = det3x3(T);
	      int detB = det3x3(B);
	      int detC = det3x3(C);
	      int detA = det3x3(A);
	      double ta = double(detT)/double(detA);
	      double ba = double(detB)/double(detA);
	      double ca = double(detC)/double(detA);
	      
	      double sm = getT(ta,ba,ca);
	      if(sm < smallest)
		{
		  smallest = sm;
		  indexOfFront = i;
		}
	      
	      for(int i = 0; i < 3; i++)
		{
		  delete []T[i];
		  delete []B[i];
		  delete []C[i];
		  delete []A[i];		  
		}
	      delete []T;
	      delete []B;
	      delete []C;
	      delete []A;
	    }

	  double ill = illuminate(fs,indexOfFront,q,p,light);
	  if(smallest < 1000000)
	    renderPixel(q,p,1.0*ill,1.0*ill,1.0*ill);
	  
	  smallest = 1000000;
	}
    }
}

void renderWire(vector<face>fs)
{
  for(int i = 0; i < fs.size(); i++)
    {
      renderLine(fs[i].a.x,fs[i].a.y,fs[i].b.x,fs[i].b.y);
      renderLine(fs[i].a.x,fs[i].a.y,fs[i].c.x,fs[i].c.y);
      renderLine(fs[i].b.x,fs[i].b.y,fs[i].c.x,fs[i].c.y);
    }
}

void renderWire2(vector<point>vs, vector<connections>ed)
{
  for(int i = 0; i < ed.size(); i++)
    {
      renderLine(vs[ed[i].v0].x,vs[ed[i].v0].y,vs[ed[i].v1].x,vs[ed[i].v1].y);
      renderLine(vs[ed[i].v0].x,vs[ed[i].v0].y,vs[ed[i].v2].x,vs[ed[i].v2].y);
      renderLine(vs[ed[i].v1].x,vs[ed[i].v1].y,vs[ed[i].v2].x,vs[ed[i].v2].y);
    }
}


//globals
vector<face> fs;
int state = 1000;
vector<point> vs;
vector<connections> ed;

int degreeX = 0;
int degreeY = 0;

void transLate(vector<point>&vs,vector<point>&cp ,int hw)
{
  for(int i = 0; i < vertices; i++)
    {
      if(hw == 0)
	{
	  vs[i].y++;
	  cp[i].y++;
	}
      else if(hw == 1)
	{
	  vs[i].x++;
	  cp[i].x++;
	}
      else if(hw == 2)
	{
	  vs[i].y--;
	  cp[i].y--;
	}
      else if(hw == 3)
	{
	  vs[i].x--;
	  cp[i].x--;
	}
    }
}

void roTate(vector<point>&vs,vector<point>cp, int degX, int degY,int hw)
{
  
  double aveX = 0;
  double aveY = 0;
  double aveZ = 0;
  for(int i = 0; i < vertices; i++)
    {
      aveX += cp[i].x;
      aveY += cp[i].y;
      aveZ += cp[i].z;
    }
  //centroid
  double centrX = aveX/(vertices);
  double centrY = aveY/(vertices);
  double centrZ = aveZ/(vertices);
  //cout << centrX << " " << centrY << " " << centrZ << endl;

  for(int i = 0; i < vertices; i++)
    {
      cp[i].x -= int(centrX+0.5);
      cp[i].y -= int(centrY+0.5);
      cp[i].z -= int(centrZ+0.5);
    }
  if(hw == 0)
    {
      for(int i = 0; i < vertices; i++)
	{ 
	  double cosTh = cos(degY*PI/180.0);
	  double sinTh = sin(degY*PI/180.0);
	  double zp = cp[i].y*sinTh + cp[i].z*cosTh;
	  double yp = cp[i].y*cosTh - cp[i].z*sinTh;
	  //cout << yp << endl;
	  cp[i].z = int(zp+0.5);
	  cp[i].y = int(yp+0.5);
	  // cout << cp[i].z << " " << cp[i].y << endl;
	}
      for(int i = 0; i < vertices; i++)
	{
      
	  cp[i].x += int(centrX+0.5);
	  cp[i].y += int(centrY+0.5);
	  cp[i].z += int(centrZ+0.5);
	}
      for(int i = 0; i < vertices; i++)
	{
	  vs[i].z = cp[i].z;
	  vs[i].y = cp[i].y;
	}
    }
  else if(hw == 1)
    {
      for(int i = 0; i < vertices; i++)
	{ 
	  double cosTh = cos(degX*PI/180.0);
	  double sinTh = sin(degX*PI/180.0);
	  double xp = cp[i].z*sinTh + cp[i].x*cosTh;
	  double zp = cp[i].z*cosTh - cp[i].x*sinTh;
	  //cout << yp << endl;
	  cp[i].z = int(zp+0.5);
	  cp[i].x = int(xp+0.5);
	  //cout << cp[i].x << " " << cp[i].z << endl;
	}
      for(int i = 0; i < vertices; i++)
	{
      
	  cp[i].x += int(centrX+0.5);
	  cp[i].y += int(centrY+0.5);
	  cp[i].z += int(centrZ+0.5);
	}
      for(int i = 0; i < vertices; i++)
	{
	  vs[i].x = cp[i].x;
	  vs[i].z = cp[i].z;
	}
    }
  else if(hw == 2)
    {
      for(int i = 0; i < vertices; i++)
	{ 
	  double cosTh = cos(degY*PI/180.0);
	  double sinTh = sin(degY*PI/180.0);
	  double zp = cp[i].y*sinTh + cp[i].z*cosTh;
	  double yp = cp[i].y*cosTh - cp[i].z*sinTh;
	  //cout << yp << endl;
	  cp[i].z = int(zp+0.5);
	  cp[i].y = int(yp+0.5);
	  //cout << cp[i].z << " " << cp[i].y << endl;
	}
      for(int i = 0; i < vertices; i++)
	{
      
	  cp[i].x += int(centrX+0.5);
	  cp[i].y += int(centrY+0.5);
	  cp[i].z += int(centrZ+0.5);
	}
      for(int i = 0; i < vertices; i++)
	{
	  vs[i].z = cp[i].z;
	  vs[i].y = cp[i].y;
	}
    }
   else if(hw == 3)
    {
      for(int i = 0; i < vertices; i++)
	{ 
	  double cosTh = cos(degX*PI/180.0);
	  double sinTh = sin(degX*PI/180.0);
	  double xp = cp[i].z*sinTh + cp[i].x*cosTh;
	  double zp = cp[i].z*cosTh - cp[i].x*sinTh;
	  //cout << yp << endl;
	  cp[i].z = int(zp+0.5);
	  cp[i].x = int(xp+0.5);
	  //cout << cp[i].x << " " << cp[i].z << endl;
	}
      for(int i = 0; i < vertices; i++)
	{
      
	  cp[i].x += int(centrX+0.5);
	  cp[i].y += int(centrY+0.5);
	  cp[i].z += int(centrZ+0.5);
	}
      for(int i = 0; i < vertices; i++)
	{
	  vs[i].x = cp[i].x;
	  vs[i].z = cp[i].z;
	}
    }
  
 
}


void keyPressed(unsigned char key, int x, int y)
{
  if(key == 't')
    state = 0;
  else if(key == 'r')
    state = 1;
  else if(key == 'i')
    state = 2;
  

  if(key == 'w' && state == 0)
    transLate(vs,cpy,0);
  else if(key == 's' && state == 0)
    transLate(vs,cpy,2);
  else if(key == 'd' && state == 0)
    transLate(vs,cpy,1);
  else if(key == 'a' && state == 0)
    transLate(vs,cpy,3);

  if(key == 'w' && state == 1)
    {
      roTate(vs,cpy,degreeX,degreeY,0);
      degreeY++;
      degreeY%=360;
    }
  else if(key == 's' && state == 1)
    {
      roTate(vs,cpy,degreeX,degreeY,2);
      degreeY--;
      degreeY%=360;
    }
  else if(key == 'd' && state == 1)
    {
      roTate(vs,cpy,degreeX,degreeY,1);
      degreeX++;
      degreeX%=360;
    }
  else if(key == 'a' && state == 1)
    {
      roTate(vs,cpy,degreeX,degreeY,3);
      degreeX--;
      degreeX%=360;
    }
  glutPostRedisplay();
}


//Output function to OpenGL Buffer
void GL_render()
{
   glClear(GL_COLOR_BUFFER_BIT);
   
   // renderLine(700,100,100,100);
   //renderLine(100,700,100,100);

   
   //cout << fs.size() << endl;
   if(state != 2)
     {
       renderWire2(vs,ed);
       //renderWire(fs);
     }
   else if(state == 2)
     rayTraceUsingEdges(vs,ed);
     //rayTraceUsingFaces(fs);
   
   //fs.clear();
   glutSwapBuffers();
}

//Initializes OpenGL attributes
void GLInit(int* argc, char** argv)
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);

	// ...
	// Complete this function
	// ...
	glutCreateWindow("CS 130 - <Insert Name Here>");

	// The default view coordinates is (-1.0, -1.0) bottom left & (1.0, 1.0) top right.
	// For the purposes of this lab, this is set to the number of pixels
	// in each dimension.
	glMatrixMode(GL_PROJECTION_MATRIX);
	glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, -1, 1);
	glutDisplayFunc(GL_render);       
}



int main(int argc, char** argv)
{	
   getVertsandFaces(vs,ed,vertices,faces,fs);
   vectr a(3,4,5);
   vectr b(1,7,6);
   point p0;
   p0.x = 1;
   p0.y = 1;
   p0.z = 5;
   cout << dotProd(a,b) << endl;
   cout << crossProd(a,b,p0).first.a1 << " " << crossProd(a,b,p0).first.a2
	<< " " << crossProd(a,b,p0).first.a3 <<endl;
   cout << crossProd(a,b,p0).second << endl;
   cout << normalize(a).a1 << " " << normalize(a).a2 << " " << normalize(a).a3
	<< endl;

   illuminate(fs,0,160,160,p0);
   point background0;
   point background1;
   point background2;
   point background3;
   point background4;
   background0.x = 0;
   background0.y = 0;
   background0.z = 0;

   background1.x = 100;
   background1.y = 200;
   background1.z = 800;

   background2.x = 700;
   background2.y = 200;
   background2.z = 800;

   background3.x = 800;
   background3.y = 0;
   background3.z = 0;

   background4.x = 400;
   background4.y = -600;
   background4.z = 400;

   vs.push_back(background0);
   vs.push_back(background1);
   vs.push_back(background2);
   vs.push_back(background3);
   vs.push_back(background4);
   connections s1;
   connections s2;
   connections s3;
   connections s4;
   connections s5;
   connections s6;
   s1.v0 = vertices;
   s1.v1 = vertices+1;
   s1.v2 = vertices+2;
   
   s2.v0 = vertices;
   s2.v1 = vertices+2;
   s2.v2 = vertices+3;

   s3.v0 = vertices+4;
   s3.v1 = vertices;
   s3.v2 = vertices+1;

   s4.v0 = vertices+4;
   s4.v1 = vertices+1;
   s4.v2 = vertices+2;

   s5.v0 = vertices+4;
   s5.v1 = vertices+2;
   s5.v2 = vertices+3;

   s6.v0 = vertices+4;
   s6.v1 = vertices;
   s6.v2 = vertices+3;
   
   ed.push_back(s1);
   ed.push_back(s2);
   ed.push_back(s3);
   ed.push_back(s4);
   ed.push_back(s5);
   ed.push_back(s6);
  
   GLInit(&argc, argv);
   glutKeyboardFunc(keyPressed);
   glutMainLoop();
	

	return 0;
}
