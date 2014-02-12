#include <cstdlib>
#include <GL/glut.h>
#include <GL/glu.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <complex>
#include <map>


#define PI 3.14159265  // Should be used from mathlib
inline float sqr(float x) { return x*x; }

using namespace std;

//////////////////////////
//----------------------//
//| :: Vector Utils :: |//
//----------------------//
//////////////////////////

//****************************************************
// Adds vector u and vector v
//***************************************************
vector<double> makeVector(double x, double y, double z){
    vector<double> out(3);
    out[0]=x;
    out[1]=y;
    out[2]=z;
    return out;
}

//****************************************************
// Adds vector u and vector v
//***************************************************
vector<double> operator+(const vector<double>& u, const vector<double>& v){
    vector<double> out(3);
    out[0]=u[0]+v[0];
    out[1]=u[1]+v[1];
    out[2]=u[2]+v[2];
    return out;
}

//****************************************************
// Subtracts vector u and vector v
//***************************************************
vector<double> operator-(const vector<double>& u, const vector<double>& v){
    vector<double> out(3);
    out[0]=u[0]-v[0];
    out[1]=u[1]-v[1];
    out[2]=u[2]-v[2];
    return out;
}

//****************************************************
// Returns the cross product of vector u and vector v
//***************************************************
vector<double> cross(const vector<double>& u, const vector<double>& v){
    vector<double> out(3);
    out[0]=u[1] * v[2] - u[2] * v[1];
    out[1]=u[2] * v[0] - u[0] * v[2];
    out[2]=u[0] * v[1] - u[1] * v[0];
    return out;
}

//****************************************************
// Takes the dot product of vector u and vector v
//***************************************************
double operator*(const vector<double>& u, const vector<double>& v){
    return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
}

//****************************************************
// Scales vector v by scalar s
//***************************************************
vector<double> operator*(const vector<double>& v, double s){
    vector<double> out(3);
    out[0]=s*v[0];
    out[1]=s*v[1];
    out[2]=s*v[2];
    return out;
}

//****************************************************
// Scales vector v by scalar s
//***************************************************
vector<double> operator*(double s, const vector<double>& v){
    return v*s;
}
vector<double> operator/(const vector<double>& v, double s){
    return v*(1.0/s);
}

//****************************************************
// Elementwise multiplication of vectors u and v
//***************************************************
vector<double> mult(vector<double> u, vector<double> v){
    vector<double> out(3);
    out[0]=u[0]*v[0];
    out[1]=u[1]*v[1];
    out[2]=u[2]*v[2];
    return out;
}

//****************************************************
// Takes the squared norm of vector v
//***************************************************
double normSq(const vector<double> v){
    return v*v;
}

//****************************************************
// Takes the norm of vector v
//***************************************************
double norm(const vector<double> v){
    return sqrt(normSq(v));
}

//****************************************************
// Normalizes vector v
//***************************************************
vector<double> normalize(const vector<double> v){
    return v/norm(v);
}

//****************************************************
// Projects vector u into vector v
//***************************************************
vector<double> project(const vector<double>& u, const vector<double>& v){
    return (((u)*(v))/normSq(v))*v;
}

//****************************************************
// Takes the perpendicular component of vector u on
// vector v.
//***************************************************
vector<double> perp(const vector<double>& u, const vector<double>& v){
    return u-project(u,v);
}

//****************************************************
// Takes the perpendicular component of vector u on
// vector v.
//***************************************************
void print(const vector<double>& u){
    cout << u[0] << " ";
    cout << u[1] << " ";
    cout << u[2] << endl;
}




/////////////////////
//-----------------//
//| :: Classes :: |//
//-----------------//
/////////////////////


typedef vector<double>(ColorFn)(double s, double t);
typedef void(renderFn)();

class Viewport; 
class Light;
class PointLight;
class Material;

class Viewport {
    public:
        int w, h; // width and height
};

Viewport viewport;
double sphereRadius;

class Light{
    public:
        Light(double x, double y, double z, vector<double> cl, bool origin_coords){
            vector<double> pos(3);
            pos[0]=x;
            pos[1]=y;
            pos[2]=z;
            this->pos = pos;
            this->cl = cl;
            this->origin_coords = origin_coords;
        }
        virtual vector<double> getPos()=0;
        virtual vector<double> getLightVector(vector<double> position)=0;
        vector<double> getCl(){
            return this->cl;
        }
    protected:
        vector<double> pos;
        vector<double> cl;
        bool origin_coords;
};

class DirectionalLight : public Light{
    public:
        DirectionalLight(double x, double y, double z, vector<double> cl, bool origin_coords) : Light(x,y,z,cl, origin_coords){}
        DirectionalLight(double x, double y, double z, vector<double> cl)  : Light(x,y,z,cl, false){}
        virtual vector<double> getPos(){
            return this->pos;
        }
        virtual vector<double> getLightVector(vector<double> position){
            return -1*normalize(this->getPos());
        }
};

class PointLight : public Light{
    public:
        PointLight(double x, double y, double z, vector<double> cl, bool origin_coords) : Light(x,y,z,cl, origin_coords){}
        PointLight(double x, double y, double z, vector<double> cl) : Light(x,y,z,cl, false){}
        virtual vector<double> getPos(){
            if (!this->origin_coords){
                return this->pos;
            }else{
                return (this->pos*sphereRadius+makeVector(viewport.w/2.0,viewport.h/2.0,0.0));
            }
        }
        virtual vector<double> getLightVector(vector<double> position){
            return normalize(this->getPos()-position);
        }
};

class Material{
    public:
        Material(){
            this->colorFn=NULL;
            this->cr = makeVector(1,1,1);
            this->init(makeVector(1.0,1.0,1.0),16.0);
        }
        Material(ColorFn cr, vector<double> cp, double p){
            this->colorFn=cr;
            this->init(cp,p);
        }
        Material(vector<double> cr, vector<double> cp, double p){
            this->colorFn=NULL;
            this->cr = cr;
            this->init(cp,p);
        }
        vector<double> getCr(double s, double t){
            if (this->colorFn){
                return (*colorFn)(s,t);
            }else{
                return this->cr;
            }
        }
        vector<double> getCp(){
            return this->cp;
        }
        double getP(){
            return this->p;
        }
        double setP(double p){
            this->p=p;
        }
    protected:
        void init(vector<double> cp, double p){
            this->cp = cp;
            this->p = p;
        }
        ColorFn *colorFn;
        vector<double> cr;
        vector<double> cp;
        double p;
};

class Sphere{
    public:
        Sphere(double x, double y, double radius, Material material, vector<Light*> lights){
            this->material=material;
            this->lights=lights;
            this->x=x;
            this->y=y;
            this->radius=radius;
        }
        Material material;
        vector<Light*> lights;
        double x;
        double y;
        double radius;
};

//****************************************************
// Global Variables
//****************************************************


struct Args{
    //Default args
    vector<Light*> lights;
    vector<double> ca;
    vector<double> cp;
    vector<double> cr;
    string texture;
    double p;
    string outFile;

    //Shadow args
    int areaLightWidth;
    int areaLightHeight;
    double areaLightZ;
    double areaLightResolution;
} args;

map<string,Material> textures;



//****************************************************
// Simple init function
//****************************************************
void initScene(){
}


//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
    viewport.w = w;
    viewport.h = h;
    sphereRadius = min(viewport.w,viewport.h)/2.5;

    glViewport (0,0,viewport.w,viewport.h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, viewport.w, 0, viewport.h);
}


/////////////////////////////////
//-----------------------------//
//| :: Rendering Functions :: |//
//-----------------------------//
/////////////////////////////////

double mandelCheck( double x, double y, int n, double threshold){
    complex<double> cOffset(x*2,y*2);
    complex<double> z(0.0,0.0);
    int i;
    for (i = 0; abs(z) < threshold && i < n; i++){
        z=z*z+cOffset;
    }
    return (double)(i)/n;
}
vector<double> mandelTexture(double s, double t){
    static vector<double> c1 = makeVector(1.0,0.0,0.0);
    static vector<double> c2 = makeVector(0.0,0.0,0.0);
    double a = pow(mandelCheck((s-.2)*1.90,t*1.90,18,2.0),1.65);
    return a*c2+(1.0-a)*c1;
}
vector<double> tigerTexture(double s, double t){
    static vector<double> c1 = makeVector(0.0,0.0,0.0);
    static vector<double> c2 = makeVector(1.0,0.2,0.0);
    return (int)round(100*((s-t)*(s+t)))%2==0?c1:c2;
}
vector<double> pokeball(double s, double t){
    static vector<double> red = makeVector(1.0,0.0,0.0);
    static vector<double> white = makeVector(1.0,1.0,1.0);
    static vector<double> black = makeVector(0.0,0.0,0.0);
    if (s*s+t*t < .025*.025){
        return white;
    }else if (s*s+t*t < .03*.03){
        return black;
    }else if (s*s+t*t < .06*.06){
        return white;
    }else if (s*s+t*t < .073*.073){
        return black;
    }else if (abs(t) < .010){
        return black;
    }
    return t>0?red:white;
}
vector<double> greatball(double s, double t){
    static vector<double> green = makeVector(0.0,1.0,1.0);
    static vector<double> red = makeVector(1.0,0.0,0.0);
    static vector<double> white = makeVector(1.0,1.0,1.0);
    static vector<double> black = makeVector(0.0,0.0,0.0);
    if (s*s+t*t < .025*.025){
        return white;
    }else if (s*s+t*t < .03*.03){
        return black;
    }else if (s*s+t*t < .06*.06){
        return white;
    }else if (s*s+t*t < .073*.073){
        return black;
    }else if (abs(t) < .010){
        return black;
    }
    if (t>0 && s*s+t*t<.22*.22){
        return green;
    }
    if (abs(.75*atan(abs(s))-t)<.04 && abs(s)>.15){
        return red;
    }

    return t>0?green:white;
}
vector<double> ultraball(double s, double t){
    static vector<double> yellow = makeVector(1.0,1.0,0.0);
    static vector<double> white = makeVector(1.0,1.0,1.0);
    static vector<double> black = makeVector(0.0,0.0,0.0);
    static vector<double> grey = makeVector(0.1,0.1,0.1);
    if (s*s+t*t < .025*.025){
        return white;
    }else if (s*s+t*t < .03*.03){
        return black;
    }else if (s*s+t*t < .06*.06){
        return white;
    }else if (s*s+t*t < .073*.073){
        return black;
    }else if (abs(t) < .010){
        return black;
    }
    if (1.5*(abs(s)-.5)*(abs(s)-.5)+(t)*(t)<.50*.40 && 1.5*(abs(s)-.5)*(abs(s)-.5)+(t)*(t)>0.30*0.30){
        return t>0?yellow:white;
    }
    return t>0?grey:white;
}
vector<double> masterball(double s, double t){
    static vector<double> purple = makeVector(0.45,0.0,1.0);
    static vector<double> pink = makeVector(1.0,0.0,1.0);
    static vector<double> white = makeVector(1.0,1.0,1.0);
    static vector<double> black = makeVector(0.0,0.0,0.0);
    if (s*s+t*t < .025*.025){
        return white;
    }else if (s*s+t*t < .03*.03){
        return black;
    }else if (s*s+t*t < .06*.06){
        return white;
    }else if (s*s+t*t < .073*.073){
        return black;
    }else if (abs(t) < .010){
        return black;
    }
    if (0.5*(s-.35)*(s-.35)+(t-.25)*(t-.25)<.16*.16 || 0.5*(s+.35)*(s+.35)+(t-.25)*(t-.25)<.16*.16){
        return pink;
    }
    if (abs(s-t+.15)<=.02 && s >= 0 && s < .07){
        return white;
    }
    if (abs(-s-t+.15)<=.02 && s <= 0 && s > -.07){
        return white;
    }
    if ((abs(-s+.09)<=.02 || abs(-s-.09)<=.02) && t>=.1 && t<.2395){
        return white;
    }

    return t>0?purple:white;
}
vector<double> premierball(double s, double t){
    static vector<double> red = makeVector(1.0,0.0,0.0);
    static vector<double> white = makeVector(1.0,1.0,1.0);
    static vector<double> black = makeVector(0.0,0.0,0.0);
    if (s*s+t*t < .025*.025){
        return white;
    }else if (s*s+t*t < .03*.03){
        return black;
    }else if (s*s+t*t < .06*.06){
        return white;
    }else if (s*s+t*t < .073*.073){
        return red;
    }else if (abs(t) < .010){
        return red;
    }
    return white;
}
vector<double> gsball(double s, double t){
    static vector<double> gold = makeVector(1.0,1.0,0.0);
    static vector<double> white = makeVector(1.0,1.0,1.0);
    static vector<double> blue = makeVector(0.0,0.0,0.5);
    if (s*s+t*t < .025*.025){
        return white;
    }else if (s*s+t*t < .03*.03){
        return blue;
    }else if (s*s+t*t < .06*.06){
        return white;
    }else if (s*s+t*t < .073*.073){
        return blue;
    }else if (abs(t) < .010){
        return blue;
    }
    return t>0?gold:white;
}
vector<double> red(double s, double t){return makeVector(1.0,0.0,0.0);}
vector<double> green(double s, double t){return makeVector(0.0,1.0,0.0);}
vector<double> blue(double s, double t){return makeVector(0.0,0.0,1.0);}
vector<double> white(double s, double t){return makeVector(1.0,1.0,1.0);}
vector<double> black(double s, double t){return makeVector(0.0,0.0,0.0);}

//****************************************************
// Renders a sphere
//***************************************************
void drawSpheres(vector<Sphere> spheres, vector<double> ca, bool shadow){
    glBegin(GL_POINTS);

    for (vector<Sphere>::iterator sphere = spheres.begin(); sphere != spheres.end(); sphere++){

        float centerX = sphere->x;
        float centerY = sphere->y;
        float radius = sphere->radius;
        Material m = sphere->material;
        vector<Light*> lights = sphere->lights;

        int minI = max(0,(int)floor(centerX-radius));
        int maxI = min(viewport.w-1,(int)ceil(centerX+radius));

        int minJ = max(0,(int)floor(centerY-radius));
        int maxJ = min(viewport.h-1,(int)ceil(centerY+radius));

        vector<double> n(3), r(3), l(3), pos(3);
        vector<double> e = makeVector(0,0,1);


        for (int i=0;i<viewport.w;i++) {
            for (int j=0;j<viewport.h;j++) {
                double x = (i+0.5-centerX);
                double y = (j+0.5-centerY);
                double distSq = x*x + y*y;
                if (distSq<=radius*radius) {
                    double z = sqrt(radius*radius-distSq);

                    double t = atan(y/sqrt(x*x+z*z))/PI;
                    double s = atan(x/z)/(PI);

                    pos = makeVector(centerX+x,centerY+y,z);
                    n = normalize(makeVector(x,y,z));


                    vector<double> clr = mult(m.getCr(s,t),ca);
                    for (vector<Light*>::iterator light = lights.begin(); light != lights.end(); light++){
                        l = (*light)->getLightVector(pos);
                        r=2*(l*n)*n-l;

                        vector<double> diffuse = mult(m.getCr(s,t),(*light)->getCl()*max(0.0,n*l));
                        vector<double> specular = mult((*light)->getCl(),m.getCp())*pow(max(0.0,e*r),m.getP());

                        if (shadow){
                            double fuzz = 1.0;
                            vector<double> otherCenter;
                            vector<double> displacement;
                            for (vector<Sphere>::iterator otherSphere = spheres.begin(); otherSphere != spheres.end(); otherSphere++){
                                if (otherSphere!=sphere){
                                    otherCenter = makeVector(otherSphere->x,otherSphere->y,0);
                                    displacement = otherCenter-pos;
                                    if (norm(displacement-project(displacement,l)) < otherSphere->radius && displacement*l>0){
                                        fuzz*=(norm(displacement-project(displacement,l)))/otherSphere->radius;
                                    }
                                }
                            }
                            //clr = clr + pow(fuzz,10)*(diffuse+pow(fuzz,5)*specular);
                            if (fuzz == 1.0){
                                clr = clr + diffuse + specular;
                            }
                        }else{
                            clr = clr + diffuse + specular;
                        }
                    }
                    glColor3f(clr[0], clr[1], clr[2]);
                    glVertex2f(i + 0.5, j + 0.5);
                }


            }
        }
    }
    glEnd();
}
void drawSpheres(vector<Sphere> spheres, vector<double> ca){
    drawSpheres(spheres, ca, false);
}



void shadowRender() {
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);	
    glLoadIdentity();

    vector<Light*> lights;
    const double OFFSET = args.areaLightResolution;
    const double Z_OFFSET = args.areaLightZ*min(viewport.w,viewport.h)/4.5;
    const double numA = args.areaLightWidth;
    const double numB = args.areaLightHeight;
    for (int i = -numA; i <= numA; i++){
        for (int j = -numB; j <= numB; j++){
            lights.push_back(new PointLight(2.0*viewport.w,viewport.h/2.0+i*OFFSET,j*OFFSET + Z_OFFSET,makeVector(1.00,1.00,1.00)/((2*numA+1)*(2*numB+1))));
        }
    }

    Material tiger(tigerTexture, makeVector(.7,.7,.7), 16.0);
    Material poke(pokeball, makeVector(.7,.7,.7), 16.0);

    vector<Sphere> spheres;
    spheres.push_back(Sphere(viewport.w*3.0/4.0,viewport.h*1.0/2.0,min(viewport.w,viewport.h)/16.0,tiger,lights));
    spheres.push_back(Sphere(viewport.w*1.0/4.0,viewport.h/2.0,min(viewport.w,viewport.h)/4.5,poke,lights));
    drawSpheres(spheres,makeVector(.1,.1,.1),true);

    glFlush();
    glutSwapBuffers();
}

void demoRender() {
    glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer
    glMatrixMode(GL_MODELVIEW);			        // indicate we are specifying camera transformations
    glLoadIdentity();				        // make sure transformation is "zero'd"

    vector<Light*> lights;

    vector<double> ca = makeVector(.15, .15, .15); //color of ambient light
    vector<double> cp = makeVector(.7, .7, .7); //specular reflectance

    lights.push_back(new PointLight(viewport.w*1.0/2.0,viewport.h*0.0/4.0,min(viewport.w,viewport.h)/4.5+100,makeVector(1.0,1.0,1.0)));
    lights.push_back(new PointLight(viewport.w*1.0/2.0,viewport.h*4.0/4.0,min(viewport.w,viewport.h)/4.5+100,makeVector(1.0,1.0,1.0)));
    lights.push_back(new PointLight(viewport.w*0.0/4.0,viewport.h*1.0/2.0,min(viewport.w,viewport.h)/4.5+100,makeVector(1.0,1.0,1.0)));
    lights.push_back(new PointLight(viewport.w*4.0/4.0,viewport.h*1.0/2.0,min(viewport.w,viewport.h)/4.5+100,makeVector(1.0,1.0,1.0)));

    Material tiger = textures["tiger"];
    Material poke = textures["pokeball"];
    Material mandel = textures["mandelbrot"];
    Material color(makeVector(0,0,1), makeVector(.7,.7,.7), 16.0);

    vector<Sphere> spheres;
    spheres.push_back(Sphere(viewport.w*1.0/4.0 , viewport.h*1.0/4.0, min(viewport.w,viewport.h)/4.5, tiger, lights));
    spheres.push_back(Sphere(viewport.w*3.0/4.0 , viewport.h*1.0/4.0, min(viewport.w,viewport.h)/4.5, poke, lights));
    spheres.push_back(Sphere(viewport.w*1.0/4.0 , viewport.h*3.0/4.0, min(viewport.w,viewport.h)/4.5, mandel, lights));
    spheres.push_back(Sphere(viewport.w*3.0/4.0 , viewport.h*3.0/4.0, min(viewport.w,viewport.h)/4.5, color, lights));
    drawSpheres(spheres,ca,false);

    glFlush();
    glutSwapBuffers();
}

void pokemansRender() {
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);	
    glLoadIdentity();

    vector<Light*> lights;
    vector<Sphere> spheres;
    lights.push_back(new PointLight(viewport.w*1.0/2.0,viewport.h*2.0/2.0,min(viewport.w,viewport.h)/4.5+600,makeVector(1.0,1.0,1.0)));
    spheres.push_back(Sphere(viewport.w*0.8/4.0 , viewport.h*2.0/3.0, min(viewport.w,viewport.h)/8.0, textures["pokeball"], lights));
    spheres.push_back(Sphere(viewport.w*2.0/4.0 , viewport.h*2.0/3.0, min(viewport.w,viewport.h)/8.0, textures["greatball"], lights));
    spheres.push_back(Sphere(viewport.w*3.2/4.0 , viewport.h*2.0/3.0, min(viewport.w,viewport.h)/8.0, textures["ultraball"], lights));
    spheres.push_back(Sphere(viewport.w*0.8/4.0 , viewport.h*1.0/3.0, min(viewport.w,viewport.h)/8.0, textures["premierball"], lights));
    spheres.push_back(Sphere(viewport.w*2.0/4.0 , viewport.h*1.0/3.0, min(viewport.w,viewport.h)/8.0, textures["masterball"], lights));
    spheres.push_back(Sphere(viewport.w*3.2/4.0 , viewport.h*1.0/3.0, min(viewport.w,viewport.h)/8.0, textures["gsball"], lights));
    drawSpheres(spheres,makeVector(.1,.1,.1),false);


    glFlush();
    glutSwapBuffers();
}

void render(){
    glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer
    glMatrixMode(GL_MODELVIEW);			        // indicate we are specifying camera transformations
    glLoadIdentity();				        // make sure transformation is "zero'd"

    vector<Sphere> spheres;
    Material m;
    if (textures.count(args.texture)){
        m=textures[args.texture];
        m.setP(args.p);
    }else{
        m=Material(args.cr,args.cp,args.p);
    }
    spheres.push_back(Sphere(viewport.w/2.0,viewport.h/2.0,sphereRadius,m,args.lights));
    drawSpheres(spheres,args.ca);

    glFlush();
    glutSwapBuffers();
}

void usage(string myName){
    cout << "Usage: " << myName << " [-ka <r> <g> <b>] [-kd <r> <g> <b>] [-ks <r> <g> <b>] [-sp <v>] [-pl <x> <y> <z> <r> <g> <b>] [-dl <x> <y> <z> <r> <g> <b>]" << endl;
    cout << "[-ka <r> <g> <b>] sets the ambient light color to rgb(<r>,<g>,<b>)" << endl;
    cout << "[-kd <r> <g> <b>] sets the sphere's diffuse color to rgb(<r>,<g>,<b>)" << endl;
    cout << "[-ks <r> <g> <b>] sets the sphere's specular color to rgb(<r>,<g>,<b>)" << endl;
    cout << "[-sp <v>] sets the sphere's specular exponent <v>" << endl;
    cout << "[-pl <x> <y> <z> <r> <g> <b>] adds a color rgb(<r>,<g>,<b>) point light at position (<x>,<y>,<z>)" << endl;
    cout << "[-dl <x> <y> <z> <r> <g> <b>] adds a color rgb(<r>,<g>,<b>) direction light pointing in direction (<x>,<y>,<z>)" << endl;
}

//****************************************************
// the usual stuff, nothing exciting here
//****************************************************
int main(int argc, char *argv[]) {
    //Init textures
    textures["tiger"]=Material(tigerTexture, makeVector(.7,.7,.7), 16.0);
    textures["pokeball"]=Material(pokeball, makeVector(.7,.7,.7), 16.0);
    textures["masterball"]=Material(masterball, makeVector(.7,.7,.7), 16.0);
    textures["greatball"]=Material(greatball, makeVector(.7,.7,.7), 16.0);
    textures["ultraball"]=Material(ultraball, makeVector(.7,.7,.7), 16.0);
    textures["premierball"]=Material(premierball, makeVector(.7,.7,.7), 16.0);
    textures["gsball"]=Material(gsball, makeVector(.7,.7,.7), 16.0);
    textures["mandelbrot"]=Material(mandelTexture, makeVector(.7,.7,.7), 16.0);

    //Parse Command Line Args, exiting if error
    renderFn *myDisplay;
    myDisplay = render;

    //Set default args
    args.ca=makeVector(0.0,0.0,0.0);
    args.cr=makeVector(0,0,0);
    args.cp=makeVector(0,0,0);
    args.p=0;
    args.outFile="";

    args.areaLightWidth=1;
    args.areaLightHeight=1;
    args.areaLightZ=-1.3;
    args.areaLightResolution=5;


    for (int i = 1; i < argc; i++){
        if ((string)argv[i]=="-h" || (string)argv[i]=="--help"){ 
            usage(argv[0]);
            return 0;
        }else if((string)argv[i]=="-o"){
            if (i+1 >= argc){
                cout << "Unexpected end of arglist." << endl;
                usage(argv[0]);
                return 1;
            }
            args.outFile = (string)argv[i+1];
            i+=1;
        }else if ((string)argv[i]=="--demo"){ 
            myDisplay = demoRender;
        }else if ((string)argv[i]=="--pokemans"){ 
            myDisplay = pokemansRender;
        }else if ((string)argv[i]=="--shadow"){ 
            myDisplay = shadowRender;
        }else if ((string)argv[i]=="--shadow-areaLight-width"){ 
            if (i+1 >= argc){
                cout << "Unexpected end of arglist." << endl;
                usage(argv[0]);
                return 1;
            }
            args.areaLightWidth=atof(argv[i+1]);
            i+=1;
        }else if ((string)argv[i]=="--shadow-areaLight-height"){ 
            if (i+1 >= argc){
                cout << "Unexpected end of arglist." << endl;
                usage(argv[0]);
                return 1;
            }
            args.areaLightHeight=atof(argv[i+1]);
            i+=1;
        }else if ((string)argv[i]=="--shadow-areaLight-z"){ 
            if (i+1 >= argc){
                cout << "Unexpected end of arglist." << endl;
                usage(argv[0]);
                return 1;
            }
            args.areaLightZ=atof(argv[i+1]);
            i+=1;
        }else if ((string)argv[i]=="--shadow-areaLight-resolution"){ 
            if (i+1 >= argc){
                cout << "Unexpected end of arglist." << endl;
                usage(argv[0]);
                return 1;
            }
            args.areaLightResolution=atof(argv[i+1]);
            i+=1;
        }else if ((string)argv[i]=="--texture" || (string)argv[i]=="-t" ){ 
            if (i+1 >= argc){
                cout << "Unexpected end of arglist." << endl;
                usage(argv[0]);
                return 1;
            }
            args.texture=(string)argv[i+1];
            i+=1;
        }else if ((string)argv[i]=="-ka"){ //Set ambient color
            if (i+3 >= argc){
                cout << "Unexpected end of arglist." << endl;
                usage(argv[0]);
                return 1;
            }
            args.ca=makeVector(atof(argv[i+1]),atof(argv[i+2]),atof(argv[i+3]));
            i+=3;
        }else if ((string)argv[i]=="-kd"){ //Set sphere diffuse color
            if (i+3 >= argc){
                cout << "Unexpected end of arglist." << endl;
                usage(argv[0]);
                return 1;
            }
            args.cr=makeVector(atof(argv[i+1]),atof(argv[i+2]),atof(argv[i+3]));
            i+=3;
        }else if ((string)argv[i]=="-ks"){ //Set sphere specular color
            if (i+3 >= argc){
                cout << "Unexpected end of arglist." << endl;
                usage(argv[0]);
                return 1;
            }
            args.cp=makeVector(atof(argv[i+1]),atof(argv[i+2]),atof(argv[i+3]));
            i+=3;
        }else if ((string)argv[i]=="-sp"){ //Set specular term
            if (i+1 >= argc){
                cout << "Unexpected end of arglist." << endl;
                usage(argv[0]);
                return 1;
            }
            args.p=atof(argv[i+1]);
            i+=1;
        }else if ((string)argv[i]=="-pl"){ //Add point light
            if (i+6 >= argc){
                cout << "Unexpected end of arglist." << endl;
                usage(argv[0]);
                return 1;
            }
            args.lights.push_back(new PointLight(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]), makeVector(atof(argv[i+4]),atof(argv[i+5]),atof(argv[i+6])),true));
            i+=6;
        }else if ((string)argv[i]=="-dl"){ //Add directional light
            if (i+6 >= argc){
                cout << "Unexpected end of arglist." << endl;
                usage(argv[0]);
                return 1;
            }
            args.lights.push_back(new DirectionalLight(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]), makeVector(atof(argv[i+4]),atof(argv[i+5]),atof(argv[i+6])),true));
            i+=6;
        }else{
            cout << "Unexpected " << argv[i] << "." << endl;
            usage(argv[0]);
            return 1;
        }
    }


    //Render
    glutInit(&argc, argv); //This initializes glut
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB); //This tells glut to use a double-buffered window with red, green, and blue channels 
    viewport.w = 400; // Initalize theviewport size
    viewport.h = 400;
    glutInitWindowSize(viewport.w, viewport.h); //The size and position of the window
    glutInitWindowPosition(0,0);
    glutCreateWindow(argv[0]);
    initScene(); // set up scene
    glutDisplayFunc(myDisplay);	// function to run when its time to draw something
    glutReshapeFunc(myReshape); // function to run when the window gets resized
    glutMainLoop(); // infinite loop that will keep drawing and resizing
    return 0;
}
