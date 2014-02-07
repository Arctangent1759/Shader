#include <cstdlib>
#include <GL/glut.h>
#include <GL/glu.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <math.h>


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
            return normalize(this->getPos());
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
        Material(double r, double g, double b, vector<double> cp, double p){
            this->colorFn=NULL;
            this->cr = makeVector(r,g,b);
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

//****************************************************
// Global Variables
//****************************************************


struct Args{
    vector<Light*> lights;
    vector<double> ca;
    vector<double> cp;
    vector<double> cr;
    double p;
    string outFile;
} args;



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
    sphereRadius = min(viewport.w,viewport.h)/4.5;

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

#include <complex>
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
    }else if (s*s+t*t < .07*.07){
        return black;
    }else if (abs(t) < .008){
        return black;
    }
    return t>0?red:white;
}
vector<double> red(double s, double t){return makeVector(1.0,0.0,0.0);}
vector<double> green(double s, double t){return makeVector(0.0,1.0,0.0);}
vector<double> blue(double s, double t){return makeVector(0.0,0.0,1.0);}
vector<double> white(double s, double t){return makeVector(1.0,1.0,1.0);}
vector<double> black(double s, double t){return makeVector(0.0,0.0,0.0);}

//****************************************************
// Renders a sphere
//***************************************************
void sphere(float centerX, float centerY, float radius, Material m, vector<double> ca, vector<Light*> lights){
    glBegin(GL_POINTS);

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
                    clr = clr + mult(m.getCr(s,t),(*light)->getCl()*max(0.0,n*l))+mult((*light)->getCl(),m.getCp())*pow(max(0.0,e*r),m.getP());
                }
                glColor3f(clr[0], clr[1], clr[2]);
                glVertex2f(i + 0.5, j + 0.5);
            }


        }
    }
    glEnd();
}


//****************************************************
// function that does the actual drawing of stuff
//***************************************************
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

    Material tiger(tigerTexture, makeVector(.7,.7,.7), 16.0);
    Material poke(pokeball, makeVector(.7,.7,.7), 16.0);
    Material mandel(mandelTexture, makeVector(.7,.7,.7), 16.0);
    Material color(0,0,1, makeVector(.7,.7,.7), 16.0);

    sphere(viewport.w*1.0/4.0 , viewport.h*1.0/4.0, min(viewport.w,viewport.h)/4.5, tiger, ca, lights);
    sphere(viewport.w*3.0/4.0 , viewport.h*1.0/4.0, min(viewport.w,viewport.h)/4.5, poke, ca, lights);
    sphere(viewport.w*1.0/4.0 , viewport.h*3.0/4.0, min(viewport.w,viewport.h)/4.5, mandel, ca, lights);
    sphere(viewport.w*3.0/4.0 , viewport.h*3.0/4.0, min(viewport.w,viewport.h)/4.5, color, ca, lights);

    glFlush();
    glutSwapBuffers();
}
void render(){
    glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer
    glMatrixMode(GL_MODELVIEW);			        // indicate we are specifying camera transformations
    glLoadIdentity();				        // make sure transformation is "zero'd"

    Material m(args.cr,args.cp,args.p);

    sphere(viewport.w/2.0,viewport.h/2.0,sphereRadius,m,args.ca,args.lights);

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
    //Parse Command Line Args, exiting if error
    renderFn *myDisplay;
    myDisplay = render;

    //Set default args
    args.ca=makeVector(0.1,0.1,0.1);
    args.cr=makeVector(0,0,0);
    args.cp=makeVector(1,1,1);
    args.p=16;
    args.outFile="";


    for (int i = 1; i < argc; i++){
        if ((string)argv[i]=="-h" || (string)argv[i]=="--help"){ 
            usage(argv[0]);
            return 0;
        }else if((string)argv[i]=="-o"){
            args.outFile = (string)argv[i+1];
            i+=1;
        }else if ((string)argv[i]=="--demo"){ 
            myDisplay = demoRender;
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
    if (args.outFile==""){
        glutCreateWindow(argv[0]);
    }
    initScene(); // set up scene

    glutDisplayFunc(myDisplay);	// function to run when its time to draw something

    glutReshapeFunc(myReshape); // function to run when the window gets resized
    if (args.outFile==""){
        glutMainLoop(); // infinite loop that will keep drawing and resizing
    }

    GLubyte *data = (GLubyte*) malloc(3 * viewport.w * viewport.h);
    glReadPixels(0, 0, viewport.w, viewport.h, GL_RGB, GL_UNSIGNED_BYTE, data);

    return 0;
}
