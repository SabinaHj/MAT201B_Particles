// Sabina Hyoju Ahn
// particles-p2
// MAT201B_2022_Assignment2-2

#include "al/app/al_App.hpp"
#include "al/app/al_GUIDomain.hpp"
#include "al/math/al_Random.hpp"

using namespace al;

#include <fstream>
#include <vector>
using namespace std;



template <typename T> T sqrt(T &&that) {
  return T(that).normalize(sqrt(that.mag()));
}

Vec3f randomVec3f(float scale) {
  return Vec3f(rnd::uniformS(), rnd::uniformS(), rnd::uniformS()) * scale; //uniform -1 ~ 1
}


string slurp(string fileName); // forward declaration //slurp: file name -> string for shader

struct AlloApp : App {
    
  Parameter pointSize{"/pointSize", "", 1.0, "", 0.0, 5.0};
  Parameter timeStep{"/timeStep", "", 1, "", 0.01, 6};
//   Parameter gravityRandomConst{"/gravityRandomConst", "", 0.01, "", 0.005, 0.2};
  // Parameter drag{"/drag", "", 0.3, "", 0.01, 0.6};
  // Parameter maxAcc{"/maxAcc", "", 0.01, "", 0.001, 0.1};

  ShaderProgram pointShader;

  //  simulation state
  Mesh mesh; // position *is inside the mesh* mesh.vertices() are the positions
  vector<Vec3f> velocity;
  vector<Vec3f> acceleration;
  vector<float> mass;
//   float gravityConst[1000];

  float force;
  // float gravityConst; 
  
  
  // Vec3f initpos;

  Vec3f unitVector;
  Vec3f dist;
  Vec3f gravity;

  bool warp_size = true;
  bool warp_distance = true;

  struct Body {
  string name;

  double mass;     // kilograms
  double radius;   // meters
  double distance; // meters
  double velocity; // meters/second
};

std::vector<Body> body = {
    {    "Sun",1989100000e21, 695508e3,        0,      0},
    {"Jupiter",   1898187e21,  69911e3,  778.3e9, 13.1e3},
    { "Saturn",    568317e21,  58232e3,   1427e9,  9.7e3},
    { "Uranus",     86813e21,  25362e3,   2871e9,  6.8e3},
    {"Neptune",    102413e21,  24622e3, 4497.1e9,  5.4e3},
    {  "Earth",    5972.4e21,   6371e3,  149.6e9, 29.8e3},
    {  "Venus",    4867.5e21,   6052e3,  108.2e9,   35e3},
    {   "Mars",    641.71e21, 3389.5e3,  227.9e9, 24.1e3},
    {"Mercury",    330.11e21, 2439.4e3,   57.9e9, 47.9e3},
    // Pluto? Moons? 'Oumuamua?
};

const double G = 6.674e-11;

  void onInit() override {
    // set up GUI
    auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
    auto &gui = GUIdomain->newGUI();
    gui.add(pointSize); // add parameter to GUI
    gui.add(timeStep);  // add parameter to GUI
    // gui.add(gravityRandomConst);
    // gui.add(drag);
    // gui.add(maxAcc);
  }

  void onCreate() override {
    navControl().useMouse(false);

    // compile shaders
    pointShader.compile(slurp("../point-vertex.glsl"),
                        slurp("../point-fragment.glsl"),
                        slurp("../point-geometry.glsl"));
 
 
//  if (slurp(../point-vertex.glsl"")==""){
//     cout << "ERROR:shader failed absent" <<endl;
//     exit(1);
//  }
// if (!pointShader.compile(slurp("../point-vertex.glsl"),
//                         slurp("../point-fragment.glsl"),
//                         slurp("../point-geometry.glsl"))){
//                           cout << "ERROR:shader failed to complie(absent)" <<endl;
//                           exit(1);
//                         }
    // set initial conditions of the simulation
    

    // c++11 "lambda" function
    auto randomColor = []() { 
      return HSV(rnd::uniform(), 1.0f, 1.0f); 
      };

    mesh.primitive(Mesh::POINTS);
    // does 1000 work on your system? how many can you make before you get a low
    // frame rate? do you need to use <1000?
    
 
    
    //  initpos = Vec3f(rnd::uniform, rnd::uniform, 0.0).normalize() * sqrt(body[_].distance / 1e10);//unitvector * body distance
    
    // initpos = Vec3f(1,0,0).normalize() * (body[_].distance);

      // mesh.vertex(initpos); // random postion btw -5~5
     
      // float m = rnd::uniform(3.0, 0.5);
    //   float m = 3 + rnd::normal() / 2;
    // float m = body[_].mass;
    //   if (m < 0.5) m = 0.5;
      

      // using a simplified volume/size relationship
       // s, t // texture example 보기

vector<Vec3f> &position(mesh.vertices());
      // separate state arrays
    //   velocity.push_back(randomVec3f(0.1));
for (int i = 0; i < body.size(); i++) {
  double t = rnd::uniform() * M_2PI; // angle = random / to make stable 0~1 random number 
  position.push_back({body[i].distance * sin(t),
                      body[i].distance * cos(t),
                      0.0});
  velocity.push_back({body[i].velocity * sin(t + M_PI / 2),
                      body[i].velocity * cos(t + M_PI / 2),
                      0.0});
  acceleration.push_back({0.0, 0.0, 0.0});
   mass.push_back(body[i].mass);
   mesh.color(randomColor());
   mesh.texCoord((body[i].radius), 0);

}


      // velocity.push_back(Vec3f(0,(body[_].velocity),0)); //y axis moving

//      acceleration.push_back(0);
    //  cout << sqrt(body[_].velocity) <<endl;
  

    nav().pos(0, 0, 70);
    //nav().near(0.1); or nav().far(0.1);

  }

  bool freeze = false;
  
  void onAnimate(double dt) override {
    if (freeze)
      return;

    // ignore the real dt and set the time step;
    dt = timeStep * 36000; //?

    
    // Calculate forces

    // XXX you put code here that calculates gravitational forces
    // These are pair-wise. Each unique pairing of two particles
    // These are equal but opposite: A exerts a force on B while B exerts that
    // same amount of force on A (but in the opposite direction!) Use a nested
    // for loop to visit each pair once The time complexity is O(n*n)
    //

    // Vec3f has lots of operations you might use...
    // • +=
    // • -=
    // • +
    // • -
    // • .normalize()  - Vec3f points in the direction as it did, but has length 1
    // • .normalize(float scale) - same but length 'scale '
    // • .mag()  - length of the Vec3f  - how long 
    // • .magSqr() - squared length of the Vec3f
    // • .dot(Vec3f f) - 
    // • .cross(Vec3f f)




    //Newton's law of universal gravitation
    vector<Vec3f> &position(mesh.vertices());
    
    for (int i = 0; i < velocity.size(); i++){
      for (int j = i + 1; j < velocity.size(); j++){
        dist = position[j] - position[i];
        force = G * body[i].mass * body[j].mass / dist.magSqr();
        unitVector = dist / dist.mag(); //to make it into vector's direction
        gravity = unitVector * force;
        acceleration[i] += gravity / mass [i]; //f = ma
        acceleration[j] -= gravity / mass [j]; // 
        

      }
      

        //  acceleration[i] -= velocity[i] * drag;

    //  if (acceleration[i].mag() > maxAcc){
    //      acceleration[i] = acceleration[i]/acceleration[i].mag() * maxAcc; //unitvertor

    //  }
     
    //  if (velocity[i].mag() > 0.1){
    //      velocity[i] = velocity[i]/velocity[i].mag() * 0.1; //unitvertor

    //  }
         velocity[i] += acceleration[i] * dt;
         position[i] += velocity[i] * dt;

// cout << position[0] <<endl;
  // cout << i << "---" << velocity[i] <<endl;
     

    }

    // cout << force << "," << dist << "," << gravity << endl;

    // clear all accelerations (IMPORTANT!!)
    for (auto &a : acceleration) a.zero();
  
///measuring distance
    for (int i = 0; i < velocity.size(); i++){
      float distance = (mesh.vertices()[i] -  nav().pos()).mag();
      // if (distance > 10) {
      //   cout << "particle:" << i << "is" <<distance << "away" << endl;
      // }
    }   
        
  }

  bool onKeyDown(const Keyboard &k) override {
    if (k.key() == ' ') {
      freeze = !freeze;
    }

    if (k.key() == '1') {
      // introduce some "random" forces
      for (int i = 0; i < velocity.size(); i++) {
        // F = ma
        acceleration[i] = randomVec3f(1) / mass[i];
      }
    }

        if (k.key() == '2') {
          warp_size = !warp_size; 

    
  
    }

  if (k.key() == '3') {
    
    warp_distance = !warp_distance;

    }


//camera setting - center
     if (k.key() == 'i') {
       Vec3f sum(0,0,0);
       // introduce some "random" forces
      for (int i = 0; i < velocity.size(); i++) {
        // F = ma
        sum+= mesh.vertices()[i];
        // acceleration[i] = randomVec3f(1) / mass[i];
      }
      sum/= velocity.size();
      nav().pos(sum);
    }
    

    return true;
  }

  void onDraw(Graphics &g) override {
    g.clear(0.3);
    g.shader(pointShader);
    g.shader().uniform("pointSize", pointSize / 100);
    g.blending(true);
    g.blendTrans();
    g.depthTesting(true);

 vector<Vec3f> &position(mesh.vertices());
    Mesh mesh;
    mesh.primitive(Mesh::POINTS);

 for (int i = 0; i < position.size(); i++){
   mesh.color(HSV((0.1666666 + float(i)/position.size())));
   if(warp_distance)
    mesh.vertex(sqrt(position[i]/ 1e10));
    else
    mesh.vertex(position[i] / 1e10);

    if (warp_size)
     mesh.texCoord(sqrt(body[i].radius / 5e4),0);
     else
     mesh.texCoord((body[i].radius / 5e4),0);
     
     cout << i << "---" << sqrt(position[i] / 1e10 )<<endl;

 }
    g.draw(mesh);

//space origin

    // mesh origin(Mesh::POINTS);
    // origin.vertex(0,0,0);
    // origin.color(1,1,1);
    // origin.texCoord(20, 0);
    // g.draw(origin);
  }
};

int main() {
AlloApp app;


app.start(); 
}

string slurp(string fileName) {
  fstream file(fileName);
  string returnValue = "";
  while (file.good()) {
    string line;
    getline(file, line);
    returnValue += line + "\n";
  }
  return returnValue;
}
