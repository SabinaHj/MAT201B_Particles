// Sabina Hyoju Ahn
// particles-p4
// MAT201B_2022_Assignment2-4

#include "al/app/al_App.hpp"
#include "al/app/al_GUIDomain.hpp"
#include "al/math/al_Random.hpp"

using namespace al;

#include <fstream>
#include <vector>
using namespace std;



Vec3f randomVec3f(float scale) {
  return Vec3f(rnd::uniformS(), rnd::uniformS(), rnd::uniformS()) * scale; //uniform -1 ~ 1
}


string slurp(string fileName); // forward declaration //slurp: file name -> string for shader

struct AlloApp : App {
    
  Parameter pointSize{"/pointSize", "", 0.5, "", 0.1, 1.1};
  Parameter timeStep{"/timeStep", "", 0.5, "", 0.1, 1};
  Parameter gravityRandomConst{"/gravityRandomConst", "", 0.01, "", 0.005, 0.2};
  Parameter drag{"/drag", "", 0.3, "", 0.01, 0.6};
  Parameter maxAcc{"/maxAcc", "", 0.1, "", 0.08, 1};

  ShaderProgram pointShader;

  //  simulation state
  Mesh mesh; // position *is inside the mesh* mesh.vertices() are the positions
  vector<Vec3f> velocity;
  vector<Vec3f> acceleration;
  vector<float> mass;
  
  float gravityConst[1000];

  float force;
  // float gravityConst;

  Vec3f unitVector;
  Vec3f dist;
  Vec3f gravity;

  void onInit() override {
    // set up GUI
    auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
    auto &gui = GUIdomain->newGUI();
    gui.add(pointSize); // add parameter to GUI
    gui.add(timeStep);  // add parameter to GUI
    gui.add(gravityRandomConst);
    gui.add(drag);
    gui.add(maxAcc);
     
  }

  void onCreate() override {
    navControl().useMouse(false);

    // compile shaders
    pointShader.compile(slurp("../point-vertex.glsl"),
                        slurp("../point-fragment.glsl"),
                        slurp("../point-geometry.glsl"));

    // set initial conditions of the simulation
    //

    // c++11 "lambda" function
    auto randomColor = []() { 
      // return HSV(rnd::uniform(), 1.0f, 1.0f);
      return HSV(rnd::uniform(), rnd::uniform(),rnd::uniform()); 
      };

    mesh.primitive(Mesh::POINTS);
    // does 1000 work on your system? how many can you make before you get a low
    // frame rate? do you need to use <1000?
    for (int _ = 0; _ < 1000; _++) {
      mesh.vertex(randomVec3f(5)); // random postion btw -5~5
      mesh.color(randomColor());

      // float m = rnd::uniform(3.0, 0.5);
      float m = 3 + rnd::normal() / 2;
      if (m < 0.5) m = 0.5;
      mass.push_back(m);

      // using a simplified volume/size relationship
      mesh.texCoord(pow(m, 1.0f / 3), 0); // s, t // texture example 보기

      // separate state arrays
      velocity.push_back(randomVec3f(0.1));
      acceleration.push_back(randomVec3f(1));
    }

    nav().pos(0, 0, 1);
    //nav().near(0.1); or nav().far(0.1);

    for (int i = 0; i < velocity.size(); i++) {
        gravityConst[i] = 0.01;
    }
   
  }

  bool freeze = false;
  
  void onAnimate(double dt) override {
    if (freeze)
      return;

    // ignore the real dt and set the time step;
    dt = timeStep;

    
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
        force = gravityConst[i] / dist.magSqr();
        unitVector = dist / dist.mag(); //to make it into vector's direction
        gravity = unitVector * force;
        acceleration[i] += gravity;
        

      }
      

         acceleration[i] -= velocity[i] * drag;

     if (acceleration[i].mag() > maxAcc){
        //  acceleration[i] = acceleration[i]/acceleration[i].mag() * maxAcc; //unitvertor
            acceleration[i]= force * (1. / 10) * maxAcc; //F = ma -> a = F/m
     }
     
    //  if (velocity[i].mag() > 0.1){
    //      velocity[i] = velocity[i]/velocity[i].mag() * 0.1; //unitvertor

    //  }
        //  velocity[i] += acceleration[i] * dt;
        //  position[i] += velocity[i] * dt;

        position[i] = clone(position[i]).rotate(M_PI / 2) * Vec3f(1, 1, -1) * 0.2;
        velocity[i] += acceleration[i] * dt * clone(position[i]).rotate(M_PI);
        //  position[i] += velocity[i] * dt;
        position[i] +=  velocity[i] * dt;
         

    }

    // cout << force << "," << dist << "," << gravity << endl;

    // clear all accelerations (IMPORTANT!!)
    for (auto &a : acceleration) a.zero();
    
     cout<< "push key 2 and 3" <<endl;


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

      for (int i = 0; i < velocity.size(); i++) {
        // F = ma
        //gravityConst[i] = ((10 * rnd::uniformS()1) * 0.01;
        gravityConst[i] = rnd::uniform(1, 10)*gravityRandomConst; //ranadom gravity
      }
      cout << gravityConst[1] << "," << gravityRandomConst << endl;
    }

  if (k.key() == '3') {
    
      for (int i = 0; i < velocity.size(); i++) {
        gravityConst[i] = 0.01;
      }
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
    g.draw(mesh);


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