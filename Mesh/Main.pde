// LecturesInGraphics: smoothing
// Template for sketches
// Author: Jarek ROSSIGNAC, last edited on September 10, 2012

//**************************** global variables ****************************
pts P = new pts();
pts S = new pts();
DelaunayMesh mesh;
float t=0, f=0;
Boolean animate=false, showAcceleration=true;

//**************************** initialization ****************************
void setup() {               // executed once at the begining 
  size(800, 600);            // window size
  frameRate(30);             // render 30 frames per second
  P.declare(); // P.declare().makeGrid(3);
  P.loadPts("data/pts");
  mesh=new DelaunayMesh(P);
  
  }

//**************************** display current frame ****************************
void draw() {      // executed at each frame
  background(white); // clear screen and paints white background
  pen(black,2); // P.IDs(); 
  //P.draw();
  //pt A=P.G[0];  pt B=P.G[1];  pt C=P.G[2]; 
  //A.show();
  mesh.showVertices();
  mesh.Delaunay();
  mesh.showTriangles();
  if(scribeText) displayTextImage();
  }  // end of draw()
  
//**************************** user actions ****************************
void keyPressed() { // executed each time a key is pressed: sets the "keyPressed" and "key" state variables, 
                    // till it is released or another key is pressed or released
  if(key=='?') scribeText=!scribeText; // toggle display of help text and authors picture
//  if(key=='!') snapPicture(); // make a picture of theas
//  if(key=='o') P.loop=!P.loop;
//  if(key=='/') animate=!animate;
//  if(key==']') P.fitToCanvas();
//  if(key=='1') P.empty();
//  if(key=='c') S.copyFrom(P);
//  if(key=='S') P.savePts();   
//  if(key=='L') P.loadPts(); 
//  if(key=='s') P.savePts("data/pts");   
//  if(key=='l') P.loadPts("data/pts"); 
//  if(key=='T') {S.tuck(.5); }
//  if(key=='U') {S.tuck(-.5); }
//  if(key=='B') {S.tuck(.5); S.tuck(-.5);}

  if(key=='Q') exit();  // quit application
  }

void mousePressed() {  // executed when the mouse is pressed
  if (keyPressed && key=='i') P.pickClosestMidEdge(Mouse()); else P.pickClosest(Mouse());
  if (keyPressed) {
     if (key=='a')  P.addPt(Mouse()); 
     if (key=='i')  P.insertPt(Mouse()); 
     if (key=='d')  P.deletePickedPt(); 
     }  
  }

void mouseDragged() {
  if (!keyPressed || (key=='a') || (key=='i')) P.dragPicked();   // drag selected point with mouse
  if (keyPressed) {
      if (key=='.') {t+=2.*float(mouseX-pmouseX)/width; f=2*sq(sin(t*PI/4));}  // adjust current frame   
      if (key=='t') P.dragAll(); // move all vertices
      if (key=='r') P.rotateAllAroundCentroid(Mouse(),Pmouse()); // turn all vertices around their center of mass
      if (key=='z') P.scaleAllAroundCentroid(Mouse(),Pmouse()); // scale all vertices with respect to their center of mass
      if (key=='1' && P.nv<256) if(P.nv<1) P.addPt(Mouse()); else if(d(Mouse(),P.G[P.nv-1])>4) P.addPt(Mouse());
      }
  }  

//**************************** text for name, title and help  ****************************
String title ="Project 2 Mesh and Cutting", name ="XX",
       menu="?:txt, t:tran, r:rot, z:zoom, ]:fit, S,s:save, L,l:load, !:pic, Q:quit",
       guide="a:add, d:delete, i:insert, o:open/closed, f:filter, b:both ways, 1:draw stroke, /:animate"; // help info



