class DelaunayMesh {
    pts P;
    ArrayList<pt> edges = new ArrayList<pt>();
    ArrayList<tri> t = new ArrayList<tri>();
    int nv=0, maxnv = 1000;  pt[] G = new pt [maxnv];                        // VERTICES
    int nt = 0, maxnt = maxnv*2; boolean[] visible = new boolean[maxnt];     // TRIANGLES
    int nc = 0; int[] V = new int [3*maxnt];   int[] O = new int [3*maxnt];  // CORNERS 
    int c=0;    // current corner that can be edited with keys
    pt[] points;  
  
  DelaunayMesh(pts p) {P=p; nv=p.nv; G=p.G; points=G;};
  void reset() {nv=0; nt=0; nc=0; c=0;}                                                  // removes all vertices and triangles
  void addVertex(pt P) { G[nv++].setTo(P); }                                             // adds a vertex to vertex table G
  void addTriangle(int i, int j, int k) {V[nc++]=i; V[nc++]=j; V[nc++]=k; nt=nc/3; }     // adds triangle (i,j,k) to V table
  void showVertices() {for (int i=0; i<nv; i++) {G[i].tag(Integer.toString(i));} }                          // shows all vertices as dots
  void showTriangles() { 
    //println("point1");
    for(int i = 0; i < edges.size(); i++) {
      //println((int)edges.get(i).x);
    pt p1 = points[(int)edges.get(i).x];
    pt p2 = points[(int)edges.get(i).y];
    line(p1.x, p1.y, p2.x, p2.y);
  }
}         // draws all triangles (edges, or filled)
//  void showShrunkTriangles() { for (int c=0; c<nc; c+=3) show(g(c), g(c+1), g(c+2),4); }  // shows triangles offset inwards by 4 (for filling)
//  void showTrianglesOrientation() { for (int c=0; c<nc; c++) arrow(cg(n(c)),cg(p(c)));}   // shows arrows as dart inside triangle
//  void showBorder() {for (int i=0; i<nc; i++) {if (bord(i)) {showEdge(i);}; }; };         // draws all border edges of mesh
//  void showCorner() {cg(c).show(3); };                                                    // draws current corner c

  int t (int c) {int r=int(c/3); return(r);}                 int t () {return t(c);}  // triangle of corner c
  int n (int c) {int r=3*int(c/3)+(c+1)%3; return(r);}       void n () {c=n(c); }    // next corner
  int p (int c) {int r=3*int(c/3)+(c+2)%3; return(r);}       void p () {c=p(c); }  // previous corner
  int v (int c) {return(V[c]);}                              int v () {return v(c);} // vertex of c
  int o (int c) {return(O[c]);}                              void o () {c=o(c);}   // opposite corner
  int l (int c) {return(o(n(c)));}                           void l () {c=l(c);} // left
  int r (int c) {return(o(p(c)));}                           void r () {c=r(c);}  // right
  pt g (int c) {return(G[V[c]]);};                           pt g () {return g(c);}  // shortcut to get the point where the vertex v(c) of corner c is located
//  pt cg(int c) {return g(c).makeOffset(g(p(c)),g(n(c)),-4);}; pt cg () {return cg(c);}  // computes offset location of point at corner c

  boolean nb(int c) {return(O[c]!=-1);};  // not a border corner
  boolean bord(int c) {return(O[c]==-1);};  // a border corner
  
  void Delaunay()
{
  edges.clear();
  t.clear();
  
  int leftMostIdx = 0;
  pt leftMost = points[0];
  for(int i = 1; i < nv; i++)
    if(points[i].x < leftMost.x)
    {
      leftMost = points[i];
      leftMostIdx = i;
    }
  
  PVector upVec = new PVector(0, 1);
  float lowestAng = 360f;
  int lowestIdx = 0;
  pt chosenPt = null;
  //Find Vector with lowest angle
  for(int i = 0; i < nv; i++)
  {
    if(points[i] == leftMost) continue;
    
    PVector v1 = new PVector(leftMost.x - points[i].x,
                              leftMost.y - points[i].y);

    //Not sure I have to normalize this
    v1.normalize();
    
    float ang = PVector.angleBetween(upVec, v1);
    if(ang < lowestAng)
    {
      chosenPt = points[i];
      lowestAng = ang;
      lowestIdx = i;
    }
  }
 
  //Add the first edges to the array
  PVector AB = new PVector(chosenPt.x - leftMost.x, chosenPt.y - leftMost.y);
  edges.add(new pt(leftMostIdx, lowestIdx)); 
  
  int closestIdx = getNearestPt(0, rightOrth(AB));
  
  edges.add(new pt(lowestIdx, closestIdx)); 
  edges.add(new pt(closestIdx, leftMostIdx));
 
  t.add(new tri(0, 1, 2)); 
  
  bulge(0);
  
}

void bulge(int triIdx)
{
  //if(triIdx > 2) return;
  println("Bulge: " + triIdx);
  tri T = t.get(triIdx);
  
  pt a = points[(int)edges.get(T.e2).x];
  pt b = points[(int)edges.get(T.e2).y];
  
  PVector dir = new PVector(b.x - a.x, b.y - a.y);
  int p1 = getNearestPt(T.e2, leftOrth(dir));
  if(p1 >= 0)
  {
  
    println("Left: " + edges.get(T.e2).x + " " + edges.get(T.e2).y);
    
    //first expand
    edges.add(new pt(edges.get(T.e2).x, p1)); 
    edges.add(new pt(p1, edges.get(T.e2).y));
    t.add(new tri(T.e2, edges.size() - 2, edges.size() - 1));
    println(edges.get(T.e2).x + " " + edges.get(T.e2).y);
    bulge(t.size() - 1);
  }
  
  //second expand
  a = points[(int)edges.get(T.e3).x];
  b = points[(int)edges.get(T.e3).y];
  dir = new PVector(b.x - a.x, b.y - a.y);
  p1 = getNearestPt(T.e3, leftOrth(dir));
  if(p1 >= 0)
  {
    println("Right: " + edges.get(T.e3).x + " " + edges.get(T.e3).y);
    edges.add(new pt(edges.get(T.e3).x, p1)); 
    edges.add(new pt(p1, edges.get(T.e3).y));
    t.add(new tri(T.e3, edges.size() - 2, edges.size() - 1));
    println(edges.get(T.e3).x + " " + edges.get(T.e3).y);
    bulge(t.size() - 1);
  }
}

PVector rightOrth(PVector orig)
{
  return new PVector(-orig.y, orig.x);
}

PVector leftOrth(PVector orig)
{
  return new PVector(orig.y, -orig.x);
}

int getNearestPt(int edge, PVector dir)
{
  boolean broken = false;
  pt a = points[(int)edges.get(edge).x];
  pt b = points[(int)edges.get(edge).y];
  //if(edges.get(edge).x == 9 && edges.get(edge).y == 7)
  {
    broken = true;
    //println("Broken edge!");
  }
  
  float smallestBulge = 10000f;
  int closestPt = -1;
  
  for(int i = 0; i < nv; i++)
  {
    if(i == edges.get(edge).x || i == edges.get(edge).y || isEdgeConnected(edge, i) || wouldCollide(edge, i))
    {
      if(broken && isEdgeConnected(edge, i)) println("Connected: " + i);
      if(broken && wouldCollide(edge, i)) println("Collide: " + i);
      continue;
    }
    
    PVector tmp = new PVector(points[i].x - (a.x + b.x) / 2, points[i].y - (a.y + b.y) / 2);
    float angl = PVector.angleBetween(dir, tmp);
    if(abs(angl) > PI / 2) continue;
    
    pt g = circumCenter(a, points[i], b);
    PVector AB = new PVector(b.x - a.x, b.y - a.y);
    PVector AG = new PVector(g.x - a.x, g.y - a.y);
    
    float r = distance(g, a);
    float d = AB.x * AG.y - AB.y * AG.x;
    d /= AB.mag();
    
    if((Math.abs(PVector.angleBetween(tmp, AG)) > PI / 2 && d < 0) ||
        Math.abs(PVector.angleBetween(tmp, AG)) < PI / 2 && d > 0)
      d *= -1;
      
    println(d);
    
    float bulge = r - d;
    println(i + " = " + bulge);
    if(bulge < smallestBulge)
    {
      smallestBulge = bulge;
      closestPt = i;
    }
  }
  
  println("Closest:" + closestPt);
  return closestPt;
}

boolean isEdgeConnected(int edgeIdx, int pointIndex)
{
  if(edges.get(edgeIdx).x == pointIndex || edges.get(edgeIdx).y == pointIndex) return true;
  
  int count = 0;
  for(pt edge : edges)
  {
    if(edge.x == pointIndex || edge.y == pointIndex)
    {
      if((edge.x == edges.get(edgeIdx).x || edge.y == edges.get(edgeIdx).x) ||
      (edge.x == edges.get(edgeIdx).y || edge.y == edges.get(edgeIdx).y)) count++;
    }
  }
  
  if(count > 1) return true;
  return false;
}

boolean wouldCollide(int edgeIdx, int pointIndex)
{
  pt myEdge = edges.get(edgeIdx);
  
  pt e1 = points[(int)myEdge.x];
  pt e2 = points[(int)myEdge.y];
  
  pt dest = points[pointIndex];
  
  for(pt edge : edges)
  {
    //if(edge == myEdge|| edge.x == pointIndex || edge.y == pointIndex) continue;
    //if(edge.x == myEdge.x || edge.x == myEdge.y || edge.y == myEdge.x || edge.y == myEdge.y) continue;
    
    //println("SDF");
    
    pt p1 = points[(int)edge.x];
    pt p2 = points[(int)edge.y];
    
    if(edge.x != myEdge.x && edge.y != myEdge.x)
      if(util.intersect(p1.x, p1.y, p2.x, p2.y, e1.x, e1.y, dest.x, dest.y) > 0) return true;
    if(edge.x != myEdge.y && edge.y != myEdge.y)
      if(util.intersect(p1.x, p1.y, p2.x, p2.y, e2.x, e2.y, dest.x, dest.y) > 0) return true;
  }
  
  return false;
}

float distance(pt A, pt B)
{
  return sqrt(pow(A.x - B.x, 2) + pow(A.y - B.y, 2));
}

//Taken from the class example
pt circumCenter (pt A, pt B, pt C) {    // computes the center of a circumscirbing circle to triangle (A,B,C)
  PVector AB =  new PVector(B.x - A.x, B.y - A.y);  float ab2 = PVector.dot(AB,AB);
  PVector AC =  new PVector(C.x - A.x, C.y - A.y); 
  float tmpX = AC.x; AC.x = -AC.y; AC.y = tmpX;
  
  float ac2 = PVector.dot(AC,AC);
  float d = 2*PVector.dot(AB,AC);
  tmpX = AB.x; AB.x = -AB.y; AB.y = tmpX;
  AB.mult(-ac2); 
  AC.mult(ab2);
  AB.add(AC);
  AB.mult(1./d);
  pt X =  new pt(A.x, A.y);
  X.x += AB.x;
  X.y += AB.y;
  return(X);
  };
  
  
//  void showCorner(int c, float r) {cg(c).show(r); };   // renders corner c as small ball
//  void showEdge(int c) {g(p(c)).to(g(n(c))); };  // draws edge of t(c) opposite to corner c

// void triangulate() {     // performs Delaunay triangulation using a quartic algorithm
//   c=0;                   // to reset current corner
//   pt X = new pt(0,0);
//   float r=1;
//   for (int i=0; i<nv-2; i++) for (int j=i+1; j<nv-1; j++) for (int k=j+1; k<nv; k++) {
//      X=circumCenter (G[i],G[j],G[k]);  r = X.disTo(G[i]);
//      boolean found=false; 
//      for (int m=0; m<nv; m++) if ((m!=i)&&(m!=j)&&(m!=k)&&(X.disTo(G[m])<=r)) found=true;  
//     if (!found) {if (ccw(G[i],G[j],G[k])) addTriangle(i,j,k); else addTriangle(i,k,j);};
//     }; 
//   }  
//
//  void computeOslow() {                      // slow method to set the O table from the V table, assumes consistent orientation of tirangles
//    for (int i=0; i<3*nt; i++) {O[i]=-1;};  // init O table to -1: has no opposite (i.e. is a border corner)
//    for (int i=0; i<3*nt; i++) {  for (int j=i+1; j<3*nt; j++) {       // for each corner i, for each other corner j
//      if( (v(n(i))==v(p(j))) && (v(p(i))==v(n(j))) ) {O[i]=j; O[j]=i;};};}; // make i and j opposite if they match 
//   }
//  
//  void computeO() {                                          // faster method fo r computin gO
//    int nIC [] = new int [maxnv];                            // number of incident corners on each vertex
//    println("COMPUTING O: nv="+nv +", nt="+nt +", nc="+nc );
//    int maxValence=0;
//    for (int c=0; c<nc; c++) {O[c]=-1;};                      // init O table to -1: has no opposite (i.e. is a border corner)
//    for (int v=0; v<nv; v++) {nIC[v]=0; };                    // init the valence value for each vertex to 0
//    for (int c=0; c<nc; c++) {nIC[v(c)]++;}                   // computes vertex valences
//    for (int v=0; v<nv; v++) {if(nIC[v]>maxValence) {maxValence=nIC[v]; };};  println(" Max valence = "+maxValence+". "); // computes and prints maximum valence 
//    int IC [][] = new int [maxnv][maxValence];                 // declares 2D table to hold incident corners (htis can be folded into a 1D table !!!!!)
//    for (int v=0; v<nv; v++) {nIC[v]=0; };                     // resets the valence of each vertex to 0 . It will be sued as a counter of incident corners.
//    for (int c=0; c<nc; c++) {IC[v(c)][nIC[v(c)]++]=c;}        // appends incident corners to corresponding vertices     
//    for (int c=0; c<nc; c++) {                                 // for each corner c
//      for (int i=0; i<nIC[v(p(c))]; i++) {                     // for each incident corner a of the vertex of the previous corner of c
//        int a = IC[v(p(c))][i];      
//        for (int j=0; j<nIC[v(n(c))]; j++) {                   // for each other corner b in the list of incident corners to the previous corner of c
//           int b = IC[v(n(c))][j];
//           if ((b==n(a))&&(c!=n(b))) {O[c]=n(b); O[n(b)]=c; };  // if a and b have matching opposite edges, make them opposite
//           };
//        };
//      };
//    } // end computeO
//
//int countloops() {                                      // returns number of bounding loops (holes) in mesh
//  boolean M[] = new boolean [3*nt];			// flag for marking visited corners
//  for (int c=0; c<3*nt; c++) M[c]=false;		// all corners not visited initially
//  int loops=0;						// counter of loops
//  for (int b=0; b<3*nt; b++) {			 	// look for not marked border corners
//     int c=b;                                           // use c to trace loop
//     if ( (!M[c]) && (o(c) == -1) ) {			// found not marked and border
//	loops++;					// new loop
//        while (!M[c]) {					// while not finished tracing
//           M[c]=true;					// mark it
//           c=n(c);					// start with the next corner
//           while ( o(c) != -1 ) c=n(o(c)); 		// keep turning until a border corner is reached 
//     }; }; };
//  return(loops);
//  }
//  
//  pt triCenter() {return triCenter(c);}                         // returns center of mass of triangle of current corner 
//  pt triCenter(int c) {return center(g(c),g(n(c)),g(p(c))); }  // returns center of mass of triangle of corner c
//  pt triCircumcenter() {return triCircumcenter(c);}                         // returns circumcenter of triangle of current corner 
//  pt triCircumcenter(int c) {return circumCenter(g(c),g(n(c)),g(p(c))); }  // returns circumcenter of triangle of corner c
//  void excludeInvisibleTriangles () {for (int b=0; b<nc; b++) {if (!visible[t(o(b))]) {O[b]=-1;};};} // sets o of opposite corners to -1 for triangles marked as invisible
//  void showDual () {for (int b=0; b<nc; b++) if (nb(b)) if (b<o(b)) triCircumcenter(b).to(triCircumcenter(o(b)));}
//
//  void compactVO() {                                                // compacts V and O tables
//    println("COMPACT TRIANGLES: nv="+nv +", nt="+nt +", nc="+nc );
//    int[] U = new int [nc];
//    int lc=-1; for (int c=0; c<nc; c++) {if (visible[t(c)]) {U[c]=++lc; }; };
//    for (int c=0; c<nc; c++) {if (nb(c)) {O[c]=U[o(c)];} else {O[c]=-1;}; };
//    int lt=0;
//    for (int t=0; t<nt; t++)  if (visible[t]) {
//        V[3*lt]=V[3*t]; V[3*lt+1]=V[3*t+1]; V[3*lt+2]=V[3*t+2]; 
//        O[3*lt]=O[3*t]; O[3*lt+1]=O[3*t+1]; O[3*lt+2]=O[3*t+2]; 
//        lt++;
//        };
//    nt=lt; nc=3*nt;
//    println("      ...  NOW: nv="+nv +", nt="+nt +", nc="+nc );
//    }
//
//  void compactV() {                                                 // compacts the G table to remove deleted trangles
//    println("COMPACT VERTICES: nv="+nv +", nt="+nt +", nc="+nc );
//    int[] U = new int [nv];
//    boolean[] deleted = new boolean [nv];
//    for (int v=0; v<nv; v++) {deleted[v]=true;};
//    for (int c=0; c<nc; c++) {deleted[v(c)]=false;};
//    int lv=-1; for (int v=0; v<nv; v++) {if (!deleted[v]) {U[v]=++lv; }; };
//    for (int c=0; c<nc; c++) {V[c]=U[v(c)]; };
//    lv=0;
//    for (int v=0; v<nv; v++) if (!deleted[v]) {G[lv].setTo(G[v]);  deleted[lv]=false; lv++; };
//    nv=lv;
//    println("      ...  NOW: nv="+nv +", nt="+nt +", nc="+nc );
//    }
//
//void flip () {flip(c);}                     // flips edge of current corner
//void flip(int c) {V[n(o(c))]=v(c); V[n(c)]=v(o(c)); int co=o(c); O[co]=r(c); O[r(c)]=co; O[c]=r(co); O[r(co)]=c; O[p(c)]=p(co); O[p(co)]=p(c); }  // fip edge opposite to corner c
//void doFlips() {for (int c=0; c<nc; c++) {if (g(n(c)).disTo(g(p(c)))>g(c).disTo(g(o(c)))) {flip(c);}; };  } // assumes manifold and flips all bad edges
//
//void collapse() {collapse(c);}  // collapses edge opposite to current corner
//void collapse(int c) {      // collapse edge opposite to corner c for simplification 
//    int b=p(c), oc=o(c), vnc=v(n(c));
//    for (int a=b; a!=n(oc); a=p(r(a))) {V[a]=vnc;}; V[p(c)]=vnc; V[n(oc)]=vnc; 
//    O[l(c)]=r(c); O[r(c)]=l(c);     O[l(oc)]=r(oc); O[r(oc)]=l(oc); 
//  }

} // end of MESH
class tri
{
  int e1, e2, e3;
  
  tri(int e1, int e2, int e3)
  {
    this.e1 = e1;
    this.e2 = e2;
    this.e3 = e3;
  }
}
