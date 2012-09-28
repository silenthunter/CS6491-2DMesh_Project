
pt[] points = new pt[10];
ArrayList<pt> edges = new ArrayList<pt>();
ArrayList<tri> t = new ArrayList<tri>();

float pointSize = 5.0f;

void setup()
{
  size(800,600);
  points[0] = new pt(100, 100);
  points[1] = new pt(50, 100);
  points[2] = new pt(120, 300);
  points[3] = new pt(300, 500);
  points[4] = new pt(400, 125);
  points[5] = new pt(350, 275);
  points[6] = new pt(475, 325);
  points[7] = new pt(250, 250);
  points[8] = new pt(425, 300);
  points[9] = new pt(100, 400);
  
  Delaunay();
}

void draw()
{
  background(255);
  for(int i = 0; i < points.length; i++)
  {
    ellipse(points[i].x, points[i].y, pointSize, pointSize);
    fill(0);
    text("" + i, points[i].x + 10, points[i].y);
  }
  
  for(int i = 0; i < edges.size(); i++)
  {
    pt p1 = points[edges.get(i).x];
    pt p2 = points[edges.get(i).y];
    line(p1.x, p1.y, p2.x, p2.y);
  }
}

void Delaunay()
{
  int leftMostIdx = 0;
  pt leftMost = points[0];
  for(int i = 1; i < points.length; i++)
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
  for(int i = 0; i < points.length; i++)
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
  
  pt a = points[edges.get(T.e2).x];
  pt b = points[edges.get(T.e2).y];
  
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
  a = points[edges.get(T.e3).x];
  b = points[edges.get(T.e3).y];
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
  pt a = points[edges.get(edge).x];
  pt b = points[edges.get(edge).y];
  
  float smallestBulge = 0000f;
  int closestPt = -1;
  
  for(int i = 0; i < points.length; i++)
  {
    if(i == edges.get(edge).x || i == edges.get(edge).y || isEdgeConnected(edge, i) || wouldCollide(edge, i)) continue;
    
    PVector tmp = new PVector(points[i].x - (a.x + b.x) / 2, points[i].y - (a.y + b.y) / 2);
    float angl = PVector.angleBetween(dir, tmp);
    if(abs(angl) > PI / 2) continue;
    
    pt g = circumCenter(a, b, points[i]);
    PVector AB = new PVector(b.x - a.x, b.y - a.y);
    PVector AG = new PVector(g.x - a.x, g.y - a.y);
    
    float r = distance(g, a);
    PVector d = AB.cross(AG);
    //println(d);
    //AB.normalize();
    d.div(AB.mag());
    //println(d);
    
    float bulge = r - d.mag();
    //println(bulge);
    if(bulge > smallestBulge)
    {
      smallestBulge = bulge;
      closestPt = i;
    }
  }
  
  return closestPt;
}

boolean isEdgeConnected(int edgeIdx, int pointIndex)
{
  if(edges.get(edgeIdx).x == pointIndex || edges.get(edgeIdx).y == pointIndex) return true;
  
  return false;
}

boolean wouldCollide(int edgeIdx, int pointIndex)
{
  pt myEdge = edges.get(edgeIdx);
  
  pt e1 = points[myEdge.x];
  pt e2 = points[myEdge.y];
  
  pt dest = points[pointIndex];
  
  for(pt edge : edges)
  {
    //if(edge == myEdge|| edge.x == pointIndex || edge.y == pointIndex) continue;
    //if(edge.x == myEdge.x || edge.x == myEdge.y || edge.y == myEdge.x || edge.y == myEdge.y) continue;
    
    //println("SDF");
    
    pt p1 = points[edge.x];
    pt p2 = points[edge.y];
    
    if(util.intersect(p1.x, p1.y, p2.x, p2.y, e1.x, e1.y, dest.x, dest.y) > 0) return true;
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

class pt
{
  public int x,y;
  
  pt(int x, int y)
  {
    this.x = x;
    this.y = y;
  }
}

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
