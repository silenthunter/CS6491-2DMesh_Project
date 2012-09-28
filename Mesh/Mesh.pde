
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
  points[4] = new pt(400, 325);
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
  
  int closestIdx = getNearestPt(0, new PVector(AB.y, AB.x));
  println(closestIdx);
  
  edges.add(new pt(leftMostIdx, closestIdx)); 
  edges.add(new pt(closestIdx, lowestIdx));
 
  t.add(new tri(0, 1, 2)); 
  
  bulge(0);
  
}

void bulge(int triIdx)
{
  println("Bulge: " + triIdx);
  tri T = t.get(triIdx);
  
  pt a = points[edges.get(T.e2).x];
  pt b = points[edges.get(T.e2).y];
  
  PVector dir = new PVector(b.y - a.y, b.x - a.x);
  int p1 = getNearestPt(T.e2, dir);
  if(p1 < 0) return;
  
  println(edges.get(T.e2).x + " " + edges.get(T.e2).y);
  
  //first expand
  edges.add(new pt(edges.get(T.e2).x, p1)); 
  edges.add(new pt(p1, edges.get(T.e2).y));
  t.add(new tri(T.e2, edges.size() - 2, edges.size() - 1));
  bulge(t.size() - 1);
  
  //second expand
  a = points[edges.get(T.e3).x];
  b = points[edges.get(T.e3).y];
  dir = new PVector(b.y - a.y, b.x - a.x);
  getNearestPt(T.e3, dir);
}

int getNearestPt(int edge, PVector dir)
{
  pt a = points[edges.get(edge).x];
  pt b = points[edges.get(edge).y];
  
  float smallestBulge = 0000f;
  int closestPt = -1;
  
  for(int i = 0; i < points.length; i++)
  {
    if(i == edges.get(edge).x || i == edges.get(edge).y || isEdgeConnected(edge, i)) continue;
    
    PVector tmp = new PVector(points[i].x - a.x, points[i].y - a.y);
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
  for(pt edge: edges)
  {
    if(edge.x == pointIndex || edge.y == pointIndex) return true;
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