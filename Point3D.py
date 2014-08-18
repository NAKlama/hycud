import math


class Point3D:
  """A ppoint in 3D space capable of basic vector mathematics"""
  def __init__(self, x=0.0, y=0.0, z=0.0, point=None):
    if point is None:
      self.x = x
      self.y = y
      self.z = z
    else:
      self.x = point.x
      self.y = point.y
      self.z = point.z

  def magnitude(self):
    return math.sqrt(self.x**2 + self.y**2 + self.z**2)

  def dotProduct(self, a):
    return a.x * self.x + a.y * self.y + a.z * self.z

  def mult(self, n):
    return Point3D( x=(self.x * n)
                  , y=(self.y * n)
                  , z=(self.z * n)
                  )

  def multTo(self, n):
    self.x *= n
    self.y *= n
    self.z *= n

  def add(self, p):
    return Point3D( x=(self.x + p.x)
                  , y=(self.y + p.y)
                  , z=(self.z + p.z)
                  )

  def addTo(self, p):
    self.x += p.x
    self.y += p.y
    self.z += p.z

  def div(self, n):
    return Point3D( x=(self.x / n)
                  , y=(self.y / n)
                  , z=(self.z / n)
                  )

  def dist(self, p):
    return math.sqrt( (self.x - p.x)**2
                    + (self.y - p.y)**2
                    + (self.z - p.z)**2
                    )
