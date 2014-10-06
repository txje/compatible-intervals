# Jeremy Wang
# Sep 12, 2014

def coreScan(left, right):
  if len(left) == 0 or not len(left) == len(right):
    print "Can't core scan because left != right or left = []"
    print "left: " + str(len(left)) + ", right: " + str(len(right))
    return None
  cores = [(left[i][0], right[i][1]) for i in xrange(len(left))]
  return cores

def maxkScan(left, right, uber):
  cores = coreScan(left, right)

  # weird thing happens in CR intervals, last uber interval needs to be one further
  uber[-1] = (uber[-1][0], cores[-1][1])

  maxk = _maxkScan(uber, cores)
  return maxk


# ----------------------------------------------------------------
# MaxK data structures and computation
# ----------------------------------------------------------------

class DPNode:
  def __init__(self, a, b):
    self.next = []
    self.a = a
    self.b = b
    self.back = None
    self.total = 0
    self.back_count = 0

  def AddNext(self, node):
    overlap = self.b + 1 - node.a
    if overlap >= 0:
      self.next.append(node)
      newTotal = self.total + overlap
      if newTotal == node.total:
        node.back = self
        self.back_count += 1
      if newTotal > node.total:
        node.total = newTotal
        node.back = self
        self.back_count = 1

class DPArray:
  def __init__(self):
    self.current = []
    self.last = None
    self.start = self.current

  def NextLevel(self):
    if len(self.current) > 0:
      self.last = self.current
      self.current = []

  def AddNode(self, a, b):
    newNode = DPNode(a, b)
    self.current.append(newNode)
    if self.last != None:
      for node in self.last:
        node.AddNext(newNode)

  def LastNode(self):
    max = -1
    lastNode = None
    for node in self.current:
      if node.total >= max:
        max = node.total
        lastNode = node

    return lastNode


def _maxkScan(uber, cores):
  minDPArray = DPArray()
  last = 0
  for i in xrange(len(cores)):
    minDPArray.NextLevel()
    for j in xrange(last, len(uber)):
      if uber[j][0] > cores[i][0]:
        break
      elif (uber[j][1] >= cores[i][1]):
        minDPArray.AddNode(uber[j][0], uber[j][1])
    last = j
  invs = []
  node = minDPArray.LastNode()
  while node != None:
    invs.append((node.a, node.b))
    node = node.back
  invs.reverse()
  return invs
