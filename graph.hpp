#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include <iostream>
#include <fstream>
#include <utility>
#include <functional>
#include <vector>
#include <string>
#include <queue>
#include <unordered_map>
#include <limits>

template <typename T>
class Graph {
 private:
  std::vector<std::unordered_map<int, T> > adjList {};
  int numVertices {};

 public:
  // empty graph with N vertices
  explicit Graph(int N);

  // construct graph from edge list in filename
  explicit Graph(const std::string& filename);

  // add an edge directed from vertex i to vertex j with given weight
  void addEdge(int i, int j, T weight);

  // removes edge from vertex i to vertex j
  void removeEdge(int i, int j);

  // is there an edge from vertex i to vertex j?
  bool isEdge(int i, int j) const;

  // return weight of edge from i to j
  // will throw an exception if there is no edge from i to j
  T getEdgeWeight(int i, int j) const;

  // returns number of vertices in the graph
  int size() const;

  // return iterator to a particular vertex
  const std::unordered_map<int, T>& neighbours(int a) const {
    return adjList.at(a);
  }
};

template <typename T>
Graph<T>::Graph(int N) : adjList(N), numVertices {N} {}

template <typename T>
Graph<T>::Graph(const std::string& inputFile) {
  std::ifstream infile {inputFile};
  if (!infile) {
    std::cerr << inputFile << " could not be opened\n";
    return;
  }
  // first line has number of vertices
  infile >> numVertices;
  adjList.resize(numVertices);
  int i {};
  int j {};
  double weight {};
  // assume each remaining line is of form
  // origin dest weight
  while (infile >> i >> j >> weight) {
    addEdge(i, j, static_cast<T>(weight));
  }
}

template <typename T>
int Graph<T>::size() const {
  return numVertices;
}

template <typename T>
void Graph<T>::addEdge(int i, int j, T weight) {
  if (i < 0 or i >= numVertices or j < 0 or j >= numVertices) {
    throw std::out_of_range("invalid vertex number");
  }
  adjList[i].insert({j, weight});
}

template <typename T>
void Graph<T>::removeEdge(int i, int j) {
  // check if i and j are valid
  if (i >= 0 && i < numVertices && j >= 0 && j < numVertices) {
    adjList[i].erase(j);
  }
}

template <typename T>
bool Graph<T>::isEdge(int i, int j) const {
  if (i >= 0 && i < numVertices && j >= 0 && j < numVertices) {
    return adjList.at(i).contains(j);
  }
  return false;
}

template <typename T>
T Graph<T>::getEdgeWeight(int i, int j) const {
  return adjList.at(i).at(j);
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const Graph<T>& G) {
  for (int i = 0; i < G.size(); ++i) {
    out << i << ':';
    for (const auto& [neighbour, weight] : G.neighbours(i)) {
      out << " (" << i << ", " << neighbour << ")[" << weight << ']';
    }
    out << '\n';
  }
  return out;
}


// APSP functions
// Use this function to return an "infinity" value
// appropriate for the type T
template <typename T>
T infinity() {
  if (std::numeric_limits<T>::has_infinity) {
    return std::numeric_limits<T>::infinity();
  } else {
    return std::numeric_limits<T>::max();
  }
}

// implement an algorithm for determining if G
// has a negative weight cycle here

/* NEGATIVE CYCLE DETECTION ALGORITHM
This uses modified Bellman-Ford to detect negative cycles anywhere in the graph.
It initialises all distances to 0 (as if there were a virtual source connected to all
vertices with weight 0), making all vertices reachable.
*/
template <typename T>
bool existsNegativeCycle(const Graph<T>& G) {
  int n = G.size();
  if (n == 0) return false;

  std::vector<T> dist(n, T(0)); // Initialises all distances to 0 instead of infinity (simulates virtual source).

  //Standard Bellman-Ford relaxation, repeats (n-1) times.
  //Ultimately, this finds shortest paths if no negative cycles exist.
  for (int i = 0; i < n - 1; i++) {
    for (int u = 0; u < n; u++) {
      for (const auto& [v, weight] : G.neighbours(u)) {
        if (dist[u] != infinity<T>() && dist[u] + weight < dist[v]) {
          dist[v] = dist[u] + weight;
        }
      }
    }
  }

  //Checks for negative cycles - if can still relax edges after (n-1) iterations, then a negative cycle is considered to exist.
  for (int u = 0; u < n; u++) {
    for (const auto& [v, weight] : G.neighbours(u)) {
      if (dist[u] != infinity<T>() && dist[u] + weight < dist[v]) {
        return true;
      }
    }
  }

  return false;
}

// implement Johnson's APSP algorithm here

/* JOHNSON'S ALL PAIRS SHORTEST PATH ALGORITHM
Combines Bellman-Ford and Dijkstra for optimal performance on sparse graphs. 
Ultimately, the edges are reweighted to make them non-negative, then uses Dijkstra from each vertex.
*/
template <typename T>
std::vector<std::vector<T> >
johnsonAPSP(const Graph<T>& G) {
  
  int n = G.size();

  // Early termination. If a negative cycle exists, no shortest paths are well-defined.
  if (existsNegativeCycle(G)) {
    return std::vector<std::vector<T> >();
  }

  else {
    /* Computes vertex potentials using Bellman-Ford.
    h[v] represents the shortest distance from virtual source to vertex v.
    */
    std::vector<T> h(n, T(0));

    // Run Bellman-Ford to find vertex potentials.
    for (int i = 0; i < n - 1; i++) {
      for (int u = 0; u < n; u++) {
        for (const auto& [v, weight] : G.neighbours(u)) {
          if (h[u] + weight < h[v]) {
            h[v] = h[u] + weight;
          }
        }
      }
    }
    // Runs Dijkstra from each vertex on reweighted graph.
    std::vector<std::vector<T> > distMatrix(n, std::vector<T>(n));

    for (int u = 0; u < n; u++) {
      // Gets the shortest distances from u in the reweighted graph.
      std::vector<T> dist = dijkstra(G, u, h);

      // Converts back to original edge weights.
      for (int v = 0; v < n; v++) {

        if (dist[v] == infinity<T>()) {
          distMatrix[u][v] = infinity<T>();
        }
        else { // Transform: reweighted_distance - h[u] + h[v] = original_distance
          distMatrix[u][v] = dist[v] - h[u] + h[v];
        }
      }
    }
    return distMatrix;
  }

}

/* DIJKSTRA'S ALGORITHM WITH REWEIGHTED EDGES
Modified Dijkstra that works with JOhnson's reweighting scheme. 
Edge weights are transformed w'(u,v) = w(u,v) + h[u] - h[v]
*/
template <typename T>
std:: vector<T> dijkstra(const Graph<T>& G, int source, const std::vector<T>& h) {
  int n = G.size();
  std::vector<T> dist(n, infinity<T>());
  std::vector<bool> visited(n, false);

  // Priority queue: pair<distance, vertex> with min-heap behaviour.
  std::priority_queue<std::pair<T, int>, std::vector<std::pair<T, int> >, std::greater<std::pair<T, int> > > pq;

  // Initialises the source vertex.
  dist[source] = T(0);
  pq.push({T(0), source});

  // Standard Dijkstra with reweighted edges.
  while (!pq.empty()) {
    auto [d, u] = pq.top();
    pq.pop();

    if (visited[u]) continue;
    visited[u] = true;

  // Relax all outgoing edges with Johnson's reweighting.
    for (const auto& [v, weight] : G.neighbours(u)) {
      T reweightedEdge = weight + h[u] - h[v]; // Apply reweighting formula
      T newDist = dist[u] + reweightedEdge;

      if (newDist < dist[v]) {
        dist[v] = newDist;
        pq.push({newDist, v});
      }
    }
  }
  return dist;
}

// implement the Floyd-Warshall APSP algorithm here
/* FLOYD-WARSHALL ALL-PAIRS SHORTEST PATH ALGORITHM
Dynamic programming approach that considers all possible intermediate vertices.
Works by gradually allowing more vertices as intermediate points in paths.
*/
template <typename T>
std::vector<std::vector<T> >
floydWarshallAPSP(const Graph<T>& G) {
  int n = G.size();
  // Initialises distance matrix where infinity is given for non-adjacent vertices.
  std::vector<std::vector<T> > distMatrix(n, std::vector<T>(n, infinity<T>()));


  // Initialises with direct edge weights.
  for (int i = 0; i < n; i++) {
    for (const auto& [j, weight] : G.neighbours(i)) {
      distMatrix[i][j] = weight;
    }
  }

  // Sets the diagonal to 0 (distance from vertex to itself)
  for (int i = 0; i < n; i++) {
    distMatrix[i][i] = 0;
  }

  // This is the main section of the Floyd-Warshall algorithm.
  // For each possible intermediate vertex k
  for (int k = 0; k < n; ++k) {
    // For each pair of vertices (i,j)
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        // Checks if the path i->k->j is shorter than the current path i->j.
        if (distMatrix[i][k] != infinity<T>() && distMatrix[k][j] != infinity<T>()) {
          T newDist = distMatrix[i][k] + distMatrix[k][j];
          if (newDist < distMatrix[i][j]) {
            distMatrix[i][j] = newDist; // Updates with the shorter path found.
          }
        }
        }
      }
    }
    return distMatrix;
}

#endif      // GRAPH_HPP_
