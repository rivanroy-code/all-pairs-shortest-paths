#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include "graph.hpp"

// *** Negative Cycle Test Cases

// negative cycle 0, 1, 0
TEST(negativeCycleTest, twoVerticesNegativeLoop) {
  Graph<int> G {2};
  G.addEdge(0, 1, 1);
  G.addEdge(1, 0, -3);
  ASSERT_TRUE(existsNegativeCycle(G));
}

// no negative cycle
TEST(negativeCycleTest, twoVerticesNoNegativeLoop) {
  Graph<int> G {2};
  G.addEdge(0, 1, 5);
  G.addEdge(1, 0, -3);
  EXPECT_FALSE(existsNegativeCycle(G));
}

// negative cycle 0, 1, 2, 0
TEST(negativeCycleTest, negativeTriangle) {
  Graph<int> G {3};
  G.addEdge(0, 1, 1);
  G.addEdge(1, 2, -3);
  G.addEdge(2, 1, 1);
  ASSERT_TRUE(existsNegativeCycle(G));
}

// no negative cycle
TEST(negativeCycleTest, zeroWeightTriangle) {
  Graph<int> G {3};
  G.addEdge(0, 1, 1);
  G.addEdge(1, 2, -2);
  G.addEdge(2, 0, 1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

// no negative cycle
TEST(negativeCycleTest, fourVerticesPositiveCycle) {
  Graph<int> G {4};
  G.addEdge(0, 1, 1);
  G.addEdge(1, 2, -2);
  G.addEdge(2, 3, 1);
  G.addEdge(3, 0, 1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

// no negative cycle
TEST(negativeCycleTest, fourVerticesZeroCycle) {
  Graph<int> G {4};
  G.addEdge(0, 1, 1);
  G.addEdge(1, 2, -3);
  G.addEdge(2, 3, 1);
  G.addEdge(3, 0, 1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

// negative cycle 0, 1, 2, 3, 0
TEST(negativeCycleTest, fourVerticesNegativeCycle) {
  Graph<int> G {4};
  G.addEdge(0, 1, 1);
  G.addEdge(1, 2, -2);
  G.addEdge(2, 3, -1);
  G.addEdge(3, 0, 1);
  ASSERT_TRUE(existsNegativeCycle(G));
}

// no negative cycle
TEST(negativeCycleTest, fourVerticesMoreEdges) {
  Graph<int> G {4};
  G.addEdge(0, 1, 1);
  G.addEdge(1, 0, 1);
  G.addEdge(1, 2, -2);
  G.addEdge(1, 3, -1);
  G.addEdge(2, 1, 3);
  G.addEdge(2, 3, -1);
  G.addEdge(3, 0, 2);
  G.addEdge(0, 3, 1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

// negative cycle 1, 2, 3, 1
TEST(negativeCycleTest, fourVerticesNegativeCycleAgain) {
  Graph<int> G {4};
  G.addEdge(0, 1, 1);
  G.addEdge(1, 0, 1);
  G.addEdge(1, 2, -2);
  G.addEdge(1, 3, -1);
  G.addEdge(3, 1, 1);
  G.addEdge(2, 1, 3);
  G.addEdge(2, 3, -1);
  G.addEdge(3, 0, 2);
  G.addEdge(0, 3, 1);
  ASSERT_TRUE(existsNegativeCycle(G));
}

// negative cycle 1, 2, 4, 3, 1
TEST(negativeCycleTest, fiveVerticesNegativeCycle) {
  Graph<int> G {5};
  G.addEdge(0, 1, 2);
  G.addEdge(1, 2, -2);
  G.addEdge(2, 4, 2);
  G.addEdge(4, 3, 1);
  G.addEdge(3, 1, -2);
  G.addEdge(1, 0, -1);
  G.addEdge(0, 4, -1);
  ASSERT_TRUE(existsNegativeCycle(G));
}

// function to create a random directed and weighted graph with N vertices
// p is the probability of an edge between any two vertices
Graph<int> createRandomGraph(int N, unsigned seed, double p = 0.5) {
  std::mt19937 mt {seed};
  std::uniform_int_distribution<int> vertexDist {0, N - 1};
  std::binomial_distribution<int> heads {1, p};
  // edge weights are randomly chosen from -1 to 100
  std::uniform_int_distribution<int> weight {-1, 100};
  Graph<int> G {N};
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (heads(mt)) {
        int edgeWeight = weight(mt);
        G.addEdge(i, j, edgeWeight);
      }
    }
  }
  return G;
}

// Random test cases on size 100 graphs
TEST(negativeCycleTest, random100_0) {
  Graph<int> G = createRandomGraph(100, 1198988931, 0.15);
  EXPECT_TRUE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_1) {
  Graph<int> G = createRandomGraph(100, 2495118394, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_2) {
  Graph<int> G = createRandomGraph(100, 340189135, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_3) {
  Graph<int> G = createRandomGraph(100, 3635050020, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_4) {
  Graph<int> G = createRandomGraph(100, 1756526730, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_5) {
  Graph<int> G = createRandomGraph(100, 4164912679, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_6) {
  Graph<int> G = createRandomGraph(100, 2449556836, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_7) {
  Graph<int> G = createRandomGraph(100, 3351942058, 0.15);
  EXPECT_TRUE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_8) {
  Graph<int> G = createRandomGraph(100, 1447474330, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_9) {
  Graph<int> G = createRandomGraph(100, 3456554122, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_10) {
  Graph<int> G = createRandomGraph(100, 2240266931, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_11) {
  Graph<int> G = createRandomGraph(100, 3463085465, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_12) {
  Graph<int> G = createRandomGraph(100, 1424997335, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_13) {
  Graph<int> G = createRandomGraph(100, 49815473, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_14) {
  Graph<int> G = createRandomGraph(100, 4048150027, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_15) {
  Graph<int> G = createRandomGraph(100, 4052540741, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_16) {
  Graph<int> G = createRandomGraph(100, 4289914535, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_17) {
  Graph<int> G = createRandomGraph(100, 3394890305, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_18) {
  Graph<int> G = createRandomGraph(100, 3832010753, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random100_19) {
  Graph<int> G = createRandomGraph(100, 931073880, 0.15);
  EXPECT_FALSE(existsNegativeCycle(G));
}

// random test cases on size 200 graphs

TEST(negativeCycleTest, random200_0) {
  Graph<int> G = createRandomGraph(200, 850177966, 0.1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_1) {
  Graph<int> G = createRandomGraph(200, 793222184, 0.1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_2) {
  Graph<int> G = createRandomGraph(200, 1283402404, 0.1);
  EXPECT_TRUE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_3) {
  Graph<int> G = createRandomGraph(200, 523281336, 0.1);
  EXPECT_TRUE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_4) {
  Graph<int> G = createRandomGraph(200, 1476350622, 0.1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_5) {
  Graph<int> G = createRandomGraph(200, 4156308090, 0.1);
  EXPECT_TRUE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_6) {
  Graph<int> G = createRandomGraph(200, 3658781174, 0.1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_7) {
  Graph<int> G = createRandomGraph(200, 2789138855, 0.1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_8) {
  Graph<int> G = createRandomGraph(200, 4062547577, 0.1);
  EXPECT_TRUE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_9) {
  Graph<int> G = createRandomGraph(200, 2266342153, 0.1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_10) {
  Graph<int> G = createRandomGraph(200, 4249418288, 0.1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_11) {
  Graph<int> G = createRandomGraph(200, 3611979686, 0.1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_12) {
  Graph<int> G = createRandomGraph(200, 2618170795, 0.1);
  EXPECT_TRUE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_13) {
  Graph<int> G = createRandomGraph(200, 917475276, 0.1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_14) {
  Graph<int> G = createRandomGraph(200, 3470194655, 0.1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_15) {
  Graph<int> G = createRandomGraph(200, 2884946528, 0.1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_16) {
  Graph<int> G = createRandomGraph(200, 1443499060, 0.1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_17) {
  Graph<int> G = createRandomGraph(200, 1432307034, 0.1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_18) {
  Graph<int> G = createRandomGraph(200, 3367014845, 0.1);
  EXPECT_FALSE(existsNegativeCycle(G));
}

TEST(negativeCycleTest, random200_19) {
  Graph<int> G = createRandomGraph(200, 4020578704, 0.1);
  EXPECT_TRUE(existsNegativeCycle(G));
}

// *** End of negative cycle tests

// some machinery for testing the APSP functions

// check that a given vector of distances is optimal
// checks that no distance can be improved via relaxation
template <typename T>
bool allEdgesRelaxed(const std::vector<T>& bestDistanceTo, const Graph<T>& G,
                      int source) {
  if (bestDistanceTo.at(source) != T {}) {
    return false;
  }
  T inf = std::numeric_limits<T>::max();
  for (int vertex = 0; vertex < G.size(); ++vertex) {
    if (bestDistanceTo.at(vertex) == inf) {
      continue;
    }
    for (const auto& [neighbour, weight] : G.neighbours(vertex)) {
      T distViaVertex = bestDistanceTo.at(vertex) + weight;
      if (bestDistanceTo.at(neighbour) > distViaVertex) {
        return false;
      }
    }
  }
  return true;
}

// type to represent johnsonAPSP or floydWarshallAPSP
using apspFunction =
  std::function<std::vector<std::vector<int> >(const Graph<int>&)>;

// function to create a random graph and test the output of f on it
void randomTest(apspFunction f, int N, unsigned seed, double p = 0.5) {
  std::mt19937 mt {seed};
  std::uniform_int_distribution<int> vertexDist {0, N - 1};
  std::binomial_distribution<int> heads {1, p};
  std::uniform_int_distribution<int> weight {-1, 100};
  Graph<int> G {N};
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (heads(mt)) {
        int edgeWeight = weight(mt);
        G.addEdge(i, j, edgeWeight);
      }
    }
  }
  ASSERT_FALSE(existsNegativeCycle(G));
  std::vector<std::vector<int> > distanceMatrix = f(G);
  for (int v = 0; v < G.size(); ++v) {
    ASSERT_TRUE(allEdgesRelaxed(distanceMatrix.at(v), G, v));
  }
}

// ***Tests of johnsonAPSP

TEST(johnsonTest, dijkstraFailure) {
  Graph<int> G {4};
  G.addEdge(0, 1, 1);
  G.addEdge(0, 2, 4);
  G.addEdge(2, 1, -4);
  G.addEdge(1, 3, 1);
  ASSERT_FALSE(existsNegativeCycle(G));
  std::vector<std::vector<int> > distanceMatrix = johnsonAPSP<int>(G);
  for (int v = 0; v < G.size(); ++v) {
    ASSERT_TRUE(allEdgesRelaxed(distanceMatrix.at(v), G, v));
  }
}

TEST(johnsonTest, dijkstraFailure2) {
  Graph<int> G {5};
  G.addEdge(0, 1, 1);
  G.addEdge(0, 2, 4);
  G.addEdge(2, 1, -4);
  G.addEdge(1, 3, 1);
  G.addEdge(3, 4, 1);
  ASSERT_FALSE(existsNegativeCycle(G));
  std::vector<std::vector<int> > distanceMatrix = johnsonAPSP<int>(G);
  for (int v = 0; v < G.size(); ++v) {
    ASSERT_TRUE(allEdgesRelaxed(distanceMatrix.at(v), G, v));
  }
}

// Test on graph in the file tinyEWD.txt
TEST(johnsonTest, tinyJohnson) {
  Graph<int> G {"tinyEWD.txt"};
  std::vector<std::vector<int> > result = johnsonAPSP(G);
  ASSERT_EQ(result.at(0).at(0), 0);
  ASSERT_EQ(result.at(0).at(1), 105);
  ASSERT_EQ(result.at(0).at(2), 26);
  ASSERT_EQ(result.at(0).at(3), 99);
  ASSERT_EQ(result.at(0).at(4), 38);
  ASSERT_EQ(result.at(0).at(5), 73);
  ASSERT_EQ(result.at(0).at(6), 151);
  ASSERT_EQ(result.at(0).at(7), 60);
  ASSERT_EQ(result.at(1).at(0), 139);
  ASSERT_EQ(result.at(1).at(1), 0);
  ASSERT_EQ(result.at(1).at(2), 121);
  ASSERT_EQ(result.at(1).at(3), 29);
  ASSERT_EQ(result.at(1).at(4), 174);
  ASSERT_EQ(result.at(1).at(5), 183);
  ASSERT_EQ(result.at(1).at(6), 81);
  ASSERT_EQ(result.at(1).at(7), 155);
  ASSERT_EQ(result.at(2).at(0), 183);
  ASSERT_EQ(result.at(2).at(1), 94);
  ASSERT_EQ(result.at(2).at(2), 0);
  ASSERT_EQ(result.at(2).at(3), 73);
  ASSERT_EQ(result.at(2).at(4), 97);
  ASSERT_EQ(result.at(2).at(5), 62);
  ASSERT_EQ(result.at(2).at(6), 125);
  ASSERT_EQ(result.at(2).at(7), 34);
  ASSERT_EQ(result.at(3).at(0), 110);
  ASSERT_EQ(result.at(3).at(1), 186);
  ASSERT_EQ(result.at(3).at(2), 92);
  ASSERT_EQ(result.at(3).at(3), 0);
  ASSERT_EQ(result.at(3).at(4), 145);
  ASSERT_EQ(result.at(3).at(5), 154);
  ASSERT_EQ(result.at(3).at(6), 52);
  ASSERT_EQ(result.at(3).at(7), 126);
  ASSERT_EQ(result.at(4).at(0), 186);
  ASSERT_EQ(result.at(4).at(1), 67);
  ASSERT_EQ(result.at(4).at(2), 168);
  ASSERT_EQ(result.at(4).at(3), 76);
  ASSERT_EQ(result.at(4).at(4), 0);
  ASSERT_EQ(result.at(4).at(5), 35);
  ASSERT_EQ(result.at(4).at(6), 128);
  ASSERT_EQ(result.at(4).at(7), 37);
  ASSERT_EQ(result.at(5).at(0), 171);
  ASSERT_EQ(result.at(5).at(1), 32);
  ASSERT_EQ(result.at(5).at(2), 153);
  ASSERT_EQ(result.at(5).at(3), 61);
  ASSERT_EQ(result.at(5).at(4), 35);
  ASSERT_EQ(result.at(5).at(5), 0);
  ASSERT_EQ(result.at(5).at(6), 113);
  ASSERT_EQ(result.at(5).at(7), 28);
  ASSERT_EQ(result.at(6).at(0), 58);
  ASSERT_EQ(result.at(6).at(1), 134);
  ASSERT_EQ(result.at(6).at(2), 40);
  ASSERT_EQ(result.at(6).at(3), 113);
  ASSERT_EQ(result.at(6).at(4), 93);
  ASSERT_EQ(result.at(6).at(5), 102);
  ASSERT_EQ(result.at(6).at(6), 0);
  ASSERT_EQ(result.at(6).at(7), 74);
  ASSERT_EQ(result.at(7).at(0), 149);
  ASSERT_EQ(result.at(7).at(1), 60);
  ASSERT_EQ(result.at(7).at(2), 131);
  ASSERT_EQ(result.at(7).at(3), 39);
  ASSERT_EQ(result.at(7).at(4), 63);
  ASSERT_EQ(result.at(7).at(5), 28);
  ASSERT_EQ(result.at(7).at(6), 91);
  ASSERT_EQ(result.at(7).at(7), 0);
}


// Test on graph in the file mediumEWD.txt
TEST(johnsonTest, mediumJohnson) {
  Graph<double> G {"mediumEWD.txt"};
  std::vector<std::vector<double> > result = johnsonAPSP(G);
  ASSERT_EQ(result.at(0).at(0), 0);
  ASSERT_EQ(result.at(0).at(1), 71178);
  ASSERT_EQ(result.at(0).at(2), 65237);
  ASSERT_EQ(result.at(0).at(3), 46262);
  ASSERT_EQ(result.at(0).at(4), 42292);
  ASSERT_EQ(result.at(0).at(5), 31728);
  ASSERT_EQ(result.at(0).at(6), 87118);
  ASSERT_EQ(result.at(0).at(7), 43429);
  ASSERT_EQ(result.at(0).at(8), 70729);
  ASSERT_EQ(result.at(0).at(9), 20353);
  ASSERT_EQ(result.at(0).at(10), 90191);
  ASSERT_EQ(result.at(0).at(11), 73843);
  ASSERT_EQ(result.at(0).at(12), 72201);
  ASSERT_EQ(result.at(0).at(13), 84085);
  ASSERT_EQ(result.at(0).at(14), 72024);
  ASSERT_EQ(result.at(0).at(15), 5719);
  ASSERT_EQ(result.at(0).at(16), 91595);
  ASSERT_EQ(result.at(0).at(17), 80442);
  ASSERT_EQ(result.at(0).at(18), 65692);
  ASSERT_EQ(result.at(0).at(19), 77088);
  ASSERT_EQ(result.at(0).at(20), 80934);
  ASSERT_EQ(result.at(0).at(21), 44164);
  ASSERT_EQ(result.at(0).at(22), 98366);
  ASSERT_EQ(result.at(0).at(23), 17027);
  ASSERT_EQ(result.at(0).at(24), 10191);
  ASSERT_EQ(result.at(0).at(25), 78521);
  ASSERT_EQ(result.at(0).at(26), 34115);
  ASSERT_EQ(result.at(0).at(27), 43330);
  ASSERT_EQ(result.at(0).at(28), 66247);
  ASSERT_EQ(result.at(0).at(29), 104001);
  ASSERT_EQ(result.at(0).at(30), 68439);
  ASSERT_EQ(result.at(0).at(31), 57558);
  ASSERT_EQ(result.at(0).at(32), 20674);
  ASSERT_EQ(result.at(0).at(33), 18203);
  ASSERT_EQ(result.at(0).at(34), 100977);
  ASSERT_EQ(result.at(0).at(35), 70615);
  ASSERT_EQ(result.at(0).at(36), 79419);
  ASSERT_EQ(result.at(0).at(37), 48071);
  ASSERT_EQ(result.at(0).at(38), 115099);
  ASSERT_EQ(result.at(0).at(39), 14770);
  ASSERT_EQ(result.at(0).at(40), 82209);
  ASSERT_EQ(result.at(0).at(41), 77104);
  ASSERT_EQ(result.at(0).at(42), 53781);
  ASSERT_EQ(result.at(0).at(43), 64780);
  ASSERT_EQ(result.at(0).at(44), 6471);
  ASSERT_EQ(result.at(0).at(45), 34360);
  ASSERT_EQ(result.at(0).at(46), 92473);
  ASSERT_EQ(result.at(0).at(47), 105899);
  ASSERT_EQ(result.at(0).at(48), 24841);
  ASSERT_EQ(result.at(0).at(49), 4849);
  ASSERT_EQ(result.at(0).at(50), 14883);
  ASSERT_EQ(result.at(0).at(51), 69018);
  ASSERT_EQ(result.at(0).at(52), 23237);
  ASSERT_EQ(result.at(0).at(53), 92925);
  ASSERT_EQ(result.at(0).at(54), 98295);
  ASSERT_EQ(result.at(0).at(55), 42859);
  ASSERT_EQ(result.at(0).at(56), 96929);
  ASSERT_EQ(result.at(0).at(57), 36678);
  ASSERT_EQ(result.at(0).at(58), 9955);
  ASSERT_EQ(result.at(0).at(59), 10657);
  ASSERT_EQ(result.at(0).at(60), 81924);
  ASSERT_EQ(result.at(0).at(61), 66473);
  ASSERT_EQ(result.at(0).at(62), 44994);
  ASSERT_EQ(result.at(0).at(63), 87887);
  ASSERT_EQ(result.at(0).at(64), 99850);
  ASSERT_EQ(result.at(0).at(65), 36353);
  ASSERT_EQ(result.at(0).at(66), 17220);
  ASSERT_EQ(result.at(0).at(67), 37012);
  ASSERT_EQ(result.at(0).at(68), 11816);
  ASSERT_EQ(result.at(0).at(69), 58054);
  ASSERT_EQ(result.at(0).at(70), 70216);
  ASSERT_EQ(result.at(0).at(71), 43196);
  ASSERT_EQ(result.at(0).at(72), 75425);
  ASSERT_EQ(result.at(0).at(73), 102677);
  ASSERT_EQ(result.at(0).at(74), 114938);
  ASSERT_EQ(result.at(0).at(75), 78037);
  ASSERT_EQ(result.at(0).at(76), 45397);
  ASSERT_EQ(result.at(0).at(77), 31559);
  ASSERT_EQ(result.at(0).at(78), 42525);
  ASSERT_EQ(result.at(0).at(79), 68775);
  ASSERT_EQ(result.at(0).at(80), 6821);
  ASSERT_EQ(result.at(0).at(81), 86205);
  ASSERT_EQ(result.at(0).at(82), 73032);
  ASSERT_EQ(result.at(0).at(83), 31454);
  ASSERT_EQ(result.at(0).at(84), 78345);
  ASSERT_EQ(result.at(0).at(85), 77376);
  ASSERT_EQ(result.at(0).at(86), 62879);
  ASSERT_EQ(result.at(0).at(87), 63312);
  ASSERT_EQ(result.at(0).at(88), 77757);
  ASSERT_EQ(result.at(0).at(89), 75523);
  ASSERT_EQ(result.at(0).at(90), 53629);
  ASSERT_EQ(result.at(0).at(91), 100157);
  ASSERT_EQ(result.at(0).at(92), 40508);
  ASSERT_EQ(result.at(0).at(93), 13264);
  ASSERT_EQ(result.at(0).at(94), 67297);
  ASSERT_EQ(result.at(0).at(95), 43436);
  ASSERT_EQ(result.at(0).at(96), 87400);
  ASSERT_EQ(result.at(0).at(97), 7705);
  ASSERT_EQ(result.at(0).at(98), 86324);
  ASSERT_EQ(result.at(0).at(99), 93526);
  ASSERT_EQ(result.at(0).at(100), 81525);
  ASSERT_EQ(result.at(0).at(101), 52599);
  ASSERT_EQ(result.at(0).at(102), 29925);
  ASSERT_EQ(result.at(0).at(103), 82407);
  ASSERT_EQ(result.at(0).at(104), 22943);
  ASSERT_EQ(result.at(0).at(105), 84396);
  ASSERT_EQ(result.at(0).at(106), 84000);
  ASSERT_EQ(result.at(0).at(107), 63694);
  ASSERT_EQ(result.at(0).at(108), 58041);
  ASSERT_EQ(result.at(0).at(109), 105854);
  ASSERT_EQ(result.at(0).at(110), 57527);
  ASSERT_EQ(result.at(0).at(111), 73212);
  ASSERT_EQ(result.at(0).at(112), 47449);
  ASSERT_EQ(result.at(0).at(113), 65457);
  ASSERT_EQ(result.at(0).at(114), 9610);
  ASSERT_EQ(result.at(0).at(115), 49055);
  ASSERT_EQ(result.at(0).at(116), 80168);
  ASSERT_EQ(result.at(0).at(117), 95939);
  ASSERT_EQ(result.at(0).at(118), 29543);
  ASSERT_EQ(result.at(0).at(119), 94430);
  ASSERT_EQ(result.at(0).at(120), 91902);
  ASSERT_EQ(result.at(0).at(121), 70836);
  ASSERT_EQ(result.at(0).at(122), 52210);
  ASSERT_EQ(result.at(0).at(123), 90520);
  ASSERT_EQ(result.at(0).at(124), 26561);
  ASSERT_EQ(result.at(0).at(125), 42078);
  ASSERT_EQ(result.at(0).at(126), 115186);
  ASSERT_EQ(result.at(0).at(127), 86205);
  ASSERT_EQ(result.at(0).at(128), 46158);
  ASSERT_EQ(result.at(0).at(129), 82761);
  ASSERT_EQ(result.at(0).at(130), 66996);
  ASSERT_EQ(result.at(0).at(131), 83033);
  ASSERT_EQ(result.at(0).at(132), 40529);
  ASSERT_EQ(result.at(0).at(133), 78673);
  ASSERT_EQ(result.at(0).at(134), 90613);
  ASSERT_EQ(result.at(0).at(135), 52534);
  ASSERT_EQ(result.at(0).at(136), 52970);
  ASSERT_EQ(result.at(0).at(137), 100909);
  ASSERT_EQ(result.at(0).at(138), 35516);
  ASSERT_EQ(result.at(0).at(139), 51857);
  ASSERT_EQ(result.at(0).at(140), 96004);
  ASSERT_EQ(result.at(0).at(141), 59227);
  ASSERT_EQ(result.at(0).at(142), 26181);
  ASSERT_EQ(result.at(0).at(143), 77179);
  ASSERT_EQ(result.at(0).at(144), 14843);
  ASSERT_EQ(result.at(0).at(145), 99791);
  ASSERT_EQ(result.at(0).at(146), 97434);
  ASSERT_EQ(result.at(0).at(147), 91292);
  ASSERT_EQ(result.at(0).at(148), 42770);
  ASSERT_EQ(result.at(0).at(149), 9659);
  ASSERT_EQ(result.at(0).at(150), 82086);
  ASSERT_EQ(result.at(0).at(151), 26057);
  ASSERT_EQ(result.at(0).at(152), 70704);
  ASSERT_EQ(result.at(0).at(153), 49792);
  ASSERT_EQ(result.at(0).at(154), 35954);
  ASSERT_EQ(result.at(0).at(155), 26385);
  ASSERT_EQ(result.at(0).at(156), 57224);
  ASSERT_EQ(result.at(0).at(157), 43451);
  ASSERT_EQ(result.at(0).at(158), 75299);
  ASSERT_EQ(result.at(0).at(159), 51273);
  ASSERT_EQ(result.at(0).at(160), 11714);
  ASSERT_EQ(result.at(0).at(161), 89158);
  ASSERT_EQ(result.at(0).at(162), 95645);
  ASSERT_EQ(result.at(0).at(163), 9368);
  ASSERT_EQ(result.at(0).at(164), 78814);
  ASSERT_EQ(result.at(0).at(165), 21384);
  ASSERT_EQ(result.at(0).at(166), 80120);
  ASSERT_EQ(result.at(0).at(167), 105366);
  ASSERT_EQ(result.at(0).at(168), 16904);
  ASSERT_EQ(result.at(0).at(169), 88030);
  ASSERT_EQ(result.at(0).at(170), 72686);
  ASSERT_EQ(result.at(0).at(171), 33069);
  ASSERT_EQ(result.at(0).at(172), 31212);
  ASSERT_EQ(result.at(0).at(173), 58106);
  ASSERT_EQ(result.at(0).at(174), 78790);
  ASSERT_EQ(result.at(0).at(175), 83778);
  ASSERT_EQ(result.at(0).at(176), 8927);
  ASSERT_EQ(result.at(0).at(177), 81859);
  ASSERT_EQ(result.at(0).at(178), 91507);
  ASSERT_EQ(result.at(0).at(179), 74538);
  ASSERT_EQ(result.at(0).at(180), 23140);
  ASSERT_EQ(result.at(0).at(181), 47275);
  ASSERT_EQ(result.at(0).at(182), 73993);
  ASSERT_EQ(result.at(0).at(183), 109748);
  ASSERT_EQ(result.at(0).at(184), 40046);
  ASSERT_EQ(result.at(0).at(185), 15484);
  ASSERT_EQ(result.at(0).at(186), 84784);
  ASSERT_EQ(result.at(0).at(187), 20904);
  ASSERT_EQ(result.at(0).at(188), 39905);
  ASSERT_EQ(result.at(0).at(189), 79131);
  ASSERT_EQ(result.at(0).at(190), 85857);
  ASSERT_EQ(result.at(0).at(191), 10711);
  ASSERT_EQ(result.at(0).at(192), 86173);
  ASSERT_EQ(result.at(0).at(193), 84044);
  ASSERT_EQ(result.at(0).at(194), 71585);
  ASSERT_EQ(result.at(0).at(195), 24938);
  ASSERT_EQ(result.at(0).at(196), 55714);
  ASSERT_EQ(result.at(0).at(197), 36445);
  ASSERT_EQ(result.at(0).at(198), 66781);
  ASSERT_EQ(result.at(0).at(199), 86300);
  ASSERT_EQ(result.at(0).at(200), 71948);
  ASSERT_EQ(result.at(0).at(201), 22440);
  ASSERT_EQ(result.at(0).at(202), 4678);
  ASSERT_EQ(result.at(0).at(203), 75288);
  ASSERT_EQ(result.at(0).at(204), 5476);
  ASSERT_EQ(result.at(0).at(205), 52136);
  ASSERT_EQ(result.at(0).at(206), 19793);
  ASSERT_EQ(result.at(0).at(207), 62742);
  ASSERT_EQ(result.at(0).at(208), 25676);
  ASSERT_EQ(result.at(0).at(209), 9511);
  ASSERT_EQ(result.at(0).at(210), 60068);
  ASSERT_EQ(result.at(0).at(211), 8438);
  ASSERT_EQ(result.at(0).at(212), 65392);
  ASSERT_EQ(result.at(0).at(213), 26387);
  ASSERT_EQ(result.at(0).at(214), 59481);
  ASSERT_EQ(result.at(0).at(215), 111024);
  ASSERT_EQ(result.at(0).at(216), 32372);
  ASSERT_EQ(result.at(0).at(217), 32511);
  ASSERT_EQ(result.at(0).at(218), 102985);
  ASSERT_EQ(result.at(0).at(219), 57293);
  ASSERT_EQ(result.at(0).at(220), 81606);
  ASSERT_EQ(result.at(0).at(221), 61015);
  ASSERT_EQ(result.at(0).at(222), 7573);
  ASSERT_EQ(result.at(0).at(223), 70122);
  ASSERT_EQ(result.at(0).at(224), 107478);
  ASSERT_EQ(result.at(0).at(225), 2383);
  ASSERT_EQ(result.at(0).at(226), 24452);
  ASSERT_EQ(result.at(0).at(227), 93793);
  ASSERT_EQ(result.at(0).at(228), 53897);
  ASSERT_EQ(result.at(0).at(229), 83610);
  ASSERT_EQ(result.at(0).at(230), 41360);
  ASSERT_EQ(result.at(0).at(231), 16855);
  ASSERT_EQ(result.at(0).at(232), 25743);
  ASSERT_EQ(result.at(0).at(233), 45114);
  ASSERT_EQ(result.at(0).at(234), 55496);
  ASSERT_EQ(result.at(0).at(235), 35692);
  ASSERT_EQ(result.at(0).at(236), 89599);
  ASSERT_EQ(result.at(0).at(237), 92849);
  ASSERT_EQ(result.at(0).at(238), 36328);
  ASSERT_EQ(result.at(0).at(239), 44590);
  ASSERT_EQ(result.at(0).at(240), 37592);
  ASSERT_EQ(result.at(0).at(241), 53286);
  ASSERT_EQ(result.at(0).at(242), 63744);
  ASSERT_EQ(result.at(0).at(243), 85566);
  ASSERT_EQ(result.at(0).at(244), 69713);
  ASSERT_EQ(result.at(0).at(245), 34224);
  ASSERT_EQ(result.at(0).at(246), 80214);
  ASSERT_EQ(result.at(0).at(247), 87806);
  ASSERT_EQ(result.at(0).at(248), 16303);
  ASSERT_EQ(result.at(0).at(249), 78283);
  ASSERT_EQ(result.at(1).at(0), 71178);
  ASSERT_EQ(result.at(1).at(1), 0);
  ASSERT_EQ(result.at(1).at(2), 56820);
  ASSERT_EQ(result.at(1).at(3), 50917);
  ASSERT_EQ(result.at(1).at(4), 29771);
  ASSERT_EQ(result.at(1).at(5), 40955);
  ASSERT_EQ(result.at(1).at(6), 55449);
  ASSERT_EQ(result.at(1).at(7), 52742);
  ASSERT_EQ(result.at(1).at(8), 82827);
  ASSERT_EQ(result.at(1).at(9), 81010);
  ASSERT_EQ(result.at(1).at(10), 99361);
  ASSERT_EQ(result.at(1).at(11), 87531);
  ASSERT_EQ(result.at(1).at(12), 32606);
  ASSERT_EQ(result.at(1).at(13), 69636);
  ASSERT_EQ(result.at(1).at(14), 56730);
  ASSERT_EQ(result.at(1).at(15), 76896);
  ASSERT_EQ(result.at(1).at(16), 57179);
  ASSERT_EQ(result.at(1).at(17), 27288);
  ASSERT_EQ(result.at(1).at(18), 49395);
  ASSERT_EQ(result.at(1).at(19), 75902);
  ASSERT_EQ(result.at(1).at(20), 17127);
  ASSERT_EQ(result.at(1).at(21), 38359);
  ASSERT_EQ(result.at(1).at(22), 33598);
  ASSERT_EQ(result.at(1).at(23), 78843);
  ASSERT_EQ(result.at(1).at(24), 81369);
  ASSERT_EQ(result.at(1).at(25), 29118);
  ASSERT_EQ(result.at(1).at(26), 38066);
  ASSERT_EQ(result.at(1).at(27), 38833);
  ASSERT_EQ(result.at(1).at(28), 29211);
  ASSERT_EQ(result.at(1).at(29), 50775);
  ASSERT_EQ(result.at(1).at(30), 78842);
  ASSERT_EQ(result.at(1).at(31), 65812);
  ASSERT_EQ(result.at(1).at(32), 51085);
  ASSERT_EQ(result.at(1).at(33), 83586);
  ASSERT_EQ(result.at(1).at(34), 36489);
  ASSERT_EQ(result.at(1).at(35), 38685);
  ASSERT_EQ(result.at(1).at(36), 37041);
  ASSERT_EQ(result.at(1).at(37), 59429);
  ASSERT_EQ(result.at(1).at(38), 61898);
  ASSERT_EQ(result.at(1).at(39), 83624);
  ASSERT_EQ(result.at(1).at(40), 13578);
  ASSERT_EQ(result.at(1).at(41), 31864);
  ASSERT_EQ(result.at(1).at(42), 54573);
  ASSERT_EQ(result.at(1).at(43), 79195);
  ASSERT_EQ(result.at(1).at(44), 64707);
  ASSERT_EQ(result.at(1).at(45), 47417);
  ASSERT_EQ(result.at(1).at(46), 23046);
  ASSERT_EQ(result.at(1).at(47), 52673);
  ASSERT_EQ(result.at(1).at(48), 57312);
  ASSERT_EQ(result.at(1).at(49), 66653);
  ASSERT_EQ(result.at(1).at(50), 64313);
}

TEST(johnsonTest, random30) {
  randomTest(johnsonAPSP<int>, 30, 3'289'211);
}

TEST(johnsonTest, random100a) {
  randomTest(johnsonAPSP<int>, 100, 987'129, 0.1);
}

TEST(johnsonTest, random100b) {
  randomTest(johnsonAPSP<int>, 100, 2'382'398, 0.1);
}

TEST(johnsonTest, random100c) {
  randomTest(johnsonAPSP<int>, 100, 329'823, 0.1);
}

TEST(johnsonTest, random100d) {
  randomTest(johnsonAPSP<int>, 100, 389'211, 0.1);
}

TEST(johnsonTest, random100e) {
  randomTest(johnsonAPSP<int>, 100, 10'918'101, 0.1);
}

TEST(johnsonTest, random150_0) {
  randomTest(johnsonAPSP<int>, 150, 1812528527, 0.1);
}

TEST(johnsonTest, random150_1) {
  randomTest(johnsonAPSP<int>, 150, 1915543033, 0.1);
}

TEST(johnsonTest, random150_2) {
  randomTest(johnsonAPSP<int>, 150, 1728632203, 0.1);
}

TEST(johnsonTest, random150_3) {
  randomTest(johnsonAPSP<int>, 150, 3400310732, 0.1);
}

TEST(johnsonTest, random150_4) {
  randomTest(johnsonAPSP<int>, 150, 1509646079, 0.1);
}

TEST(johnsonTest, random150_5) {
  randomTest(johnsonAPSP<int>, 150, 877500092, 0.1);
}

TEST(johnsonTest, random150_6) {
  randomTest(johnsonAPSP<int>, 150, 307693294, 0.1);
}

TEST(johnsonTest, random150_7) {
  randomTest(johnsonAPSP<int>, 150, 1740768105, 0.1);
}

TEST(johnsonTest, random150_8) {
  randomTest(johnsonAPSP<int>, 150, 2445968672, 0.1);
}

TEST(johnsonTest, random150_9) {
  randomTest(johnsonAPSP<int>, 150, 55419902, 0.1);
}

TEST(johnsonTest, random150_10) {
  randomTest(johnsonAPSP<int>, 150, 449020997, 0.1);
}

TEST(johnsonTest, random150_11) {
  randomTest(johnsonAPSP<int>, 150, 1922416908, 0.1);
}

TEST(johnsonTest, random200a) {
  randomTest(johnsonAPSP<int>, 200, 1'238'923, 0.1);
}

TEST(johnsonTest, random200b) {
  randomTest(johnsonAPSP<int>, 200, 108'181'745, 0.1);
}

TEST(johnsonTest, random200c) {
  randomTest(johnsonAPSP<int>, 200, 3'281'918, 0.1);
}

TEST(johnsonTest, random200_1) {
  randomTest(johnsonAPSP<int>, 200, 445192270, 0.1);
}

TEST(johnsonTest, random200_2) {
  randomTest(johnsonAPSP<int>, 200, 1609628408, 0.1);
}

TEST(johnsonTest, random200_3) {
  randomTest(johnsonAPSP<int>, 200, 1801774926, 0.1);
}

TEST(johnsonTest, random200_4) {
  randomTest(johnsonAPSP<int>, 200, 2002738901, 0.1);
}

TEST(johnsonTest, random200_5) {
  randomTest(johnsonAPSP<int>, 200, 1138962375, 0.1);
}

TEST(johnsonTest, random200_6) {
  randomTest(johnsonAPSP<int>, 200, 1326754919, 0.1);
}

TEST(johnsonTest, random200_7) {
  randomTest(johnsonAPSP<int>, 200, 2178494701, 0.1);
}

TEST(johnsonTest, random200_8) {
  randomTest(johnsonAPSP<int>, 200, 895276633, 0.1);
}

TEST(johnsonTest, random200_9) {
  randomTest(johnsonAPSP<int>, 200, 419725955, 0.1);
}

TEST(johnsonTest, random200_10) {
  randomTest(johnsonAPSP<int>, 200, 2031438742, 0.1);
}

TEST(johnsonTest, random200_11) {
  randomTest(johnsonAPSP<int>, 200, 3188665175, 0.1);
}

TEST(johnsonTest, random200_12) {
  randomTest(johnsonAPSP<int>, 200, 4014034417, 0.1);
}

TEST(johnsonTest, random200_13) {
  randomTest(johnsonAPSP<int>, 200, 2553198683, 0.1);
}

TEST(johnsonTest, random200_14) {
  randomTest(johnsonAPSP<int>, 200, 1186300218, 0.1);
}

TEST(johnsonTest, random250a) {
  randomTest(johnsonAPSP<int>, 200, 23'892'389, 0.1);
}

TEST(johnsonTest, random250b) {
  randomTest(johnsonAPSP<int>, 200, 34'894'398, 0.1);
}

TEST(johnsonTest, random250c) {
  randomTest(johnsonAPSP<int>, 200, 890'119'109, 0.1);
}

TEST(johnsonTest, random300a) {
  randomTest(johnsonAPSP<int>, 300, 98'982, 0.05);
}

TEST(johnsonTest, random300b) {
  randomTest(johnsonAPSP<int>, 300, 32'892, 0.05);
}

TEST(johnsonTest, random300c) {
  randomTest(johnsonAPSP<int>, 300, 238'923, 0.05);
}

TEST(johnsonTest, random400) {
  randomTest(johnsonAPSP<int>, 400, 89'239'823, 0.01);
}

TEST(johnsonTest, random500) {
  randomTest(johnsonAPSP<int>, 500, 2'829'211, 0.01);
}

TEST(johnsonTest, random1000a) {
  randomTest(johnsonAPSP<int>, 1000, 32'218'119, 0.005);
}

TEST(johnsonTest, random1000b) {
  randomTest(johnsonAPSP<int>, 1000, 98'239'283, 0.005);
}

TEST(johnsonTest, random1000c) {
  randomTest(johnsonAPSP<int>, 1'000, 2'389'239, 0.005);
}

// *** End of tests of johnsonAPSP

// *** Tests of floydWarshallAPSP
TEST(FWTest, dijkstraFailure) {
  Graph<int> G {4};
  G.addEdge(0, 1, 1);
  G.addEdge(0, 2, 4);
  G.addEdge(2, 1, -4);
  G.addEdge(1, 3, 1);
  ASSERT_FALSE(existsNegativeCycle(G));
  std::vector<std::vector<int> > distanceMatrix = floydWarshallAPSP<int>(G);
  for (int v = 0; v < G.size(); ++v) {
    ASSERT_TRUE(allEdgesRelaxed(distanceMatrix.at(v), G, v));
  }
}

TEST(FWTest, dijkstraFailure2) {
  Graph<int> G {5};
  G.addEdge(0, 1, 1);
  G.addEdge(0, 2, 4);
  G.addEdge(2, 1, -4);
  G.addEdge(1, 3, 1);
  G.addEdge(3, 4, 1);
  ASSERT_FALSE(existsNegativeCycle(G));
  std::vector<std::vector<int> > distanceMatrix = floydWarshallAPSP<int>(G);
  for (int v = 0; v < G.size(); ++v) {
    ASSERT_TRUE(allEdgesRelaxed(distanceMatrix.at(v), G, v));
  }
}

// Test of FloydWarshall on graph in the file tinyEWD.txt
TEST(FWTest, tinyFloydWarshall) {
  Graph<double> G {"tinyEWD.txt"};
  std::vector<std::vector<double> > result = floydWarshallAPSP(G);
  ASSERT_EQ(result.at(0).at(0), 0);
  ASSERT_EQ(result.at(0).at(1), 105);
  ASSERT_EQ(result.at(0).at(2), 26);
  ASSERT_EQ(result.at(0).at(3), 99);
  ASSERT_EQ(result.at(0).at(4), 38);
  ASSERT_EQ(result.at(0).at(5), 73);
  ASSERT_EQ(result.at(0).at(6), 151);
  ASSERT_EQ(result.at(0).at(7), 60);
  ASSERT_EQ(result.at(1).at(0), 139);
  ASSERT_EQ(result.at(1).at(1), 0);
  ASSERT_EQ(result.at(1).at(2), 121);
  ASSERT_EQ(result.at(1).at(3), 29);
  ASSERT_EQ(result.at(1).at(4), 174);
  ASSERT_EQ(result.at(1).at(5), 183);
  ASSERT_EQ(result.at(1).at(6), 81);
  ASSERT_EQ(result.at(1).at(7), 155);
  ASSERT_EQ(result.at(2).at(0), 183);
  ASSERT_EQ(result.at(2).at(1), 94);
  ASSERT_EQ(result.at(2).at(2), 0);
  ASSERT_EQ(result.at(2).at(3), 73);
  ASSERT_EQ(result.at(2).at(4), 97);
  ASSERT_EQ(result.at(2).at(5), 62);
  ASSERT_EQ(result.at(2).at(6), 125);
  ASSERT_EQ(result.at(2).at(7), 34);
  ASSERT_EQ(result.at(3).at(0), 110);
  ASSERT_EQ(result.at(3).at(1), 186);
  ASSERT_EQ(result.at(3).at(2), 92);
  ASSERT_EQ(result.at(3).at(3), 0);
  ASSERT_EQ(result.at(3).at(4), 145);
  ASSERT_EQ(result.at(3).at(5), 154);
  ASSERT_EQ(result.at(3).at(6), 52);
  ASSERT_EQ(result.at(3).at(7), 126);
  ASSERT_EQ(result.at(4).at(0), 186);
  ASSERT_EQ(result.at(4).at(1), 67);
  ASSERT_EQ(result.at(4).at(2), 168);
  ASSERT_EQ(result.at(4).at(3), 76);
  ASSERT_EQ(result.at(4).at(4), 0);
  ASSERT_EQ(result.at(4).at(5), 35);
  ASSERT_EQ(result.at(4).at(6), 128);
  ASSERT_EQ(result.at(4).at(7), 37);
  ASSERT_EQ(result.at(5).at(0), 171);
  ASSERT_EQ(result.at(5).at(1), 32);
  ASSERT_EQ(result.at(5).at(2), 153);
  ASSERT_EQ(result.at(5).at(3), 61);
  ASSERT_EQ(result.at(5).at(4), 35);
  ASSERT_EQ(result.at(5).at(5), 0);
  ASSERT_EQ(result.at(5).at(6), 113);
  ASSERT_EQ(result.at(5).at(7), 28);
  ASSERT_EQ(result.at(6).at(0), 58);
  ASSERT_EQ(result.at(6).at(1), 134);
  ASSERT_EQ(result.at(6).at(2), 40);
  ASSERT_EQ(result.at(6).at(3), 113);
  ASSERT_EQ(result.at(6).at(4), 93);
  ASSERT_EQ(result.at(6).at(5), 102);
  ASSERT_EQ(result.at(6).at(6), 0);
  ASSERT_EQ(result.at(6).at(7), 74);
  ASSERT_EQ(result.at(7).at(0), 149);
  ASSERT_EQ(result.at(7).at(1), 60);
  ASSERT_EQ(result.at(7).at(2), 131);
  ASSERT_EQ(result.at(7).at(3), 39);
  ASSERT_EQ(result.at(7).at(4), 63);
  ASSERT_EQ(result.at(7).at(5), 28);
  ASSERT_EQ(result.at(7).at(6), 91);
  ASSERT_EQ(result.at(7).at(7), 0);
}

// Test of FloydWarshall on graph in the file mediumEWD.txt
TEST(FWTest, mediumFloydWarshall) {
  Graph<double> G {"mediumEWD.txt"};
  std::vector<std::vector<double> > result = floydWarshallAPSP(G);
  ASSERT_EQ(result.at(0).at(0), 0);
  ASSERT_EQ(result.at(0).at(1), 71178);
  ASSERT_EQ(result.at(0).at(2), 65237);
  ASSERT_EQ(result.at(0).at(3), 46262);
  ASSERT_EQ(result.at(0).at(4), 42292);
  ASSERT_EQ(result.at(0).at(5), 31728);
  ASSERT_EQ(result.at(0).at(6), 87118);
  ASSERT_EQ(result.at(0).at(7), 43429);
  ASSERT_EQ(result.at(0).at(8), 70729);
  ASSERT_EQ(result.at(0).at(9), 20353);
  ASSERT_EQ(result.at(0).at(10), 90191);
  ASSERT_EQ(result.at(0).at(11), 73843);
  ASSERT_EQ(result.at(0).at(12), 72201);
  ASSERT_EQ(result.at(0).at(13), 84085);
  ASSERT_EQ(result.at(0).at(14), 72024);
  ASSERT_EQ(result.at(0).at(15), 5719);
  ASSERT_EQ(result.at(0).at(16), 91595);
  ASSERT_EQ(result.at(0).at(17), 80442);
  ASSERT_EQ(result.at(0).at(18), 65692);
  ASSERT_EQ(result.at(0).at(19), 77088);
  ASSERT_EQ(result.at(0).at(20), 80934);
  ASSERT_EQ(result.at(0).at(21), 44164);
  ASSERT_EQ(result.at(0).at(22), 98366);
  ASSERT_EQ(result.at(0).at(23), 17027);
  ASSERT_EQ(result.at(0).at(24), 10191);
  ASSERT_EQ(result.at(0).at(25), 78521);
  ASSERT_EQ(result.at(0).at(26), 34115);
  ASSERT_EQ(result.at(0).at(27), 43330);
  ASSERT_EQ(result.at(0).at(28), 66247);
  ASSERT_EQ(result.at(0).at(29), 104001);
  ASSERT_EQ(result.at(0).at(30), 68439);
  ASSERT_EQ(result.at(0).at(31), 57558);
  ASSERT_EQ(result.at(0).at(32), 20674);
  ASSERT_EQ(result.at(0).at(33), 18203);
  ASSERT_EQ(result.at(0).at(34), 100977);
  ASSERT_EQ(result.at(0).at(35), 70615);
  ASSERT_EQ(result.at(0).at(36), 79419);
  ASSERT_EQ(result.at(0).at(37), 48071);
  ASSERT_EQ(result.at(0).at(38), 115099);
  ASSERT_EQ(result.at(0).at(39), 14770);
  ASSERT_EQ(result.at(0).at(40), 82209);
  ASSERT_EQ(result.at(0).at(41), 77104);
  ASSERT_EQ(result.at(0).at(42), 53781);
  ASSERT_EQ(result.at(0).at(43), 64780);
  ASSERT_EQ(result.at(0).at(44), 6471);
  ASSERT_EQ(result.at(0).at(45), 34360);
  ASSERT_EQ(result.at(0).at(46), 92473);
  ASSERT_EQ(result.at(0).at(47), 105899);
  ASSERT_EQ(result.at(0).at(48), 24841);
  ASSERT_EQ(result.at(0).at(49), 4849);
  ASSERT_EQ(result.at(0).at(50), 14883);
  ASSERT_EQ(result.at(0).at(51), 69018);
  ASSERT_EQ(result.at(0).at(52), 23237);
  ASSERT_EQ(result.at(0).at(53), 92925);
  ASSERT_EQ(result.at(0).at(54), 98295);
  ASSERT_EQ(result.at(0).at(55), 42859);
  ASSERT_EQ(result.at(0).at(56), 96929);
  ASSERT_EQ(result.at(0).at(57), 36678);
  ASSERT_EQ(result.at(0).at(58), 9955);
  ASSERT_EQ(result.at(0).at(59), 10657);
  ASSERT_EQ(result.at(0).at(60), 81924);
  ASSERT_EQ(result.at(0).at(61), 66473);
  ASSERT_EQ(result.at(0).at(62), 44994);
  ASSERT_EQ(result.at(0).at(63), 87887);
  ASSERT_EQ(result.at(0).at(64), 99850);
  ASSERT_EQ(result.at(0).at(65), 36353);
  ASSERT_EQ(result.at(0).at(66), 17220);
  ASSERT_EQ(result.at(0).at(67), 37012);
  ASSERT_EQ(result.at(0).at(68), 11816);
  ASSERT_EQ(result.at(0).at(69), 58054);
  ASSERT_EQ(result.at(0).at(70), 70216);
  ASSERT_EQ(result.at(0).at(71), 43196);
  ASSERT_EQ(result.at(0).at(72), 75425);
  ASSERT_EQ(result.at(0).at(73), 102677);
  ASSERT_EQ(result.at(0).at(74), 114938);
  ASSERT_EQ(result.at(0).at(75), 78037);
  ASSERT_EQ(result.at(0).at(76), 45397);
  ASSERT_EQ(result.at(0).at(77), 31559);
  ASSERT_EQ(result.at(0).at(78), 42525);
  ASSERT_EQ(result.at(0).at(79), 68775);
  ASSERT_EQ(result.at(0).at(80), 6821);
  ASSERT_EQ(result.at(0).at(81), 86205);
  ASSERT_EQ(result.at(0).at(82), 73032);
  ASSERT_EQ(result.at(0).at(83), 31454);
  ASSERT_EQ(result.at(0).at(84), 78345);
  ASSERT_EQ(result.at(0).at(85), 77376);
  ASSERT_EQ(result.at(0).at(86), 62879);
  ASSERT_EQ(result.at(0).at(87), 63312);
  ASSERT_EQ(result.at(0).at(88), 77757);
  ASSERT_EQ(result.at(0).at(89), 75523);
  ASSERT_EQ(result.at(0).at(90), 53629);
  ASSERT_EQ(result.at(0).at(91), 100157);
  ASSERT_EQ(result.at(0).at(92), 40508);
  ASSERT_EQ(result.at(0).at(93), 13264);
  ASSERT_EQ(result.at(0).at(94), 67297);
  ASSERT_EQ(result.at(0).at(95), 43436);
  ASSERT_EQ(result.at(0).at(96), 87400);
  ASSERT_EQ(result.at(0).at(97), 7705);
  ASSERT_EQ(result.at(0).at(98), 86324);
  ASSERT_EQ(result.at(0).at(99), 93526);
  ASSERT_EQ(result.at(0).at(100), 81525);
  ASSERT_EQ(result.at(0).at(101), 52599);
  ASSERT_EQ(result.at(0).at(102), 29925);
  ASSERT_EQ(result.at(0).at(103), 82407);
  ASSERT_EQ(result.at(0).at(104), 22943);
  ASSERT_EQ(result.at(0).at(105), 84396);
  ASSERT_EQ(result.at(0).at(106), 84000);
  ASSERT_EQ(result.at(0).at(107), 63694);
  ASSERT_EQ(result.at(0).at(108), 58041);
  ASSERT_EQ(result.at(0).at(109), 105854);
  ASSERT_EQ(result.at(0).at(110), 57527);
  ASSERT_EQ(result.at(0).at(111), 73212);
  ASSERT_EQ(result.at(0).at(112), 47449);
  ASSERT_EQ(result.at(0).at(113), 65457);
  ASSERT_EQ(result.at(0).at(114), 9610);
  ASSERT_EQ(result.at(0).at(115), 49055);
  ASSERT_EQ(result.at(0).at(116), 80168);
  ASSERT_EQ(result.at(0).at(117), 95939);
  ASSERT_EQ(result.at(0).at(118), 29543);
  ASSERT_EQ(result.at(0).at(119), 94430);
  ASSERT_EQ(result.at(0).at(120), 91902);
  ASSERT_EQ(result.at(0).at(121), 70836);
  ASSERT_EQ(result.at(0).at(122), 52210);
  ASSERT_EQ(result.at(0).at(123), 90520);
  ASSERT_EQ(result.at(0).at(124), 26561);
  ASSERT_EQ(result.at(0).at(125), 42078);
  ASSERT_EQ(result.at(0).at(126), 115186);
  ASSERT_EQ(result.at(0).at(127), 86205);
  ASSERT_EQ(result.at(0).at(128), 46158);
  ASSERT_EQ(result.at(0).at(129), 82761);
  ASSERT_EQ(result.at(0).at(130), 66996);
  ASSERT_EQ(result.at(0).at(131), 83033);
  ASSERT_EQ(result.at(0).at(132), 40529);
  ASSERT_EQ(result.at(0).at(133), 78673);
  ASSERT_EQ(result.at(0).at(134), 90613);
  ASSERT_EQ(result.at(0).at(135), 52534);
  ASSERT_EQ(result.at(0).at(136), 52970);
  ASSERT_EQ(result.at(0).at(137), 100909);
  ASSERT_EQ(result.at(0).at(138), 35516);
  ASSERT_EQ(result.at(0).at(139), 51857);
  ASSERT_EQ(result.at(0).at(140), 96004);
  ASSERT_EQ(result.at(0).at(141), 59227);
  ASSERT_EQ(result.at(0).at(142), 26181);
  ASSERT_EQ(result.at(0).at(143), 77179);
  ASSERT_EQ(result.at(0).at(144), 14843);
  ASSERT_EQ(result.at(0).at(145), 99791);
  ASSERT_EQ(result.at(0).at(146), 97434);
  ASSERT_EQ(result.at(0).at(147), 91292);
  ASSERT_EQ(result.at(0).at(148), 42770);
  ASSERT_EQ(result.at(0).at(149), 9659);
  ASSERT_EQ(result.at(0).at(150), 82086);
  ASSERT_EQ(result.at(0).at(151), 26057);
  ASSERT_EQ(result.at(0).at(152), 70704);
  ASSERT_EQ(result.at(0).at(153), 49792);
  ASSERT_EQ(result.at(0).at(154), 35954);
  ASSERT_EQ(result.at(0).at(155), 26385);
  ASSERT_EQ(result.at(0).at(156), 57224);
  ASSERT_EQ(result.at(0).at(157), 43451);
  ASSERT_EQ(result.at(0).at(158), 75299);
  ASSERT_EQ(result.at(0).at(159), 51273);
  ASSERT_EQ(result.at(0).at(160), 11714);
  ASSERT_EQ(result.at(0).at(161), 89158);
  ASSERT_EQ(result.at(0).at(162), 95645);
  ASSERT_EQ(result.at(0).at(163), 9368);
  ASSERT_EQ(result.at(0).at(164), 78814);
  ASSERT_EQ(result.at(0).at(165), 21384);
  ASSERT_EQ(result.at(0).at(166), 80120);
  ASSERT_EQ(result.at(0).at(167), 105366);
  ASSERT_EQ(result.at(0).at(168), 16904);
  ASSERT_EQ(result.at(0).at(169), 88030);
  ASSERT_EQ(result.at(0).at(170), 72686);
  ASSERT_EQ(result.at(0).at(171), 33069);
  ASSERT_EQ(result.at(0).at(172), 31212);
  ASSERT_EQ(result.at(0).at(173), 58106);
  ASSERT_EQ(result.at(0).at(174), 78790);
  ASSERT_EQ(result.at(0).at(175), 83778);
  ASSERT_EQ(result.at(0).at(176), 8927);
  ASSERT_EQ(result.at(0).at(177), 81859);
  ASSERT_EQ(result.at(0).at(178), 91507);
  ASSERT_EQ(result.at(0).at(179), 74538);
  ASSERT_EQ(result.at(0).at(180), 23140);
  ASSERT_EQ(result.at(0).at(181), 47275);
  ASSERT_EQ(result.at(0).at(182), 73993);
  ASSERT_EQ(result.at(0).at(183), 109748);
  ASSERT_EQ(result.at(0).at(184), 40046);
  ASSERT_EQ(result.at(0).at(185), 15484);
  ASSERT_EQ(result.at(0).at(186), 84784);
  ASSERT_EQ(result.at(0).at(187), 20904);
  ASSERT_EQ(result.at(0).at(188), 39905);
  ASSERT_EQ(result.at(0).at(189), 79131);
  ASSERT_EQ(result.at(0).at(190), 85857);
  ASSERT_EQ(result.at(0).at(191), 10711);
  ASSERT_EQ(result.at(0).at(192), 86173);
  ASSERT_EQ(result.at(0).at(193), 84044);
  ASSERT_EQ(result.at(0).at(194), 71585);
  ASSERT_EQ(result.at(0).at(195), 24938);
  ASSERT_EQ(result.at(0).at(196), 55714);
  ASSERT_EQ(result.at(0).at(197), 36445);
  ASSERT_EQ(result.at(0).at(198), 66781);
  ASSERT_EQ(result.at(0).at(199), 86300);
  ASSERT_EQ(result.at(0).at(200), 71948);
  ASSERT_EQ(result.at(0).at(201), 22440);
  ASSERT_EQ(result.at(0).at(202), 4678);
  ASSERT_EQ(result.at(0).at(203), 75288);
  ASSERT_EQ(result.at(0).at(204), 5476);
  ASSERT_EQ(result.at(0).at(205), 52136);
  ASSERT_EQ(result.at(0).at(206), 19793);
  ASSERT_EQ(result.at(0).at(207), 62742);
  ASSERT_EQ(result.at(0).at(208), 25676);
  ASSERT_EQ(result.at(0).at(209), 9511);
  ASSERT_EQ(result.at(0).at(210), 60068);
  ASSERT_EQ(result.at(0).at(211), 8438);
  ASSERT_EQ(result.at(0).at(212), 65392);
  ASSERT_EQ(result.at(0).at(213), 26387);
  ASSERT_EQ(result.at(0).at(214), 59481);
  ASSERT_EQ(result.at(0).at(215), 111024);
  ASSERT_EQ(result.at(0).at(216), 32372);
  ASSERT_EQ(result.at(0).at(217), 32511);
  ASSERT_EQ(result.at(0).at(218), 102985);
  ASSERT_EQ(result.at(0).at(219), 57293);
  ASSERT_EQ(result.at(0).at(220), 81606);
  ASSERT_EQ(result.at(0).at(221), 61015);
  ASSERT_EQ(result.at(0).at(222), 7573);
  ASSERT_EQ(result.at(0).at(223), 70122);
  ASSERT_EQ(result.at(0).at(224), 107478);
  ASSERT_EQ(result.at(0).at(225), 2383);
  ASSERT_EQ(result.at(0).at(226), 24452);
  ASSERT_EQ(result.at(0).at(227), 93793);
  ASSERT_EQ(result.at(0).at(228), 53897);
  ASSERT_EQ(result.at(0).at(229), 83610);
  ASSERT_EQ(result.at(0).at(230), 41360);
  ASSERT_EQ(result.at(0).at(231), 16855);
  ASSERT_EQ(result.at(0).at(232), 25743);
  ASSERT_EQ(result.at(0).at(233), 45114);
  ASSERT_EQ(result.at(0).at(234), 55496);
  ASSERT_EQ(result.at(0).at(235), 35692);
  ASSERT_EQ(result.at(0).at(236), 89599);
  ASSERT_EQ(result.at(0).at(237), 92849);
  ASSERT_EQ(result.at(0).at(238), 36328);
  ASSERT_EQ(result.at(0).at(239), 44590);
  ASSERT_EQ(result.at(0).at(240), 37592);
  ASSERT_EQ(result.at(0).at(241), 53286);
  ASSERT_EQ(result.at(0).at(242), 63744);
  ASSERT_EQ(result.at(0).at(243), 85566);
  ASSERT_EQ(result.at(0).at(244), 69713);
  ASSERT_EQ(result.at(0).at(245), 34224);
  ASSERT_EQ(result.at(0).at(246), 80214);
  ASSERT_EQ(result.at(0).at(247), 87806);
  ASSERT_EQ(result.at(0).at(248), 16303);
  ASSERT_EQ(result.at(0).at(249), 78283);
  ASSERT_EQ(result.at(1).at(0), 71178);
  ASSERT_EQ(result.at(1).at(1), 0);
  ASSERT_EQ(result.at(1).at(2), 56820);
  ASSERT_EQ(result.at(1).at(3), 50917);
  ASSERT_EQ(result.at(1).at(4), 29771);
  ASSERT_EQ(result.at(1).at(5), 40955);
  ASSERT_EQ(result.at(1).at(6), 55449);
  ASSERT_EQ(result.at(1).at(7), 52742);
  ASSERT_EQ(result.at(1).at(8), 82827);
  ASSERT_EQ(result.at(1).at(9), 81010);
  ASSERT_EQ(result.at(1).at(10), 99361);
  ASSERT_EQ(result.at(1).at(11), 87531);
  ASSERT_EQ(result.at(1).at(12), 32606);
  ASSERT_EQ(result.at(1).at(13), 69636);
  ASSERT_EQ(result.at(1).at(14), 56730);
  ASSERT_EQ(result.at(1).at(15), 76896);
  ASSERT_EQ(result.at(1).at(16), 57179);
  ASSERT_EQ(result.at(1).at(17), 27288);
  ASSERT_EQ(result.at(1).at(18), 49395);
  ASSERT_EQ(result.at(1).at(19), 75902);
  ASSERT_EQ(result.at(1).at(20), 17127);
  ASSERT_EQ(result.at(1).at(21), 38359);
  ASSERT_EQ(result.at(1).at(22), 33598);
  ASSERT_EQ(result.at(1).at(23), 78843);
  ASSERT_EQ(result.at(1).at(24), 81369);
  ASSERT_EQ(result.at(1).at(25), 29118);
  ASSERT_EQ(result.at(1).at(26), 38066);
  ASSERT_EQ(result.at(1).at(27), 38833);
  ASSERT_EQ(result.at(1).at(28), 29211);
  ASSERT_EQ(result.at(1).at(29), 50775);
  ASSERT_EQ(result.at(1).at(30), 78842);
  ASSERT_EQ(result.at(1).at(31), 65812);
  ASSERT_EQ(result.at(1).at(32), 51085);
  ASSERT_EQ(result.at(1).at(33), 83586);
  ASSERT_EQ(result.at(1).at(34), 36489);
  ASSERT_EQ(result.at(1).at(35), 38685);
  ASSERT_EQ(result.at(1).at(36), 37041);
  ASSERT_EQ(result.at(1).at(37), 59429);
  ASSERT_EQ(result.at(1).at(38), 61898);
  ASSERT_EQ(result.at(1).at(39), 83624);
  ASSERT_EQ(result.at(1).at(40), 13578);
  ASSERT_EQ(result.at(1).at(41), 31864);
  ASSERT_EQ(result.at(1).at(42), 54573);
  ASSERT_EQ(result.at(1).at(43), 79195);
  ASSERT_EQ(result.at(1).at(44), 64707);
  ASSERT_EQ(result.at(1).at(45), 47417);
  ASSERT_EQ(result.at(1).at(46), 23046);
  ASSERT_EQ(result.at(1).at(47), 52673);
  ASSERT_EQ(result.at(1).at(48), 57312);
  ASSERT_EQ(result.at(1).at(49), 66653);
  ASSERT_EQ(result.at(1).at(50), 64313);
}

TEST(FWTest, random30) {
  randomTest(floydWarshallAPSP<int>, 30, 3'289'211);
}

TEST(FWTest, random100a) {
  randomTest(floydWarshallAPSP<int>, 100, 987'129, 0.1);
}

TEST(FWTest, random100b) {
  randomTest(floydWarshallAPSP<int>, 100, 2'382'398, 0.1);
}

TEST(FWTest, random100c) {
  randomTest(floydWarshallAPSP<int>, 100, 329'823, 0.1);
}

TEST(FWTest, random100d) {
  randomTest(floydWarshallAPSP<int>, 100, 389'211, 0.1);
}

TEST(FWTest, random100e) {
  randomTest(floydWarshallAPSP<int>, 100, 10'918'101, 0.1);
}

TEST(FWTest, random150_0) {
  randomTest(floydWarshallAPSP<int>, 150, 1812528527, 0.1);
}

TEST(FWTest, random150_1) {
  randomTest(floydWarshallAPSP<int>, 150, 1915543033, 0.1);
}

TEST(FWTest, random150_2) {
  randomTest(floydWarshallAPSP<int>, 150, 1728632203, 0.1);
}

TEST(FWTest, random150_3) {
  randomTest(floydWarshallAPSP<int>, 150, 3400310732, 0.1);
}

TEST(FWTest, random150_4) {
  randomTest(floydWarshallAPSP<int>, 150, 1509646079, 0.1);
}

TEST(FWTest, random150_5) {
  randomTest(floydWarshallAPSP<int>, 150, 877500092, 0.1);
}

TEST(FWTest, random150_6) {
  randomTest(floydWarshallAPSP<int>, 150, 307693294, 0.1);
}

TEST(FWTest, random150_7) {
  randomTest(floydWarshallAPSP<int>, 150, 1740768105, 0.1);
}

TEST(FWTest, random150_8) {
  randomTest(floydWarshallAPSP<int>, 150, 2445968672, 0.1);
}

TEST(FWTest, random150_9) {
  randomTest(floydWarshallAPSP<int>, 150, 55419902, 0.1);
}

TEST(FWTest, random150_10) {
  randomTest(floydWarshallAPSP<int>, 150, 449020997, 0.1);
}

TEST(FWTest, random150_11) {
  randomTest(floydWarshallAPSP<int>, 150, 1922416908, 0.1);
}

TEST(FWTest, random150_12) {
  randomTest(floydWarshallAPSP<int>, 150, 1851447901, 0.1);
}

TEST(FWTest, random200a) {
  randomTest(floydWarshallAPSP<int>, 200, 1'238'923, 0.1);
}

TEST(FWTest, random200b) {
  randomTest(floydWarshallAPSP<int>, 200, 108'181'745, 0.1);
}

TEST(FWTest, random200c) {
  randomTest(floydWarshallAPSP<int>, 200, 3'281'918, 0.1);
}

TEST(FWTest, random200_1) {
  randomTest(floydWarshallAPSP<int>, 200, 445192270, 0.1);
}

TEST(FWTest, random200_2) {
  randomTest(floydWarshallAPSP<int>, 200, 1609628408, 0.1);
}

TEST(FWTest, random200_3) {
  randomTest(floydWarshallAPSP<int>, 200, 1801774926, 0.1);
}

TEST(FWTest, random200_4) {
  randomTest(floydWarshallAPSP<int>, 200, 2002738901, 0.1);
}

TEST(FWTest, random200_5) {
  randomTest(floydWarshallAPSP<int>, 200, 1138962375, 0.1);
}

TEST(FWTest, random200_6) {
  randomTest(floydWarshallAPSP<int>, 200, 1326754919, 0.1);
}

TEST(FWTest, random200_7) {
  randomTest(floydWarshallAPSP<int>, 200, 2178494701, 0.1);
}

TEST(FWTest, random200_8) {
  randomTest(floydWarshallAPSP<int>, 200, 895276633, 0.1);
}

TEST(FWTest, random200_9) {
  randomTest(floydWarshallAPSP<int>, 200, 419725955, 0.1);
}

TEST(FWTest, random200_10) {
  randomTest(floydWarshallAPSP<int>, 200, 2031438742, 0.1);
}

TEST(FWTest, random200_11) {
  randomTest(floydWarshallAPSP<int>, 200, 3188665175, 0.1);
}

TEST(FWTest, random200_12) {
  randomTest(floydWarshallAPSP<int>, 200, 4014034417, 0.1);
}

TEST(FWTest, random200_13) {
  randomTest(floydWarshallAPSP<int>, 200, 2553198683, 0.1);
}

TEST(FWTest, random200_14) {
  randomTest(floydWarshallAPSP<int>, 200, 1186300218, 0.1);
}

TEST(FWTest, random200_15) {
  randomTest(floydWarshallAPSP<int>, 200, 3395108580, 0.1);
}

TEST(FWTest, random200_16) {
  randomTest(floydWarshallAPSP<int>, 200, 148678007, 0.1);
}

TEST(FWTest, random250a) {
  randomTest(floydWarshallAPSP<int>, 200, 23'892'389, 0.1);
}

TEST(FWTest, random250b) {
  randomTest(floydWarshallAPSP<int>, 200, 34'894'398, 0.1);
}

TEST(FWTest, random250c) {
  randomTest(floydWarshallAPSP<int>, 200, 890'119'109, 0.1);
}

TEST(FWTest, random300a) {
  randomTest(floydWarshallAPSP<int>, 300, 98'982, 0.05);
}

TEST(FWTest, random300b) {
  randomTest(floydWarshallAPSP<int>, 300, 32'892, 0.05);
}

TEST(FWTest, random300c) {
  randomTest(floydWarshallAPSP<int>, 300, 238'923, 0.05);
}

TEST(FWTest, random400) {
  randomTest(floydWarshallAPSP<int>, 400, 89'239'823, 0.01);
}

TEST(FWTest, random500) {
  randomTest(floydWarshallAPSP<int>, 500, 2'829'211, 0.01);
}

// TEST(FWTest, random1000a) {
//   randomTest(floydWarshallAPSP<int>, 1000, 32'218'119, 0.005);
// }
//
// TEST(FWTest, random1000b) {
//   randomTest(floydWarshallAPSP<int>, 1000, 98'239'283, 0.005);
// }
//
// TEST(FWTest, random1000c) {
//   randomTest(floydWarshallAPSP<int>, 1'000, 2'389'239, 0.005);
// }

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
