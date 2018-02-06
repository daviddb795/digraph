// Digraph.hpp
//
// ICS 46 Spring 2017
// Project #4: Rock and Roll Stops the Traffic
//
// This header file declares a class template called Digraph, which is
// intended to implement a generic directed graph.  The implementation
// uses the adjacency lists technique, so each vertex stores a linked
// list of its outgoing edges.
//
// Along with the Digraph class template is a class DigraphException
// and a couple of utility structs that aren't generally useful outside
// of this header file.
//
// In general, directed graphs are all the same, except in the sense
// that they store different kinds of information about each vertex and
// about each edge; these two types are the type parameters to the
// Digraph class template.

#ifndef DIGRAPH_HPP
#define DIGRAPH_HPP

#include <functional>
#include <list>
#include <map>
#include <utility>
#include <vector>
#include <iostream> 
#include <queue>



// DigraphExceptions are thrown from some of the member functions in the
// Digraph class template, so that exception is declared here, so it
// will be available to any code that includes this header file.

class DigraphException
{
public:
    DigraphException(const std::string& reason): reason_{reason} { }

    std::string reason() const { return reason_; }

private:
    std::string reason_;
};



// A DigraphEdge lists a "from vertex" (the number of the vertex from which
// the edge points), a "to vertex" (the number of the vertex to which the
// edge points), and an EdgeInfo object.  Because different kinds of Digraphs
// store different kinds of edge information, DigraphEdge is a template
// struct.

template <typename EdgeInfo>
struct DigraphEdge
{
    int fromVertex;
    int toVertex;
    EdgeInfo einfo;
};



// A DigraphVertex includes two things: a VertexInfo object and a list of
// its outgoing edges.  Because different kinds of Digraphs store different
// kinds of vertex and edge information, DigraphVertex is a template struct.

template <typename VertexInfo, typename EdgeInfo>
struct DigraphVertex
{
    VertexInfo vinfo;
    std::list<DigraphEdge<EdgeInfo>> edges;
};



// Digraph is a class template that represents a directed graph implemented
// using adjacency lists.  It takes two type parameters:
//
// * VertexInfo, which specifies the kind of object stored for each vertex
// * EdgeInfo, which specifies the kind of object stored for each edge
//
// You'll need to implement the member functions declared here; each has a
// comment detailing how it is intended to work.
//
// Each vertex in a Digraph is identified uniquely by a "vertex number".
// Vertex numbers are not necessarily sequential and they are not necessarily
// zero- or one-based.

template <typename VertexInfo, typename EdgeInfo>
class Digraph
{
public:
    // The default constructor initializes a new, empty Digraph so that
    // contains no vertices and no edges.
    Digraph();

    // The copy constructor initializes a new Digraph to be a deep copy
    // of another one (i.e., any change to the copy will not affect the
    // original).
    Digraph(const Digraph& d);

    // The move constructor initializes a new Digraph from an expiring one.
    Digraph(Digraph&& d);

    // The destructor deallocates any memory associated with the Digraph.
    ~Digraph();

    // The assignment operator assigns the contents of the given Digraph
    // into "this" Digraph, with "this" Digraph becoming a separate, deep
    // copy of the contents of the given one (i.e., any change made to
    // "this" Digraph afterward will not affect the other).
    Digraph& operator=(const Digraph& d);

    // The move assignment operator assigns the contents of an expiring
    // Digraph into "this" Digraph.
    Digraph& operator=(Digraph&& d);

    // vertices() returns a std::vector containing the vertex numbers of
    // every vertex in this Digraph.
    std::vector<int> vertices() const;

    // edges() returns a std::vector of std::pairs, in which each pair
    // contains the "from" and "to" vertex numbers of an edge in this
    // Digraph.  All edges are included in the std::vector.
    std::vector<std::pair<int, int>> edges() const;

    // This overload of edges() returns a std::vector of std::pairs, in
    // which each pair contains the "from" and "to" vertex numbers of an
    // edge in this Digraph.  Only edges outgoing from the given vertex
    // number are included in the std::vector.  If the given vertex does
    // not exist, a DigraphException is thrown instead.
    std::vector<std::pair<int, int>> edges(int vertex) const;

    // vertexInfo() returns the VertexInfo object belonging to the vertex
    // with the given vertex number.  If that vertex does not exist, a
    // DigraphException is thrown instead.
    VertexInfo vertexInfo(int vertex) const;

    // edgeInfo() returns the EdgeInfo object belonging to the edge
    // with the given "from" and "to" vertex numbers.  If either of those
    // vertices does not exist *or* if the edge does not exist, a
    // DigraphException is thrown instead.
    EdgeInfo edgeInfo(int fromVertex, int toVertex) const;

    // addVertex() adds a vertex to the Digraph with the given vertex
    // number and VertexInfo object.  If there is already a vertex in
    // the graph with the given vertex number, a DigraphException is
    // thrown instead.
    void addVertex(int vertex, const VertexInfo& vinfo);

    // addEdge() adds an edge to the Digraph pointing from the given
    // "from" vertex number to the given "to" vertex number, and
    // associates with the given EdgeInfo object with it.  If one
    // of the vertices does not exist *or* if the same edge is already
    // present in the graph, a DigraphException is thrown instead.
    void addEdge(int fromVertex, int toVertex, const EdgeInfo& einfo);

    // removeVertex() removes the vertex (and all of its incoming
    // and outgoing edges) with the given vertex number from the
    // Digraph.  If the vertex does not exist already, a DigraphException
    // is thrown instead.
    void removeVertex(int vertex);

    // removeEdge() removes the edge pointing from the given "from"
    // vertex number to the given "to" vertex number from the Digraph.
    // If either of these vertices does not exist *or* if the edge
    // is not already present in the graph, a DigraphException is
    // thrown instead.
    void removeEdge(int fromVertex, int toVertex);

    // vertexCount() returns the number of vertices in the graph.
    int vertexCount() const;

    // edgeCount() returns the total number of edges in the graph,
    // counting edges outgoing from all vertices.
    int edgeCount() const;

    // This overload of edgeCount() returns the number of edges in
    // the graph that are outgoing from the given vertex number.
    // If the given vertex does not exist, a DigraphException is
    // thrown instead.
    int edgeCount(int vertex) const;

    // isStronglyConnected() returns true if the Digraph is strongly
    // connected (i.e., every vertex is reachable from every other),
    // false otherwise.
    bool isStronglyConnected() const;

    // findShortestPaths() takes a start vertex number and a function
    // that takes an EdgeInfo object and determines an edge weight.
    // It uses Dijkstra's Shortest Path Algorithm to determine the
    // shortest paths from the start vertex to every other vertex
    // in the graph.  The result is returned as a std::map<int, int>
    // where the keys are vertex numbers and the value associated
    // with each key k is the precedessor of that vertex chosen by
    // the algorithm.  For any vertex without a predecessor (e.g.,
    // a vertex that was never reached, or the start vertex itself),
    // the value is simply a copy of the key.
    std::map<int, int> findShortestPaths(
        int startVertex,
        std::function<double(const EdgeInfo&)> edgeWeightFunc) const;

    class CompareClass
    {
    public:
        bool operator()(std::pair<int, double> p1, std::pair<int, double> p2)
        {
            if(p1.second < p2.second)
            return true;
        else
            return false;
        }
    };
    


private:

    std::map<int, DigraphVertex<VertexInfo, EdgeInfo>> theGraph;
    std::map<int, int> shortestPath;
    int numEdges;
    int numVertices;
    // Add whatever member variables you think you need here.  One
    // possibility is a std::map where the keys are vertex numbers
    // and the values are DigraphVertex<VertexInfo, EdgeInfo> objects.


    // You can also feel free to add any additional member functions
    // you'd like (public or private), so long as you don't remove or
    // change the signatures of the ones that already exist.
};


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

    // The default constructor initializes a new, empty Digraph so that
    // contains no vertices and no edges.
    template <typename VertexInfo, typename EdgeInfo>
    Digraph<VertexInfo, EdgeInfo>::Digraph() : numEdges(0), numVertices(0)
    {
    }

    // The copy constructor initializes a new Digraph to be a deep copy
    // of another one (i.e., any change to the copy will not affect the
    // original).
    template <typename VertexInfo, typename EdgeInfo>
    Digraph<VertexInfo, EdgeInfo>::Digraph(const Digraph& d)
    {

        std::map<int, DigraphVertex<VertexInfo, EdgeInfo>> newGraph;

        for( auto it = theGraph.begin(); it != theGraph.end(); it++ )
        {
            newGraph[it->first] = theGraph[it->first];
        }

    }

    // The move constructor initializes a new Digraph from an expiring one.
    template <typename VertexInfo, typename EdgeInfo>
    Digraph<VertexInfo, EdgeInfo>::Digraph(Digraph&& d)
    {
        std::map<int, DigraphVertex<VertexInfo, EdgeInfo>> newGraph;

        theGraph.swap(newGraph);

    }

    // The destructor deallocates any memory associated with the Digraph.
    template <typename VertexInfo, typename EdgeInfo>
    Digraph<VertexInfo, EdgeInfo>::~Digraph()
    {
    }

    // The assignment operator assigns the contents of the given Digraph
    // into "this" Digraph, with "this" Digraph becoming a separate, deep
    // copy of the contents of the given one (i.e., any change made to
    // "this" Digraph afterward will not affect the other).
    template <typename VertexInfo, typename EdgeInfo>
   Digraph<VertexInfo, EdgeInfo>&  Digraph<VertexInfo, EdgeInfo>::operator=(const Digraph& d)
    {
        if(this != &d)
        {
            theGraph.clear();
            numVertices = d.numVertices;
            numEdges = d.numEdges;

            for( auto it = d.theGraph.begin(); it != d.theGraph.end(); it++ )
            {
                theGraph[it->first] = d.theGraph.at(it->first);
            }
        }
        return *this;

    }

    // The move assignment operator assigns the contents of an expiring
    // Digraph into "this" Digraph.
    template <typename VertexInfo, typename EdgeInfo>
    Digraph<VertexInfo, EdgeInfo>& Digraph<VertexInfo, EdgeInfo>::operator=(Digraph&& d)
    {
        theGraph = std::move(d.theGraph);
        return *this;
    }

    // vertices() returns a std::vector containing the vertex numbers of
    // every vertex in this Digraph.
    template <typename VertexInfo, typename EdgeInfo>
    std::vector<int> Digraph<VertexInfo, EdgeInfo>::vertices() const
    {
        std::vector<int> vertices;

        for(auto it = theGraph.begin(); it != theGraph.end(); ++it)
        {
            vertices.push_back(it -> first);
        }
        
        return vertices;
    }

    // edges() returns a std::vector of std::pairs, in which each pair
    // contains the "from" and "to" vertex numbers of an edge in this
    // Digraph.  All edges are included in the std::vector.
    template <typename VertexInfo, typename EdgeInfo>
    std::vector<std::pair<int, int>> Digraph<VertexInfo, EdgeInfo>::edges() const
    {
        std::vector<std::pair<int, int>> edgeFinder;

        for( int i = 0; i < theGraph.size(); ++i)
        {
            for(auto it = theGraph.at(i).edges.begin(); it != theGraph.at(i).edges.end(); ++it)
            {   
                edgeFinder.push_back(std::make_pair(it->fromVertex, it->toVertex));
            }
        }
        
        return edgeFinder;

    }

    // This overload of edges() returns a std::vector of std::pairs, in
    // which each pair contains the "from" and "to" vertex numbers of an
    // edge in this Digraph.  Only edges outgoing from the given vertex
    // number are included in the std::vector.  If the given vertex does
    // not exist, a DigraphException is thrown instead.
    template <typename VertexInfo, typename EdgeInfo>
    std::vector<std::pair<int, int>> Digraph<VertexInfo, EdgeInfo>::edges(int vertex) const
    {
         std::vector<std::pair<int, int>> edgeFinder;
         
         if(theGraph.find(vertex) == theGraph.end()) 
        {
            throw new DigraphException("Vertex already exists");
        }

        else
        {

         for(auto it = theGraph.at(vertex).edges.begin(); it != theGraph.at(vertex).edges.end(); ++it)
            {   
                edgeFinder.push_back(std::make_pair(it->fromVertex, it->toVertex));
            }

            return edgeFinder;
        }
    }

    // vertexInfo() returns the VertexInfo object belonging to the vertex
    // with the given vertex number.  If that vertex does not exist, a
    // DigraphException is thrown instead.
    template <typename VertexInfo, typename EdgeInfo>
    VertexInfo Digraph<VertexInfo, EdgeInfo>::vertexInfo(int vertex) const
    {
        if(theGraph.find(vertex) == theGraph.end()) 
        {
            throw new DigraphException("Vertex already exists");
        }
        else
        return theGraph.at(vertex).vinfo;
    }

    // edgeInfo() returns the EdgeInfo object belonging to the edge
    // with the given "from" and "to" vertex numbers.  If either of those
    // vertices does not exist *or* if the edge does not exist, a
    // DigraphException is thrown instead.
    template <typename VertexInfo, typename EdgeInfo>
    EdgeInfo Digraph<VertexInfo, EdgeInfo>::edgeInfo(int fromVertex, int toVertex) const
    {
        EdgeInfo edgeInfo;

        if((theGraph.find(toVertex) == theGraph.end()) || (theGraph.find(fromVertex) == theGraph.end()))
        {
            throw new DigraphException("Vertex does not exist");
        }
        else
        {

        for(auto it = theGraph.at(fromVertex).edges.begin(); it != theGraph.at(fromVertex).edges.end(); ++it)
            {
                if((it->fromVertex == fromVertex) && (it->toVertex == toVertex))
                {
                    edgeInfo = it->einfo;
                    break;
                }
            }
            //throw new DigraphException("Edge does not exist");
        }    
        return edgeInfo;
    }

    // addVertex() adds a vertex to the Digraph with the given vertex
    // number and VertexInfo object.  If there is already a vertex in
    // the graph with the given vertex number, a DigraphException is
    // thrown instead.
    template <typename VertexInfo, typename EdgeInfo>
    void Digraph<VertexInfo, EdgeInfo>::addVertex(int vertex, const VertexInfo& vinfo)
    {
        if(theGraph.find(vertex) == theGraph.end()) 
        {
            DigraphVertex<VertexInfo, EdgeInfo> newVertex;
            newVertex.vinfo = vinfo;
            theGraph.emplace(vertex, newVertex);
            ++numVertices;
        }

        else
        {
            throw new DigraphException("Vertex already exists");
        }

    }

    // addEdge() adds an edge to the Digraph pointing from the given
    // "from" vertex number to the given "to" vertex number, and
    // associates with the given EdgeInfo object with it.  If one
    // of the vertices does not exist *or* if the same edge is already
    // present in the graph, a DigraphException is thrown instead.
    template <typename VertexInfo, typename EdgeInfo>
    void Digraph<VertexInfo, EdgeInfo>::addEdge(int fromVertex, int toVertex, const EdgeInfo& einfo)
    {
        if(theGraph.find(fromVertex) == theGraph.end())
        {
            throw new DigraphException("Vertex not found in graph");
        }
        else
        {
            DigraphEdge<EdgeInfo> newEdge;
            newEdge.fromVertex = fromVertex;
            newEdge.toVertex = toVertex;
            newEdge.einfo = einfo;

            for(auto it = theGraph.at(fromVertex).edges.begin(); it != theGraph.at(fromVertex).edges.end(); ++it)
            {
                if((it->fromVertex == fromVertex) && (it->toVertex == toVertex))
                {
                    throw new DigraphException("Edge already exist");
                    return;
                }
            }
            theGraph.at(fromVertex).edges.push_back(newEdge);
            ++numEdges;
        }

    }

    // removeVertex() removes the vertex (and all of its incoming
    // and outgoing edges) with the given vertex number from the
    // Digraph.  If the vertex does not exist already, a DigraphException
    // is thrown instead.
    template <typename VertexInfo, typename EdgeInfo>
    void Digraph<VertexInfo, EdgeInfo>::removeVertex(int vertex)
    {
        if(theGraph.find(vertex) == theGraph.end()) 
        {
            throw new DigraphException("Vertex does not exist");
        }
        else
        {
            theGraph.erase(vertex);
            --numVertices;
        }

    }

    // removeEdge() removes the edge pointing from the given "from"
    // vertex number to the given "to" vertex number from the Digraph.
    // If either of these vertices does not exist *or* if the edge
    // is not already present in the graph, a DigraphException is
    // thrown instead.
    template <typename VertexInfo, typename EdgeInfo>
    void Digraph<VertexInfo, EdgeInfo>::removeEdge(int fromVertex, int toVertex)
    {
        if((theGraph.find(toVertex) == theGraph.end()) || (theGraph.find(toVertex) == theGraph.end()))
        {
            throw new DigraphException("Vertex does not exist");
        }

        else
        {

        for(auto it = theGraph.at(fromVertex).edges.begin(); it != theGraph.at(fromVertex).edges.end(); ++it)
            {
                if((it->fromVertex == fromVertex) && (it->toVertex == toVertex))
                {
                    
                    theGraph.at(fromVertex).edges.erase(it);
                    --numEdges;
                    return;
                }
            }
            
            throw new DigraphException("Edge does not exist");
        }

    }

    // vertexCount() returns the number of vertices in the graph.
    template <typename VertexInfo, typename EdgeInfo>
    int Digraph<VertexInfo, EdgeInfo>::vertexCount() const
    {
        return numVertices;
    }

    // edgeCount() returns the total number of edges in the graph,
    // counting edges outgoing from all vertices.
    template <typename VertexInfo, typename EdgeInfo>
    int Digraph<VertexInfo, EdgeInfo>::edgeCount() const
    {
        return numEdges;
    }

    // This overload of edgeCount() returns the number of edges in
    // the graph that are outgoing from the given vertex number.
    // If the given vertex does not exist, a DigraphException is
    // thrown instead.
    template <typename VertexInfo, typename EdgeInfo>
    int Digraph<VertexInfo, EdgeInfo>::edgeCount(int vertex) const
    {
        int tempEdgeNum = 0;
        for(auto it = theGraph.at(vertex).edges.begin(); it != theGraph.at(vertex).edges.end(); ++it)
            {
                ++tempEdgeNum;
            }

            return tempEdgeNum;
    }

    // isStronglyConnected() returns true if the Digraph is strongly
    // connected (i.e., every vertex is reachable from every other),
    // false otherwise.
    template <typename VertexInfo, typename EdgeInfo>
    bool Digraph<VertexInfo, EdgeInfo>::isStronglyConnected() const
    {


    }
    
    // findShortestPaths() takes a start vertex number and a function
    // that takes an EdgeInfo object and determines an edge weight.
    // It uses Dijkstra's Shortest Path Algorithm to determine the
    // shortest paths from the start vertex to every other vertex
    // in the graph.  The result is returned as a std::map<int, int>
    // where the keys are vertex numbers and the value associated
    // with each key k is the precedessor of that vertex chosen by
    // the algorithm.  For any vertex without a predecessor (e.g.,
    // a vertex that was never reached, or the start vertex itself),
    // the value is simply a copy of the key.
    template <typename VertexInfo, typename EdgeInfo>
    std::map<int, int> Digraph<VertexInfo, EdgeInfo>::findShortestPaths(
        int startVertex,
        std::function<double(const EdgeInfo&)> edgeWeightFunc) const
    {
        const int INF = std::numeric_limits<int>::infinity();
        
        std::map<int, bool> shortestPathKnown;
        std::map<int, int>  shortestPathGraph;
        std::map<int, double>  distanceOfPath;
        std::pair<int,double> vertexEdgeWeight;
        std::vector<int> allVertices = vertices();

        std::priority_queue<std::pair<int, double>, std::vector<std::pair<int,int>>, CompareClass> pQueue;

        for(int i = 0; i < allVertices.size(); ++i)
        {
            shortestPathKnown.emplace(allVertices[i], false);

            if(allVertices[i] == startVertex)
            {
                shortestPathGraph.emplace(allVertices[i], allVertices[i]);
                distanceOfPath.emplace(allVertices[i], allVertices[i]);
            }
            shortestPathGraph.emplace(allVertices[i], 0);
            distanceOfPath.emplace(allVertices[i], INF);
        }
        // for(auto it = theGraph.begin(); it != theGraph.end(); ++it)
        // {
        //     //shortestPathKnown.emplace(it->first, false) = theGraph.at(it->first);
        //     shortestPathKnown.at(it->second) = false;

        //     shortestPathGraph.at(it->first) = theGraph.at(it->first);
        //     if(theGraph.at(it->first) == startVertex)
        //     {
        //         shortestPathGraph.at(it->first) = theGraph.at(it->first);
        //         shortestPathGraph.at(it->second) = theGraph.at(it->first);
        //     }
        //     shortestPathGraph.at(it->second) = 0;

        //     distanceOfPath.at(it->first) = theGraph.at(it->first);
        //     if(theGraph.at(it->first) == startVertex)
        //     {
        //         distanceOfPath.at(it->first) = theGraph.at(it->first);
        //         distanceOfPath.at(it->second) = 0;
        //     }
        //     distanceOfPath.at(it->second) = INF;

        // }

        vertexEdgeWeight.first = startVertex;
        vertexEdgeWeight.second = 0;

        while(!pQueue.empty())
        {
            int vTex = pQueue.top().first;
            pQueue.pop();

            if(shortestPathKnown.at(vTex) == false)
            {
                shortestPathKnown.at(vTex) = true; 
            
                 for(auto it = theGraph.at(vTex).edges.begin(); it != theGraph.at(vTex).edges.end(); ++it)
                {   

                    for(int i = 0; i < theGraph.size(); ++i)
                    {
                        if((allVertices[i] == it->toVertex) && (shortestPathGraph.at(allVertices[i] == false)))
                        {
                            if(distanceOfPath.at(allVertices[i]) > edgeWeightFunc(it->einfo))
                            {
                                distanceOfPath.at(allVertices[i]) = edgeWeightFunc(it->einfo);
                                shortestPathGraph.at(allVertices[i]) = vTex;

                                vertexEdgeWeight.first = allVertices[i];
                                vertexEdgeWeight.second = distanceOfPath.at(allVertices[i]);
                                pQueue.push(vertexEdgeWeight);
                            }
                        }
                    }




                }
            } 
        }
        return shortestPathGraph;
    }


#endif // DIGRAPH_HPP

