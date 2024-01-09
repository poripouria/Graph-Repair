//POURIA ALIMORADPOR 9912035
//Data Structure (Dr.Akbari) Final Project

#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>
#include <exception>
#include <stdexcept>
using namespace std ;

typedef pair<int, char> neighbor_weight ;
typedef pair<int, int> sorc_dest ;

class WeightedGraph
{
public:
    WeightedGraph(unsigned int) ;
    ~WeightedGraph() {delete[] adjList ;}
    void addDirectedEdge(int, int, char) ;
    void DFS_Repair(char, char, char) ;
    void BFS_Repair(char, char, char) ;
    void repair(char, char, char) ;

    int getvertexNum() const {return vertexNum ;}
    int getedgeNum() const {return edgeNum ;}
    void printGraphInfo() const ;
private:
    unsigned int vertexNum ;
    unsigned int edgeNum ;
    list<neighbor_weight> *adjList ;
    list< pair<char, sorc_dest> > edges ;

    void DFS_Repair_Visit(int, bool *, bool *, char, char, char) ;
    void repairCheck(int, bool *, char, char, char) ;
};


int main()
{
    try
    {
        unsigned vertexNum, firstV, secondV ;
        char edgeWeight ;
        char firstE, secondE, newE ;
        unsigned int addedEdges = 0 ;

        fstream rulesFile ;
        fstream graphFile ;
        rulesFile.open("TestCases/simple_grammar.txt", ios::in) ;
        graphFile.open("TestCases/simple_graph.txt", ios::in) ;
        if(!rulesFile || !graphFile)
        {
            cerr << "Something wrong during opening files!" << endl;
            return 1 ;
        }

        graphFile >> vertexNum ;
        WeightedGraph graph(vertexNum) ;
        while(!graphFile.eof())
        {
            graphFile >> firstV >> secondV >> edgeWeight ;
            graph.addDirectedEdge(firstV, secondV, edgeWeight) ;
        }
        int edgeNumBeforRepair = graph.getedgeNum() ;
        //cout << "Graph info BEFORE repairing it :" << endl ;
        //graph.printGraphInfo() ;
        while(!rulesFile.eof())
        {
            rulesFile >> newE >> firstE >> secondE ;
            graph.repair(firstE, secondE, newE) ;
            //graph.DFS_Repair(firstE, secondE, newE) ;
            //ظgraph.BFS_Repair(firstE, secondE, newE) ;
        }
        int edgeNumAfterRepair = graph.getedgeNum() ;
        //cout << "Graph info AFTER repairing it :" << endl ;
        //graph.printGraphInfo() ;


        addedEdges = edgeNumAfterRepair - edgeNumBeforRepair ;
        cout << addedEdges ;

        rulesFile.close() ;
        graphFile.close() ;
    }
    catch(const exception &err)
    {
        cout << err.what() << endl ;
    }

    return 0 ;
}


WeightedGraph::WeightedGraph(unsigned int v)
    :vertexNum(v), edgeNum(0)
{
    adjList = new list<neighbor_weight>[vertexNum] ;
}

void WeightedGraph::addDirectedEdge(int u, int v, char weight)
{
    if(adjList[u].back() != make_pair(v, weight))
    {
        adjList[u].push_back({v, weight}) ;
        edges.push_back({weight, {u, v}}) ;
        ++edgeNum ;
    }
}

void WeightedGraph::repair(char firstEdge, char secondEdge, char newEdge)
{
    for(auto it1=edges.begin(); it1 != edges.end(); ++it1)
        for(auto it2=edges.begin(); it2 != edges.end(); ++it2)
            if(it1->first == firstEdge && it2->first == secondEdge)
                if(it1->second.second == it2->second.first)
                    addDirectedEdge(it1->second.first, it2->second.second, newEdge) ;

}

void WeightedGraph::DFS_Repair(char firstEdge, char secondEdge, char newEdge)
{
    bool *isVisited = new bool[vertexNum] ;
    bool *isRepaired = new bool[vertexNum] ;
    fill(isVisited, isVisited+vertexNum, false) ;
    fill(isRepaired, isRepaired+vertexNum, false) ;
    for(unsigned int i=0; i<vertexNum; ++i)
        if(!isVisited[i])
            DFS_Repair_Visit(i, isVisited, isRepaired, firstEdge, secondEdge, newEdge) ;

    delete[] isVisited ;
    delete[] isRepaired ;
}

void WeightedGraph::DFS_Repair_Visit(int s, bool *visited, bool *repaired, char firstEdge, char secondEdge, char newEdge)
{
    repairCheck(s, repaired, firstEdge, secondEdge, newEdge) ;
    if(!repaired[s])
        visited[s] = true ;

    for(auto it=adjList[s].begin(); it != adjList[s].end(); ++it)
        if(!visited[it->first])
            DFS_Repair_Visit(it->first, visited, repaired, firstEdge, secondEdge, newEdge) ;
}

void WeightedGraph::BFS_Repair(char firstEdge, char secondEdge, char newEdge)
{
    bool *isVisited = new bool[vertexNum] ;
    bool *isRepaired = new bool[vertexNum] ;
    fill(isVisited, isVisited+vertexNum, false) ;
    fill(isRepaired, isRepaired+vertexNum, false) ;
    list<int> Queue ;
    int s=0 ;

    repairCheck(s, isRepaired, firstEdge, secondEdge, newEdge) ;
    if(!isRepaired[s])
        isVisited[s] = true ;
    Queue.push_back(s) ;
    while(!Queue.empty())
    {
        s = Queue.front() ;
        Queue.pop_front() ;
        repairCheck(s, isRepaired, firstEdge, secondEdge, newEdge) ;
        for(auto i=adjList[s].begin(); i != adjList[s].end(); ++i)
            if(!isVisited[i->first])
            {
                if(!isRepaired[i->first])
                    isVisited[i->first] = true ;
                Queue.push_back(i->first) ;
            }
    }

    delete[] isVisited ;
    delete[] isRepaired ;
}

void WeightedGraph::repairCheck(int s, bool *isRepaired, char firstEdge, char secondEdge, char newEdge)
{
    int v1=s, v2=adjList[v1].begin()->first, v3=adjList[v2].begin()->first ;
    char w1=adjList[v1].begin()->second, w2=adjList[v2].begin()->second ;

    if(w1 == firstEdge && w2 == secondEdge)
    {
        addDirectedEdge(v1, v3, newEdge) ;
        isRepaired[s] = true ;
    }
    else
        isRepaired[s] = false ;
}

void WeightedGraph::printGraphInfo() const
{
    cout << "Number of Vertexes : " << vertexNum << endl ;
    cout << "Number of Edges : " << edgeNum << endl ;

    cout << "Graph Edges : " << endl ;
    for(auto it=edges.begin(); it != edges.end(); ++it)
        cout << it->second.first << '-' << it->second.second << " : " << it->first << endl ;

    cout << endl ;
}
