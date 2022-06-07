#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include<vector>
#include<set>
#include <memory>
#include <queue>
#define UNREFERENCED_PARAMETER(P) (P)
using namespace std;



struct BaseGraph
{
    bool isVisited = false;
    vector<long long int> connections;
    bool isBottomNode = false;
    bool isaddedEdge = false;
    ~BaseGraph() { }
};

struct DeepVertex
{
    int vertex;
    int level;
    DeepVertex(int vertex, int level)
    {
        this->vertex = vertex;
        this->level = level;
    }
};

struct Edge
{
    int vertexA;
    int vertexB;

    Edge(int A, int B)
    {
        this->vertexA = A;
        this->vertexB = B;
    }

};


class Alghotitm
{
public:
    shared_ptr<BaseGraph> graph;
    int nodes = 0;
    int stepsTOFinish = 0;
    int r = 0;

    Alghotitm() {}
    virtual void buildGraph() = 0;
    virtual int calcNodes(int nodes, int r = 0) = 0;

    void bfs(int vertex)
    {

        int sum = 0;
        auto tmpGraph = this->graph.get();
        int tmpAmountOFVertex = calcNodes(vertex, this->r);
        for (int i = 0; i < tmpAmountOFVertex; i++)
        {

            queue<DeepVertex> queue;
            queue.push(DeepVertex(i, 0));


            tmpGraph[i].isVisited = 1;

            while (!queue.empty())
            {
                DeepVertex tmpVertex = queue.front();
                queue.pop();

                for (int j = 0; j < (int)tmpGraph[tmpVertex.vertex].connections.size(); j++)
                {

                    if (!tmpGraph[tmpGraph[tmpVertex.vertex].connections[j]].isVisited && tmpGraph[tmpVertex.vertex].connections[j] < tmpAmountOFVertex)
                    {
                        queue.push(DeepVertex((int)tmpGraph[tmpVertex.vertex].connections[j], tmpVertex.level + 1));
                        tmpGraph[tmpGraph[tmpVertex.vertex].connections[j]].isVisited = 1; //oznaczenie, ¿e dodano do kolejki

                        if (i < tmpGraph[tmpVertex.vertex].connections[j])
                        {
                            sum += tmpVertex.level + 1;
                        }
                    }
                }

            }
            this->setUnvisitedVetrexInGraph(0);
        }
        cout << sum << endl;
    }


    void getGraph()
    {
        cout << "get graph" << endl;
        auto tmpGraph = this->graph.get();
        for (int i = 0; i < this->nodes; i++)
        {
            cout << "wierzcholek " << i << " polaczono z ";
            for (int j = 0; j < (int)tmpGraph[i].connections.size(); j++)
            {
                cout << tmpGraph[i].connections[j] << " ";
            }
            cout << endl;
        }
    }

    void setUnvisitedVetrexInGraph(int i)
    {
        auto tmpGraph = this->graph.get();
        for (int j = i; j < this->nodes; j++)
            tmpGraph[j].isVisited = 0;
    }

};

// Model wzrostowo-iteracyjny
class MWI : public Alghotitm
{
    int edges;
public:

    MWI(int k, int r)
    {
        this->stepsTOFinish = k;
        this->r = r;
        this->edges = MWI::calcEdges(k, r);
        this->nodes = this->calcNodes(k, r);
        this->graph = shared_ptr<BaseGraph>(new BaseGraph[this->nodes]);

    }
    static int calcEdges(int k, int r2)
    {
        return ((int)pow((2 * r2), k + 1) - 1) / ((2 * r2) - 1);
    }

    int calcNodes (int k, int r2) override
    {
        return ((r2 * (int)pow(2 * r2, k)) + (3 * r2) - 2) / ((2 * r2) - 1);
    }

    void buildGraph() override
    {

        // krok 0
        auto tmpGraph = this->graph.get();
        tmpGraph[0].connections.push_back(1);
        tmpGraph[1].connections.push_back(0);

        queue<Edge> edge;
        edge.push(Edge(0, 1));

        int tmpNodes = nodes;

        vector<vector<int>> matrix(nodes, vector<int>(nodes, 0));

        matrix[0][1] = 1;
        matrix[1][0] = 1;

        // tymczasowo dodana krawedz

        for (int i = 1; i <= this->stepsTOFinish; i++)
        {

            int edgeAmount = edge.size();
            int vertexAmountBefore = calcNodes(i - 1, this->r) - 1;
            int firstVertexIndex = vertexAmountBefore + 1;

            // ile krawedzi sci¹gn¹æ
            for (int tmpEdge = 0; tmpEdge < edgeAmount; tmpEdge++)
            {

                auto pulledEdge = edge.front();

                edge.pop();
                // zebraæ dodane krawêdzie dla tego kroku
                // dodawanie wezlow do sciagnietej krawedzi
                for (int j = 0; j < this->r; j++)
                {

                    // pierwszy wierzcholek lacze z dodanym wierzcholkiem
                    tmpGraph[pulledEdge.vertexA].connections.push_back(firstVertexIndex);
                    // dodany wierzcholek lacze z pierwszym
                    tmpGraph[firstVertexIndex].connections.push_back(pulledEdge.vertexA);
                    // drugi wierzcholek lacze z dodanym wierzcholkiem
                    tmpGraph[pulledEdge.vertexB].connections.push_back(firstVertexIndex);
                    // dodany wierzcholek lacze z drugim
                    tmpGraph[firstVertexIndex].connections.push_back(pulledEdge.vertexB);

                    edge.push(Edge(pulledEdge.vertexA, firstVertexIndex));
                    edge.push(Edge(pulledEdge.vertexB, firstVertexIndex));

                    // przy dodawaniu
                    matrix[pulledEdge.vertexA][firstVertexIndex] = 1;
                    matrix[firstVertexIndex][pulledEdge.vertexA] = 1;

                    matrix[pulledEdge.vertexB][firstVertexIndex] = 1;
                    matrix[firstVertexIndex][pulledEdge.vertexB] = 1;


                    MWI::addDistance(matrix, pulledEdge.vertexA, pulledEdge.vertexB, firstVertexIndex);

                    firstVertexIndex++;
                }
            }
        }

        MWI::countDistance(matrix, tmpNodes);
    }

    static void addDistance(vector<vector<int>>& matrix, int rowA, int rowB, int addedVertex)
    {
        for (int i = 0; i < addedVertex; i++)
        {
            // jesli ktorys z dodanych mam polaczony z jakims innym wierzcholkiem to mode go polazyc z dodanym zwierkszajac wartosc polaczenie o 1 z krotszej odleglosci pomiedzy nimi
            if ((matrix[rowA][i] > 0 || matrix[rowB][i] > 0) && i != addedVertex && matrix[addedVertex][i] == 0)
            {
                if (matrix[rowA][i] <= matrix[rowB][i])
                {
                    matrix[addedVertex][i] = matrix[rowA][i] + 1;
                    matrix[i][addedVertex] = matrix[rowA][i] + 1;
                }
                else
                {
                    matrix[addedVertex][i] = matrix[rowB][i] + 1;
                    matrix[i][addedVertex] = matrix[rowB][i] + 1;
                }
            }
        }
    };

    static void countDistance(vector<vector<int>>& matrix, int tmpNodes)
    {
        int sum = 0;
        for (int i = 0; i < tmpNodes; i++)
        {
            for (int j = 0; j < tmpNodes; j++)
            {
                sum += matrix[i][j];
            }
        }
        cout << sum / 2 << endl;
    }

};
// Barabasi-Ravasz-Vicsek
class BRV : public Alghotitm
{
public:
    /*int stepsTOFinish;*/
    explicit BRV(int k)
    {
        this->stepsTOFinish = k;
        this->nodes = this->calcNodes(k);
        this->graph = shared_ptr<BaseGraph>(new BaseGraph[this->nodes]);
    };
    int calcNodes(int k, int r2 = 0) override
    {
        this->stepsTOFinish = k;
        UNREFERENCED_PARAMETER(r2);
        return (int)pow(3, k);
    }
    void buildGraph() override
    {
        // krok 1
        auto tmpGraph = this->graph.get();
        tmpGraph[0].connections.push_back(1);
        tmpGraph[0].connections.push_back(2);

        tmpGraph[1].connections.push_back(0);
        tmpGraph[2].connections.push_back(0);

        tmpGraph[1].isBottomNode = 1;
        tmpGraph[2].isBottomNode = 1;

        int step = 2;

        for (int i = step; i <= this->stepsTOFinish; i++)
        {
            int startNode = (int)pow(3, i) / 3;
            int endNode = (startNode * 3) - 1;
            int startCounting = (int)pow(3, i - 1) - startNode;
            int lowLevelNodesLevelUp = (int)pow(2, i - 1);
            int increaseOFVal = (int)pow(3, i - 1);

            int indexFindedLoweLewelNodes = 0;
            for (int j = startNode; j <= endNode; j++)
            {

                // wezly do odznaczenia
                if (tmpGraph[startCounting].isBottomNode == true && indexFindedLoweLewelNodes < lowLevelNodesLevelUp)
                {
                    // odznaczam ze nie sa juz dolne
                    tmpGraph[startCounting].isBottomNode = 0;

                    // nowy dolny wezel
                    tmpGraph[j].isBottomNode = 1;

                    // lacze werzcholek z root
                    tmpGraph[j].connections.push_back(0);
                    // lacze root z wezlem
                    tmpGraph[0].connections.push_back(j);

                    indexFindedLoweLewelNodes++;

                }
                else if (tmpGraph[startCounting].isBottomNode == true && indexFindedLoweLewelNodes >= lowLevelNodesLevelUp)
                {
                    // nie odzanczam
                    tmpGraph[j].isBottomNode = 1;
                    // lacze werzcholek z root
                    tmpGraph[j].connections.push_back(0);
                    // lacze root z wezlem
                    tmpGraph[0].connections.push_back(j);
                    indexFindedLoweLewelNodes++;
                }

                for (int w = 0; w < (int)tmpGraph[startCounting].connections.size(); w++)
                {
                    if (indexFindedLoweLewelNodes <= lowLevelNodesLevelUp)
                        tmpGraph[j].connections.push_back(tmpGraph[startCounting].connections[w] + increaseOFVal);
                    else if (tmpGraph[startCounting].connections[w] != 0)
                        tmpGraph[j].connections.push_back(tmpGraph[startCounting].connections[w] + increaseOFVal);
                }
                startCounting++;
            }

        }
    }
};

//MODEL-DCN
class DCN : public Alghotitm
{
public:

    explicit DCN(int k)
    {
        this->stepsTOFinish = k;
        this->nodes = this->calcNodes(k);
        this->graph = shared_ptr<BaseGraph>(new BaseGraph[this->nodes]);
    }

    static int calcNodes2(int k, int r2 = 0)
    {
        // -1 bo od zera
        UNREFERENCED_PARAMETER(r2);
        return  ((int)pow(2, k + 1) - 1) - 1;
    }

    int calcNodes(int k, int r2 = 0) override
    {
        // -1 bo od zera
        UNREFERENCED_PARAMETER(r2);
        return  k;
    }

    void buildGraph() override
    {
        // krok 1
        auto tmpGraph = this->graph.get();
        tmpGraph[0].connections.push_back(1);
        tmpGraph[0].connections.push_back(2);

        tmpGraph[1].connections.push_back(0);
        tmpGraph[2].connections.push_back(0);
        int step = 2;

        int intVertexForSteps;
        do {

            intVertexForSteps = DCN::calcNodes2(step);

            int firstIndexOnLevel = DCN::calcNodes2(step - 1) + 1;
            int lastIndexOnLevel = DCN::calcNodes2(step);
            int firstIndexLevelUp = DCN::calcNodes2(step - 2) + 1;

            lastIndexOnLevel = lastIndexOnLevel < this->nodes - 1 ? lastIndexOnLevel : this->nodes - 1;

            int increaseUpIndeks = 1;

            for (firstIndexOnLevel; firstIndexOnLevel <= lastIndexOnLevel; firstIndexOnLevel++)
            {
                if (increaseUpIndeks > 2)
                {
                    firstIndexLevelUp++;
                    increaseUpIndeks = 1;
                }

                // lacze wierzcholek z rodzicem
                tmpGraph[firstIndexOnLevel].connections.push_back(firstIndexLevelUp);


                for (int j = 0; j < (int)tmpGraph[firstIndexLevelUp].connections.size(); j++)
                {
                    // wszystkie polaczenia rodzica
                    if (tmpGraph[firstIndexLevelUp].connections[j] <= firstIndexLevelUp)
                    {
                        tmpGraph[firstIndexOnLevel].connections.push_back(tmpGraph[firstIndexLevelUp].connections[j]);
                        tmpGraph[tmpGraph[firstIndexLevelUp].connections[j]].connections.push_back(firstIndexOnLevel);
                    }
                }

                // lacze rodzica z wierzcholkiem
                tmpGraph[firstIndexLevelUp].connections.push_back(firstIndexOnLevel);

                // polaczyc wiercholek ze wszystkimi przodkami w zaleznosci od ilosci krokow
                increaseUpIndeks++;
            }

            step++;

        } while (intVertexForSteps < this->nodes);
    }
};


// Lu-Su-Go
class LuSuGuo : public Alghotitm
{
public:
    explicit LuSuGuo(int k)
    {
        this->nodes = this->calcNodes(k);
        this->graph = shared_ptr<BaseGraph>(new BaseGraph[this->nodes]);
    };
    int calcNodes(int nodes2,int r2 = 0) override
    {
        UNREFERENCED_PARAMETER(r2);
        return nodes2;
    }
    void buildGraph() override
    {

        // krok 1
        auto tmpGraph = this->graph.get();
        tmpGraph[0].connections.push_back(1);
        tmpGraph[0].connections.push_back(2);

        tmpGraph[1].connections.push_back(0);
        tmpGraph[1].connections.push_back(2);

        tmpGraph[2].connections.push_back(0);
        tmpGraph[2].connections.push_back(1);

        int step = 2;

        // wyliczyæ iloœæ niezbêdnych wywo³añ
        int tmp = (int)ceil(sqrt(nodes));
        int endLoop = tmp - step;

        for (int i = step; i <= step + endLoop; i++)
        {
            long long int startParent = (long long int)pow(2, i - 1) - 1;
            long long int startVertex = (long long int)pow(2, i) - 1;
            long long int stopVertex = startVertex * 2 > nodes - 1 ? nodes - 1 : startVertex * 2;

            // indeks dziecka
            int childIndex = 1;
            // polaczenia dla kazdego wierzcholka
            for (long long int vertex = startVertex; vertex <= stopVertex; vertex++)
            {
                // polaczenie z sasiadem
                long long int neighbour;

                if (vertex % 2 == 0)
                    neighbour = vertex - 1;
                else
                    neighbour = vertex + 1 <= stopVertex ? vertex + 1 : stopVertex;

                // nie lacze z samym soba
                if (neighbour != vertex)
                    tmpGraph[vertex].connections.push_back(neighbour);

                long long int parent = startParent;

                // dodanie dziecka do rodzica
                tmpGraph[parent].connections.push_back(vertex);
                // dodanie rodzica do dziecka
                tmpGraph[vertex].connections.push_back(parent);

                childIndex++;
                if (childIndex > 2)
                {
                    childIndex = 1;
                    startParent++;
                }

                long long int treeIndex = (long long int)floor(((vertex - 1) / 2));
                long long int ancient = treeIndex % (i - 1);

                //lacze przodka z wierzcholkiem
                tmpGraph[ancient].connections.push_back(vertex);
                // lacze wierzcholek z przodkiem
                tmpGraph[vertex].connections.push_back(ancient);

            }
        }

    }
};

class GMS
{
public:
    int sumVertecise = 0;
    shared_ptr<Alghotitm> alghoritm;

    explicit GMS(shared_ptr<Alghotitm> alghoritm)
    {
        this->alghoritm = alghoritm;
    };

    void buildGraph()
    {
        alghoritm.get()->buildGraph();
    }

    void getGraph()
    {
        alghoritm.get()->getGraph();
    }

    void bfs(int max)
    {
        alghoritm.get()->bfs(max);
    }

};


int main()
{
    bool fistBuildBRV = true;
    bool fistBuildLUSUGO = true;
    bool fistBuildDCN = true;
 
    int brvMax = 8;
    GMS* gms0 = nullptr;


    int lusugoMax = 1022;
    GMS* gms1 = nullptr;

    int dcnMax = 600;
    GMS* gms3 = nullptr;

    string line;
    while (getline(cin, line)) {


        if (line.empty())
            continue;

        std::istringstream ss(line);
        string model, k, r;
        ss >> model;
        ss >> k;
        if (stoi(model) == 3)
            ss >> r;
        else
            r = "0";

        switch (stoi(model))
        {

        case 0:
        {
            if (fistBuildBRV)
            {
                fistBuildBRV = false;
                gms0 = new GMS(shared_ptr<BRV>(new BRV(brvMax)));
                gms0->buildGraph();
            }
            if (gms0 != nullptr)
                gms0->bfs(stoi(k));
            break;
        }

        case 1:
        {
            if (fistBuildLUSUGO)
            {
                fistBuildLUSUGO = false;
                gms1 = new GMS(shared_ptr<LuSuGuo>(new LuSuGuo(lusugoMax)));
                gms1->buildGraph();
            }

            if (gms1 != nullptr)
                gms1->bfs(stoi(k));
            break;
        }

        case 2:
            break;
        case 3:
        {
            GMS* gms2 = new GMS(shared_ptr<MWI>(new MWI(stoi(k), stoi(r))));
            gms2->buildGraph();
            break;

        }
        case 4:
        {
            if (fistBuildDCN)
            {
                fistBuildDCN = false;
                gms3 = new GMS(shared_ptr<DCN>(new DCN(dcnMax)));
                gms3->buildGraph();
            }

            if (gms3 != nullptr)
                gms3->bfs(stoi(k));
        }

        default:
            break;
        }
    }
}
