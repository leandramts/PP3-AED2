
#include <iostream>
#include <list>
#include <vector>
#include <stdexcept>
#include <limits>
#include <utility> 

using uint = unsigned int;
using Vertex = unsigned int;
using Weight = unsigned int;
using VertexWeightPair = std::pair<Vertex, Weight>;

class WeightedGraphAL 
{
    private:
    uint num_vertices;
    uint num_edges;
    std::list <VertexWeightPair> *adj;

    public:
    WeightedGraphAL(uint num_vertices):
    num_vertices(num_vertices),
    num_edges(0)
    {
        adj = new std::list<VertexWeightPair>[num_vertices];
    }
    ~WeightedGraphAL()
    {
        delete []adj;
        adj = nullptr;
    }

    void add_edge(const Vertex& u, const Vertex& v, const Weight& w)
    {
        if (u < 0 || v < 0 || u > num_vertices || v > num_vertices || u == v)
        throw std::invalid_argument("Vertices invalidos");
        
            auto pair1 = std::make_pair(v,w);
            
            adj[u].push_back(pair1);

            auto pair2 = std::make_pair(u,w);

            adj[v].push_back(pair2);

            num_edges++;

    }

    const std::list <VertexWeightPair> get_adj(const Vertex& u) const
    {
        if ( u < 0 || u > num_vertices)
        throw std:: invalid_argument("Vertice invalido");

        return adj[u];
    }

    uint get_num_vertices ( ) const
    {
        return num_vertices;
    }

    uint get_num_edges() const
    {
        return num_edges;
    }

    void print_adjacency_list(const WeightedGraphAL& g) const
    {
        uint n = g.get_num_vertices();
        uint e = g.get_num_edges();

        std:: cout << "num_vertices: " << n << std::endl;
        std:: cout << "num_edges: " << e << std::endl;

        for (uint u = 0; u < n; u++)
        {
            std:: cout << u << ": ";
            std:: list<VertexWeightPair> l = get_adj(u);
            for(const auto& item: l)
            {
                std:: cout << "(" << item.first << ", " << item.second << "), ";
            }
            std:: cout <<std::endl;
        }
    }
};

class MinimumPriorityQueue
{
    private:
        std::vector<std::pair<Weight, Vertex>> MinHeap;

        int parent(int i) 
        {
            return (i-1)/2;
        }

        int left(int i)
        {
            return 2*i+1;
        }

        int right(int i)
        {
            return 2*i + 2;
        }

        void min_heapify(int i)
        {
            int l = left(i);
            int r = right(i);
            int smallest = i; 

            int tam = MinHeap.size();

            if (l < tam && MinHeap[l].first < MinHeap[i].first) //define se o menor é o original, esquerda ou direita
                smallest = l; 

            if (r < tam && MinHeap[r].first < MinHeap[smallest].first) 
            {
                smallest = r;
            }

            if (smallest!=i) //se o menor for esquerda ou direita, troca de lugar e recursivamente repete o processo
            {
                std::pair<Weight, Vertex> aux = MinHeap[i];
                MinHeap[i] = MinHeap[smallest];
                MinHeap[smallest] = aux;
                min_heapify(smallest);
            }
        }
    
    public:
        MinimumPriorityQueue() {}
   
    bool empty()
    {
        return MinHeap.empty();
    }

    std::vector<std::pair<Weight, Vertex>> getHeap() { return MinHeap; }


    void build_min_heap() //constroi a heap binaria minima 
    {
        int tam = MinHeap.size();
        for (int i = tam/2-1; i>=0; i--) //repete o processo de encontrar o menor ate a raiz
        {
            min_heapify(i);
        }
    }

    void insert(std::pair<Weight, Vertex> key)
    {
        MinHeap.push_back(key); //insere chaves na heap
        int i = MinHeap.size() - 1;

        while (i>0 && MinHeap[i].first< MinHeap[parent(i)].first) //enquanto a chave nao for a raiz e for menor que a pai, eles mudam de lugar
        {
            std::pair<Weight, Vertex> aux = MinHeap[i];
            MinHeap[i] = MinHeap[parent(i)];
            MinHeap[parent(i)] = aux;

            i = parent(i);

        }
    }

    std::pair<Weight, Vertex> extract_min() //remove o elemento de menor chave
    {
        if (MinHeap.size() == 0)
        {
            throw std::runtime_error("Fila vazia");
        }

        std::pair<Weight, Vertex> min = MinHeap[0];
        MinHeap[0] = MinHeap.back();
        MinHeap.pop_back(); //vai pro fim da fila pra ser removido
        
        if (!MinHeap.empty()) 
            min_heapify(0);

        return min;
    }
};

class AlgorithmDijkstra
{
    private:
        WeightedGraphAL g;
        std::vector<uint> dist;
        std::vector<int> pred;
    
    public:
         AlgorithmDijkstra(WeightedGraphAL& graph) : 
         g(graph) 

        {
            dist.resize(g.get_num_vertices());
            pred.resize(g.get_num_vertices());
        }

        void initialize(Vertex s)

        {
            uint n = g.get_num_vertices();
            for (Vertex u = 0; u < n; u++)
            {
                if (u == s)
                {
                    dist[u] = 0;
                    pred[s] = -1;
                }
                else
                {
                    dist[u] = std::numeric_limits<uint>::max(); // considera primeiramente distancia infinita
                    pred[u] = -1; 
                }
            }

        }

        void relax(Vertex u, Vertex v, Weight w)
        {

            if (dist[v] > dist[u] + w)
            {
                dist[v] = dist[u] + w;
                pred[v] = u;
            }

        }

        void Dijkstra(Vertex s)
        {
            initialize(s);
            MinimumPriorityQueue mq;
            uint n = g.get_num_vertices();

            for (Vertex u = 0; u < n; u++)
            {
                mq.insert({dist[u], u});
            }

            while (!mq.empty())
            {
                auto [d, u] = mq.extract_min();

                for (const auto& adj_vertice: g.get_adj(u))
                {
                    Vertex v = adj_vertice.first;
                    Weight w = adj_vertice.second;
                    if (dist[v] > dist[u] + w)
                    {
                        relax(u, v, w);
                        mq.insert({dist[v], v});
                    }

                }
            }

        }

        //funcoes para testes de depuracoes (retirar depois)

        void print_distances() 
        {
            for (uint u = 0; u < g.get_num_vertices(); u++) 
            {
                std::cout << "dist[" << u << "] = ";
                if (dist[u] == std::numeric_limits<uint>::max())
                    std::cout << "INF\n";
                else
                    std::cout << dist[u] << "\n";
            }
        }

        const std::vector<uint>& getDistances() const 
        {
        return dist;
        }
};


class ArmiesAttack
{
private:
    WeightedGraphAL g;
    int N;


    void armies_paths(WeightedGraphAL& g, int N)
    { 
    int L_moves[8][2] =
    { {-1,-2}, {1, -2}, {-2,-1}, {2, -1},  //matriz com todas as possibilidades de movimento em L de um exercito 
        {-2,1}, {2,1}, {-1,2}, {1, 2}, 
    };

        for (int initial_line = 0; initial_line < N; initial_line++) //percorre toda a linha
        {
            for (int initial_row = 0; initial_row < N; initial_row++) // percorre toda coluna
            {
                Vertex initial_vertice = initial_line * N + initial_row; //atribui um valor para o vertice baseado na linha e coluna

                for (const auto& possible_moves : L_moves)   
                {
                    int new_line = initial_line + possible_moves[1];
                    int new_row  = initial_row + possible_moves[0];

                    if (new_line >= 0 && new_line < N && new_row >= 0 && new_row < N) 
                    {
                        Vertex new_vertice = new_line * N + new_row;
                        int ascii_initial_row = 'a' + initial_row; //codigo ascii começando por "a" -> 97
                        int ascii_new_row = 'a' + new_row;

                        Weight w_edge = (ascii_initial_row*(initial_line+1) + ascii_new_row*(new_line+1)) % 19; //formula para determinar o peso

                        if (initial_vertice < new_vertice) 
                        {
                            g.add_edge(initial_vertice, new_vertice, w_edge);
                        }
                    }
                }
            }
        }
    }


public:
    ArmiesAttack(int N) : g(N*N), N(N)
    {
        armies_paths(g, N);
    }

    void print_graph() const //funcao pra teste de depuracao 
    {
        g.print_adjacency_list(g);
    }

    WeightedGraphAL getGraph() const //funcao pra teste de depuracao 
    {
        return g;
    }

     Vertex position_to_vertice(std::string position)
    {
        int initial_row = position[0] - 'a';
        int initial_line = position[1] - '1';

        return initial_line * N + initial_row;
    }

};

int main() {

    //testes do AlgorithmDijkstra (aparentemente coerente com os testes)

    // int N;
    // std::cout << "Digite o tamanho do tabuleiro NxN (8 a 15): ";
    // std::cin >> N;

    // ArmiesAttack battle(N);
    // WeightedGraphAL g_battle = battle.getGraph();

    // std::string army_color, army_pos;
    // std::cout << "Digite a cor do exército: ";
    // std::cin >> army_color;
    // std::cout << "Digite a posição inicial do exército (ex: b8): ";
    // std::cin >> army_pos;

    // Vertex start = battle.position_to_vertice(army_pos);

    // std::string castle_pos;
    // std::cout << "Digite a posição do castelo de Hunnus: ";
    // std::cin >> castle_pos;
    // Vertex castle = battle.position_to_vertice(castle_pos);

    // Algorithm_Dijkstra dijkstra(g_battle);
    // dijkstra.Dijkstra(start);

    // const auto& dist = dijkstra.getDistances(); // precisa do getter em Algorithm_Dijkstra

    // std::cout << "\nExército " << army_color 
    //           << " até o castelo " << castle_pos << ": ";

    // if (dist[castle] == std::numeric_limits<uint>::max())
    //     std::cout << "INACESSÍVEL\n";
    // else
    //     std::cout << "distância mínima = " << dist[castle] << "\n";

    return 0;
}
