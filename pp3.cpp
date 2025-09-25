
#include <iostream>
#include <list>
#include <vector>
#include <stdexcept>
#include <limits>
#include <utility>
#include <string>

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

            if (l < tam && MinHeap[l].first < MinHeap[i].first) //define qual o menor: o original, esquerda ou direita
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
        const WeightedGraphAL& g;
        std::vector<uint> dist;
        std::vector<int> pred;
    
    public:
         AlgorithmDijkstra(const WeightedGraphAL& graph) : 
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
      
    std::vector<Vertex> get_path(Vertex final) const
    {
        std::vector<Vertex> path_inverse;

        if (dist[final] == std::numeric_limits<uint>::max())
        {
            return path_inverse; //quando o vertice final nao é alcancanvel
        }

        Vertex current_vertex = final; //comecamos o processo pelo vertice final

        while (current_vertex != -1) //a rota do caminho minimo pelo predecessores
        {
            path_inverse.push_back(current_vertex);
            current_vertex = pred[current_vertex];
        }

        std::vector <Vertex> path;

        for (int i = path_inverse.size() - 1; i >= 0; i--) //botar o caminho na ordem correnta
        {
            path.push_back(path_inverse[i]);
        }

        return path;
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

    struct Army
    {
        std::string cor;
        std::string posicao;
        std::vector<std::string> inimigos;
        int rodadas;
    };

    std::vector<Army> armies_list;

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
                        int ascii_initial_row = 'a' + initial_row; //codigo ascii comecando por "a" -> 97
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

    void add_army(const std::string& color, const std::string& position, const std::vector<std::string>& enemies)
    {
        Army a;
        a.cor = color;
        a.posicao = position;
        a.inimigos = enemies;
        armies_list.push_back(a);

         // falta adicionar a logica de introducao do exercito no grafo

    }

    const std::vector<Army>& get_armies() const
    {
        return armies_list;
    }

    uint army_djk_distance_to_castle(std::string army_position, std::string castle_position)
    {
        Vertex a = position_to_vertice(army_position);
        Vertex c = position_to_vertice(castle_position);
        AlgorithmDijkstra ad(g);
        ad.Dijkstra(a); //aplica o algoritmo de Dijkstra com o vertice do exercito sendo a origem
        const auto& array_dist = ad.getDistances(); //pega os vetores das distancias minimas
        uint min_dist = array_dist[c];

        if (min_dist == std::numeric_limits<uint>::max())
        {
            throw std::runtime_error("Nao tem caminho ate castelo");
        }

        return min_dist; //retorna a distancia minima segundo o algoritmo ate o castelo
    }
  
    void print_graph() const //funcao pra teste de depuracao 
    {
        g.print_adjacency_list(g);
    }

    const WeightedGraphAL& getGraph() const //funcao pra teste de depuracao 
    {
        return g;
    }

     Vertex position_to_vertice(std::string position) const
    {
        int initial_row = position[0] - 'a';
        int initial_line = std::stoi(position.substr(1)) - 1;  //funcao que converte a substring da position a partir do index 1 para inteiro
        return initial_line * N + initial_row;
    }

    bool detect_enemies(Vertex next_vertex, const Army& current_army) const
    {
            for (const auto& enemy_color: current_army.inimigos)
            {
                for (const auto& other_army: armies_list)
                {
                if (other_army.cor == enemy_color)
                {
                    Vertex enemy_vertex = position_to_vertice(other_army.posicao);
    
                    if (next_vertex == enemy_vertex)
                    {
                        return true;
                    }
                }
        
                }
            }
             return false;
    }

};

std::string get_next_token(const std::string& str, int& pos) {
    while (pos < str.length() && str[pos] == ' ') {
        pos++;
    }

    int start = pos;
    while (pos < str.length() && str[pos] != ' ') {
        pos++;
    }

    if (start < str.length()) {
        return str.substr(start, pos - start); // Retorna a sub-string que eh a palavra
    }
    return ""; // Retorna string vazia se nao houver mais palavras
}

int main()
{
    int N;
    std:: cin >> N;
    ArmiesAttack at(N);
   // at.print_graph();
    
    int num_royal_armies;
    std::cin >> num_royal_armies;
    std::cin.ignore();

    for (int i = 0; i < num_royal_armies; ++i)
    {
        std::string line;
        std::getline(std::cin, line);
        int pos = 0;

        std::string color = get_next_token(line, pos);
        std::string position = get_next_token(line, pos);

        std::vector<std::string> enemies;
        std::string enemy;
        while ((enemy = get_next_token(line, pos)) != "")
        {
            enemies.push_back(enemy);
        }
        
        at.add_army(color, position, enemies);

        std::cout << "Inimigos:";
        for (auto &e : enemies)
        {
            std::cout << " " << e;
        }
        std::cout << "\n";
        
    }

    std::string castle_position;
    std::cin >> castle_position;

    int num_tormentas;
    std::cin >> num_tormentas;

    std::vector<std::string> tormenta_positions;
    for (int i = 0; i < num_tormentas; ++i)
    {
        std::string tormenta_pos;
        std::cin >> tormenta_pos;
        tormenta_positions.push_back(tormenta_pos);
    }
    
    //testes de depuracao 
    const auto& all_armies = at.get_armies();


    for (const auto& army : all_armies) 
    {
        try
        {

            uint distance = at.army_djk_distance_to_castle(army.posicao, castle_position);
            std::cout << "Exercito " << army.cor << " (" << army.posicao << "): "
                      << "Distancia minima = " << distance << std::endl;
        }
        catch (const std::runtime_error& e)
        {
            std::cout << "Exercito " << army.cor << " (" << army.posicao << "): "
                      << "INACESSIVEL" << std::endl;
        }
    }



}
