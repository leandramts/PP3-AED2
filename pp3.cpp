
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

    std::list <VertexWeightPair> get_adj(const Vertex& u) const
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
    };

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

    Vertex position_to_vertice(std::string position)
    {
        int initial_row = position[0] - 'a';
        int initial_line = position[1] - '1';

        return initial_line * N + initial_row;
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

        Vertex v = position_to_vertice(position);

        // falta adicionar a logica de introducao do exercito no grafo
    }

    void print_graph() const
    {
        g.print_adjacency_list(g);
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
        return str.substr(start, pos - start); // Retorna a sub-string que é a palavra
    }
    return ""; // Retorna string vazia se não houver mais palavras
}

int main()
{
    int N;
    std:: cin >> N;
    ArmiesAttack at(N);
    at.print_graph();
    
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

    // std::string castle_position;
    // std::cin >> castle_position;

    // int num_tormentas;
    // std::cin >> num_tormentas;

    // std::vector<std::string> tormenta_positions;
    // for (int i = 0; i < num_tormentas; ++i)
    // {
    //     std::string tormenta_pos;
    //     std::cin >> tormenta_pos;
    //     tormenta_positions.push_back(tormenta_pos);
    // }
    
    return 0;
}
