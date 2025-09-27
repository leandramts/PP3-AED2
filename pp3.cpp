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
    std::list<VertexWeightPair> *adj;

public:
    WeightedGraphAL(uint num_vertices) : num_vertices(num_vertices),
                                         num_edges(0)
    {
        adj = new std::list<VertexWeightPair>[num_vertices];
    }

    ~WeightedGraphAL()
    {
        delete[] adj;
        adj = nullptr;
    }

    void add_edge(const Vertex &u, const Vertex &v, const Weight &w)
    {
        if (u < 0 || v < 0 || u > num_vertices || v > num_vertices || u == v)
            throw std::invalid_argument("Vertices invalidos");

        auto pair1 = std::make_pair(v, w);

        adj[u].push_back(pair1);

        auto pair2 = std::make_pair(u, w);

        adj[v].push_back(pair2);

        num_edges++;
    }

    const std::list<VertexWeightPair> get_adj(const Vertex &u) const
    {
        if (u < 0 || u > num_vertices)
            throw std::invalid_argument("Vertice invalido");

        return adj[u];
    }

    uint get_num_vertices() const
    {
        return num_vertices;
    }
};

class MinimumPriorityQueue
{
private:
    std::vector<std::pair<Weight, Vertex>> MinHeap;

    int parent(int i)
    {
        return (i - 1) / 2;
    }

    int left(int i)
    {
        return 2 * i + 1;
    }

    int right(int i)
    {
        return 2 * i + 2;
    }

    void min_heapify(int i)
    {
        int l = left(i);
        int r = right(i);
        int smallest = i;

        int tam = MinHeap.size();

        if (l < tam && MinHeap[l].first < MinHeap[i].first) // define qual o menor: o original, esquerda ou direita
            smallest = l;

        if (r < tam && MinHeap[r].first < MinHeap[smallest].first)
        {
            smallest = r;
        }

        if (smallest != i) // se o menor for esquerda ou direita, troca de lugar e recursivamente repete o processo
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

    void build_min_heap() // constroi a heap binaria minima
    {
        int tam = MinHeap.size();
        for (int i = tam / 2 - 1; i >= 0; i--) // repete o processo de encontrar o menor ate a raiz
        {
            min_heapify(i);
        }
    }

    void insert(std::pair<Weight, Vertex> key)
    {
        MinHeap.push_back(key); // insere chaves na heap
        int i = MinHeap.size() - 1;

        while (i > 0 && MinHeap[i].first < MinHeap[parent(i)].first) // enquanto a chave nao for a raiz e for menor que a pai, eles mudam de lugar
        {
            std::pair<Weight, Vertex> aux = MinHeap[i];
            MinHeap[i] = MinHeap[parent(i)];
            MinHeap[parent(i)] = aux;

            i = parent(i);
        }
    }

    std::pair<Weight, Vertex> extract_min() // remove o elemento de menor chave
    {
        if (MinHeap.size() == 0)
        {
            throw std::runtime_error("Fila vazia");
        }

        std::pair<Weight, Vertex> min = MinHeap[0];
        MinHeap[0] = MinHeap.back();
        MinHeap.pop_back(); // vai pro fim da fila pra ser removido

        if (!MinHeap.empty())
            min_heapify(0);

        return min;
    }
};

class AlgorithmDijkstra
{
private:
    const WeightedGraphAL &g;
    std::vector<uint> dist;
    std::vector<int> pred;

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

public:
    AlgorithmDijkstra(const WeightedGraphAL &graph) : g(graph)
    {
        dist.resize(g.get_num_vertices());
        pred.resize(g.get_num_vertices());
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

            for (const auto &adj_vertice : g.get_adj(u))
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
            return path_inverse; // quando o vertice final nao e alcancanvel
        }

        Vertex current_vertex = final; // comecamos o processo pelo vertice final

        while (current_vertex != (uint)-1) // a rota do caminho minimo pelo predecessores
        {
            path_inverse.push_back(current_vertex);
            current_vertex = pred[current_vertex];
        }

        std::vector<Vertex> path;

        for (int i = path_inverse.size() - 1; i >= 0; i--) // botar o caminho na ordem correnta
        {
            path.push_back(path_inverse[i]);
        }

        return path;
    }

    const std::vector<uint> &getDistances() const
    {
        return dist;
    }
};

class ArmiesAttack
{
private:
    struct Army
    {
        // dados de entrada
        std::string color;
        std::string position;
        std::vector<std::string> enemies;

        // dados processados
        std::vector<Vertex> path_to_castle;
        uint turns_to_castle;
        uint min_distance_to_castle;
        std::vector<std::string> allies;
    };
    uint N;            // tamanho do tabuleiro NxN
    WeightedGraphAL g; // grafo que representa o tabuleiro
    std::vector<Army> armies_list;
    std::string castle_position;
    std::vector<std::string> tormenta_positions;

    // funcoes principais
    void armies_paths()
    {
        int L_moves[8][2] =
            {
                {-1, -2}, {1, -2}, {-2, -1}, {2, -1}, // matriz com todas as possibilidades de movimento em L de um exercito
                {-2, 1}, {2, 1}, {-1, 2}, {1, 2},
            };

        for (uint initial_line = 0; initial_line < N; initial_line++) // percorre toda a linha
        {
            for (uint initial_row = 0; initial_row < N; initial_row++) // percorre toda coluna
            {
                Vertex initial_vertice = initial_line * N + initial_row; // atribui um valor para o vertice baseado na linha e coluna

                for (const auto &possible_moves : L_moves)
                {
                    uint new_line = initial_line + possible_moves[1];
                    uint new_row = initial_row + possible_moves[0];

                    if (new_line >= 0 && new_line < N && new_row >= 0 && new_row < N)
                    {
                        Vertex new_vertice = new_line * N + new_row;
                        uint ascii_initial_row = 'a' + initial_row;
                        uint ascii_new_row = 'a' + new_row;

                        Weight w_edge = (ascii_initial_row * (initial_line + 1) + ascii_new_row * (new_line + 1)) % 19; // formula para determinar o peso

                        if (initial_vertice < new_vertice)
                        {
                            g.add_edge(initial_vertice, new_vertice, w_edge);
                        }
                    }
                }
            }
        }
    }

    void process_data()
    {
        for (auto &army : armies_list)
        {
            Vertex a = position_to_vertice(army.position);
            Vertex c = position_to_vertice(castle_position);

            AlgorithmDijkstra ad(g);
            ad.Dijkstra(a);

            uint min_dist = ad.getDistances()[c]; // pega a distancia minima ate o castelo

            army.path_to_castle = ad.get_path(c);
            army.turns_to_castle = army.path_to_castle.size() - 1;
            army.min_distance_to_castle = min_dist;
        }
    }

    void input_data()
    {
        uint num_royal_armies;
        std::cin >> num_royal_armies;
        std::cin.ignore();

        for (uint i = 0; i < num_royal_armies; ++i)
        {
            std::string line;
            std::getline(std::cin, line);
            uint pos = 0;

            std::string color = get_next_token(line, pos);
            std::string position = get_next_token(line, pos);

            std::vector<std::string> enemies;
            std::string enemy;
            while ((enemy = get_next_token(line, pos)) != "")
            {
                enemies.push_back(enemy);
            }

            add_army(color, position, enemies);
        }

        std::cin >> castle_position;

        uint num_tormentas;
        std::cin >> num_tormentas;

        for (uint i = 0; i < num_tormentas; ++i)
        {
            std::string tormenta_pos;
            std::cin >> tormenta_pos;
            tormenta_positions.push_back(tormenta_pos);
        }
    }

    void game_logic()
    {
        int max_turns = 0;
        for (auto &army : armies_list)
        {
            if ((int)army.path_to_castle.size() > max_turns)
                max_turns = army.path_to_castle.size(); //qtd de turno maximo eh a quantidade de passos pra chegar ao castelo
        }

        // simula rodada a rodada
        for (int round = 0; round < max_turns - 1; ++round)
        {
            for (auto &army : armies_list)
            {
                // se o exercito ja chegou ao castelo, ignora
                if (round >= (int)army.path_to_castle.size() - 1)
                    continue;

                // incrementa turns_to_castle se houver inimigos ou tormentas sem alianca, pois esta bloqueado uma rodada
                if (detect_enemies(army, round) || (detect_tormentas(army, round) && !detect_alliance(army, round)))
                {
                    army.turns_to_castle++;
                }

                // atualiza a posicao do army
                army.position = vertice_to_position(army.path_to_castle[round + 1]);
            }
        }
    }

      void find_best_army()
    {
        // encontra o menor numero de turns_to_castle
        uint min_turns = armies_list[0].turns_to_castle;
        for (auto &army : armies_list)
        {
            if (army.turns_to_castle < min_turns)
            {
                min_turns = army.turns_to_castle;
            }
        }

        // cria um vetor com os armies que tem o menor numero de turns_to_castle
        std::vector<ArmiesAttack::Army> best_armies;
        for (auto &army : armies_list)
        {
            if (army.turns_to_castle == min_turns)
            {
                best_armies.push_back(army);
            }
        }

        // ordena em ordem alfabetica (Selection Sort)
        for (uint i = 0; i < best_armies.size(); ++i)
        {
            uint min_idx = i;
            for (uint j = i + 1; j < best_armies.size(); ++j)
            {
                if (best_armies[j].color < best_armies[min_idx].color)
                {
                    min_idx = j;
                }
            }
            if (min_idx != i)
            {
                ArmiesAttack::Army temp = best_armies[i];
                best_armies[i] = best_armies[min_idx];
                best_armies[min_idx] = temp;
            }
        }

        // imprime os melhores armies
        for (auto &army : best_armies)
        {
            std::cout << army.color << " " << army.turns_to_castle << " " << army.min_distance_to_castle << " ";
        }
    }
    // funcoes auxiliares
    void add_army(const std::string &color, const std::string &position, const std::vector<std::string> &enemies)
    {
        Army a;
        a.color = color;
        a.position = position;
        a.enemies = enemies;
        armies_list.push_back(a);
    }

    bool detect_enemies(const Army &current_army, int current_round)
    {
        if (current_round + 1 >= (int)current_army.path_to_castle.size())
            return false; // ja chegou ao castelo

        Vertex next_vertex = current_army.path_to_castle[current_round + 1];
        std:: string next_position = vertice_to_position(next_vertex);
        
        for (auto &army : armies_list)
        {
            for (const auto &enemy_color : current_army.enemies)
            {
                if (army.color == enemy_color && army.position == next_position) //verifica se o proximo movimento tem um inimigo
                {
                    return true;
                }
            }
        }

        return false;
    }

    bool detect_alliance(Army & current_army, int current_round)
    {
        if (current_round >= (int)current_army.path_to_castle.size())
            return false; //ja chegou no castelo

        for (auto &other_army: armies_list)
        {
            if (&other_army == &current_army)
                continue;
            
            bool is_enemy = false;
            for (const auto &e_army: current_army.enemies)
            {
                if(other_army.color == e_army)
                {
                    is_enemy = true;
                    break;
                }
            }

            if (is_enemy)
            {
                continue;
            }

            if (other_army.position == current_army.position) //alia-se a exercito nao inimigo quando estao na mesma posicao
            {
                current_army.allies.push_back(other_army.color);
                other_army.allies.push_back(current_army.color);
                return true;
            }
        }

        return false;

    }

    bool detect_tormentas(const Army &current_army, int current_round)
    {
        if (current_round >= (int)current_army.path_to_castle.size())
                return false;

        for (auto it = tormenta_positions.begin(); it != tormenta_positions.end();)
        {
            if (position_to_vertice(*it) == position_to_vertice(current_army.position)) 
            {
                it = tormenta_positions.erase(it); //remove a tormenta do tabuleiro
                return true;
            }
            else
            {
                ++it;
            }
        }

        return false;
    }

    std::string get_next_token(const std::string &str, uint &pos)
    {
        while (pos < str.length() && str[pos] == ' ')
        {
            pos++;
        }

        uint start = pos;
        while (pos < str.length() && str[pos] != ' ')
        {
            pos++;
        }

        if (start < str.length())
        {
            return str.substr(start, pos - start); // retorna a sub-string que eh a palavra
        }
        return ""; // retorna string vazia se nao houver mais palavras
    }

    Vertex position_to_vertice(std::string position) const
    {
        uint initial_row = position[0] - 'a';
        uint initial_line = std::stoi(position.substr(1)) - 1; // funcao que converte a substring da position a partir do index 1 para inteiro
        return initial_line * N + initial_row;
    }

    std::string vertice_to_position(Vertex v) const
    {
        uint row = v % N;
        uint line = v / N;
        char col_letter = 'a' + row;
        return std::string(1, col_letter) + std::to_string(line + 1);
    }

public:
    ArmiesAttack(uint N) : N(N), g(N * N)
    {
        armies_paths();
    }

    void run()
    {
        input_data();
        process_data();
        game_logic();
        find_best_army();
    }
};

int main()
{
    uint N; // tamanho do tabuleiro NxN
    std::cin >> N;
    ArmiesAttack at(N);
    at.run();

}