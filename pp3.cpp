
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

//tirar duvida de começar em 1 conforme os slides ou 0 conforme std::vector
//esta comecando pelo indice 0 atualmente
class MinHeap
{
    private:
        std::vector<Weight> heap;

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
    
    public:
        MinHeap() {}
    
    std::vector<Weight> getHeap()
    {
        return heap;
    }


    void min_heapify(int i)
    {
        int l = left(i);
        int r = right(i);
        int smallest;

        int tam = heap.size();

        if (l < tam && heap[l] < heap[i])
            smallest = l;
        else
        {
            smallest = i;
        }

        if (r < tam && heap[r] < heap[smallest])
        {
            smallest = r;
        }

        if (smallest!=i)
        {
            Weight aux;
            aux = heap[i];
            heap[i] = heap[smallest];
            heap[smallest] = aux;
            min_heapify(smallest);
        }
    }

    void build_min_heap()
    {
        int tam = heap.size();
        for (int i = tam/2-1; i>=0; i--)
        {
            min_heapify(i);
        }
    }

    void insert(Weight key)
    {
        heap.push_back(key);
        int i = heap.size() - 1;

        while (i>0 and heap[i]< heap[parent(i)])
        {
            Weight aux = heap[i];
            heap[i] = heap[parent(i)];
            heap[parent(i)] = aux;
            i = parent(i);

        }
    }

    Weight extract_min()
    {
        if (heap.size() == 0)
        {
            throw std::runtime_error("Fila vazia");
        }

        Weight min = heap[0];
        heap[0] = heap.back();
        heap.pop_back();
        min_heapify(0);

        return min;
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

    void print_graph() const
    {
        g.print_adjacency_list(g);
    }

};

// int main()
// {
//     int N;
//     std:: cin >> N;
//     ArmiesAttack at(N);
//     //at.print_graph();
    
//     return 0;
// }



// int main() {
//     MinHeap heap;

//     // Teste 1: Inserções
//     std::cout << "Inserindo elementos:\n";
//     heap.insert(5);
//     heap.insert(3);
//     heap.insert(8);
//     heap.insert(1);
//     heap.insert(6);
//     heap.insert(10);
//     heap.insert(0);
//     heap.insert(11);
//     heap.insert(7);

//     // Mostrar conteúdo interno do heap
//     std::cout << "Heap interno apos insercoes: ";
//     // (precisa de um método para mostrar, por exemplo getHeap())
//     for (auto x : heap.getHeap()) {
//         std::cout << x << " ";
//     }
//     std::cout << "\n";

//     // Teste 2: Extração do mínimo
//     std::cout << "Extraindo min: " << heap.extract_min() << "\n";
//     std::cout << "Heap apos extract_min: ";
//     for (auto x : heap.getHeap()) {
//         std::cout << x << " ";
//     }
//     std::cout << "\n";

//     // Teste 3: Várias extrações até esvaziar
//     std::cout << "Extraindo todos:\n";
//     while (true) {
//         try {
//             std::cout << heap.extract_min() << " ";
//         } catch (std::runtime_error& e) {
//             std::cout << "\n[Fila vazia]\n";
//             break;
//         }
//     }

//     return 0;
// }