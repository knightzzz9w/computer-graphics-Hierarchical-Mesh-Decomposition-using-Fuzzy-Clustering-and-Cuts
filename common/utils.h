#ifndef UTILS_H
#define UTILS_H

#include <chrono>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <Eigen/Dense>
#include <chrono>
#include <unordered_map>
#include <utility>
#include <map>
#include <cmath>
#include <queue>


// 声明全局变量
extern double deta = 0.8;
extern double eta = 0.2;
extern int max_iter = 10;
extern double  fuzzy_thresh = 0.04;
extern double  stop_dis_ratio = 0.03;
extern double  stop_ang_ratio = 0.05;
extern double  const_s = 0.7;
extern double  const_v = 0.8;

std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);

    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }

    return tokens;
}



// TicToc 类声明
class TicToc {
public:
    TicToc();
    void tic();
    double toc(bool display = false);

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
};

class PlyReader {
public:
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::vector<unsigned int>> faces;
    int add_property = 0;

    PlyReader(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error opening file." << std::endl;
            return;
        }

        unsigned int numVertices = 0, numFaces = 0;
        readHeader(file, numVertices, numFaces);
        readVertices(file, numVertices);
        readFaces(file, numFaces);

        file.close();
    }

    std::vector<Eigen::Vector3d> getVertices() const {
        return vertices;
    }

    std::vector<std::vector<unsigned int>> getFaces() const {
        return faces;
    }

private:
   
    void readHeader(std::ifstream& file, unsigned int &numVertices, unsigned int &numFaces) {
        std::string line; bool meet_face = false; int property_num = 0;
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string word;
            ss >> word;
            if (word == "element") {
                ss >> word;
                if (word == "vertex") {
                    ss >> numVertices;
                } else if (word == "face") {
                    meet_face = true;
                    ss >> numFaces;
                }
            } 
            else if (word == "property") {
                if(!meet_face)
                {
                    property_num += 1;
                    continue;
                }
            }
            else if (line == "end_header") {
                break;
            }
        }
        add_property = property_num -3;
    }

    void readVertices(std::ifstream& file, unsigned int numVertices) {
        std::vector<double> add(add_property);
        for (unsigned int i = 0; i < numVertices; ++i) {
            Eigen::Vector3d v;
            file >> v[0] >> v[1] >> v[2];
            for(int j = 0 ; j < add_property ; j++)
            {
                file >> add[j];
            }
            vertices.push_back(v);
        }
    }

    void readFaces(std::ifstream& file, unsigned int numFaces) {
        for (unsigned int i = 0; i < numFaces; ++i) {
            unsigned int n;
            file >> n;
            std::vector<unsigned int> face;
            unsigned int index;
            for (unsigned int j = 0; j < n; ++j) {
                file >> index;
                face.push_back(index);
            }
            faces.push_back(face);
        }
    }
};

struct RGB {
    double r;  // 范围: [0, 1]
    double g;  // 范围: [0, 1]
    double b;  // 范围: [0, 1]
};

struct HSV {
    double h;  // 范围: [0, 360]
    double s;  // 范围: [0, 1]
    double v;  // 范围: [0, 1]

    HSV(double hue, double saturation, double value) 
        : h(hue), s(saturation), v(value) {
        // Ensure the values are within the appropriate ranges
        h = std::max(0.0, std::min(h, 360.0));
        s = std::max(0.0, std::min(s, 1.0));
        v = std::max(0.0, std::min(v, 1.0));
    }

    HSV() 
        : h(0.), s(0.), v(0.) {
    }
};

RGB hsv_to_rgb(HSV hsv);

template <typename T>
std::vector<std::reference_wrapper<T>> extractElementsWithIndices(
    std::vector<T>& originalVector, const std::vector<int>& indices);


class Neighbor {
public:
    Neighbor(int face_id, std::pair<int, int> edge_id, Eigen::Matrix<double, 2, 3> edge, double arc,  double dihedral_angle , double arc_dis, double geo_dis);

    int face_id;
    std::pair<int, int> edge_id;
    Eigen::Matrix<double, 2, 3> edge;
    double arc;
    double dihedral_angle;
    double ang_dis;
    double geo_dis;
    double weight;
};

// Face 类声明
class Face {
public:
    Face(int id, const std::vector<Eigen::Vector3d>& vertices, const std::vector<unsigned int>& vertices_ids);

    void add_neighbor(Face& neighbor_face, std::pair<int, int> edge_id, const Eigen::Matrix<double, 2, 3>& edge);
    int id;
    std::vector<unsigned int> vertices_ids;
    Eigen::Vector3d center;
    Eigen::Vector3d norm;
    std::vector<Neighbor> nbrs;
    std::vector<int> refine_nbrids;
    RGB rgb ; 
    HSV hsv;
};

// 函数声明
std::tuple<double, Eigen::Vector3d> get_DisPoint(const Eigen::Vector3d& P, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2);
void find_neighbors(const std::vector<Eigen::Vector3d>& vertices, std::vector<Face>& faces);

void init_DistMatrix(const std::vector<Face>& faces ,  Eigen::MatrixXd& Dist);

void globalDijkstra(Eigen::MatrixXd &distances , const Eigen::MatrixXd &graph   , std::vector<std::vector<int>>& neighbor_matrix) ;

std::vector<int> find_allREP(const Eigen::MatrixXd &distances ,  int k ) ;  

Eigen::MatrixXd get_ProbMatrix(const std::vector<int>& allREP ,  const Eigen::MatrixXd& Dist) ;

std::vector<Face> seg(const Eigen::MatrixXd& Dist ,   std::vector<Face>& faces ,  std::pair<double , double> h_section  , double global_max_dist);

void write_ply(std::stringstream& header_ss , std::vector<Eigen::Vector3d> & mesh_vertices , std::vector<Face>&  mesh_faces , std::string filename);


#endif // GEOMETRY_H
