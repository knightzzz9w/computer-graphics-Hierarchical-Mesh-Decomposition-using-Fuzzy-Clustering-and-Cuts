#include <utils.h>
#include <iostream>
#include <string>
#include <iomanip> 
using namespace std;

//string filename = "./models/bunny/reconstruction/bun_zipper_res3.ply";
//string filename = "./models/bunny/reconstruction/bun_zipper_res2.ply";
//string filename = "./models/dino.ply";
string filename = "./models/horse.ply";


string header_name = "./standard_header";


int main(int argc, char** argv) {

    

    PlyReader plyReader(filename);

    std::stringstream ply_header ;

    std::vector<Eigen::Vector3d> mesh_vertices = plyReader.getVertices();

    std::vector<std::vector<unsigned int>> mesh_faces = plyReader.getFaces();

    std::cout << "Mesh vertices: " << mesh_vertices.size() << std::endl;
    std::cout << "Mesh faces: " << mesh_faces.size() << std::endl;

    std::ifstream plyFile(filename);
    char delimiter = '/';

    std::vector<std::string> tokens = split(filename, delimiter);
    auto& ply_name =  *(tokens.end() - 1);
    char delimiter2 = '.';
    std::vector<std::string> ply_tokens = split(ply_name, delimiter2);

    std::stringstream stream1 , stream2 , stream3 ;
    stream1 << std::fixed << std::setprecision(2) << fuzzy_thresh;  // 设置为固定的小数点格式，并保留2位小数
    stream2 << std::fixed << std::setprecision(2) << stop_dis_ratio; 
    stream3 << std::fixed << std::setprecision(2) << stop_ang_ratio; 
    std::string str1 = stream1.str();
    std::string str2 = stream2.str();
    std::string str3 = stream3.str();


    auto ply_name2 =  *(ply_tokens.begin()) + "_seg_"+    str1  + "_" + str2  + "_" + str3  + "_" +".ply";
    std::string refine_filename =   ply_name2;
    std::cout << refine_filename << std::endl;

    std::ifstream headerFile(header_name); std::string line;
    while (std::getline(headerFile , line))
    {
        ply_header << line << std::endl;
    }
    


    // for(int i = 0; i < mesh_vertices.size(); i++)
    // {
    //     std::cout << "i is "<<  i <<  "  " << mesh_vertices[i] << std::endl;
    // }

    // for(int i = 0; i < mesh_faces.size(); i++)
    // {
    //     std::cout << mesh_faces[i][0] << mesh_faces[i][1] << mesh_faces[i][2] << std::endl;
    // }


    std::vector<Face> faces;
    for (int i = 0; i < mesh_faces.size(); i++)
    {
        std::vector<unsigned int> cur_vertices_id = mesh_faces[i];
        std::vector<Eigen::Vector3d> vertices_cur;
        for (int j = 0; j < cur_vertices_id.size(); j++)
        {
            vertices_cur.push_back(mesh_vertices[cur_vertices_id[j]]);
        }
        faces.push_back(Face(i , vertices_cur , mesh_faces[i])) ;
    }
    std::cout << "continue 1" << std::endl;

    find_neighbors(mesh_vertices , faces  );

    unsigned int id = 0. ; double avg_geod = 0. ; double avg_ang = 0. ; 
    for (auto &face : faces)
    {
        for (auto &nbr : face.nbrs)
        {
            id += 1 ; 
            avg_geod = (id - 1)*avg_geod/id + nbr.geo_dis/id ;
            avg_ang = (id - 1)*avg_ang/id + nbr.ang_dis/id ;
            //std::cout << "id is " << id << " and cur geo is " << nbr.geo_dis << " and avg geo is " << avg_geod << std::endl ; 
        }
    }

    std::cout << "avg geod is " << avg_geod << " and avg ang is " << avg_ang << std::endl;

    int N_faces = faces.size();

    std::vector<std::vector<int>> neighbor_matrix(N_faces);

    for (auto &face : faces)
    {
        for (auto &nbr : face.nbrs)
        {
            nbr.weight = deta * nbr.geo_dis / avg_geod + (1 - deta) * nbr.ang_dis / avg_ang ;
            neighbor_matrix[face.id].push_back(nbr.face_id); 
        }
            
    }


    //std::cout << std::numeric_limits<double>::max() ; 

    Eigen::MatrixXd graph(N_faces , N_faces);

    graph.setConstant(std::numeric_limits<double>::max());  //先设置为最大值

    init_DistMatrix( faces , graph);

    Eigen::MatrixXd Dist(N_faces , N_faces);  //the final matrix

    Dist.setConstant(std::numeric_limits<double>::max()); 

    globalDijkstra(Dist ,  graph , neighbor_matrix );

    int n = Dist.rows() ;

    // for ( int i = 0 ; i < n ; i++)
    //     for(int j = 0 ; j < n ; j++)
    //         assert(  (std::abs((i,j) - Dist(j,i)) < 1e-4) && !std::isnan(Dist(i,j)))  ;

    //std::vector<int> allREP  = find_allREP(Dist ,  max_iter );

    auto h_section = std::make_pair<double , double >( 0 , 360) ; 

    auto refined_faces = seg(Dist ,   faces, h_section ,  Dist.maxCoeff()) ;

    std::cout << "all refined faces size is " << refined_faces.size()  << std::endl;
    

    //std::cout << "Prob one and two is  " << Prob.row(0) << std::endl << Prob.row(1) <<std::endl;
    // for(auto &id : allREP)
    //     std::cout << id << "," ;

    //find the top 2 prob for each id


    write_ply(ply_header , mesh_vertices , refined_faces , refine_filename);
    return 0;
}