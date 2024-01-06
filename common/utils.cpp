#include "utils.h"


void TicToc::tic() {
    start_time = std::chrono::high_resolution_clock::now();
}

double TicToc::toc(bool display) {

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    if (display) {
        std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
    }

    return elapsed.count();
}


RGB hsv_to_rgb(HSV hsv) {
    RGB rgb;
    double h = hsv.h, s = hsv.s, v = hsv.v;

    double c = v * s;  // 色度
    double x = c * (1 - abs(fmod(h / 60.0, 2) - 1));
    double m = v - c;

    double r_temp, g_temp, b_temp;

    if(h >= 0 && h < 60) {
        r_temp = c, g_temp = x, b_temp = 0;
    } else if(h >= 60 && h < 120) {
        r_temp = x, g_temp = c, b_temp = 0;
    } else if(h >= 120 && h < 180) {
        r_temp = 0, g_temp = c, b_temp = x;
    } else if(h >= 180 && h < 240) {
        r_temp = 0, g_temp = x, b_temp = c;
    } else if(h >= 240 && h < 300) {
        r_temp = x, g_temp = 0, b_temp = c;
    } else {
        r_temp = c, g_temp = 0, b_temp = x;
    }

    rgb.r = r_temp + m;
    rgb.g = g_temp + m;
    rgb.b = b_temp + m;

    return rgb;
}

template <typename T>
std::vector<std::reference_wrapper<T>> extractElementsWithIndices(
    std::vector<T>& originalVector, const std::vector<int>& indices) {

    std::vector<std::reference_wrapper<T>> refVector;

    for (auto index : indices) {
        if (index < originalVector.size()) {
            refVector.push_back(originalVector[index]);
        }
        // 可以选择处理索引越界的情况
    }

    return refVector;
}

// get_DisPoint 函数定义
std::tuple<double, Eigen::Vector3d> get_DisPoint(const Eigen::Vector3d& P, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) {
    Eigen::Vector3d v = v1 - v2;
    Eigen::Vector3d PA = P - v1;
    double t = PA.dot(v) / v.dot(v);
    Eigen::Vector3d Q = v1 + t * v; // 垂点
    double dis = (P - Q).norm();    // 距离

    return {dis, Q}; // 返回距离和垂点
}

// Neighbor 类定义-
Neighbor::Neighbor(int face_id, std::pair<int, int> edge_id, Eigen::Matrix<double, 2, 3> edge, double arc, double dihedral_angle,  double arc_dis, double geo_dis)
    : face_id(face_id), edge_id(edge_id), edge(edge), arc(arc), dihedral_angle(dihedral_angle) , ang_dis(arc_dis), geo_dis(geo_dis), weight(0.0) {}

// Face 类定义
Face::Face(int id, const std::vector<Eigen::Vector3d>& vertices, const std::vector<unsigned int>& vertices_ids)
    : id(id), vertices_ids(vertices_ids) {
    
    center = (vertices[0] + vertices[1] + vertices[2]) / 3.0;

    Eigen::Vector3d init_norm = (vertices[1] - vertices[0]).cross(vertices[2] - vertices[0]);
    norm = init_norm / init_norm.norm();
}

void Face::add_neighbor(Face& neighbor_face, std::pair<int, int> edge_id, const Eigen::Matrix<double , 2 , 3>& edge) {
    int face_id = neighbor_face.id;

    double cos_value = std::max(std::min(norm.dot(neighbor_face.norm), 1.0), -1.0);
    double p_cos_value = std::min( std::abs(norm.dot(neighbor_face.norm)) , 1.);
    double angle = std::acos(cos_value);
    double dihedral_angle = std::acos(p_cos_value);
    bool is_convex = norm.dot(neighbor_face.center - center) < 1e-12;
    double coef = is_convex ? eta : 1.0;

    //double cos_alpha = norm.dot(neighbor_face.norm);
    double ang_dis = coef * (1 - cos_value);

    auto [h0, Point0] = get_DisPoint(center, edge.row(0), edge.row(1));
    auto [h1, Point1] = get_DisPoint(neighbor_face.center, edge.row(0), edge.row(1));
    
    double height = h1 + h0;
    Eigen::Vector3d point_vec = Point1 - Point0;
    double geo_dis = sqrt(height * height + point_vec.dot(point_vec));
    // std::cout << "face id is " << face_id << " and edge id is " << edge_id.first << "," << edge_id.second << 
    // " and ang dis is " << ang_dis << " and geo dis is " << geo_dis << std::endl;
    // std::cout << "h0 is" << h0 <<   "h1 is " << h1 << "neighborcenter is " << neighbor_face.center << " and Point0 is " << Point0 << " and Point1 is " << Point1 << std::endl;
    nbrs.push_back(Neighbor(face_id, edge_id, edge, angle,  dihedral_angle  , ang_dis, geo_dis));
    refine_nbrids.push_back(face_id);
}


void find_neighbors(const std::vector<Eigen::Vector3d>& vertices, std::vector<Face>& faces) {
    std::map<std::pair<int, int>, std::vector<int>> edge_to_faces;

    // 记录每条边对应的面
    for (int face_id = 0; face_id < faces.size(); ++face_id) {
        std::vector<unsigned int> sorted_vertices_id = faces[face_id].vertices_ids;
        std::sort(sorted_vertices_id.begin(), sorted_vertices_id.end());

        for (int i = 0; i < 3; ++i) {
            std::pair<unsigned int, unsigned int> edge_id = std::minmax(sorted_vertices_id[i], sorted_vertices_id[(i + 1) % 3]);
            //std::cout << "edge id is " << edge_id.first << "," <<   edge_id.second << " and  face_id is " << face_id << std::endl;
            edge_to_faces[edge_id].push_back(face_id);
        }
    }

    // 找出邻居面
    for (const auto& [edge_id, face_ids] : edge_to_faces) {
        if (face_ids.size() < 2) {
            continue;
        }

        for (int face_id : face_ids) {
            for (int nbr_face_id : face_ids) {
                if (face_id != nbr_face_id) {
                    Eigen::Matrix<double, 2, 3> edge;
                    edge.row(0) = vertices[edge_id.first];
                    edge.row(1) = vertices[edge_id.second];

                    // 这里假设 Face 类有一个 add_neighbor 方法
                    // 需要根据实际情况来定义 Neighbor 类型及其构造函数
                    faces[face_id].add_neighbor(faces[nbr_face_id] , edge_id , edge);
                    // std::cout << "face id is " << face_id <<" and nbr face id is  " << nbr_face_id << 
                    // " and edge  is " << edge  << std::endl;
                }
            }
        }
    }
}

// 用于优先队列的比较函数

void init_DistMatrix(const std::vector<Face>& faces ,  Eigen::MatrixXd& Dist)
{
    int N  = Dist.rows(); assert(faces.size() == N);

    for (int i = 0; i < N; i++)
    {
        Dist(i , i) = 0.;
        for(auto &nbr: faces[i].nbrs)  //for each nbr
        {
            int nbr_id = nbr.face_id;
            Dist(i , nbr_id) = nbr.weight;
        }
    }

}

// Dijkstra 算法

struct CompareDist {
    bool operator()(const std::pair<int, double>& p1, const std::pair<int, double>& p2) {
        return p1.second > p2.second;
    }
};

void dijkstra(Eigen::MatrixXd &distances, int start , const Eigen::MatrixXd& graph , std::vector<std::vector<int>>& neighbor_matrix) {
    int V = distances.rows();
    std::vector<bool> visited(V, false);
    std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, CompareDist> pq;  //优先是一个最小堆

    distances(start, start) = 0.0;
    pq.push(std::make_pair(start , 0.0)); //初始的节点在S中

    while (!pq.empty()) {
        int u = pq.top().first;
        //std::cout << "cur s u is " << start  << ","  << u << " and dist (s , u) is " << distances(start , u) << std::endl;
        pq.pop();

        if (visited[u]) continue;
        visited[u] = true;

        std::vector<int> cur_nbrid = neighbor_matrix[u];   //这个邻居

        for(auto &v : cur_nbrid) 
        {
            assert(graph(u,v) != std::numeric_limits<double>::max());
            if(!visited[v])
            {
                double new_distance = distances(start , u) + graph(u,v);
                if(new_distance < distances(start , v))
                {
                    distances(start, v) = new_distance;
                    pq.push(std::make_pair(v, new_distance));
                }
            }
        }
    }
}

void globalDijkstra(Eigen::MatrixXd &distances , const Eigen::MatrixXd& graph   , std::vector<std::vector<int>>& neighbor_matrix) {
    int V = distances.rows();
    for (int i = 0; i < V; ++i) {
        dijkstra(distances, i ,  graph , neighbor_matrix);
    }
    std::cout <<    "mean of the distances is  "<<  distances.mean() <<    " and max of the distances is  "  <<    distances.maxCoeff()   <<std::endl;
}



int find_firstREP(const Eigen::MatrixXd &distances) //min distances to all neighbors
{
    int N = distances.rows(); double max_sum = std::numeric_limits<double>::max(); int index = 0;
    for(int i = 0 ; i < N; ++i) 
    {
        double cur_sum = distances.row(i).sum();
        if(cur_sum < max_sum)
        {
            max_sum = cur_sum;
            index = i;
        }

    }
    return index;
}

std::vector<int> find_allREP(const Eigen::MatrixXd &distances ,  int k ) // find k , k stands for the max iter time
{
    int N = distances.rows();

    std::vector<int> allREP ;

    int min_index = find_firstREP(distances);

    allREP.push_back(min_index);

    // find k - 1 rep, each time is the max distance of  the min distance to all rep
    for(int i = 0 ; i < k - 1; ++i)
    {
        double max_min_dis = 0.0; int max_min_index = 0;
        for(int j = 0 ; j < N; ++j)
        {
            double cur_min_dis = std::numeric_limits<double>::max();
            for(auto &rep : allREP)
            {
                cur_min_dis = std::min(cur_min_dis , distances(j , rep));
            }
            if(cur_min_dis > max_min_dis)
            {
                max_min_dis = cur_min_dis;
                max_min_index = j;
            }
        }
        allREP.push_back(max_min_index);
    }

    //find the max derivate of the distance to all rep

    int max_derivate_index = 0; double  last_G = 0.; double max_derivate = 0.0;
    for(int i = 1 ; i < k; ++i)
    {
        double G_i = std::numeric_limits<double>::max();
        for(int j = 0 ; j < i ; j++)
        {
            G_i = std::min(G_i , distances(allREP[i] , allREP[j]));
        }

        if(i == 1)
        {
            last_G = G_i;
            continue; 
        }

        double derivate = last_G - G_i;
        //std::cout << "derivate is " << derivate << std::endl;
        if(derivate > max_derivate)
        {
            max_derivate = derivate;
            max_derivate_index = i - 1;
        }

        last_G = G_i ; 

    }

    return std::vector<int>(allREP.begin() , allREP.begin() + max_derivate_index + 1);

}


Eigen::MatrixXd get_ProbMatrix(const std::vector<int>& allREP ,  const Eigen::MatrixXd& Dist) 
{
    int N = Dist.rows(); int K = allREP.size();

    Eigen::MatrixXd Prob(N , K);

    for(int i = 0 ; i < N; ++i)
    {
        double sum = 0.0;
        auto it  = std::find(allREP.begin() , allREP.end() , i); 
        if(it != allREP.end())  //if i is in the allREP
        {
            int find_place = it - allREP.begin();
            for(int j = 0 ; j < K; ++j)
            {
                if(j == find_place)
                    Prob(i , j) = 1.0;
                else
                    Prob(i , j) = 0.0;
            }
            continue;
        }
        for(int j = 0 ; j < K; ++j)
        {
            Prob(i , j) = 1.0 / Dist(i , allREP[j]);
            sum += Prob(i , j);
        }
        Prob.row(i) /= sum;
    }

    return Prob;
}

 std::vector<Eigen::VectorXd>  get_topprob(const Eigen::MatrixXd& Prob ,
 std::vector<std::vector<int>>&  clear_part)
{
    int N = Prob.rows(); int K = Prob.cols();  
    
    std::vector<std::vector<double>> topProbs_vector(K);

    for(int i = 0 ; i < N ; i++)
    {
        double maxVal = -1;
        int maxIndex = -1;
        std::vector<int> cur_probs;

        for (int j = 0; j < K; ++j) {
            double val = Prob(i, j);
            // 更新最大值和次大值及其索引
            if (val > maxVal) {
                maxVal = val;
                maxIndex = j;
            } 
        }
        //double probDiff = (maxVal - secondMaxVal)/(maxVal + secondMaxVal);
        // std::cout << "continue 2" << std::endl;
        // std::cout << "Prob all"    <<  Prob.row(i)  <<  std::endl;
        // std::cout << "MAX index  and second :"   << maxIndex  <<  "  "  << secondMaxIndex <<  std::endl;
        // std::cout <<   "MAX val  and second :"   <<  maxVal  <<  "  "  << secondMaxVal <<  std::endl;
        clear_part[maxIndex].push_back(i);

       // std::cout << "continue 3" << std::endl;
        topProbs_vector[maxIndex].push_back(maxVal) ;
        //std::cout << "continue 4" << std::endl;

    }

    std::vector<Eigen::VectorXd> topProbs(K);

    for(int i = 0 ; i < K ; i++)
    {
        int partSize = topProbs_vector[i].size();
        topProbs[i] = Eigen::VectorXd(partSize);
        for(int j = 0 ; j < partSize ; j++)
        {
            topProbs[i](j) = topProbs_vector[i][j];
        }
    }

    return topProbs;
    
}

Eigen::MatrixXd extractSubMatrice(const Eigen::MatrixXd& Dist, const std::vector<int>& clear_part) {

    int partSize = clear_part.size();
    if (partSize == 0) return Eigen::MatrixXd(0,0);
    //std::cout << "partSize is " <<  partSize << std::endl;
    // 创建一个新的子矩阵
    Eigen::MatrixXd subMatrix(partSize, partSize);

    // 填充子矩阵
    for (int i = 0; i < partSize; ++i) {
        for (int j = 0; j < partSize; ++j) {
            subMatrix(i, j) = Dist(clear_part[i], clear_part[j]);
        }
    }
    return subMatrix;
}



std::vector<Eigen::MatrixXd> extractSubMatrices(const Eigen::MatrixXd& Dist, const std::vector<std::vector<int>>& clear_part) {
    std::vector<Eigen::MatrixXd> sub_matrices;

    for (const auto& part : clear_part) {
        sub_matrices.push_back(extractSubMatrice(Dist , part));
    }
    return sub_matrices;
}




void refine_REP(const Eigen::MatrixXd& Dist ,  const std::vector<std::vector<int>>& clear_part , 
const std::vector<Eigen::VectorXd>& topprobs , std::vector<int>& allREP)    //all REP是我即将要改变的值
{
    int K = allREP.size();
    auto sub_matrices = extractSubMatrices(Dist , clear_part);
    //find the weighted sub matrices and find the min point
    for(int k = 0 ; k < K ; k++)
    {
        int partSize = clear_part[k].size();
        if (partSize == 0) continue;
        Eigen::MatrixXd subMatrix = sub_matrices[k];
        Eigen::VectorXd topProb = topprobs[k];
        Eigen::VectorXd weightedProb = subMatrix* topProb;
        //std::cout <<"continue 1" << std::endl;
        int minIndex = -1 ;
        double minVal = std::numeric_limits<double>::max();
        for (int i = 0; i < partSize; ++i) {
            if (weightedProb(i) < minVal) {
                minVal = weightedProb(i);
                minIndex = i;
            }
        }
        //std::cout << "min val in "  << k  << " is "  << minVal;
        allREP[k] = clear_part[k][minIndex];
    }

}

void get_top2probandfuzzy(const Eigen::MatrixXd& Prob ,
 std::vector<std::vector<int>>&  clear_part , std::vector<std::vector<std::vector<int>>> & fuzzy_part )
{
    int N = Prob.rows(); int K = Prob.cols();  
    
    for(int i = 0 ; i < N ; i++)
    {
        double maxVal = -1, secondMaxVal = -1;
        int maxIndex = -1, secondMaxIndex = -1;

        for (int j = 0; j < K; ++j) {
            double val = Prob(i, j);
            // 更新最大值和次大值及其索引
            if (val > maxVal) {
                secondMaxVal = maxVal;
                secondMaxIndex = maxIndex;
                maxVal = val;
                maxIndex = j;
            } else if (val > secondMaxVal) {
                secondMaxVal = val;
                secondMaxIndex = j;
            }
        }
        assert(maxIndex != secondMaxIndex) ; 
        double probDiff = (maxVal - secondMaxVal)/(maxVal + secondMaxVal);
        // std::cout << "continue 2" << std::endl;
        // std::cout << "Prob all"    <<  Prob.row(i)  <<  std::endl;
        // std::cout << "MAX index  and second :"   << maxIndex  <<  "  "  << secondMaxIndex <<  std::endl;
        // std::cout <<   "MAX val  and second :"   <<  maxVal  <<  "  "  << secondMaxVal <<  std::endl;
        if (probDiff > fuzzy_thresh) {
                // 明显大于第二个概率，属于一个中心
                clear_part[maxIndex].push_back(i);  //also the relative ones
        }
        else {
                // 属于两个中心的模糊区域
                int first_index = std::min(maxIndex, secondMaxIndex);
                int second_index = std::max(maxIndex, secondMaxIndex);
                fuzzy_part[first_index][second_index].push_back(i);    //注意这里的顺序
        }
       // std::cout << "continue 3" << std::endl;
        //std::cout << "continue 4" << std::endl;

    }
}



bool bfs(const Eigen::MatrixXd &rGraph, int s, int t, std::vector<int> &parent) {
    int V = rGraph.rows();
    std::vector<bool> visited(V, false);
    std::queue<int> queue;

    queue.push(s);
    visited[s] = true;
    parent[s] = -1;

    while (!queue.empty()) {
        int u = queue.front();
        queue.pop();

        for (int v = 0; v < V; v++) {
            if (visited[v] == false && rGraph(u, v) > 0) {
                if (v == t) {
                    parent[v] = u;
                    return true;
                }
                queue.push(v);
                parent[v] = u;
                visited[v] = true;
            }
        }
    }

    return false;
}

double fordFulkerson(Eigen::MatrixXd &graph, int s, int t, std::vector<int> &minCut) {
    int V = graph.rows();
    Eigen::MatrixXd rGraph = graph; // 创建残余网络

    std::vector<int> parent(V); // 存储路径
    double max_flow = 0;

    // 增广路径存在时更新流
    while (bfs(rGraph, s, t, parent)) {
        double path_flow = std::numeric_limits<double>::infinity();
        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            path_flow = std::min(path_flow, rGraph(u, v));
        }

        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            rGraph(u, v) -= path_flow;
            rGraph(v, u) += path_flow;
        }

        max_flow += path_flow;
    }

    // 找出从源点可以到达的所有节点
    for (int i = 0; i < V; i++)
        if (bfs(rGraph, s, i, parent))
            minCut.push_back(i);

    return max_flow;
}

std::pair<std::vector<double> , std::vector< std::pair<double , double> >    > get_allHandpart(std::pair<double , double> h_section , int K)
{
    double h_min = h_section.first; double h_max = h_section.second;
    double h_step = (h_max - h_min) / K ;
    std::vector<double> allH(K);
    std::vector<std::pair<double , double>> allsection(K);
    for(int i = 0 ; i < K ; i++)
    {
        allH[i] = h_min + (i + 0.5)*h_step;
        allsection[i] = std::make_pair(h_min + i*h_step , h_min + (i+1)*h_step);
    }
    return std::make_pair (allH , allsection);
}

bool stop_condition(std::vector<int>& allREP , std::vector<Face>& faces , const Eigen::MatrixXd& Dist , double global_max_dist) //
{
    // if(faces.size() < 3851)
    //     return true;
    double max_Dist = Dist(allREP[0] ,  allREP[1]);
    if(max_Dist < stop_dis_ratio*global_max_dist)
    {
        std::cout <<"meet condition 1    Stop!" << std::endl;
        return true;
    }
        

    unsigned int id = 0. ; double min_dihedral_angle = std::numeric_limits<double>::max(); double max_dihedral_angle = -1.;
    for (auto &face : faces) //update the geo and ang
    {
        for (auto &nbr : face.nbrs)
        {
            id += 1 ; 
            max_dihedral_angle = std::max(max_dihedral_angle , nbr.dihedral_angle);
            min_dihedral_angle = std::min(min_dihedral_angle , nbr.dihedral_angle);
            //std::cout << "id is " << id << " and cur geo is " << nbr.geo_dis << " and avg geo is " << avg_geod << std::endl ; 
        }
    }

    if( max_dihedral_angle - min_dihedral_angle < stop_ang_ratio*max_dihedral_angle )
    {
        std::cout <<"meet condition 2   Stop!" << std::endl;
        return true;
    }
    
    return false;
}


void removeDuplicates(std::vector<int>& vec) {
    // 首先，对向量进行排序
    std::sort(vec.begin(), vec.end());

    // 使用 std::unique 重排向量并删除重复项。
    // 这个函数并不实际删除元素，而是将重复的元素移到向量的末尾，
    // 然后返回一个指向新不重复序列末尾的迭代器。
    auto last = std::unique(vec.begin(), vec.end());

    // 删除向量中所有重复的元素
    vec.erase(last, vec.end());
}

std::vector<Face> seg(const Eigen::MatrixXd& Dist ,  std::vector<Face>& faces  , std::pair<double , double> h_section  , double global_max_dist)   //only seg for once
{  //在face这边就做判断即可
    //find  the top 2 prob that each face belong to

    assert(h_section.first >= 0 && h_section.second <= 360);

    std::vector<int> allREP  = find_allREP(Dist ,  max_iter );  //this is the relative id

    int N = Dist.rows(); int K = allREP.size();  

    if (stop_condition(allREP , faces , Dist , global_max_dist ))
    {
        //add color to faces
        double face_H = (h_section.first + h_section.second)/2;
        HSV face_hsv(face_H ,  const_s , const_v);
        RGB face_color = hsv_to_rgb(face_hsv);
        for(auto& face:faces)
        {
            face.rgb = face_color ; face.hsv = face_hsv;
        }
        return faces;
    }  //constuct a projection from the origin face id to the relative face id

    std::unordered_map<int , int> face_id_to_relative_id;  //the relative id is the id in the allREP

    for(int p = 0; p < faces.size() ; p++)  //consturct  face_id_to_relative_id
    {
        int face_id = faces[p].id;
        face_id_to_relative_id[face_id] = p;

    }

    auto [allH , allsection] = get_allHandpart(h_section , K);

    for(int iter = 0 ; iter < max_iter ; iter++)
    {
        auto Prob =  get_ProbMatrix(allREP ,   Dist) ;
        std::vector<std::vector<int>> clear_part(K); //K vecotr , each contains the label for sure
        //std::cout << "continue 1" << std::endl;
        auto topProbs  = get_topprob(Prob , clear_part);
        // for(int i = 0 ; i < K ; i++)
        //     std::cout << "clear part size " << i << " is " << clear_part[i].size() << std::endl;
        // for(int i = 0 ; i < K ; i++)
        //     for(int j = 0 ; j < K; j++)
        //         std::cout << "fuzzy part size " <<  i  <<   " and " << j << " is " << fuzzy_part[i][j].size() << std::endl;
        
        refine_REP(Dist ,   clear_part , topProbs ,allREP);
    }
    
    auto Prob =  get_ProbMatrix(allREP ,   Dist) ;
    std::vector<std::vector<int>> clear_part(K);
    std::vector<std::vector<std::vector<int>>> fuzzy_part(K , std::vector<std::vector<int>>(K));  //a matrix , each item is the labels of the fuzzy part
    get_top2probandfuzzy(Prob ,  clear_part ,  fuzzy_part )  ;  //还需要有face nbr的信息判断连接性
    
    // for(int i = 0 ; i < K ; i++)   
    // {
    //     std::cout << "clear part size " << i << " is " << clear_part[i].size() << std::endl;
    //     for(int j = i + 1 ; j < K; j++)
    //         std::cout << "fuzzy part size " <<  i  <<   " and " << j << " is " << fuzzy_part[i][j].size() << std::endl;
    // }
    //find the connected neighbors of the fuzzy part
    std::vector<std::vector<std::vector<int>>> fuzzy_part_nbr_left(K , std::vector<std::vector<int>>(K));
    std::vector<std::vector<std::vector<int>>> fuzzy_part_nbr_right(K , std::vector<std::vector<int>>(K));

    for(int i = 0 ; i < K ; i++)      //clear part i
    {
        for(int j = i + 1 ; j < K; j++)  //clear part j 
        {
            std::vector<int>& cur_fuzzy_part = fuzzy_part[i][j];
            //std::cout << "i j is "<< i << " , " << j << " and Part size is " << cur_fuzzy_part.size() << std::endl; 
            std::vector<int>& cur_fuzzy_part_nbr_left = fuzzy_part_nbr_left[i][j] ; 
            std::vector<int>& cur_fuzzy_part_nbr_right = fuzzy_part_nbr_right[i][j] ; 
            for(auto &relative_face_id : cur_fuzzy_part)
            {
                for(auto &relative_nbr_id : faces[relative_face_id].refine_nbrids)
                {    //check if the nbr in the clear part
                    //int nbr_id = nbr.face_id;
                    auto it = std::find(clear_part[i].begin() , clear_part[i].end() , relative_nbr_id);
                    if(it != clear_part[i].end())  //与左边相邻的顶点
                    {
                        cur_fuzzy_part_nbr_left.push_back(relative_nbr_id);
                        continue;
                    }
                    it = std::find(clear_part[j].begin() , clear_part[j].end() , relative_nbr_id);  //与右边相邻的顶点
                    if(it != clear_part[j].end())
                    {
                        cur_fuzzy_part_nbr_right.push_back(relative_nbr_id);
                        continue;
                    }
                }
            }
            //cat cur_fuzzy_part and cur_fuzzy_part_nbr
        }
    }

    for(int i = 0 ; i < K ; i++)      //clear part i
    {
        for(int j = i + 1 ; j < K; j++)  //clear part j 
        {
            std::vector<int>& cur_fuzzy_part_nbr_left = fuzzy_part_nbr_left[i][j] ; 
            std::vector<int>& cur_fuzzy_part_nbr_right = fuzzy_part_nbr_right[i][j] ; 
            removeDuplicates(cur_fuzzy_part_nbr_left);
            removeDuplicates(cur_fuzzy_part_nbr_right);
        }
    }

    std::cout << "Begin new cluster" << std::endl;

    for(int i = 0 ; i < K ; i++)      //clear part i
    {
        for(int j = i + 1 ; j < K; j++)  //clear part j 
        {
            std::vector<int>& cur_fuzzy_part = fuzzy_part[i][j];  //这个是相对应的face的id是多少
            if(cur_fuzzy_part.size() == 0)
                continue;
            std::vector<int>& cur_fuzzy_part_nbr_left = fuzzy_part_nbr_left[i][j] ; 
            std::vector<int>& cur_fuzzy_part_nbr_right = fuzzy_part_nbr_right[i][j] ; 

            std::vector<int> all_part;
            all_part.insert(all_part.end() , cur_fuzzy_part_nbr_left.begin() , cur_fuzzy_part_nbr_left.end());
            all_part.insert(all_part.end() , cur_fuzzy_part.begin() , cur_fuzzy_part.end());
            all_part.insert(all_part.end() , cur_fuzzy_part_nbr_right.begin() , cur_fuzzy_part_nbr_right.end());

            int left_size = cur_fuzzy_part_nbr_left.size();
            int mid_size =  cur_fuzzy_part.size() ;
            int right_size = cur_fuzzy_part_nbr_right.size();

            if(mid_size == 119)
            {
                int xxx = 1;
            }

            int PartSize =  left_size + mid_size + right_size;
            //std::cout << "i j is "<< i << " , " << j << " and Part size is " << PartSize << std::endl; 
            int graph_size = PartSize + 2;
            Eigen::MatrixXd cur_capgraph = Eigen::MatrixXd::Zero(graph_size , graph_size) ;

            double avg_angdis = 0.; int avg_id = 0;
            //find all the edge in EC the fuzzy part
            std::cout << "all  part size is " << all_part.size() << std::endl;
            std::cout << "left size is " << left_size << " and mid size is " << mid_size << " and right size is " << right_size << std::endl;

            for(int k1 = 0 ; k1 < cur_fuzzy_part.size() ; k1++)
            {
                int face_id1 = cur_fuzzy_part[k1];  //relative
                for(auto &relative_nbr_id : faces[face_id1].refine_nbrids )
                {
                    int face_id2 = relative_nbr_id;     //找出邻居节点在哪个位置
                    auto it = std::find(cur_fuzzy_part.begin() , cur_fuzzy_part.end() , face_id2 );
                    if(it != cur_fuzzy_part.end())
                    {
                        avg_id += 1; int matrix_pose1 = 1 + left_size + k1; int matrix_pose2 = 1 + left_size + it - cur_fuzzy_part.begin();
                         double cos_value = std::max(std::min(faces[face_id1].norm.dot(faces[face_id2].norm) , 1.0 ), -1.0);
                        //double angle = std::acos(cos_value);
                        bool is_convex = faces[face_id1].norm.dot(faces[face_id1].center - faces[face_id2].center) < 1e-12;
                        double coef = is_convex ? eta : 1.0;
                        double ang_dis = coef * (1 - cos_value);
                        cur_capgraph(matrix_pose1 , matrix_pose2) = ang_dis;
                        cur_capgraph(matrix_pose2 , matrix_pose1) = ang_dis;
                        avg_angdis = (avg_id - 1)*avg_angdis/avg_id + ang_dis/avg_id;  avg_id += 1;
                        avg_angdis = (avg_id - 1)*avg_angdis/avg_id + ang_dis/avg_id;
                        continue;
                    }

                    auto it2 = std::find(cur_fuzzy_part_nbr_left.begin() , cur_fuzzy_part_nbr_left.end() , face_id2 );
                    if( it2 != cur_fuzzy_part_nbr_left.begin())
                    {
                        avg_id += 1; int matrix_pose1 = 1 + left_size + k1; int matrix_pose2 = 1 + (it2 - cur_fuzzy_part_nbr_left.begin());
                         double cos_value = std::max(std::min(faces[face_id1].norm.dot(faces[face_id2].norm) , 1.0 ), -1.0);
                        //double angle = std::acos(cos_value);
                        bool is_convex = faces[face_id1].norm.dot(faces[face_id1].center - faces[face_id2].center) < 1e-12;
                        double coef = is_convex ? eta : 1.0;
                        double ang_dis = coef * (1 - cos_value);
                        cur_capgraph(matrix_pose1 , matrix_pose2) = ang_dis;
                        cur_capgraph(matrix_pose2 , matrix_pose1) = ang_dis;
                        avg_angdis = (avg_id - 1)*avg_angdis/avg_id + ang_dis/avg_id;  avg_id += 1;
                        avg_angdis = (avg_id - 1)*avg_angdis/avg_id + ang_dis/avg_id;
                        continue;
                    }

                    auto it3 = std::find(cur_fuzzy_part_nbr_right.begin() , cur_fuzzy_part_nbr_right.end() , face_id2 );
                    if( it3 != cur_fuzzy_part_nbr_right.begin())
                    {
                        avg_id += 1; int matrix_pose1 = 1 + left_size + k1; int matrix_pose2 = 1 + left_size + mid_size + (it3 - cur_fuzzy_part_nbr_right.begin());
                         double cos_value = std::max(std::min(faces[face_id1].norm.dot(faces[face_id2].norm) , 1.0 ), -1.0);
                        //double angle = std::acos(cos_value);
                        bool is_convex = faces[face_id1].norm.dot(faces[face_id1].center - faces[face_id2].center) < 1e-12;
                        double coef = is_convex ? eta : 1.0;
                        double ang_dis = coef * (1 - cos_value);
                        cur_capgraph(matrix_pose1 , matrix_pose2) = ang_dis;
                        cur_capgraph(matrix_pose2 , matrix_pose1) = ang_dis;
                        avg_angdis = (avg_id - 1)*avg_angdis/avg_id + ang_dis/avg_id;  avg_id += 1;
                        avg_angdis = (avg_id - 1)*avg_angdis/avg_id + ang_dis/avg_id;
                        continue;
                    }
                }
            }
            std::cout << "cur avg id is " << avg_id << " and avg ang dis is " << avg_angdis << std::endl;

            for(int k1 = 1 ; k1 < graph_size - 1 ; k1++)
            {
                for(int k2 = 1 ; k2 < graph_size - 1; k2++)
                {
                    double cur_value = cur_capgraph(k1 , k2); 
                    assert(cur_value != std::numeric_limits<double>::max());
                    if(  std::abs(cur_value - 0. ) > 1e-12 )
                    {
                        cur_capgraph(k1 , k2) = 1./(1. + cur_value/avg_angdis);
                    }
                }
            }


            for(int k2 = 0 ; k2 < cur_fuzzy_part_nbr_left.size() ; k2++)
            {
                cur_capgraph(0 ,  k2 + 1) = std::numeric_limits<double>::infinity();;   //S to the left part is inf
                cur_capgraph(k2 + 1 ,  0) = std::numeric_limits<double>::infinity();; 
            }

            for(int k3 = 0 ; k3 < cur_fuzzy_part_nbr_right.size() ; k3++)
            {
                int matrix_pose2 = 1 + left_size + mid_size + k3;
                cur_capgraph(graph_size - 1 ,  matrix_pose2) = std::numeric_limits<double>::infinity();;   //T to the right part is inf
                cur_capgraph(matrix_pose2 ,  graph_size - 1) = std::numeric_limits<double>::infinity();; 
            }

            std::vector<int> minCut;
            double max_flow = fordFulkerson(cur_capgraph, 0 ,  graph_size -1 ,  minCut);

            std::cout << "Maximum flow: " << max_flow << std::endl;
            std::cout << "Min cut: ";
            std::sort(minCut.begin() , minCut.end());
            std::for_each(minCut.begin() , minCut.end() , [](int& n){n -= 1;});
            for (int v : minCut)
                std::cout << v << " ";   //
            std::cout << std::endl;
            
            std::vector<int> left_part; std::vector<int> right_part; int cut_id = 0;
            for(int v_id = 0 ; v_id < PartSize ; v_id ++ )
            {
                if(cut_id ==  minCut.size())
                {   
                    int k = v_id ;
                    while(k < PartSize)
                    {
                        right_part.push_back(all_part[k]);
                        k++;
                    }
                    break;
                }
                if(v_id == minCut[cut_id])
                {
                    left_part.push_back(all_part[v_id]);
                    cut_id += 1;
                    continue;
                }
                right_part.push_back(all_part[v_id]);

            }

            // std::cout << "left part size is " << left_part.size() << " right part size is " << right_part.size() << std::endl;
            // assert(left_part.size() + right_part.size() == PartSize);

            cut_id = 0;
            int last_isize = clear_part[i].size(); int last_jsize = clear_part[j].size();
            for(int v_id = 0 ; v_id < PartSize ; v_id ++ )
            {
                if(cut_id ==  minCut.size())
                {   
                    int k = v_id ; 
                    if(k >= left_size + mid_size + right_size )
                        break;
                    while(k < left_size )
                        k++;
                    assert(k >= left_size) ; 
                    while(k < left_size + mid_size)
                    {
                        clear_part[j].push_back(cur_fuzzy_part[k - left_size]);
                        k++;
                    }
                    break;
                }
                if(v_id == minCut[cut_id])
                {
                    if(v_id < left_size + mid_size && v_id >= left_size)
                        clear_part[i].push_back(cur_fuzzy_part[v_id - left_size]);
                    cut_id += 1;
                    continue;
                }
                if(v_id < left_size + mid_size && v_id >= left_size)
                    clear_part[j].push_back(cur_fuzzy_part[v_id - left_size]);
            }
            
            int cur_isize = clear_part[i].size(); int cur_jsize = clear_part[j].size();
            // std::cout << "cur size is" << cur_isize + cur_jsize << std::endl;
            //  std::cout << "last size is" << last_isize + last_jsize << std::endl;
            assert(cur_isize + cur_jsize - last_isize - last_jsize == mid_size);

            // for(int k = 0; k < K ; k++)
            // {

            // }
            
            // for(int p = 0 ; p < faces.size() ; p++)  //origin face id 
            // {
            //     if(std::find(all_part.begin()     ,  all_part.end()  , p) != all_part.end()) //in the fuzzy and neighbour part
            //     {
            //         int relative_face_id = p;
            //         if(std::find(left_part.begin()     ,  left_part.end()  , relative_face_id) != left_part.end())  //face in the left 
            //         {
            //             face.nbrs.erase(std::remove_if(face.nbrs.begin(), face.nbrs.end(), 
            //                   [right_part](const Neighbor& nbr) { return std::find(right_part.begin()     ,  right_part.end()  , nbr.face_id) != right_part.end(); }),
            //    face.nbrs.end());
            //         }
            //         if(std::find(right_part.begin()     ,  right_part.end()  , face_id) != right_part.end())  //face in the left 
            //         {
            //             face.nbrs.erase(std::remove_if(face.nbrs.begin(), face.nbrs.end(), 
            //                   [left_part](const Neighbor& nbr) { return std::find(left_part.begin()     ,  left_part.end()  , nbr.face_id) != left_part.end(); }),
            //    face.nbrs.end());
            //         }
            //     }
            // }
        }
    }

    for(int k = 0; k < K ; k++)
    {
        auto& part_k = clear_part[k];
        for(auto part_id : part_k)  //relative
        {
            faces[part_id].refine_nbrids.erase(std::remove_if(faces[part_id].refine_nbrids.begin(), faces[part_id].refine_nbrids.end(), 
                              [part_k](const int& refine_nbrid) { return std::find(part_k.begin()     ,  part_k.end()  , refine_nbrid) == part_k.end(); }),
                faces[part_id].refine_nbrids.end())  ;

            faces[part_id].nbrs.erase(std::remove_if(faces[part_id].nbrs.begin(), faces[part_id].nbrs.end(), 
                              [part_k , &face_id_to_relative_id](const Neighbor& nbr) { return std::find(part_k.begin()     ,  part_k.end()  , face_id_to_relative_id[nbr.face_id]) == part_k.end(); }),
                faces[part_id].nbrs.end())  ;
            assert(faces[part_id].nbrs.size() == faces[part_id].refine_nbrids.size());
        }

    }

    std::vector<Face> refine_faces;

    std::cout <<"part size is " << K << std::endl;

    if(clear_part[0].size() == 51 && clear_part[1].size() == 22)
    {
        int xxx = 1;
    }

    for(int k = 0 ; k < K; k ++)
    {
        std::cout << "clear part size " << k << " is " << clear_part[k].size() << std::endl;
    }

    for(int k = 0 ; k < K ; k++)    //construct a map that project the face id to the relative id
    {
        //std::cout << "clear part size " << k << " is " << clear_part[k].size() << std::endl;
        std::vector<Face> cur_segfaces;
        std::unordered_map<int , int> cur_face_id_to_relative_id;
        for(int i = 0 ; i < clear_part[k].size() ; i++)
        {
            int face_id = clear_part[k][i];
            cur_face_id_to_relative_id[face_id] = i;
        }
        for(int i = 0 ; i < clear_part[k].size() ; i++)
        {
            int face_id = clear_part[k][i];
            for (int j = 0; j < faces[face_id].refine_nbrids.size() ; j++)
            {
                faces[face_id].refine_nbrids[j]  =  cur_face_id_to_relative_id[faces[face_id].refine_nbrids[j]];
            }
            cur_segfaces.push_back(faces[face_id]);
        }
        auto submatrix = extractSubMatrice(Dist , clear_part[k]);
        auto cur_refine_faces = seg(submatrix , cur_segfaces , allsection[k] , global_max_dist);
        refine_faces.insert(refine_faces.end() , cur_refine_faces.begin() , cur_refine_faces.end());
    }

    return refine_faces;
}



void write_ply(std::stringstream& header_ss, std::vector<Eigen::Vector3d>& mesh_vertices, std::vector<Face>& mesh_faces, std::string filename) {
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }

    // 写入头部信息
    std::string line;
    while (std::getline(header_ss, line)) {
        if (line.find("element vertex") != std::string::npos) {
            outFile << "element vertex " << mesh_vertices.size() << std::endl;
        } else if (line.find("element face") != std::string::npos) {
            outFile << "element face " << mesh_faces.size() << std::endl;
        } else {
            outFile << line << std::endl;
        }
    }

    // 写入顶点数据
    for (const auto& vertex : mesh_vertices) {
        outFile << vertex[0] << " " << vertex[1] << " " << vertex[2] << std::endl;
    }

    // 写入面数据
    for (const auto& face : mesh_faces) {
        outFile << face.vertices_ids.size();
        for (auto id : face.vertices_ids) {
            outFile << " " << id;
        }

        // 写入 RGB 数据
        int r = static_cast<int>(face.rgb.r * 255);
        int g = static_cast<int>(face.rgb.g * 255);
        int b = static_cast<int>(face.rgb.b * 255);
        outFile << " " << r << " " << g << " " << b << std::endl;
    }

    outFile.close();
}






