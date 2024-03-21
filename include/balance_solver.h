#include <pybind11/pybind11.h>
#include <iostream>
#include <Eigen/Dense>
#include <pybind11/eigen.h>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <list>
#include <ctime>
using namespace std;
namespace py=pybind11;
using namespace Eigen;
#define GravityAcceleration 9.8
//reference https://zhuanlan.zhihu.com/p/441301638
//          Anderson Acceleration for Geometry Optimization and Physics Simulation

#if 0
typedef float scalar;
typedef MatrixX3f MatrixX3s;
typedef VectorXf VectorXs;
typedef Vector3f Vector3s;
typedef Vector2f Vector2s;
typedef MatrixXf MatrixXs;
typedef Matrix3f Matrix3s;
#else
typedef double scalar;
typedef MatrixX3d MatrixX3s;
typedef VectorXd VectorXs;
typedef Vector3d Vector3s;
typedef MatrixXd MatrixXs;
typedef Matrix3d Matrix3s;
#endif

typedef SparseMatrix<scalar> spma;
typedef std::vector<Triplet<scalar>> vtp;

class balance_solver
{
public:
    balance_solver(){}
    ~balance_solver(){}
public:
    void set_info(VectorXi &spring_, VectorXs &stiffness_, VectorXs &v_vertices,
                    VectorXi &rope_index_, VectorXi &stick_index_, 
                    scalar h_interval, scalar m_mass, int a_attach)
    {
        stiffness=stiffness_;
        fixed_stiffness=stiffness(0);
        v=v_vertices; v_new=v; y=v;
        rope_index=rope_index_; stick_index=stick_index_;
        h=h_interval; m=m_mass*3.0/v.size(); 

        
        spring.resize(spring_.size()/2,2);
        for(int i=0;i<spring.rows();++i)
        {
            spring(i,0)=spring_(i*2);
            spring(i,1)=spring_(i*2+1);
        }
        spring_len.resize(spring.rows());
        for(int i=0;i<spring.rows();++i)
        {
            int id0=spring(i,0)*3;
            int id1=spring(i,1)*3;
            spring_len(i)=Vector3s(v(id0)-v(id1),
                        v(id0+1)-v(id1+1),v(id0+2)-v(id1+2)).norm();
        }
        // a.resize(2+stick_index.size());
        // a(0)=a_attach*3; a(1)=a_attach*3+1;
        // for(int i=2;i<stick_index.size()+2;++i)
        // {
        //     a(i)=stick_index(i-2)*3+2;
        // }
        // a.resize(2+stick_index.size()*3);
        // a(0)=a_attach*3; a(1)=a_attach*3+1;
        // for(int i=0;i<stick_index.size();++i)
        // {
        //     for(int j=0;j<3;++j)
        //     {
        //         a(2+i*3+j)=stick_index(i)*3+j;
        //     }
        // }
        a.resize(stick_index.size()*3);
        for(int i=0;i<stick_index.size();++i)
        {
            for(int j=0;j<3;++j)
            {
                a(i*3+j)=stick_index(i)*3+j;
            }
        }
        a_value.resize(a.size());
        for(int i=0;i<a_value.size();++i)
        {
            a_value(i)=v(a(i));
        }

        predecomposition();
        d.resize(J.cols());
        f.resize(v.size());
    }
    void set_forces(VectorXs &f_forces)
    {
        f.setZero();
        for(int i=0;i<v.size()/3;++i)
        {
            f(i*3+2)=-m*GravityAcceleration;
        }
        for(int i=0;i<rope_index.size();++i)
        {
            for(int j=0;j<3;++j)
            {
                f(rope_index(i)*3+j)+=f_forces(i*3+j);
            }
        }
    }

    MatrixX2i spring;
    VectorXs spring_len;
    VectorXs stiffness;
    VectorXs v;
    VectorXs v_new;
    VectorXs y;
    VectorXs f;
    VectorXi rope_index;
    VectorXi stick_index;

    scalar h;
    scalar m;
    VectorXi a;
    VectorXs a_value;
    scalar fixed_stiffness;

    //pd
    SimplicialLDLT<spma> pd_solver;
    VectorXs v_aa;
    VectorXs d;
    spma J;
    std::list<VectorXs> v_list;
    scalar energy;
    long unsigned int anderson_length=7;
    int pd_times=10;
    void predecomposition();
    void local_step(VectorXs &pos);
    void global_step();
    void pd();
    void Anderson();
    void get_energy(VectorXs &pos, scalar &out_energy,bool with_mass=true);

    //newton-raphson
    SparseLU<spma> newton_raphson_solver;
    bool newton_raphson_solver_info=false;
    void newton_raphson();
    bool newton_raphson_linesearch(VectorXs &pos,VectorXs &dir,VectorXs &new_pos);

    //newton
    SimplicialLDLT<spma> newton_solver;
    bool newton_solver_info=false;
    void newton();
    bool newton_linesearch(VectorXs &pos, VectorXs &dir, VectorXs &new_pos);

    //balance
    bool x_invariant=false;
    scalar balance_result;
    scalar compute_balance_times;
    scalar balance_threshold;
    scalar newton_rate;
    void balance_state(VectorXs &pos, VectorXs &force_on_vertex);
    void set_compute_balance_parameter(int compute_balance_times_,scalar balance_threshold_,scalar newton_rate_)
    {
        compute_balance_times=compute_balance_times_;
        balance_threshold=balance_threshold_;
        newton_rate=newton_rate_;
    }
    void compute_balance();

    //jacobi
    VectorXi jacobirow;
    VectorXi jacobicol;
    VectorXs jacobival;
    MatrixXs jacobiright;
    void compute_csr();
    void compute_csr_right();



#pragma region tool functions
    scalar loc_rad;
    scalar dir_rad;
    scalar start_rad;
    scalar end_rad;
    int sampling_num;
    VectorXd sampling_lambda;
    void set_sampling_parameter(scalar loc,scalar dir,scalar start,scalar end,int num)
    {
        loc_rad=loc;dir_rad=dir;start_rad=start;end_rad=end;sampling_num=num;
    }
    void cubic_sampling();

#pragma endregion
};




    
#include <queue>
class stickLooker
{
public:
    stickLooker(){}
    ~stickLooker(){}
public:
    struct Node
    {
        Node(int id_) { prev=nullptr; suc=nullptr; start=nullptr; count=0; id=id_; }
        Node* prev;
        Node* suc;
        Node* start;
        int count;
        int id;
    };
    typedef std::vector<Node> nl;
    nl loop;
    int loopsize;
    Node* entry;
    VectorXi stick_index;

    struct NodeDist
    {
        NodeDist(){}
        NodeDist(int id_,int count_,scalar dist_) { id=id_; count=count_; dist=dist_;}
        int id; int count; scalar dist;
        bool operator>(const NodeDist right) const { return dist>right.dist; }
    };
    std::priority_queue<NodeDist,std::vector<NodeDist>,std::greater<NodeDist>> pq;

    void init(int n)
    {
        for(int i=0;i<n;++i) { loop.emplace_back(i); }
        for(int i=0;i<n-1;++i) { loop[i].suc=&loop[i+1]; loop[i+1].prev=&loop[i];}
        loop[n-1].suc=&loop[0];
        loop[0].prev=&loop[n-1];
        loopsize=n;
        entry=&loop[0];
    }
    void delete_node(int i)
    {
        if(i<0) { std::cout<<"i<0 in node index"<<std::endl; exit(0);}
        Node* node=&loop[i];
        if(entry==node)
        {
            entry=node->suc;
        }
        node->prev->suc=node->suc;
        node->suc->prev=node->prev;
        --loopsize;
    }
    void get_stick_index()
    {
        std::vector<int> idx;
        Node* itr=entry;
        stick_index.resize(loopsize);
        int count=0;
        while(count!=loopsize)
        {
            stick_index(count)=itr->id;
            ++count;
            itr=itr->suc;
        }
    }
    scalar compute_hanging_distance(MatrixX3s &point,Node* node);
    void stick_locating(MatrixX3s &point,int stick_num);
};



PYBIND11_MODULE(balance_solver, m) {
  
  py::class_<balance_solver>(m,"balance_solver")
  .def(py::init<>())
  .def_readwrite("v",&balance_solver::v)
  .def_readwrite("balance_result",&balance_solver::balance_result)
  .def_readwrite("x_invariant",&balance_solver::x_invariant)
  .def_readwrite("jacobirow",&balance_solver::jacobirow)
  .def_readwrite("jacobicol",&balance_solver::jacobicol)
  .def_readwrite("jacobival",&balance_solver::jacobival)
  .def_readwrite("jacobiright",&balance_solver::jacobiright)
  .def_readwrite("sampling_lambda",&balance_solver::sampling_lambda)
  .def("set_info",&balance_solver::set_info)
  .def("set_forces",&balance_solver::set_forces)
  .def("set_compute_balance_parameter",&balance_solver::set_compute_balance_parameter)
  .def("compute_balance",&balance_solver::compute_balance)
  .def("compute_csr",&balance_solver::compute_csr)
  .def("compute_csr_right",&balance_solver::compute_csr_right)
  .def("set_sampling_parameter",&balance_solver::set_sampling_parameter)
  .def("cubic_sampling",&balance_solver::cubic_sampling);
  
  py::class_<stickLooker>(m,"stickLooker")
  .def(py::init<>())
  .def_readwrite("stick_index",&stickLooker::stick_index)
  .def("stick_locating",&stickLooker::stick_locating);
}
//cmake .. -DCMAKE_BUILD_TYPE=Release
//make;rm /home/BA22001030/program/TarpDesign/algorithm/balance_solver.cpython-39-x86_64-linux-gnu.so;mv balance_solver.cpython-39-x86_64-linux-gnu.so /home/BA22001030/program/TarpDesign/algorithm/
//make;rm /home/yanyisheshou/Program/TarpDesign/algorithm/balance_solver.cpython-39-x86_64-linux-gnu.so;mv balance_solver.cpython-39-x86_64-linux-gnu.so /home/yanyisheshou/Program/TarpDesign/algorithm/