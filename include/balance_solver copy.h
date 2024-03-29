#include <pybind11/pybind11.h>
#include <iostream>
#include <Eigen/Dense>
#include <pybind11/eigen.h>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
//#include <Eigen/IterativeLinearSolvers>
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

//typedef Matrix<scalar, 3, 3, 0, 3 ,3> Matrix3s;
typedef SparseMatrix<scalar> spma;
typedef std::vector<Triplet<scalar>> vtp;

class Constraint
{
public:
    Constraint(scalar stiffness):s(stiffness){}
    virtual ~Constraint(){}

    virtual void  EvaluateEnergy(VectorXs& x, scalar &energy) {}
    virtual void  EvaluateGradient(VectorXs& x, VectorXs& gradient) {}
    virtual void  EvaluateHessian(VectorXs& x, vtp& hessian_triplets) {}
    virtual void  EvaluateWeightedLaplacian(vtp& laplacian_triplets) {}

    // for pd
    virtual void EvaluateDVector(unsigned int index, VectorXs& x, VectorXs& d){}
    virtual void EvaluateJMatrix(unsigned int index, vtp& J_triplets) {}

    scalar s;
};

class SpringConstraint : public Constraint
{
public:
    SpringConstraint(scalar stiffness);
    SpringConstraint(scalar stiffness, unsigned int p1, unsigned int p2, scalar length);
    ~SpringConstraint();

    void  EvaluatePotentialEnergy(VectorXs& x, scalar &energy)
    {
        Vector3s xij(x(m_p1*3)-x(m_p2*3),x(m_p1*3+1)-x(m_p2*3+1),x(m_p1*3+2)-x(m_p2*3+2));
        scalar ld=xij.norm()-m_rest_length;
        energy=0.5*s*ld*ld;
    }
    void  EvaluateGradient(VectorXs& x, VectorXs& gradient)
    {
        Vector3s xij(x(m_p1*3)-x(m_p2*3),x(m_p1*3+1)-x(m_p2*3+1),x(m_p1*3+2)-x(m_p2*3+2));
        scalar norm=xij.norm();
        xij=s*(norm-m_rest_length)*xij/norm;
        for(int i=0;i<3;++i)
        {
            gradient(m_p1*3+i)+=xij(i);
            gradient(m_p2*3+i)-=xij(i);
        }
    }
    void  EvaluateHessian(VectorXs& x, vtp& hessian_triplets)
    {
        Vector3s xij(x(m_p1*3)-x(m_p2*3),x(m_p1*3+1)-x(m_p2*3+1),x(m_p1*3+2)-x(m_p2*3+2));
        scalar norm=xij.norm();
        Matrix3s mat3=s*(Matrix3s::Identity()-m_rest_length/norm*(Matrix3s::Identity()-xij*xij.transpose()/(norm*norm)));
        for(int i=0;i<3;++i)
        {
            for(int j=0;j<3;++j)
            {
                scalar value=mat3(i,j);
                hessian_triplets.emplace_back(m_p1*3+i,m_p1*3+j,value);
                hessian_triplets.emplace_back(m_p1*3+i,m_p2*3+j,-value);
                hessian_triplets.emplace_back(m_p2*3+i,m_p1*3+j,-value);
                hessian_triplets.emplace_back(m_p2*3+i,m_p2*3+j,value);
            }
        }
    }
    void  EvaluateWeightedLaplacian(vtp& laplacian_triplets)
    {
        for(int i=0;i<3;++i)
        {
            laplacian_triplets.emplace_back(m_p1*3+i,m_p1*3+i,s);
            laplacian_triplets.emplace_back(m_p1*3+i,m_p2*3+i,-s);
            laplacian_triplets.emplace_back(m_p2*3+i,m_p1*3+i,-s);
            laplacian_triplets.emplace_back(m_p2*3+i,m_p2*3+i,s);
        }
    }
    void EvaluateDVector(unsigned int index, VectorXs& x, VectorXs& d)
    {
        Vector3s xij(x(m_p1*3)-x(m_p2*3),x(m_p1*3+1)-x(m_p2*3+1),x(m_p1*3+2)-x(m_p2*3+2));
        xij.normalized();
        for(int i=0;i<3;++i)
        {
            d(index*3+i)=xij(i)*m_rest_length;
        }
    }
    void EvaluateJMatrix(unsigned int index, vtp& J_triplets)
    {
        for(int i=0;i<3;++i)
        {
            J_triplets.emplace_back(m_p1*3+i,index*3+i,s);
            J_triplets.emplace_back(m_p2*3+i,index*3+i,-s);
        }
    }

protected:
    unsigned int m_p1, m_p2;
    scalar m_rest_length;
};


class balance_solver
{
public:
    balance_solver(){}
    ~balance_solver(){}

    void set_info(MatrixX3i &t_triangles,VectorXs &v_vertices,VectorXi &in_index,
    scalar s_stiffness,scalar h_interval,scalar m_mass,int a_attach)
    {
        t=t_triangles;v=v_vertices;in=in_index;
        v_new=v;y=v;
        s=s_stiffness;h=h_interval;m=m_mass*3.0/v.size();a=a_attach;
        for(int i=0;i<3;++i)
        {
            a_pos(i)=v(a*3+i);
        }
        calc_edgelist();
        predecomposition();
        d.resize(e.rows()*3);
        for(int i=0;i<in.size();++i)
        {
            if(in(i)==a)
            {
                std::cout<<"\nfatal error!\nfixed point is in forcing!\n";
            }
        }
    }
    void set_forces(VectorXs &f_forces)
    {
        f.resize(v.size());
        f.setZero();
        for(int i=0;i<v.size()/3;++i)
        {
            f(i*3+2)=-m*GravityAcceleration;
        }
        for(int i=0;i<in.size();++i)
        {
            for(int j=0;j<3;++j)
            {
                f(in(i)*3+j)+=f_forces(i*3+j);
                //f(in(i)*3+j)+=f_forces(i*6+j)+f_forces(i*6+3+j);
            }
        }
    }
    void set_balance_parameter(int compute_balance_times_,scalar balance_threshold_,scalar newton_rate_,int anderson_length_,
                               int pd_itertimes_,int newton_itertimes_);
    void set_compute_balance_parameter(int compute_balance_times_,scalar balance_threshold_,scalar newton_rate_)
    {
        compute_balance_times=compute_balance_times_;
        balance_threshold=balance_threshold_;
        newton_rate=newton_rate_;
    }
    void calc_edgelist();
    void predecomposition();
    void local_step(VectorXs &pos);
    void global_step();
    void Anderson();
    void pd();
    void newton();
    void newton_raphson();
    void compute_balance();
    void trans_fix();
    bool newton_linesearch(VectorXs &pos, VectorXs &dir, VectorXs &new_pos);
    bool newton_raphson_linesearch(VectorXs &pos, VectorXs &dir, VectorXs &new_pos);
    void balance_state(VectorXs &pos, VectorXs &force_on_vertex);
    scalar compute_balance_result();
    void get_energy(VectorXs &pos, scalar &out_energy, bool with_mass_matrix=true, bool with_fixed_energy=false);
    void compute_jacobi();
    void compute_csr();
    void compute_csr_right();
    
    int method_times[3];

    MatrixX3i t;
    scalar s;
    scalar h;
    scalar m;
    scalar energy;
    int a;
    Vector3s a_pos;
    
    std::vector<Constraint> constraints;
    clock_t runtime[6];
    VectorXs v;
    VectorXs v_new;
    VectorXs v_aa;
    VectorXs y;
    VectorXs f;
    VectorXs d;
    MatrixX2i e;
    VectorXs el;

    scalar compute_balance_times;
    scalar balance_threshold;
    scalar newton_rate;
    scalar balance_result;
    bool x_invariant=false;

    std::list<VectorXs> v_list;
    long unsigned int anderson_length=7;
    int pd_times=10;

    
    //pd
    SimplicialLDLT<spma> solver;
    spma mass_matrix;
    spma J;
    //newton
    SimplicialLDLT<spma> newton_solver;
    bool newton_solver_info=false;
    //newton_raphson
    SparseLU<spma> newton_raphson_solver;
    bool newton_raphson_solver_info=false;
    //jacobi
    VectorXi in;
    MatrixXs jacobi;
    MatrixXs delta;
    SimplicialLDLT<spma> v_solver;
    bool v_solver_info=false;
    VectorXi jacobirow;
    VectorXi jacobicol;
    VectorXs jacobival;
    MatrixXs jacobiright;

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
  .def_readwrite("jacobi",&balance_solver::jacobi)
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
  .def("compute_balance_result",&balance_solver::compute_balance_result)
  .def("compute_jacobi",&balance_solver::compute_jacobi)
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