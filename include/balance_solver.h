#include <pybind11/pybind11.h>
#include <iostream>
#include <Eigen/Dense>
#include <pybind11/eigen.h>
#include <Eigen/SparseCholesky>
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
typedef MatrixXf MatrixXs;
//typedef Matrix3f Matrix3s;
#else
typedef double scalar;
typedef MatrixX3d MatrixX3s;
typedef VectorXd VectorXs;
typedef Vector3d Vector3s;
typedef MatrixXd MatrixXs;
//typedef Matirx3d Matrix3s;
#endif


typedef SparseMatrix<scalar> spma;
typedef std::vector<Triplet<scalar>> vtp;


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
        calc_edgelist();
        init_adj();
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
    void init_adj();
    void predecomposition();
    void local_step(VectorXs &pos);
    void global_step();
    void Anderson();
    void pd();
    void newton();
    void compute_balance();
    void trans_fix();
    bool linesearch(VectorXs &pos, VectorXs &dir, VectorXs &new_pos);
    void balance_state(VectorXs &force_on_vertex);
    void get_energy(VectorXs &pos, scalar &out_energy);
    void compute_jacobi();
    


    MatrixX3i t;
    scalar s;
    scalar h;
    scalar m;
    scalar energy;
    int a;
    

    clock_t runtime[6];
    VectorXs v;
    VectorXs v_new;
    VectorXs v_aa;
    VectorXs y;
    VectorXs f;
    VectorXs d;
    MatrixX2i e;
    VectorXs el;
    std::vector<std::vector<int>> adj;
    std::vector<std::vector<scalar>> adj_len;

    scalar compute_balance_times;
    scalar balance_threshold;
    scalar newton_rate;
    scalar balance_result;

    std::list<VectorXs> v_list;
    long unsigned int anderson_length=7;
    int pd_times=10;

    
    //pd
    SimplicialLDLT<spma> solver;
    spma Mt2L;
    spma J;
    //newton
    SimplicialLDLT<spma,Eigen::Upper> hessian_solver;
    bool hessian_solver_info=false;
    //jacobi
    VectorXi in;
    MatrixXs jacobi;
    MatrixXs delta;
    SimplicialLDLT<spma> v_solver;
    bool v_solver_info=false;
};


PYBIND11_MODULE(balance_solver, m) {
  
  py::class_<balance_solver>(m,"balance_solver")
  .def(py::init<>())
  .def_readwrite("v",&balance_solver::v)
  .def_readwrite("jacobi",&balance_solver::jacobi)
  .def_readwrite("balance_result",&balance_solver::balance_result)
  .def("set_info",&balance_solver::set_info)
  .def("set_forces",&balance_solver::set_forces)
  .def("set_compute_balance_parameter",&balance_solver::set_compute_balance_parameter)
  .def("compute_balance",&balance_solver::compute_balance)
  .def("compute_jacobi",&balance_solver::compute_jacobi);
}
//cmake .. -DCMAKE_BUILD_TYPE=Release
//make;rm /home/yanyisheshou/Program/TarpDesign/algorithm/balance_solver.cpython-39-x86_64-linux-gnu.so;mv balance_solver.cpython-39-x86_64-linux-gnu.so /home/yanyisheshou/Program/TarpDesign/algorithm/