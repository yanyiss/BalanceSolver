#include <pybind11/pybind11.h>
#include <iostream>
#include <Eigen/Dense>
#include <pybind11/eigen.h>
#include <Eigen/SparseCholesky>
#include <list>
#include <omp.h>
using namespace std;
namespace py=pybind11;
using namespace Eigen;
#define GravityAcceleration 9.8
//reference https://zhuanlan.zhihu.com/p/441301638
//          Anderson Acceleration for Geometry Optimization and Physics Simulation


class DiffSimulation
{
public:
    DiffSimulation(){}
    ~DiffSimulation(){}
public:
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

    void set_info(MatrixX3i &t_triangles,VectorXs &v_vertices,VectorXi &in_index,
    scalar s_stiffness,scalar h_interval,scalar m_mass,int a_attach)
    {
        t=t_triangles;v=v_vertices;in=in_index;
        v_new=v;y=v;
        s=s_stiffness;h=h_interval;m=m_mass*3.0/v.size();a=a_attach;
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
        /* std::cout<<"triangles:\n";
        std::cout<<0<<"    "<<t(0,0)<<" "<<t(0,1)<<" "<<t(0,2)<<std::endl;
        std::cout<<1<<"    "<<t(1,0)<<" "<<t(1,1)<<" "<<t(1,2)<<std::endl;
        std::cout<<"positions:\n";
        std::cout<<0<<"    "<<v(0)<<" "<<v(1)<<" "<<v(2)<<std::endl;
        std::cout<<1<<"    "<<v(3)<<" "<<v(4)<<" "<<v(5)<<std::endl;
        std::cout<<"stiffness: "<<s<<std::endl;
        std::cout<<"interval: "<<h<<std::endl;
        std::cout<<"mass: "<<m<<std::endl; */
    }
    void set_vertices(VectorXs &v_vertices) { v=v_vertices; }
    void set_forces(VectorXs &f_forces)
    {
        f.resize(v.size());
        f.setZero();
        for(int i=0;i<v.size()/3;++i)
        {
            f(i*3+2)=-m*GravityAcceleration;
        }
#if 0
        Vector3s sum_force(0.0,0.0,0.0);
        for(int i=0;i<in.size()-1;++i)
        {
            for(int j=0;j<3;++j)
            {
                f(in(i)*3+j)+=f_forces(i*3+j);
                sum_force(j)+=f_forces(i*3+j);
            }
        }
        int f_last=in(in.size()-1)*3;
        f(f_last)+=-sum_force(0);
        f(f_last+1)+=-sum_force(1);
        f(f_last+2)+=m*v.size()/3.0*GravityAcceleration-sum_force(2);
#else
        for(int i=0;i<in.size();++i)
        {
            for(int j=0;j<3;++j)
            {
                f(in(i)*3+j)+=f_forces(i*3+j);
            }
        }
#endif

        /* float vm=4.5f*9.8f/(v.size()/3);
        f.resize(v.size());f.setZero();
        for(int i=0;i<v.size()/3;++i)
        {
            f(i*3+2)=-vm;
        }
        std::cout<<"fesfsfs"<<f.sum()<<std::endl;
        f(0)=-70;f(1)=-70;f(2)+=-20;
        f(3)=70;f(4)=-70;f(5)+=-20;
        f(12)=-70;f(13)=70;f(14)+=-20;
        f(15)=70;f(16)=70;f(17)+=-20;
        f(168)=0;f(169)=100;f(170)+=-20;
        f(468)=0;f(469)=-100;f(470)+=-20;
        f(6)=100;f(7)=0;f(8)+=82.05;
        f(9)=-100;f(10)=0;f(11)+=82.05;
        std::cout<<f.sum()<<std::endl; */
    }
    void calc_edgelist();
    void predecomposition();
    void local_step(VectorXs &pos);
    void global_step();
    void linesearch(VectorXs &pos, VectorXs &dir, VectorXs &new_pos);
    void newton();
    void Anderson();
    void Opt();
    bool print_balance_info(scalar cof);
    void get_energy(VectorXs &pos, scalar &energy);
    void compute_jacobi();
    


    MatrixX3i t;
    scalar s;
    scalar h;
    scalar m;
    scalar energy;
    int a;

    VectorXs v;
    VectorXs v_new;
    VectorXs v_aa;
    VectorXs y;
    VectorXs f;
    VectorXs d;
    MatrixX3s e;

    std::list<VectorXs> v_list;
    int max_num=7;

    typedef SparseMatrix<scalar> spma;
    typedef std::vector<Triplet<scalar>> vtp;
    
    SimplicialLDLT<spma> solver;
    spma Mt2L;
    spma J;
    
    VectorXi in;
    MatrixXs jacobi;
    MatrixXs delta;
    SimplicialLDLT<spma> v_solver;
    MatrixXs leftmat;
    bool v_solver_info=false;

    SimplicialLDLT<spma> hessian_solver;
    bool hessian_solver_info=false;
};


PYBIND11_MODULE(pd_cpp, m) {
  
  py::class_<DiffSimulation>(m,"DiffSimulation")
  .def(py::init<>())
  .def_readwrite("v",&DiffSimulation::v)
  .def_readwrite("jacobi",&DiffSimulation::jacobi)
  .def_readwrite("leftmat",&DiffSimulation::leftmat)
  .def("set_info",&DiffSimulation::set_info)
  .def("set_vertices",&DiffSimulation::set_vertices)
  .def("set_forces",&DiffSimulation::set_forces)
  .def("newton",&DiffSimulation::newton)
  .def("Opt",&DiffSimulation::Opt)
  .def("print_balance_info",&DiffSimulation::print_balance_info)
  .def("compute_jacobi",&DiffSimulation::compute_jacobi);
}