
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <vector>
using namespace std;
using namespace Eigen;

#if 0
typedef float scalar;
typedef MatrixX3f MatrixX3s;
typedef VectorXf VectorXs;
typedef Vector3f Vector3s;
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
//#define block_vector(a) block<3,1>(3*(a), 0)

class Constraint
{
public:
    Constraint(){}
    Constraint(scalar stiffness):s(stiffness){}
    virtual ~Constraint(){}

    virtual void  EvaluateEnergy(VectorXs& x, scalar &energy) {std::cout<<"?";}
    virtual void  EvaluateGradient(VectorXs& x, VectorXs& gradient) {std::cout<<"?";}
    virtual void  EvaluateHessian(VectorXs& x, vtp& hessian_triplets) {std::cout<<"?";}
    virtual void  EvaluateWeightedLaplacian(vtp& laplacian_triplets) {std::cout<<"?";}
    virtual void EvaluateDVector(unsigned int index, VectorXs& x, VectorXs& d){std::cout<<"?";}
    virtual void EvaluateJMatrix(unsigned int index, vtp& J_triplets) {std::cout<<"?";}
public:
    scalar s;
};

class AttachmentConstraint : public Constraint
{
public:
    AttachmentConstraint(scalar stiffness){}
    AttachmentConstraint(scalar stiffness, unsigned int p0, Vector3s pos):
                        m_p0(p0),m_fixed_point(pos){s=stiffness;}
    ~AttachmentConstraint(){}

    void  EvaluateEnergy(VectorXs& x, scalar &energy)
    {
        Vector3s xij(x(m_p0*3)-m_fixed_point(0),x(m_p0*3+1)-m_fixed_point(1),x(m_p0*3+2)-m_fixed_point(2));
        energy+=0.5*s*xij.squaredNorm();
    }
    void  EvaluateGradient(VectorXs& x, VectorXs& gradient)
    {
        Vector3s xij(x(m_p0*3)-m_fixed_point(0),x(m_p0*3+1)-m_fixed_point(1),x(m_p0*3+2)-m_fixed_point(2));
        for(int i=0;i<3;++i)
        {
            gradient(m_p0*3+i)+=s*xij(i);
        }
        //gradient.block_vector(m_p0)+=s*xij;
    }
    void  EvaluateHessian(VectorXs& x, vtp& hessian_triplets)
    {
        EvaluateWeightedLaplacian(hessian_triplets);
    }
    void  EvaluateWeightedLaplacian(vtp& laplacian_triplets)
    {
        for(int i=0;i<3;++i)
        {
            laplacian_triplets.emplace_back(m_p0*3+i,m_p0*3+i,s);
        }
    }
    void EvaluateDVector(unsigned int index, VectorXs& x, VectorXs& d)
    {
        //d.block_vector(index)=m_fixed_point;
        for(int i=0;i<3;++i)
        {
            d(index*3+i)=m_fixed_point(i);
        }
    }
    void EvaluateJMatrix(unsigned int index, vtp& J_triplets)
    {
        for(int i=0;i<3;++i)
        {
            J_triplets.emplace_back(m_p0*3+i,index*3+i,s);
        }
    }

public:
    unsigned int m_p0;
    Vector3s m_fixed_point;
};

class SpringConstraint : public Constraint
{
public:
    SpringConstraint(scalar stiffness){}
    SpringConstraint(scalar stiffness, unsigned int p1, unsigned int p2, scalar length):
                     m_p1(p1),m_p2(p2),m_rest_length(length){s=stiffness;}
    ~SpringConstraint(){}
public:
    void EvaluateEnergy(VectorXs& x, scalar &energy)
    {
        Vector3s xij(x(m_p1*3)-x(m_p2*3),x(m_p1*3+1)-x(m_p2*3+1),x(m_p1*3+2)-x(m_p2*3+2));
        scalar ld=xij.norm()-m_rest_length;
        energy+=0.5*s*ld*ld;
    }
    void EvaluateGradient(VectorXs& x, VectorXs& gradient)
    {
        Vector3s xij(x(m_p1*3)-x(m_p2*3),x(m_p1*3+1)-x(m_p2*3+1),x(m_p1*3+2)-x(m_p2*3+2));
        scalar norm=xij.norm();
        xij=s*(norm-m_rest_length)*xij/norm;
        //gradient.block_vector(m_p1)+=xij;
        //gradient.block_vector(m_p2)-=xij;
        for(int i=0;i<3;++i)
        {
            gradient(m_p1*3+i)+=xij(i);
            gradient(m_p2*3+i)-=xij(i);
        }
    }
    void EvaluateHessian(VectorXs& x, vtp& hessian_triplets)
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
    void EvaluateWeightedLaplacian(vtp& laplacian_triplets)
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
        xij.normalize();
        //d.block_vector(index)=xij*m_rest_length;
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

public:
    unsigned int m_p1, m_p2;
    scalar m_rest_length;
};