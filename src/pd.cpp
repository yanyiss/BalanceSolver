#include "pd.h"
#include <set>
void DiffSimulation::calc_edgelist()
{
    std::set<std::pair<int,int>> edgelist;
    for(int i=0;i<t.rows();++i)
    {
        for(int j=0;j<3;++j)
        {
            int id0=t(i,j);
            int id1=t(i,(j+1)%3);
            if(id0<id1)
            {
                edgelist.insert(std::pair<int,int>(id0,id1));
            }
            else
            {
                edgelist.insert(std::pair<int,int>(id1,id0));
            }
        }
    }
    e.resize(edgelist.size(),3);
    int i=0;
    for(auto edge:edgelist)
    {
        e(i,0)=edge.first;
        e(i,1)=edge.second;
        e(i,2)=0;
        for(int j=0;j<3;++j)
            e(i,2)+=(v(edge.first*3+j)-v(edge.second*3+j))
            *(v(edge.first*3+j)-v(edge.second*3+j));
        e(i,2)=std::sqrt(e(i,2));
        //e(i,2)=(v.row(edge.first)-v.row(edge.second)).norm();
        //std::cout<<i<<" "<<e(i,0)<<" "<<e(i,1)<<std::endl;
        ++i;
    }
    //std::cout<<"edge num: "<<e.rows()<<std::endl;
}

void DiffSimulation::predecomposition()
{
    int n=v.size();
    vtp tri;
    tri.reserve(n+e.rows()*12);
    //inertia
    for(int i=0;i<n;++i)
    {
        tri.emplace_back(i,i,m);
    }
    scalar h2_s=h*h*s;
    //spring force
    for(int i=0;i<e.rows();++i)
    {
        int id0=e(i,0)*3;
        int id1=e(i,1)*3;
        for(int j=0;j<3;++j)
        {
            tri.emplace_back(id0+j,id0+j,h2_s);
            tri.emplace_back(id1+j,id1+j,h2_s);
            tri.emplace_back(id0+j,id1+j,-h2_s);
            tri.emplace_back(id1+j,id0+j,-h2_s);
        }
    }
    //fixed point
    /* for(int i=0;i<3;++i)
    {
        tri.emplace_back(a*3+i,a*3+i,h2_s);
    } */

    Mt2L.resize(n,n);
    Mt2L.setFromTriplets(tri.begin(),tri.end());
    solver.compute(Mt2L);
    

    tri.clear();
    //spring force
    for(int i=0;i<e.rows();++i)
    {
        int id0=e(i,0)*3;
        int id1=e(i,1)*3;
        for(int j=0;j<3;++j)
        {
            tri.emplace_back(id0+j,i*3+j,s);
            tri.emplace_back(id1+j,i*3+j,-s);
        }
    }
    //fixed point
    /* for(int i=0;i<3;++i)
    {
        tri.emplace_back(a*3+i,e.rows()*3+i,s);
    } */
    J.resize(n,e.rows()*3);
    J.setFromTriplets(tri.begin(),tri.end());
}

void DiffSimulation::local_step(VectorXs &pos)
{
    for(int i=0;i<e.rows();++i)
    {
        int id0=e(i,0)*3;
        int id1=e(i,1)*3;
        Vector3s vij;
        vij<<pos(id0)-pos(id1),pos(id0+1)-pos(id1+1),
        pos(id0+2)-pos(id1+2);
        vij.normalize();
        for(int j=0;j<3;++j)
        {
            d(i*3+j)=vij(j)*e(i,2);
        }
        /* std::cout<<id0/3<<id1/3<<std::endl;
        std::cout<<vij<<std::endl;
        std::cout<<std::endl; */
    }
}

void DiffSimulation::global_step()
{
    v_list.push_back(std::move(v_new));
    v_new=solver.solve(m*y+h*h*(J*d+f));
}

void DiffSimulation::Anderson()
{
    if(v_list.size()>max_num)
        v_list.pop_front();
    if(v_list.size()<max_num)
        return;

    int n=v_list.size();
    MatrixXs A;
    A.resize(v.size(),n-1);
    VectorXs b=v_new-v_list.back();

    int i=0;
    auto v_begin=v_list.begin();
    ++v_begin;
    for(;v_begin!=v_list.end();++v_begin)
    {
        auto itr=v_begin;--itr;
        A.col(i)=*v_begin-*itr-b;
        ++i;
    }
    MatrixXs B;
    B.resize(n-1,n-1);
    for(int i=0;i<n-1;++i)
    {
        for(int j=0;j<n-1;++j)
        {
            B(i,j)=A.col(i).transpose()*A.col(j);
        }
    }
    VectorXs alpha=B.ldlt().solve(-A.transpose()*b);
    

    v_aa=(1-alpha.sum())*v_new;
    i=0;
    v_begin=v_list.begin();
    ++v_begin;
    for(;v_begin!=v_list.end();++v_begin)
    {
        v_aa+=alpha(i)*(*v_begin);
        ++i;
    }
}

void DiffSimulation::Opt()
{
    for(int i=0;i<max_num;++i)
    {
        local_step(v_new);
        global_step();
        //std::cout<<i<<" "<<(v_new-v_list.back()).norm()<<std::endl;
        //get_energy(v_new,energy);
        //std::cout<<"energy: "<<energy<<std::endl;
    }
    get_energy(v_new,energy);
    for(int i=0;i<10;++i)
    {
        Anderson();
        local_step(v_aa);
        scalar aa_energy;
        get_energy(v_aa,aa_energy);
        if(aa_energy>energy)
        {
            local_step(v_new);
            //get_energy(v_new,energy);
        }
        else
        {
            v_new=std::move(v_aa);
            energy=aa_energy;
        }
        
        global_step();
        //Anderson();
        //std::cout<<i<<" "<<(v_new-v_list.back()).norm()<<std::endl;
        //std::cout<<"energy: "<<energy<<std::endl;
    }

    v_list.clear();
    //translate the mesh
    Vector3s trans(v(a*3)-v_new(a*3),v(a*3+1)-v_new(a*3+1),v(a*3+2)-v_new(a*3+2));
    for(int i=0;i<v.size()/3;++i)
    {
        for(int j=0;j<3;++j)
        {
            v_new(i*3+j)+=trans(j);
        }
    }
    //bool balance=print_balance_info();
    
    y=v_new+(v_new-v)*0.9;
    v=v_new;
}

bool DiffSimulation::print_balance_info()
{
    MatrixX3s balance;
    balance.resize(f.size()/3,3);
    for(int i=0;i<f.size()/3;++i)
    {
        for(int j=0;j<3;++j)
        {
            balance(i,j)=f(i*3+j);
        }
    }
    Vector3s vij;
    for(int i=0;i<e.rows();++i)
    {
        int id0=e(i,0)*3;
        int id1=e(i,1)*3;
        vij<<v_new(id0)-v_new(id1),v_new(id0+1)-v_new(id1+1),
        v_new(id0+2)-v_new(id1+2);
        vij*=s*(1-e(i,2)/vij.norm());
        balance.row(e(i,0))-=vij;
        balance.row(e(i,1))+=vij;
    }
    scalar meanNorm=balance.rowwise().norm().sum()/balance.rows();
    std::cout<<"\nBalance Info:\n";
    std::cout<<"Mean Norm: "<<meanNorm<<std::endl;
    if(meanNorm<9.8*m)
        return true;
    else
        return false;
    /* std::cout<<"balance info:\n";
    std::cout<<"resultant force: "<<f.sum()<<std::endl;
    std::cout<<"Frobenius Norm: "<<balance.norm()<<std::endl;
    VectorXs norm=balance.rowwise().norm();
    std::cout<<"Mean Norm: "<<norm.sum()/norm.size()<<std::endl;
    std::cout<<"Max Norm: "<<norm.maxCoeff(); */
    /* for(int i=0;i<norm.size();++i)
    {
        std::cout<<i<<" "<<norm(i)<<std::endl;
    } */
}

void DiffSimulation::get_energy(VectorXs &pos,scalar &energy)
{
    energy= 0.5*pos.transpose()*(Mt2L*pos);
    scalar e_temp=pos.transpose()*(J*d+f);
    energy-=e_temp*h*h;
    e_temp=-pos.transpose()*y;
    energy-=e_temp*m;
    energy+=0.5*s*h*h*d.squaredNorm();
}

void DiffSimulation::compute_jacobi()
{
    vtp tri;
    tri.reserve(e.rows()*36+3);
    MatrixX3s vjk;
    vjk.resize(3,3);
    for(int i=0;i<e.rows();++i)
    {
        int id0=e(i,0);
        int id1=e(i,1);
        Vector3s v01(v(id1*3)-v(id0*3),v(id1*3+1)-v(id0*3+1),v(id1*3+2)-v(id0*3+2));
        scalar len_inv=1.0/v01.norm();
        scalar len_inv3=len_inv*len_inv*len_inv;
        vjk<<1-len_inv+v01(0)*v01(0)*len_inv3,v01(0)*v01(1)*len_inv3,v01(0)*v01(2)*len_inv3,
             v01(1)*v01(0)*len_inv3,1-len_inv+v01(1)*v01(1)*len_inv3,v01(1)*v01(2)*len_inv3,
             v01(2)*v01(0)*len_inv3,v01(2)*v01(1)*len_inv3,1-len_inv+v01(2)*v01(2)*len_inv3;
        for(int j=0;j<3;++j)
        {
            for(int k=0;k<3;++k)
            {
                tri.emplace_back(id0*3+j,id0*3+k,vjk(j,k));
                tri.emplace_back(id0*3+j,id1*3+k,-vjk(j,k));
                tri.emplace_back(id1*3+j,id0*3+k,-vjk(j,k));
                tri.emplace_back(id1*3+j,id1*3+k,vjk(j,k));
            }
        }
    }
    for(int i=0;i<3;++i)
    {
        tri.emplace_back(a*3+i,a*3+i,1.0);
    }
    spma v_mat;
    v_mat.resize(v.size(),v.size());
    v_mat.setFromTriplets(tri.begin(),tri.end());
    if(!v_solver_info)
    {
        delta.resize(v.size(),in.size()*3);
        delta.setZero();
        scalar s_inv=1.0/s;
        for(int i=0;i<in.size();++i)
        {
            for(int j=0;j<3;++j)
            {
                delta(in(i)*3+j,i*3+j)=s_inv;
            }
        }
        v_solver.analyzePattern(v_mat);
        v_solver_info=true;
    }
    v_solver.factorize(v_mat);
    jacobi=v_solver.solve(delta);
}
