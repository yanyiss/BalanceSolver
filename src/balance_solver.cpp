#include "balance_solver.h"
#include <set>
#include <iomanip>
#include <fstream>
#include <ctime>

void balance_solver::calc_edgelist()
{
    //get stretch spring
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
    for(auto &edge:edgelist)
    {
        int i1=edge.first;
        int i2=edge.second;
        SpringConstraint* sc=new SpringConstraint(s,i1,i2,Vector3s(v(i1*3)-v(i2*3),v(i1*3+1)-v(i2*3+1),v(i1*3+2)-v(i2*3+2)).norm());
        constraints.push_back(sc);
    }

    /* std::vector<std::vector<std::pair<int,int>>> triangle_list;
    triangle_list.resize(v.size()/3);
    for(int i=0;i<t.rows();++i)
    {
        for(int j=0;j<3;++j)
        {
            triangle_list[t(i,j)].emplace_back(t(i,(j+1)%3),t(i,(j+2)%3));
        }
    }
    for(int i=0;i<triangle_list.size();++i)
    {
        for(int j=0;j<triangle_list[i].size();++j)
        {
            int i1=triangle_list[i][j].first;
            int i2=triangle_list[i][j].second;
            if(i1<i2)
            {
                constraints.push_back(std::move(SpringConstraint(s,i1,i2,Vector3s(v(i1*3)-v(i2*3),v(i1*3+1)-v(i2*3+1),v(i1*3+2)-v(i2*3+2)).norm())));
            }
            if(i<i1)
            {
                constraints.push_back(std::move(SpringConstraint(s,i,i1,Vector3s(v(i*3)-v(i1*3),v(i*3+1)-v(i1*3+1),v(i*3+2)-v(i1*3+2)).norm())));
                int id=0;
                for(;id<triangle_list[i1].size();++id)
                {
                    if(triangle_list[i1][id].first==i)
                        break;
                }
                if(id==triangle_list[i1].size())
                    continue;
                int id0=triangle_list[i][j].second;
                int id1=triangle_list[i1][id].second;
                constraints.push_back(std::move(SpringConstraint(s*0.1,id0,id1,Vector3s(v(id0*3)-v(id1*3),v(id0*3+1)-v(id1*3+1),v(id0*3+2)-v(id1*3+2)).norm())));
            }
        }
    } */
    
    AttachmentConstraint* ac=new AttachmentConstraint(s,a,Vector3s(v(a*3),v(a*3+1),v(a*3+2)));
    constraints.push_back(ac);
}


void balance_solver::predecomposition()
{
    int n=v.size();
    vtp tri;
    tri.reserve(constraints.size()*12);
    for(int i=0;i<constraints.size();++i)
    {
        constraints[i]->EvaluateWeightedLaplacian(tri);
    }
    std::cout<<constraints.size();
    /* for(int i=0;i<3;++i)
    {
        tri.emplace_back(a*3+i,a*3+i,s);
    }
 */
    spma Mt2L;
    Mt2L.resize(n,n);
    Mt2L.setFromTriplets(tri.begin(),tri.end());
    pd_solver.compute(Mt2L);
    if(pd_solver.info()!=0)
    {
        std::cout<<"pd_solver factorize failed"<<std::endl;
        exit(0);
    }
    else
    {
        std::cout<<"pd_solver factorize succeed"<<std::endl;
    }
    
    tri.clear();
    for(int i=0;i<constraints.size();++i)
    {
        constraints[i]->EvaluateJMatrix(i,tri);
    }
    J.resize(n,constraints.size()*3);
    J.setFromTriplets(tri.begin(),tri.end());
    d.resize(constraints.size()*3);
}

void balance_solver::local_step(VectorXs &pos)
{
    for(int i=0;i<constraints.size();++i)
    {
        constraints[i]->EvaluateDVector(i,pos,d);
    }
}

void balance_solver::global_step()
{
    v_list.push_back(std::move(v_new));
    v_new=pd_solver.solve(J*d+f);
}

void balance_solver::Anderson()
{
    if(v_list.size()>anderson_length)
        v_list.pop_front();
    if(v_list.size()<anderson_length)
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

void balance_solver::pd()
{
    v_new=v;
    for(long unsigned int i=0;i<anderson_length;++i)
    {
        local_step(v_new);
        global_step();
    }
    scalar energy;
    get_energy(v_new,energy);
    for(int i=0;i<pd_times;++i)
    {
        Anderson();
        scalar aa_energy;
        get_energy(v_aa,aa_energy);
        if(aa_energy>energy)
        {
            local_step(v_new);
        }
        else
        {
            local_step(v_aa);
            v_new=std::move(v_aa);
        }
        global_step();
        get_energy(v_new,energy);


        /* Anderson();
        local_step(v_aa);
        scalar aa_energy;
        get_energy(v_aa,aa_energy);
        if(aa_energy>energy)
        {
            local_step(v_new);
        }
        else
        {
            v_new=std::move(v_aa);
            energy=aa_energy;
        }
        //local_step(v_new);
        global_step(); */
    }
    std::cout<<"pd once"<<std::endl;
    v=v_new;
    v_list.clear();


    VectorXs force_on_vertex;
    balance_state(v_new,force_on_vertex);
    balance_result=force_on_vertex.sum()/(force_on_vertex.size()*m);
    std::cout<<"Mean Norm Rate: "<<balance_result<<std::endl;
}

bool balance_solver::newton_linesearch(VectorXs &pos, VectorXs &dir)
{
    scalar current_energy;
    get_energy(pos,current_energy);

    scalar times=2.0;
    VectorXs v_tdv;
    scalar step_energy;
    do
    {
        times*=0.5;
        v_tdv=pos+times*dir;
        get_energy(v_tdv,step_energy);
    }while(times>1.0e-6&&step_energy>current_energy);
    if(times>1.0e-6)
    {
        v=std::move(v_tdv);
        return true;
    }
    else
    {
        return false;
    }
}

bool balance_solver::newton()
{
    VectorXs gradient=-f;
    for(int i=0;i<constraints.size();++i)
    {
        constraints[i]->EvaluateGradient(v,gradient);
    }

    vtp tri;
    tri.reserve(constraints.size()*36);
    for(int i=0;i<constraints.size();++i)
    {
        constraints[i]->EvaluateHessian(v,tri);
    }
    spma hessian;
    hessian.resize(v.size(),v.size());
    hessian.setFromTriplets(tri.begin(),tri.end());

    if(!newton_solver_info)
    {
        newton_solver.analyzePattern(hessian);
        newton_solver_info=true;
    }
    newton_solver.factorize(hessian);
    //hessian_solver.compute(hessian);
    if(newton_solver.info()!=0)
    {
        std::cout<<"newton hessian factorize failed"<<std::endl;
        return false;
    }
    VectorXs dir=-newton_solver.solve(gradient);
    return newton_linesearch(v,dir);
}

void balance_solver::compute_balance()
{
    VectorXs force_on_vertex;
    int itertimes=0;
    do
    {
        balance_state(v,force_on_vertex);
        balance_result=force_on_vertex.sum()/(force_on_vertex.size()*m);
        ++itertimes;
        if(itertimes>compute_balance_times||balance_result<balance_threshold)
        {
            break;
        }
        if(balance_result>newton_rate)
        {
            pd();
        }
        else
        {
            if(!newton())
            {
                pd();
            }
        }
    }while(1);
    std::cout<<"\nitertimes: "<<itertimes;
    std::cout<<"\nbalance info:\n";
    std::cout<<"Mean Norm: "<<balance_result<<std::endl;
    std::cout<<"Max Norm: "<<force_on_vertex.maxCoeff()/m<<std::endl;
}

void balance_solver::balance_state(VectorXs &pos, VectorXs &force_on_vertex)
{
    VectorXs force=-f;
    for(int i=0;i<constraints.size();++i)
    {
        constraints[i]->EvaluateGradient(pos,force);
    }
    force_on_vertex.resize(f.size()/3);
    for(int i=0;i<f.size()/3;++i)
    {
        scalar sum=0;
        for(int j=0;j<3;++j)
        {
            sum+=force(i*3+j)*force(i*3+j);
        }
        force_on_vertex(i)=sqrt(sum);
    }
}

void balance_solver::get_energy(VectorXs &pos,scalar &out_energy)
{
    out_energy=-pos.dot(f);
    for(int i=0;i<constraints.size();++i)
    {
        constraints[i]->EvaluateEnergy(pos,out_energy);
    }
}

void balance_solver::compute_jacobi()
{
    vtp tri;
    tri.reserve(constraints.size()*36);
    for(int i=0;i<constraints.size();++i)
    {
        constraints[i]->EvaluateHessian(v,tri);
    }
    for(auto &tr:tri)
    {
        if((tr.row()==a)!=(tr.col()==a))
        {
            std::cout<<"fesf";
            exit(0);
        }
    }
    spma v_mat;
    v_mat.resize(v.size(),v.size());
    v_mat.setFromTriplets(tri.begin(),tri.end());
    
    if(!v_solver_info)
    {
        delta.resize(v.size(),in.size()*3-3);
        delta.setZero();
        scalar s_inv=1.0/s;
#if 0
        for(int i=0;i<in.size()-1;++i)
        {
            for(int j=0;j<3;++j)
            {
                delta(in(i)*3+j,i*3+j)=s_inv;
                delta(in(in.size()-1)*3+j,i*3+j)=-s_inv;
            }
        }
#else
        for(int i=0;i<in.size()-1;++i)
        {
            for(int j=0;j<3;++j)
            {
                delta(in(i)*3+j,i*3+j)=s_inv;
                delta(in(i+1)*3+j,i*3+j)=-s_inv;
            }
        }
#endif
        v_solver.analyzePattern(v_mat);
        v_solver_info=true;
    }
    v_solver.factorize(v_mat);
    
    jacobi.resize(v.size(),in.size()*3);
    jacobi=v_solver.solve(delta);
}
