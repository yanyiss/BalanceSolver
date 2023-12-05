#include "balance_solver.h"
#include <set>
#include <iomanip>
#include <ctime>

void balance_solver::calc_edgelist()
{
#if 0
    std::vector<std::vector<std::pair<int,int>>> triangle_list;
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
    }
#else
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
    e.resize(edgelist.size(),2);
    el.resize(edgelist.size());
    int i=0;
    for(auto edge:edgelist)
    {
        e(i,0)=edge.first;
        e(i,1)=edge.second;
        el(i)=0;
        for(int j=0;j<3;++j)
            el(i)+=(v(edge.first*3+j)-v(edge.second*3+j))
                   *(v(edge.first*3+j)-v(edge.second*3+j));
        el(i)=std::sqrt(el(i));
        ++i;
    }
#endif
}


void balance_solver::predecomposition()
{
#if 0
    int n=v.size();
    vtp tri;
    tri.reserve(n);
    //inertia
    for(int i=0;i<n;++i)
    {
        tri.emplace_back(i,i,m);
    }
    mass_matrix.resize(n,n);
    mass_matrix.setFromTriplets(tri.begin(),tri.end());
    //spring force
    tri.clear();
    tri.reserve(constraint.size()*12);
    for(int i=0;i<constraints.size();++i)
    {
        constraints.EvaluateWeightedLaplacian(tri);
    }
    Mt2L.resize(n,n);
    Mt2L.setFromTriplets(tri.begin(),tri.end());
    solver.compute(Mt2L*h*h+mass_matrix);
    
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
    J.resize(n,e.rows()*3);
    J.setFromTriplets(tri.begin(),tri.end());
#else
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
    spma Mt2L;
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
    J.resize(n,e.rows()*3);
    J.setFromTriplets(tri.begin(),tri.end());
#endif
}

void balance_solver::local_step(VectorXs &pos)
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
            d(i*3+j)=vij(j)*el(i);
        }
    }
}

void balance_solver::global_step()
{
    v_list.push_back(std::move(v_new));
    v_new=solver.solve(m*y+h*h*(J*d+f));
}

bool balance_solver::newton_linesearch(VectorXs &pos, VectorXs &dir, VectorXs &new_pos)
{
    clock_t start=clock();
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
    runtime[4]=clock()-start;
    if(times>1.0e-6)
    {
        energy=step_energy;
        new_pos=std::move(v_tdv);
        return true;
    }
    else
    {
        energy=current_energy;
        return false;
    }
}
bool balance_solver::newton_raphson_linesearch(VectorXs &pos, VectorXs &dir, VectorXs &new_pos)
{
    clock_t start=clock();
    scalar current_energy;
    get_energy(pos,current_energy,false);

    scalar times=2.0;
    VectorXs v_tdv;
    scalar step_energy;
    do
    {
        times*=0.5;
        v_tdv=pos+times*dir;
        get_energy(v_tdv,step_energy,false);
    }while(times>1.0e-6&&step_energy>current_energy);
    runtime[4]=clock()-start;
    if(times>1.0e-6)
    {
        energy=step_energy;
        new_pos=std::move(v_tdv);
        return true;
    }
    else
    {
        energy=current_energy;
        return false;
    }
}

void balance_solver::trans_fix()
{
    Vector3s trans(v(a*3)-v_new(a*3),v(a*3+1)-v_new(a*3+1),v(a*3+2)-v_new(a*3+2));
    for(int i=0;i<v.size()/3;++i)
    {
        for(int j=0;j<3;++j)
        {
            v_new(i*3+j)+=trans(j);
        }
    }
    y=v_new+(v_new-v)*0.1;
    v=v_new;
}

void balance_solver::newton()
{
    clock_t start=clock();
    //compute gradient and hessian
    //for(int iter=0;iter<5;++iter)
    //{
    VectorXs gradient;
    gradient.resize(v.size());
    gradient.setZero();
    vtp tri;
    Vector3s v01;
    MatrixX3s vijk;
    vijk.resize(3,3);
    scalar h2s=h*h*s;
    for(int i=0;i<e.rows();++i)
    {
        int id0=e(i,0);
        int id1=e(i,1);
        v01<<v(id0*3)-v(id1*3),v(id0*3+1)-v(id1*3+1),v(id0*3+2)-v(id1*3+2);
        Vector3s g01=s*(v01.norm()-el(i))*v01.normalized();
        for(int j=0;j<3;++j)
        {
            gradient(id0*3+j)+=g01(j);
            gradient(id1*3+j)-=g01(j);
        }

        scalar len_inv=1.0/v01.norm();
        scalar len_inv3=el(i)*len_inv*len_inv*len_inv;
        len_inv*=el(i);
        vijk<<1.0-len_inv+v01(0)*v01(0)*len_inv3,             v01(0)*v01(1)*len_inv3,             v01(0)*v01(2)*len_inv3,
                          v01(1)*v01(0)*len_inv3, 1.0-len_inv+v01(1)*v01(1)*len_inv3,             v01(1)*v01(2)*len_inv3,
                          v01(2)*v01(0)*len_inv3,             v01(2)*v01(1)*len_inv3, 1.0-len_inv+v01(2)*v01(2)*len_inv3;
        for(int j=0;j<3;++j)
        {
            for(int k=0;k<3;++k)
            {
                tri.emplace_back(id0*3+j,id0*3+k,h2s*vijk(j,k));
                tri.emplace_back(id0*3+j,id1*3+k,-h2s*vijk(j,k));
                tri.emplace_back(id1*3+j,id0*3+k,-h2s*vijk(j,k));
                tri.emplace_back(id1*3+j,id1*3+k,h2s*vijk(j,k));
            }
        }
    }

    gradient-=f;
    gradient=m*(v-y)+h*h*gradient;

    for(int i=0;i<v.size();++i)
    {
        tri.emplace_back(i,i,m);
    }
    runtime[0]+=clock()-start;
    start=clock();
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
        for(int i=0;i<5;++i)
        {
            pd();
        }
        return;
    }
    runtime[1]+=clock()-start;
    start=clock();

    VectorXs dir=-newton_solver.solve(gradient);
    runtime[2]+=clock()-start;
    start=clock();
    if(newton_linesearch(v,dir,v_new))
    {
        ++method_times[1];
        trans_fix();
    }
    else
    {
        ++method_times[0];
        start=clock();
        std::cout<<"from pd"<<std::endl;
        for(int i=0;i<5;++i)
        {
            pd();
        }
        runtime[5]+=clock()-start;
    }
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
    get_energy(v_new,energy);
    for(int i=0;i<pd_times;++i)
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
    }

    v_list.clear();
    //translate the mesh
    trans_fix();
}

void balance_solver::newton_raphson()
{
    VectorXs gradient;
    gradient.resize(v.size());
    gradient.setZero();
    vtp tri;
    Vector3s v01;
    MatrixX3s vijk;
    vijk.resize(3,3);
    for(int i=0;i<e.rows();++i)
    {
        int id0=e(i,0);
        int id1=e(i,1);
        v01<<v(id0*3)-v(id1*3),v(id0*3+1)-v(id1*3+1),v(id0*3+2)-v(id1*3+2);
        Vector3s g01=s*(v01.norm()-el(i))*v01.normalized();
        for(int j=0;j<3;++j)
        {
            gradient(id0*3+j)+=g01(j);
            gradient(id1*3+j)-=g01(j);
        }

        scalar len_inv=1.0/v01.norm();
        scalar len_inv3=el(i)*len_inv*len_inv*len_inv;
        len_inv*=el(i);
        vijk<<1.0-len_inv+v01(0)*v01(0)*len_inv3,             v01(0)*v01(1)*len_inv3,             v01(0)*v01(2)*len_inv3,
                          v01(1)*v01(0)*len_inv3, 1.0-len_inv+v01(1)*v01(1)*len_inv3,             v01(1)*v01(2)*len_inv3,
                          v01(2)*v01(0)*len_inv3,             v01(2)*v01(1)*len_inv3, 1.0-len_inv+v01(2)*v01(2)*len_inv3;
        for(int j=0;j<3;++j)
        {
            for(int k=0;k<3;++k)
            {
                tri.emplace_back(id0*3+j,id0*3+k,s*vijk(j,k));
                tri.emplace_back(id0*3+j,id1*3+k,-s*vijk(j,k));
                tri.emplace_back(id1*3+j,id0*3+k,-s*vijk(j,k));
                tri.emplace_back(id1*3+j,id1*3+k,s*vijk(j,k));
            }
        }
    }

    gradient-=f;

    spma hessian;
    hessian.resize(v.size(),v.size());
    hessian.setFromTriplets(tri.begin(),tri.end());
    if(!newton_raphson_solver_info)
    {
        newton_raphson_solver.analyzePattern(hessian);
        newton_raphson_solver_info=true;
    }
    newton_raphson_solver.factorize(hessian);
    //hessian_solver.compute(hessian);
    if(newton_raphson_solver.info()!=0)
    {
        std::cout<<"newton raphson hessian factorize failed"<<std::endl;
        for(int i=0;i<5;++i)
        {
            pd();
        }
        return;
        //exit(0);
    }

    VectorXs dir=-newton_raphson_solver.solve(gradient);


    if(newton_raphson_linesearch(v,dir,v_new))
    {
        ++method_times[2];
        trans_fix();
    }
    else
    {
        ++method_times[1];
        newton();
        /* std::cout<<"to pd"<<std::endl;
        for(int i=0;i<5;++i)
        {
            pd();
        } */
    }
}

void balance_solver::compute_balance()
{
    VectorXs force_on_vertex;
    int itertimes=0;
    for(int i=0;i<6;++i)
    {
        runtime[i]=0;
    }
    for(int i=0;i<3;++i)
    {
        method_times[i]=0;
    }
    do
    {
        clock_t start=clock();
        balance_state(v_new,force_on_vertex);
        runtime[3]+=clock()-start;
        balance_result=force_on_vertex.sum()/(force_on_vertex.size()*m);
        ++itertimes;
        if(itertimes>compute_balance_times||balance_result<balance_threshold)
        { 
            std::cout<<"\nitertimes: "<<itertimes;
            std::cout<<"\nbalance info:\n";
            std::cout<<"energy: "<<energy<<std::endl;
            std::cout<<"Mean Norm Rate: "<<balance_result<<std::endl;
            std::cout<<"Max Norm Rate: "<<force_on_vertex.maxCoeff()/m<<std::endl;

            clock_t sum=0;
            for(int i=0;i<6;++i)
            {
                std::cout<<runtime[i]<<" ";
                sum+=runtime[i];
            }
            std::cout<<runtime[1]*1.0/sum<<std::endl;
            std::cout<<"sum: "<<sum<<std::endl;
            std::cout<<std::endl;
            std::cout<<method_times[0]<<" "<<method_times[1]<<" "<<method_times[2]<<std::endl;
            return;
        }
        if(balance_result>newton_rate)
        {
            newton();
        }
        else
        {
            newton_raphson();
        }
    }while(1);
    std::cout<<"\nitertimes: "<<itertimes;
    std::cout<<"\nbalance info:\n";
    std::cout<<"energy: "<<energy<<std::endl;
    std::cout<<"Mean Norm: "<<balance_result<<std::endl;
    std::cout<<"Max Norm: "<<force_on_vertex.maxCoeff()/m<<std::endl;
    clock_t sum=0;
            for(int i=0;i<6;++i)
            {
                std::cout<<runtime[i]<<" ";
                sum+=runtime[i];
            }
    std::cout<<runtime[1]*1.0/sum<<std::endl;
    std::cout<<"sum: "<<sum<<std::endl;
    std::cout<<std::endl;
    std::cout<<method_times[0]<<" "<<method_times[1]<<" "<<method_times[2]<<std::endl;
}

void balance_solver::balance_state(VectorXs &pos, VectorXs &force_on_vertex)
{
    MatrixXs force;
    force.resize(f.size()/3,3);
    for(int i=0;i<f.size()/3;++i)
    {
        for(int j=0;j<3;++j)
        {
            force(i,j)=f(i*3+j);
        }
    }
    Vector3s vij;
    for(int i=0;i<e.rows();++i)
    {
        int id0=e(i,0)*3;
        int id1=e(i,1)*3;
        vij<<v_new(id0)-v_new(id1),v_new(id0+1)-v_new(id1+1),
        v_new(id0+2)-v_new(id1+2);
        vij*=s*(1-el(i)/vij.norm());
        force.row(e(i,0))-=vij;
        force.row(e(i,1))+=vij;
    }
    force_on_vertex=force.rowwise().norm();
}

void balance_solver::get_energy(VectorXs &pos,scalar &out_energy, bool with_mass_matrix)
{
    if (with_mass_matrix)
    {
        out_energy=(pos-y).squaredNorm()*m*0.5;
    }
    else
    {
        out_energy=0;
    }
    scalar h2s=h*h*s;
    Vector3s vij;
    for(int i=0;i<e.rows();++i)
    {
        int id0=e(i,0);
        int id1=e(i,1);
        vij<<pos(id0*3)-pos(id1*3),pos(id0*3+1)-pos(id1*3+1),pos(id0*3+2)-pos(id1*3+2);
        scalar delta=vij.norm()-el(i);
        out_energy+=0.5*h2s*delta*delta;
    }
    scalar Work=pos.transpose()*f;
    out_energy-=Work*h*h;
}

void balance_solver::compute_jacobi()
{
    vtp tri;
    tri.reserve(e.rows()*36+3);
    MatrixX3s vjk;
    vjk.resize(3,3);
    for(int i=0;i<e.rows();++i)
    {
        int id0=e(i,0);
        int id1=e(i,1);
        //std::cout<<i<<" "<<id0<<" "<<id1<<std::endl;
        Vector3s v01(v(id1*3)-v(id0*3),v(id1*3+1)-v(id0*3+1),v(id1*3+2)-v(id0*3+2));
        scalar len_inv=1.0/v01.norm();
        scalar len_inv3=el(i)*len_inv*len_inv*len_inv;
        len_inv*=el(i);
        vjk<<1.0-len_inv+v01(0)*v01(0)*len_inv3,             v01(0)*v01(1)*len_inv3,             v01(0)*v01(2)*len_inv3,
                         v01(1)*v01(0)*len_inv3, 1.0-len_inv+v01(1)*v01(1)*len_inv3,             v01(1)*v01(2)*len_inv3,
                         v01(2)*v01(0)*len_inv3,             v01(2)*v01(1)*len_inv3, 1.0-len_inv+v01(2)*v01(2)*len_inv3;
        auto zer=[&](int i0,int i1){if(i0==a||i1==a) return 0.0;else return 1.0;};
        for(int j=0;j<3;++j)
        {
            for(int k=0;k<3;++k)
            {
                tri.emplace_back(id0*3+j,id0*3+k,vjk(j,k)*zer(id0,id0));
                tri.emplace_back(id0*3+j,id1*3+k,-vjk(j,k)*zer(id0,id1));
                tri.emplace_back(id1*3+j,id0*3+k,-vjk(j,k)*zer(id1,id0));
                tri.emplace_back(id1*3+j,id1*3+k,vjk(j,k)*zer(id1,id1));
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
