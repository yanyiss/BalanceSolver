#include "balance_solver.h"
#include <set>
#include <iomanip>
#include <omp.h>
#include <ctime>
#include <fstream>


void balance_solver::predecomposition()
{
    int n=v.size();
    vtp tri;
    tri.reserve(n+spring.rows()*12+a.size());
    //inertia
    for(int i=0;i<n;++i)
    {
        tri.emplace_back(i,i,m);
    }
    //spring force
    for(int i=0;i<spring.rows();++i)
    {
        scalar h2_s=h*h*stiffness(i);
        int id0=spring(i,0)*3;
        int id1=spring(i,1)*3;
        for(int j=0;j<3;++j)
        {
            tri.emplace_back(id0+j,id0+j,h2_s);
            tri.emplace_back(id1+j,id1+j,h2_s);
            tri.emplace_back(id0+j,id1+j,-h2_s);
            tri.emplace_back(id1+j,id0+j,-h2_s);
        }
    }
    //fixed coordinate
    for(int i=0;i<a.size();++i)
    {
        tri.emplace_back(a(i),a(i),h*h*fixed_stiffness);
    }
    spma Mt2L;
    Mt2L.resize(n,n);
    Mt2L.setFromTriplets(tri.begin(),tri.end());
    pd_solver.compute(Mt2L);

    tri.clear();
    for(int i=0;i<spring.rows();++i)
    {
        int id0=spring(i,0)*3;
        int id1=spring(i,1)*3;
        scalar s=stiffness(i);
        for(int j=0;j<3;++j)
        {
            tri.emplace_back(id0+j,i*3+j,s);
            tri.emplace_back(id1+j,i*3+j,-s);
        }
    }
    for(int i=0;i<a.size();++i)
    {
        tri.emplace_back(a(i),spring.rows()*3+i,fixed_stiffness);
    }
    J.resize(n,spring.rows()*3+a.size());
    J.setFromTriplets(tri.begin(),tri.end());
}

void balance_solver::local_step(VectorXs &pos)
{
    for(int i=0;i<spring.rows();++i)
    {
        int id0=spring(i,0)*3;
        int id1=spring(i,1)*3;
        Vector3s vij;
        vij<<pos(id0)-pos(id1),pos(id0+1)-pos(id1+1),
        pos(id0+2)-pos(id1+2);
        vij.normalize();
        for(int j=0;j<3;++j)
        {
            d(i*3+j)=vij(j)*spring_len(i);
        }
    }
    for(int i=0;i<a.size();++i)
    {
        d(i+spring.rows()*3)=a_value(i);
    }
}

void balance_solver::global_step()
{
    v_list.push_back(std::move(v_new));
    v_new=pd_solver.solve(m*y+h*h*(J*d+f));
}


void balance_solver::get_energy(VectorXs &pos,scalar &out_energy,bool with_mass)
{
    if(with_mass)
    {
        out_energy=(pos-y).squaredNorm()*m*0.5;
    }
    else
    {
        out_energy=0.0;
    }
    Vector3s vij;
    for(int i=0;i<spring.rows();++i)
    {
        int id0=spring(i,0)*3;
        int id1=spring(i,1)*3;
        vij<<pos(id0)-pos(id1),pos(id0+1)-pos(id1+1),pos(id0+2)-pos(id1+2);
        scalar delta=vij.norm()-spring_len(i);
        out_energy+=0.5*h*h*stiffness(i)*delta*delta;
    }
    for(int i=0;i<a.size();++i)
    {
        scalar delta=pos(a(i))-a_value(i);
        out_energy+=0.5*h*h*fixed_stiffness*delta*delta;
    }
    scalar work=pos.transpose()*f;
    out_energy-=work*h*h;
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
        }
        else
        {
            v_new=std::move(v_aa);
            energy=aa_energy;
        }
        
        global_step();
    }
    v_list.clear();
    y=v_new+(v_new-v)*0.1;
    v=v_new;
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

bool balance_solver::newton_raphson_linesearch(VectorXs &pos, VectorXs &dir, VectorXs &new_pos)
{
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
    }while(times>1.0e-4&&step_energy>current_energy);
    if(times>1.0e-4)
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
    }while(times>1.0e-4&&step_energy>current_energy);
    if(times>1.0e-4)
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


void balance_solver::newton_raphson()
{
    VectorXs gradient;
    gradient.resize(v.size());
    gradient.setZero();
    vtp tri;
    Vector3s v01;
    MatrixX3s vijk;
    vijk.resize(3,3);
    for(int i=0;i<spring.rows();++i)
    {
        int id0=spring(i,0);
        int id1=spring(i,1);
        v01<<v(id0*3)-v(id1*3),v(id0*3+1)-v(id1*3+1),v(id0*3+2)-v(id1*3+2);
        scalar s=stiffness(i);
        Vector3s g01=s*(v01.norm()-spring_len(i))*v01.normalized();
        for(int j=0;j<3;++j)
        {
            gradient(id0*3+j)+=g01(j);
            gradient(id1*3+j)-=g01(j);
        }

        scalar len_inv=1.0/v01.norm();
        scalar len_inv3=spring_len(i)*len_inv*len_inv*len_inv;
        len_inv*=spring_len(i);
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
    //fixed point
    for(int i=0;i<a.size();++i)
    {
        gradient(a(i))+=fixed_stiffness*(v(a(i))-a_value(i));
        tri.emplace_back(a(i),a(i),fixed_stiffness);
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
    if(newton_raphson_solver.info()!=0)
    {
        std::cout<<"newton raphson hessian factorize failed"<<std::endl;
        for(int i=0;i<5;++i)
        {
            pd();
        }
        return;
    }

    VectorXs dir=-newton_raphson_solver.solve(gradient);
    if(newton_raphson_linesearch(v,dir,v_new))
    {
        y=v_new+(v_new-v)*0.1;
        v=v_new;
    }
    else
    {
        for(int i=0;i<5;++i)
        {
            pd();
        }
    }
}


void balance_solver::newton()
{
    VectorXs gradient;
    gradient.resize(v.size());
    gradient.setZero();
    vtp tri;
    Vector3s v01;
    MatrixX3s vijk;
    vijk.resize(3,3);
    for(int i=0;i<spring.rows();++i)
    {
        scalar h2s=h*h*stiffness(i);
        int id0=spring(i,0);
        int id1=spring(i,1);
        v01<<v(id0*3)-v(id1*3),v(id0*3+1)-v(id1*3+1),v(id0*3+2)-v(id1*3+2);
        Vector3s g01=stiffness(i)*(v01.norm()-spring_len(i))*v01.normalized();
        for(int j=0;j<3;++j)
        {
            gradient(id0*3+j)+=g01(j);
            gradient(id1*3+j)-=g01(j);
        }

        scalar len_inv=1.0/v01.norm();
        scalar len_inv3=spring_len(i)*len_inv*len_inv*len_inv;
        len_inv*=spring_len(i);
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
    //fixed point
    scalar h2s=h*h*fixed_stiffness;
    for(int i=0;i<a.size();++i)
    {
        gradient(a(i))+=fixed_stiffness*(v(a(i))-a_value(i));
        tri.emplace_back(a(i),a(i),h2s);
    }
    gradient-=f;
    gradient=m*(v-y)+h*h*gradient;

    for(int i=0;i<v.size();++i)
    {
        tri.emplace_back(i,i,m);
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
    if(newton_solver.info()!=0)
    {
        std::cout<<"newton hessian factorize failed"<<std::endl;
        for(int i=0;i<5;++i)
        {
            pd();
        }
        return;
    }

    VectorXs dir=-newton_solver.solve(gradient);
    if(newton_linesearch(v,dir,v_new))
    {
        y=v_new+(v_new-v)*0.1;
        v=v_new;
    }
    else
    {
        for(int i=0;i<10;++i)
        {
            pd();
        }
    }
}

void balance_solver::compute_balance()
{
    VectorXs force_on_vertex;
    int itertimes=0;
    do
    {
        clock_t start=clock();
        
        balance_state(v_new,force_on_vertex);
        balance_result=force_on_vertex.maxCoeff();
        //std::cout<<"balance result: "<<balance_result<<std::endl;
        ++itertimes;
        if(itertimes>compute_balance_times||balance_result<balance_threshold)
        { 
            if(itertimes==1)
            {
                x_invariant=true;
            }
            std::cout<<"balance result: "<<balance_result<<std::endl;
            return;
        }
        x_invariant=false;
        if(balance_result>newton_rate)
        {
            pd();
        }
        else
        {
            newton();
        }
    }while(1);
}

void balance_solver::balance_state(VectorXs &pos, VectorXs &force_on_vertex)
{
    MatrixXs force;
    int n=f.size()/3;
    force.resize(n,3);
    for(int i=0;i<n;++i)
    {
        for(int j=0;j<3;++j)
        {
            force(i,j)=f(i*3+j);
        }
    }
    Vector3s vij;
    for(int i=0;i<spring.rows();++i)
    {
        int id0=spring(i,0)*3;
        int id1=spring(i,1)*3;
        vij<<pos(id0)-pos(id1),pos(id0+1)-pos(id1+1),
        pos(id0+2)-pos(id1+2);
        vij*=stiffness(i)*(1.0-spring_len(i)/vij.norm());
        force.row(spring(i,0))-=vij;
        force.row(spring(i,1))+=vij;
    }
    for(int i=0;i<a.size();++i)
    {
        force(a(i)/3,a(i)%3)+=fixed_stiffness*(a_value(i)-pos(a(i)));
    }
    // for(int i=0;i<n;++i)
    // {
    //     std::cout<<i<<" "<<force(i,0)<<" "<<force(i,1)<<" "<<force(i,2)<<std::endl;
    // }
    force_on_vertex=force.rowwise().squaredNorm();
}


// void balance_solver::compute_csr()
// {
//     jacobirow.resize(spring.rows()*36+a.size());
//     jacobicol.resize(spring.rows()*36+a.size());
//     jacobival.resize(spring.rows()*36+a.size());
//     MatrixX3s vjk;
//     vjk.resize(3,3);
//     auto assemble_tri=[&](int id,int r,int c,scalar v)
//     {
//         jacobirow(id)=r;jacobicol(id)=c;jacobival(id)=v;
//     };
//     int count=0;
//     for(int i=0;i<spring.rows();++i)
//     {
//         int id0=spring(i,0);
//         int id1=spring(i,1);
//         //std::cout<<i<<" "<<id0<<" "<<id1<<std::endl;
//         Vector3s v01(v(id1*3)-v(id0*3),v(id1*3+1)-v(id0*3+1),v(id1*3+2)-v(id0*3+2));
//         scalar len_inv=1.0/v01.norm();
//         scalar len_inv3=spring_len(i)*len_inv*len_inv*len_inv;
//         len_inv*=spring_len(i);
//         vjk<<1.0-len_inv+v01(0)*v01(0)*len_inv3,             v01(0)*v01(1)*len_inv3,             v01(0)*v01(2)*len_inv3,
//                          v01(1)*v01(0)*len_inv3, 1.0-len_inv+v01(1)*v01(1)*len_inv3,             v01(1)*v01(2)*len_inv3,
//                          v01(2)*v01(0)*len_inv3,             v01(2)*v01(1)*len_inv3, 1.0-len_inv+v01(2)*v01(2)*len_inv3;
//         vjk*=stiffness(i)/fixed_stiffness;
//         auto is_in_a=[&](int i0,int i1)
//         {
//             for(int i=0;i<a.size();++i)
//             {
//                 if(i0==a(i)||i1==a(i))
//                 {
//                     return 0.0;
//                 }
//             }
//             return 1.0;
//         };
//         for(int j=0;j<3;++j)
//         {
//             for(int k=0;k<3;++k)
//             {
//                 assemble_tri(count++,id0*3+j,id0*3+k,vjk(j,k)*is_in_a(id0*3+j,id0*3+k));
//                 assemble_tri(count++,id0*3+j,id1*3+k,-vjk(j,k)*is_in_a(id0*3+j,id1*3+k));
//                 assemble_tri(count++,id1*3+j,id0*3+k,-vjk(j,k)*is_in_a(id1*3+j,id0*3+k));
//                 assemble_tri(count++,id1*3+j,id1*3+k,vjk(j,k)*is_in_a(id1*3+j,id1*3+k));
//             }
//         }
//     }
//     for(int i=0;i<a.size();++i)
//     {
//         assemble_tri(count++,a(i),a(i),1.0);
//     }
// }


void balance_solver::compute_csr()
{
    jacobirow.resize(spring.rows()*36+a.size());
    jacobicol.resize(spring.rows()*36+a.size());
    jacobival.resize(spring.rows()*36+a.size());
    MatrixX3s vjk;
    vjk.resize(3,3);
    auto assemble_tri=[&](int id,int r,int c,scalar v)
    {
        jacobirow(id)=r;jacobicol(id)=c;jacobival(id)=v;
    };
    int count=0;
    for(int i=0;i<spring.rows();++i)
    {
        int id0=spring(i,0);
        int id1=spring(i,1);
        //std::cout<<i<<" "<<id0<<" "<<id1<<std::endl;
        Vector3s v01(v(id1*3)-v(id0*3),v(id1*3+1)-v(id0*3+1),v(id1*3+2)-v(id0*3+2));
        scalar len_inv=1.0/v01.norm();
        scalar len_inv3=spring_len(i)*len_inv*len_inv*len_inv;
        len_inv*=spring_len(i);
        vjk<<1.0-len_inv+v01(0)*v01(0)*len_inv3,             v01(0)*v01(1)*len_inv3,             v01(0)*v01(2)*len_inv3,
                         v01(1)*v01(0)*len_inv3, 1.0-len_inv+v01(1)*v01(1)*len_inv3,             v01(1)*v01(2)*len_inv3,
                         v01(2)*v01(0)*len_inv3,             v01(2)*v01(1)*len_inv3, 1.0-len_inv+v01(2)*v01(2)*len_inv3;
        vjk*=stiffness(i)/fixed_stiffness;
        for(int j=0;j<3;++j)
        {
            for(int k=0;k<3;++k)
            {
                assemble_tri(count++,id0*3+j,id0*3+k,vjk(j,k));
                assemble_tri(count++,id0*3+j,id1*3+k,-vjk(j,k));
                assemble_tri(count++,id1*3+j,id0*3+k,-vjk(j,k));
                assemble_tri(count++,id1*3+j,id1*3+k,vjk(j,k));
            }
        }
    }
    for(int i=0;i<a.size();++i)
    {
        assemble_tri(count++,a(i),a(i),1.0);
    }
}

void balance_solver::compute_csr_right()
{
    jacobiright.resize(v.size(),rope_index.size()*3);
    jacobiright.setZero();
    scalar s_inv=1.0/fixed_stiffness;
    for(int i=0;i<rope_index.size();++i)
    {
        for(int j=0;j<3;++j)
        {
            jacobiright(rope_index(i)*3+j,i*3+j)=s_inv;
        }
    }
}

void balance_solver::cubic_sampling()
{
    scalar sinloc_cosdir=sin(loc_rad)*cos(dir_rad);
    scalar cosloc_sindir=cos(loc_rad)*sin(dir_rad);
    scalar cosloc_cosdir=cos(loc_rad)*cos(dir_rad);
    scalar sinloc_sindir=sin(loc_rad)*sin(dir_rad);
    scalar cosdir=cos(dir_rad);

    int n=10000;
    scalar cubic_len=0;
    VectorXs sampling_len;sampling_len.resize(n);
    VectorXs x;x.resize(n+1);
    VectorXs y;y.resize(n+1);
    x(0)=(sinloc_cosdir*cos(start_rad)-cosloc_sindir)/(cosloc_cosdir*cos(start_rad)+sinloc_sindir);
    y(0)=cosdir*sin(start_rad)/(cosloc_cosdir*cos(start_rad)+sinloc_sindir);
    scalar step=(end_rad-start_rad)/n;
    for(int i=1;i<=n;++i)
    {
        scalar lambda=i*step+start_rad;
        x(i)=(sinloc_cosdir*cos(lambda)-cosloc_sindir)/(cosloc_cosdir*cos(lambda)+sinloc_sindir);
        y(i)=cosdir*sin(lambda)/(cosloc_cosdir*cos(lambda)+sinloc_sindir);
        scalar xi=x(i)-x(i-1);
        scalar yi=y(i)-y(i-1);
        sampling_len(i-1)=sqrt(xi*xi+yi*yi);
        cubic_len+=sampling_len(i-1);
    }

    int count=1;
    //step=(end_rad-start_rad)/(sampling_num-1);
    scalar acc_len=0;
    scalar avg_len=cubic_len/(sampling_num-1);
    sampling_lambda.resize(sampling_num);
    sampling_lambda(0)=start_rad;
    for(int i=1;i<=n;++i)
    {
        acc_len+=sampling_len(i-1);
        scalar t=acc_len-count*avg_len;
        if(t>0)
        {
            t/=sampling_len(i-1);
            sampling_lambda(count)=t*((i-1)*step+start_rad)+(1-t)*(i*step+start_rad);
            //sampling_lambda(count)=compute_lambda(t*x(i-1)+(1-t)*x(i),t*y(i-1)+(1-t)*y(i));
            ++count;
        }
    }
    sampling_lambda(sampling_num-1)=end_rad;
}

scalar stickLooker::compute_hanging_distance(MatrixX3s &point,Node* node)
{
    int prev=node->prev->id;
    int suc=node->suc->id;
    int np=point.rows();
    if(prev>suc) suc+=np;
    
 
    std::vector<int> support;
    support.push_back(prev);
    for(int i=prev;i<suc;)
    {
        int j=i+1;
        int k=j+1;
        scalar theta=0.0;
        scalar max_theta=0.0;
        int mark=j;
        auto calc_angle=[&]()
        {
            Vector3s vij=point.row(j%np)-point.row(i%np);
            Vector3s vik=point.row(k%np)-point.row(i%np);
            return std::atan2(vij.cross(vik)(2),vij.dot(vik));
        };
        while(k<=suc)
        {
            theta+=calc_angle();
            if(theta>max_theta)
            {
                max_theta=theta;
                mark=k;
            }
            ++j;++k;
        }
        support.push_back(mark);
        i=mark;
    }
    scalar max_dis=0.0;
    for(int i=0;i<support.size()-1;++i)
    {
        int j=support[i];
        int k=support[i+1];
        Vector3s vjk=(point.row(k%np)-point.row(j%np)).normalized();
        Vector3s vjt;
        for(int t=j+1;t<k;++t)
        {
            vjt=point.row(t%np)-point.row(j%np);
            max_dis=std::max(max_dis,std::fabs(vjt.cross(vjk)(2)));
        }
    }
    return max_dis;
}

void stickLooker::stick_locating(MatrixX3s &point,int stick_num)
{
    int np=point.rows();
    init(np);
    for(int i=0;i<np;++i)
    {
        pq.emplace(i,++loop[i].count,compute_hanging_distance(point,&loop[i]));
    }
    while(loopsize>stick_num)
    {
        NodeDist nd;
        do
        {
            nd=pq.top();
            pq.pop();
        } while (nd.count!=loop[nd.id].count);
        Node* node=&loop[nd.id];
        delete_node(nd.id);
        pq.emplace(node->prev->id,++node->prev->count,compute_hanging_distance(point, node->prev));
        pq.emplace(node->suc->id,++node->suc->count,compute_hanging_distance(point, node->suc));
    }
    get_stick_index();
}

