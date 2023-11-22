#include "pd.h"
#include <set>
#include <iomanip>

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
            d(i*3+j)=vij(j)*el(i);
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

bool DiffSimulation::linesearch(VectorXs &pos, VectorXs &dir, VectorXs &new_pos)
{
#if 1
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
    }while(times>1.0e-8&&step_energy>current_energy);
    
    if(times>1.0e-8)
    {
        new_pos=std::move(v_tdv);
        return true;
    }
    else
    {
        return false;
    }
#else
    scalar current_rate=0;
    print_balance_info(balance_cof);
    current_rate=balance_rate;

    scalar times=2.0;
    VectorXs v_tdv;
    scalar step_rate;
    do
    {
        times*=0.5;
        v_tdv=pos+times*dir;
        print_balance_info(balance_cof);
        step_rate=balance_rate;
    }while(times>1.0e-8&&step_rate>current_rate);
    if(times>1.0e-8)
    {
        new_pos=std::move(v_tdv);
        return true;
    }
    else
    {
        return false;
    }
#endif
}

void DiffSimulation::trans_fix()
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

void DiffSimulation::newton()
{
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

    for(int i=0;i<3;++i)
    {
        tri.emplace_back(a*3+i,a*3+i,h2s);
    }
    for(int i=0;i<v.size();++i)
    {
        tri.emplace_back(i,i,m);
    }
    spma hessian;
    hessian.resize(v.size(),v.size());
    hessian.setFromTriplets(tri.begin(),tri.end());
    if(!hessian_solver_info)
    {
        hessian_solver.analyzePattern(hessian);
        hessian_solver_info=true;
    }
    hessian_solver.factorize(hessian);
    if(hessian_solver.info()!=0)
    {
        std::cout<<"hessian factorize failed"<<std::endl;
        exit(0);
    }
    VectorXs dir=-hessian_solver.solve(gradient);
    if(linesearch(v,dir,v_new))
    {
        trans_fix();
    }
    else
    {
        for(int i=0;i<10;++i)
        {
            Opt();
        }
    }
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
    v_new=v;
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
    trans_fix();
}

void DiffSimulation::print_balance_info(scalar cof)
{
    balance_cof=cof;
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
        vij*=s*(1-el(i)/vij.norm());
        balance.row(e(i,0))-=vij;
        balance.row(e(i,1))+=vij;
    }
    VectorXs Norm=balance.rowwise().norm();
    std::cout<<"\nBalance Info:\n";
    scalar meanNorm=Norm.sum()/Norm.size();
    std::cout<<"Mean Norm: "<<meanNorm<<std::endl;
    std::cout<<"Max Norm: "<<Norm.maxCoeff()<<std::endl;
    /* std::cout<<"balance info:\n";
    std::cout<<"resultant force: "<<f.sum()<<std::endl;
    std::cout<<"Frobenius Norm: "<<balance.norm()<<std::endl;
    VectorXs norm=balance.rowwise().norm();
    std::cout<<"Mean Norm: "<<norm.sum()/norm.size()<<std::endl;
    std::cout<<"Max Norm: "<<norm.maxCoeff(); */
    /* VectorXs mean1=balance.rowwise().norm();
    for(int i=0;i<mean1.size();++i)
    {
        std::cout<<i<<" "<<mean1(i)<<std::endl;
    } */
    balance_rate=meanNorm/(9.8*m);
    /* return meanNorm/(9.8*m);
    if(meanNorm<9.8*m*cof)
    {
        return true;
    }
    else
    {
        return false;
    } */
    /* for(int i=0;i<norm.size();++i)
    {
        std::cout<<i<<" "<<norm(i)<<std::endl;
    } */
}

void DiffSimulation::get_energy(VectorXs &pos,scalar &out_energy)
{
    /* energy= 0.5*pos.transpose()*(Mt2L*pos);
    scalar e_temp=pos.transpose()*(J*d+f);
    energy-=e_temp*h*h;
    e_temp=-pos.transpose()*y;
    energy-=e_temp*m;
    energy+=0.5*s*h*h*d.squaredNorm(); */
    out_energy=(pos-y).squaredNorm()*m*0.5;
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

void DiffSimulation::get_energy_with_fixed(VectorXs &pos,scalar &out_energy)
{
    get_energy(pos,out_energy);
    Vector3s v01;
    v01<<pos(a*3)-a_pos(0),pos(a*3+1)-a_pos(1),pos(a*3+2)-a_pos(2);
    out_energy+=0.5*h*h*s*v01.squaredNorm();
}

void DiffSimulation::get_gradient_with_fixed(VectorXs &pos,std::vector<scalar> &gradient)
{
    gradient.resize(v.size(),0);
    scalar h2s=h*h*s;
    Vector3s v01;
    for(int i=0;i<e.rows();++i)
    {
        int id0=e(i,0);
        int id1=e(i,1);
        v01<<pos(id0*3)-pos(id1*3),pos(id0*3+1)-pos(id1*3+1),pos(id0*3+2)-pos(id1*3+2);
        Vector3s g01=h2s*(v01.norm()-el(i))*v01.normalized();
        for(int j=0;j<3;++j)
        {
            gradient[id0*3+j]+=g01(j);
            gradient[id1*3+j]-=g01(j);
        }
    }
    for(int i=0;i<f.size();++i)
    {
        gradient[i]+=-h*h*f(i)+m*(pos(i)-y(i));
    }
    for(int i=0;i<3;++i)
    {
        gradient[a*3+i]+=h2s*(pos(a*3+i)-a_pos(i));
    }
}

void DiffSimulation::energy_and_gradient_evaluation(const std::vector<scalar>& x_evaluation, scalar& energy_evaluation, 
                                   std::vector<scalar>& gradient_evaluation)
{
    VectorXs pos;
    pos.resize(v.size());
    for(int i=0;i<v.size();++i)
    {
        pos(i)=x_evaluation[i];
    }
    get_energy_with_fixed(pos,energy_evaluation);
    get_gradient_with_fixed(pos,gradient_evaluation);
}

void DiffSimulation::LBFGS()
{
    HLBFGS hlbfgs_solver;
	hlbfgs_solver.set_number_of_variables(v.size());
	hlbfgs_solver.set_verbose(true);
	hlbfgs_solver.set_func_callback(func_callback, 0, 0, 0, 0);
    v_new=v;
	hlbfgs_solver.optimize_without_constraints(&v_new(0), 100, this);
    trans_fix();
}

void DiffSimulation::compute_jacobi()
{
#if 0
    std::cout<<"info\n"<<std::endl;
    std::cout<<f.norm()<<" "<<y.norm()<<" "<<v.norm()<<std::endl;
    jacobi.resize(v.size(),in.size()*3-3);
    jacobi.setZero();
    VectorXs f_temp=f;
    VectorXs v_stable=v;
    VectorXs y_temp=y;
    scalar dif=1e-3;
    int inback=in(in.size()-1);
    //for(int i=0;i<in.size()-1;++i)
    for(int i=0;i<3;++i)
    {
        for(int j=0;j<3;++j)
        {
            f(in(i)*3+j)+=dif;
            f(inback*3+j)-=dif;
            do
            {
#if 0
                Opt();
#else
                if(balance_rate>1e-4)
                Opt();
                else
                newton();
#endif
                print_balance_info(balance_cof);
            }while(balance_rate>balance_cof);

            VectorXs v_dif=(v-v_stable)/dif;
            for(int k=0;k<v.size();++k)
            {
                jacobi(k,i*3+j)=v_dif(k);
            }

            f=f_temp;
            v=v_stable;
            y=y_temp;
        }
    }
    std::cout<<"jacobi\n";
    std::cout<<jacobi.norm()<<std::endl;
    return;
 #endif

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
        /* for(int j=0;j<3;++j)
        {
            tri.emplace_back(a*3+i,a*3+j,1.0);
        } */
        tri.emplace_back(a*3+i,a*3+i,1.0);
    }
    spma v_mat;
    v_mat.resize(v.size(),v.size());
    v_mat.setFromTriplets(tri.begin(),tri.end());
    /* std::cout<<"v_mat!!!!!!!!!!!!!!!!!!!!!\n";
    for(int i=0;i<v_mat.outerSize();++i)
    {
        for(spma::InnerIterator it(v_mat,i);it;++it)
        {
            std::cout<<it.row()<<" "<<it.col()<<" "<<it.value()<<std::endl;
        }
    } */
    
    if(!v_solver_info)
    {
#if 1
        delta.resize(v.size(),in.size()*3-3);
        delta.setZero();
        scalar s_inv=1.0/s;
        for(int i=0;i<in.size()-1;++i)
        {
            for(int j=0;j<3;++j)
            {
                delta(in(i)*3+j,i*3+j)=s_inv;
                delta(in(in.size()-1)*3+j,i*3+j)=-s_inv;
            }
        }
        /* for(int i=0;i<3;++i)
        {
            delta(in(0)*3+i,0*3+i)=s_inv;
            for(int j=1;j<in.size()-1;++j)
            {
                delta(in(j-1)*3+i,(j-1)*3+i)=-s_inv;
                delta(in(j)*3+i,j*3+i)=s_inv;
            }
            delta(in(in.size()-2)*3+i,(in.size()-2)*3+i)=-s_inv;
        } */
#else
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
#endif
        v_solver.analyzePattern(v_mat);
        v_solver_info=true;
    }
    v_solver.factorize(v_mat);
    /* std::cout<<v_solver.info()<<std::endl;
    std::cout<<v_solver.determinant()<<std::endl; */
    
    jacobi.resize(v.size(),in.size()*3);
/* #pragma omp parallel for num_threads(6)
    for(int i=0;i<in.size()*3;++i)
    {
        jacobi.col(i)=v_solver.solve(delta.col(i));
    } */
    jacobi=v_solver.solve(delta);

    ///std::cout<<"vmat\n";
    //std::cout<<v_mat*jacobi-delta<<std::endl;
    leftmat=v_mat;

    /* int adj[3]={1,4,6};
    scalar len[3]={2.7,1.684399,1.9885};
    MatrixX3s sum;
    sum.resize(3,3);
    sum.setZero();
    MatrixX3s id;
    id.resize(3,3);
    id.setZero();
    id(0,0)=1.0;id(1,1)=1.0;id(2,2)=1.0;
    for(int i=0;i<3;++i)
    {
        Vector3s v01(v(adj[i]*3)-v(0),v(adj[i]*3+1)-v(1),v(adj[i]*3+2)-v(2));
        scalar norm=v01.norm();
        MatrixX3s mat=id*(1-len[i]/norm)+len[i]*v01*v01.transpose()/(norm*norm*norm);
        sum+=mat*(jacobi.block(adj[i]*3,0,3,3)-jacobi.block(0,0,3,3));
        std::cout<<"\n"<<i<<std::endl;
        std::cout<<mat<<std::endl;
        std::cout<<jacobi.block(adj[i]*3,0,3,3)-jacobi.block(0,0,3,3)<<std::endl;
    }
    std::cout<<sum<<std::endl; */
}

void func_callback(const size_t N, const std::vector<scalar>& x_evaluation, scalar& energy_evaluation, 
                                   std::vector<scalar>& gradient_evaluation, void* user_supply)
{
	auto ptr_this = static_cast<DiffSimulation*>(user_supply);
	ptr_this->energy_and_gradient_evaluation(x_evaluation, energy_evaluation, gradient_evaluation);
}