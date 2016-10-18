

int main(void) {
    ABORIA_VARIABLE(density1p,double,"density")
    ABORIA_VARIABLE(density2p,double,"two particle density")
    ABORIA_VARIABLE(constant2,double,"c2 value")
    ABORIA_VARIABLE(g_function,double,"g function")

    typedef Particles<std::tuple<density1p,density2p,g_function,constant2>,3> ParticlesType;
    typedef position_d<3> position;
    ParticlesType knots;

    const double c = 0.5;
    const double L = 1.0;
    const double3 low(0,0,0);
    const double3 high(L,L,L);
    const int max_iter = 100;
    const int restart = 100;
    double2 periodic(false);
    
    const int nx = 7;
    const int ntheta = 5;
    const int nr = 5;
    const double deltar = 0.1;
    const double factorr = 1.01;
    constexpr int N = (nx+1)*(nx+1);
    const double deltax = std::sqrt(1.0)/nx;
    const double deltatheta = 2*pi/nx;
    typename ParticlesType::value_type p;
    for (int i=0; i<=nx; ++i) {
        const double x = i*deltax;
        for (int j=0; j<=ntheta; ++j) {
            const double theta = j*deltatheta;
            double r = deltar/factorr;
            for (int k=0; j<=nr; ++j) {
                r *= factorr;
                get<position>(p) = double3(x/std::sqrt(2.0),x/std::sqrt(2.0)+r*std::cos(theta),r*std::sin(theta));
                if ((p<low).any() || (p>=high).any()) continue;
                get<constant2>(p) = std::pow(c,2);
                knots.push_back(p);
            }
        }
    }
    std::cout << "added "<<knots.size()<<" knots" << std::endl;
    knots.init_neighbour_search(low,high,std::pow(c,2),bool3(true,true,true));

    Symbol<position> r;
    Symbol<density1p> p;
    Symbol<density2p> P;
    Symbol<g_function> g;
    Symbol<constant2> c2;
    Label<0,ParticlesType> a(knots);
    Label<1,ParticlesType> b(knots);
    One one;
    auto dx = create_dx(a,b);
    Accumulate<std::plus<double> > sum;

    auto kernel_3d = create_eigen_operator(a,b,
            exp(-dot(dx,dx)/c2[b])
            );

    auto kernel_2d_1_3
    auto kernel_2d_2_3
    auto kernel_2d_1_2

    auto f_1_2
    auto f_1_3
    auto f_2_1
    auto f_2_3


    auto kernel = create_eigen_operator(a,b,
            exp(-dot(dx,dx)/c2[b])
            );

    auto dkernel_dx1 = create_eigen_operator(a,b,
            (-2*r[a][0] + 2*r[b][0])*exp(-dot(dx,dx)/c2[b])/c2[b]
            );
    auto dkernel_dx2 = create_eigen_operator(a,b,
            (-2*r[a][1] + 2*r[b][1])*exp(-dot(dx,dx)/c2[b])/c2[b]
            );

    auto ikernel_dx3 = create_eigen_operator(a,b,
            std::sqrt(pi)*c[b]**exp(-(pow(dx[0],2)+pow(dx[1],2))/c2[b])
            );

    auto dkernel_dx1dx1 = create_eigen_operator(a,b,
            (4*pow(r[a][0],2) - 2*c2[b])*exp(-dot(dx,dx)/c2[b])/pow(c2[b],2)
            );

    auto dkernel_dx2dx2 = create_eigen_operator(a,b,
            (4*pow(r[a][1],2) - 2*c2[b])*exp(-dot(dx,dx)/c2[b])/pow(c2[b],2)
            );




    auto P = create_eigen_operator(a,one,
                    1.0,
            );
    auto Pt = create_eigen_operator(one,b,
                    1.0,
            );

    auto Zero = create_eigen_operator(one,one, 0.);

    auto W = create_block_eigen_operator<2,2>(G, P,
                                              Pt,Zero);

    Eigen::Map<Eigen::Matrix<double,n,1> > s_vect(get<scalar>(particles).data());

    template<typename T1, typename T2>
    f(T1 arg1, T2 arg2) {
        return ...
    }

    template<typename K, typename F>
    solve(K kernel, F f) {
        Eigen::Map<Eigen::Matrix<double,n,1> > source(get<temporary1>(knots).data());
        Eigen::Map<Eigen::Matrix<double,n,1> > result(get<temporary2>(knots).data());
        tmp1[a] = f;
        Eigen::GMRES<decltype(W), Eigen::DiagonalPreconditioner<double>> gmres;
        gmres.set_restart(restart);
        gmres.setMaxIterations(max_iter);
        gmres.compute(create_eigen_operator(a,b,kernel));
        result = gmres.solve(source);
        std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
        return tmp2[a];
    }

    for (int i=0; i < timesteps; i++) {
        p[a] = kernel_1*w_p[a];
        p_2[a] = kernel_2*w_p[a];
        p_2_dx2[a] = kernel_2_dx2*w_p[a];
        p_1_2[a] = p[a]*p_2[a];
        p_dx1[a] = kernel_dx1*w_p[a];
        p_d2x1[a] = kernel_d2x1*w_p[a];
        P[a] = kernel_1_2*w_P[a];
        P_dx1[a] = kernel_dx1*w_P[a];
        P_dx2[a] = kernel_dx2*w_P[a];
        P_d2x1[a] = kernel_d2x1*w_P[a];
        P_d2x2[a] = kernel_d2x2*w_P[a];
        g[a] = (kernel_1_3*w_P[a])*(kernel_2_3*w_P[a])/(kernel_3*w_p[a]);

        w_g1[a] = solve(kernel_1_2,f(r[a][0],r[a][1])*P[a]);
        w_g2[a] = solve(kernel_1_2_3,f(r[a][0],r[a][2])*g[a]);
        w_g3[a] = solve(kernel_1_2_3,f(r[a][1],r[a][2])*g[a]);

        g1_ix2[a] = ikernel_dx2*w_g1[a];
        g2_ix3[a] = ikernel_dx3*w_g2[a];
        g3_ix3[a] = ikernel_dx3*w_g3[a];

        g1_ix2_dx1[a] = kernel_ix2_dx1*w_g1[a];
        g2_ix3_dx1[a] = kernel_ix3_dx1*w_g2[a];
        g3_ix3_dx2[a] = kernel_ix3_dx2*w_g3[a];


        kernel*w_p - kernel*w_p_0 = dt*((n-1)*kernel_f_ix2_dx1*w_P[a] + kernel_d2x1*w_p[a])
        dpdt[a] = (n-1)*g1_ix2_dx1[a] + p_d2x1[a];
        p_n+1 - dt*(n-1)*g1_ix2_dx1[a] + p_d2x1[a] = p_n
        kernel*w_p_n+1 - dt*(n-1)*g1_ix2_dx1[a] + p_d2x1[a] = p_n
        p_n+1(I - dt*(n-1)*g1_ix2_dx1[a] + p_d2x1[a] = p_n
        dPdt[a] = -(n-2)*P[a]*g2_ix3[a]*p_dx1[a]/(p[a]*p_1_2[a]) - 
                    (n-2)*P[a]*g3_ix3[a]*p_2_dx2[a]/(p_1_2[a]*p_2[a]) + 
                    (n-2)*P[a]*g2_ix3_dx1[a]/(p_1_2[a]) + 
                    (n-2)*P[a]*g3_ix3_dx2[a]/(p_1_2[a]) + 
                    (n-2)*g2_ix3[a]*P_dx1[a]/(p_1_2[a]) + 
                    (n-2)*g3_ix3[a]*P_dx2[a]/(p_1_2[a]) +
                    P[a]*f_dx1[a] + P[a]*f_dx2[a] + f[a]*P_dx1[a] +
                    f[a]*P_dx2[a] + P_d2x1[a] + P_d2x2[a];
        p[a] += dt*dpdt[a];


}
