

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


    auto P = create_eigen_operator(a,one,
                    1.0,
            );
    auto Pt = create_eigen_operator(one,b,
                    1.0,
            );

    auto Zero = create_eigen_operator(one,one, 0.);

    auto W = create_block_eigen_operator<2,2>(G, P,
                                              Pt,Zero);


    Eigen::VectorXd phi(knots.size()+1), gamma;
    for (int i=0; i<knots.size(); ++i) {
        const double x = get<position>(knots[i])[0];
        const double y = get<position>(knots[i])[1];
        if (get<boundary>(knots[i])) {
            phi[i] = funct(x,y);
        } else {
            phi[i] = laplace_funct(x,y);
        }
    }
    phi[knots.size()] = 0;

    std::cout << std::endl;
   

    Eigen::GMRES<decltype(W), Eigen::DiagonalPreconditioner<double>> gmres;
    gmres.set_restart(restart);
    gmres.setMaxIterations(max_iter);
    gmres.compute(W);
    gamma = gmres.solve(phi);
    std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;


    phi = W*gamma;
    for (int i=0; i<knots.size(); ++i) {
        const double x = get<position>(knots[i])[0];
        const double y = get<position>(knots[i])[1];
        if (get<boundary>(knots[i])) {
            TS_ASSERT_DELTA(phi[i],funct(x,y),2e-3); 
        } else {
            TS_ASSERT_DELTA(phi[i],laplace_funct(x,y),2e-3); 
        }
    }
    TS_ASSERT_DELTA(phi[knots.size()],0,2e-3); 


    // This could be more intuitive....
    Eigen::Map<Eigen::Matrix<double,N,1>> alpha_wrap(get<alpha>(knots).data());
    alpha_wrap = gamma.segment<N>(0);

    const double beta = gamma(knots.size());

    interp[a] = sum(b,true,al[b]*kernel) + beta;

    double rms_error = 0;
    double scale = 0;
    for (int i=0; i<knots.size(); ++i) {
        const double x = get<position>(knots[i])[0];
        const double y = get<position>(knots[i])[1];
        const double truth = funct(x,y);
        const double eval_value = get<interpolated>(knots[i]);
        rms_error += std::pow(eval_value-truth,2);
        scale += std::pow(truth,2);
        TS_ASSERT_DELTA(eval_value,truth,1e-2); 
    }
    std::cout << "rms_error for global support, at centers  = "<<std::sqrt(rms_error/scale)<<std::endl;
    TS_ASSERT_LESS_THAN(std::sqrt(rms_error/scale),1e-3);

}
