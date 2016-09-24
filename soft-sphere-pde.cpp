

int main(void) {
    ABORIA_VARIABLE(density1p,double,"density")
    ABORIA_VARIABLE(density2p,double,"two particle density")
        ABORIA_VARIABLE(constant2,double,"c2 value")

    typedef Particles<std::tuple<density1p,density2p,constant2>,2> ParticlesType;
    typedef position_d<2> position;
    ParticlesType knots;

    const double c = 0.5;
    const int max_iter = 100;
    const int restart = 100;
    double2 periodic(false);
    
    const int nx = 7;
    constexpr int N = (nx+1)*(nx+1);
    const double delta = 1.0/nx;
    typename ParticlesType::value_type p;
    for (int i=0; i<=nx; ++i) {
        for (int j=0; j<=nx; ++j) {
            get<position>(p) = double2(i*delta,j*delta);
            if ((i==0)||(i==nx)||(j==0)||(j==nx)) {
                get<boundary>(p) = true;
            } else {
                get<boundary>(p) = false;
            }
            get<constant2>(p) = std::pow(c,2);
            knots.push_back(p);
        }
    }

    Symbol<position> r;
    Symbol<density1p> p;
    Symbol<density2p> P;
    Symbol<constant2> c2;
    Label<0,ParticlesType> a(knots);
    Label<1,ParticlesType> b(knots);
    One one;
    auto dx = create_dx(a,b);
    Accumulate<std::plus<double> > sum;

    auto kernel = deep_copy(
            exp(-pow(norm(dx),2)/c2[b])
            );

    auto laplace_kernel = deep_copy(
            //(2*c2[b] + pow(norm(dx),2)) / pow(pow(norm(dx),2) + c2[b],1.5)
            4*(pow(norm(dx),2) - c2[b]) * exp(-pow(norm(dx),2)/c2[b])/pow(c2[a],2)
            );

    auto G = create_eigen_operator(a,b, 
                if_else(is_b[a],
                    kernel,
                    laplace_kernel
                )
            );
    auto P = create_eigen_operator(a,one,
                if_else(is_b[a],
                    1.0,
                    0.0
                )
            );
    auto Pt = create_eigen_operator(one,b,
                if_else(is_b[b],
                    1.0,
                    0.0
                )
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
