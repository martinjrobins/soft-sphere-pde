

int main(void) {
    ABORIA_VARIABLE(density1p,double,"density")
    ABORIA_VARIABLE(density2p,double,"two particle density")
    ABORIA_VARIABLE(constant,double,"constant for basis functions")
    ABORIA_VARIABLE(g_function,double,"g function")

    typedef Particles<std::tuple<density1p,density2p,g_function,constant>,3> ParticlesType;
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
                get<constant>(p) = 1.0/std::pow(c,2);
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
    Symbol<constant> c;
    Label<0,ParticlesType> a(knots);
    Label<1,ParticlesType> b(knots);
    One one;
    auto dx = create_dx(a,b);
    Accumulate<std::plus<double> > sum;

    auto k = deep_copy(
            exp(-dot(dx,dx)*c[b])
            );
    auto f = deep_copy(
            if_else(r[a][0] > r[a][1],
                -epsilon*exp(-epsilon*(r[a][0]-r[a][1])),
                epsilon*exp(-epsilon*(r[a][1]-r[a][0])),
                )
            );


    auto suberf1 = deep_copy(erf(0.5*pow(c[b],-0.5)*(epsilon-2*c[b]*(r[a][0]-r[b][1])))-1.0);
    auto suberf2 = deep_copy(erf(0.5*pow(c[b],-0.5)*(epsilon+2*c[b]*(r[a][0]-r[b][1])))-1.0);
    auto subexp1 = deep_copy(exp(0.25*pow(epsilon-2*c[b]*(r[a][0]-r[b][1]),2)/c[b]));
    auto subexp2 = deep_copy(exp(0.25*pow(epsilon+2*c[b]*(r[a][0]-r[b][1]),2)/c[b]));

    auto K1 = deep_copy(std::sqrt(pi)*epsilon*pow(c2[b],-0.5)*0.5*(
                -(epsilon+2*c[b]*dx[0])*suberf1*subexp1
                -(epsilon-2*c[b]*dx[0])*suberf2*subexp2)
            *exp(-c[b]*(pow(dx[0],2)+pow(r[a][0]-r[b][1],2)));

    auto K2 = deep_copy(2*c[b]*(2*c[b]*pow(dx[0],2)-1)*exp(-c[b]*pow(dx[0],2)));
    auto K3 = deep_copy(exp(-c[b]*(pow(dx[0],2)+pow(dx[1],2))));
    auto K4 = deep_copy(exp(-c[b]*(pow(dx[0],2))));
    auto K5 = deep_copy(exp(-c[b]*(pow(r[b][0]-r[a][1],2))));
    
    auto suberf3 = deep_copy(erf(0.5*pow(c[b],-0.5)*(epsilon-2*c[b]*(r[a][0]-r[b][2])))-1);
    auto suberf4 = deep_copy(erf(0.5*pow(c[b],-0.5)*(epsilon+2*c[b]*(r[a][0]-r[b][2])))-1);
    auto subexp3 = deep_copy(exp(0.25*pow(epsilon-2*c[b]*(r[a][0]-r[b][2]),2)/c[b]));
    auto subexp4 = deep_copy(exp(0.25*pow(epsilon+2*c[b]*(r[a][0]-r[b][2]),2)/c[b]));
    auto subbigexp1 = deep_copy(exp(-c[b]*(pow(dx[0],2)+pow(dx[1],2)+pow(r[a][0]-r[b][2],2))));
    auto subbigexp2 = deep_copy(exp(-c[b]*(pow(dx[0],2)+pow(dx[1],2)+pow(r[a][1]-r[b][2],2))));

    auto K6 = deep_copy(std::sqrt(pi)*epsilon*0.5*pow(c[b],-0.5)*(
                -(epsilon+2*c[b]*dx[0])*suberf3*subexp3
                -(epsilon-2*c[b]*dx[0])*suberf4*subexp4)
            *subbigexp1);

    auto suberf5 = deep_copy(erf(0.5*pow(c[b],-0.5)*(epsilon-2*c[b]*(r[a][1]-r[b][2])))-1);
    auto suberf6 = deep_copy(erf(0.5*pow(c[b],-0.5)*(epsilon+2*c[b]*(r[a][1]-r[b][2])))-1);
    auto subexp5 = deep_copy(exp(0.25*pow(epsilon-2*c[b]*(r[a][1]-r[b][2]),2)/c[b]));
    auto subexp6 = deep_copy(exp(0.25*pow(epsilon+2*c[b]*(r[a][1]-r[b][2]),2)/c[b]));

    auto K7 = deep_copy(std::sqrt(pi)*epsilon*sqrt(c[b])*
        (-dx[0]*(suberf5*subexp5-suberf6*subexp6))*subbigexp2);

    auto K8 = deep_copy(2*c[b]*(r[b][0]-r[a][1])*exp(-c[b]*pow(r[b][0]-r[a][1],2)));

    auto K9 = deep_copy(std::sqrt(pi)*epsilon*pow(c[b],-0.5)*0.5*
        (suberf5*subexp5-suberf6*subexp6)*subbigexp2);

    auto K10 = deep_copy(std::sqrt(pi)*epsilon*pow(c[b],-0.5)*0.5*
        (suberf3*subexp3-suberf4*subexp4)*subbigexp1);

    auto K11 = deep_copy(-2*c[b]*dx[0]*exp(-c[b]*(pow(dx[0],2)+pow(dx[1],2))));
    auto K12 = deep_copy(-2*c[b]*dx[1]*exp(-c[b]*(pow(dx[0],2)+pow(dx[1],2))));

    auto subdiff = deep_copy(c[b]*(pow(dx[0],2)+pow(dx[1],2)));
    auto K13 = deep_copy(4.0*c[b]*(subdiff-1)*exp(-subdiff));

    auto K14 = deep_copy(exp(-c[b]*(pow(dx[0],2)+pow(dx[1],2)+pow(dx[2],2))));
    auto K15 = deep_copy(exp(-c[b]*(pow(dx[0],2)+pow(r[b][1]-r[a][2],2))));
    auto K16 = deep_copy(exp(-c[b]*(pow(r[b][0]-r[a][1],2)+pow(r[b][1]-r[a][2],2))));
    auto K17 = deep_copy(exp(-c[b]*(pow(r[b][0]-r[a][2],2))));



    Eigen::Map<Eigen::Matrix<double,n,1> > s_vect(get<scalar>(particles).data());


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
        for (int i=0; i < max_iter; i++) {
            //solve for J(x_n) (x_{n+1} - x_n) = -F(x_n)
            //terminate if ||x_{n+1}-x_n|| < tol
        }
    }
}

