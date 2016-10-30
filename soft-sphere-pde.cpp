#include <math.h>




int main(void) {
    ABORIA_VARIABLE(density1p,double,"density")
    ABORIA_VARIABLE(density2p,double,"two particle density")
    ABORIA_VARIABLE(g_function,double,"g function")

    ABORIA_VARIABLE(density1p_weights,double,"density weights")
    ABORIA_VARIABLE(density2p_weights,double,"two particle density weights")
    ABORIA_VARIABLE(g_function_weights,double,"g function weights")

    ABORIA_VARIABLE(density1p_weights0,double,"n-1 density weights")
    ABORIA_VARIABLE(density2p_weights0,double,"n-1 two particle density weights")
    ABORIA_VARIABLE(g_function_weights0,double,"n-1 g function weights")

    ABORIA_VARIABLE(constant,double,"constant for basis functions")

    typedef Particles<std::tuple<density1p,density2p,g_function,density1p_weights,density2p_weights,g_function_weights,,density1p_weights,density2p_weights,g_function_weights,constant>,3> ParticlesType;
    typedef position_d<3> position;
    ParticlesType knots;

    const double c = 0.5;
    const double L = 1.0;
    const double3 low(0,0,0);
    const double3 high(L,L,L);
    const int max_iter_gmres = 100;
    const int max_iter_newton = 100;
    const int restart_gmres = 100;
    double2 periodic(false);
    const double epsilon = 0.05;
    const double alpha = 30.0;
    const double theta = 0.2;
    const double Tf = 0.02;
    const double dt_aim = 0.001;
    const int timesteps = Tf/dt;
    const double dt = Tf/timesteps;
    const int N = 16;
    
    const int nx = 7;
    const int ntheta = 5;
    const int nr = 5;
    const int n = nx*ntheta*nr;
    const double deltar = 0.1;
    const double factorr = 1.01;
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
                get<density1p>(p) = 0.5*(tanh(alpha*(x-theta))+tanh(alpha*(1-theta-x)));
                knots.push_back(p);
            }
        }
    }
    std::cout << "added "<<knots.size()<<" knots" << std::endl;
    std::cout << "expected "<<n<<" knots" << std::endl;
    knots.init_neighbour_search(low,high,std::pow(c,2),bool3(true,true,true));

    Symbol<position> r;
    Symbol<density1p> p;
    Symbol<density2p> P;
    Symbol<g_function> g;

    Symbol<density1p_weights> wp;
    Symbol<density2p_weights> wP;
    Symbol<g_function_weights> wg;

    Symbol<density1p_weights0> wp0;
    Symbol<density2p_weights0> wP0;
    Symbol<g_function_weights0> wg0;

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

    auto fp = deep_copy((N-1)*sum(b,true,K1*wP[b]) + sum(b,true,K2*wp[b]));
    auto fP = deep_copy((N-2)*(
                (sum(b,true,K3*wP[b])*(sum(b,true,K6*wg[b]) 
                                      + sum(b,true,K7*wg[b]) 
                                      - sum(b,true,K8*wp[b])*sum(b,true,K9*wg[b])/sum(b,true,K5*wp[b])
                                      - sum(b,true,K11*wp[b])*sum(b,true,K10*wg[b])/sum(b,true,K4*wp[b]))
                + sum(b,true,K11*wP[b])*sum(b,true,K10*wg[b])
                + sum(b,true,K12*wP[b])*sum(b,true,K9*wg[b]))/(sum(b,true,K4*wp[b])*sum(b,true,K5*wp[b]))
                + sum(b,true,K13*wP[b])));
    auto fg = deep_copy(sum(b,true,K15*wP[b])*sum(b,true,K16*wP[b])/sum(b,true,K17*wp[b]) - sum(b,true,K14*wg[b]));


    auto F1 = 
    auto J = create_block_eigen_operator<3,3>(A,B,
                                                     C,Zero);
    typedef Eigen::Matrix<double,n,1> vector_type;

    template<typename Source, typename Result,typename Kernel>
    solve(Kernel &kernel, vector_type &result, vector_type &source) {
        Eigen::GMRES<decltype(W), Eigen::DiagonalPreconditioner<double>> gmres;
        gmres.set_restart(restart_gmres);
        gmres.setMaxIterations(max_iter_gmres);
        gmres.compute(create_eigen_operator(a,b,kernel));
        result = gmres.solve(source);
        std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
    }

    vector_type w0; 
    vector_type dw; 

    //solve for w0 weights from initial condition
    //
    //solve for p weights
    solve<density1p,density1p_weights>(K4);
    
    //estimate P and g
    P[a] = sum(b,true,K4*wp)*sum(b,true,K5*wp);
    solve<density2p,density2p_weights>(K3);

    g[a] = sum(b,true,K15*wP)*sum(b,true,K16*wP)/sum(b,true,K17*wP);
    solve<g_function,g_function_weights>(K14);

    wp0[a] = wp[a];
    wP0[a] = wP[a];
    wg0[a] = wg[a];

    for (int i=0; i < timesteps; i++) {
        for (int i=0; i < max_iter_newton; i++) {
            //solve for J(x_n) (x_{n+1} - x_n) = -F(x_n)
            solve(
            
            //terminate if ||x_{n+1}-x_n|| < tol
        }
    }
}

