#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_VECTOR_SIZE 40
#include "Aboria.h"
using namespace Aboria;

#include <boost/math/constants/constants.hpp>
#include "boost/program_options.hpp" 
namespace po = boost::program_options;

const double PI = boost::math::constants::pi<double>();
enum linear_solver {CG, BiCGSTAB, GMRES};

template<typename Kernel,typename VectorType>
void solve(Kernel &&kernel, VectorType &&result, VectorType &&source, size_t max_iter=10, size_t restart=10, linear_solver solver=CG) {
    switch (solver) {
        case CG: {
            Eigen::ConjugateGradient<Kernel, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg;
            cg.setMaxIterations(max_iter);
            cg.compute(kernel);
            result = cg.solve(source);
            std::cout << "CG:    #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
            break;
                 }
        case BiCGSTAB: {
            Eigen::BiCGSTAB<Kernel, Eigen::DiagonalPreconditioner<double>> bicg;
            bicg.setMaxIterations(max_iter);
            bicg.compute(kernel);
            result = bicg.solve(source);
            std::cout << "BiCGSTAB:    #iterations: " << bicg.iterations() << ", estimated error: " << bicg.error() << std::endl;
            break;
               }
        case GMRES: {
            Eigen::GMRES<Kernel, Eigen::DiagonalPreconditioner<double>> gmres;
            gmres.set_restart(restart);
            gmres.setMaxIterations(max_iter);
            gmres.compute(kernel);
            result = gmres.solve(source);
            std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
            break;
                    }
    }

}

int main(int argc, char **argv) {

    unsigned int nout,max_iter_linear,restart_linear,nx,nr;
    double dt_aim,c0,factorr;
    unsigned int solver_in;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("max_iter_linear", po::value<unsigned int>(&max_iter_linear)->default_value(20), "maximum iterations for linear solve")
        ("restart_linear", po::value<unsigned int>(&restart_linear)->default_value(20), "iterations until restart for linear solve")
        ("linear_solver", po::value<unsigned int>(&solver_in)->default_value(0), "linear solver")
        ("nout", po::value<unsigned int>(&nout)->default_value(10), "number of output points")
        ("c0", po::value<double>(&c0)->default_value(1.0), "kernel constant")
        ("nx", po::value<unsigned int>(&nx)->default_value(10), "nx")
        ("nr", po::value<unsigned int>(&nr)->default_value(10), "nr")
        ("factorr", po::value<double>(&factorr)->default_value(1.5), "factorr")
        ("dt", po::value<double>(&dt_aim)->default_value(0.0001), "timestep")
    ;
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);  

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    

    ABORIA_VARIABLE(density1p,double,"density")
    ABORIA_VARIABLE(density2p,double,"two particle density")
    ABORIA_VARIABLE(g_function,double,"g function")

    ABORIA_VARIABLE(density1p_weights,double,"density weights")
    ABORIA_VARIABLE(density2p_weights,double,"two particle density weights")
    ABORIA_VARIABLE(g_function_weights,double,"g function weights")

    ABORIA_VARIABLE(density1p_weights0,double,"n-1 density weights")
    ABORIA_VARIABLE(density2p_weights0,double,"n-1 two particle density weights")

    ABORIA_VARIABLE(FP_vector,double,"FP_vector")
    ABORIA_VARIABLE(Fp_vector,double,"Fp_vector")
    ABORIA_VARIABLE(Fg_vector,double,"Fg_vector")

    ABORIA_VARIABLE(constant,double,"constant for basis functions")

    ABORIA_VARIABLE(K1wP_vector,double,"K1wP")
    ABORIA_VARIABLE(K2wp_vector,double,"K2wp")
    ABORIA_VARIABLE(K3wP_vector,double,"K3wP")
    ABORIA_VARIABLE(K4wp_vector,double,"K4wp")
    ABORIA_VARIABLE(K5wp_vector,double,"K5wp")
    ABORIA_VARIABLE(K6wg_vector,double,"K6wg")
    ABORIA_VARIABLE(K7wg_vector,double,"K7wg")
    ABORIA_VARIABLE(K8wp_vector,double,"K8wp")
    ABORIA_VARIABLE(K9wg_vector,double,"K9wg")
    ABORIA_VARIABLE(K10wg_vector,double,"K10wg")
    ABORIA_VARIABLE(K11wP_vector,double,"K11wP")
    ABORIA_VARIABLE(K12wP_vector,double,"K12wp")
    ABORIA_VARIABLE(K13wP_vector,double,"K13wP")
    ABORIA_VARIABLE(K14wg_vector,double,"K14wg")
    ABORIA_VARIABLE(K15wP_vector,double,"K15wP")
    ABORIA_VARIABLE(K16wP_vector,double,"K16wP")
    ABORIA_VARIABLE(K17wp_vector,double,"K17wP")
    ABORIA_VARIABLE(K18wp_vector,double,"K18wp")

    typedef Particles<std::tuple<
            density1p,density2p,g_function,
            density1p_weights,density2p_weights,g_function_weights,
            density1p_weights0,density2p_weights0,
            FP_vector, Fp_vector, Fg_vector,
            constant,
            K1wP_vector,K2wp_vector,K3wP_vector,K4wp_vector,K5wp_vector,
            K6wg_vector,K7wg_vector,K8wp_vector,K9wg_vector,K10wg_vector,
            K11wP_vector,K12wP_vector,K13wP_vector,K14wg_vector,
            K15wP_vector,K16wP_vector,K17wp_vector,K18wp_vector
            >,3> ParticlesType;
    typedef position_d<3> position;
    ParticlesType knots;

    const double L = 1.0;
    const double3 low(0,0,0);
    const double3 high(L,L,L);
    const int max_iter_newton = 10;
    const double tol2_newton = std::pow(1e-6,2);
    const double tol2_gmres = std::pow(1e-6,2);
    double2 periodic(false);
    const double epsilon = 0.05;
    const double theta = 0.2;
    const double Tf = 0.02;
    const int timesteps = Tf/dt_aim;
    const double dt = Tf/timesteps;
    const int N = 16;
    const double id_theta = 0.2;
    const double id_alpha = 30.0;

    
    int n = nx*nr*nx;
    const double deltax = std::sqrt(2.0)/nx;
    typename ParticlesType::value_type particle;
    for (int i=0; i<=nx; ++i) {
        const double x = i*deltax;
        double deltar = L*(1-factorr)/(1-std::pow(factorr,nr+1));
        double r = 0;
        for (int k=0; k<=nr; ++k) {
            r += deltar;
            deltar *= factorr;
            double deltatheta = deltar/r;
            const double ntheta = 2*PI/deltatheta;
            const double start_theta = deltatheta*(k%2)/2.0;
            for (int j=0; j<=ntheta; ++j) {
                const double theta = j*deltatheta + start_theta;
                //std::cout << "x, r, theta = "<<x<<","<<r<<","<<theta<<std::endl;
//get<position>(particle) = double3((x-r*cos(theta))/std::sqrt(2.0),(x+r*cos(theta))/std::sqrt(2.0),r*std::sin(theta));
                get<position>(particle) = double3(
                            x/std::sqrt(3) - std::sqrt(6)*r*sin(theta+PI/3.0)/3.0,
                            x/std::sqrt(3) + std::sqrt(6)*r*cos(theta+PI/3.0)/3.0,
                            x/std::sqrt(3) + std::sqrt(6)*r*sin(theta)/3.0
                            );
                if ((get<position>(particle)<low).any() || (get<position>(particle)>=high).any()) continue;
                get<constant>(particle) = c0;
                get<density1p>(particle) = 0.5*(tanh(id_alpha*(get<position>(particle)[0]-id_theta))+tanh(id_alpha*(1-id_theta-get<position>(particle)[0])));
                knots.push_back(particle);
            }
        }
    }
    std::cout << "added "<<knots.size()<<" knots" << std::endl;
    std::cout << "expected "<<n<<" knots" << std::endl;
    knots.init_neighbour_search(low,high,L/10,bool3(true,true,true));
    n = knots.size();

#ifdef HAVE_VTK
    vtkWriteGrid("knots_init",0,knots.get_grid(true));
#endif
 

    Symbol<position> r;
    Symbol<density1p> p;
    Symbol<density2p> P;
    Symbol<g_function> g;

    Symbol<density1p_weights> wp;
    Symbol<density2p_weights> wP;
    Symbol<g_function_weights> wg;

    Symbol<density1p_weights0> wp0;
    Symbol<density2p_weights0> wP0;

    Symbol<FP_vector> FP;
    Symbol<Fp_vector> Fp;
    Symbol<Fg_vector> Fg;

    Symbol<K1wP_vector> K1wP;
    Symbol<K2wp_vector> K2wp;
    Symbol<K3wP_vector> K3wP;
    Symbol<K4wp_vector> K4wp;
    Symbol<K5wp_vector> K5wp;
    Symbol<K6wg_vector> K6wg;
    Symbol<K7wg_vector> K7wg;
    Symbol<K8wp_vector> K8wp;
    Symbol<K9wg_vector> K9wg;
    Symbol<K10wg_vector> K10wg;
    Symbol<K11wP_vector> K11wP;
    Symbol<K12wP_vector> K12wP;
    Symbol<K13wP_vector> K13wP;
    Symbol<K14wg_vector> K14wg;
    Symbol<K15wP_vector> K15wP;
    Symbol<K16wP_vector> K16wP;
    Symbol<K17wp_vector> K17wp;
    Symbol<K18wp_vector> K18wp;

    Symbol<constant> c;
    Label<0,ParticlesType> a(knots);
    Label<1,ParticlesType> b(knots);
    auto dx = create_dx(a,b);
    Accumulate<std::plus<double> > sum;

    auto k = deep_copy(
            exp(-dot(dx,dx)*c[b])
            );
    auto f = deep_copy(
            if_else(r[a][0] > r[a][1],
                -epsilon*exp(-epsilon*(r[a][0]-r[a][1])),
                epsilon*exp(-epsilon*(r[a][1]-r[a][0]))
                )
            );

    auto suberf1 = deep_copy(erfc(0.5*pow(c[b],-0.5)*(epsilon-2*c[b]*(r[a][0]-r[b][1]))));
    auto suberf2 = deep_copy(erfc(0.5*pow(c[b],-0.5)*(epsilon+2*c[b]*(r[a][0]-r[b][1]))));
    auto subexp1 = deep_copy(exp(0.25*pow(epsilon-2*c[b]*(r[a][0]-r[b][1]),2)/c[b]));
    auto subexp2 = deep_copy(exp(0.25*pow(epsilon+2*c[b]*(r[a][0]-r[b][1]),2)/c[b]));

    auto K1 = deep_copy(std::sqrt(PI)*epsilon*pow(c[b],-0.5)*0.5*(
                (epsilon+2*c[b]*dx[0])*suberf1*subexp1
                +(epsilon-2*c[b]*dx[0])*suberf2*subexp2)
            *exp(-c[b]*(pow(dx[0],2)+pow(r[a][0]-r[b][1],2))));

    auto K2 = deep_copy(2*c[b]*(2*c[b]*pow(dx[0],2)-1)*exp(-c[b]*pow(dx[0],2)));
    auto K3 = deep_copy(exp(-c[b]*(pow(dx[0],2)+pow(dx[1],2))));
    auto K4 = deep_copy(exp(-c[b]*(pow(dx[0],2))));
    auto K5 = deep_copy(exp(-c[b]*(pow(r[b][0]-r[a][1],2))));
    
    auto suberf3 = deep_copy(erfc(0.5*pow(c[b],-0.5)*(epsilon-2*c[b]*(r[a][0]-r[b][2]))));
    auto suberf4 = deep_copy(erfc(0.5*pow(c[b],-0.5)*(epsilon+2*c[b]*(r[a][0]-r[b][2]))));
    auto subexp3 = deep_copy(exp(0.25*pow(epsilon-2*c[b]*(r[a][0]-r[b][2]),2)/c[b]));
    auto subexp4 = deep_copy(exp(0.25*pow(epsilon+2*c[b]*(r[a][0]-r[b][2]),2)/c[b]));
    auto subbigexp1 = deep_copy(exp(-c[b]*(pow(dx[0],2)+pow(dx[1],2)+pow(r[a][0]-r[b][2],2))));
    auto subbigexp2 = deep_copy(exp(-c[b]*(pow(dx[0],2)+pow(dx[1],2)+pow(r[a][1]-r[b][2],2))));

    auto K6 = deep_copy(std::sqrt(PI)*epsilon*0.5*pow(c[b],-0.5)*(
                (epsilon+2*c[b]*dx[0])*suberf3*subexp3
                +(epsilon-2*c[b]*dx[0])*suberf4*subexp4)
            *subbigexp1);

    auto suberf5 = deep_copy(erfc(0.5*pow(c[b],-0.5)*(epsilon-2*c[b]*(r[a][1]-r[b][2]))));
    auto suberf6 = deep_copy(erfc(0.5*pow(c[b],-0.5)*(epsilon+2*c[b]*(r[a][1]-r[b][2]))));
    auto subexp5 = deep_copy(exp(0.25*pow(epsilon-2*c[b]*(r[a][1]-r[b][2]),2)/c[b]));
    auto subexp6 = deep_copy(exp(0.25*pow(epsilon+2*c[b]*(r[a][1]-r[b][2]),2)/c[b]));

    auto K7 = deep_copy(std::sqrt(PI)*epsilon*sqrt(c[b])*
        (dx[0]*(suberf5*subexp5-suberf6*subexp6))*subbigexp2);

    auto K8 = deep_copy(2*c[b]*(r[b][0]-r[a][1])*exp(-c[b]*pow(r[b][0]-r[a][1],2)));
    auto K18 = deep_copy(-2*c[b]*(r[a][0]-r[b][0])*exp(-c[b]*pow(r[a][0]-r[b][0],2)));

    auto K9 = deep_copy(std::sqrt(PI)*epsilon*pow(c[b],-0.5)*0.5*
        (suberf6*subexp6-suberf5*subexp5)*subbigexp2);

    auto K10 = deep_copy(std::sqrt(PI)*epsilon*pow(c[b],-0.5)*0.5*
        (suberf4*subexp4-suberf3*subexp3)*subbigexp1);

    auto K11 = deep_copy(-2*c[b]*dx[0]*exp(-c[b]*(pow(dx[0],2)+pow(dx[1],2))));
    auto K12 = deep_copy(-2*c[b]*dx[1]*exp(-c[b]*(pow(dx[0],2)+pow(dx[1],2))));

    auto subdiff = deep_copy(c[b]*(pow(dx[0],2)+pow(dx[1],2)));
    auto K13 = deep_copy(4.0*c[b]*(subdiff-1)*exp(epsilon*abs(r[a][0]-r[a][1])+subdiff)*
                                              exp(-epsilon*abs(r[a][0]-r[a][1])-subdiff));

    auto K14 = deep_copy(exp(-c[b]*(pow(dx[0],2)+pow(dx[1],2)+pow(dx[2],2))));
    auto K15 = deep_copy(exp(-c[b]*(pow(dx[0],2)+pow(r[b][1]-r[a][2],2))));
    auto K16 = deep_copy(exp(-c[b]*(pow(r[b][0]-r[a][1],2)+pow(r[b][1]-r[a][2],2))));
    auto K17 = deep_copy(exp(-c[b]*(pow(r[b][0]-r[a][2],2))));


    auto fp = deep_copy((N-1)*K1wP[a] + K2wp[a]);
    auto fP = deep_copy((N-2)*(
                (K3wP[a]*(K6wg[a] + K7wg[a] 
                          - K8wp[a]*K9wg[a]/K5wp[a]
                          - K18wp[a]*K10wg[a]/K4wp[a])
                + K11wP[a]*K10wg[a]
                + K12wP[a]*K9wg[a])/(K4wp[a]*K5wp[a])
                + K13wP[a]));
    auto fg = deep_copy(K15wP[a]*K16wP[a]/K17wp[a]);


    typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_type;
    typedef Eigen::Map<vector_type> map_type;

    

    solve(create_eigen_operator(a,b,K4),
            map_type(get<density1p_weights>(knots).data(),n),
            map_type(get<density1p>(knots).data(),n),max_iter_linear,restart_linear,(linear_solver)solver_in);
    
    //estimate P and solve for weights 
    P[a] = sum(b,true,K4*wp[b])*sum(b,true,K5*wp[b]);
    solve(create_eigen_operator(a,b,K3),
            map_type(get<density2p_weights>(knots).data(),n),
            map_type(get<density2p>(knots).data(),n),max_iter_linear,restart_linear,(linear_solver)solver_in);
    
    //estimate g and solve for weights 
    g[a] = sum(b,true,K15*wP[b])*sum(b,true,K16*wP[b])/sum(b,true,K17*wp[b]);
    solve(create_eigen_operator(a,b,K14),
            map_type(get<g_function_weights>(knots).data(),n),
            map_type(get<g_function>(knots).data(),n),max_iter_linear,restart_linear,(linear_solver)solver_in);

    //init last iteration weights
    wp0[a] = wp[a];
    wP0[a] = wP[a];


#ifdef HAVE_VTK
    vtkWriteGrid("knots_init",1,knots.get_grid(true));
#endif

    vector_type source(3*n);
    vector_type result(3*n);
    for (int i=0; i < timesteps; i++) {
        //evaluate temporaries
        K3wP[a] = sum(b,true,K3*wP[b]);
        K4wp[a] = sum(b,true,K4*wp[b]);
        K5wp[a] = sum(b,true,K5*wp[b]);
        K6wg[a] = sum(b,true,K6*wg[b]);
        K7wg[a] = sum(b,true,K7*wg[b]);
        K8wp[a] = sum(b,true,K8*wp[b]);
        K9wg[a] = sum(b,true,K9*wg[b]);
        K10wg[a] = sum(b,true,K10*wg[b]);
        K11wP[a] = sum(b,true,K11*wP[b]);
        K12wP[a] = sum(b,true,K12*wP[b]);
        K13wP[a] = sum(b,true,K13*wP[b]);
        K14wg[a] = sum(b,true,K14*wg[b]);
        K18wp[a] = sum(b,true,K18*wp[b]);

        //explicit euler step
        P[a] = dt*fP + K3wP[a];
        p[a] = dt*fp + K4wp[a];
        solve(create_eigen_operator(a,b,K4),
            map_type(get<density1p_weights>(knots).data(),n),
            map_type(get<density1p>(knots).data(),n),max_iter_linear,restart_linear,(linear_solver)solver_in);
        solve(create_eigen_operator(a,b,K3),
            map_type(get<density2p_weights>(knots).data(),n),
            map_type(get<density2p>(knots).data(),n),max_iter_linear,restart_linear,(linear_solver)solver_in);
        
        K15wP[a] = sum(b,true,K15*wP[b]);
        K16wP[a] = sum(b,true,K16*wP[b]);
        K17wp[a] = sum(b,true,K17*wp[b]);
        g[a] = fg;
        solve(create_eigen_operator(a,b,K14),
            map_type(get<g_function_weights>(knots).data(),n),
            map_type(get<g_function>(knots).data(),n),max_iter_linear,restart_linear,(linear_solver)solver_in);

#ifdef HAVE_VTK
        vtkWriteGrid("knots_explicit",i,knots.get_grid(true));
#endif
        
        /*
        int ii=0
        for (; ii < max_iter_newton; ii++) {
            
            //evaluate temporaries
            K3wP[a] = sum(b,K3*wP[b]);
            K4wp[a] = sum(b,K4*wp[b]);
            ...

            FP[a] = dt*fP - sum(b,K3*(wp[b]-wp0[b]));
            Fp[a] = dt*fp - sum(b,K4*(wp[b]-wp0[b]));
            Fg[a] = fg - sum(b,true,K14*wg[b]);

            //copy F* to b
            source.segment(0,n) = -map_type(get<FP_vector>(knots).data(),n);
            source.segment(n,2*n) = -map_type(get<Fp_vector>(knots).data(),n);
            source.segment(2*n,3*n) = -map_type(get<Fg_vector>(knots).data(),n);

            //create J operator
            auto JPP = create_eigen_operator(a,b,dt*dfPdP - K3);
            auto JPp = create_eigen_operator(a,b,dt*dfPdp);
            auto JPg = create_eigen_operator(a,b,dt*dfPdg);
            auto JpP = create_eigen_operator(a,b,dt*dfpdP);
            auto Jpp = create_eigen_operator(a,b,dt*dfpdp - K4);
            auto Jpg = create_eigen_operator(a,b,dt*dfpdg);
            auto JgP = create_eigen_operator(a,b,dt*dfgdP);
            auto Jgp = create_eigen_operator(a,b,dt*dfgdp);
            auto Jgg = create_eigen_operator(a,b,dt*dfgdg);

            //solve for J(x_n) (x_{n+1} - x_n) = -F(x_n)
            solve(create_block_eigen_operator<3,3>(
                        JPP,JPp,JPg,
                        JpP,Jpp,Jpg,
                        JgP,Jgp,Jgg,
                        )
                    result,source);

            
            //increment weights by result
            map_type(get<density2p_weights>(knots).data(),n) += result.segment(0,n);
            map_type(get<density1p_weights>(knots).data(),n) += result.segment(n,2*n);
            map_type(get<g_function_weights>(knots).data(),n) += result.segment(2*n,3*n);
                        
            //terminate if ||x_{n+1}-x_n|| < tol
            if (result.squaredNorm() < tol2_gmres) {
                break;
            }
        }

        std::cout << "finished newton iterations, ii = "<<ii<<"/"<<max_iter_newton<<std::endl;
        */

    }
}

