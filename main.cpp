#include<fstream>
#include<iostream>
#include<sstream>
#include<cmath>
#include<vector>

using namespace std;

/*
Main simulation file

If manual execution:
This script needs to executed after creating the Plots/rho Plots/star ans results directory. 
It handles the setup and calculation of the given density distribution 
and saves snapshots of the system in timestep intervals of savetime = 0.1.
*/

//parameters
double K = 0.1; //polytropic constant
double gam = 2.0; //polytropic exponant
double nu = 0.1; //dampening coeff
double lambda = 2.012032860; //linear acc
double h = 0.2;  //kernel size
double delta_t = 0.01; //fixed time step
double tend = 20.0; //simulation


typedef vector<double> Vec;
typedef vector<vector<double> > Mat;


/// @brief class that handles setup and interaction of particles with additional i/o capabilities
class particles{
    public:
    int num;
    Mat x;
    Mat v;
    Mat a;
    Vec rho;
    Vec m;
    Vec p;
    vector<vector<int> > nn;

    //setup and i/o:
    particles(string);
    particles(Vec, Vec);
    void print_particle(unsigned int);
    void print_to_file(int);

    //helper functions:
    double dist(int, int);
    vector<int> calc_nn(int);
    void calc_nn();
    double calc_W(double);
    Vec calc_grad_W(double, int, int);

    //calculation of physical quantities:
    void calc_rho();
    void calc_acc();

    //one function to rule them all:
    void run();

};


/// @brief a function to split a string seperated via delimiters into multiple substrings
/// @param str string to split
/// @param delim char on which to split
/// @param out pointer to vector<string> object, in which to store seperate substrings
void tokenize(std::string const &str, const char delim,
            std::vector<std::string> &out)
{
    size_t start;
    size_t end = 0;
 
    while ((start = str.find_first_not_of(delim, end)) != std::string::npos)
    {
        end = str.find(delim, start);
        out.push_back(str.substr(start, end - start));
    }
}


/// @brief constructor for particles class: initializes object with particle parameters specified in sourcefile
/// @param sourcefile name of txt file containing particle parameters in form x y z vx vy vz m
particles::particles(string sourcefile){
    ifstream file;
    file.open(sourcefile);
    string line;
    int n = 0;
    if (file.is_open()){
        while (getline(file, line)){
            Vec pos;
            Vec vel;
            vector<string> s;
            tokenize(line, ' ', s);
            pos.push_back(stod(s[0]));
            pos.push_back(stod(s[1]));
            pos.push_back(stod(s[2]));
            x.push_back(pos);
            vel.push_back(stod(s[3]));
            vel.push_back(stod(s[4]));
            vel.push_back(stod(s[5]));
            v.push_back(vel);
            m.push_back(stod(s[6]));
            n++;
        }
        num = n;
        a = Mat(num, Vec(3, 0.0));
        rho = Vec(num, 0.0);
        p = Vec(num, 0.0);
        nn = vector<vector<int>>(num);
    }

}


/// @brief constructor for 2 particles. Only for debugging.
/// @param a x, y, z, vx, vy, vz, m particle 1
/// @param b x, y, z, vx, vy, vz, m particle 2
particles::particles(Vec a, Vec b){
    h = 1.0;
    num = 2;
    this -> x = Mat(2, Vec(3, 0.0));
    this -> v = Mat(2, Vec(3, 0.0));
    this -> a = Mat(2, Vec(3, 0.0));
    this -> rho = Vec(2, 0.0);
    this -> m = Vec(2, 0.0);
    this -> p = Vec(2, 0.0);

    for (size_t i = 0; i != x[0].size(); ++i){
        x[0][i] = a[i];
        v[0][i] = a[i+3];
    }
    m[0] = a[-1];
    for (size_t i=0; i!=x[1].size(); ++i){
        x[1][i] = b[i];
        v[1][i] = b[i+3];
    }
    m[1] = b[-1];
}


/// @brief prints attributes of single particle
/// in form i  x y z  vx vy vz  ax ay az    rho p\n
/// @param i index of particle
void particles::print_particle(unsigned int i){
    cout << i << "  ";
    for (size_t j = 0; j != x[i].size(); j++){
        cout << x[i][j] << " ";
    }
    cout << " ";
    for (size_t j = 0; j != v[i].size(); ++j){
        cout << v[i][j] << " ";
    }
    cout << " ";
    for (size_t j = 0; j != a[i].size(); ++j){
        cout << a[i][j] << " ";
    }
    cout << "    " << rho[i] << " " << p[i] << endl;
}


/// @brief gives the distance between to particles
/// @param i index particle 1
/// @param j index particle 2
/// @return distance as double
inline double particles::dist(int i, int j){
    //bugchecked
    return sqrt(pow((x[i][0]-x[j][0]), 2.0) + pow((x[i][1]-x[j][1]), 2.0) + pow((x[i][2]-x[j][2]), 2.0));
}


/// @brief --!--OUTDATED DUE TO EFFICIENCY--!-- 
/// calculates the neighbors of a particle
/// @param i index of particle for which neighbors need to be found
/// @return vector<int> vector of indices of neighbors, includes self
vector<int> particles::calc_nn(int i){
    //slow version
    vector<int> neighbors;
    for (int j = 0; j < num; ++j){
        if (dist(i, j)/h < 1.0){
            neighbors.push_back(j);
        }
    }
    return neighbors;
}


/// @brief calculates the list of nearest neighbors for each particle
void particles::calc_nn(){
    for (int i = 0; i < num; ++i){
        //clear out prev:
        nn[i].clear();
        //add itself to list:
        nn[i].push_back(i);
        //using symmetry r_ij = r_ji,
        //therefore i is nn to j implies j is nn to i
        if (i != 0){
            for (int j = 0; j < i; ++j){
                if (dist(i, j) < h){
                    //adding i to j nn list and j to i nn list:
                    nn[i].push_back(j);
                    nn[j].push_back(i);
                }
            }
        }
    }
}


/// @brief calculates the cubic b-spline kernel function
/// @param r distance on which to evaluate W as double
/// @return W(r) as double
inline double particles::calc_W(double r){
    double temp = r/h;
    double ans = 0.0;
    if (temp<0.5){
        ans = 6.0*pow(temp, 3.0)-6.0*pow(temp, 2.0)+1.0;
    }
    else
    {
        ans = 2.0*pow((1.0-temp), 3.0);
    }
    return ans*8.0/(M_PI*pow(h, 3.0));
}


/// @brief calculates the gradient of the cubic b-spline kernel function
/// @param r distance between particles
/// @param i index of particle one
/// @param j index of particle two
/// @return gradient vector of W as vector<double>
inline Vec particles::calc_grad_W(double r, int i, int j){
    double dW = 0.0;
    double temp = r/h;
    //calculate del W/del r
    if (temp<0.5){
        dW = 3.0*pow(temp, 2.0)-2*temp;
    }
    else{
        dW = -pow((1-temp), 2.0);
    }
    dW = dW*6.0*8.0/(M_PI*pow(h, 4.0)); //Eq (10)
    // calculate grad W
    vector<double> gradW(3, 0.0);
    if (r != 0.0){
        for (int k = 0; k<3; ++k){
            gradW[k] = (x[i][k] - x[j][k])/r *dW;
        }
    }
    return gradW;
}


/// @brief calculates the density of each particle
void particles::calc_rho(){
    for (int i=0; i<num; ++i){
        double sum = 0.0;
        for (int j=0; j<nn[i].size(); ++j){
            sum += m[nn[i][j]]*calc_W(dist(i, nn[i][j]));
        }
        rho[i] = sum;
    }
}


/// @brief calculates the acceleration vector of each particle
void particles::calc_acc(){
    for (int i = 0; i<num; ++i){
        // first pressure acceleration:
        Vec sum = {0.0, 0.0, 0.0};
        for (int j=0; j<nn[i].size(); ++j){
            Vec grad = calc_grad_W(dist(i, nn[i][j]), i, nn[i][j]);
            double scalars = m[nn[i][j]]*(K*pow(rho[i], gam)/pow(rho[i], 2.0) + K*pow(rho[nn[i][j]], gam)/pow(rho[nn[i][j]], 2.0));
            sum[0] -= scalars * grad[0];
            sum[1] -= scalars * grad[1];
            sum[2] -= scalars * grad[2];
        }
        a[i] = sum;

        //linear acceleration and dampening force:
        for (int l = 0; l<3; ++l){
            a[i][l] -= (lambda * x[i][l] + nu*v[i][l]);
        }
    }
}


/// @brief Runs the simulation until time reaches tend and saves snapshots of system
void particles::run(){
    double  t = 0.0;
    int n_savestep = 0;
    double savetime = 0.0;
    Mat v_at_t(num, Vec(3));
    while (t<tend)
    {
        calc_nn();
        calc_rho();
        //print to file after rho is calculated
        if (t >= savetime){
            print_to_file(n_savestep);
            savetime += 0.1;
            n_savestep += 1;
        }
        calc_acc();
        //update x, v for half step:
        for (int i = 0; i < num; ++i){
            for (int l = 0; l<3; ++l){
                x[i][l] += delta_t/2.0*v[i][l]; //calc r(t+1/2)
                v_at_t[i][l] = v[i][l]; //preserve v(t)
                v[i][l] += delta_t/2.0*a[i][l]; //calc v(t+1/2)
            }
        }
        //calc all needed quantities for t+1/2:
        calc_nn();
        calc_rho();
        calc_acc();
        //update x, v to get x(t+1) and v(t+1)
        for (int i = 0; i < num; ++i){
            for (int l = 0; l<3; ++l){
                v[i][l] = v_at_t[i][l] + delta_t*a[i][l]; //calc v(t+1)
                x[i][l] += delta_t/2.0*v[i][l]; //calc r(t+1)    
            }
        }
        t += delta_t;
    }
    
}


/// @brief prints lines containing the x y z rho of all particles to out{num_save}.txt
/// @param num_save int variable to individualize outputfiles
void particles::print_to_file(int num_save){
    ofstream outfile;
    string filename = "results/out" + to_string(num_save) + ".txt";
    outfile.open(filename);
    for (int i=0; i<num; ++i){
        outfile << x[i][0] << ' ' << x[i][1] << ' ' << x[i][2] << ' ';
        //outfile << v[i][0] << ' ' << v[i][1] << ' ' << v[i][2] << endl;
        outfile << rho[i] << endl;
    }
    outfile.close();
}


/// @brief test function for particles class, 
void debug(){
    //particles test(x, y);
    particles test("random_distribution.dat");
    cout << "Particles: " << endl;
    test.print_particle(0);
    test.print_particle(1);
    
    //test nearest neighbors and acc
    cout << "Neighbor lists:" << endl;
    test.calc_nn();
    for (int k = 0; k<test.num; ++k){
        for (int i = 0; i<test.nn[k].size(); ++i){
            cout << test.nn[k][i] << endl;
        }
        cout << endl;
    }
    test.calc_rho();
    test.calc_acc();
    test.print_particle(0);
    test.print_particle(1);
    test.print_particle(2);

    //test integration:
    test.run();
    test.print_particle(0);
    test.print_particle(1);
    test.print_particle(2);

}


int main(){
    //debug();
    particles ToyStar("random_distribution.dat");//init instance of particle class
    ToyStar.run(); //run simulation
}