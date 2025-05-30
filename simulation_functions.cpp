/*

OPERATING SYSTEM
NAME="Ubuntu"
VERSION="22.04.1 LTS (Jammy Jellyfish)"

COMPILER VERSION
g++ version 11.4.0 

COMPILATION COMMAND 
g++ simulation_functions.cpp -o simulation_functions

COMMAND LINE ARGUMENTS:
input1 (double) = relative mass of the sun (solar mass)
input2 (double) = relative mass of moon (lunar mass)
input3 (string) = name of csv file to be saved

EXAMPLE COMMAND
./assessment_6 1 1 orbit_data.csv

*/



#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <omp.h>
#include <string>
#include <fstream>
using namespace std;

template <typename T>
void print_vector(T vector);
double calc_mod(double x, double y);
vector<double> divide_vector_by_double(vector<double> v, double d);
vector<double> multiply_vector_by_double(vector<double> v, double d);
vector<double> add_vectors(vector<double> v1, vector<double> v2);
vector<vector<double>> find_other_positions(vector<vector<double>> all_current_positions, int index);
void print_len_sim(int max_time_in_hours);

class celestial_body {
    // The following variables and functions are private, 
    // This means they cannot be accessed from within the main script - only from functions within this class.
    private:
        vector<double> previous_acceleration = {0,0};
        vector<double> mod_r = {0,0};
        vector<double> current_acceleration = {0,0};

        double calc_single_force(double r_component, double single_mod_r, double other_mass) 
        /*
        Function to calculate the force component in 1d acting on a body from another body.

        inputs:
        r_component: component of vector in required dimension
        single_mod_r: scalar distance between the two bodies. 
        other_mass: mass of body causing force

        output:
        F_component: component of force in required dimension
        */
        {
            double G = 6.6743E-11;
            // double r_hat_component;
            double F_component;
            // r_hat_component = r_component / single_mod_r; 
            F_component = G*this->mass*other_mass*r_component / pow(single_mod_r,3);
            return F_component;
        }

        double calc_mod(int body_num)
        /*
        calculate distance between two points

        inputs:
        x: distance between x values of the two points
        y: distance between y values of the two points

        output:
        modulus: distance between two points
        */
        {
            double sum = 0;
            for (int i=0; i<num_dims;i++)
            {
                double component = this->other_positions[body_num][i] - this->current_position[i];
                sum += pow(component, 2);
            }
            // maybe use std::accumulate for more dimensions
            double modulus = sqrt( sum );
            return modulus;
        }


        void calc_mod_r()
        /*
        calculate vector of distances between body of interest and all other bodies.
        Automatically updates new mod_r of celestial body when calculated. 
        */
        {
            // vector<double> current_position = all_position[index];
            vector<double> mod_rs;
            double single_mod_r;
            for (int i=0; i<num_other_bods; i++)
            {
                single_mod_r = this->calc_mod(i);
                mod_rs.push_back(single_mod_r);
            }
            this->mod_r = mod_rs;
        }

        vector<double> calc_total_force()
        /*
        Function to calculate the total force component acting on a body.

        output:
        force_vector: vector of the components of the total force acting on a body. 
        */
        {
            int num_other = this->num_other_bods;
            vector<double> force_vector(this->num_dims);
            this->calc_mod_r();


            // loop through each number of dimensions
            for (int i=0; i<this->num_dims;i++)
            {
                // Loop through number of other bodies
                double total_force=0;
                for (int j=0; j<num_other;j++)
                {
                    double other_mass = this->other_masses[j];
                    double r_component = this->other_positions[j][i] - this->current_position[i];
                    double single_force = this->calc_single_force(r_component, this->mod_r[j], other_mass);           
                    total_force += single_force;
                }
                force_vector[i] = total_force;
            }
            return force_vector;
        }

        
    public:
        string name;
        int num_dims;
        ///////////////////////////// change num_bodies
        int num_other_bods;
        
        // mass
        double mass; // once set doesnt change
        vector<double> other_masses; // 1xnum_bodies

        // positions
        vector<vector<double>> all_positions; // vector of all positions in each dimensions ((xxx...)(yyy...))
        vector<double> current_position; // (xy)
        vector<vector<double>> other_positions; // (xy) (xy)
        
        // velocity
        vector<double> current_velocity; // (dx dy)

        double timestep;

        celestial_body(string name, int num_dims, int num_other_bods, double mass, vector<double> other_masses, 
                        vector<vector<double>> all_positions, vector<double> current_position,
                        vector<vector<double>> other_positions, vector<double> current_velocity, double timestep)
        {
            this->name = name;
            this->num_dims = num_dims;
            this->num_other_bods = num_other_bods;
            this->mass = mass;
            this->other_masses = other_masses; 
            this->all_positions = all_positions;
            this->current_position = current_position;
            this->other_positions = other_positions;
            this->current_velocity = current_velocity;
            this->timestep = timestep;
        }

        string describe() 
        {
            return this->name + " is a planet with mass " + to_string(this->mass) ;
        }

        
        void calc_acceleration()
        /*
        Calculate acceleration of a body.

        Automatically updates class acceleration values  
        */
        {
            // Update old acceleration
            this->previous_acceleration = this->current_acceleration;
            // Calculate new force
            vector<double> force_vector = this->calc_total_force();
            // Calculate acceleration
            vector<double> acceleration = divide_vector_by_double(force_vector, this->mass);
            // Update acceleration
            this->current_acceleration = acceleration;
        }

        void update_position()
        /*
        Calculate acceleration of a body.

        Automatically updates class position values  
        Also appends new current_position value to all_positions
        */
        {    
            // get each component of position sum
            vector<double> position_component = this->current_position;
            vector<double> velocity_component = multiply_vector_by_double(this->current_velocity,this->timestep);
            vector<double> acceleration_component = multiply_vector_by_double(this->current_acceleration,0.5*pow(this->timestep, 2));
            // Add each component of sum
            vector<double> pv_component = add_vectors(position_component, velocity_component);
            vector<double> new_position = add_vectors(pv_component, acceleration_component);
            // Update current position
            this->current_position = new_position;

            // Update all_positions for each dimension. Each dimension has to be updated separately. 
            for (int i=0; i<this->num_dims;i++)
            {
                this->all_positions[i].push_back(new_position[i]);
            }
        }

        void update_velocity()
        /*
        Calculate velocity of a body.

        Automatically updates class position values  
        Also appends new current value to 
        
        */
        {
            // get each component of position sum
            vector<double> velocity_component = this->current_velocity;
            vector<double> acceleration_component = multiply_vector_by_double(add_vectors(this->current_acceleration, this->previous_acceleration), (0.5*this->timestep));
            // Add each component of sum
            vector<double> new_velocity = add_vectors(velocity_component, acceleration_component);
            // Update velocity
            this->current_velocity = new_velocity;
        }
};

// These functions need the celestial_body class to be defined first. 
void write_to_file(string filename, int num_iters, int num_bodies, int num_dims, vector<celestial_body> celestial_bodies);
vector<vector<double>> update_all_positions(vector<celestial_body> celestial_bodies, int num_dims);


int main(int argc, char *argv[] )
{
    int num_inputs = 3;

    if (argc!=num_inputs+1)
    {
        throw invalid_argument(to_string(num_inputs) + " inputs expected");
        return 1;
    }

    // atof converts string to double
    double relative_sun_mass = atof(argv[1]);
    double relative_moon_mass = atof(argv[2]);
    string csv_filename = argv[3];


    ///// SETUP OF INITIAL CONDITIONS

    // State number of dimensions of simulation and number of bodies in simulation
    int num_dims=2;
    int num_bodies=3;

    // Set inital masses
    cout << "Relative mass of the earth: " << 1 << endl;
    double earth_mass = 5.9772E24;
    cout << "Relative mass of the sun: " << relative_sun_mass << endl;
    double sun_mass = 1.98847E30*relative_sun_mass; 
    cout << "Relative mass of the moon: " <<  relative_moon_mass << endl;
    double moon_mass = 7.3476E22*relative_moon_mass;

    vector<double> all_masses =  {sun_mass, earth_mass, moon_mass};

    // Initial conditions for the sun
    // Position
    // vector of vectors because all the x and y positions are going to be stored: ((x0, x1, x2...), (y0, y1, y1...))
    vector<vector<double>> sun_positions = {{0}, {0}};
    // Velocity
    vector<double> sun_vel_init = {0, 0};

    // Initial conditions for the earth
    // Position
    vector<vector<double>> earth_positions = {{ 1.495E11}, {0.0}};
    // Velocity
    vector<double> earth_vel_init = {0,  2.978E4};

    // Initial conditions for the moon
    // Position
    vector<vector<double>> moon_positions = {{ 1.495E11}, {3.84E8}};
    // Velocity
    vector<double> moon_vel_init = {-1.022E3, 2.978E4};    


    // Set timestep and maximum time to run simulation
    double timestep_in_hours = 1;
    // Convert to seconds
    double timestep = timestep_in_hours * 60 * 60;
    double max_time_in_years = 3;
    // Convert to hours
    double max_time_in_hours = max_time_in_years * 365 * 24 * (1/timestep_in_hours);

    // Vector of celestial_body classes
    vector<celestial_body> celestial_bodies;  
    vector<string> names = {"Sun", "Earth", "Moon"};

    // Vector of all x and y positions for each body
    vector<vector<vector<double>>> all_body_positions = {sun_positions, earth_positions, moon_positions};  // ((x)(y)) ((x)(y)) ((x)(y))

    // Initial velocities of each body
    vector<vector<double>> all_init_vel = {sun_vel_init, earth_vel_init, moon_vel_init}; // (xy) (xy) (xy)

    // Initialise each celestial body class
    for (int i=0; i<num_bodies; i++)  
    {
        // Determine initial position vector of body 
        vector<double> initial_position;
        for (int j=0; j<num_dims; j++)
        {
            double initial_position_component = all_body_positions[i][j][0];
            initial_position.push_back(initial_position_component);
        }

        // Get vector of position vectors of other bodies.  
        vector<vector<double>> other_initial_positions;
        for (int k=0;k<num_bodies;k++)
        {
            // Only include other bodies
            if (k!=i)
            {
                // position vector of one other body 
                vector<double> body_positions;
                for (int j=0;j<num_dims;j++)
                {
                    body_positions.push_back(all_body_positions[k][j][0]);
                }
                other_initial_positions.push_back(body_positions);
            }
        }

        // Calculate vector of masses of other bodies by removing mass of current body 
        vector<double> other_masses = all_masses;
        other_masses.erase(other_masses.begin()+i);

        // fill in celestial body class
        auto single_celestial_body = celestial_body(names[i], num_dims, num_bodies-1, all_masses[i], other_masses, all_body_positions[i],
                        initial_position, other_initial_positions, all_init_vel[i], timestep);

        // Add celestial body to vector of celestial bodies
        celestial_bodies.push_back(single_celestial_body);
    }           
    ////// END OF SETUP OF INITIAL CONDITIONS


    // Calc initial acceleration for each body 
    omp_set_num_threads(num_bodies);
    # pragma omp parallel
    {
        # pragma omp for
        for (int i=0; i<num_bodies; i++)
        {
            celestial_bodies[i].calc_acceleration();
        }
    }
    

    // Loop through each iteration
    for (int step=0; step<max_time_in_hours; step++)
    {
        // calculate new position of each body
        # pragma omp parallel
        {
            # pragma omp for
            for (int i=0; i<num_bodies; i++)
            {
                celestial_bodies[i].update_position();
            }
        }

        // Create vector of position vectors of all bodies 
        vector<vector<double>> all_current_positions = update_all_positions(celestial_bodies, num_dims); // (xy xy xy)

        // Use this vector to find positions of other bodies for each body 
        for (int i=0; i<num_bodies; i++)
        {
            celestial_bodies[i].other_positions = find_other_positions(all_current_positions, i); 
        }

        // calculate new acceleration
        # pragma omp parallel
        {
            # pragma omp for
            for (int i=0; i<num_bodies; i++)
            {
                celestial_bodies[i].calc_acceleration();
            }
        }
        
        // calculate new velocity
        # pragma omp parallel
        {
            # pragma omp for
            for (int i=0; i<num_bodies; i++)
            {
                celestial_bodies[i].update_velocity();
            }
        }
    }

    // Print length of simulation
    print_len_sim(max_time_in_hours);

    // Write positions to file
    write_to_file(csv_filename, max_time_in_hours, num_bodies, num_dims, celestial_bodies);
    
    return EXIT_SUCCESS;
}



// FUNCTIONS

vector<double> divide_vector_by_double(vector<double> v, double d)
/*
Function to divide every element in a vector by a double

inputs:
v: vector to be divided
d: number to divide each element in vector by

output:
new_vector: vector with each element divided by d.

*/
{
    vector<double> new_vector;
    if (d==0)
    {
        throw invalid_argument( "Can't divide by zero.");
    }
    else
    {
        double new_value;
        for (int i=0; i<v.size(); i++)
        {
            new_value = v[i] / d;
            new_vector.push_back(new_value);
        }
    }

    return new_vector;
}

vector<double> multiply_vector_by_double(vector<double> v, double d)
/*
Function to multiply every element in a vector by a double

inputs:
v: vector to be multiplied
d: number to multiply each element in vector by

output:
new_vector: vector with each element multiplied by d.

*/
{
    vector<double> new_vector;
    double new_value;
    for (int i=0; i<v.size(); i++)
    {
        new_value = v[i] * d;
        new_vector.push_back(new_value);
    }
    return new_vector;
}

vector<double> add_vectors(vector<double> v1, vector<double> v2)
/*
Function to add each corresponding element of two vectors together

inputs:
v1: first vector to be added
v2: second vector to be added

output:
new_vector: vector with each element added from two vectors.
*/
{
    vector<double> new_vector;
    if (v1.size() != v2.size())
    {
        throw invalid_argument("Vectors must be same length");
    }
    else
    {
        
        double new_value;
        for (int i=0; i<v1.size(); i++)
        {
            new_value = v1[i] + v2[i];
            new_vector.push_back(new_value);
        }
    }
    return new_vector;
}

double calc_mod(double x, double y)
/*
calculate distance between two points

inputs:
x: distance between x values of the two points
y: distance between y values of the two points

output:
modulus: distance between two points
*/
{
    // maybe use std::accumulate for more dimensions
    double modulus = sqrt( pow(x,2) + pow(y,2));
    return modulus;
}


void write_to_file(string filename, int num_iters, int num_bodies, int num_dims, vector<celestial_body> celestial_bodies)
/*
Function to write positions of all bodies to a csv file.

inputs
filename: name of csv to data write to
num_iters: number of measurements of position taken 
num_bodies: number of celestial bodies in the system
num_dims: number of dimensions in this system
celestial_bodies: vector of all the celestial body classes. 
*/
{
    cout << "writing to csv file: " << filename << endl;

    ofstream file_out (filename);

    vector<string> dimensions = {"x", "y", "z"};
    for (int i=-1; i<num_iters; i++)
    {
        for (int j=0; j<num_bodies; j++)
        {
            for (int k=0; k<num_dims; k++)
            {
                if (i==-1)
                {
                    file_out << dimensions[k] << to_string(j+1);
                }
                else
                {
                    file_out << celestial_bodies[j].all_positions[k][i];  
                }

                if (k!=num_dims-1 || j!=num_bodies-1)
                {
                    file_out << ",";
                }
            }
        }
        file_out << endl;
    }

}

template <typename T>
void print_vector(T vector)
/*
print a vector

Inputs
vector: vector to be inputted
*/
{
    for (auto i: vector) 
    {
        cout << i << ", ";
    }
    cout << endl;
}

vector<vector<double>> update_all_positions(vector<celestial_body> celestial_bodies, int num_dims)
/*
puts current position vectors of each celestial body into a single vector

Inputs
celestial_bodies: vector of all
num_dims: number of dimensions

Output
all_positions: vector of position vectors for each celestial body
*/
{   
    int num_bodies = celestial_bodies.size();
    vector<vector<double>> all_positions;
    for (int i=0;i<num_bodies;i++)
    {
        // positions of all bodies in eg x dimension
        vector<double> positions;
        for (int j=0;j<num_dims;j++)
        {
            positions.push_back(celestial_bodies[i].current_position[j]);
        }
        all_positions.push_back(positions);
    }
    return all_positions;
}

vector<vector<double>> find_other_positions(vector<vector<double>> all_current_positions, int index)
/*
Removes position of body at the provided index, to leave only the positions of the other bodies

Inputs:
all_current_positions: vector of position vectors for each celestial body
index: index in all_current_positions of body to be remove

Output:
other_positions: vector of position vectors of of each other celestial body. 

*/
{
    // Copy all current positions. 
    vector<vector<double>> other_positions = all_current_positions;
    // Remove current body from other positions 
    other_positions.erase(other_positions.begin()+index);
    return other_positions;    
}


void print_len_sim(int max_time_in_hours)
{

    double max_time_in_years = max_time_in_hours / (365*24);
    if (max_time_in_years > 1)
    {
        cout << "Number of simulated years: " << max_time_in_years << endl;
    }
    else if (max_time_in_hours < 24)
    {
        cout << "Number of simulated hours: " << max_time_in_hours << endl;
    }
    else
    {
        cout << "Number of simulated days: " << max_time_in_hours*365.25 << endl;
    }
}