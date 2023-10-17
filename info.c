#include <stdio.h>
#include <math.h>

float gravity(float mass, float distance){
    float earthmass  = 5.9722*pow(10,24); //kg
    float G = 6.674*pow(10, -11);
    float dis = 6378137 + distance;
    float grav = G*((earthmass*mass)/(pow(dis,2)));
    return(grav);
}

// Function to calculate air resistance
float den(float y){
    float density = 0;
    float tem = 0;
    float presure = 0;
    if(y<=11000){
        tem = 15.4 - (0.00649*y);
        presure = 101.29 * pow(((tem + 273.1)/288.08), 5.256);
        return density = presure/(0.2869*(tem + 2773.1));
    }
    else if(11000<y && y<25000){
        tem = -56.46;
        presure = 22.655 * exp(1.73-0.000157*y);
        return density = presure/(0.2869*(tem + 2773.1));
    }
    else if (y >= 25000) {
        tem = -131.21 + (0.00299*y);
        presure = 2.488 * pow(((tem + 273.1)/216.6), -11.388);
        return density = presure/(0.2869*(tem + 2773.1));
    }
}

float airResistance(float y, float velocity, float dragCoefficient) {
    float density = den(y);
    return -0.5 * dragCoefficient * density * velocity * fabs(velocity);
}

// Function to calculate the desired trajectory
float desiredTrajectory(float x) {
    return -1000 * exp(-0.01 * x);
}

int main() {
    // Create a file
    FILE *fpt;
    fpt = fopen("2D_Simulation.csv", "w+");
    fprintf(fpt, "time, fuel, x, y, velocity_x, velocity_y, velocity,acceleration_x, acceleration_y, gravChange\n");

    // Initial conditions
    float mass = 198400; // kg
    float emptymass = 34500;  // kg
    float fuelMass = 163900; // kg
    float thrust = 3746000; // N
    float fuelRate = 1330; // kg/s
    float dragCoefficient = 0.23; // Example value for drag coefficient
    float timeStep = 0.1; // Time step (adjust as needed)

    // Initialize variables
    float x = 0; // Initial x position
    float y = 0; // Initial y position
    float velocity_x = 0; // Initial x velocity
    float velocity_y = 0; // Initial y velocity
    float acceleration_x = 0;
    float acceleration_y = 0;
    float gravChange = 0;

    for (float t = 0; t < 720; t += timeStep) {
        // Fuel consumption
        if (fuelMass > 0) {
            fuelMass -= fuelRate * timeStep;
            if (fuelMass < 0){
                fuelMass = 0 ;
                thrust = 0;
            }
        }

        // Gravity
        float trueWeight = emptymass + fuelMass;
        float grav = (gravity(trueWeight, y))/trueWeight;

        // Air resistance
        float velocity = sqrt(velocity_x * velocity_x + velocity_y * velocity_y);
        float airResistance_x = airResistance(y, velocity_x, dragCoefficient);
        float airResistance_y = airResistance(y, velocity_y, dragCoefficient);

        // Calculate acceleration
        acceleration_x = (thrust - airResistance_x) / trueWeight;
        
        // Desired y-acceleration to follow the trajectory
        float desiredY = desiredTrajectory(x);
        if(y<1000){
            acceleration_y = ((thrust - airResistance_y) / trueWeight ) - grav;            
        }
        else{
        acceleration_y = ((thrust - airResistance_y) / trueWeight ) - grav - desiredY;
        }
        // Update velocity
        velocity_x += acceleration_x * timeStep;
        velocity_y += acceleration_y * timeStep;
        velocity = (velocity_x*velocity_x)+(velocity_y*velocity_y)

        // Update position
        x += velocity_x * timeStep;
        y += velocity_y * timeStep;


        fprintf(fpt, "%.1f, %.1f , %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n", t, fuelMass, x, y, velocity_x, velocity_y, velocity,acceleration_x, acceleration_y, grav);
    }

    fclose(fpt);
    return 0;
}
