#include <iostream>
#include <ceres/ceres.h>
#include <cmath>

/*

Exercise 3.1 Solution of Numerical Optimization (2nd Edition) By Nocedal & Wright

Program the steepest descent and Newton algorithms using the backtracking line search, Algorithm 3.1. 
Use them to minimize the Rosenbrock function (2.22). 
Set the initial step length alpha = 1 and print the step length used by each method at each iteration.
First try the initial point x0 = (1.2, 1.2)^T and then the more difficult starting point x0 = (−1.2, 1)^T 

*/

class CostFunctor {
public:
  CostFunctor(){}

  template <typename T>
  bool operator()(const T* const x, T* residual) const {
    //This is the rosenbrock function (2.22) in the book.
    residual[0] = ceres::pow((x[1] - ceres::pow(x[0],2)),2)*100.0 + ceres::pow((-x[0] + 1.0),2);
    return true;
  }
};

struct resultsStruct{
    double steepestDescentStep;
    double newtonStep;
    resultsStruct(){
        steepestDescentStep = 0;
        newtonStep = 0;
    }
};

resultsStruct backtrackingLineSearch(
                        double alpha,
                        double rho,
                        double c,
                        const double* const* parameters,
                        ceres::CostFunction* costFunctor){


    std::cout << "---Backtracking Line Search---" << std::endl;
    resultsStruct result = resultsStruct();

    double jacobian[2]; 
    double* jacobians[1] = {jacobian};
    double residuals[1];

    /*
    
    Ceres Solver has a neat AutoDiff implementation that unfortunately only
    works up to the Jacobian out of the box. I use it here for the steepest
    descent method since I think it's cool.

    */

    costFunctor->Evaluate(parameters, residuals, jacobians);
    double ogResidual = residuals[0];
    double ogAlpha = alpha;

    double steepestPk[2] = {-(jacobians[0][0]),-(jacobians[0][1])};
    double directionalDerivative = (jacobians[0][0]*steepestPk[0] + jacobians[0][1] * steepestPk[1]);

    double xkplus1[2] = {parameters[0][0] + alpha * steepestPk[0],parameters[0][1] + alpha * steepestPk[1]};
    const double* xkplus_1[1] = {xkplus1};
    const double* const* nextIterate = xkplus_1;

    costFunctor->Evaluate(nextIterate, residuals, nullptr);

    double leftHandSide = residuals[0];
    double rightHandSide = ogResidual + c * alpha * directionalDerivative;

    while (leftHandSide >= rightHandSide){
        alpha *= rho;
        rightHandSide = ogResidual + c * alpha * directionalDerivative;
        xkplus1[0] = parameters[0][0] + alpha * steepestPk[0];
        xkplus1[1] = parameters[0][1] + alpha * steepestPk[1];
        costFunctor->Evaluate(nextIterate, residuals, nullptr);
        leftHandSide = residuals[0];
        std::cout << "Steepest Descent Stepsize: " << alpha << std::endl;
    }
    result.steepestDescentStep = alpha;

    /*
    
    After grappling with Ceres Solver's Jets and trying to nest them to compute the exact Hessian,
    I fell into templating and operator overload hell. Ceres Solver's Jet overloaded operators are only meant for
    Jet and scalar operations, and I tried doing nested Jet<Jet> and scalar ops. Maybe I could revisit this in the future.
    It kind of sucks because it's so easy to compute the exact Hessian using duals. Oh well...
    For simplicity, I compute the Hessian analytically.

    */
    
    alpha = ogAlpha;

    double h00 = -400*(parameters[0][1]) + 1200*(pow(parameters[0][0],2)) + 2;
    double hDiag = -400*(parameters[0][0]);
    double h11 = 200;
    
    double inverseCoeff = 1.0/ ((h00 * h11) - (hDiag*hDiag));
    double inverseH00 = h11 * inverseCoeff;
    double inverseH11 = h00 * inverseCoeff;
    double inverseHDiag = -hDiag * inverseCoeff;
    double newtonPk[2] = {-(jacobians[0][0]* inverseH00 + jacobians[0][1] * inverseHDiag),-(jacobians[0][0]* inverseHDiag + jacobians[0][1] * inverseH11)};

    xkplus1[0] = parameters[0][0] + alpha * newtonPk[0];
    xkplus1[1] = parameters[0][1] + alpha * newtonPk[1];

    directionalDerivative = (jacobians[0][0]*newtonPk[0] + jacobians[0][1] * newtonPk[1]);

    costFunctor->Evaluate(nextIterate, residuals, nullptr);

    leftHandSide = residuals[0];
    rightHandSide = ogResidual + c * alpha * directionalDerivative;

    while(leftHandSide >= rightHandSide){
        alpha *= rho;
        rightHandSide = ogResidual + c * alpha * directionalDerivative;
        xkplus1[0] = parameters[0][0] + alpha * newtonPk[0];
        xkplus1[1] = parameters[0][1] + alpha * newtonPk[1];
        costFunctor->Evaluate(nextIterate, residuals, nullptr);
        leftHandSide = residuals[0];
        std::cout << "Newton Stepsize: " << alpha << std::endl;
    }
    result.newtonStep = alpha;

    return result;
}

int main(){
    double x1 = 0;
    double x2 = 0;
    double alpha = 0;
    double c = 0;
    double rho = 0;
    std::string input = "";
    //No numerical constraints on the inputs. Try some weird inputs like a negative alpha and see what happens
    try{
        std::cout << "Enter value for x1: ";
        std::cin >> input;
        x1 = stod(input);

        std::cout << "Enter value for x2: ";
        std::cin >> input;
        x2 = stod(input);

        std::cout << "Enter value for alpha: ";
        std::cin >> input;
        alpha = stod(input);

        std::cout << "Enter value for rho: ";
        std::cin >> input;
        rho = stod(input);

        std::cout << "Enter value for c: ";
        std::cin >> input;
        c = stod(input);

    } catch (std::invalid_argument) {
        std::cerr << "Invalid argument, exiting." << std::endl;
        std::exit(1);
    }
    
    double x[2] = {x1, x2};
    const double* params[1] = {x};
    const double* const* parameters = params;

    auto rosenbrock = std::make_unique<ceres::AutoDiffCostFunction<CostFunctor,1,2>>(new CostFunctor());

    resultsStruct results = backtrackingLineSearch(alpha,rho,c,params, rosenbrock.get());

    std::cout << "---Final Step sizes---" << std::endl;
    std::cout << "Steepest Descent Method: " << results.steepestDescentStep << std::endl;
    std::cout << "Newton Method: " << results.newtonStep << std::endl;

    return 0;
}