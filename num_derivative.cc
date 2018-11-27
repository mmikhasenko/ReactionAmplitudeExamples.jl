#include <iostream>
#include <functional>

std::function<double(double)> 
  deriv(std::function<double(double)> original) {

  auto func1 = [original](double x0)->double{
    double delta = 1e-5;
    return (original(x0+delta) - original(x0)) / delta;
  };

  return func1;
}

std::function<double(double,double)> 
deriv(std::function<double(double,double)> original) {

  auto func1 = [original](double x0, double y)->double{
    double delta = 1e-5;
    return (original(x0+delta,y) - original(x0,y)) / delta;
  };

  return func1;
}


int main() {

  auto func1 = [](double x)->double{return x*x*x*x;};
  auto dfunc1 = deriv(func1);

  std::cout << dfunc1(1.0) << "\n";

  auto func2 = [](double x, double y)->double{return x*x*x*x*y;};
  auto dfunc2 = deriv(func2);

  std::cout << dfunc2(1.0, 3.0) << "\n";

  return 1;
}
