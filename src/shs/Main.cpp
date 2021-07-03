#include "SSModels.h"

#include <optional>
#include <stdio.h>

class A
{
public:
  A(std::optional<double> _t, std::optional<int> _b ) : t(_t), b(_b)
  {
    if (t.has_value())
    {
      printf("%f\n", *t );
    }
    else
    {
      printf("not set\n");
    }
  };

  std::optional<double> t;
  std::optional<int> b;

};

int main()
{
    A A(1.0, std::nullopt);

    // Define dynamics of stochastic process
    // Mode q1
    arma::mat Aq1 = {{0.43, 0.52}, {0.65, 0.12}};
    arma::mat Gq1 = {{1, 0.1}, {0, 0.1}};
    ssmodels_t modelq1(Aq1, Gq1);

    std::optional<arma::mat> A_matrix = modelq1.getAMatrix();

    if (A_matrix.has_value())
    {
      std::cout << "A matrix: " << *A_matrix << std::endl;
    }
    else
    {
      std::cout << "Incorrect init" << std::endl;
    }

    modelq1.checkModelStructure();

    return 0;
}
