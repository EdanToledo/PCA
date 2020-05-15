//TLDEDA001
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace TLDEDA001
{

class PCA
{

    int NumberOfDimensions;

    //Number of Dimensions of DATASET
    std::vector<std::vector<double>> dimensions;
    //Dimensional Means
    std::vector<double> means;

    Eigen::MatrixXd CovarianceMatrix;
    Eigen::MatrixXd EigenValues;
    Eigen::MatrixXd EigenVectors;

    double TotalVariance;

    //Inner Functions - not needed by User
    void SubtractTheMean();

    const Eigen::MatrixXd CalculateCovarianceMatrix() const;

    const Eigen::MatrixXd CalculateEigenVectors() const;

    const Eigen::VectorXd CalculateEigenValues() const;

public:
    PCA(std::vector<std::vector<double>> dimensions, int NumberOfDimensions);

    const std::vector<Eigen::MatrixXd> CreateFeatureVector() const;

    void TransformDataSet();

    const double VarianceProportions(int PCAindex)const;

    const Eigen::MatrixXd getCovarianceMatrix() const;

    const Eigen::MatrixXd getEigenVectors() const;

    const Eigen::VectorXd getEigenValues() const;

    const double getTotalVariance() const;

    friend std::ostream & operator<<(std::ostream & os , const PCA & instance);
};

} // namespace TLDEDA001