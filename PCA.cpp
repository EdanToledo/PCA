//TLDEDA001

#include "PCA.h"

TLDEDA001::PCA::PCA(std::vector<std::vector<double>> dimensions, int NumberOfDimensions)
{
    this->dimensions = dimensions;
    this->NumberOfDimensions = NumberOfDimensions;
    for (int i = 0; i < this->dimensions.size(); i++)
    {
        float mean = 0;
        for (int j = 0; j < dimensions[i].size(); j++)
        {
            mean += dimensions[i][j];
        }

        mean /= dimensions[i].size();
        means.push_back(mean);
    }

    SubtractTheMean();
    CovarianceMatrix = CalculateCovarianceMatrix();
    EigenValues = CalculateEigenValues().reverse();

    EigenVectors = Eigen::MatrixXd(NumberOfDimensions, NumberOfDimensions);
    auto Temp = CalculateEigenVectors();
    EigenVectors.col(0) = Temp.col(1);
    EigenVectors.col(1) = Temp.col(0);

    TotalVariance = CovarianceMatrix.coeff(0, 0) + CovarianceMatrix.coeff(1, 1);
}

void TLDEDA001::PCA::SubtractTheMean()
{

    for (int i = 0; i < this->dimensions.size(); i++)
    {

        for (int j = 0; j < dimensions[i].size(); j++)
        {
            dimensions[i][j] - means[i];
        }
    }
}

const Eigen::MatrixXd TLDEDA001::PCA::CalculateCovarianceMatrix() const
{

    Eigen::MatrixXd cov(NumberOfDimensions, NumberOfDimensions);

    for (int i = 0; i < NumberOfDimensions; i++)
    {
        for (int j = 0; j < NumberOfDimensions; j++)
        {
            double val = 0;
            for (int k = 0; k < dimensions[i].size(); k++)
            {
                val += (dimensions[i][k] - means[i]) * (dimensions[j][k] - means[j]);
            }
            val /= (dimensions[i].size() - 1);
            cov(i, j) = val;
        }
    }
    return cov;
}

const Eigen::MatrixXd TLDEDA001::PCA::CalculateEigenVectors() const
{

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> Solver(CalculateCovarianceMatrix());

    return Solver.eigenvectors();
}

const Eigen::VectorXd TLDEDA001::PCA::CalculateEigenValues() const
{

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> Solver(CalculateCovarianceMatrix());

    return Solver.eigenvalues();
}

const std::vector<Eigen::MatrixXd> TLDEDA001::PCA::CreateFeatureVector() const
{

    Eigen::MatrixXd EigenVec = CalculateEigenVectors();

    std::vector<Eigen::MatrixXd> Feature;
    for (int i = 1; i <= NumberOfDimensions; i++)
    {
        Eigen::MatrixXd Temp(1, NumberOfDimensions);
        for (int j = 0; j < NumberOfDimensions; j++)
        {
            Temp(0, j) = EigenVec.coeff(j, NumberOfDimensions - i);
        }

        Feature.push_back(Temp);
    }
}

void TLDEDA001::PCA::TransformDataSet()
{

    std::vector<Eigen::MatrixXd> Feature = CreateFeatureVector();

    for (int i = 0; i < dimensions[0].size(); i++)
    {
        Eigen::MatrixXd Temp(2, 1);
        Temp(0, 0) = dimensions[0][i];
        Temp(1, 0) = dimensions[1][i];

        Eigen::MatrixXd Result1 = Feature[0] * Temp;
        Eigen::MatrixXd Result2 = Feature[1] * Temp;
        dimensions[0][i] = Result1.coeff(0, 0);
        dimensions[1][i] = Result2.coeff(0, 0);
    }
}

const double TLDEDA001::PCA::VarianceProportions(int PCAindex) const
{
    return this->EigenValues.coeff(PCAindex) / this->TotalVariance;
}

const Eigen::MatrixXd TLDEDA001::PCA::getCovarianceMatrix() const { return this->CovarianceMatrix; }

const Eigen::MatrixXd TLDEDA001::PCA::getEigenVectors() const { return this->EigenVectors; }

const Eigen::VectorXd TLDEDA001::PCA::getEigenValues() const { return this->EigenValues; }

const double TLDEDA001::PCA::getTotalVariance() const
{
    return this->TotalVariance;
}



std::ostream & TLDEDA001::operator<<(std::ostream & os , const TLDEDA001::PCA & PCAInstance){
    
    os << std::endl;

   os << "Results:" << std::endl;

    os<< std::endl;

    os << "Eigen Values: " << std::endl;
   os<< PCAInstance.getEigenValues() << std::endl;

    os << std::endl;

    os << "Eigen Vectors: " << std::endl;
    os << PCAInstance.getEigenVectors() << std::endl;

    os << std::endl;

    os << "Covariance Matrix: " << std::endl;
    os << PCAInstance.getCovarianceMatrix() << std::endl;

    os << std::endl;

    os << "Total Variance: " << std::endl;
    os << PCAInstance.getTotalVariance() << std::endl;

    os << std::endl;

    os << "Principal Component 1 Explains :" << PCAInstance.VarianceProportions(0)*100 <<"%" <<std::endl;
    os << std::endl;
    os << "Principal Component 2 Explains :" << PCAInstance.VarianceProportions(1)*100<<"%" << std::endl;

}