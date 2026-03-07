#include <vector>
#include <stdexcept>
#include <ostream>
#include <map>
class Vector{
    std::vector<double> v_;
public:
    Vector(const std::vector<double>& v) : v_(v){}
    const std::vector<double>& data() const {return v_;}
    
    Vector operator+(const Vector& o) const {
        if(v_.size() != o.v_.size()) {throw std::invalid_argument("Impossible to sum vectors");}
        std::vector<double> res;
        res.reserve(v_.size());
        double sum = 0;
        for(std::size_t i = 0; i < v_.size(); ++i){
            sum = v_[i] + o.v_[i];
            res.push_back(sum);
        }
        return {res};
    }
    Vector operator-(const Vector& o) const {
        if(v_.size() != o.v_.size()) {throw std::invalid_argument("Impossible to sum vectors");}
        std::vector<double> res;
        res.reserve(v_.size());
        double dif = 0;
        for(std::size_t i = 0; i < v_.size(); ++i){
            dif = v_[i] - o.v_[i];
            res.push_back(dif);
        }
        return {res};
    }
    Vector operator*(double a) const {
        std::vector<double> res;
        res.reserve(v_.size());
        for(const auto& item : v_){
            res.push_back(a * item);
        }
        return {res};
    }
    inline double operator*(const Vector& o) const;

    double& operator[](std::size_t i) {return v_[i];}

    const double& operator[](std::size_t i) const {return v_[i];}
};
inline std::ostream& operator<<(std::ostream& os, const Vector& v) {
    for(const auto& item : v.data()){
        os << item << " " << std::endl;
    }
    return os;
}
double Vector::operator*(const Vector& o) const {
    if(v_.size() != o.v_.size()) {throw std::invalid_argument("Impossible to mult vectors");}
    double res = 0;
    for(std::size_t i = 0; i < v_.size(); ++i){
        res += v_[i] * o.v_[i];
    }
    return res;
}
/*nx - число столбцов, ny - число строк матрицы*/
class Matrix{
    std::vector<double> mtx_;
public:
    const std::size_t nx_;
    const std::size_t ny_;
    Matrix(std::size_t nx, std::size_t ny, const std::vector<double>& mtx) : nx_(nx), ny_(ny){
        if (mtx.empty() or mtx.size() != nx_ * ny_){
            throw std::invalid_argument("Matrix is empty");
        }
        mtx_ = mtx;
    }
    const double& operator()(std::size_t i, std::size_t j) const {return mtx_[i * nx_ + j];}

    double& operator()(std::size_t i, std::size_t j) {return mtx_[i * nx_ + j];}

    const std::vector<double>& data() const {return mtx_;}

    Vector operator*(const Vector& v) const {
        std::size_t v_size = v.data().size();
        if(nx_ != v_size) {throw std::invalid_argument("Impossible to mult on vector");}
        std::vector<double> res;
        res.reserve(ny_);
        for(std::size_t i = 0; i < ny_; ++i){
            double current = 0;
            for(std::size_t k = 0; k < v_size; ++k){
                current += mtx_[i * nx_ + k] * v[k];
            }
            res.push_back(current);
        }
        return {res};
    }
};
using SparseMtx = std::map<std::pair<std::size_t, std::size_t>, double>;
class CsrMatrix{
    std::vector<double> values_;
    std::vector<std::size_t> cols_, rows_;
public:
    const std::size_t nx_;
    const std::size_t ny_;
    
    CsrMatrix(std::size_t nx, std::size_t ny, const SparseMtx& mtx) : nx_(nx), ny_(ny){
        /*rows(0) = 0*/
        rows_.push_back(0);
        std::size_t non_zero = 0;
        for(std::size_t i = 0; i < ny_; ++i){
            for(std::size_t j = 0; j < nx_; ++j){
                if(mtx.count(std::make_pair(i, j))){
                    values_.push_back(mtx.at(std::make_pair(i, j)));
                    cols_.push_back(j);
                    non_zero++;
                }
            }
            rows_.push_back(non_zero);
        }
    };
    double operator()(std::size_t i, std::size_t j) const {
        std::size_t start = rows_[i];
        std::size_t end = rows_[i + 1];
        double res = 0;
        for(std::size_t k = start; k < end; ++k){
            if(cols_[k] == j){
                res = values_[k];
                return res;
            }
        }
        return res;
    }

    Vector operator*(const Vector& v) const {
        std::size_t v_size = v.data().size();
        if(nx_ != v_size) {throw std::invalid_argument("Impossible to mult on vector");}
        std::vector<double> res;
        res.reserve(ny_);
        for(std::size_t i = 0; i < ny_; ++i){
            double current = 0;
            for(std::size_t k = rows_[i]; k < rows_[i + 1]; ++k){
                std::size_t q = cols_[k];
                current += values_[k] * v[q];
            }
            res.push_back(current);
        }
        return {res};
    }
};
