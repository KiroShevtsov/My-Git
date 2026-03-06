#include <vector>
#include <stdexcept>
#include <ostream>
#include <map>
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

    Matrix operator+(const Matrix& o) const {
        if(nx_ != o.nx_ or ny_ != o.ny_){throw std::invalid_argument("Impossible to sum");}
        std::vector<double> res = mtx_;
        for(std::size_t i = 0; i < mtx_.size(); ++i){
            res[i] = mtx_[i] + o.mtx_[i];
        }
        return {nx_, ny_, res};
    }
    Matrix operator*(const Matrix& o) const {
        if(nx_ != o.ny_) throw std::invalid_argument("Impossible to multiply");
        std::size_t s = o.nx_ * ny_;
        std::vector<double> res(s, 0);
        for(std::size_t i = 0; i < ny_; ++i){
            for(std::size_t j = 0; j < o.nx_; ++j){
                double current = 0;
                for(std::size_t r = 0; r < nx_; ++r){
                    current += mtx_[i * nx_ + r] * o.mtx_[r * o.nx_ + j];
                }
                res[i * o.nx_ + j] = current;
            }
        }
        return {o.nx_, ny_, res};
    }
    Matrix operator*(double a) const {
        std::vector<double> res = mtx_;
        for(auto& item : res) {
            item *= a;
        }
        return {nx_, ny_, res};
    }
};
inline std::ostream& operator<<(std::ostream& os, const Matrix& mtx){
    for(std::size_t i = 0; i < mtx.ny_; ++i){
        for(std::size_t j = 0; j < mtx.nx_; ++j){
            os << mtx(i, j) << " ";
        }
        os << std::endl;
    }
    return os;
}
class CsrMatrix{
    std::vector<double> values_, cols_, rows_;
    using SparseMtx = std::map<std::pair<std::size_t, std::size_t>, double>;
    SparseMtx mtx_;
public:
    const std::size_t nx_;
    const std::size_t ny_;
    
    CsrMatrix(std::size_t nx, std::size_t ny, const SparseMtx& mtx) : nx_(nx), ny_(ny), mtx_(mtx){
        if(mtx.empty()) {throw std::invalid_argument("Mtx is empty");}
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
    CsrMatrix operator*(const CsrMatrix& o) const {
        if(nx_ != o.ny_){throw std::invalid_argument("Impossible to mult");}
        SparseMtx res;
        for(std::size_t i = 0; i < ny_; ++i){
            /*p - позиции ненулевых элементов в i-ой строке исходной матрицы*/
            for(std::size_t p = rows_[i]; p < rows_[i + 1]; ++p){
                std::size_t k = cols_[p];
                /*ненулевой элемент лежащий на пересечении i-строки и k-столбца исх.матрицы, a[i][k]*/
                double nnz_a = values_[p];
                /*q - позиции ненулевых элементов k - ой строки матрицы o_*/
                for(std::size_t q = o.rows_[k]; q < o.rows_[k + 1]; ++q){
                    std::size_t l = o.cols_[q];
                    /*ненулевой элемент лежащий на пересечении k-строки, l-столбца, b[k][l]*/                    
                    double nnz_b = o.values_[q];
                    res[std::make_pair(i, l)] += nnz_a * nnz_b;
                }
            }
        }
        return {o.nx_, ny_, res};
    }
    CsrMatrix operator+(const CsrMatrix& o) const {
        if(nx_ != o.nx_ or ny_ != o.ny_){throw std::invalid_argument("Impossible to sum");}
        SparseMtx res = mtx_;
        for(const auto& [key, value] : o.mtx_){
            res[key] += value;
        }
        for(const auto& [key, value] : res){
            if(value == 0){
                res.erase(key);
            }
        }
        return {nx_, ny_, res};
    }
    CsrMatrix operator*(double a) const {
        SparseMtx res = mtx_;
        for(auto& [key, value] : res){
            value *= a;
        }
        return {nx_, ny_, res};
    }
};