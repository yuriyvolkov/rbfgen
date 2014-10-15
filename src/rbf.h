#ifndef _RBF_NET_H_
#define _RBF_NET_H_
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <float.h>
// symbolic C++ computations headers
#include <vector.h>
#include <matrix.h>

//genetic algorithm headers
#include <ga/ga.h>
#include <rbfgagenome.h>

namespace ainet {

// Implementation of RBF artificial neural network
// see http://en.wikipedia.org/wiki/Radial_basis_function_network
template <typename T> class RBFnet  {
    size_t _input_size, _hidden_size;
    Matrix<T> _centers;
    Vector<T> _windows, _weights;
    const Matrix<T> *_vec_in;
    const Vector<T> *_vec_out;
    T fitness(GAGenome&);
public:
    RBFnet(size_t in_size, size_t hidden_size=0)
        :_input_size(in_size), _hidden_size(hidden_size),
        _centers(hidden_size,in_size),
        _windows(hidden_size),
        _weights(hidden_size),
        _vec_in(0),
        _vec_out(0)
        {};

    virtual ~RBFnet(){};

    //calculate the network output for the input vector v
    virtual T operator()(Vector<T>&v){
        Vector<T> hidden_layer_out(_windows.size(),0);
        for (size_t i=0;i < _centers.rows();i++){
            hidden_layer_out[i]=activation(v,_centers[i],_windows[i]);
        }
        return hidden_layer_out|_weights;
    };

    //calculate the network output for the input value v
    virtual T operator()(const T&v){
        if (_input_size > 1) {
            // throw smth here
        }
        Vector<T> in(1,v);
        Vector<T> hidden_layer_out(_windows.size(),0);
        for (size_t i=0;i < _centers.rows();i++){
            hidden_layer_out[i]=activation(in,_centers[i],_windows[i]);
        }
        return hidden_layer_out|_weights;
    };

    //radial basis activation function to calculate an output the hidden layer neuron
    //http://en.wikipedia.org/wiki/Radial_basis_function
    static T activation(const Vector<T>& x,const Vector<T>& c,const T& w){
        Vector<T> temp (x-c);
        T inner_product=temp|temp;
        return exp(-inner_product/pow(w,2));
    }

    //serialize network params to string
    std::string serialize() const;
    //load network params from file
    void load(const char* );

};

template<typename T> T RBFnet<T>::fitness(GAGenome&gene){
    RBFGenome& g= dynamic_cast<RBFGenome&> (gene);
    float score=FLT_MAX; double sm=0.01;

    //decode genotype into phenotype
    list <int> index;
    for(int i=0;i<g.binstr().length();i++)
        if(g.binstr().gene(i)&&( g.real().gene(i)>sm))
            index.push_back(i);

    if (index.size()>1){
        _hidden_size = index.size();
        _centers.resize(_hidden_size, _centers.cols());
        _windows.resize(_hidden_size);
        int    j=0;
        //aply phenotype to real object
        for    (list<int>::iterator i=index.begin();i!=index.end();i++){
            //centres
            _centers[j]=*_vec_in[*i];
            //window_width
            _windows[j]=g.real().gene(*i);
            ++j;
        }
        //create an approximation matrix
        Matrix<double> aMatr(_vec_in->rows(),_hidden_size);
        for (int i=0;i<aMatr.rows();i++)
            for(int j=0;j<aMatr.cols();j++)
                aMatr[i][j]=activation(*_vec_in[i],_centers[j],_windows[j]);
        //calculate G matrix
        Matrix<double> G=aMatr.transpose()*aMatr;
        //calculate b vector
        Vector<double> b=aMatr.transpose()*(*_vec_out);
        //calculate weights
        _weights=G.inverse()*b;

        //calculate the mean square error
        for (int i=0;i<_vec_out->size();i++){
            score+=pow(*this(*_vec_in[i])-*_vec_out[i],2);
        }
        score=pow(double(score),0.5);
    }

    return score;
}

template<typename T> std::string RBFnet<T>::serialize() const {
    std::ostringstream s;

    s<<_input_size<<std::endl<<_hidden_size<<std::endl;
    s<<_windows<<std::endl;
    s<<_weights<<std::endl;
    s<<_centers<<std::endl;
    /*
     for(int i=0;i<_centers.size();i++){
        s<<_windows[i]<<" ";
        for(int j=0;j<_input_size;j++)
            s<<_centers[i][j]<<" ";
        s<<_weights[i]<<std::endl;
     }
     */

     return s.str();
}

template<typename T> void RBFnet<T>::load(const char* fname){
    std::ifstream s(fname); int hidden_count=0;
    s>>_input_size; s>>_hidden_size;
    _centers.resize(_hidden_size, _input_size);
    _windows.resize(_hidden_size);
    _weights.resize(_hidden_size);
    //FIXME: add actual data load
    /*
    for (int i = 0; i <hidden_count; i++) {
        double v;
        s>>v;
        _windows[i]=v;
        _centers[i].resize(_input_size);

        for (int j=0;j<_input_size;j++){
            s>>v;
            _centers[i][j]=v;
        }
        s>>v;
        _weights[i]=v;
    }
    */
}

template<typename T> std::ostream& operator<<(std::ostream& s, const RBFnet<T>& net) {
    return s << net.serialize();
}

}
#endif
