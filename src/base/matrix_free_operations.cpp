// -----------------------------------------------------------------------------
// Matrix-free operations for tensor-product and Khatri–Rao structures
// Implementation in Rcpp (C++)
// -----------------------------------------------------------------

#include <Rcpp.h>
#include <numeric>
#include <iostream>
#include <math.h>

using namespace Rcpp;

// -----------------------------------------------------------------------------
// Computes the matrix–vector product (I_left ⊗ A ⊗ I_right) * x
// [[Rcpp::export]]
NumericVector mvp_normalfactor(
    const NumericMatrix& A, 
    int left, 
    int right, 
    const NumericVector& x
){
  
  int J = A.nrow();
  //int left = nf_size[1];
  //int right = nf_size[2];
  int base = 0;
  NumericVector v(x.size());

  for(int l=0; l<left; l++){
    for(int r=0; r<right; r++){

      int index = base+r;
      
      NumericVector z_in(J);
      for(int j=0; j<J; j++){
        z_in[j] = x[index];
        index += right;
      }
      
      NumericVector z_out(J);
      for(int nrow=0; nrow<J; nrow++){
       for(int ncol=0; ncol<J; ncol++){
         z_out[nrow] += A(nrow,ncol)*z_in[ncol];
       }
      }
      
      index = base+r;

      for(int j=0; j<J; j++){
        v[index] = z_out[j];
        index += right;
      }
      
    }
    
    base += (right*J);
    
  }
  
  return v;
}

// -----------------------------------------------------------------------------
// Computes the Khatri–Rao matrix–vector product (A_1 ⊙ ... ⊙ A_P) * x
// [[Rcpp::export]]
NumericVector mvp_khatrirao(
    const List& A_list, 
    const NumericVector& x
){
  
  int P=A_list.size();
  int n=x.size(); 
  int m=1;
  NumericVector m_p(P);
  for(int p=0; p<P; p++){
    const NumericMatrix& A_p = A_list[p];
    m_p[p]=A_p.nrow();
    m*=m_p[p];
  }
  NumericVector res(m);
  if(P==1){
    const NumericMatrix& A_0=A_list[0];
    for(int m_0=0; m_0<m; m_0++){
      for(int i=0; i<n; i++){
        res[m_0] += A_0(m_0,i)*x[i];
      }
    }
    return res;
  }
  else if(P==2){
    const NumericMatrix& A_0=A_list[0];
    const NumericMatrix& A_1=A_list[1];
    for(int i=0; i<n; i++){
      int j=0;
      for(int m_0=0; m_0<m_p[0]; m_0++){
        const double a0_m0_i=A_0(m_0,i);
        if(a0_m0_i==0){
          j+=m_p[1];
        }
        else{
          for(int m_1=0; m_1<m_p[1]; m_1++){
            res[j++] += a0_m0_i*A_1(m_1,i)*x[i];
          }
        }
      }
    }
    return res;
  }
  else if(P==3){
    const NumericMatrix& A_0=A_list[0];
    const NumericMatrix& A_1=A_list[1];
    const NumericMatrix& A_2=A_list[2];
    for(int i=0; i<n; i++){
      int j=0;
      for(int m_0=0; m_0<m_p[0]; m_0++){
        const double a0_m0_i=A_0(m_0,i);
        if(a0_m0_i==0){
          j+=m_p[1]*m_p[2];
        }
        else{
          for(int m_1=0; m_1<m_p[1]; m_1++){
            const double a1_m1_i=A_1(m_1,i);
            if(a1_m1_i==0){
              j+=m_p[2];
            } else{
              for(int m_2=0; m_2<m_p[2]; m_2++){
                res[j++] += a0_m0_i*a1_m1_i*A_2(m_2,i)*x[i];
              }
            }
          }
        }
      }
    }
    return res;
  }
  else if(P==4){
    const NumericMatrix& A_0=A_list[0];
    const NumericMatrix& A_1=A_list[1];
    const NumericMatrix& A_2=A_list[2];
    const NumericMatrix& A_3=A_list[3];
    for(int i=0; i<n; i++){
      int j=0;
      for(int m_0=0; m_0<m_p[0]; m_0++){
        const double a0_m0_i=A_0(m_0,i);
        if(a0_m0_i==0){
          j+=m_p[1]*m_p[2]*m_p[3];
        }
        else{
          for(int m_1=0; m_1<m_p[1]; m_1++){
            const double a1_m1_i=A_1(m_1,i);
            if(a1_m1_i==0){
              j+=m_p[2]*m_p[3];
            } else{
              for(int m_2=0; m_2<m_p[2]; m_2++){
                const double a2_m2_i=A_2(m_2,i);
                if(a2_m2_i==0){
                  j+=m_p[3];
                } else{
                  for(int m_3=0; m_3<m_p[3]; m_3++){
                    res[j++] += a0_m0_i*a1_m1_i*a2_m2_i*A_3(m_3,i)*x[i];
                  }
                }
              }
            }
          }
        }
      }
    }
    return res;
  }
  else{
    std::cout<<"P too large"<<std::endl;
    return 0;
  }
}

// -----------------------------------------------------------------------------
// Computes the transposed Khatri–Rao matrix-vector product (A_1 ⊙ ... ⊙ A_P)^T * y
// [[Rcpp::export]]
NumericVector mvp_transposed_khatrirao(
    const List& A_list, 
    const NumericVector& y
){
  
  int P=A_list.size();
  const NumericMatrix& A_0 = A_list[0];
  int n = A_0.ncol();
  NumericVector res(n);
  //std::vector<int> m_p(P);
  NumericVector m_p(P);
  //int m=1;
  for(int p=0; p<P; p++){
    const NumericMatrix& A_p = A_list[p];
    m_p[p]=A_p.nrow();
    //m*=m_p[p];
  }
  
  if(P==1){
    const NumericMatrix& A_0=A_list[0];
    for(int m_0=0; m_0<m_p[0]; m_0++){
      for(int i=0; i<n; i++){
        res[i] += A_0(m_0,i)*y[m_0];
      }
    }
    return res;
  }
  else if(P==2){
    const NumericMatrix& A_0=A_list[0];
    const NumericMatrix& A_1=A_list[1];
    for(int i=0; i<n; i++){
      int j=0;
      for(int m_0=0; m_0<m_p[0]; m_0++){
        const double a0_m0_i=A_0(m_0,i);
        if(a0_m0_i==0){
          j+=m_p[1];
        }
        else{
          for(int m_1=0; m_1<m_p[1]; m_1++){
            res[i] += a0_m0_i*A_1(m_1,i)*y[j++];
          }
        }
      }
    }
    return res;
  }
  else if(P==3){
    const NumericMatrix& A_0=A_list[0];
    const NumericMatrix& A_1=A_list[1];
    const NumericMatrix& A_2=A_list[2];
    for(int i=0; i<n; i++){
      int j=0;
      for(int m_0=0; m_0<m_p[0]; m_0++){
        const double a0_m0_i=A_0(m_0,i);
        if(a0_m0_i==0){
          j+=m_p[1]*m_p[2];
        }
        else{
          for(int m_1=0; m_1<m_p[1]; m_1++){
            const double a1_m1_i=A_1(m_1,i);
            if(a1_m1_i==0){
              j+=m_p[2];
            }
            else{
              for(int m_2=0; m_2<m_p[2]; m_2++){
                res[i] += a0_m0_i*a1_m1_i*A_2(m_2,i)*y[j++];
              }
            }
          }
        }
      }
    }
    return res;
  }
  else if(P==4){
    const NumericMatrix& A_0=A_list[0];
    const NumericMatrix& A_1=A_list[1];
    const NumericMatrix& A_2=A_list[2];
    const NumericMatrix& A_3=A_list[3];
    for(int i=0; i<n; i++){
      int j=0;
      for(int m_0=0; m_0<m_p[0]; m_0++){
        const double a0_m0_i=A_0(m_0,i);
        if(a0_m0_i==0){
          j+=m_p[1]*m_p[2]*m_p[3];
        }
        else{
          for(int m_1=0; m_1<m_p[1]; m_1++){
            const double a1_m1_i=A_1(m_1,i);
            if(a1_m1_i==0){
              j+=m_p[2]*m_p[3];
            }
            else{
              for(int m_2=0; m_2<m_p[2]; m_2++){
                const double a2_m2_i=A_2(m_2,i);
                if(a2_m2_i==0){
                  j+=m_p[3];
                } else{
                  for(int m_3=0; m_3<m_p[3]; m_3++){
                    res[i] += a0_m0_i*a1_m1_i*a2_m2_i*A_3(m_3,i)*y[j++];
                  }
                }
              }
            }
          }
        }
      }
    }
    return res;
  }
  else{
    std::cout<<"P too large"<<std::endl;
    return 0;
  }
}

// -----------------------------------------------------------------------------
// Computes the gram matrix-vector product (A Aᵀ) * x with A = A_1 ⊙ ... ⊙ A_P
// [[Rcpp::export]]
NumericVector mvp_gram_khatrirao(
    const List& A_list,
    const NumericVector& x
){
  NumericVector v = mvp_transposed_khatrirao(A_list, x);
  return mvp_khatrirao(A_list, v);
}

// -----------------------------------------------------------------------------
// Computes the diagonal of the gram matrix: diag(A Aᵀ) with A = A_1 ⊙ ... ⊙ A_P
// [[Rcpp::export]]
NumericVector diag_gram_khatrirao(
    const List& A_list
){
  
  int P=A_list.size();
  const NumericMatrix& A_0=A_list[0];
  int n=A_0.ncol(); 
  int K=1;
  std::vector<int> K_p(P);
  for(int p=0; p<P; p++){
    const NumericMatrix& A_p = A_list[p];
    K_p[p]=A_p.nrow();
    K*=K_p[p];
  }
  NumericVector diag(K);
  
  if(P==1){
    for(int i=0; i<n; i++){
      int j=0;
      const NumericMatrix& A_0 = A_list[0];
      for(int m_0=0; m_0<K_p[0]; m_0++){
        const double w_i_j =A_0(m_0,i);
        diag[j++] += w_i_j*w_i_j;
      }
    }
    return diag;
  }
  else if(P==2){
    for(int i=0; i<n; i++){
      int j=0;
      const NumericMatrix& A_0 = A_list[0];
      const NumericMatrix& A_1 = A_list[1];
      for(int m_0=0; m_0<K_p[0]; m_0++){
        const double a0_m0_i = A_0(m_0,i);
        for(int m_1=0; m_1<K_p[1]; m_1++){
          const double w_i_j = a0_m0_i*A_1(m_1,i);
          diag[j++] += w_i_j*w_i_j;
        }
      }
    }
    return diag;
  } 
  else if(P==3){
    for(int i=0; i<n; i++){
      int j=0;
      const NumericMatrix& A_0 = A_list[0];
      const NumericMatrix& A_1 = A_list[1];
      const NumericMatrix& A_2 =A_list[2];
      for(int m_0=0; m_0<K_p[0]; m_0++){
        const double a0_m0_i = A_0(m_0,i);
        for(int m_1=0; m_1<K_p[1]; m_1++){
          const double a1_m1_i = A_1(m_1,i);
          for(int m_2=0; m_2<K_p[2]; m_2++){
            const double w_i_j = a0_m0_i*a1_m1_i*A_2(m_2,i);
            diag[j++] += w_i_j*w_i_j;
          }
        }
      }
    }
    return diag;
  }
  else if(P==4){
    for(int i=0; i<n; i++){
      int j=0;
      const NumericMatrix& A_0 = A_list[0];
      const NumericMatrix& A_1 = A_list[1];
      const NumericMatrix& A_2 = A_list[2];
      const NumericMatrix& A_3 = A_list[3];
      for(int m_0=0; m_0<K_p[0]; m_0++){
        const double a0_m0_i = A_0(m_0,i);
        for(int m_1=0; m_1<K_p[1]; m_1++){
          const double a1_m1_i = A_1(m_1,i);
          for(int m_2=0; m_2<K_p[2]; m_2++){
            const double a2_m2_i = A_2(m_2,i);
            for(int m_3=0; m_3<K_p[3]; m_3++){
              const double w_i_j = a0_m0_i*a1_m1_i*a2_m2_i*A_3(m_3,i);
              diag[j++] += w_i_j*w_i_j;
            }
          }
        }
      }
    }
    return diag;
  } 
  else{
    std::cout<<"P too large"<<std::endl;
  }
  return diag;
  
}
