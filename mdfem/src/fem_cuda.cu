#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/extrema.h>
#include <thrust/inner_product.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>

#include "fem.h"
//#define NeoHooken
#define __NTHREADS 256

void FEMFrame::cuda_memory_alloc() {
  assert(allocmultiple == 1);
  assert(_NP > 0);
  gpuErrchk(cudaMalloc( &_d_map23,         2 * sizeof(int)) );
  gpuErrchk( cudaMalloc( &_d_SR,           _NP*sizeof(G_Vector3)) );
  gpuErrchk( cudaMalloc( &_d_SRref,        _NP*sizeof(G_Vector3)) );
  gpuErrchk( cudaMalloc( &_d_Rref,         _NP*sizeof(G_Vector3)) );
  gpuErrchk( cudaMalloc( &_d_F,            _NP*sizeof(G_Vector3)) );
  gpuErrchk( cudaMalloc( &_d_fixed,        _NP*sizeof(int)) );
  gpuErrchk( cudaMalloc( &_d_EPOT_IND,     _NP*sizeof(double)) );
  gpuErrchk( cudaMalloc( &_d_EPOT_RMV,     _NP*sizeof(double)) );

  gpuErrchk(cudaMalloc(&_d_H_element, 3*3*sizeof(double))); 
  gpuErrchk( cudaMalloc( &_d_inv_elements, _NP*_MAX_NELEM_SHARE_NODE*sizeof(int)) );

#ifdef DEBUG_USECUDA
  Realloc( _h_d_map23,    int,             2);
  Realloc( _h_d_SR,       G_Vector3,      _NP);
  Realloc( _h_d_SRref,    G_Vector3,      _NP);
  Realloc( _h_d_Rref,     G_Vector3,      _NP);
  Realloc( _h_d_F,        G_Vector3,      _NP);
  Realloc( _h_d_fixed,    int,            _NP);
  Realloc( _h_d_EPOT_IND, double,         _NP);
  Realloc( _h_d_EPOT_RMV, double,         _NP);
  Realloc( _h_d_EPOT,     double,         _NP);
  Realloc( _h_d_H_element,  double,         3*3); 
  Realloc( _h_d_inv_elements,int,         _NP*_MAX_NELEM_SHARE_NODE); //>>>>>>
#endif

}

void FEMFrame::cuda_memory_alloc_elements() {
  assert(_NELE > 0);
  gpuErrchk( cudaMalloc( &_d_EPOT,             _NELE*sizeof(double)) );
  gpuErrchk( cudaMalloc( &_d_elements,         _NELE*_NNODE_PER_ELEMENT*sizeof(int)) );

#ifdef DEBUG_USECUDA
  Realloc( _h_d_EPOT,        double,  _NELE);
  Realloc( _h_d_elements,    int,     _NELE*_NNODE_PER_ELEMENT);
#endif
}

void FEMFrame::cuda_memory_alloc_element_coeff() {
  assert(_NDIM > 0 && _NELE > 0 && _NNODE_PER_ELEMENT > 0 && _NINT_PER_ELEMENT > 0);

  int size =  _NINT_PER_ELEMENT*_NDIM*_NDIM*_NDIM*_NNODE_PER_ELEMENT*_NELE;
  gpuErrchk( cudaMalloc( &_d_gauss_weight, _NINT_PER_ELEMENT*sizeof(double)) );
  gpuErrchk( cudaMalloc( &_d_dFdu,         size*sizeof(double)) );
  gpuErrchk( cudaMalloc( &_d_F_padding,    _NELE*_NNODE_PER_ELEMENT*sizeof(G_Vector3)) );

#ifdef DEBUG_USECUDA
  Realloc( _h_d_gauss_weight,double,     _NINT_PER_ELEMENT);
  Realloc( _h_d_dFdu,        double,     size);
  Realloc( _h_d_F_padding,   G_Vector3,  _NELE*_NNODE_PER_ELEMENT);
  Realloc( _h_EPOT,          double,     _NELE);
#endif
}

void FEMFrame::cuda_memcpy_all() {
  assert(sizeof(G_Vector3) == sizeof(Vector3));
  assert(sizeof(G_Matrix33) == sizeof(Matrix33));
  assert(_NELE > 0); assert(_NNODE_PER_ELEMENT>0);assert(_NINT_PER_ELEMENT>0);
  assert(_H[0][0]>0 && _H[1][1]>0 && _H[2][2]>0);
  int size1 = _NELE*_NNODE_PER_ELEMENT;
  int size2 = _NP*_MAX_NELEM_SHARE_NODE;
  int size3 =  _NINT_PER_ELEMENT*_NDIM*_NDIM*_NDIM*_NNODE_PER_ELEMENT*_NELE;

  gpuErrchk( cudaMemcpy( _d_map23,    map23,      2*sizeof(int),         cudaMemcpyHostToDevice) );
  gpuErrchk(  cudaMemcpy( _d_SR,      _SR,        _NP*sizeof(G_Vector3), cudaMemcpyHostToDevice) );
  gpuErrchk(  cudaMemcpy( _d_SRref,   _SRref,     _NP*sizeof(G_Vector3), cudaMemcpyHostToDevice) );
  gpuErrchk(  cudaMemcpy( _d_Rref,    _Rref,      _NP*sizeof(G_Vector3), cudaMemcpyHostToDevice) );
  gpuErrchk(  cudaMemcpy( _d_F,       _F,         _NP*sizeof(G_Vector3), cudaMemcpyHostToDevice) );
  gpuErrchk(  cudaMemcpy( _d_fixed,   fixed,      _NP*sizeof(int),       cudaMemcpyHostToDevice) );
  gpuErrchk(  cudaMemcpy( _d_H_element, _H.element, 3*3*sizeof(double),  cudaMemcpyHostToDevice) );
  gpuErrchk(  cudaMemcpy(_d_elements,   elements,    size1*sizeof(int),cudaMemcpyHostToDevice) );
  create_inverse_connectivities_matrix();
  gpuErrchk(  cudaMemcpy(_d_inv_elements,inv_elements,size2*sizeof(int),cudaMemcpyHostToDevice) );
  gpuErrchk(  cudaMemcpy(_d_gauss_weight,gauss_weight,_NINT_PER_ELEMENT*sizeof(double),cudaMemcpyHostToDevice) );
  gpuErrchk(  cudaMemcpy(_d_dFdu,        dFdu,       size3*sizeof(double),cudaMemcpyHostToDevice) );
}

__device__ G_Matrix33 getEigenF(G_Vector3 p, G_Matrix33 Fdef, 
                                double y_eigen_zbound_max, double y_eigen_zbound_min,
                                double x_eigen_zbound_max, double x_eigen_zbound_min,
                                double y_eigen_strain, double x_eigen_strain ) {
  G_Matrix33 I; I.eye();    
   if (p[2] <= y_eigen_zbound_max && p[2] >= y_eigen_zbound_min ){ 
    I[1][1] = y_eigen_strain;
  }
  if (p[2] <= x_eigen_zbound_max && p[2] >= x_eigen_zbound_min ){ 
    I[0][0] = x_eigen_strain;
  }
  return I;
}


__global__ void kernel_beam_fem_energy_force(int _NDIM, int _NELE, int _NINT_PER_ELEMENT, int _NNODE_PER_ELEMENT, int *_d_elements, int *_d_inv_elements, int *_d_map23, int *_d_fixed, double *_d_gauss_weight, double *_d_dFdu, double *_d_EPOT, G_Vector3 *_d_SR, G_Vector3 *_d_SRref, G_Vector3 *_d_Rref, G_Vector3 *_d_F, G_Vector3 *_d_F_padding, double* _d_H_element, double __V0 ) {

  for (int iele = blockDim.x * blockIdx.x + threadIdx.x;iele < _NELE; iele+= blockDim.x * gridDim.x) { 
    int i,j,jpt,iA, in, ip, iq, ind;
    G_Vector3 dsj, drj, elem_center;
    G_Matrix33 Fe, Fdef, B, C, Cinv, FCinvtran, dEdF, hinv;
    G_Matrix33 eigF, invEigF;
    double Eele, Eint, Jet, trace, detB, J2;
    G_Matrix33 E, E2, I, pk2, temp1, temp2;    

    G_Matrix33 _d_H(_d_H_element);
    I.eye();
    /* energy of this element */
    hinv = _d_H.inv();
    Eele = 0;

    /* center of the element */
    elem_center.clear();elem_center[0] = elem_center[1]=elem_center[2]= 0;
    for(j=0;j<_NNODE_PER_ELEMENT;j++) {
      jpt=_d_elements[iele*_NNODE_PER_ELEMENT+j];	
      for (i = 0;i<_NDIM;i++)  {
        elem_center[i] += 1.0/_NNODE_PER_ELEMENT *_d_Rref[jpt][i];
      }
    }
    for(j=0;j<_NNODE_PER_ELEMENT;j++) {
      jpt=iele*_NNODE_PER_ELEMENT+j;
      _d_F_padding[jpt].clear();
    }

    for(iA=0;iA<_NINT_PER_ELEMENT;iA++) {
      /* energy of this Gauss integration point */
      Eint = 0;
      /* deformation gradient */
      Fdef.clear(); Fdef[0][0] = Fdef[1][1] = Fdef[2][2] = 1.0;
      for(j=0;j<_NNODE_PER_ELEMENT;j++) {
        jpt=_d_elements[iele*_NNODE_PER_ELEMENT+j];
        _d_SRref[jpt ] = hinv*_d_Rref[jpt];
        dsj = _d_SR[jpt] - _d_SRref[jpt];
        dsj.subint();
        _d_H.multiply(dsj, drj);
        /* Add contribution from drj to F using dFdu */
        for(ip=0;ip<_NDIM;ip++) {
          for(iq=0;iq<_NDIM;iq++) {
            for(in=0;in<_NDIM;in++) {
              ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
              if (_NDIM == 2)
                Fdef[ip][iq] += _d_dFdu[ind]*drj[_d_map23[in]]; 
              else //_NDIM == 3
                Fdef[ip][iq] += _d_dFdu[ind]*drj[in];
        } } }
      }

      E = Fdef.tran()*Fdef-I;
      B = Fdef*Fdef.tran();
      C = Fdef.tran()*Fdef;
      Cinv = C.inv();
	    FCinvtran = Fdef*Cinv.tran();
      Jet = Fdef.det();  J2 = Jet*Jet;
      detB = B.det(); 
      Eint = 0;

      /* Add contribution from F to Eint */ 
      if (_NDIM == 2) {
	      double MU = 1;
	      trace = B[0][0] + B[1][1];
	      Eint = 0.5*(trace + 1.0/detB - 3); /* multiply MU and V0 later */
	      dEdF = FCinvtran*(-1.0/J2) + Fdef; /* multiply MU later */
	      Eint *= MU;
	      dEdF *= MU;

        /* Add contribution from drj to F using dFdu */
        for(j=0;j<_NNODE_PER_ELEMENT;j++) {
          jpt=iele*_NNODE_PER_ELEMENT+j;
          for(ip=0;ip<_NDIM;ip++) {
            for(iq=0;iq<_NDIM;iq++) {
    	      for(in=0;in<_NDIM;in++) {
    	        ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
              if (1) {// _d_fixed[jpt] == 0) { 
    	          _d_F_padding[jpt][_d_map23[in] ] -= _d_gauss_weight[iA] *dEdF[ip][iq]*_d_dFdu[ind] * __V0; 
//                printf("_d_Fpadding[%d][%d] = %g\n", jpt, _d_map23[in], _d_F_padding[jpt][_d_map23[in] ]);
              }
    	} } } }
      }
      
      Eele += _d_gauss_weight[iA]*Eint *__V0;
    }
    _d_EPOT[iele] = Eele;
  }
}


__global__ void kernel_snap_fem_energy_force(int _NDIM, int _NELE, int _NINT_PER_ELEMENT, int _NNODE_PER_ELEMENT, int *_d_elements, int *_d_inv_elements, int *_d_map23, int *_d_fixed, double *_d_gauss_weight, double *_d_dFdu, double *_d_EPOT, G_Vector3 *_d_SR, G_Vector3 *_d_SRref, G_Vector3 *_d_Rref, G_Vector3 *_d_F, G_Vector3 *_d_F_padding, double* _d_H_element, double y_eigen_zbound_max, double y_eigen_zbound_min, double x_eigen_zbound_max, double x_eigen_zbound_min, double y_eigen_strain, double x_eigen_strain, double __V0 ) {

  for (int iele = blockDim.x * blockIdx.x + threadIdx.x;iele < _NELE; iele+= blockDim.x * gridDim.x) { 
    int i,j,jpt,iA, in, ip, iq, ind, p, q, r;
    G_Vector3 dsj, drj, elem_center;
    G_Matrix33 Fe, Fdef, B, C, Cinv, FCinvtran, dEdF, hinv;
    G_Matrix33 eigF, invEigF;
    double Eele, Eint, Jet, trace, detB, J2;
    double lambda, mu, temp;
    G_Matrix33 E, E2, I, pk2, temp1, temp2;    
    mu = 1; lambda = 1.95; 

    G_Matrix33 _d_H(_d_H_element);
    I.eye();
    /* energy of this element */
    hinv = _d_H.inv();
    Eele = 0;

    /* center of the element */
    elem_center.clear();elem_center[0] = elem_center[1]=elem_center[2]= 0;
    for(j=0;j<_NNODE_PER_ELEMENT;j++) {
      jpt=_d_elements[iele*_NNODE_PER_ELEMENT+j];	
      for (i = 0;i<_NDIM;i++)  {
        elem_center[i] += 1.0/_NNODE_PER_ELEMENT *_d_Rref[jpt][i];
      }
    }
    for(j=0;j<_NNODE_PER_ELEMENT;j++) {
      jpt=iele*_NNODE_PER_ELEMENT+j;
      _d_F_padding[jpt].clear();
    }

    for(iA=0;iA<_NINT_PER_ELEMENT;iA++) {
      /* energy of this Gauss integration point */
      Eint = 0;
      /* deformation gradient */
      Fdef.clear(); Fdef[0][0] = Fdef[1][1] = Fdef[2][2] = 1.0;
      for(j=0;j<_NNODE_PER_ELEMENT;j++) {
        jpt=_d_elements[iele*_NNODE_PER_ELEMENT+j];
        _d_SRref[jpt ] = hinv*_d_Rref[jpt];
        dsj = _d_SR[jpt] - _d_SRref[jpt];
        dsj.subint();
        _d_H.multiply(dsj, drj);
        /* Add contribution from drj to F using dFdu */
        for(ip=0;ip<_NDIM;ip++) {
          for(iq=0;iq<_NDIM;iq++) {
            for(in=0;in<_NDIM;in++) {
              ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
              if (_NDIM == 2)
                Fdef[ip][iq] += _d_dFdu[ind]*drj[_d_map23[in]]; 
              else //_NDIM == 3
                Fdef[ip][iq] += _d_dFdu[ind]*drj[in];
        } } }
      }

      eigF = getEigenF(elem_center, Fdef, y_eigen_zbound_max, y_eigen_zbound_min, x_eigen_zbound_max, x_eigen_zbound_min, y_eigen_strain, x_eigen_strain );

      invEigF = eigF.inv();
      Fe = Fdef*invEigF;     
      E = Fe.tran()*Fe-I;
      B = Fe*Fe.tran();
      C = Fe.tran()*Fe;
      Cinv = C.inv();
      Jet = Fdef.det();  J2 = Jet*Jet;
      detB = B.det(); 
      Eint = 0;

      /* Add contribution from F to Eint */ 
      if (_NDIM == 2) {
        trace = B[0][0] + B[1][1];
        Eint = 0.5*(trace + 1.0/detB - 3); /* multiply MU and V0 later */
        dEdF = FCinvtran*(-1.0/J2) + Fdef; /* multiply MU later */
        /* Add contribution from drj to F using dFdu */
        for(j=0;j<_NNODE_PER_ELEMENT;j++) {
          jpt=iele*_NNODE_PER_ELEMENT+j;
          for(ip=0;ip<_NDIM;ip++) {
            for(iq=0;iq<_NDIM;iq++) {
    	      for(in=0;in<_NDIM;in++) {
    	        ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
    	        _d_F_padding[jpt][_d_map23[in] ] -= _d_gauss_weight[iA] * dEdF[ip][iq]*_d_dFdu[ind] * __V0; 
    	} } } }
      }
      
      /* Add contribution from F to Eint */ 
      if (_NDIM == 3) {
		
#if defined NeoHooken
        double C10 = 1; double D = 1e-1;
        double jet23 = pow(Jet, -2.0/3.0);
        double Ibar = B.trace() * jet23;
        Eint = C10*(Ibar-3) + 1.0/D *(Jet-1)*(Jet-1);
        dEdF = (Fdef*(invEigF*(invEigF.tran())) * 2.0*jet23 - ((Fdef.inv()).tran())*2.0/3.0 * Ibar)*C10 + ((Fdef.inv()).tran())*2.0/D *(Jet-1)*Jet;

#else
        E2 = E*E;
        dEdF.clear();
        for (i = 0;i<_NDIM;i++)
          for (j = 0;j<_NDIM;j++)
            for (p = 0;p<_NDIM;p++)
              for (r = 0;r<_NDIM;r++) {
                temp = 0;
                for (q = 0;q<_NDIM;q++)
                  temp += 2*mu*invEigF[j][p]*Fdef[i][r]*invEigF[r][q]*E[p][q] + 
                          2*mu*invEigF[r][p]*Fdef[i][r]*invEigF[j][q]*E[p][q];
                  dEdF[i][j] += 0.5*(2*lambda*E.trace()*invEigF[j][p]*Fdef[i][r] * invEigF[r][p] + temp);
              }
        Eint = 0.5*lambda*(E.trace())*(E.trace()) + mu*(E2.trace());
#endif
        /* Add contribution from drj to F using dFdu */

    	for(j=0;j<_NNODE_PER_ELEMENT;j++) {
    	  jpt=iele*_NNODE_PER_ELEMENT+j;
    	  for(ip=0;ip<_NDIM;ip++) {
    	    for(iq=0;iq<_NDIM;iq++) {
    	      for(in=0;in<_NDIM;in++) {
    		ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
    		_d_F_padding[jpt][in ] -= _d_gauss_weight[iA] * dEdF[ip][iq]*_d_dFdu[ind] * __V0; 				
    	} } } }
      }
      /* Add contribution from Eint to Eele */
      Eele += _d_gauss_weight[iA]*Eint *__V0;
    }
    _d_EPOT[iele] = Eele;
  }
}

__global__ void kernel_assemble_back_force(int _NP, int *_d_fixed, int *_d_inv_elements, G_Vector3 *_d_F, G_Vector3 *_d_F_padding) { 
  for(int ind = blockDim.x * blockIdx.x + threadIdx.x; ind<_NP;ind+=blockDim.x*gridDim.x) {
    _d_F[ind].clear();  
    for (int k = 0; k<_MAX_NELEM_SHARE_NODE; k++) { 
      int indice = ind *_MAX_NELEM_SHARE_NODE+k;
      if (_d_inv_elements[ indice ] >= 0 && _d_fixed[ind] == 0) { 
         _d_F[ind] += _d_F_padding[ _d_inv_elements[indice] ];
      }
    }
  }
}

void FEMFrame::cuda_beam_fem_energy_force() {
  DUMP("FEM");
  gpuErrchk( cudaMemcpy( _d_SR,        _SR,        _NP*sizeof(G_Vector3), cudaMemcpyHostToDevice));
  gpuErrchk( cudaMemcpy( _d_H_element, _H.element, 3*3*sizeof(double),    cudaMemcpyHostToDevice));
  gpuErrchk( cudaMemcpy( _d_fixed,     fixed,      _NP*sizeof(int),       cudaMemcpyHostToDevice));

  _EPOT=0; //put this back
  for(int i=0;i<_NP;i++) {
    _EPOT_IND[i]=0; _EPOT_RMV[i]=0; _VIRIAL_IND[i].clear();
  }

  _VIRIAL.clear();

#ifdef DEBUG_USECUDA
  assert(check_host_device_memory_transfer() == 0);
#endif

  kernel_beam_fem_energy_force<<<(_NELE+255)/256,__NTHREADS>>>(_NDIM, _NELE, _NINT_PER_ELEMENT, _NNODE_PER_ELEMENT, _d_elements, _d_inv_elements, _d_map23, _d_fixed, _d_gauss_weight, _d_dFdu, _d_EPOT, _d_SR, _d_SRref, _d_Rref, _d_F, _d_F_padding, _d_H_element, __V0);

  //cudaDeviceSynchronize();

#ifdef DEBUG_USECUDA
  assert(check_host_device_memory_transfer() == 0);
#endif

 kernel_assemble_back_force<<<(_NP+255)/256,__NTHREADS>>>(_NP, _d_fixed, _d_inv_elements, _d_F, _d_F_padding);

  // cudaDeviceSynchronize();

#ifdef DEBUG_USECUDA
  gpuErrchk(cudaMemcpy(_h_EPOT,_d_EPOT, _NELE*sizeof(double), cudaMemcpyDeviceToHost));
  for (int i = 0; i<_NELE; i++) { printf("_h_EPOT[%d]= %g, _NELE = %d\n", i,_h_EPOT[i],_NELE); _EPOT += _h_EPOT[i]; } 
#else
  /* Reduce potential energy to CPU for relax function to call */ 
  thrust::device_ptr<double> t_EPOT = thrust::device_pointer_cast(_d_EPOT);
  _EPOT = thrust::reduce(t_EPOT,t_EPOT+_NELE); 
#endif

  /* Copy force (Vector3 *) back to CPU for relax function to call */ 
  cudaMemcpy( _F,       _d_F,         _NP*sizeof(G_Vector3), cudaMemcpyDeviceToHost);
}

void FEMFrame::cuda_snap_fem_energy_force() {
  DUMP("FEM");
  gpuErrchk( cudaMemcpy( _d_SR,        _SR,        _NP*sizeof(G_Vector3), cudaMemcpyHostToDevice));
  gpuErrchk( cudaMemcpy( _d_H_element, _H.element, 3*3*sizeof(double),    cudaMemcpyHostToDevice));
  gpuErrchk( cudaMemcpy( _d_fixed,     fixed,      _NP*sizeof(int),       cudaMemcpyHostToDevice));

  _EPOT=0; //put this back
  for(int i=0;i<_NP;i++) {
    _EPOT_IND[i]=0; _EPOT_RMV[i]=0; _VIRIAL_IND[i].clear();
  }

  _VIRIAL.clear();

#ifdef DEBUG_USECUDA
  assert(check_host_device_memory_transfer() == 0);
#endif

  kernel_snap_fem_energy_force<<<(_NELE+255)/256,__NTHREADS>>>(_NDIM, _NELE, _NINT_PER_ELEMENT, _NNODE_PER_ELEMENT, _d_elements, _d_inv_elements, _d_map23, _d_fixed, _d_gauss_weight, _d_dFdu, _d_EPOT, _d_SR, _d_SRref, _d_Rref, _d_F, _d_F_padding, _d_H_element,y_eigen_zbound_max, y_eigen_zbound_min,  x_eigen_zbound_max, x_eigen_zbound_min, y_eigen_strain,     x_eigen_strain , __V0);

  //cudaDeviceSynchronize();

#ifdef DEBUG_USECUDA
  assert(check_host_device_memory_transfer() == 0);
#endif

 kernel_assemble_back_force<<<(_NP+255)/256,__NTHREADS>>>(_NP, _d_fixed, _d_inv_elements, _d_F, _d_F_padding);

  // cudaDeviceSynchronize();

#ifdef DEBUG_USECUDA
  gpuErrchk(cudaMemcpy(_h_EPOT,_d_EPOT, _NELE*sizeof(double), cudaMemcpyDeviceToHost));
  for (int i = 0; i<_NELE; i++) { printf("_h_EPOT[%d]= %g, _NELE = %d\n", i,_h_EPOT[i],_NELE); _EPOT += _h_EPOT[i]; } 
#else
  /* Reduce potential energy to CPU for relax function to call */ 
  thrust::device_ptr<double> t_EPOT = thrust::device_pointer_cast(_d_EPOT);
  _EPOT = thrust::reduce(t_EPOT,t_EPOT+_NELE); 
#endif

  /* Copy force (Vector3 *) back to CPU for relax function to call */ 
  cudaMemcpy( _F,       _d_F,         _NP*sizeof(G_Vector3), cudaMemcpyDeviceToHost);
}

void FEMFrame::free_device_ptr() {
  cudaFree(_d_elements);
  cudaFree(_d_inv_elements);
  cudaFree(_d_map23);
  cudaFree(_d_fixed);
  cudaFree(_d_colorids);
  cudaFree(_d_gauss_weight); 
  cudaFree(_d_dFdu);
  cudaFree(_d_EPOT);
  cudaFree(_d_EPOT_IND);
  cudaFree(_d_EPOT_RMV);
  cudaFree(_d_SR);
  cudaFree(_d_SRref);
  cudaFree(_d_Rref);        
  cudaFree(_d_F);
  cudaFree(_d_F_padding);
}

/* This is a simple test for GPU. run the function to see if maxErro == 0. If not, GPU device is not set correctly */
__global__ void saxpy(int n, float a, const float *x, float *y)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) y[i] = a*x[i] + y[i];
}
int FEMFrame::test_saxpy(void)
{
  int N = 1<<20;
  float *x, *y, *d_x, *d_y;
  x = (float*)malloc(N*sizeof(float));
  y = (float*)malloc(N*sizeof(float));
  cudaMalloc(&d_x, N*sizeof(float));
  cudaMalloc(&d_y, N*sizeof(float));
  for (int i = 0; i < N; i++) { x[i] = 1.0f; y[i] = 2.0f; }
  cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, y, N*sizeof(float), cudaMemcpyHostToDevice);

  // Perform SAXPY on 1M elements
  saxpy<<<(N+255)/256, __NTHREADS>>>(N, 2.0f, d_x, d_y);
  cudaMemcpy(y, d_y, N*sizeof(float), cudaMemcpyDeviceToHost);
  float maxError = 0.0f;
  for (int i = 0; i < N; i++) maxError = max(maxError, abs(y[i]-4.0f));
  INFO_Printf("Max error: %f\n", maxError);
  cudaFree(d_x);
  cudaFree(d_y);
  free(x);
  free(y);
  return 0;
}

#ifdef DEBUG_USECUDA
int FEMFrame::check_host_device_memory_transfer() 
{
  INFO_Printf("I am in check_host_device memory transfer\n");
  assert(sizeof(G_Vector3) == sizeof(Vector3));
  assert(sizeof(G_Matrix33) == sizeof(Matrix33));
  assert(_NP>0);
  assert(_NELE > 0); assert(_NNODE_PER_ELEMENT>0);assert(_NINT_PER_ELEMENT>0);
  assert(_H[0][0]>0 && _H[1][1]>0 && _H[2][2]>0);
  int size1 = _NELE*_NNODE_PER_ELEMENT;
  int size2 = _NP*_MAX_NELEM_SHARE_NODE;

  gpuErrchk( cudaMemcpy( _h_d_map23,   _d_map23,      2*sizeof(int),        cudaMemcpyDeviceToHost));
  gpuErrchk( cudaMemcpy( _h_d_SR,      _d_SR,        _NP*sizeof(G_Vector3), cudaMemcpyDeviceToHost));
  gpuErrchk( cudaMemcpy( _h_d_SRref,   _d_SRref,     _NP*sizeof(G_Vector3), cudaMemcpyDeviceToHost));
  gpuErrchk( cudaMemcpy( _h_d_Rref,    _d_Rref,      _NP*sizeof(G_Vector3), cudaMemcpyDeviceToHost));
  gpuErrchk( cudaMemcpy( _h_d_F,       _d_F,         _NP*sizeof(G_Vector3), cudaMemcpyDeviceToHost));
  gpuErrchk( cudaMemcpy( _h_d_fixed,   _d_fixed,      _NP*sizeof(int),      cudaMemcpyDeviceToHost));
  gpuErrchk( cudaMemcpy( _h_d_H_element,_d_H_element, 3*3*sizeof(double),    cudaMemcpyDeviceToHost));
  gpuErrchk(  cudaMemcpy(_h_d_elements,  _d_elements,    size1*sizeof(int),     cudaMemcpyDeviceToHost));
  gpuErrchk( cudaMemcpy(_h_d_inv_elements, _d_inv_elements,size2*sizeof(int),  cudaMemcpyDeviceToHost));
  gpuErrchk( cudaMemcpy(_h_d_gauss_weight, _d_gauss_weight,_NINT_PER_ELEMENT*sizeof(double),cudaMemcpyDeviceToHost));
  gpuErrchk( cudaMemcpy(_h_EPOT, _d_EPOT,_NELE*sizeof(double),cudaMemcpyDeviceToHost));

  for (int i = 0;i<2;i++)   assert(map23[i]==_h_d_map23[i]);
  for (int i = 0;i<_NP;i++){
  //	  printf("_SR[i]=%g,%g,%g, _h_d_SR[i]=%g,%g,%g, diff = %g,%g,%g\n", _SR[i][0], _SR[i][1],  _SR[i][2], _h_d_SR[i][0], _h_d_SR[i][1], _h_d_SR[i][2], _SR[i][0]-_h_d_SR[i][0],  _SR[i][1]-_h_d_SR[i][1], _SR[i][2]-_h_d_SR[i][2]);
    assert((G_Vector3(_SR[i])==G_Vector3(_h_d_SR[i])));
  }
  for (int i = 0;i<_NP;i++)   assert((G_Vector3(_SRref[i])==G_Vector3(_h_d_SRref[i])));
  for (int i = 0;i<_NP;i++)   assert((G_Vector3(_Rref[i])==G_Vector3(_h_d_Rref[i])));
  for (int i = 0;i<_NP;i++)   assert((G_Vector3(_F[i])==G_Vector3(_h_d_F[i])));
  for (int i = 0;i<_NP;i++)   assert(fixed[i]==_h_d_fixed[i]);
  for (int i = 0;i<size1;i++) assert(elements[i]==_h_d_elements[i]);
  for (int i = 0;i<size2;i++) assert(inv_elements[i]==_h_d_inv_elements[i]);  
  for (int i = 0;i<_NINT_PER_ELEMENT;i++) assert(fabs(gauss_weight[i]-_h_d_gauss_weight[i])<1e-15);

  INFO_Printf("I am about to get out of check_host_device memory transfer\n");
  return 0;
}
#endif
