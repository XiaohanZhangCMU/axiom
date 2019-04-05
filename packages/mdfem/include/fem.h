/*
  fem.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Fri Mar 11 21:07:00 2016

  FUNCTION  :  FEM model of hyperelastic material with brick elements
  NOTE :       All elements must be the same shape and size for now.
*/

#ifndef _FEM_H
#define _FEM_H
#include <fstream>      // std::ofstream
#include <cassert>
#include "md.h"
#include <iostream>
#include <vector>

//#define _USECUDA
#ifdef _USECUDA
//#define DEBUG_USECUDA
#include <cuda_runtime.h>
#include "linalg3_cu.h"
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code,const char *file,int line,bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %dn",cudaGetErrorString(code),file,line);
    if (abort) { getchar(); exit(code); } } } 
#endif
#define _MAX_NELEM_SHARE_NODE 10

class FEMFrame : public virtual MDFrame /* FEM with brick elements */
{
 protected:
    /* FEM parameters */
    int _NELE;              /* number of elements */
    int _NNODE_PER_ELEMENT; /* number of nodes per element */
    int _NINT_PER_ELEMENT;  /* number of Gauss integration points per element */
    int _NFACE_PER_ELEMENT;
    int _NDIM;              /* dimension of elements */
    int _EquationType;
    int *elements;
    int *inv_elements;
    int *colorids;
    int *map23;

#ifdef _USECUDA
    int *_d_elements;
    int *_d_inv_elements;
    int *_d_map23;
    int *_d_fixed;
    int *_d_colorids;

    double *_d_VIRIAL_IND;
    double * _d_H_element;
    double *_d_gauss_weight; 
    double *_d_dFdu;        
    double *_d_EPOT;
    double *_d_EPOT_IND;
    double* _d_EPOT_RMV;
    double *_h_EPOT;        /* used for host reduction only */
    G_Vector3 *_d_SR;
    G_Vector3 *_d_SRref;       
    G_Vector3 *_d_Rref;        
    G_Vector3 *_d_F;
    G_Vector3 *_d_F_padding; /* padding forces to avoid race problem */
#endif

#ifdef DEBUG_USECUDA
    int *_h_d_elements;
    int *_h_d_inv_elements;
    int *_h_d_map23;
    int *_h_d_fixed;
    int *_h_d_colorids;
    class G_Matrix33 *_h_d_VIRIAL_IND;
    double *_h_d_H_element;
    double *_h_d_gauss_weight; 
    double *_h_d_dFdu;        
    double *_h_d_EPOT;
    double *_h_d_EPOT_IND;
    double *_h_d_EPOT_RMV;
    G_Vector3 *_h_d_SR;
    G_Vector3 *_h_d_SRref;       
    G_Vector3 *_h_d_Rref;        
    G_Vector3 *_h_d_F;
    G_Vector3 *_h_d_F_padding; //for padding forces of elements
    int check_host_device_memory_transfer(); 
#endif

    char ELEMENT_TYPE[100], ELEMENT_INT_TYPE[100];
    char elements_file[MAXFILENAMELEN];
    char fem_coeff_file[MAXFILENAMELEN];

    double y_eigen_strain;
    double x_eigen_strain;
    double x_eigen_zbound_min;
    double x_eigen_zbound_max;
    double y_eigen_zbound_min;
    double y_eigen_zbound_max;

    double __V0;
    double *gauss_weight;   /* weight of each Gauss integration point */
    double *dFdu;           /* derivative of F (deformation gradient) wrt nodal displacement */
    Vector3 *_SRref;        /* reference nodal position (scaled coordinates) */
    Vector3 *_Rref;         /* reference nodal position (real coordinates) */

public:
    FEMFrame():_NELE(0),_NNODE_PER_ELEMENT(0),_NINT_PER_ELEMENT(0),
      elements(0),gauss_weight(0),dFdu(0),_SRref(0),_Rref(0), n_bdy_nodes(0),n_bdy_tags(0), n_existed_atoms(_NP), nfemNodes(0), _ReadColorID(0) {};
    virtual ~FEMFrame() {
      delete(win);
#ifdef _USECUDA
     free_device_ptr();
#endif

    } 
    void RtoRref();
    int read_elements(const char*);
    int read_fem_coeff(const char*);
    int read_uniform_fem_coeff(const char*);
    int read_element_wise_fem_coeff(const char*, int);
    int read_1stfem_to_Alloc(const char*);
    Matrix33 getEigenF(Vector3 p, Matrix33 Fdef);

    void WriteStressCoord();

    void beam_fem_energy_force();
    void create_inverse_connectivities_matrix();
    void islam_fem_energy_force();
    void snap_fem_energy_force();

    virtual void potential();
    void Alloc_Elements();
    void Alloc_Element_Coeff();
    
    virtual void plot();
    virtual void Alloc();

    
public:
  char fem_bdy_nodes_file[MAXFILENAMELEN];
  char contact_file[MAXFILENAMELEN];
  int *bNds, *bTags;

  int *cNds;
  double *bXs, *bYs, *bZs;
  int n_bdy_nodes;
  int n_bdy_tags;  
  int n_existed_atoms;
  int nfemNodes;

  int n_contact_nodes;

  int read_bdy_nodes(const char*);
  int read_contact(const char*);
  int _ReadColorID;
#ifdef _PARALLEL
  void Broadcast_FEM_Param();
#endif

  int read2cn();
  void shift_fem_node_id(int);

#ifdef _USECUDA
  int test_saxpy();
  void cuda_memcpy_all(void);
  void cuda_memory_alloc(void);
  void cuda_memory_alloc_elements(void);
  void cuda_memory_alloc_element_coeff(void);
  void free_device_ptr(void);
  void cuda_snap_fem_energy_force();
  void cuda_beam_fem_energy_force();
#endif

/* Expose fields to python through pybind11 */
  int test_function(int m);

  std::string get_contact_file();
  void set_contact_file(std::string);
  std::vector<int> get_bNds();
  void set_bNds(std::vector<int>);
  std::vector<int> get_bTags();
  void set_bTags(std::vector<int>);


};

#endif // _FEM_H

