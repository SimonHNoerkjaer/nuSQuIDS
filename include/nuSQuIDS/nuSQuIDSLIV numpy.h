#ifndef nusquidslv_H
#define nusquidslv_H

#include <vector>
#include <iostream>
#include <string>
#include <nuSQuIDS/nuSQuIDS.h>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <Python.h>
#include <numpy/arrayobject.h>



// Printing matrix function
void print_gsl_matrix(gsl_matrix_complex* matrix) {
  for (size_t i = 0; i < matrix->size1; i++) {
    for (size_t j = 0; j < matrix->size2; j++) {
      //std::cout << std::fixed << std::setprecision(2) << GSL_REAL(gsl_matrix_complex_get(matrix, i, j)) << " + " << GSL_IMAG(gsl_matrix_complex_get(matrix, i, j)) << "i   "; 
      std::cout << std::scientific << GSL_REAL(gsl_matrix_complex_get(matrix, i, j)) << " + " << GSL_IMAG(gsl_matrix_complex_get(matrix, i, j)) << "i   "; 
    }
    std::cout << std::endl;
  }
}







namespace nusquids {

class nuSQUIDSLIV: public nuSQUIDS {

  private:

    // Sidereal vs isotropic
    bool directional = false;

    // The parameters of the LIV model
    int LIV_n;
    double LIV_cft;


    double ra_rad;
    double dec_rad;
    double sme_a_x;
    double sme_a_y;
    double sme_a_z;
    double sme_c_tx;
    double sme_c_ty;
    double sme_c_tz;

    
    


    //
    //Coordinates
    //

    // celestial colatitude and longitude
    double theta = M_pi/2 + dec_rad;
    double phi = ra 

    
    // spherical coordinates unit propagation vectors
    
    // r vector
    double NX = sin(theta) * cos(phi)
    double NY = sin(theta) * sin(phi)
    double = cos(theta)
    
    // theta vector
    double ThetaX = cos(theta) * cos(phi)
    double ThetaY = cos(theta) * sin(phi)
    double ThetaZ = -sin(theta)  
    
    // phi vector
    double PhiX = -sin(phi)
    double PhiY = cos(phi)
    double PhiZ = 0

    //
    //  Mass independent operators
    //
    
    // Amplitude to be multiplied with sin(omega_sid L)
    double As0 = -NY * a_eV_x + NX * a_eV_y
    double As1 = + 2 * NY * c_tx - 2 * NX * c_ty
    double As = As0 + E * As1

    // Amplitude to be multiplied with cos(omega_sid L)
    double Ac0 = - NX * a_eV_x - NY * a_eV_y
    double Ac1 = 2 * NX * c_tx + 2 * NY * c_ty
    double Ac = Ac0 + E * Ac1



    // double LIV_a_eV;
    // double LIV_c;

//TODO: Fix if-statement for sidereal
    // The LIV potential
    // void def_LIVP() {
    //   if (directional) {
    //     squids::SU_vector LIVPX;
    //     squids::SU_vector LIVPY;
    //     squids::SU_vector LIVPZ;
    //     std::vector<squids::SU_vector> LIVPX_evol;
    //     std::vector<squids::SU_vector> LIVPY_evol;
    //     std::vector<squids::SU_vector> LIVPZ_evol;
    //   } else {
    //     squids::SU_vector LIVP;
    //     std::vector<squids::SU_vector> LIVP_evol;
    //   }
    // }  

    squids::SU_vector LIVP;
    std::vector<squids::SU_vector> LIVP_evol;



    void AddToPreDerive(double x){
      for(unsigned int ei = 0; ei < ne; ei++){
        // asumming same mass hamiltonian for neutrinos/antineutrinos
        squids::SU_vector h0 = H0( E_range[ei], 0 );
        // LIVP_evol[ei] = LIVP.Evolve( h0, (x-Get_t_initial()) );

        // //TODO: Fix if-statement for sidereal
        if (directional) {
          std::cout << "Sidereal not yet implemented. Changed to isotropic" << std::endl;
        //   LIVPX_evol[ei] = LIVPX.Evolve( h0, (x-Get_t_initial()) );
        //   LIVPY_evol[ei] = LIVPY.Evolve( h0, (x-Get_t_initial()) );
        //   LIVPZ_evol[ei] = LIVPZ.Evolve( h0, (x-Get_t_initial()) );
          LIVP_evol[ei] = LIVP.Evolve( h0, (x-Get_t_initial()) );
        } else {
          LIVP_evol[ei] = LIVP.Evolve( h0, (x-Get_t_initial()) );
        }

        
      }
    }


    // squids::SU_vector H0(double E, unsigned int irho) const {
    //   squids::SU_vector h0 = nuSQUIDS::H0( E, irho );
    //   squids::SU_vector h0_liv = h0 + ( LIV_cft * pow(E, LIV_n) );
    //   return h0_liv;
    // }


    squids:: SU_vector HI(unsigned int ie,unsigned int irho) const {
      //std::cout << "+++ HI : start" << std::endl;
      squids::SU_vector potential = nuSQUIDS::HI(ie, irho);
      double sign = 1;
      // if ((irho == 1 and NT==both) or NT==antineutrino){
      //     // antineutrino matter potential flips sign
      //     sign*=(-1);
      // }
      // ================= HERE WE ADD THE NEW PHYSICS ===================
      potential += sign*pow(E_range[ie], LIV_n) * LIVP_evol[ie]; // <- super important line here is where all the physics is set
      // ================= HERE WE ADD THE NEW PHYSICS ===================
      //std::cout << "+++ HI : end" << std::endl;
      return potential;
    }

  void init() {
    /*
      Common init function, called by various constructors
    */

    //std::cout << "+++ nuSQUIDSLIV::init" << std::endl;

    // Init LIV params to null values
    LIV_n = 0;
    LIV_cft = 0;

    // Allocate some matrices
    LIVP_evol.resize(ne);
    for(unsigned int ei = 0; ei < ne; ei++){
      LIVP_evol[ei] = squids::SU_vector(nsun);
    }

  }


  public:


    //
    // Constructors
    //

    //Default void constructor (never used but required by boost python, specifically RegisterBasicNuSQuIDSPythonBindings, not sure why).
    nuSQUIDSLIV() {}

    /// Multiple energy mode constructor. Basically just passing down to base class constructor and calling init.
    nuSQUIDSLIV(marray<double,1> E_vector,unsigned int numneu, NeutrinoType NT = both, bool iinteraction = false, std::shared_ptr<CrossSectionLibrary> ncs = nullptr)
      : nuSQUIDS(E_vector, numneu, NT, iinteraction, ncs)
      // : nuSQUIDS(E_vector, numneu, NT, false, nullptr)
    { init(); }

    /// Single energy mode constructor. Basically just passing down to base class constructor and calling init.
    nuSQUIDSLIV(unsigned int numneu, NeutrinoType NT = neutrino)
      : nuSQUIDS(numneu, NT)
    { init(); }

    // Load from serialized file (not yet implemented, but required to exist for pybindings)
    nuSQUIDSLIV(std::string hdf5_filename, std::string grp = "/", std::shared_ptr<InteractionStructure> int_struct = nullptr)
      // : nuSQUIDS(hdf5_filename, grp, int_struct)
    {
      throw std::runtime_error("nuSQUIDSLIV::Error::HDF5 serialization not implemented");
    }






// 
//  FUNCTIONS FOR IMPORTING NUMPY ARRAYS AS GSL MATRICES
// 

    void convert_numpy_to_gsl(
      PyObject* a_eV_x_np,
      PyObject* a_eV_y_np,
      PyObject* a_eV_z_np,
      PyObject* c_tx_np,
      PyObject* c_ty_np,
      PyObject* c_tz_np) 
      {

      Py_Initialize();                                           // Initialize the Python interpreter
      PyObject* numpy = PyImport_ImportModule("numpy");          // Import the NumPy module
      import_array();                                            // Initialize the NumPy C API


      // Convert NumPy arrays to GSL matrices
      gsl_matrix_complex* a_eV_x = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* a_eV_y = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* a_eV_z = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* c_tx = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* c_ty = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* c_tz = gsl_matrix_complex_calloc(3, 3);

      for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            gsl_complex a_eV_x_value = gsl_complex_rect(
                PyArray_GETPTR2(a_eV_x_np, i, j)->real,
                PyArray_GETPTR2(a_eV_x_np, i, j)->imag
            );
            gsl_matrix_complex_set(a_eV_x, i, j, a_eV_x_value);

            gsl_complex a_eV_y_value = gsl_complex_rect(
                PyArray_GETPTR2(a_eV_y_np, i, j)->real,
                PyArray_GETPTR2(a_eV_y_np, i, j)->imag
            );
            gsl_matrix_complex_set(a_eV_y, i, j, a_eV_y_value);

            gsl_complex a_eV_z_value = gsl_complex_rect(
                PyArray_GETPTR2(a_eV_z_np, i, j)->real,
                PyArray_GETPTR2(a_eV_z_np, i, j)->imag
            );
            gsl_matrix_complex_set(a_eV_z, i, j, a_eV_z_value);

            gsl_complex c_tx_value = gsl_complex_rect(
                PyArray_GETPTR2(c_tx_np, i, j)->real,
                PyArray_GETPTR2(c_tx_np, i, j)->imag
            );
            gsl_matrix_complex_set(c_tx, i, j, c_tx_value);

            gsl_complex c_ty_value = gsl_complex_rect(
                PyArray_GETPTR2(c_ty_np, i, j)->real,
                PyArray_GETPTR2(c_ty_np, i, j)->imag
            );
            gsl_matrix_complex_set(c_ty, i, j, c_ty_value);

            gsl_complex c_tz_value = gsl_complex_rect(
                PyArray_GETPTR2(c_tz_np, i, j)->real,
                PyArray_GETPTR2(c_tz_np, i, j)->imag
            );
            gsl_matrix_complex_set(c_tz, i, j, c_tz_value);
        }
    }

      // Finalize the Python interpreter (optional)
      Py_Finalize();
    return 0;
}












    //
    // FUNCTION FOR SETTING SME param getters/setters
    //

    void  Set_Directional_SME(
      double ra_rad,
      double dec_rad,
      numpy_matrix a_eV_x,
      numpy_matrix a_eV_y,
      numpy_matrix a_eV_z,
      numpy_matrix c_tx,
      numpy_matrix c_ty,
      numpy_matrix c_tz,
      double E){
        
      //
      // Convert numpy to GSL matrices 
      //

      gsl_matrix_complex* a_eV_x = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* a_eV_y = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* a_eV_z = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* c_tx = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* c_ty = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* c_tz = gsl_matrix_complex_calloc(3, 3);

      // set gsl matrices to numpy matrices
      // for (size_t i = 0; i < 3; i++) {
      //   for (size_t j = 0; j < 3; j++) {
      //     gsl_matrix_complex_set(a_eV_x, i, j, a_eV_x[i][j]);
      //     gsl_matrix_complex_set(a_eV_y, i, j, a_eV_y[i][j]);
      //     gsl_matrix_complex_set(a_eV_z, i, j, a_eV_z[i][j]);
      //     gsl_matrix_complex_set(c_tx, i, j, c_tx[i][j]);
      //     gsl_matrix_complex_set(c_ty, i, j, c_ty[i][j]);
      //     gsl_matrix_complex_set(c_tz, i, j, c_tz[i][j]);
      //   }
      // }

      // Create copies of matrices to be used in calculations
      gsl_matrix_complex* a_eV_x_copy1 = gsl_matrix_complex_alloc(3, 3);
      gsl_matrix_complex* a_eV_x_copy2 = gsl_matrix_complex_alloc(3, 3);
      gsl_matrix_complex_memcpy(a_eV_x_copy1, a_eV_x);
      gsl_matrix_complex_memcpy(a_eV_x_copy2, a_eV_x);
      gsl_matrix_complex* a_eV_y_copy1 = gsl_matrix_complex_alloc(3, 3);
      gsl_matrix_complex* a_eV_y_copy2 = gsl_matrix_complex_alloc(3, 3);
      gsl_matrix_complex_memcpy(a_eV_y_copy1, a_eV_y);
      gsl_matrix_complex_memcpy(a_eV_y_copy2, a_eV_y);
      gsl_matrix_complex* a_eV_z_copy1 = gsl_matrix_complex_alloc(3, 3);
      gsl_matrix_complex* a_eV_z_copy2 = gsl_matrix_complex_alloc(3, 3);
      gsl_matrix_complex_memcpy(a_eV_z_copy1, a_eV_z);
      gsl_matrix_complex_memcpy(a_eV_z_copy2, a_eV_z);
      gsl_matrix_complex* c_tx_copy1 = gsl_matrix_complex_alloc(3, 3);
      gsl_matrix_complex* c_tx_copy2 = gsl_matrix_complex_alloc(3, 3);
      gsl_matrix_complex_memcpy(c_tx_copy1, c_tx);
      gsl_matrix_complex_memcpy(c_tx_copy2, c_tx);
      gsl_matrix_complex* c_ty_copy1 = gsl_matrix_complex_alloc(3, 3);
      gsl_matrix_complex* c_ty_copy2 = gsl_matrix_complex_alloc(3, 3);
      gsl_matrix_complex_memcpy(c_ty_copy1, c_ty);
      gsl_matrix_complex_memcpy(c_ty_copy2, c_ty);
      gsl_matrix_complex* c_tz_copy1 = gsl_matrix_complex_alloc(3, 3);
      gsl_matrix_complex* c_tz_copy2 = gsl_matrix_complex_alloc(3, 3);
      gsl_matrix_complex_memcpy(c_tz_copy1, c_tz);
      gsl_matrix_complex_memcpy(c_tz_copy2, c_tz);

      



      //
      //Coordinates
      //

      // celestial colatitude and longitude
      double theta = M_pi/2 + dec_rad;
      double phi = ra;
      // spherical coordinates unit propagation vectors
      // r vector
      double NX = sin(theta) * cos(phi);
      double NY = sin(theta) * sin(phi);
      double NZ = cos(theta);
      


      //
      //  Mass independent operators
      //
      
      // // Amplitude to be multiplied with sin(omega_sid L)
      // As0 = -NY * a_eV_x + NX * a_eV_y;
      // As1 = + 2 * NY * c_tx - 2 * NX * c_ty;
      // As = As0 + E * As1;

      gsl_matrix_complex* As0 = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* As1 = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* As = gsl_matrix_complex_calloc(3, 3);

      // Perform the operation: As0 = -NY * a_eV_x + NX * a_eV_y
      gsl_matrix_complex_scale(a_eV_x_copy1, gsl_complex_rect(-NY, 0));  // Scale a_eV_x by -NY
      gsl_matrix_complex_scale(a_eV_y_copy1, gsl_complex_rect(NX, 0));  // Scale a_eV_y by NX
      gsl_matrix_complex_add(As0, a_eV_x_copy1, a_eV_y_copy1);  // Add scaled matrices to get As0

      // Perform the operation: As1 = + 2 * NY * c_tx - 2 * NX * c_ty
      gsl_matrix_complex_scale(c_tx_copy1, gsl_complex_rect(2 * NY, 0));  // Scale c_tx by 2 * NY
      gsl_matrix_complex_scale(c_ty_copy1, gsl_complex_rect(-2 * NX, 0));  // Scale c_ty by -2 * NX
      gsl_matrix_complex_add(As1, c_tx_copy1, c_ty_copy1);  // Add scaled matrices to get As1

      // Perform the operation: As = As0 + E * As1
      gsl_matrix_complex_scale(As1, gsl_complex_rect(E, 0));  // Scale As1 by E
      gsl_matrix_complex_add(As, As0, As1);  // Add matrices to get As


  



     
      // Amplitude to be multiplied with cos(omega_sid L)
      // Ac0 = - NX * a_eV_x - NY * a_eV_y;
      // Ac1 = 2 * NX * c_tx + 2 * NY * c_ty;
      // Ac = Ac0 + E * Ac1;

      gsl_matrix_complex* Ac0 = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* Ac1 = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* Ac = gsl_matrix_complex_calloc(3, 3);

      // Perform the operation: Ac0 = - NX * a_eV_x - NY * a_eV_y
      gsl_matrix_complex_scale(a_eV_x_copy2, gsl_complex_rect(-NX, 0));  // Scale a_eV_x by -NX
      gsl_matrix_complex_scale(a_eV_y_copy2, gsl_complex_rect(-NY, 0));  // Scale a_eV_y by -NY
      gsl_matrix_complex_add(Ac0, a_eV_x_copy2, a_eV_y_copy2);  // Add scaled matrices to get Ac0

      // Perform the operation: Ac1 = 2 * NX * c_tx + 2 * NY * c_ty
      gsl_matrix_complex_scale(c_tx_copy2, gsl_complex_rect(2 * NX, 0));  // Scale c_tx by 2 * NX
      gsl_matrix_complex_scale(c_ty_copy2, gsl_complex_rect(2 * NY, 0));  // Scale c_ty by 2 * NY
      gsl_matrix_complex_add(Ac1, c_tx_copy2, c_ty_copy2);  // Add scaled matrices to get Ac1

      // Perform the operation: Ac = Ac0 + E * Ac1
      gsl_matrix_complex_scale(Ac1, gsl_complex_rect(E, 0));  // Scale Ac1 by E
      gsl_matrix_complex_add(Ac, Ac0, Ac1);  // Add matrices to get Ac


      

      // write as gsl: Const = NZ * a_eV_z  
      gsl_matrix_complex* Const = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex_scale(a_eV_z, gsl_complex_rect(NZ, 0));  // Scale a_eV_z by NZ
      gsl_matrix_complex_memcpy(Const, a_eV_z);  // Copy a_eV_z to Const
      gsl_matrix_complex_free(a_eV_z);

      // write as gsl: H_eff = As + Const
      gsl_matrix_complex* H_eff = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex_add(H_eff, As, Const);  // Add matrices to get H_eff
      gsl_matrix_complex_free(As);
      gsl_matrix_complex_free(Const);

      
      LIVP = squids::SU_vector(H_eff);
      gsl_matrix_complex_free(H_eff);



      // free allocated matrices
      gsl_matrix_complex_free(a_eV_x_copy1);
      gsl_matrix_complex_free(a_eV_x_copy2);
      gsl_matrix_complex_free(a_eV_y_copy1);
      gsl_matrix_complex_free(a_eV_y_copy2);
      gsl_matrix_complex_free(a_eV_z_copy1);
      gsl_matrix_complex_free(a_eV_z_copy2);
      gsl_matrix_complex_free(c_tx_copy1);
      gsl_matrix_complex_free(c_tx_copy2);
      gsl_matrix_complex_free(c_ty_copy1);
      gsl_matrix_complex_free(c_ty_copy2);
      gsl_matrix_complex_free(c_tz_copy1);
      gsl_matrix_complex_free(c_tz_copy2);
      gsl_matrix_complex_free(a_eV_x);
      gsl_matrix_complex_free(a_eV_y);
      gsl_matrix_complex_free(a_eV_z);
      gsl_matrix_complex_free(c_tx);
      gsl_matrix_complex_free(c_ty);
      gsl_matrix_complex_free(c_tz);

    }





    void Set_LIVCoefficient(double cft){
      
      
      gsl_complex c{{ cft , 0.0 }}; //Only using real part right now


       // defining a complex matrix M which will contain our flavor
       // violating flavor structure.
       gsl_matrix_complex * M = gsl_matrix_complex_calloc(3,3); //TODO check num nu 
       gsl_matrix_complex_set(M, 2, 2, c);


       LIVP = squids::SU_vector(M);

       std::cout << "Matrix:" << std::endl;
       print_gsl_matrix(M);

       // rotate from flavor to mass basis
       // LIVP.RotateToB1(params);

       // free allocated matrix
       gsl_matrix_complex_free(M);

    }

    void Set_LIVEnergyDependence(int n){
      LIV_n = n;
    }

};

} // close nusquids namespace

#endif //nusquidslv_h

